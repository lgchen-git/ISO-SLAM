#!/usr/bin/env python3
"""Shared long-read gene assignment helpers for LR-SLAM callers."""

from __future__ import annotations

import bisect
import gzip
from collections import defaultdict
from dataclasses import dataclass

import pysam


BIN_SIZE = 20000


@dataclass
class TranscriptModel:
    transcript_id: str
    gene_id: str
    gene_name: str
    chrom: str
    strand: str
    start: int
    end: int
    exons: list[tuple[int, int]]
    introns: list[tuple[int, int]]


@dataclass
class AssignmentResult:
    gene_id: str
    gene_name: str
    transcript_id: str
    strand: str
    region_class: str
    assignment_class: str
    assignment_level: str
    assignment_support_tx_count: int
    assignment_support_gene_count: int
    assignment_resolution_method: str


def attrs_to_dict(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in text.rstrip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        pieces = part.split(" ", 1)
        if len(pieces) != 2:
            continue
        out[pieces[0]] = pieces[1].strip().strip('"')
    return out


def load_interval_annotations(gtf: str, feature_type: str):
    annotations = defaultdict(lambda: defaultdict(list))
    opener = gzip.open if gtf.endswith(".gz") else open
    with opener(gtf, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != feature_type:
                continue
            attrs = attrs_to_dict(fields[8])
            gene_id = attrs.get("gene_id", ".").split(".")[0]
            gene_name = attrs.get("gene_name", gene_id)
            interval = {
                "start": int(fields[3]) - 1,
                "end": int(fields[4]),
                "gene_id": gene_id,
                "gene_name": gene_name,
                "strand": fields[6],
            }
            for bin_id in range(interval["start"] // BIN_SIZE, (interval["end"] - 1) // BIN_SIZE + 1):
                annotations[fields[0]][bin_id].append(interval)
    return annotations


def load_exon_annotations(gtf: str):
    return load_interval_annotations(gtf, "exon")


def load_gene_annotations(gtf: str):
    return load_interval_annotations(gtf, "gene")


def build_introns(exons: list[tuple[int, int]]) -> list[tuple[int, int]]:
    introns = []
    for idx in range(len(exons) - 1):
        intron = (exons[idx][1], exons[idx + 1][0])
        if intron[0] < intron[1]:
            introns.append(intron)
    return introns


def load_transcript_models(gtf: str) -> tuple[dict[str, TranscriptModel], dict[str, list[TranscriptModel]], dict[str, list[int]]]:
    transcript_rows: dict[str, dict[str, object]] = {}
    opener = gzip.open if gtf.endswith(".gz") else open
    with opener(gtf, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            attrs = attrs_to_dict(fields[8])
            transcript_id = attrs.get("transcript_id", "").split(".")[0]
            gene_id = attrs.get("gene_id", ".").split(".")[0]
            if not transcript_id or gene_id == ".":
                continue
            gene_name = attrs.get("gene_name", gene_id)
            entry = transcript_rows.setdefault(
                transcript_id,
                {
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "chrom": fields[0],
                    "strand": fields[6],
                    "exons": [],
                },
            )
            entry["exons"].append((int(fields[3]) - 1, int(fields[4])))
    transcripts: dict[str, TranscriptModel] = {}
    chrom_models: dict[str, list[TranscriptModel]] = defaultdict(list)
    for transcript_id, entry in transcript_rows.items():
        exons = sorted(set(entry["exons"]))
        if not exons:
            continue
        model = TranscriptModel(
            transcript_id=transcript_id,
            gene_id=str(entry["gene_id"]),
            gene_name=str(entry["gene_name"]),
            chrom=str(entry["chrom"]),
            strand=str(entry["strand"]),
            start=exons[0][0],
            end=exons[-1][1],
            exons=exons,
            introns=build_introns(exons),
        )
        transcripts[transcript_id] = model
        chrom_models[model.chrom].append(model)
    model_lists: dict[str, list[TranscriptModel]] = {}
    start_lists: dict[str, list[int]] = {}
    for chrom, models in chrom_models.items():
        sorted_models = sorted(models, key=lambda m: m.start)
        model_lists[chrom] = sorted_models
        start_lists[chrom] = [m.start for m in sorted_models]
    return transcripts, model_lists, start_lists


def overlap_score(blocks, interval) -> int:
    score = 0
    for start, end in blocks:
        ov_start = max(start, interval["start"])
        ov_end = min(end, interval["end"])
        if ov_end > ov_start:
            score += ov_end - ov_start
    return score


def assign_best_interval(chrom: str, blocks, annotations):
    chrom_bins = annotations.get(chrom)
    if not chrom_bins:
        return None, False
    candidates = {}
    for start, end in blocks:
        for bin_id in range(start // BIN_SIZE, (end - 1) // BIN_SIZE + 1):
            for interval in chrom_bins.get(bin_id, []):
                candidates[id(interval)] = interval
    scores = defaultdict(lambda: {"score": 0, "gene_name": "", "strand": ""})
    for interval in candidates.values():
        score = overlap_score(blocks, interval)
        if score <= 0:
            continue
        gene = scores[interval["gene_id"]]
        gene["score"] += score
        gene["gene_name"] = interval["gene_name"]
        gene["strand"] = interval["strand"]
    if not scores:
        return None, False
    best_gene = ""
    best_payload = None
    best_score = -1
    second_score = -1
    for gene_id, payload in scores.items():
        score = payload["score"]
        if score > best_score:
            second_score = best_score
            best_gene = gene_id
            best_payload = payload
            best_score = score
        elif score > second_score:
            second_score = score
    tied = second_score == best_score
    return (best_gene, best_payload["gene_name"], best_payload["strand"]), tied


def block_overlap(blocks: list[tuple[int, int]], intervals: list[tuple[int, int]]) -> int:
    score = 0
    for start, end in blocks:
        for interval_start, interval_end in intervals:
            overlap_start = max(start, interval_start)
            overlap_end = min(end, interval_end)
            if overlap_end > overlap_start:
                score += overlap_end - overlap_start
    return score


def read_strand_from_alignment(read: pysam.AlignedSegment) -> str:
    return "-" if read.is_reverse else "+"


def effective_read_strand(read: pysam.AlignedSegment, library_strand_mode: str) -> str:
    """Return the transcript strand implied by the configured library model."""
    observed = read_strand_from_alignment(read)
    if library_strand_mode == "reverse":
        return "+" if observed == "-" else "-"
    return observed


def requires_strand_filter(args) -> bool:
    mode = getattr(args, "library_strand_mode", "legacy")
    if mode in {"forward", "reverse"}:
        return True
    if mode == "unstranded":
        return False
    return bool(getattr(args, "assignment_require_same_strand", False))


def collect_overlapping_transcripts(
    chrom: str, blocks, model_lists: dict[str, list[TranscriptModel]], start_lists: dict[str, list[int]]
) -> list[TranscriptModel]:
    models = model_lists.get(chrom)
    starts = start_lists.get(chrom)
    if not models:
        return []
    read_start = blocks[0][0]
    read_end = blocks[-1][1]
    idx = bisect.bisect_right(starts, read_end)
    candidates: list[TranscriptModel] = []
    for i in range(idx - 1, -1, -1):
        model = models[i]
        if model.end < read_start:
            break
        candidates.append(model)
    return candidates


def transcript_compatibility(
    blocks: list[tuple[int, int]],
    read_introns: list[tuple[int, int]],
    model: TranscriptModel,
    min_exon_overlap_bp: int,
) -> tuple[bool, int, int, int]:
    exon_overlap = block_overlap(blocks, model.exons)
    if exon_overlap < min_exon_overlap_bp:
        return False, exon_overlap, 0, 0
    if read_introns:
        transcript_introns = set(model.introns)
        matched_introns = sum(1 for intron in read_introns if intron in transcript_introns)
        if matched_introns != len(read_introns):
            return False, exon_overlap, matched_introns, len(read_introns)
    else:
        matched_introns = 0
    end_distance = abs(blocks[0][0] - model.start) + abs(blocks[-1][1] - model.end)
    return True, exon_overlap, matched_introns, end_distance


def assign_transcript_aware(
    read: pysam.AlignedSegment,
    gene_annotations,
    model_lists: dict[str, list[TranscriptModel]],
    start_lists: dict[str, list[int]],
    args,
) -> AssignmentResult:
    blocks = read.get_blocks()
    chrom = read.reference_name
    library_strand_mode = getattr(args, "library_strand_mode", "legacy")
    observed_strand = effective_read_strand(read, library_strand_mode)
    strand_filter = requires_strand_filter(args)
    read_introns = build_introns(blocks)
    overlapping_transcripts = collect_overlapping_transcripts(chrom, blocks, model_lists, start_lists)
    if not overlapping_transcripts:
        return AssignmentResult(".", ".", ".", ".", "intergenic", "intergenic", "unresolved", 0, 0, "no_overlap")

    same_strand_models = [model for model in overlapping_transcripts if model.strand == observed_strand]
    opposite_strand_models = [model for model in overlapping_transcripts if model.strand != observed_strand]
    candidate_models = same_strand_models if strand_filter else overlapping_transcripts

    if not candidate_models:
        gene_assigned, gene_tied = assign_best_interval(chrom, blocks, gene_annotations)
        if gene_assigned is not None and not gene_tied:
            gene_id, gene_name, gene_strand = gene_assigned
            return AssignmentResult(gene_id, gene_name, ".", gene_strand, "antisense", "antisense_gene_salvaged", "gene", 0, 1, "strand_gene_salvage")
        if gene_assigned is not None and gene_tied:
            opp_gene_ids = {model.gene_id for model in opposite_strand_models}
            return AssignmentResult(".", ".", ".", ".", "antisense", "shared_antisense_ambiguous", "unresolved", 0, len(opp_gene_ids) if opp_gene_ids else 2, "strand_filter_multi_gene")
        return AssignmentResult(".", ".", ".", ".", "antisense_only", "antisense_only", "unresolved", 0, len({model.gene_id for model in opposite_strand_models}), "strand_filter")

    compatible = []
    for model in candidate_models:
        is_compatible, exon_overlap, matched_introns, end_distance = transcript_compatibility(
            blocks,
            read_introns,
            model,
            args.assignment_min_exon_overlap_bp,
        )
        if is_compatible:
            compatible.append(
                {
                    "model": model,
                    "exon_overlap": exon_overlap,
                    "matched_introns": matched_introns,
                    "end_distance": end_distance,
                    "score": (matched_introns, exon_overlap),
                }
            )
    if compatible:
        best_score = max(item["score"] for item in compatible)
        best = [item for item in compatible if item["score"] == best_score]
        gene_ids = {item["model"].gene_id for item in best}
        if len(best) == 1:
            model = best[0]["model"]
            return AssignmentResult(model.gene_id, model.gene_name, model.transcript_id, model.strand, "exon", "tx_unique", "transcript", 1, 1, "unique_transcript")
        if len(gene_ids) == 1 and args.assignment_allow_gene_level_salvage:
            model = sorted(best, key=lambda item: item["model"].transcript_id)[0]["model"]
            return AssignmentResult(model.gene_id, model.gene_name, ".", model.strand, "exon", "gene_unique_multi_tx", "gene", len(best), 1, "shared_transcripts_same_gene")
        if args.assignment_end_distance_tiebreak:
            best_end_distance = min(item["end_distance"] for item in best)
            end_best = [item for item in best if item["end_distance"] == best_end_distance]
            end_gene_ids = {item["model"].gene_id for item in end_best}
            if len(end_best) == 1 and len(end_gene_ids) == 1 and args.assignment_count_shared_resolved_reads:
                model = end_best[0]["model"]
                return AssignmentResult(model.gene_id, model.gene_name, ".", model.strand, "exon", "shared_gene_resolved_by_ends", "gene", len(best), len(gene_ids), "transcript_ends_tiebreak")
        return AssignmentResult(".", ".", ".", ".", "ambiguous", "shared_gene_ambiguous", "unresolved", len(best), len(gene_ids), "multi_gene_compatible")

    gene_assigned, gene_tied = assign_best_interval(chrom, blocks, gene_annotations)
    if gene_assigned is None:
        if opposite_strand_models and strand_filter:
            return AssignmentResult(".", ".", ".", ".", "antisense_only", "antisense_only", "unresolved", 0, len({model.gene_id for model in opposite_strand_models}), "gene_strand_filter")
        return AssignmentResult(".", ".", ".", ".", "intergenic", "intergenic", "unresolved", 0, 0, "no_gene_overlap")
    gene_id, gene_name, gene_strand = gene_assigned
    if strand_filter and gene_strand != observed_strand:
        return AssignmentResult(gene_id, gene_name, ".", gene_strand, "antisense", "antisense_gene_salvaged", "gene", 0, 1, "gene_strand_salvage")
    if gene_tied:
        return AssignmentResult(".", ".", ".", ".", "ambiguous", "shared_gene_ambiguous", "unresolved", 0, 2, "multi_gene_overlap")
    exon_overlap_same_gene = 0
    for model in same_strand_models if strand_filter else overlapping_transcripts:
        if model.gene_id == gene_id:
            exon_overlap_same_gene = max(exon_overlap_same_gene, block_overlap(blocks, model.exons))
    if exon_overlap_same_gene > 0:
        return AssignmentResult(gene_id, gene_name, ".", gene_strand, "intron", "gene_overlap_unresolved", "gene", 0, 1, "gene_overlap_without_tx_compat")
    return AssignmentResult(gene_id, gene_name, ".", gene_strand, "intron", "intron_only_unique_gene", "gene", 0, 1, "unique_gene_span")


def assign_region(read: pysam.AlignedSegment, exon_annotations, gene_annotations, model_lists, start_lists, args) -> AssignmentResult:
    chrom = read.reference_name
    blocks = read.get_blocks()
    if args.assignment_mode == "transcript_aware_gene":
        return assign_transcript_aware(read, gene_annotations, model_lists, start_lists, args)
    if exon_annotations is None:
        gene_assigned, gene_tied = assign_best_interval(chrom, blocks, gene_annotations)
        if gene_assigned is not None and not gene_tied:
            gene_id, gene_name, strand = gene_assigned
            return AssignmentResult(gene_id, gene_name, ".", strand, "gene", "gene", "gene", 0, 1, "gene_span_overlap")
        if gene_tied:
            return AssignmentResult(".", ".", ".", ".", "ambiguous", "shared_gene_ambiguous", "unresolved", 0, 2, "multi_gene_overlap")
        return AssignmentResult(".", ".", ".", ".", "intergenic", "intergenic", "unresolved", 0, 0, "no_gene_overlap")
    exon_assigned, exon_tied = assign_best_interval(chrom, blocks, exon_annotations)
    if exon_assigned is not None and not exon_tied:
        gene_id, gene_name, strand = exon_assigned
        return AssignmentResult(gene_id, gene_name, ".", strand, "exon", "exon", "gene", 0, 1, "unique_exon_overlap")
    gene_assigned, gene_tied = assign_best_interval(chrom, blocks, gene_annotations)
    if gene_assigned is not None and not gene_tied:
        gene_id, gene_name, strand = gene_assigned
        return AssignmentResult(gene_id, gene_name, ".", strand, "intron", "intron_only_unique_gene", "gene", 0, 1, "unique_gene_span")
    if exon_tied or gene_tied:
        return AssignmentResult(".", ".", ".", ".", "ambiguous", "shared_gene_ambiguous", "unresolved", 0, 2, "multi_gene_overlap")
    return AssignmentResult(".", ".", ".", ".", "intergenic", "intergenic", "unresolved", 0, 0, "no_gene_overlap")
