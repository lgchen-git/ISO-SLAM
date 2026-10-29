#!/usr/bin/env python3
"""Validate and render conserved bioinformatics flow graphs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import os
import re
import sys
from collections import defaultdict, deque
from pathlib import Path
from typing import Any


NODE_FIELDS = ["view_id", "node_id", "stage", "label", "value", "unit", "scope", "node_kind"]
LINK_FIELDS = ["view_id", "link_id", "source", "target", "value", "reason"]
VIEW_FIELDS = ["view_id", "display_name", "order", "is_default"]
TRANSFORM_FIELDS = [
    "transform_id", "source_view_id", "source_node_id", "target_view_id", "target_node_id",
    "relation", "label", "reason",
]
NODE_KINDS = {"retained", "loss", "branch"}
TRANSFORM_RELATIONS = {"one_to_one", "aggregate"}
NODE_LAYOUTS = {"snap", "main_top_losses"}
TRUE_VALUES = {"1", "true", "yes", "y"}
PACKED_PAD_FRACTION = 12.0 / 460.0
# Keep loss nodes close to their source stage so their right-side labels do not
# collide with the retained node in the next stage.
LOSS_GUTTER_FRACTION = 0.12
NPG = [
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488",
    "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
    "#7E6148", "#B09C85", "#6F99AD", "#8A7AA9",
    "#C9A227", "#6E7F80", "#D4537E", "#378ADD",
]


class ContractError(RuntimeError):
    """Raised when a flow violates the declared contract."""


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nodes", required=True, type=Path)
    parser.add_argument("--links", required=True, type=Path)
    parser.add_argument("--views", type=Path)
    parser.add_argument("--transformations", type=Path)
    parser.add_argument("--default-mode", choices=("connected", "view"), default="view")
    parser.add_argument("--node-layout", choices=sorted(NODE_LAYOUTS), default="snap")
    parser.add_argument("--output-dir", type=Path, help="Required unless --validate-only is used")
    parser.add_argument("--title", default="Conserved bioinformatics flow")
    parser.add_argument("--plotly-js", type=Path)
    parser.add_argument("--validate-only", action="store_true")
    return parser.parse_args()


def read_tsv(path: Path, required: list[str]) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise ContractError(f"Missing or empty TSV: {path}")
    with path.open(encoding="utf-8", errors="strict", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ContractError(f"Missing header: {path}")
        missing = [field for field in required if field not in reader.fieldnames]
        if missing:
            raise ContractError(f"Missing columns in {path}: {', '.join(missing)}")
        rows = list(reader)
    if not rows:
        raise ContractError(f"No data rows: {path}")
    return rows


def strict_integer(value: str, field: str, row_name: str) -> int:
    text = str(value).strip()
    if not re.fullmatch(r"[0-9]+", text):
        raise ContractError(f"{row_name}: {field} must be a non-negative integer, got {value!r}")
    return int(text)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def color_for(kind: str, stage_index: int) -> str:
    if kind == "loss":
        return "#D9E0E7"
    if kind == "branch":
        return "#00A087"
    return NPG[stage_index % len(NPG)]


def parse_nodes(path: Path) -> tuple[list[dict[str, Any]], list[str]]:
    raw = read_tsv(path, NODE_FIELDS)
    seen: set[tuple[str, str]] = set()
    view_order: list[str] = []
    nodes: list[dict[str, Any]] = []
    for line, row in enumerate(raw, start=2):
        view_id = row["view_id"].strip()
        node_id = row["node_id"].strip()
        label = row["label"].strip()
        unit = row["unit"].strip()
        scope = row["scope"].strip()
        kind = row["node_kind"].strip()
        name = f"{path}:{line}"
        if not all((view_id, node_id, label, unit, scope)):
            raise ContractError(f"{name}: view_id, node_id, label, unit, and scope must be non-empty")
        key = (view_id, node_id)
        if key in seen:
            raise ContractError(f"{name}: duplicate node identifier {view_id}/{node_id}")
        if kind not in NODE_KINDS:
            raise ContractError(f"{name}: node_kind must be retained, loss, or branch")
        seen.add(key)
        if view_id not in view_order:
            view_order.append(view_id)
        nodes.append({
            "view_id": view_id,
            "node_id": node_id,
            "stage": strict_integer(row["stage"], "stage", name),
            "label": label,
            "value": strict_integer(row["value"], "value", name),
            "unit": unit,
            "scope": scope,
            "node_kind": kind,
        })
    return nodes, view_order


def parse_links(path: Path) -> list[dict[str, Any]]:
    raw = read_tsv(path, LINK_FIELDS)
    seen: set[tuple[str, str]] = set()
    links: list[dict[str, Any]] = []
    for line, row in enumerate(raw, start=2):
        view_id = row["view_id"].strip()
        link_id = row["link_id"].strip()
        source = row["source"].strip()
        target = row["target"].strip()
        reason = row["reason"].strip()
        name = f"{path}:{line}"
        if not all((view_id, link_id, source, target, reason)):
            raise ContractError(f"{name}: identifiers and reason must be non-empty")
        key = (view_id, link_id)
        if key in seen:
            raise ContractError(f"{name}: duplicate link identifier {view_id}/{link_id}")
        seen.add(key)
        links.append({
            "view_id": view_id,
            "link_id": link_id,
            "source": source,
            "target": target,
            "value": strict_integer(row["value"], "value", name),
            "reason": reason,
        })
    return links


def parse_transformations(path: Path | None) -> list[dict[str, str]]:
    if path is None:
        return []
    raw = read_tsv(path, TRANSFORM_FIELDS)
    seen: set[str] = set()
    transformations: list[dict[str, str]] = []
    for line, row in enumerate(raw, start=2):
        name = f"{path}:{line}"
        item = {field: row[field].strip() for field in TRANSFORM_FIELDS}
        if not all(item.values()):
            raise ContractError(f"{name}: transformation fields must be non-empty")
        if item["transform_id"] in seen:
            raise ContractError(f"{name}: duplicate transform_id {item['transform_id']}")
        if item["relation"] not in TRANSFORM_RELATIONS:
            raise ContractError(f"{name}: relation must be one_to_one or aggregate")
        seen.add(item["transform_id"])
        transformations.append(item)
    return transformations


def validate_transformations(
    transformations: list[dict[str, str]], nodes: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    node_lookup = {(node["view_id"], node["node_id"]): node for node in nodes}
    enriched: list[dict[str, Any]] = []
    for item in transformations:
        source_key = (item["source_view_id"], item["source_node_id"])
        target_key = (item["target_view_id"], item["target_node_id"])
        if source_key not in node_lookup or target_key not in node_lookup:
            raise ContractError(f"{item['transform_id']}: transformation endpoint does not exist")
        source = node_lookup[source_key]
        target = node_lookup[target_key]
        if source["unit"] == target["unit"]:
            raise ContractError(f"{item['transform_id']}: transformations must connect different units")
        if item["relation"] == "one_to_one" and source["value"] != target["value"]:
            raise ContractError(f"{item['transform_id']}: one_to_one endpoints must have equal values")
        if item["relation"] == "aggregate" and source["value"] < target["value"]:
            raise ContractError(f"{item['transform_id']}: aggregate source cannot be smaller than target")
        enriched.append({
            **item,
            "source_label": source["label"],
            "source_value": source["value"],
            "source_unit": source["unit"],
            "target_label": target["label"],
            "target_value": target["value"],
            "target_unit": target["unit"],
        })
    return enriched


def parse_views(path: Path | None, graph_views: list[str]) -> tuple[list[dict[str, Any]], str]:
    if path is None:
        views = [
            {"view_id": view_id, "display_name": view_id, "order": index, "is_default": index == 0}
            for index, view_id in enumerate(graph_views)
        ]
        return views, graph_views[0]
    raw = read_tsv(path, VIEW_FIELDS)
    seen: set[str] = set()
    defaults: list[str] = []
    views: list[dict[str, Any]] = []
    for line, row in enumerate(raw, start=2):
        view_id = row["view_id"].strip()
        display_name = row["display_name"].strip()
        name = f"{path}:{line}"
        if not view_id or not display_name:
            raise ContractError(f"{name}: view_id and display_name must be non-empty")
        if view_id in seen:
            raise ContractError(f"{name}: duplicate view_id {view_id}")
        seen.add(view_id)
        is_default = row["is_default"].strip().lower() in TRUE_VALUES
        if is_default:
            defaults.append(view_id)
        views.append({
            "view_id": view_id,
            "display_name": display_name,
            "order": strict_integer(row["order"], "order", name),
            "is_default": is_default,
        })
    missing = sorted(set(graph_views) - seen)
    extra = sorted(seen - set(graph_views))
    if missing or extra:
        raise ContractError(f"Views file mismatch: missing={missing}, extra={extra}")
    if len(defaults) > 1:
        raise ContractError(f"Views file declares multiple defaults: {defaults}")
    views.sort(key=lambda row: (row["order"], row["view_id"]))
    default_view = defaults[0] if defaults else views[0]["view_id"]
    return views, default_view


def validate_graph(
    nodes: list[dict[str, Any]], links: list[dict[str, Any]], graph_views: list[str],
) -> dict[str, dict[str, Any]]:
    nodes_by_view: dict[str, dict[str, dict[str, Any]]] = defaultdict(dict)
    links_by_view: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for node in nodes:
        nodes_by_view[node["view_id"]][node["node_id"]] = node
    for link in links:
        links_by_view[link["view_id"]].append(link)
    unknown_link_views = sorted(set(links_by_view) - set(graph_views))
    if unknown_link_views:
        raise ContractError(f"Links reference unknown views: {unknown_link_views}")

    summaries: dict[str, dict[str, Any]] = {}
    for view_id in graph_views:
        view_nodes = nodes_by_view[view_id]
        view_links = links_by_view.get(view_id, [])
        if not view_links:
            raise ContractError(f"{view_id}: at least one link is required")
        units = sorted({str(node["unit"]) for node in view_nodes.values()})
        if len(units) != 1:
            raise ContractError(f"{view_id}: mixed units are not allowed: {units}")
        incoming: dict[str, int] = defaultdict(int)
        outgoing: dict[str, int] = defaultdict(int)
        adjacency: dict[str, list[str]] = defaultdict(list)
        indegree: dict[str, int] = {node_id: 0 for node_id in view_nodes}
        for link in view_links:
            source = link["source"]
            target = link["target"]
            if source not in view_nodes or target not in view_nodes:
                raise ContractError(f"{view_id}/{link['link_id']}: unknown source or target")
            source_stage = int(view_nodes[source]["stage"])
            target_stage = int(view_nodes[target]["stage"])
            if target_stage != source_stage + 1:
                raise ContractError(
                    f"{view_id}/{link['link_id']}: link must advance exactly one stage "
                    f"({source_stage} -> {target_stage})"
                )
            incoming[target] += int(link["value"])
            outgoing[source] += int(link["value"])
            adjacency[source].append(target)
            indegree[target] += 1
        roots = [node_id for node_id, degree in indegree.items() if degree == 0]
        if len(roots) != 1:
            raise ContractError(f"{view_id}: expected exactly one root, found {roots}")
        root = roots[0]
        for node_id, node in view_nodes.items():
            value = int(node["value"])
            if node_id != root and incoming[node_id] != value:
                raise ContractError(
                    f"{view_id}/{node_id}: incoming {incoming[node_id]} does not equal node value {value}"
                )
            if node_id in outgoing and outgoing[node_id] != value:
                raise ContractError(
                    f"{view_id}/{node_id}: outgoing {outgoing[node_id]} does not equal node value {value}"
                )
        queue = deque([node_id for node_id, degree in indegree.items() if degree == 0])
        remaining = dict(indegree)
        visited: list[str] = []
        while queue:
            node_id = queue.popleft()
            visited.append(node_id)
            for target in adjacency.get(node_id, []):
                remaining[target] -= 1
                if remaining[target] == 0:
                    queue.append(target)
        if len(visited) != len(view_nodes):
            raise ContractError(f"{view_id}: graph contains a cycle")
        reachable: set[str] = set()
        stack = [root]
        while stack:
            node_id = stack.pop()
            if node_id in reachable:
                continue
            reachable.add(node_id)
            stack.extend(adjacency.get(node_id, []))
        if reachable != set(view_nodes):
            raise ContractError(f"{view_id}: nodes are not all reachable from root {root}")
        summaries[view_id] = {
            "root_node_id": root,
            "root_value": int(view_nodes[root]["value"]),
            "unit": units[0],
            "nodes": len(view_nodes),
            "links": len(view_links),
        }
    return summaries


def enrich(
    nodes: list[dict[str, Any]], links: list[dict[str, Any]], views: list[dict[str, Any]],
    summaries: dict[str, dict[str, Any]], transformations: list[dict[str, Any]], node_layout: str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    nodes_by_view: dict[str, list[dict[str, Any]]] = defaultdict(list)
    links_by_view: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for node in nodes:
        nodes_by_view[node["view_id"]].append(dict(node))
    for link in links:
        links_by_view[link["view_id"]].append(dict(link))

    all_nodes: list[dict[str, Any]] = []
    all_links: list[dict[str, Any]] = []
    payload_views: dict[str, Any] = {}
    kind_order = {"retained": 0, "branch": 1, "loss": 2}
    for view in views:
        view_id = view["view_id"]
        current_nodes = nodes_by_view[view_id]
        stages = sorted({int(node["stage"]) for node in current_nodes})
        stage_position = {stage: index for index, stage in enumerate(stages)}
        stage_x = {
            stage: (0.04 if len(stages) == 1 else 0.04 + 0.90 * index / (len(stages) - 1))
            for index, stage in enumerate(stages)
        }
        stage_nodes: dict[int, list[dict[str, Any]]] = defaultdict(list)
        for node in current_nodes:
            stage_nodes[int(node["stage"])].append(node)
        packed_scale = 0.0
        if node_layout == "main_top_losses":
            candidates: list[float] = []
            for stage in stages:
                positive = [node for node in stage_nodes[stage] if int(node["value"]) > 0]
                total = sum(int(node["value"]) for node in positive)
                usable = 1.0 - PACKED_PAD_FRACTION * max(0, len(positive) - 1)
                if total and usable <= 0:
                    raise ContractError(f"{view_id}: too many visible nodes for packed layout")
                if total:
                    candidates.append(usable / total)
            packed_scale = min(candidates) if candidates else 0.0
        ordered_nodes: list[dict[str, Any]] = []
        for stage in stages:
            items = sorted(
                stage_nodes[stage],
                key=lambda row: (
                    kind_order[row["node_kind"]],
                    -int(row["value"]) if row["node_kind"] == "loss" else 0,
                    row["label"],
                    row["node_id"],
                ),
            )
            cursor = 0.0
            for position, node in enumerate(items):
                stage_index = stage_position[stage]
                x = stage_x[stage]
                if node_layout == "main_top_losses" and node["node_kind"] == "loss" and stage_index > 0:
                    previous_x = stage_x[stages[stage_index - 1]]
                    x = previous_x + LOSS_GUTTER_FRACTION * (stage_x[stage] - previous_x)
                if node_layout == "main_top_losses":
                    node_fraction = int(node["value"]) * packed_scale
                    y = cursor + node_fraction / 2.0 if int(node["value"]) > 0 else min(cursor, 1.0)
                    if int(node["value"]) > 0:
                        cursor += node_fraction + PACKED_PAD_FRACTION
                else:
                    count = len(items)
                    y = 0.05 if count == 1 else 0.05 + 0.85 * position / (count - 1)
                node["pct_root"] = (
                    100.0 * int(node["value"]) / summaries[view_id]["root_value"]
                    if summaries[view_id]["root_value"] else 0.0
                )
                node["color"] = color_for(node["node_kind"], stage_index)
                node["x"] = x
                node["y"] = y
                node["layout_order"] = position
                ordered_nodes.append(node)
        index_by_id = {node["node_id"]: index for index, node in enumerate(ordered_nodes)}
        lookup = {node["node_id"]: node for node in ordered_nodes}
        enriched_links: list[dict[str, Any]] = []
        for link in links_by_view[view_id]:
            source = lookup[link["source"]]
            target = lookup[link["target"]]
            parent = int(source["value"])
            root = int(summaries[view_id]["root_value"])
            enriched_links.append({
                **link,
                "source_label": source["label"],
                "target_label": target["label"],
                "source_stage": source["stage"],
                "target_stage": target["stage"],
                "target_kind": target["node_kind"],
                "unit": source["unit"],
                "parent_count": parent,
                "pct_parent": 100.0 * int(link["value"]) / parent if parent else 0.0,
                "pct_root": 100.0 * int(link["value"]) / root if root else 0.0,
                "source_index": index_by_id[link["source"]],
                "target_index": index_by_id[link["target"]],
            })
        all_nodes.extend(ordered_nodes)
        all_links.extend(enriched_links)
        payload_views[view_id] = {
            "display_name": view["display_name"],
            "unit": summaries[view_id]["unit"],
            "root_value": summaries[view_id]["root_value"],
            "node_layout": node_layout,
            "nodes": ordered_nodes,
            "links": enriched_links,
        }
    payload = {
        "node_layout": node_layout,
        "view_order": [view["view_id"] for view in views],
        "views": payload_views,
        "transformations": transformations,
    }
    return all_nodes, all_links, payload


def safe_stem(view_id: str, used: set[str]) -> str:
    base = re.sub(r"[^A-Za-z0-9._-]+", "_", view_id).strip("._-") or "view"
    stem = base
    if stem in used:
        suffix = hashlib.sha1(view_id.encode("utf-8")).hexdigest()[:8]
        stem = f"{base}-{suffix}"
    used.add(stem)
    return stem


def compact_count(value: int) -> str:
    number = int(value)
    for divisor, suffix in ((1_000_000_000, "B"), (1_000_000, "M"), (1_000, "K")):
        if abs(number) >= divisor:
            text = f"{number / divisor:.2f}".rstrip("0").rstrip(".")
            return f"{text}{suffix}"
    return f"{number:,}"


def packed_svg_geometry(
    view: dict[str, Any], width: int, flow_top: float, flow_height: float, node_pad: float = 12.0,
) -> tuple[dict[str, dict[str, float]], float]:
    nodes = [node for node in view["nodes"] if int(node["value"]) > 0]
    if view.get("node_layout") == "main_top_losses":
        stages = sorted({int(node["stage"]) for node in view["nodes"]})
        stage_nodes: dict[int, list[dict[str, Any]]] = defaultdict(list)
        for node in nodes:
            stage_nodes[int(node["stage"])].append(node)
        candidates = [
            (1.0 - PACKED_PAD_FRACTION * max(0, len(stage_nodes[stage]) - 1))
            / sum(int(node["value"]) for node in stage_nodes[stage])
            for stage in stages if stage_nodes[stage]
        ]
        normalized_scale = min(candidates) if candidates else 0.0
        scale = normalized_scale * flow_height
        geometry = {
            node["node_id"]: {
                "x": 55 + float(node["x"]) * (width - 235),
                "y0": flow_top + (float(node["y"]) - int(node["value"]) * normalized_scale / 2.0) * flow_height,
                "y1": flow_top + (float(node["y"]) + int(node["value"]) * normalized_scale / 2.0) * flow_height,
                "height": int(node["value"]) * scale,
            }
            for node in nodes
        }
        return geometry, scale
    stages = sorted({int(node["stage"]) for node in view["nodes"]})
    stage_nodes: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for node in nodes:
        stage_nodes[int(node["stage"])].append(node)
    candidates: list[float] = []
    for stage in stages:
        items = stage_nodes.get(stage, [])
        total = sum(int(node["value"]) for node in items)
        if total:
            usable = flow_height - node_pad * max(0, len(items) - 1)
            if usable <= 0:
                raise ContractError(f"{view['display_name']}: SVG height is too small for node padding")
            candidates.append(usable / total)
    scale = min(candidates) if candidates else 0.0
    geometry: dict[str, dict[str, float]] = {}
    for stage in stages:
        items = sorted(
            stage_nodes.get(stage, []),
            key=lambda row: (row["node_kind"] != "retained", -int(row["value"]), row["label"]),
        )
        cursor = flow_top
        for node in items:
            height = int(node["value"]) * scale
            geometry[node["node_id"]] = {
                "x": 55 + float(node["x"]) * (width - 235),
                "y0": cursor,
                "y1": cursor + height,
                "height": height,
            }
            cursor += height + node_pad
    return geometry, scale


def render_svg(view: dict[str, Any], title: str) -> str:
    nodes = view["nodes"]
    links = view["links"]
    stages = sorted({int(node["stage"]) for node in nodes})
    visible_nodes = [node for node in nodes if int(node["value"]) > 0]
    zero_nodes = [node for node in nodes if int(node["value"]) == 0]
    width = max(1300, 300 * len(stages))
    flow_top = 82.0
    flow_height = 460.0
    zero_height = 42 if zero_nodes else 0
    height = int(flow_top + flow_height + 44 + zero_height)
    geometry, scale = packed_svg_geometry(view, width, flow_top, flow_height)
    node_width = 20.0
    parts = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}' "
        f"data-node-layout='{html.escape(str(view.get('node_layout', 'snap')), quote=True)}' role='img' aria-label='{html.escape(title, quote=True)}'>",
        "<style>text{font-family:Arial,sans-serif}.title{font-size:20px;font-weight:700;fill:#0f172a}.stage{font-size:12px;font-weight:700;fill:#475569}.label{font-size:11px;fill:#1f2937;paint-order:stroke fill;stroke:#fff;stroke-width:3px;stroke-linejoin:round}.link{fill-opacity:.45}</style>",
        f"<rect width='{width}' height='{height}' fill='#ffffff'/>",
        f"<text class='title' x='24' y='32'>{html.escape(title)}</text>",
    ]
    for stage_index, stage in enumerate(stages):
        fraction = 0.04 if len(stages) == 1 else 0.04 + 0.90 * stage_index / (len(stages) - 1)
        x = 55 + fraction * (width - 235)
        parts.append(f"<text class='stage' x='{x:.1f}' y='58'>Stage {stage}</text>")
    outgoing_cursor = {node_id: item["y0"] for node_id, item in geometry.items()}
    incoming_cursor = {node_id: item["y0"] for node_id, item in geometry.items()}
    ordered_links = sorted(
        (link for link in links if int(link["value"]) > 0),
        key=lambda row: (
            int(row["source_stage"]),
            geometry[row["source"]]["y0"],
            geometry[row["target"]]["y0"],
        ),
    )
    for link in ordered_links:
        if int(link["value"]) == 0:
            continue
        source_geo = geometry[link["source"]]
        target_geo = geometry[link["target"]]
        sx = source_geo["x"] + node_width
        tx = target_geo["x"]
        band_height = int(link["value"]) * scale
        sy0 = outgoing_cursor[link["source"]]
        sy1 = sy0 + band_height
        ty0 = incoming_cursor[link["target"]]
        ty1 = ty0 + band_height
        outgoing_cursor[link["source"]] = sy1
        incoming_cursor[link["target"]] = ty1
        color = "#94A3B8" if link["target_kind"] == "loss" else ("#00A087" if link["target_kind"] == "branch" else "#4DBBD5")
        title_text = html.escape(
            f"{link['source_label']} -> {link['target_label']}: {link['value']:,} {link['unit']}; "
            f"{link['pct_parent']:.2f}% parent; {link['pct_root']:.2f}% root; {link['reason']}"
        )
        control = max(55.0, (tx - sx) * 0.42)
        parts.append(
            f"<path class='link' data-link-id='{html.escape(str(link['link_id']), quote=True)}' "
            f"data-value='{int(link['value'])}' data-width='{band_height:.12f}' fill='{color}' "
            f"d='M {sx:.12f} {sy0:.12f} C {sx+control:.12f} {sy0:.12f}, {tx-control:.12f} {ty0:.12f}, {tx:.12f} {ty0:.12f} "
            f"L {tx:.12f} {ty1:.12f} C {tx-control:.12f} {ty1:.12f}, {sx+control:.12f} {sy1:.12f}, {sx:.12f} {sy1:.12f} Z'>"
            f"<title>{title_text}</title></path>"
        )
    for node in visible_nodes:
        item = geometry[node["node_id"]]
        x = item["x"]
        y0 = item["y0"]
        node_height = item["height"]
        label_y = y0 + max(9.0, node_height / 2 + 4.0)
        full = html.escape(f"{node['label']}: {node['value']:,} {node['unit']} ({node['pct_root']:.2f}% root)")
        shown = f"{node['label']} · {compact_count(int(node['value']))}"
        parts.append(
            f"<rect data-node-id='{html.escape(str(node['node_id']), quote=True)}' data-value='{int(node['value'])}' "
            f"data-node-kind='{node['node_kind']}' data-layout-order='{int(node['layout_order'])}' "
            f"data-layout-x='{float(node['x']):.6f}' data-layout-y='{float(node['y']):.6f}' "
            f"data-height='{node_height:.12f}' x='{x:.12f}' y='{y0:.12f}' width='{node_width:.12f}' height='{node_height:.12f}' "
            f"rx='2' fill='{node['color']}' stroke='#ffffff'><title>{full}</title></rect>"
            f"<text class='label' x='{x+node_width+6:.2f}' y='{label_y:.2f}'>{html.escape(shown)}<title>{full}</title></text>"
        )
    if zero_nodes:
        labels = " · ".join(str(node["label"]) for node in zero_nodes)
        parts.append(
            f"<text class='stage' x='24' y='{height-18}'>0-count fates: {html.escape(labels)}"
            f"<title>{html.escape(labels)}</title></text>"
        )
    parts.append("</svg>")
    return "".join(parts)


BASE_REPORT_CSS = r"""
:root{--ink:#172033;--muted:#64748b;--line:#dbe3ea;--panel:#fff;--wash:#f6f8fa;--blue:#3C5488;--cyan:#4DBBD5}
*{box-sizing:border-box}body{margin:0;background:var(--wash);color:var(--ink);font:14px/1.55 Arial,sans-serif}main{max-width:1480px;margin:auto;padding:24px}.bio-flow-hero{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:22px;margin-bottom:18px}.bio-flow-hero h1{margin:0 0 8px;font-size:25px}.bio-flow-muted{color:var(--muted)}footer{color:var(--muted);text-align:center;padding:10px}
@media(max-width:800px){main{padding:10px}.bio-flow-hero{padding:15px}}
@media print{body{background:#fff}.bio-flow-hero{border:0;padding:0}}
"""


MODULE_CSS = r"""
#bio-read-flow{--bf-ink:#172033;--bf-muted:#64748b;--bf-line:#dbe3ea;--bf-soft:#f5f9fc;--bf-cyan:#4DBBD5;--bf-green:#00A087;color:var(--bf-ink);background:#fff;border:1px solid var(--bf-line);border-radius:10px;padding:18px;min-width:0}
#bio-read-flow *{box-sizing:border-box}.bio-flow-module-head{display:flex;gap:16px;justify-content:space-between;align-items:flex-start;flex-wrap:wrap}.bio-flow-module-head h3{font-size:18px;margin:0 0 5px}.bio-flow-module-head p{margin:0;color:var(--bf-muted)}
.bio-flow-notice{padding:10px 12px;border-left:4px solid var(--bf-cyan);background:#f2fbfd;margin:12px 0}.bio-flow-controls{display:flex;gap:10px;align-items:center;flex-wrap:wrap;margin:12px 0}.bio-flow-view-label{display:flex;gap:8px;align-items:center;font-weight:700;color:#334155}.bio-flow-view-select,.bio-flow-expand,.bio-flow-pager button{font:inherit;border:1px solid #b8c6d3;border-radius:7px;background:#fff;padding:7px 11px;color:#334155}.bio-flow-view-select{min-width:210px}.bio-flow-expand{margin-left:auto}
.bio-flow-figure{border:1px solid var(--bf-line);border-radius:9px;padding:14px;min-width:0}.bio-flow-unit-block{border:1px solid #e3eaf0;border-radius:9px;background:#fff;padding:12px;margin-bottom:10px}.bio-flow-unit-head{display:flex;justify-content:space-between;gap:12px;align-items:baseline;flex-wrap:wrap}.bio-flow-unit-head h4{margin:0;font-size:15px}.bio-flow-unit-head span{color:var(--bf-muted);font-size:12px}.bio-flow-plot-shell{width:100%;overflow-x:auto}.bio-flow-plot{width:100%;min-width:0}
.bio-flow-kpis{display:grid;grid-template-columns:repeat(3,minmax(130px,1fr));gap:8px;margin:10px 0}.bio-flow-kpi{border:1px solid #e1e9ef;background:var(--bf-soft);border-radius:7px;padding:9px 10px}.bio-flow-kpi b{display:block;font-size:15px}.bio-flow-kpi span{display:block;color:var(--bf-muted);font-size:11px;margin-top:2px}
.bio-flow-bridge{display:grid;grid-template-columns:minmax(0,1fr) 64px minmax(0,1fr);gap:10px;align-items:center;border:1px dashed #7c93a8;border-radius:10px;background:#f8fbfd;padding:12px 16px;margin:12px 4%}.bio-flow-bridge-node{text-align:center}.bio-flow-bridge-node b{display:block;font-size:16px}.bio-flow-bridge-node span{display:block;color:var(--bf-muted);font-size:12px}.bio-flow-bridge-arrow{text-align:center;color:#2f6f8f;font-size:26px;font-weight:700}.bio-flow-bridge-copy{grid-column:1/-1;text-align:center;color:#4d6376;font-size:12px}.bio-flow-zero{display:flex;gap:6px;align-items:center;flex-wrap:wrap;margin:5px 0 0;color:var(--bf-muted);font-size:11px}.bio-flow-zero b{color:#475569}.bio-flow-zero span{border:1px solid #e1e7ec;border-radius:999px;padding:2px 7px;background:#fafcfd}
.bio-flow-details{margin-top:12px;border:1px solid var(--bf-line);border-radius:8px;background:#fbfdff;padding:10px 12px}.bio-flow-details summary{font-weight:700;cursor:pointer}.bio-flow-table-wrap{overflow:auto;margin-top:10px}.bio-flow-table-wrap table{border-collapse:collapse;width:100%}.bio-flow-table-wrap th,.bio-flow-table-wrap td{border-bottom:1px solid var(--bf-line);padding:8px 9px;text-align:left;vertical-align:top}.bio-flow-table-wrap th{background:#eef3f6;position:sticky;top:0}.bio-flow-transform-row{background:#f2f8fb}.bio-flow-pager{display:flex;gap:9px;align-items:center;justify-content:flex-end;margin-top:10px}.bio-flow-expanded{position:fixed!important;inset:18px!important;z-index:1000!important;overflow:auto!important;background:#fff!important;box-shadow:0 0 0 9999px #111827cc!important}.bio-flow-expanded .bio-flow-plot{min-width:900px}.bio-flow-view-select:focus-visible,.bio-flow-expand:focus-visible,.bio-flow-pager button:focus-visible,.bio-flow-details summary:focus-visible{outline:3px solid #4DBBD566;outline-offset:2px}#bio-read-flow .sankey-node text{paint-order:stroke fill;stroke:#fff;stroke-width:3px;stroke-linejoin:round}
@media(max-width:800px){#bio-read-flow{padding:12px}.bio-flow-kpis{grid-template-columns:1fr}.bio-flow-bridge{grid-template-columns:1fr}.bio-flow-bridge-arrow{transform:rotate(90deg)}.bio-flow-bridge-copy{grid-column:1}.bio-flow-plot{min-width:720px}.bio-flow-expand{margin-left:0}}
@media print{#bio-read-flow{border:0;padding:0}.bio-flow-controls,.bio-flow-pager,.bio-flow-expand{display:none!important}.bio-flow-figure,.bio-flow-unit-block{border:0;padding:0}.bio-flow-plot{min-width:0!important}.bio-flow-details{display:none}.bio-flow-bridge{break-inside:avoid}.bio-flow-unit-block{break-inside:avoid}}
"""


def report_section(views: list[dict[str, Any]], transformations: list[dict[str, Any]]) -> str:
    if len(views) == 1 and not transformations:
        heading = "数据处理"
        subtitle = f"展示 {views[0]['display_name']}；所有带宽均为同一单位。"
        notice = "实线带仅表示本样本 reads 的互斥守恒分支；本视图不包含跨单位转换。"
        controls = (
            "<div class='bio-flow-controls'>"
            "<button type='button' class='bio-flow-expand' id='bio-flow-expand'>Expand</button>"
            "</div>"
        )
    else:
        heading = "数据守恒与单位转换"
        subtitle = "通过下拉框切换各守恒视图；不同单位分开展示。"
        notice = "实线带仅用于同单位守恒分支；跨单位转换单独标注，转换两侧数量不可按带宽直接比较。"
        controls = (
            "<div class='bio-flow-controls'>"
            "<label class='bio-flow-view-label' for='bio-flow-view'>流程视图 <select class='bio-flow-view-select' id='bio-flow-view' aria-label='流程视图'></select></label>"
            "<button type='button' class='bio-flow-expand' id='bio-flow-expand'>Expand</button>"
            "</div>"
        )
    return f"""
<section class='bio-flow-section' id='bio-read-flow'>
  <div class='bio-flow-module-head'><div><h3>{html.escape(heading)}</h3><p>{html.escape(subtitle)}</p></div></div>
  <p class='bio-flow-notice' id='bio-flow-conversion'>{html.escape(notice)}</p>
  {controls}
  <article class='bio-flow-figure' id='bio-flow-card'>
    <noscript>交互式流程图需要 JavaScript；静态结果请查看生成的 SVG。</noscript>
    <div id='bio-flow-detail'><div class='bio-flow-kpis' id='bio-flow-detail-kpis'></div><div class='bio-flow-plot-shell'><div class='bio-flow-plot' id='bio-flow-plot' role='img' aria-label='守恒流程 Sankey 图'></div></div><div class='bio-flow-zero' id='bio-flow-detail-zero'></div></div>
  </article>
  <details class='bio-flow-details' id='bio-flow-details'><summary>查看精确分支与转换明细</summary><div class='bio-flow-table-wrap'><table id='bio-flow-table'>
    <thead><tr><th>单位</th><th>From</th><th>To / fate</th><th>Count / conversion</th><th>Parent %</th><th>Root %</th><th>Reason</th></tr></thead>
    <tbody id='bio-flow-body'></tbody>
  </table></div><div class='bio-flow-pager'><button type='button' id='bio-flow-prev'>Previous</button><span id='bio-flow-status'></span><button type='button' id='bio-flow-next'>Next</button></div></details>
</section>
""".strip()


REPORT_JS = r"""
(function(){
'use strict';
const data=__PAYLOAD__,defaultView=__DEFAULT_VIEW__,defaultMode=__DEFAULT_MODE__,pageSize=10;
let page=1,active=defaultView;
const fmt=value=>Number(value).toLocaleString();
function compact(value){const number=Number(value);for(const [divisor,suffix] of [[1e9,'B'],[1e6,'M'],[1e3,'K']]){if(Math.abs(number)>=divisor)return (number/divisor).toFixed(2).replace(/\.00$/,'').replace(/(\.\d)0$/,'$1')+suffix}return fmt(number)}
const flowShortLabels={raw:'Raw',insert:'Insert',structure_loss:'结构未识别',primary_mapped:'Mapped',unmapped:'Unmapped',primary_q20:'MAPQ≥20',low_mapq:'MAPQ<20',gene_assigned:'Assigned',intergenic:'Intergenic',shared_antisense_ambiguous:'Shared antisense',shared_gene_ambiguous:'Shared gene',antisense_only:'Antisense only'};
function flowLabel(node){return (flowShortLabels[node.node_id]||node.label.replace(/ reads$/i,''))+' · '+compact(node.value)}
function visibleGraph(view){
  const nodes=view.nodes.filter(node=>Number(node.value)>0),ids=new Set(nodes.map(node=>node.node_id));
  const index={};nodes.forEach((node,i)=>index[node.node_id]=i);
  const links=view.links.filter(link=>Number(link.value)>0&&ids.has(link.source)&&ids.has(link.target)).map(link=>({...link,source_index:index[link.source],target_index:index[link.target]}));
  return {nodes,links};
}
function plotHeight(view){const counts={};view.nodes.filter(n=>Number(n.value)>0).forEach(n=>counts[n.stage]=(counts[n.stage]||0)+1);return Math.max(330,210+62*Math.max(1,...Object.values(counts)))}
function raiseSankeyLabels(el){el.querySelectorAll('.sankey').forEach(group=>{const nodes=group.querySelector('.sankey-node-set');if(nodes&&nodes.parentNode)nodes.parentNode.appendChild(nodes)})}
function plotView(el,view){
  const graph=visibleGraph(view),nodes=graph.nodes,links=graph.links,height=plotHeight(view);el.style.height=height+'px';
  if(!window.Plotly){el.textContent='Plotly.js 不可用；请查看静态 SVG。';return}
  const arrangement=data.node_layout==='main_top_losses'?'fixed':'snap';
  const nodeSpec={label:nodes.map(flowLabel),color:nodes.map(n=>n.color),x:nodes.map(n=>n.x),pad:12,thickness:19,line:{color:'#fff',width:1.2},customdata:nodes.map(n=>[n.label,n.value,n.unit,n.pct_root,n.scope,n.node_kind]),hovertemplate:'%{customdata[0]}<br>%{customdata[1]:,} %{customdata[2]}<br>%{customdata[3]:.2f}% of root<br>scope: %{customdata[4]}<br>%{customdata[5]}<extra></extra>'};
  if(arrangement==='fixed')nodeSpec.y=nodes.map(n=>n.y);
  Plotly.react(el,[{type:'sankey',orientation:'h',arrangement,node:nodeSpec,link:{source:links.map(l=>l.source_index),target:links.map(l=>l.target_index),value:links.map(l=>l.value),color:links.map(l=>l.target_kind==='loss'?'rgba(148,163,184,.38)':'rgba(77,187,213,.55)'),customdata:links.map(l=>[l.unit,l.reason,l.pct_parent,l.pct_root]),hovertemplate:'%{source.label} → %{target.label}<br>%{value:,} %{customdata[0]}<br>%{customdata[2]:.2f}% of parent<br>%{customdata[3]:.2f}% of root<br>%{customdata[1]}<extra></extra>'}}],{height,margin:{t:25,r:130,b:24,l:100},font:{family:'Arial,"Noto Sans SC",sans-serif',size:12,color:'#334155'},paper_bgcolor:'#fff',plot_bgcolor:'#fff'},{responsive:true,displaylogo:false,modeBarButtonsToRemove:['sendDataToCloud','lasso2d','select2d']}).then(()=>raiseSankeyLabels(el));
}
function zeroFates(view){return view.nodes.filter(node=>Number(node.value)===0).map(node=>node.label)}
function zeroBox(view){const wrap=document.createElement('div');wrap.className='bio-flow-zero';const zeros=zeroFates(view);if(!zeros.length)return wrap;wrap.innerHTML='<b>0-count fates</b>';zeros.forEach(label=>{const chip=document.createElement('span');chip.textContent=label;wrap.appendChild(chip)});return wrap}
function kpis(view){const retained=[...view.nodes].filter(n=>n.node_kind==='retained').sort((a,b)=>b.stage-a.stage),losses=[...view.nodes].filter(n=>n.node_kind==='loss').sort((a,b)=>b.value-a.value),terminal=retained[0],largest=losses[0];return [[fmt(view.root_value)+' '+view.unit,'Root denominator'],[terminal?fmt(terminal.value)+' '+terminal.unit:'—',terminal?terminal.label:'Final retained'],[largest?fmt(largest.value)+' '+largest.unit:'0',largest?'Largest loss: '+largest.label:'No loss']]}
function kpiBox(view){const wrap=document.createElement('div');wrap.className='bio-flow-kpis';kpis(view).forEach(([value,label])=>{const card=document.createElement('div');card.className='bio-flow-kpi';card.innerHTML='<b></b><span></span>';card.querySelector('b').textContent=value;card.querySelector('span').textContent=label;wrap.appendChild(card)});return wrap}
function rowsFor(mode){const rows=[];data.transformations.filter(item=>item.target_view_id===mode).forEach(item=>rows.push({kind:'transform',unit:item.source_unit+' → '+item.target_unit,from:item.source_label,to:item.target_label,count:fmt(item.source_value)+' → '+fmt(item.target_value),parent:'—',root:'—',reason:item.label+'; '+item.reason}));data.views[mode].links.forEach(link=>rows.push({kind:'link',unit:link.unit,from:link.source_label,to:link.target_label,count:fmt(link.value)+' '+link.unit,parent:Number(link.pct_parent).toFixed(2)+'%',root:Number(link.pct_root).toFixed(2)+'%',reason:link.reason}));return rows}
function table(mode){
  const body=document.getElementById('bio-flow-body');body.replaceChildren();
  const rows=rowsFor(mode),total=Math.max(1,Math.ceil(rows.length/pageSize));page=Math.min(Math.max(1,page),total);
  rows.forEach((item,index)=>{const row=document.createElement('tr');if(item.kind==='transform')row.className='bio-flow-transform-row';row.dataset.pageHidden=(index<(page-1)*pageSize||index>=page*pageSize)?'true':'false';if(row.dataset.pageHidden==='true')row.style.display='none';[item.unit,item.from,item.to,item.count,item.parent,item.root,item.reason].forEach(value=>{const cell=document.createElement('td');cell.textContent=value;row.appendChild(cell)});body.appendChild(row)});
  document.getElementById('bio-flow-status').textContent='Page '+page+' of '+total+' · '+rows.length+' rows';
  document.getElementById('bio-flow-prev').disabled=page<=1;document.getElementById('bio-flow-next').disabled=page>=total;
}
function renderDetail(viewId){const view=data.views[viewId],kpi=document.getElementById('bio-flow-detail-kpis');kpi.replaceChildren(...kpiBox(view).children);const zero=document.getElementById('bio-flow-detail-zero');zero.replaceChildren(...zeroBox(view).childNodes);plotView(document.getElementById('bio-flow-plot'),view)}
function selectView(mode){active=mode;page=1;const selector=document.getElementById('bio-flow-view');if(selector)selector.value=mode;renderDetail(mode);table(mode)}
function init(){const selector=document.getElementById('bio-flow-view');if(selector){data.view_order.forEach(key=>{const option=document.createElement('option');option.value=key;option.textContent=data.views[key].display_name;selector.appendChild(option)});selector.onchange=()=>selectView(selector.value)}document.getElementById('bio-flow-prev').onclick=()=>{page--;table(active)};document.getElementById('bio-flow-next').onclick=()=>{page++;table(active)};document.getElementById('bio-flow-expand').onclick=()=>{const card=document.getElementById('bio-flow-card');card.classList.toggle('bio-flow-expanded');document.getElementById('bio-flow-expand').textContent=card.classList.contains('bio-flow-expanded')?'Close':'Expand';if(window.Plotly)document.querySelectorAll('#bio-flow-detail .bio-flow-plot.js-plotly-plot').forEach(el=>{Plotly.Plots.resize(el);raiseSankeyLabels(el)})};document.addEventListener('keydown',event=>{if(event.key==='Escape'){document.getElementById('bio-flow-card').classList.remove('bio-flow-expanded');document.getElementById('bio-flow-expand').textContent='Expand'}});selectView(active)}
if(document.readyState==='loading')document.addEventListener('DOMContentLoaded',init);else init();
})();
"""


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_outputs(
    output_dir: Path, title: str, nodes: list[dict[str, Any]], links: list[dict[str, Any]],
    payload: dict[str, Any], views: list[dict[str, Any]], default_view: str,
    default_mode: str, transformations: list[dict[str, Any]], validation: dict[str, Any], plotly_path: Path,
) -> None:
    tables = output_dir / "tables"
    figures = output_dir / "figures" / "read_flow_sankey"
    reports = output_dir / "reports"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)
    node_fields = NODE_FIELDS + ["pct_root", "color", "x", "y", "layout_order"]
    link_fields = LINK_FIELDS + [
        "source_label", "target_label", "source_stage", "target_stage", "target_kind", "unit",
        "parent_count", "pct_parent", "pct_root", "source_index", "target_index",
    ]
    nodes_path = tables / "read_flow_nodes.tsv"
    links_path = tables / "read_flow_links.tsv"
    payload_path = tables / "read_flow_payload.json"
    transforms_path = tables / "read_flow_transformations.tsv"
    write_tsv(nodes_path, nodes, node_fields)
    write_tsv(links_path, links, link_fields)
    if transformations:
        write_tsv(transforms_path, transformations, TRANSFORM_FIELDS + [
            "source_label", "source_value", "source_unit", "target_label", "target_value", "target_unit",
        ])
    payload_with_meta = {"title": title, "default_view": default_view, "default_mode": default_mode, **payload}
    payload_path.write_text(json.dumps(payload_with_meta, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    used_stems: set[str] = set()
    svg_paths: list[tuple[str, Path]] = []
    for view in views:
        view_id = view["view_id"]
        stem = safe_stem(view_id, used_stems)
        path = figures / f"{stem}.svg"
        path.write_text(render_svg(payload["views"][view_id], view["display_name"]), encoding="utf-8")
        svg_paths.append((view_id, path))
    safe_payload = json.dumps(payload, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")
    script = (
        REPORT_JS.replace("__DEFAULT_VIEW__", json.dumps(default_view), 1)
        .replace("__DEFAULT_MODE__", json.dumps(default_mode), 1)
        .replace("__PAYLOAD__", safe_payload, 1)
    )
    section = report_section(views, transformations)
    fragment_path = reports / "read_flow_sankey_section.html"
    fragment_path.write_text(f"<style>{MODULE_CSS}</style>\n{section}\n<script>{script}</script>\n", encoding="utf-8")
    plotly = plotly_path.read_text(encoding="utf-8", errors="strict")
    report_path = reports / "read_flow_sankey.html"
    full = (
        "<!doctype html><html><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>"
        f"<title>{html.escape(title)}</title><style>{BASE_REPORT_CSS}{MODULE_CSS}</style><script>{plotly}</script></head><body><main>"
        f"<header class='bio-flow-hero'><h1>{html.escape(title)}</h1><p class='bio-flow-muted'>经过校验、分单位守恒并可离线查看的完整流程。</p></header>"
        f"{section}<footer>Conserved bioinformatics flow - offline artifact</footer></main><script>{script}</script></body></html>"
    )
    report_path.write_text(full, encoding="utf-8")
    validation_path = reports / "sankey_validation.json"
    validation_path.write_text(json.dumps(validation, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    def manifest_target(path: Path) -> str:
        return Path(os.path.relpath(path, reports)).as_posix()

    manifest_rows: list[dict[str, str]] = [
        {"artifact": "read_flow_nodes", "description": "canonical conserved flow nodes", "target": manifest_target(nodes_path), "status": "up_to_date"},
        {"artifact": "read_flow_links", "description": "canonical conserved flow links", "target": manifest_target(links_path), "status": "up_to_date"},
        {"artifact": "read_flow_payload", "description": "shared static and interactive payload", "target": manifest_target(payload_path), "status": "up_to_date"},
    ]
    if transformations:
        manifest_rows.append({"artifact": "read_flow_transformations", "description": "audited cross-unit transformations", "target": manifest_target(transforms_path), "status": "up_to_date"})
    manifest_rows.extend({
        "artifact": f"read_flow_svg_{view_id}", "description": "static conserved flow SVG",
        "target": manifest_target(path), "status": "up_to_date",
    } for view_id, path in svg_paths)
    manifest_rows.extend([
        {"artifact": "read_flow_section", "description": "embeddable report section without Plotly bundle", "target": manifest_target(fragment_path), "status": "up_to_date"},
        {"artifact": "read_flow_standalone", "description": "self-contained offline interactive Sankey", "target": manifest_target(report_path), "status": "up_to_date"},
        {"artifact": "read_flow_validation", "description": "machine-readable flow validation", "target": manifest_target(validation_path), "status": "up_to_date"},
    ])
    write_tsv(reports / "sankey_manifest.tsv", manifest_rows, ["artifact", "description", "target", "status"])


def main() -> int:
    cli = arguments()
    try:
        nodes, graph_views = parse_nodes(cli.nodes)
        links = parse_links(cli.links)
        transformations = parse_transformations(cli.transformations)
        views, default_view = parse_views(cli.views, graph_views)
        summaries = validate_graph(nodes, links, graph_views)
        enriched_transformations = validate_transformations(transformations, nodes)
        enriched_nodes, enriched_links, payload = enrich(
            nodes, links, views, summaries, enriched_transformations, cli.node_layout,
        )
        script_dir = Path(__file__).resolve().parent
        bundled_plotly = script_dir / "read_flow_assets" / "plotly-3.5.1.min.js"
        if not bundled_plotly.is_file():
            bundled_plotly = script_dir.parent / "assets" / "plotly-3.5.1.min.js"
        plotly_path = cli.plotly_js or bundled_plotly
        graph_checks = [
            "non_negative_integer_counts", "unique_identifiers", "single_unit_per_view",
            "single_root", "acyclic", "stage_adjacent_links", "incoming_conservation",
            "outgoing_conservation", "root_reachability",
        ]
        default_mode = "connected" if cli.default_mode == "connected" and enriched_transformations else "view"
        if default_mode == "connected":
            graph_checks.append("cross_unit_transformations")
        if cli.node_layout == "main_top_losses":
            graph_checks.append("main_top_local_losses_layout")
        validation = {
            "ok": True,
            "title": cli.title,
            "default_view": default_view,
            "default_mode": default_mode,
            "node_layout": cli.node_layout,
            "view_order": [view["view_id"] for view in views],
            "views": summaries,
            "node_count": len(nodes),
            "link_count": len(links),
            "transformation_count": len(enriched_transformations),
            "checks": graph_checks,
            "presentation_checked": not cli.validate_only,
            "inputs": {
                "nodes": {"path": str(cli.nodes.resolve()), "sha256": sha256(cli.nodes)},
                "links": {"path": str(cli.links.resolve()), "sha256": sha256(cli.links)},
                "views": ({"path": str(cli.views.resolve()), "sha256": sha256(cli.views)} if cli.views else None),
                "transformations": ({"path": str(cli.transformations.resolve()), "sha256": sha256(cli.transformations)} if cli.transformations else None),
                "plotly": None,
            },
        }
        if cli.validate_only:
            print(json.dumps(validation, ensure_ascii=False, indent=2))
            return 0
        if cli.output_dir is None:
            raise ContractError("--output-dir is required unless --validate-only is used")
        if not plotly_path.is_file() or plotly_path.stat().st_size == 0:
            raise ContractError(f"Missing local Plotly.js asset: {plotly_path}")
        validation["checks"].append("offline_plotly_asset")
        validation["inputs"]["plotly"] = {"path": str(plotly_path.resolve()), "sha256": sha256(plotly_path)}
        build_outputs(
            cli.output_dir, cli.title, enriched_nodes, enriched_links, payload, views,
            default_view, default_mode, enriched_transformations, validation, plotly_path,
        )
        print(json.dumps({"ok": True, "output_dir": str(cli.output_dir.resolve()), "views": len(views)}, ensure_ascii=False))
        return 0
    except (ContractError, OSError, UnicodeError) as exc:
        print(json.dumps({"ok": False, "error": str(exc)}, ensure_ascii=False), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
