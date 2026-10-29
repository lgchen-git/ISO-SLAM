# 06 isoform dynamics

Per-read T>C sufficient statistics are joined to filtered IsoQuant assignments. V4.1 records caller/IsoQuant mismatches as `Assignment_conflict`; those reads are excluded from isoform-level posterior summaries. Consistency-qualified reads are classified as `New` (P>=0.50), `Old` (P<=0.30), `Ambiguous` (middle with covered T sites), or `No_T_sites` (no covered T). Expected New/Old counts, deterministic technical bootstrap intervals, and within-sample isoform-usage shifts are produced without biological p-values. AS/APA is not generated in the active V4 pipeline.
