# 04_slam_calling

Stage wrapper directory. The reusable implementation remains in the top-level scripts directory and is invoked by `workflow/run.sh`.

V3 runs `audit_lr_library_strand.py` before transcript-aware assignment, records a project-wide orientation selection, and writes exact `assignment_summary.tsv`, `substitution_summary.tsv`, `site_call_summary.tsv`, and `ntr_sufficient_stats.tsv` inputs. `grand_slam_compatibility.json` is a hard readiness record; missing software or incomplete MD/NH tags must not be silently ignored.
