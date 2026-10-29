# Delivery package

`build_iso_slam_delivery_package.py` creates the customer-facing copy-only package from a completed configured result root. It copies only rendered reports and compact derived TSV tables, refuses symbolic-link sources and targets, emits SHA-256 metadata, and promotes a validated staging package atomically.

Run the configured workflow through the unified entry point. The customer package
is built after the project report and is included in the final completion target:

```bash
bash workflow/run.sh config/example_sample_5group.yaml
```
