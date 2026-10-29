# ISO-SLAM

Docker-delivered Snakemake workflow for bulk long-read RNA labeling analysis.
The repository contains reusable workflow code, environments, report builders,
tests, and layout contracts.

Production sample sheets, raw reads, references, runtime images, result
archives, and server handoff material are excluded from the public release.

~~~bash
make lint
make smoke
bash workflow/run.sh --help
~~~

Use a local configuration for real data and keep inputs and references outside
the repository.
