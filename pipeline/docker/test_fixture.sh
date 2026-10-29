#!/usr/bin/env bash
set -euo pipefail
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"
export TMPDIR="$repo_root/tmp"
export PYTHONDONTWRITEBYTECODE=1
image='iso-slam:20260930-v1.0.0'
fixture_name="${ISO_SLAM_FIXTURE_NAME:-docker_fixture_$(date -u +%Y%m%dT%H%M%S)}"
[[ "$fixture_name" =~ ^docker_fixture_[A-Za-z0-9_.-]+$ ]] || { printf 'Invalid fixture name\n' >&2; exit 2; }
fixture_root="$repo_root/tmp/$fixture_name"
[[ ! -e "$fixture_root" ]] || { printf 'Fixture already exists: %s\n' "$fixture_root" >&2; exit 2; }
docker run --rm --network none --user "$(id -u):$(id -g)" \
    --env "ISO_SLAM_FIXTURE_NAME=$fixture_name" \
    --mount "type=bind,src=$repo_root/config,dst=/work/config,readonly" \
    --mount "type=bind,src=$repo_root/pipeline,dst=/work/pipeline,readonly" \
    --mount "type=bind,src=$repo_root/tmp,dst=/work/tmp" \
    --workdir /work --entrypoint /opt/project/projects/ISO-SLAM/.envs/iso_slam_v4/bin/python \
    "$image" pipeline/scripts/create_fixture.py
cd "$fixture_root"
ISO_SLAM_CORES=16 ISO_SLAM_MEMORY_MB=131072 bash pipeline/docker/run_snakemake.sh config/fixture.yaml --dry-run
ISO_SLAM_CORES=16 ISO_SLAM_MEMORY_MB=131072 bash pipeline/docker/run_snakemake.sh config/fixture.yaml
test -s outputs/synthetic_fixture/_project/RUN_COMPLETE.tsv
