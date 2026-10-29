#!/usr/bin/env bash
set -euo pipefail
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"
image='iso-slam:20260930-v1.0.0'
build_root="$repo_root/tmp/docker_build_20260930"
mkdir -p "$build_root/runtime" outputs/docker log/docker tmp/docker_client
export DOCKER_CONFIG="$repo_root/tmp/docker_client"
export TMPDIR="$build_root"
if [[ ! -f "$build_root/runtime.complete" ]]; then
    python3 pipeline/docker/stage_runtime.py "$build_root/runtime"
    touch "$build_root/runtime.complete"
fi
cp pipeline/docker/Dockerfile "$build_root/Dockerfile"
docker build --network=none --pull=false --tag "$image" "$build_root" >log/docker/build.log 2>&1
docker save --output outputs/docker/iso-slam_20260930_v1.0.0.tar "$image"
printf 'Built and exported %s\n' "$image"
