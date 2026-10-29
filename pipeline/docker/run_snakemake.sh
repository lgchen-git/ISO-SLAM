#!/usr/bin/env bash
set -euo pipefail
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"
config="${1:-config/run.yaml}"
shift $(( $# > 0 ? 1 : 0 ))
[[ "$config" == config/*.yaml && -f "$config" ]] || { printf 'Configuration must exist under config/\n' >&2; exit 2; }
for name in config inputs log outputs pipeline reference tmp; do
    [[ -d "$name" ]] || { printf 'Missing required root: %s\n' "$name" >&2; exit 2; }
done
python3 pipeline/scripts/validate_run_config.py --config "$config"
python3 pipeline/scripts/record_run_launch.py --config "$config"
mkdir -p tmp/docker_client tmp/container_runtime tmp/snakemake log/snakemake
export DOCKER_CONFIG="$repo_root/tmp/docker_client"
cores="${ISO_SLAM_CORES:-32}"
memory_mb="${ISO_SLAM_MEMORY_MB:-196608}"
image='iso-slam:20260930-v1.0.0'
mounts=(--mount "type=bind,src=$repo_root/tmp/snakemake,dst=/work/.snakemake" --mount "type=bind,src=$repo_root/log/snakemake,dst=/work/.snakemake/log")
for name in config inputs pipeline; do mounts+=(--mount "type=bind,src=$repo_root/$name,dst=/work/$name,readonly"); done
for name in log outputs tmp; do mounts+=(--mount "type=bind,src=$repo_root/$name,dst=/work/$name"); done
mounts+=(--mount "type=bind,src=$repo_root/reference,dst=/work/reference,readonly")
while IFS= read -r -d '' link; do
    source="$(readlink -f "$link")"
    mounts+=(--mount "type=bind,src=$source,dst=$source,readonly")
done < <(find reference -maxdepth 1 -type l -print0)
exec docker run --rm --network none --user "$(id -u):$(id -g)" \
    --cpus "$cores" --memory "${memory_mb}m" --memory-swap "${memory_mb}m" \
    --read-only --workdir /work "${mounts[@]}" \
    --env HOME=/work/tmp/container_runtime --env TMPDIR=/work/tmp/container_runtime \
    --env TMP=/work/tmp/container_runtime --env TEMP=/work/tmp/container_runtime \
    --env NUMBA_CACHE_DIR=/work/tmp/container_runtime/numba \
    --env XDG_CACHE_HOME=/work/tmp/container_runtime/cache \
    --env MPLCONFIGDIR=/work/tmp/container_runtime/mplconfig \
    "$image" --snakefile pipeline/Snakefile --configfile "$config" \
    --cores "$cores" --resources "mem_mb=$memory_mb" sample_slots=4 \
    --directory /work --rerun-incomplete --printshellcmds \
    --default-resources tmpdir=/work/tmp/container_runtime "$@"
