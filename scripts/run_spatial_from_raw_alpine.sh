#!/usr/bin/env bash
# Alpine execution adapter for numbered workflow 07. All data destinations are
# defined in config.R (cfg$paths$spatial); scientific settings remain separate.
# TIFF sources and copy instructions: docs/spatial/README.md.
set -euo pipefail
[[ $# -ge 2 && $# -le 4 ]] || { echo 'Usage: bash scripts/run_spatial_from_raw_alpine.sh UPDATED_GOAL_DIR PHENOTYPE_CSV [IMAGE_IDS.txt] [NEW_RUN_NAME]' >&2; exit 2; }
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo 'First obtain a compute allocation: acompile --ntasks=1' >&2; exit 2; }
command -v apptainer >/dev/null
repo_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
absolute_file() { local parent; parent=$(cd -- "$(dirname -- "$1")" && pwd -P); printf '%s/%s\n' "$parent" "$(basename -- "$1")"; }
goals=$(cd -- "$1" && pwd -P)
phenotypes=$(absolute_file "$2")
args=("$goals" "$phenotypes")
binds=(--bind /projects --bind "$repo_dir" --bind "$goals" --bind "$(dirname -- "$phenotypes")")
if [[ $# -ge 3 ]]; then
  ids=$(absolute_file "$3")
  args+=("$ids")
  binds+=(--bind "$(dirname -- "$ids")")
fi
[[ $# -lt 4 ]] || args+=("$4")
container=${SPATIAL_CONTAINER:-${CURC_CONTAINER_DIR_OOD:?}/rstudio-server-4.4.1.sif}
overlay=${SPATIAL_OVERLAY:-/projects/${USER:?}/.rstudioserver/rstudio-4.4.1/rstudio-server-4.4.1_overlay.img}
library_dir=${SPATIAL_R_LIBRARY:-/projects/${USER}/Rstudio_libs/4.4.1}
[[ -r "$container" && -r "$overlay" && -d "$library_dir" ]] || { echo 'Existing container, overlay or library unavailable.' >&2; exit 2; }
binds+=(--bind "$library_dir" --bind "$(dirname -- "$container")")
[[ -z "${SCRATCHDIR:-}" ]] || binds+=(--bind "$SCRATCHDIR")
# R's system Renviron hard-codes the projects tmp_data path. Redirect it inside
# the container only; host temporary data stays inside this repository checkout.
spatial_tmp=$(mktemp -d "$repo_dir/.r_tmp.XXXXXX")
export TMPDIR="$spatial_tmp"
export APPTAINERENV_TMPDIR=/tmp APPTAINERENV_TMP=/tmp APPTAINERENV_TEMP=/tmp
export APPTAINERENV_R_LIBS_USER="$library_dir"
export APPTAINERENV_SPATIAL_ALLOW_UNVALIDATED=0
# Marks the container child so workflow 07 does not relaunch itself.
export APPTAINERENV_SPATIAL_ALPINE_WORKFLOW=1
binds+=(--bind "$spatial_tmp:/tmp" --bind "$spatial_tmp:/projects/${USER}/.rstudioserver/rstudio-4.4.1/tmp_data")
cd -- "$repo_dir"
apptainer exec "${binds[@]}" --overlay "${overlay}:ro" "$container" \
  Rscript "$repo_dir/scripts/07_compute_spatial_data.R" "${args[@]}"
