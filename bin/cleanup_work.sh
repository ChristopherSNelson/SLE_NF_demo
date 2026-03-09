#!/usr/bin/env bash
# Clean up Nextflow work directory, keeping conda envs cached.
# Run after a successful pipeline completion to reclaim disk space.
#
# Usage: bash bin/cleanup_work.sh [work_dir]

set -euo pipefail

WORK_DIR="${1:-work}"

if [ ! -d "$WORK_DIR" ]; then
    echo "Work directory not found: $WORK_DIR"
    exit 1
fi

# Size before cleanup
echo "Work directory: $WORK_DIR"
du -sh "$WORK_DIR" 2>/dev/null || true

# Keep conda envs (they're reusable across runs)
CONDA_DIR="$WORK_DIR/conda"

# Remove task directories (the bulk of disk usage)
echo "Removing task directories..."
find "$WORK_DIR" -maxdepth 2 -mindepth 2 -type d ! -path "*/conda/*" -exec rm -rf {} + 2>/dev/null || true

# Remove stale lock files
find "$WORK_DIR" -name ".lock" -delete 2>/dev/null || true

# Size after cleanup
echo "After cleanup:"
du -sh "$WORK_DIR" 2>/dev/null || true

if [ -d "$CONDA_DIR" ]; then
    echo "Conda envs preserved:"
    du -sh "$CONDA_DIR" 2>/dev/null || true
fi

echo "Done. Published results in --outdir are untouched."
