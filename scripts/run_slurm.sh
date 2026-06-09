#!/bin/bash
#SBATCH --job-name=merge_bedmethyl
#SBATCH --output=merge_bedmethyl_%j.out
#SBATCH --error=merge_bedmethyl_%j.err
#SBATCH --partition=ngen-ko
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
BIN="${PROJECT_ROOT}/build/merge_bedmethyl"

if [[ ! -x "${BIN}" ]]; then
    echo "Build first: cd ${PROJECT_ROOT} && make" >&2
    exit 1
fi

MODKIT_DIR="${MODKIT_DIR:-../results/modkit}"
OUTPUT="${OUTPUT:-merged_methylation.tsv}"
MIN_COV="${MIN_COV:-3}"
MIN_SAMPLES="${MIN_SAMPLES:-}"

if [[ -z "${SAMPLES:-}" ]]; then
    echo "Set SAMPLES (space-separated IDs) and optionally MODKIT_DIR, OUTPUT, MIN_COV, MIN_SAMPLES" >&2
    exit 1
fi

args=(--hap -c "${MIN_COV}")
if [[ -n "${MIN_SAMPLES}" ]]; then
    args+=(-s "${MIN_SAMPLES}")
fi
args+=("${OUTPUT}")

for id in ${SAMPLES}; do
    args+=("${id}" "${MODKIT_DIR}/${id}.bed/${id}_hp1.bedmethyl"
            "${MODKIT_DIR}/${id}.bed/${id}_hp2.bedmethyl")
done

exec "${BIN}" "${args[@]}"
