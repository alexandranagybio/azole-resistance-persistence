#!/usr/bin/env bash
set -euo pipefail

# Usage:
# SLIM_BIN="$HOME/SLiM/build/slim" \
#   bash scripts/run_neutral_two_markers.sh [output_directory]

if (( $# > 1 )); then
  printf 'Usage: %s [output_directory]\n' "$0" >&2
  exit 1
fi

project_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$project_root"

slim_bin="${SLIM_BIN:-slim}"
output_root="${1:-data/neutral_two_markers_100}"
model="models/14_neutral_two_marker_population.slim"

if ! command -v "$slim_bin" >/dev/null 2>&1; then
  printf 'SLiM executable not found: %s\n' "$slim_bin" >&2
  exit 1
fi

if [[ ! -f "$model" ]]; then
  printf 'Model not found: %s\n' "$model" >&2
  exit 1
fi

if [[ -e "$output_root" || -L "$output_root" ]]; then
  printf 'Destination already exists: %s\nChoose a fresh directory.\n' \
    "$output_root" >&2
  exit 1
fi

# OUTFILE is passed as an Eidos string literal.
if [[ "$output_root" == *\"* ||
      "$output_root" == *\\* ||
      "$output_root" == *$'\n'* ||
      "$output_root" == *$'\r'* ]]; then
  printf 'Output path cannot contain quotes, backslashes or line breaks.\n' >&2
  exit 1
fi

mkdir -p -- "$(dirname -- "$output_root")"
mkdir -- "$output_root"

printf 'sexual_fraction,recombination_fraction,seed,csv\n' \
  > "$output_root/manifest.csv"

for fraction in 0.0 0.01 0.1 1.0; do
  destination="$output_root/sexual_$fraction"
  mkdir -- "$destination"

  for seed in {1001..1100}; do
    csv="$destination/seed_$seed.csv"
    log="$destination/seed_$seed.log"

    if ! "$slim_bin" \
      -s "$seed" \
      -d "SEXUAL_FRAC=$fraction" \
      -d RECOMB_FRAC=0.00192 \
      -d "OUTFILE=\"$csv\"" \
      "$model" > "$log" 2>&1; then
      printf 'Simulation failed. See: %s\n' "$log" >&2
      exit 1
    fi

    if [[ ! -s "$csv" ]]; then
      printf 'Missing or empty output: %s\n' "$csv" >&2
      exit 1
    fi

    printf '%s,0.00192,%s,sexual_%s/seed_%s.csv\n' \
      "$fraction" "$seed" "$fraction" "$seed" \
      >> "$output_root/manifest.csv"
  done

  printf 'Completed two-parent fraction %s: 100 runs\n' "$fraction"
done

printf 'Finished all 400 simulations in %s\n' "$output_root"
