#!/usr/bin/env bash
set -euo pipefail

# Find the project root regardless of where this script is called.
project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$project_root"

slim_bin="${SLIM_BIN:-slim}"
command -v "$slim_bin" >/dev/null 2>&1 || {
  printf 'SLiM executable not found: %s\n' "$slim_bin" >&2
  exit 1
}

# Replicate count is the second argument; defaults to 20.
replicates="${2:-20}"
if ! [[ "$replicates" =~ ^[1-9][0-9]*$ ]]; then
  printf 'Replicate count must be a positive integer.\n' >&2
  exit 1
fi
last_seed=$((1000 + replicates))

# Refuse an existing destination to protect previous results.
output_root="${1:-data/reproduced}"
mkdir -p "$(dirname "$output_root")"
mkdir "$output_root"

run_scenario() {
  local model="$1"
  local folder="$2"
  shift 2

  local destination="$output_root/$folder"
  mkdir "$destination"

  for ((seed=1001; seed<=last_seed; seed++)); do
    "$slim_bin" \
      -s "$seed" \
      "$@" \
      -d "OUTFILE=\"$destination/seed_$seed.csv\"" \
      "$model" \
      > "$destination/seed_$seed.log" 2>&1

    test -s "$destination/seed_$seed.csv"
  done

  printf 'Completed %s: %s runs\n' "$folder" "$replicates"
}

for scenario in "000 0.0" "001 0.01" "005 0.05"; do
  read -r label cost <<< "$scenario"
  run_scenario \
    models/07_TR34_L98H_persistence.slim \
    "replicates_cost_$label" \
    -d "COST=$cost"
done

for cost in 0.0 0.01 0.05; do
  run_scenario \
    models/08_TR34_L98H_long_persistence.slim \
    "long_replicates_cost_$cost" \
    -d "COST=$cost"
done

for exposure in 25 50 100; do
  run_scenario \
    models/09_TR34_L98H_exposure_duration.slim \
    "exposure_duration_$exposure" \
    -d COST=0.0 \
    -d "EXPOSURE_END=$exposure"
done

printf 'Finished all %s simulations in %s\n' \
  "$((9 * replicates))" "$output_root"
