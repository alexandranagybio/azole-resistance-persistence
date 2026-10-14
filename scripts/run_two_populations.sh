#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$project_root"

slim_bin="${SLIM_BIN:-slim}"
output_root="${1:-data/two_population_dispersal_100}"

if (( $# > 1 )); then
  printf 'Usage: bash scripts/run_two_populations.sh [output_directory]\n' >&2
  exit 1
fi

if ! command -v "$slim_bin" >/dev/null 2>&1; then
  printf 'SLiM executable not found: %s\n' "$slim_bin" >&2
  exit 1
fi

# Refuse an existing destination to protect previous results.
mkdir -p "$(dirname "$output_root")"
mkdir "$output_root"

for rate in 0.0 0.001 0.01; do
  destination="$output_root/migration_$rate"
  mkdir "$destination"

  for seed in {1001..1100}; do
    "$slim_bin" \
      -s "$seed" \
      -d COST=0.0 \
      -d EXPOSURE_END=50 \
      -d "MIGRATION_RATE=$rate" \
      -d "OUTFILE=\"$destination/seed_$seed.csv\"" \
      models/10_TR34_L98H_two_populations.slim \
      > "$destination/seed_$seed.log" 2>&1

    test -s "$destination/seed_$seed.csv"
  done

  printf 'Completed migration rate %s: 100 runs\n' "$rate"
done

printf 'Finished all 300 simulations in %s\n' "$output_root"
