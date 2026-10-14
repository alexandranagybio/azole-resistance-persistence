# Neutral two-marker population model

## Purpose

Script 14 checks clonal and two-parent inheritance in a population
before adding selection or dispersal.

## Model

- One population of 500 haploid individuals.
- Initially 450 wild type and 50 double-marker carriers.
- All four genotypes have equal fitness.
- No new mutation or migration.
- Cycle 1 records the founders.
- Cycles 2–500 each produce exactly 500 offspring.
- Parents die after reproduction; generations do not overlap.

SEXUAL_FRAC is the hypothetical probability of two-parent inheritance.
Otherwise, an offspring is a clone of one uniformly sampled parent.

Two-parent events sample two distinct individuals uniformly.
RECOMB_FRAC specifies the probability of recombinant inheritance
between the markers within such an event.

Mating types, tetrads, and the complete fungal sexual cycle are not
represented. Marker coordinates are abstract.

The default RECOMB_FRAC of 0.00192 follows the cross-specific total
recombinant fraction inferred by Auxier et al. (2023):
https://doi.org/10.1371/journal.pbio.3002278

It is not an estimate of how frequently sex occurs in nature.
Applying it to other parental backgrounds is a modelling assumption.

## Diagnostic runs

All three runs use seed 1001.

| Scenario | SEXUAL_FRAC | RECOMB_FRAC |
|---|---:|---:|
| Clonal | 0 | 0.00192 |
| Two-parent, no recombination | 1 | 0 |
| Two-parent, with recombination | 1 | 0.00192 |

These settings isolate inheritance mechanisms; they are not estimates
of natural reproductive behaviour.

Neither control produced single-marker genotypes.

With recombination, a TR34-only individual first appeared at cycle 4.
At that cycle there were 444 wild-type individuals, one TR34-only
individual, no L98H-only individuals, and 55 double carriers.

Both markers were lost by cycle 500 in all three diagnostic runs.
These single-run outcomes do not estimate loss probabilities.

Offspring are independent inheritance draws. Reciprocal recombinant
genotypes therefore need not appear together in the same cycle.

## Validation

The model checks population size, removal of older parents, genotype
frequency totals, and absence of single-marker genotypes in controls.

The saved R checker additionally checks all 500 cycle records,
initial conditions, frequency ranges and increments, marker-frequency
calculations, and absorbing loss and fixation of each marker.

Recombinant appearance is reported, not required for a run to pass.
A valid stochastic run can lose the markers before recombinants arise.

## Reproduce

Run from the project root. Output files must not already exist.

```bash
mkdir -p data/checks/neutral_two_markers

~/SLiM/build/slim -s 1001 \
  -d SEXUAL_FRAC=0.0 \
  -d RECOMB_FRAC=0.00192 \
  -d 'OUTFILE="data/checks/neutral_two_markers/clonal_seed_1001.csv"' \
  models/14_neutral_two_marker_population.slim

~/SLiM/build/slim -s 1001 \
  -d SEXUAL_FRAC=1.0 \
  -d RECOMB_FRAC=0.0 \
  -d 'OUTFILE="data/checks/neutral_two_markers/two_parent_no_recomb_seed_1001.csv"' \
  models/14_neutral_two_marker_population.slim

~/SLiM/build/slim -s 1001 \
  -d SEXUAL_FRAC=1.0 \
  -d RECOMB_FRAC=0.00192 \
  -d 'OUTFILE="data/checks/neutral_two_markers/two_parent_recomb_seed_1001.csv"' \
  models/14_neutral_two_marker_population.slim

Rscript analysis/check_neutral_two_markers.R
```

Generated diagnostic CSVs are ignored by Git.
