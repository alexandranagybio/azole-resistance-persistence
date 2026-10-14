# Two-marker inheritance checks

These isolated checks examine inheritance of two neutral markers.
They do not simulate a complete fungal life cycle.

## Controlled inheritance: script 12

The parents are wild type and a double-marker carrier.

Without a crossover, the two constructed offspring retain the
parental genotypes: one wild type and one double carrier.

With a forced crossover between the markers, the two constructed
offspring are TR34-only and L98H-only.

The script checks these genotype counts directly.

## Stochastic inheritance: script 13

The parents are TR34-only and L98H-only. Each offspring is an
independent haploid inheritance draw.

The configured recombinant fraction is 0.00192. Expected offspring
probabilities are:

| Genotype | Probability |
|---|---:|
| Wild type | 0.00096 |
| TR34 only | 0.49904 |
| L98H only | 0.49904 |
| Double carrier | 0.00096 |

Wild type and double carrier are the recombinant classes in this
cross. The implementation checks that their combined count equals
the number of recombinant draws.

With seed 1001 and 100,000 offspring, the observed counts were
91 wild type, 49,826 TR34-only, 49,976 L98H-only and 107 double carriers.

The observed recombinant fraction was 0.198%, with a 95% exact
binomial interval of approximately 0.171–0.228%.

This is consistent with the configured probability. It checks
implementation, not independent biological validity.

## Scope

Sexual reproduction frequency and recombination fraction are
different parameters.

Sexual reproduction frequency determines how often offspring arise
through sexual reproduction. The recombination fraction describes
inheritance between the markers conditional on that reproduction.

This check addresses inheritance only. It does not estimate how
often sexual reproduction occurs in nature.

Independent offspring draws do not explicitly represent tetrads,
mating-type compatibility, or the complete fungal life cycle.
Marker positions are illustrative.

## Reproduce

Run from the project root:

```bash
~/SLiM/build/slim -s 1001 \
  models/12_TR34_L98H_inheritance_check.slim

mkdir -p data/checks/stochastic_cross

~/SLiM/build/slim -s 1001 \
  -d 'OUTFILE="data/checks/stochastic_cross/seed_1001.csv"' \
  models/13_TR34_L98H_stochastic_cross.slim

Rscript analysis/check_stochastic_cross.R
```

Script 13 refuses an existing output file. For another run, choose
a fresh OUTFILE and pass that same path to the R checker.
