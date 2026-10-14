# Run after models/13_TR34_L98H_stochastic_cross.slim.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1)
  stop("Usage: Rscript analysis/check_stochastic_cross.R [csv_path]")

path <- if (length(args)) args[1] else
  "data/checks/stochastic_cross/seed_1001.csv"

d <- read.csv(path)

stopifnot(
  identical(names(d), c(
    "genotype", "count", "frequency", "expected_frequency"
  )),
  identical(d$genotype, c(
    "wild_type", "TR34_only", "L98H_only", "double"
  ))
)

values <- as.matrix(d[, -1])
stopifnot(
  is.numeric(values),
  all(is.finite(values)),
  all(d$count >= 0),
  all(d$count == round(d$count))
)

n <- sum(d$count)
stopifnot(
  n > 0,
  all(d$frequency >= 0 & d$frequency <= 1),
  all(d$expected_frequency >= 0 & d$expected_frequency <= 1),
  abs(sum(d$expected_frequency) - 1) < 1e-6,
  all(abs(d$frequency - d$count / n) < 1e-6)
)

recombinant <- d$genotype %in% c("wild_type", "double")
k <- sum(d$count[recombinant])
configured <- sum(d$expected_frequency[recombinant])
ci <- binom.test(k, n)$conf.int

d$expected_count <- n * d$expected_frequency
print(d, row.names = FALSE)

cat("\nCSV checks passed.\n")
cat("Offspring:", n, "\n")
cat("Observed recombinant percentage:", 100 * k / n, "\n")
cat("Configured recombinant percentage:", 100 * configured, "\n")
cat("95% binomial interval (%):",
    100 * ci[1], "to", 100 * ci[2], "\n")

# Interval coverage is descriptive, not a pass/fail requirement:
# a correct stochastic implementation can produce an interval
# that excludes its configured probability.
