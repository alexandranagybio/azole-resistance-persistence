# Run from the project root:
# Rscript analysis/check_neutral_two_markers.R

root <- "data/checks/neutral_two_markers"

files <- c(
  clonal = "clonal_seed_1001.csv",
  two_parent_no_recomb = "two_parent_no_recomb_seed_1001.csv",
  two_parent_recomb = "two_parent_recomb_seed_1001.csv"
)

expected_columns <- c(
  "generation", "population_size",
  "wt", "tr34_only", "l98h_only", "double",
  "tr34_frequency", "l98h_frequency"
)

results <- lapply(names(files), function(scenario) {
  d <- read.csv(file.path(root, files[[scenario]]))
  fail <- function(message)
    stop(scenario, ": ", message, call. = FALSE)

  if (!identical(names(d), expected_columns))
    fail("unexpected columns")

  if (!all(vapply(d, is.numeric, logical(1))))
    fail("non-numeric values")

  if (nrow(d) != 500L || any(!is.finite(as.matrix(d))) ||
      any(d$generation != seq_len(500)))
    fail("invalid cycle records")

  if (any(d$population_size != 500))
    fail("population size differs from 500")

  frequencies <- as.matrix(d[, 3:8])
  if (any(frequencies < 0 | frequencies > 1))
    fail("frequency outside [0, 1]")

  if (any(abs(frequencies * 500 -
              round(frequencies * 500)) > 1e-6))
    fail("frequencies inconsistent with 500 haploids")

  genotypes <- as.matrix(d[, c(
    "wt", "tr34_only", "l98h_only", "double"
  )])

  if (any(abs(rowSums(genotypes) - 1) > 1e-6))
    fail("genotype frequencies do not sum to one")

  if (any(abs(genotypes[1, ] - c(0.9, 0, 0, 0.1)) > 1e-6))
    fail("incorrect initial genotypes")

  if (any(abs(d$tr34_frequency -
              (d$tr34_only + d$double)) > 1e-6) ||
      any(abs(d$l98h_frequency -
              (d$l98h_only + d$double)) > 1e-6))
    fail("marker frequencies disagree with genotype frequencies")

  # Without mutation or migration, an individual marker cannot
  # return after loss, or disappear after fixation.
  for (column in c("tr34_frequency", "l98h_frequency")) {
    x <- d[[column]]

    for (boundary in c(0, 1)) {
      hits <- which(x == boundary)
      if (length(hits)) {
        later <- seq.int(hits[1], length(x))
        if (any(x[later] != boundary))
          fail(paste(column, "left an absorbing boundary"))
      }
    }
  }

  singles <- d$tr34_only > 0 | d$l98h_only > 0

  if (scenario != "two_parent_recomb" && any(singles))
    fail("single-marker genotype appeared in a control")

  # Appearance is descriptive, not a required stochastic outcome.
  first <- if (any(singles))
    d$generation[which(singles)[1]] else NA_real_

  data.frame(
    scenario = scenario,
    first_single_marker_cycle = first,
    final_tr34_frequency = tail(d$tr34_frequency, 1),
    final_l98h_frequency = tail(d$l98h_frequency, 1)
  )
})

cat("All three neutral trajectories passed validation.\n\n")
print(do.call(rbind, results), row.names = FALSE)
