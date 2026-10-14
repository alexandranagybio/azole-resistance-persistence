# Run from the project root:
# Rscript analysis/analyse_neutral_two_markers.R

input_root <- "data/neutral_two_markers_100"
output_root <- "data/summaries/neutral_two_markers_100"
fractions <- c("0.0", "0.01", "0.1", "1.0")
seeds <- 1001:1100

columns <- c(
  "generation", "population_size",
  "wt", "tr34_only", "l98h_only", "double",
  "tr34_frequency", "l98h_frequency"
)

expected_files <- unlist(lapply(fractions, function(q) {
  file.path(paste0("sexual_", q), paste0("seed_", seeds, ".csv"))
}))

actual_files <- list.files(
  input_root, pattern = "^seed_.*\\.csv$", recursive = TRUE
)

stopifnot(
  length(actual_files) == 400L,
  setequal(actual_files, expected_files)
)

validate <- function(d, q, label) {
  fail <- function(message) stop(label, ": ", message, call. = FALSE)

  if (!identical(names(d), columns))
    fail("unexpected columns")
  if (!all(vapply(d, is.numeric, logical(1))))
    fail("non-numeric values")
  if (nrow(d) != 500L || any(!is.finite(as.matrix(d))) ||
      any(d$generation != seq_len(500)))
    fail("invalid cycle records")
  if (any(d$population_size != 500))
    fail("population size differs from 500")

  f <- as.matrix(d[, 3:8])
  if (any(f < 0 | f > 1))
    fail("frequency outside [0, 1]")
  if (any(abs(f * 500 - round(f * 500)) > 1e-6))
    fail("frequency inconsistent with 500 haploids")

  g <- as.matrix(d[, c("wt", "tr34_only", "l98h_only", "double")])
  if (any(abs(rowSums(g) - 1) > 1e-6))
    fail("genotypes do not sum to one")
  if (any(abs(g[1, ] - c(0.9, 0, 0, 0.1)) > 1e-6))
    fail("incorrect initial conditions")

  if (any(abs(d$tr34_frequency -
              d$tr34_only - d$double) > 1e-6) ||
      any(abs(d$l98h_frequency -
              d$l98h_only - d$double) > 1e-6))
    fail("marker frequencies disagree with genotypes")

  if (q == 0 && any(d$tr34_only > 0 | d$l98h_only > 0))
    fail("single-marker genotype in clonal control")

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
}

outcomes <- list()
trajectories <- list()
index <- 0L

for (q in fractions) {
  for (seed in seeds) {
    relative <- file.path(
      paste0("sexual_", q), paste0("seed_", seed, ".csv")
    )
    d <- read.csv(file.path(input_root, relative))
    validate(d, as.numeric(q), relative)

    single <- d$tr34_only > 0 | d$l98h_only > 0
    last <- d[nrow(d), ]
    association <- d$double - d$tr34_frequency * d$l98h_frequency

    index <- index + 1L
    outcomes[[index]] <- data.frame(
      sexual_fraction = as.numeric(q),
      seed = seed,
      singles_ever = any(single),
      first_single_cycle = if (any(single))
        d$generation[which(single)[1]] else NA_real_,
      final_tr34 = last$tr34_frequency,
      final_l98h = last$l98h_frequency,
      final_double = last$double,
      final_single_frequency = last$tr34_only + last$l98h_only,
      final_D = tail(association, 1)
    )

    trajectories[[index]] <- data.frame(
      sexual_fraction = as.numeric(q),
      generation = d$generation,
      tr34_frequency = d$tr34_frequency,
      l98h_frequency = d$l98h_frequency,
      double_frequency = d$double,
      single_frequency = d$tr34_only + d$l98h_only,
      D = association
    )
  }
}

runs <- do.call(rbind, outcomes)

wilson <- function(k, n) {
  z <- qnorm(0.975)
  p <- k / n
  denominator <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denominator
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) /
    denominator
  100 * c(max(0, centre - half), min(1, centre + half))
}

summary <- do.call(rbind, lapply(as.numeric(fractions), function(q) {
  r <- runs[runs$sexual_fraction == q, ]
  k <- sum(r$singles_ever)
  ci <- wilson(k, nrow(r))

  data.frame(
    sexual_fraction = q,
    runs = nrow(r),
    singles_ever = k,
    singles_ever_percent = 100 * k / nrow(r),
    ci_lower_percent = ci[1],
    ci_upper_percent = ci[2],
    median_first_cycle_if_observed = if (k > 0)
      median(r$first_single_cycle, na.rm = TRUE) else NA_real_,
    tr34_present_at_500 = sum(r$final_tr34 > 0),
    l98h_present_at_500 = sum(r$final_l98h > 0),
    both_markers_present_at_500 = sum(
      r$final_tr34 > 0 & r$final_l98h > 0
    ),
    double_carriers_present_at_500 = sum(r$final_double > 0),
    mean_final_tr34 = mean(r$final_tr34),
    mean_final_l98h = mean(r$final_l98h),
    mean_final_D = mean(r$final_D)
  )
}))

mean_trajectories <- aggregate(
  cbind(tr34_frequency, l98h_frequency,
        double_frequency, single_frequency, D) ~
    sexual_fraction + generation,
  data = do.call(rbind, trajectories),
  FUN = mean
)

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

write.csv(runs, file.path(output_root, "per_run_outcomes.csv"),
          row.names = FALSE)
write.csv(summary, file.path(output_root, "summary.csv"),
          row.names = FALSE)
write.csv(mean_trajectories,
          file.path(output_root, "mean_trajectories.csv"),
          row.names = FALSE)

cat("All 400 neutral trajectories passed validation.\n\n")
print(summary, row.names = FALSE, digits = 4)
cat("\nSaved summaries in", output_root, "\n")
