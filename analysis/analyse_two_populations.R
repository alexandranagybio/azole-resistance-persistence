# Analyse the two-patch, zero-cost dispersal experiment.
# Run from the project root with:
# Rscript analysis/analyse_two_populations.R

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2)
  stop("Usage: Rscript analysis/analyse_two_populations.R [input_root] [output_root]")

input_root <- if (length(args) >= 1) args[1] else
  "data/two_population_dispersal_100"
output_root <- if (length(args) >= 2) args[2] else
  "data/summaries/two_population_dispersal_100"
rates <- c("0.0", "0.001", "0.01")
seeds <- 1001:1100

expected_names <- c(
  "generation", "frequency_p1", "frequency_p2", "frequency_total"
)

validate_run <- function(d, rate, label) {
  fail <- function(message) stop(label, ": ", message, call. = FALSE)

  if (!identical(names(d), expected_names))
    fail("unexpected CSV columns")
  if (!identical(d$generation, 1:500))
    fail("expected cycles 1 through 500")
  if (!all(vapply(d, is.numeric, logical(1))))
    fail("non-numeric values")

  f <- as.matrix(d[, expected_names[-1]])
  if (any(!is.finite(f)) || any(f < 0 | f > 1))
    fail("invalid frequency")

  if (any(abs(f[1, ] - c(0.1, 0, 0.05)) > 1e-8))
    fail("incorrect starting frequencies")

  if (any(abs(d$frequency_total -
              (d$frequency_p1 + d$frequency_p2) / 2) > 1e-6))
    fail("overall frequency does not match patch frequencies")

  # Each patch contains 500 haploids: frequencies change in units of 1/500.
  for (column in c("frequency_p1", "frequency_p2")) {
    counts <- d[[column]] * 500
    if (any(abs(counts - round(counts)) > 1e-6))
      fail("patch frequencies are inconsistent with 500 haploids")
  }

  if (rate == 0 && any(d$frequency_p2 != 0))
    fail("resistance reached p2 without dispersal")

  # Global loss/fixation is absorbing; local loss/fixation need not be.
  for (boundary in c(0, 1)) {
    hits <- which(d$frequency_total == boundary)
    if (length(hits) &&
        any(d$frequency_total[hits[1]:nrow(d)] != boundary))
      fail("global loss or fixation was reversed")
  }
}

rows <- list()
i <- 0L

for (rate in rates) {
  folder <- file.path(input_root, paste0("migration_", rate))
  expected_files <- paste0("seed_", seeds, ".csv")
  actual_files <- list.files(folder, pattern = "^seed_.*\\.csv$")

  if (!setequal(actual_files, expected_files))
    stop("Unexpected or missing seed CSVs in ", folder)

  for (seed in seeds) {
    path <- file.path(folder, paste0("seed_", seed, ".csv"))
    d <- read.csv(path)
    validate_run(d, as.numeric(rate), path)

    arrival <- d$generation[d$frequency_p2 > 0]
    final <- d[nrow(d), ]

    i <- i + 1L
    rows[[i]] <- data.frame(
      migration_rate = as.numeric(rate),
      seed = seed,
      arrived_p2 = length(arrival) > 0,
      first_arrival_cycle = if (length(arrival)) min(arrival) else NA_integer_,
      final_p1 = final$frequency_p1,
      final_p2 = final$frequency_p2,
      final_total = final$frequency_total
    )
  }
}

runs <- do.call(rbind, rows)

wilson <- function(successes, n) {
  z <- qnorm(0.975)
  p <- successes / n
  denominator <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denominator
  half_width <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) /
    denominator
  100 * c(max(0, centre - half_width), min(1, centre + half_width))
}

summaries <- list()
endpoints <- list()

for (rate in as.numeric(rates)) {
  r <- runs[runs$migration_rate == rate, ]
  n <- nrow(r)

  outcomes <- list(
    arrival_p2 = r$arrived_p2,
    persistence_p1 = r$final_p1 > 0,
    persistence_p2 = r$final_p2 > 0,
    persistence_anywhere = r$final_total > 0
  )

  for (metric in names(outcomes)) {
    successes <- sum(outcomes[[metric]])
    ci <- wilson(successes, n)
    summaries[[length(summaries) + 1L]] <- data.frame(
      migration_rate = rate,
      metric = metric,
      total_runs = n,
      successes = successes,
      percent = 100 * successes / n,
      ci_lower_percent = ci[1],
      ci_upper_percent = ci[2]
    )
  }

  endpoints[[length(endpoints) + 1L]] <- data.frame(
    migration_rate = rate,
    total_runs = n,
    globally_lost = sum(r$final_total == 0),
    globally_segregating = sum(r$final_total > 0 & r$final_total < 1),
    globally_fixed = sum(r$final_total == 1),
    mean_final_p1 = mean(r$final_p1),
    mean_final_p2 = mean(r$final_p2),
    mean_final_total = mean(r$final_total)
  )
}

summary <- do.call(rbind, summaries)
endpoint_summary <- do.call(rbind, endpoints)

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
write.csv(runs, file.path(output_root, "per_run_outcomes.csv"),
          row.names = FALSE)
write.csv(summary, file.path(output_root, "arrival_and_persistence.csv"),
          row.names = FALSE)
write.csv(endpoint_summary, file.path(output_root, "final_frequencies.csv"),
          row.names = FALSE)


# First arrival by withdrawal, after withdrawal, or never.
arrival_timing <- do.call(rbind, lapply(as.numeric(rates), function(rate) {
  r <- runs[runs$migration_rate == rate, ]
  arrival <- r$first_arrival_cycle
  arrived <- !is.na(arrival)

  before <- sum(arrived & arrival <= 50, na.rm = TRUE)
  after <- sum(arrived & arrival > 50, na.rm = TRUE)
  never <- sum(!arrived)

  stopifnot(
    before + after + never == nrow(r),
    all(r$arrived_p2 == arrived)
  )

  data.frame(
    migration_rate = rate,
    total_runs = nrow(r),
    first_arrival_by_cycle_50 = before,
    first_arrival_after_cycle_50 = after,
    never_arrived = never,
    median_first_arrival_if_arrived =
      if (any(arrived)) median(arrival[arrived]) else NA_real_
  )
}))

write.csv(arrival_timing,
          file.path(output_root, "arrival_timing.csv"),
          row.names = FALSE)

cat("All", nrow(runs), "trajectories passed validation.\n\n")
cat("Arrival and persistence (95% Wilson intervals):\n")
print(summary, row.names = FALSE, digits = 4)
cat("\nFinal outcomes at cycle 500:\n")
print(endpoint_summary, row.names = FALSE, digits = 4)

cat("\nFirst-arrival timing (withdrawal at cycle 50):\n")
print(arrival_timing, row.names = FALSE)
