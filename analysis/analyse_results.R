# Run from the project root:
# Rscript analysis/analyse_results.R
# Reads saved simulations; does not launch SLiM.

args <- commandArgs(trailingOnly = TRUE)
raw_dir <- if (length(args) >= 1) args[1] else file.path("data", "raw")
n_runs <- if (length(args) >= 2) as.numeric(args[2]) else 20

if (length(n_runs) != 1 || !is.finite(n_runs) ||
    n_runs < 1 || n_runs != floor(n_runs))
  stop("Replicate count must be a positive integer.")

seeds <- 1000L + seq_len(n_runs)
summary_dir <- if (length(args) >= 3) args[3] else file.path("data", "summaries")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(raw_dir))
  stop("Run this script from the project root.")

read_runs <- function(folder, end_generation) {
  paths <- file.path(raw_dir, folder, paste0("seed_", seeds, ".csv"))

  if (!all(file.exists(paths)))
    stop("Missing seed files in: ", folder)

  runs <- lapply(paths, read.csv)

  for (i in seq_along(runs)) {
    d <- runs[[i]]

    valid <- identical(names(d), c("generation", "frequency")) &&
      is.numeric(d$generation) &&
      is.numeric(d$frequency) &&
      identical(d$generation, seq_len(end_generation)) &&
      all(is.finite(d$frequency)) &&
      all(d$frequency >= 0 & d$frequency <= 1) &&
      d$frequency[1] == 0.1

    if (!valid)
      stop("Invalid trajectory: ", paths[i])
  }

  x <- do.call(cbind, lapply(runs, function(d) d$frequency))
  colnames(x) <- seeds
  x
}

short <- list(
  "0%" = read_runs("replicates_cost_000", 100),
  "1%" = read_runs("replicates_cost_001", 100),
  "5%" = read_runs("replicates_cost_005", 100)
)

long <- list(
  "0%" = read_runs("long_replicates_cost_0.0", 500),
  "1%" = read_runs("long_replicates_cost_0.01", 500),
  "5%" = read_runs("long_replicates_cost_0.05", 500)
)

exposure <- list(
  "25" = read_runs("exposure_duration_25", 475),
  "50" = read_runs("exposure_duration_50", 500),
  "100" = read_runs("exposure_duration_100", 550)
)

# Same seeds must reproduce the short runs when extended.
for (label in names(short)) {
  stopifnot(all(short[[label]] == long[[label]][1:100, ]))
  stopifnot(all(long[[label]][1:50, ] == long[["0%"]][1:50, ]))
}

# Exposure scenarios must match before the earliest withdrawal.
stopifnot(
  all(exposure[["25"]][1:25, ] == exposure[["50"]][1:25, ]),
  all(exposure[["25"]][1:25, ] == exposure[["100"]][1:25, ]),
  all(exposure[["50"]] == long[["0%"]])
)

summarise_endpoint <- function(x) {
  final <- x[nrow(x), ]
  data.frame(
    total_runs = length(final),
    lost = sum(final == 0),
    segregating = sum(final > 0 & final < 1),
    fixed = sum(final == 1),
    persistence_percent = 100 * mean(final > 0),
    mean_final_frequency = mean(final)
  )
}

summarise_costs <- function(runs) {
  do.call(rbind, lapply(names(runs), function(label) {
    data.frame(cost = label, summarise_endpoint(runs[[label]]))
  }))
}

summary100 <- summarise_costs(short)
summary500 <- summarise_costs(long)

summary_exposure <- do.call(rbind, lapply(names(exposure), function(label) {
  withdrawal <- as.integer(label)
  x <- exposure[[label]]
  data.frame(
    withdrawal_generation = withdrawal,
    mean_frequency_at_withdrawal = mean(x[withdrawal, ]),
    summarise_endpoint(x)
  )
}))

# Wilson 95% intervals for the proportion of runs retaining resistance.
add_persistence_ci <- function(d) {
  n <- d$total_runs
  persisted <- d$segregating + d$fixed
  p <- persisted / n
  z <- qnorm(0.975)
  denominator <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denominator
  half_width <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) /
    denominator

  d$persisted <- persisted
  d$persistence_ci_lower_percent <- 100 * pmax(0, centre - half_width)
  d$persistence_ci_upper_percent <- 100 * pmin(1, centre + half_width)
  d
}

summary100 <- add_persistence_ci(summary100)
summary500 <- add_persistence_ci(summary500)
summary_exposure <- add_persistence_ci(summary_exposure)

# Confirm that the original seeds reproduce yesterday's trajectories.
if (normalizePath(raw_dir) != normalizePath("data/raw") && n_runs >= 20) {
  original_paths <- list.files(
    "data/raw", pattern = "^seed_[0-9]+\\.csv$", recursive = TRUE
  )
  stopifnot(length(original_paths) == 180L)

  for (relative_path in original_paths) {
    original <- read.csv(file.path("data/raw", relative_path))
    current <- read.csv(file.path(raw_dir, relative_path))
    if (!identical(original, current))
      stop("Original seed differs: ", relative_path)
  }
  cat("All 180 original trajectories reproduced exactly.\n")
}

# Keep the original summary files; save regenerated versions separately.
write.csv(summary100,
          file.path(summary_dir, "recomputed_generation_100.csv"),
          row.names = FALSE)
write.csv(summary500,
          file.path(summary_dir, "recomputed_generation_500.csv"),
          row.names = FALSE)
write.csv(summary_exposure,
          file.path(summary_dir, "recomputed_exposure_duration.csv"),
          row.names = FALSE)

cat("All", 9 * n_runs, "trajectories passed validation.\n")
cat("\nCost comparison at generation 100:\n")
print(summary100, row.names = FALSE)
cat("\nCost comparison at generation 500:\n")
print(summary500, row.names = FALSE)
cat("\nExposure-duration comparison:\n")
print(summary_exposure, row.names = FALSE)
