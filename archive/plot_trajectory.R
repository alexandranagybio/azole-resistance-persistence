trajectory <- read.csv("~/SLiM/changing_exposure.csv")

plot(
  trajectory$generation,
  trajectory$frequency,
  type = "l",
  ylim = c(0, 1),
  xlab = "Generation",
  ylab = "Allele frequency"
)

abline(v = 50, lty = 2, col = "red")

trajectory <- read.csv("~/SLiM/changing_exposure.csv")

# Check that all 100 generations were saved
nrow(trajectory)
head(trajectory)
tail(trajectory)

plot(
  trajectory$generation,
  trajectory$frequency,
  type = "l",
  lwd = 2,
  ylim = c(0, 1),
  xlab = "Generation",
  ylab = "Allele frequency"
)

abline(v = 50, lty = 2, col = "red")
file.copy(
  "~/SLiM/changing_exposure.csv",
  "~/SLiM/changing_exposure_run01.csv",
  overwrite = FALSE
)

file.copy(
  "~/SLiM/changing_exposure.csv",
  "~/SLiM/changing_exposure_run02.csv",
  overwrite = FALSE
)

run1 <- read.csv("~/SLiM/changing_exposure_run01.csv")
run2 <- read.csv("~/SLiM/changing_exposure_run02.csv")

plot(
  run1$generation, run1$frequency,
  type = "l", lwd = 2, col = "steelblue",
  ylim = c(0, 1),
  xlab = "Generation", ylab = "Allele frequency"
)

lines(
  run2$generation, run2$frequency,
  lwd = 2, col = "darkorange"
)

abline(v = 50, lty = 2, col = "red")

legend(
  "topright",
  legend = c("Run 1", "Run 2"),
  col = c("steelblue", "darkorange"),
  lwd = 2, bty = "n"
)

slim <- path.expand("~/SLiM/build/slim")
model <- path.expand("~/SLiM/05_changing_exposure.slim")
csv <- path.expand("~/SLiM/changing_exposure.csv")
outdir <- path.expand("~/SLiM/replicates")

stopifnot(file.exists(slim), file.exists(model))
dir.create(outdir, showWarnings = FALSE)

slim <- path.expand("~/SLiM/build/slim")
csv <- path.expand("~/SLiM/changing_exposure.csv")
model <- path.expand("~/SLiM/06_exposure_no_cost.slim")
outdir <- path.expand("~/SLiM/replicates_no_cost")
dir.create(outdir, showWarnings = FALSE)

# Recorded seeds make each run reproducible
seeds <- 1001:1020
runs <- vector("list", length(seeds))

for (i in seq_along(seeds)) {
  destination <- file.path(
    outdir, sprintf("run_%02d.csv", i)
  )
  logfile <- file.path(
    outdir, sprintf("run_%02d.log", i)
  )

  # Avoid replacing an existing set of results
  if (file.exists(destination) || file.exists(logfile))
    stop("Results already exist for run ", i)

  # Remove the previous temporary CSV
  if (file.exists(csv)) unlink(csv)

  status <- system2(
    slim,
    args = c("-s", seeds[i], shQuote(model)),
    stdout = logfile,
    stderr = logfile
  )

  if (status != 0 || !file.exists(csv))
    stop("Run ", i, " failed. Check: ", logfile)

  dat <- read.csv(csv)

  stopifnot(
    nrow(dat) == 100,
    all(dat$generation == 1:100),
    all(is.finite(dat$frequency)),
    all(dat$frequency >= 0 & dat$frequency <= 1),
    dat$frequency[1] == 0.1
  )

  if (!file.copy(csv, destination, overwrite = FALSE))
    stop("Could not save run ", i)

  dat$replicate <- i
  dat$seed <- seeds[i]
  runs[[i]] <- dat

  message("Finished run ", i, " of ", length(seeds))
}

results <- do.call(rbind, runs)

plot(
  NA, xlim = c(1, 100), ylim = c(0, 1),
  xlab = "Generation", ylab = "Allele frequency"
)

for (dat in runs) {
  lines(
    dat$generation, dat$frequency,
    col = adjustcolor("steelblue", alpha.f = 0.35)
  )
}

abline(v = 50, lty = 2, col = "red")

final <- subset(results, generation == 100)

c(
  persisted = sum(final$frequency > 0),
  lost = sum(final$frequency == 0),
  fixed = sum(final$frequency == 1)
)

write.csv(
  results,
  file.path(outdir, "all_trajectories.csv"),
  row.names = FALSE
)

write.csv(
  final,
  file.path(outdir, "final_outcomes.csv"),
  row.names = FALSE
)

cost <- read.csv(
  "~/SLiM/replicates/final_outcomes.csv"
)

no_cost <- read.csv(
  "~/SLiM/replicates_no_cost/final_outcomes.csv"
)

data.frame(
  scenario = c("With cost", "No cost"),
  mean_frequency = c(
    mean(cost$frequency),
    mean(no_cost$frequency)
  ),
  median_frequency = c(
    median(cost$frequency),
    median(no_cost$frequency)
  )
)

boxplot(
  list(
    "With cost" = cost$frequency,
    "No cost" = no_cost$frequency
  ),
  ylim = c(0, 1),
  ylab = "Allele frequency at generation 100",
  col = c("salmon", "skyblue")
)

stripchart(
  list(
    "With cost" = cost$frequency,
    "No cost" = no_cost$frequency
  ),
  vertical = TRUE,
  method = "stack",
  add = TRUE,
  pch = 16
)

dev.copy(
  png,
  filename = path.expand("~/SLiM/final_frequency_comparison.png"),
  width = 1600,
  height = 1200,
  res = 180
)
dev.off()

cost_all <- read.csv(
  "~/SLiM/replicates/all_trajectories.csv"
)

no_cost_all <- read.csv(
  "~/SLiM/replicates_no_cost/all_trajectories.csv"
)

before_cost <- subset(cost_all, generation <= 50)
before_no_cost <- subset(no_cost_all, generation <= 50)

# Match rows by seed and generation
before_cost <- before_cost[
  order(before_cost$seed, before_cost$generation), ]

before_no_cost <- before_no_cost[
  order(before_no_cost$seed, before_no_cost$generation), ]

identical(before_cost, before_no_cost)
mean_cost <- aggregate(
  frequency ~ generation, cost_all, mean
)

mean_cost <- aggregate(
  frequency ~ generation, cost_all, mean
)

mean_no_cost <- aggregate(
  frequency ~ generation, no_cost_all, mean
)

plot(
  mean_cost$generation, mean_cost$frequency,
  type = "l", col = "firebrick", lwd = 3,
  ylim = c(0, 1),
  xlab = "Generation",
  ylab = "Mean allele frequency"
)

lines(
  mean_no_cost$generation, mean_no_cost$frequency,
  col = "steelblue", lwd = 3
)

abline(v = 50, lty = 2)

legend(
  "topleft",
  legend = c("Cost after exposure", "No cost after exposure"),
  col = c("firebrick", "steelblue"),
  lwd = 3, bty = "n"
)

dev.copy(
  png,
  filename = path.expand("~/SLiM/mean_trajectory_comparison.png"),
  width = 1600,
  height = 1200,
  res = 180
)
dev.off()

files <- c(
  "persistence_cost_0_seed_1001.csv",
  "persistence_cost_001_seed_1001.csv",
  "persistence_cost_005_seed_1001.csv"
)

runs <- lapply(file.path(path.expand("~/SLiM"), files), read.csv)
colours <- c("#0072B2", "#E69F00", "#CC79A7")

plot(
  runs[[1]]$generation, runs[[1]]$frequency,
  type = "n", ylim = c(0, 0.2),
  xlab = "Generation",
  ylab = "Resistance frequency",
  main = "Resistance persistence: seed 1001"
)

abline(v = 50, lty = 2, col = "grey50")

for (i in seq_along(runs)) {
  lines(
    runs[[i]]$generation, runs[[i]]$frequency,
    col = colours[i], lwd = 2
  )
}

legend(
  "topright",
  legend = c("Cost: 0%", "Cost: 1%", "Cost: 5%", "Exposure ends"),
  col = c(colours, "grey50"),
  lty = c(1, 1, 1, 2),
  lwd = c(2, 2, 2, 1),
  bty = "n"
)

seed1001 <- read.csv(
  path.expand("~/SLiM/persistence_cost_001_seed_1001.csv")
)
seed1002 <- read.csv(
  path.expand("~/SLiM/persistence_cost_001_seed_1002.csv")
)

plot(
  seed1001$generation, seed1001$frequency,
  type = "l", col = "#0072B2", lwd = 2,
  ylim = c(0, 1),
  xlab = "Generation", ylab = "Resistance frequency",
  main = "Same resistance cost (1%), different seeds"
)

lines(
  seed1002$generation, seed1002$frequency,
  col = "#D55E00", lwd = 2
)

abline(v = 50, lty = 2, col = "grey50")

legend(
  "topleft",
  legend = c("Seed 1001", "Seed 1002"),
  col = c("#0072B2", "#D55E00"),
  lty = 1, lwd = 2, bty = "n"
)

folder <- path.expand("~/SLiM/replicates_cost_001")
files <- file.path(folder, paste0("seed_", 1001:1020, ".csv"))

runs <- lapply(files, read.csv)

# Check that every run contains generations 1–100.
stopifnot(all(vapply(runs, function(d) {
  identical(as.integer(d$generation), 1:100) &&
    all(is.finite(d$frequency)) &&
    all(d$frequency >= 0 & d$frequency <= 1)
}, logical(1))))

frequencies <- sapply(runs, function(d) d$frequency)

# Thin lines show individual runs; thick line shows their mean.
matplot(
  1:100, frequencies,
  type = "l", lty = 1,
  col = adjustcolor("#0072B2", alpha.f = 0.35),
  ylim = c(0, 1),
  xlab = "Generation",
  ylab = "Resistance frequency",
  main = "20 replicate runs: resistance cost 1%"
)

lines(1:100, rowMeans(frequencies), col = "black", lwd = 3)
abline(v = 50, col = "grey50", lty = 2)

legend(
  "topleft",
  legend = c("Individual runs", "Mean frequency", "Exposure ends"),
  col = c("#0072B2", "black", "grey50"),
  lty = c(1, 1, 2),
  lwd = c(1, 3, 1),
  bty = "n"
)

# Persistence means at least one carrier remains at generation 100.
final <- frequencies[100, ]

results <- data.frame(
  seed = 1001:1020,
  frequency_50 = frequencies[50, ],
  frequency_100 = final,
  carriers_100 = round(500 * final),
  persisted = final > 0
)

print(results, row.names = FALSE)

cat("\nPersisted:", sum(final > 0), "of", length(final), "runs\n")
cat("Persistence percentage:", 100 * mean(final > 0), "%\n")
cat("Fixed at generation 100:", sum(final == 1), "runs\n")
cat("Mean final frequency:", mean(final), "\n")


read_runs <- function(folder) {
  files <- file.path(
    path.expand(folder),
    paste0("seed_", 1001:1020, ".csv")
  )
  runs <- lapply(files, read.csv)

  stopifnot(all(vapply(runs, function(d) {
    identical(as.integer(d$generation), 1:100) &&
      all(is.finite(d$frequency)) &&
      all(d$frequency >= 0 & d$frequency <= 1)
  }, logical(1))))

  sapply(runs, function(d) d$frequency)
}

cost1 <- read_runs("~/SLiM/replicates_cost_001")
cost5 <- read_runs("~/SLiM/replicates_cost_005")

stopifnot(all(cost1[1:50, ] == cost5[1:50, ]))

summarise_cost <- function(x, label) {
  final <- x[100, ]
  data.frame(
    cost = label,
    persisted = sum(final > 0),
    total_runs = length(final),
    persistence_percent = 100 * mean(final > 0),
    fixed = sum(final == 1),
    mean_final_frequency = mean(final)
  )
}

comparison <- rbind(
  summarise_cost(cost1, "1%"),
  summarise_cost(cost5, "5%")
)

print(comparison, row.names = FALSE)

# Compare mean trajectories, including runs that lost resistance.
plot(
  1:100, rowMeans(cost1),
  type = "l", col = "#0072B2", lwd = 3,
  ylim = c(0, 1),
  xlab = "Generation", ylab = "Mean resistance frequency",
  main = "Resistance cost: 20 runs per scenario"
)

lines(1:100, rowMeans(cost5), col = "#D55E00", lwd = 3)
abline(v = 50, lty = 2, col = "grey50")

legend(
  "topleft",
  legend = c("Cost: 1%", "Cost: 5%"),
  col = c("#0072B2", "#D55E00"),
  lty = 1, lwd = 3, bty = "n"
)

cost0 <- read_runs("~/SLiM/replicates_cost_000")

stopifnot(all(cost0[1:50, ] == cost1[1:50, ]))

comparison <- rbind(
  summarise_cost(cost0, "0%"),
  summarise_cost(cost1, "1%"),
  summarise_cost(cost5, "5%")
)

print(comparison, row.names = FALSE)

means <- cbind(
  rowMeans(cost0),
  rowMeans(cost1),
  rowMeans(cost5)
)

colours <- c("#009E73", "#0072B2", "#D55E00")

matplot(
  1:100, means,
  type = "l", lty = 1, lwd = 3,
  col = colours, ylim = c(0, 1),
  xlab = "Generation",
  ylab = "Mean resistance frequency",
  main = "Resistance after azole withdrawal",
  sub = "20 replicates per scenario; illustrative selection coefficients"
)

abline(v = 50, lty = 2, col = "grey50")

legend(
  "topleft",
  legend = c("Cost: 0%", "Cost: 1%", "Cost: 5%"),
  col = colours, lty = 1, lwd = 3, bty = "n"
)

long <- read.csv(
  path.expand("~/SLiM/long_cost_001_seed_1002.csv")
)

# Confirm the original trajectory is reproduced.
stopifnot(long$frequency[long$generation == 100] == 0.936)

# Find when fixation first occurred.
fixed_generations <- long$generation[long$frequency == 1]
if (length(fixed_generations) > 0) {
  cat("First fixation at generation:", min(fixed_generations), "\n")
}

plot(
  long$generation, long$frequency,
  type = "l", col = "#0072B2", lwd = 2,
  ylim = c(0, 1),
  xlab = "Generation", ylab = "Resistance frequency",
  main = "Long-term persistence: 1% cost, seed 1002"
)

abline(v = 50, lty = 2, col = "grey50")

read_long_runs <- function(cost) {
  folder <- path.expand(
    paste0("~/SLiM/long_replicates_cost_", cost)
  )
  files <- file.path(folder, paste0("seed_", 1001:1020, ".csv"))
  runs <- lapply(files, read.csv)

  stopifnot(all(vapply(runs, function(d) {
    identical(as.integer(d$generation), 1:500) &&
      all(is.finite(d$frequency)) &&
      all(d$frequency >= 0 & d$frequency <= 1)
  }, logical(1))))

  sapply(runs, function(d) d$frequency)
}

long0 <- read_long_runs("0.0")
long1 <- read_long_runs("0.01")
long5 <- read_long_runs("0.05")

# Extended runs must reproduce the original first 100 generations.
stopifnot(
  all(long0[1:100, ] == cost0),
  all(long1[1:100, ] == cost1),
  all(long5[1:100, ] == cost5)
)

summarise_long <- function(x, label) {
  final <- x[500, ]

  data.frame(
    cost = label,
    lost = sum(final == 0),
    segregating = sum(final > 0 & final < 1),
    fixed = sum(final == 1),
    persistence_percent = 100 * mean(final > 0),
    mean_final_frequency = mean(final)
  )
}

comparison500 <- rbind(
  summarise_long(long0, "0%"),
  summarise_long(long1, "1%"),
  summarise_long(long5, "5%")
)

print(comparison500, row.names = FALSE)


matplot(
  1:500,
  cbind(rowMeans(long0), rowMeans(long1), rowMeans(long5)),
  type = "l", lty = 1, lwd = 3,
  col = c("#009E73", "#0072B2", "#D55E00"),
  ylim = c(0, 1),
  xlab = "Generation",
  ylab = "Mean resistance frequency",
  main = "Resistance persistence over 500 generations",
  sub = "20 replicates per scenario"
)

abline(v = 50, lty = 2, col = "grey50")

legend(
  "topright",
  legend = c("Cost: 0%", "Cost: 1%", "Cost: 5%"),
  col = c("#009E73", "#0072B2", "#D55E00"),
  lty = 1, lwd = 3, bty = "n"
)

write.csv(
  comparison,
  path.expand("~/SLiM/summary_generation_100.csv"),
  row.names = FALSE
)

write.csv(
  comparison500,
  path.expand("~/SLiM/summary_generation_500.csv"),
  row.names = FALSE
)


read_exposure <- function(exposure) {
  folder <- path.expand(
    paste0("~/SLiM/exposure_duration_", exposure)
  )
  files <- file.path(folder, paste0("seed_", 1001:1020, ".csv"))
  runs <- lapply(files, read.csv)
  end <- exposure + 450L

  stopifnot(all(vapply(runs, function(d) {
    identical(as.integer(d$generation), seq_len(end)) &&
      all(is.finite(d$frequency)) &&
      all(d$frequency >= 0 & d$frequency <= 1)
  }, logical(1))))

  sapply(runs, function(d) d$frequency)
}

exposure25 <- read_exposure(25)
exposure50 <- read_exposure(50)
exposure100 <- read_exposure(100)

# All scenarios should match before the earliest withdrawal.
stopifnot(
  all(exposure25[1:25, ] == exposure50[1:25, ]),
  all(exposure25[1:25, ] == exposure100[1:25, ])
)

summarise_exposure <- function(x, exposure) {
  final <- x[nrow(x), ]

  data.frame(
    withdrawal_generation = exposure,
    mean_frequency_at_withdrawal = mean(x[exposure, ]),
    lost = sum(final == 0),
    segregating = sum(final > 0 & final < 1),
    fixed = sum(final == 1),
    persistence_percent = 100 * mean(final > 0),
    mean_final_frequency = mean(final)
  )
}

duration_comparison <- rbind(
  summarise_exposure(exposure25, 25),
  summarise_exposure(exposure50, 50),
  summarise_exposure(exposure100, 100)
)

print(duration_comparison, row.names = FALSE)

# Align trajectories at withdrawal:
# cycle 0 = withdrawal; cycle 450 = end of follow-up.
post_means <- cbind(
  rowMeans(exposure25[25:475, ]),
  rowMeans(exposure50[50:500, ]),
  rowMeans(exposure100[100:550, ])
)

matplot(
  0:450, post_means,
  type = "l", lty = 1, lwd = 3,
  col = c("#009E73", "#0072B2", "#D55E00"),
  ylim = c(0, 1),
  xlab = "Cycles after azole withdrawal",
  ylab = "Mean resistance frequency",
  main = "Exposure duration and resistance persistence",
  sub = "20 replicates per scenario; no cost after withdrawal"
)

legend(
  "topright",
  legend = c(
    "Withdrawal at generation 25",
    "Withdrawal at generation 50",
    "Withdrawal at generation 100"
  ),
  col = c("#009E73", "#0072B2", "#D55E00"),
  lty = 1, lwd = 3, bty = "n"
)

old_par <- par(no.readonly = TRUE)
par(mar = c(7, 4, 4, 2) + 0.1)

matplot(
  0:450, post_means,
  type = "l", lty = 1, lwd = 3,
  col = c("#009E73", "#0072B2", "#D55E00"),
  ylim = c(0, 1),
  xlab = "Cycles after azole withdrawal",
  ylab = "Mean resistance frequency",
  main = "Exposure duration and resistance persistence",
  sub = "20 replicates per scenario; no cost after withdrawal"
)

legend(
  "bottom",
  inset = c(0, -0.30),
  legend = c(
    "Withdrawal: generation 25",
    "Withdrawal: generation 50",
    "Withdrawal: generation 100"
  ),
  col = c("#009E73", "#0072B2", "#D55E00"),
  lty = 1, lwd = 3,
  horiz = TRUE,
  bty = "n",
  cex = 0.8,
  xpd = NA
)

par(old_par)

write.csv(
  duration_comparison,
  path.expand("~/SLiM/summary_exposure_duration.csv"),
  row.names = FALSE
)
