# Run from the project root:
# Rscript analysis/plot_two_populations.R

input_root <- "data/two_population_dispersal_100"
figure_root <- "figures/two_population_dispersal_100"
rates <- c("0.0", "0.001", "0.01")
rate_labels <- c("0%", "0.1%", "1%")
columns <- c("frequency_p1", "frequency_p2", "frequency_total")
colours <- c("#0072B2", "#D55E00", "#222222")

dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)

runs <- lapply(rates, function(rate) {
  lapply(1001:1100, function(seed) {
    read.csv(file.path(
      input_root, paste0("migration_", rate),
      paste0("seed_", seed, ".csv")
    ))
  })
})

# Rows are cycles; columns are replicate runs.
matrices <- lapply(runs, function(scenario) {
  setNames(lapply(columns, function(column) {
    vapply(scenario, function(d) d[[column]], numeric(500))
  }), columns)
})

save_figure <- function(name, draw) {
  png(file.path(figure_root, paste0(name, ".png")),
      width = 3300, height = 1800, res = 300)
  draw()
  dev.off()

  pdf(file.path(figure_root, paste0(name, ".pdf")),
      width = 11, height = 6)
  draw()
  dev.off()
}

draw_trajectories <- function() {
  layout(matrix(c(1, 2, 3, 4, 4, 4), nrow = 2, byrow = TRUE),
         heights = c(5, 1))
  par(mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

  for (i in seq_along(rates)) {
    m <- matrices[[i]]

    plot(NA, xlim = c(1, 500), ylim = c(0, 1),
         xlab = "Reproductive cycle",
         ylab = "Resistance frequency",
         main = paste("Dispersal:", rate_labels[i]))

    # Grey lines show overall frequency in each individual run.
    matlines(1:500, m$frequency_total,
             col = adjustcolor("grey40", alpha.f = 0.18),
             lty = 1, lwd = 0.6)

    abline(v = 50, lty = 2, col = "grey40")

    for (j in seq_along(columns)) {
      lines(1:500, rowMeans(m[[columns[j]]]),
            col = colours[j], lwd = 2)
    }
  }

  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center",
         legend = c("Mean p1", "Mean p2", "Mean overall",
                    "Individual runs: overall", "Withdrawal"),
         col = c(colours, "grey65", "grey40"),
         lty = c(1, 1, 1, 1, 2),
         lwd = c(2, 2, 2, 1, 1),
         horiz = TRUE, bty = "n", cex = 0.85)

  mtext("Two environmental patches | zero cost | 100 runs per scenario",
        side = 3, outer = TRUE, cex = 1)
}

draw_withdrawal <- function() {
  par(mfrow = c(1, 3), mar = c(4, 4, 3, 1),
      oma = c(2, 0, 3, 0))
  set.seed(42)

  titles <- c("Exposed patch: p1", "Unexposed patch: p2", "Overall")

  for (j in seq_along(columns)) {
    values <- lapply(matrices, function(m) m[[columns[j]]][50, ])

    boxplot(values, names = rate_labels, ylim = c(0, 1),
            outline = FALSE, col = "grey95",
            xlab = "Dispersal per cycle",
            ylab = "Resistance frequency",
            main = titles[j], cex.main = 0.95)

    for (i in seq_along(values)) {
      points(jitter(rep(i, length(values[[i]])), amount = 0.12),
             values[[i]], pch = 16, cex = 0.5,
             col = adjustcolor(colours[j], alpha.f = 0.45))
      points(i, mean(values[[i]]), pch = 23,
             bg = "white", col = "black", cex = 1.2, lwd = 1.4)
    }
  }

  mtext("Resistance at withdrawal: cycle 50",
        side = 3, outer = TRUE, line = 1, cex = 1.1)
  mtext("Dots: individual runs | White diamonds: means | Box centre: median",
        side = 1, outer = TRUE, cex = 0.85)
}

save_figure("dispersal_trajectories", draw_trajectories)
save_figure("frequencies_at_withdrawal", draw_withdrawal)

cat("Saved both figures as PNG and PDF in", figure_root, "\n")
