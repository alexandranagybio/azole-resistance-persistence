# Run from the project root:
# Rscript analysis/make_figures.R

source("analysis/analyse_results.R")
figure_dir <- if (n_runs == 20) "figures" else
  file.path("figures", paste0("replicates_", n_runs))
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

colours <- c("#009E73", "#0072B2", "#D55E00")

plot_cost <- function() {
  par(mar = c(5, 4.5, 4, 1))
  layout(matrix(c(1, 2), ncol = 1), heights = c(5, 1))

  matplot(
    1:500,
    do.call(cbind, lapply(long, rowMeans)),
    type = "l", lty = 1, lwd = 3,
    col = colours, ylim = c(0, 1),
    xlab = "Model generation",
    ylab = "Mean resistance frequency",
    main = "Resistance persistence after azole withdrawal",
    sub = paste(n_runs, "runs per scenario; illustrative selection coefficients")
  )

  abline(v = 50, lty = 2, col = "grey40")

  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend(
    "center",
    legend = c("Cost: 0%", "Cost: 1%", "Cost: 5%", "Withdrawal: generation 50"),
    col = c(colours, "grey40"),
    lty = c(1, 1, 1, 2),
    lwd = c(3, 3, 3, 1),
    ncol = 2, bty = "n", cex = 0.9
  )
}

plot_exposure <- function() {
  post_means <- do.call(cbind, lapply(names(exposure), function(label) {
    withdrawal <- as.integer(label)
    x <- exposure[[label]]
    rowMeans(x[withdrawal:nrow(x), , drop = FALSE])
  }))

  par(mar = c(5, 4.5, 4, 1))
  layout(matrix(c(1, 2), ncol = 1), heights = c(5, 1))

  matplot(
    0:450, post_means,
    type = "l", lty = 1, lwd = 3,
    col = colours, ylim = c(0, 1),
    xlab = "Cycles after azole withdrawal",
    ylab = "Mean resistance frequency",
    main = "Exposure duration and resistance persistence",
    sub = paste(n_runs, "runs per scenario; zero cost after withdrawal")
  )

  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend(
    "center",
    legend = paste("Withdrawal: generation", names(exposure)),
    col = colours, lty = 1, lwd = 3,
    ncol = 1, bty = "n", cex = 0.9
  )
}

export_figure <- function(name, draw) {
  png(file.path(figure_dir, paste0(name, ".png")),
      width = 2400, height = 1800, res = 300)
  tryCatch(draw(), finally = dev.off())

  pdf(file.path(figure_dir, paste0(name, ".pdf")),
      width = 8, height = 6)
  tryCatch(draw(), finally = dev.off())
}

export_figure("cost_comparison_500", plot_cost)
export_figure("exposure_duration_comparison", plot_exposure)

cat("\nSaved trajectory figures in", figure_dir, "\n")


# Points show endpoint persistence; bars show 95% Wilson intervals.
plot_persistence <- function() {
  par(mfrow = c(1, 2), mar = c(5, 4.5, 4, 1), oma = c(0, 0, 2, 0))

  draw_panel <- function(d, labels, title, x_label) {
    x <- seq_len(nrow(d))
    plot(
      x, d$persistence_percent,
      type = "n", xaxt = "n",
      xlim = c(0.5, length(x) + 0.5), ylim = c(0, 100),
      xlab = x_label, ylab = "Runs retaining resistance (%)",
      main = title, cex.main = 0.85
    )
    axis(1, at = x, labels = labels)
    abline(h = seq(0, 100, 20), col = "grey90")
    arrows(
      x, d$persistence_ci_lower_percent,
      x, d$persistence_ci_upper_percent,
      angle = 90, code = 3, length = 0.06,
      col = colours, lwd = 2
    )
    points(x, d$persistence_percent, pch = 19, cex = 1.4, col = colours)
  }

  draw_panel(
    summary500, summary500$cost,
    "Cost: generation 500",
    "Post-withdrawal cost"
  )

  draw_panel(
    summary_exposure, summary_exposure$withdrawal_generation,
    "Exposure: 450 cycles after withdrawal",
    "Withdrawal generation"
  )

  mtext(
    paste(n_runs, "runs per scenario; 95% Wilson intervals"),
    outer = TRUE, side = 3, line = 0.3
  )
}

export_figure("persistence_confidence_intervals", plot_persistence)
cat("Saved persistence figure in", figure_dir, "\n")
