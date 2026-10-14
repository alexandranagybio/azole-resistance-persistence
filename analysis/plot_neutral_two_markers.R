# Run from the project root:
# Rscript analysis/plot_neutral_two_markers.R

input_root <- "data/summaries/neutral_two_markers_100"
figure_root <- "figures/neutral_two_markers_100"

s <- read.csv(file.path(input_root, "summary.csv"))
t <- read.csv(file.path(input_root, "mean_trajectories.csv"))

s <- s[order(s$sexual_fraction), ]
fractions <- c(0, 0.01, 0.1, 1)
labels <- c("0%", "1%", "10%", "100%")

stopifnot(
  nrow(s) == 4,
  all(s$sexual_fraction == fractions),
  all(s$runs == 100),
  nrow(t) == 2000
)

# Recover mean genotype frequencies from the saved mean marker
# and double-carrier frequencies. These are linear identities.
t$tr34_only <- t$tr34_frequency - t$double_frequency
t$l98h_only <- t$l98h_frequency - t$double_frequency
t$wt <- 1 - t$tr34_frequency - t$l98h_frequency +
  t$double_frequency

genotype_columns <- c(
  "wt", "tr34_only", "l98h_only", "double_frequency"
)
g <- as.matrix(t[, genotype_columns])

stopifnot(
  all(is.finite(g)),
  all(g >= -1e-8 & g <= 1 + 1e-8),
  all(abs(rowSums(g) - 1) < 1e-8)
)

dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)

save_figure <- function(name, draw, width, height) {
  for (extension in c("png", "pdf")) {
    path <- file.path(figure_root, paste0(name, ".", extension))

    if (extension == "png") {
      png(path, width = width, height = height,
          units = "in", res = 180)
    } else {
      cairo_pdf(path, width = width, height = height,
                family = "sans")
    }

    tryCatch(draw(), finally = dev.off())
  }
}

draw_appearance <- function() {
  par(mar = c(7, 5, 4, 2))
  x <- seq_along(fractions)

  plot(
    x, s$singles_ever_percent,
    type = "n", xaxt = "n",
    xlim = c(0.6, 4.4), ylim = c(0, 110),
    xlab = "Hypothetical two-parent fraction",
    ylab = "Runs with single-marker appearance (%)",
    main = "Single-marker appearance by cycle 500",
    sub = "100 runs per scenario; conditional recombination fraction = 0.00192"
  )

  axis(1, at = x, labels = labels)
  abline(h = seq(0, 100, 20), col = "grey90")

  arrows(
    x, s$ci_lower_percent,
    x, s$ci_upper_percent,
    angle = 90, code = 3, length = 0.06,
    lwd = 1.5, col = "#0072B2"
  )

  points(x, s$singles_ever_percent,
         pch = 19, cex = 1.4, col = "#0072B2")

  text(
    x, s$ci_upper_percent + 4,
    labels = paste0(s$singles_ever, "/100"),
    cex = 0.9
  )

  mtext("Points: observed proportions; bars: 95% Wilson intervals",
        side = 3, line = 0.4, cex = 0.85)
}

draw_trajectories <- function() {
  colours <- c("#666666", "#0072B2", "#D55E00", "#009E73")
  line_types <- c(1, 2, 3, 1)

  layout(
    matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE),
    heights = c(1, 1, 0.22)
  )
  par(oma = c(0, 0, 3, 0))

  for (q in fractions) {
    d <- t[t$sexual_fraction == q, ]
    d <- d[order(d$generation), ]

    stopifnot(
      nrow(d) == 500,
      all(d$generation == seq_len(500))
    )

    par(mar = c(4, 4, 3, 1))

    matplot(
      d$generation,
      as.matrix(d[, genotype_columns]),
      type = "l", lty = line_types, lwd = 2,
      col = colours, ylim = c(0, 1),
      xlab = "Reproductive cycle",
      ylab = "Mean genotype frequency",
      main = paste0("Two-parent fraction: ", 100 * q, "%")
    )
  }

  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend(
    "center",
    legend = c("Wild type", "TR34 only", "L98H only", "Double carrier"),
    col = colours, lty = line_types, lwd = 2,
    horiz = TRUE, bty = "n", cex = 0.9
  )

  mtext(
    "Neutral inheritance: means across all 100 runs per scenario",
    outer = TRUE, side = 3, line = 1, cex = 1.1
  )
}

save_figure(
  "single_marker_appearance",
  draw_appearance, width = 8, height = 6
)

save_figure(
  "mean_genotype_trajectories",
  draw_trajectories, width = 10, height = 8
)

cat("Saved both figures as PNG and PDF in", figure_root, "\n")
