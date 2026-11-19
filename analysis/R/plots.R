# Plot helpers. Colour roles:
#   - GC content and length are magnitudes -> one blue ramp, light to dark
#     (starting at step 250 so the lightest line still reads on white).
#   - Models A/B/C are identities -> fixed categorical slots, in that order.
#   - Residual sign is polarity -> blue/red diverging around a neutral grey.

blue_ramp <- c("#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
               "#2a78d6", "#1c5cab", "#104281", "#0d366b")
model_colours <- c(A = "#2a78d6", B = "#eb6834", C = "#1baf7a")
ink <- "#0b0b0b"
ink_muted <- "#52514e"
neutral_mid <- "#f0efec"

# Lengths shown when a plot needs a handful of representative lines.
selected_lengths <- c(50, 100, 250, 500, 1000, 2000, 3000, 4000)

theme_rnaexpected <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#e6e5e1", linewidth = 0.3),
      axis.title = element_text(colour = ink_muted),
      axis.text = element_text(colour = ink_muted),
      plot.title = element_text(colour = ink, face = "bold"),
      plot.caption = element_text(colour = ink_muted, hjust = 0),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      legend.title = element_text(colour = ink_muted)
    )
}

scale_colour_gc <- function(...) {
  scale_colour_gradientn(colours = blue_ramp, limits = c(0, 100),
                         name = "GC (%)", ...)
}

scale_colour_length <- function(...) {
  scale_colour_manual(values = set_names(blue_ramp, selected_lengths),
                      name = "Length (nt)", ...)
}

# Mean observed MFE per design cell (length x target GC, 10 replicates each),
# with the expected MFE at the cell's mean measured GC.
cell_means <- function(data, model = "C") {
  data |>
    group_by(length, target_GC) |>
    summarise(GC_content = mean(GC_content), MFE = mean(MFE), .groups = "drop") |>
    mutate(expected = expected_mfe(GC_content, length, model = model),
           residual = MFE - expected)
}

save_figure <- function(plot, name, width = 8, height = 5,
                        dir = "results/figures", formats = c("pdf", "png")) {
  for (ext in formats) {
    ggsave(file.path(dir, paste0(name, ".", ext)), plot,
           width = width, height = height, dpi = 200, bg = "white")
  }
}

# --- Key figures --------------------------------------------------------------

plot_mfe_vs_length <- function(cells, model = "C", overlay_gc = seq(0, 100, 10)) {
  ggplot(cells, aes(length, MFE, colour = target_GC, group = target_GC)) +
    geom_line(linewidth = 0.45) +
    geom_line(data = filter(cells, target_GC %in% overlay_gc),
              aes(y = expected), colour = ink, linewidth = 0.3, linetype = "22") +
    scale_colour_gc() +
    scale_x_continuous(breaks = seq(0, 4000, 500), expand = expansion(0.01)) +
    labs(title = "Mean MFE vs sequence length, by GC content",
         x = "Sequence length (nt)", y = "Mean MFE (kcal/mol)",
         caption = paste0("Coloured: mean of 10 RNAfold replicates per length and target GC (2% steps). ",
                          "Dashed: model ", model, " at GC 0, 10, ..., 100%.")) +
    theme_rnaexpected()
}

plot_mfe_vs_gc <- function(cells, model = "C", lengths = selected_lengths) {
  cells |>
    filter(length %in% lengths) |>
    mutate(length = factor(length, levels = lengths)) |>
    ggplot(aes(GC_content * 100, colour = length)) +
    geom_line(aes(y = expected), linewidth = 0.5) +
    geom_point(aes(y = MFE), size = 1.4) +
    scale_colour_length() +
    scale_x_continuous(breaks = seq(0, 100, 10), expand = expansion(0.01)) +
    labs(title = "Mean MFE vs GC content, by sequence length",
         x = "Measured GC content (%)", y = "Mean MFE (kcal/mol)",
         caption = paste0("Points: mean of 10 RNAfold replicates per target GC. ",
                          "Lines: model ", model, ".")) +
    theme_rnaexpected()
}

# --- Diagnostics ----------------------------------------------------------------

# Five-panel residual diagnostics for one model, on design-cell means so that
# replicate noise is averaged out and systematic misfit is visible. `limits`
# fixes the residual axis so models can be compared panel for panel.
plot_diagnostics <- function(cells, model, limits) {
  resid_lab <- "Residual (kcal/mol)"
  zero <- geom_hline(yintercept = 0, colour = ink_muted, linewidth = 0.3)

  heat <- ggplot(cells, aes(target_GC, length, fill = residual)) +
    geom_raster() +
    scale_fill_gradient2(low = "#1c5cab", mid = neutral_mid, high = "#c23b3a",
                         midpoint = 0, limits = limits, oob = scales::squish,
                         name = "Observed -\nexpected") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = "Residual by GC and length", x = "Target GC (%)",
         y = "Length (nt)") +
    theme_rnaexpected()

  by_gc <- cells |>
    filter(length %in% selected_lengths) |>
    mutate(length = factor(length, levels = selected_lengths)) |>
    ggplot(aes(GC_content * 100, residual, colour = length)) +
    zero + geom_line(linewidth = 0.4) +
    scale_colour_length() +
    coord_cartesian(ylim = limits) +
    labs(title = "Residual vs GC", x = "Measured GC (%)", y = resid_lab) +
    theme_rnaexpected()

  by_length <- ggplot(cells, aes(length, residual, colour = target_GC,
                                 group = target_GC)) +
    zero + geom_line(linewidth = 0.3) +
    scale_colour_gc() +
    coord_cartesian(ylim = limits) +
    labs(title = "Residual vs length", x = "Length (nt)", y = resid_lab) +
    theme_rnaexpected()

  vs_fitted <- ggplot(cells, aes(expected, residual)) +
    zero + geom_point(size = 0.4, alpha = 0.25, colour = ink_muted) +
    coord_cartesian(ylim = limits) +
    labs(title = "Residual vs expected", x = "Expected MFE (kcal/mol)",
         y = resid_lab) +
    theme_rnaexpected()

  qq <- ggplot(cells, aes(sample = residual)) +
    stat_qq(size = 0.4, alpha = 0.4, colour = ink_muted) +
    stat_qq_line(colour = model_colours[[model]], linewidth = 0.5) +
    labs(title = "Normal Q-Q", x = "Theoretical quantile",
         y = "Sample quantile") +
    theme_rnaexpected()

  (heat | by_gc) / (by_length | vs_fitted | qq) +
    plot_annotation(
      title = paste0("Model ", model, ": ", mfe_models[[model]]$description),
      subtitle = deparse1(mfe_models[[model]]$formula),
      caption = "Residuals are observed minus expected MFE, averaged over the 10 replicates in each design cell.",
      theme = theme_rnaexpected()
    )
}

plot_rmse_by_length <- function(by_length) {
  long <- by_length |>
    pivot_longer(c(rmse, cell_mean_rmse), names_to = "metric") |>
    mutate(metric = recode(metric,
                           rmse = "Individual sequences",
                           cell_mean_rmse = "Design-cell means (systematic error)"),
           metric = factor(metric, levels = c("Individual sequences",
                                              "Design-cell means (systematic error)")),
           band = factor(band, levels = unique(by_length$band)))
  floor <- by_length |>
    distinct(band, replicate_sd) |>
    mutate(band = factor(band, levels = levels(long$band)),
           metric = "Individual sequences")

  ggplot(long, aes(band, value, colour = model, group = model)) +
    geom_line(data = floor, mapping = aes(band, replicate_sd, group = 1),
              colour = ink_muted, linetype = "22", inherit.aes = FALSE) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2.2) +
    facet_wrap(~metric) +
    scale_colour_manual(values = model_colours, name = "Model") +
    labs(title = "Prediction error by length band",
         x = "Length band (nt)", y = "RMSE (kcal/mol)",
         caption = "Dashed grey: pooled SD between the 10 replicates of each design cell (approximate noise floor).") +
    theme_rnaexpected() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}
