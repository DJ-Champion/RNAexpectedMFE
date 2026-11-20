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
      plot.subtitle = element_text(colour = ink_muted),
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
#
# Both figures use the same encoding: observed means in colour, the model as a
# black dashed line. Each has a linear and a log version. MFE is negative, so
# the log version plots it on a log scale of its magnitude while keeping the
# sign and orientation: -1, -10, -100, -1000 run down the axis. Zero MFEs
# (short, GC-poor sequences, and the capped model there) cannot be shown on
# that scale and are dropped.

key_figure_size <- 16
key_figure_gc <- seq(0, 100, 10)

neg_log10_trans <- scales::new_transform(
  "neg_log10",
  transform = function(x) -log10(-x),
  inverse   = function(x) -10^(-x),
  domain    = c(-Inf, 0)
)

scale_y_mfe <- function(log = FALSE) {
  labels <- scales::label_number(big.mark = ",")
  if (log) {
    scale_y_continuous(transform = neg_log10_trans, labels = labels,
                       breaks = -c(1, 3, 10, 30, 100, 300, 1000, 3000))
  } else {
    scale_y_continuous(labels = labels)
  }
}

model_key <- function(model) {
  scale_linetype_manual(values = set_names("22", paste("Model", model)),
                        name = NULL)
}

plot_mfe_vs_length <- function(cells, model = "C", gc_levels = key_figure_gc,
                               log = FALSE) {
  d <- filter(cells, target_GC %in% gc_levels) |>
    mutate(gc = factor(target_GC, levels = gc_levels))
  if (log) d <- filter(d, MFE < 0)

  p <- ggplot(d, aes(length, colour = gc, group = gc)) +
    geom_line(aes(y = MFE), linewidth = 0.7) +
    geom_line(data = if (log) filter(d, expected < 0) else d,
              aes(y = expected, linetype = paste("Model", model)),
              colour = ink, linewidth = 0.45) +
    scale_colour_manual(values = colorRampPalette(blue_ramp)(length(gc_levels)),
                        name = "GC (%)") +
    model_key(model) +
    scale_y_mfe(log) +
    guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
    labs(title = paste0("Mean MFE vs sequence length, by GC content",
                        if (log) " (log scale)"),
         x = "Sequence length (nt)", y = "Mean MFE (kcal/mol)") +
    theme_rnaexpected(key_figure_size)

  if (log) {
    p + scale_x_log10(breaks = c(20, 50, 100, 200, 500, 1000, 2000, 4000))
  } else {
    p + scale_x_continuous(breaks = seq(0, 4000, 1000), expand = expansion(0.01))
  }
}

plot_mfe_vs_gc <- function(cells, model = "C", lengths = selected_lengths,
                           log = FALSE) {
  d <- filter(cells, length %in% lengths) |>
    mutate(length = factor(length, levels = lengths))
  if (log) d <- filter(d, MFE < 0)

  ggplot(d, aes(GC_content * 100, colour = length, group = length)) +
    geom_point(aes(y = MFE), size = 2) +
    geom_line(data = if (log) filter(d, expected < 0) else d,
              aes(y = expected, linetype = paste("Model", model)),
              colour = ink, linewidth = 0.45) +
    scale_colour_length() +
    model_key(model) +
    scale_y_mfe(log) +
    scale_x_continuous(breaks = seq(0, 100, 20), expand = expansion(0.01)) +
    guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
    labs(title = paste0("Mean MFE vs GC content, by sequence length",
                        if (log) " (log scale)"),
         x = "Measured GC content (%)", y = "Mean MFE (kcal/mol)") +
    theme_rnaexpected(key_figure_size)
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

# Relative error of each model on the synthetic data along the length axis,
# above the length distributions of real human transcript regions, so it reads
# as "how far off is the expected MFE at lengths like my regions?". Titles make
# clear the error is measured on synthetic sequences, not on the transcripts.
plot_error_vs_length <- function(rel, regions) {
  x_limits <- c(20, 40000)
  x_breaks <- c(20, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 30000)
  beyond <- annotate("rect", xmin = fitted_range$length[2], xmax = Inf,
                     ymin = -Inf, ymax = Inf, fill = neutral_mid)
  scale_x <- scale_x_log10(limits = x_limits, breaks = x_breaks,
                           labels = scales::label_number(big.mark = ","),
                           expand = c(0, 0))

  error <- ggplot(rel, aes(length, rel_error, colour = model)) +
    beyond +
    annotate("text", x = 12000, y = Inf, vjust = 1.5, size = 4.5,
             label = "Beyond training\ndata (> 4,000 nt)", colour = ink_muted) +
    geom_hline(yintercept = 0, colour = ink_muted, linewidth = 0.3) +
    geom_line(linewidth = 0.7) +
    scale_colour_manual(values = model_colours, name = "Model") +
    scale_x +
    # Zoomed to the useful range; below ~40 nt errors run off towards -100%.
    scale_y_continuous(breaks = seq(-50, 10, 10)) +
    coord_cartesian(ylim = c(-50, 10)) +
    labs(title = "Expected-MFE model error on synthetic random sequences",
         subtitle = "Median relative error of replicate means vs sequence length",
         x = NULL,
         y = paste0("Relative error (%)\nGC ", rel_error_gc[1], "-",
                    rel_error_gc[2], "%")) +
    theme_rnaexpected(key_figure_size) +
    theme(axis.text.x = element_blank())

  lengths <- regions |>
    filter(measure == "length") |>
    mutate(region = factor(region, levels = rev(unique(region))))
  ranges <- ggplot(lengths, aes(y = region)) +
    beyond +
    geom_linerange(aes(xmin = min, xmax = max), colour = ink_muted,
                   linewidth = 0.5) +
    geom_linerange(aes(xmin = q1, xmax = q3), colour = ink_muted,
                   linewidth = 5) +
    geom_point(aes(x = median), shape = 21, size = 3.5, fill = "white",
               colour = ink, stroke = 1) +
    scale_x +
    labs(subtitle = "For reference: lengths of human transcript regions (~13,000 transcripts)",
         x = "Sequence length (nt)", y = NULL) +
    theme_rnaexpected(key_figure_size) +
    theme(panel.grid.major.y = element_blank())

  error / ranges + plot_layout(heights = c(3, 1.3))
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
