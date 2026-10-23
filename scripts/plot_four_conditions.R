#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[[1]] else file.path("results", "calibrated")
conditions <- c("Pre_VE", "PaE", "mTE", "pTE")

dataset_labels <- c(
  "Natural E4.5 holdout" = "Natural E4.5 holdout",
  "Jorgensen_iGata4_blastoids" = "Jorgensen iGata4",
  "Xu_TPS_blastoids" = "Xu TPS",
  "Li_EPSC_blastoids" = "Li EPSC",
  "Bing_TBLC_blastoids" = "Bing TBLC",
  "Sozen_EPSC_TSC" = "Sozen EPSC + TSC",
  "Rivron_ESC_TSC_blastoids" = "Rivron ESC + TSC"
)

dataset_order <- unname(dataset_labels)

dataset_colors <- c(
  "Natural E4.5 holdout" = "#0072B2",
  "Jorgensen iGata4" = "#009E73",
  "Xu TPS" = "#CC79A7",
  "Li EPSC" = "#D55E00",
  "Bing TBLC" = "#E69F00",
  "Sozen EPSC + TSC" = "#56B4E9",
  "Rivron ESC + TSC" = "#7A5195"
)

bootstrap <- read.csv(
  file.path(input_dir, "condition_bootstrap.csv"),
  stringsAsFactors = FALSE
)
summary <- read.csv(
  file.path(input_dir, "condition_and_lineage_distances.csv"),
  stringsAsFactors = FALSE
)

bootstrap <- bootstrap[bootstrap$condition %in% conditions, ]
summary <- summary[summary$level == "condition" & summary$condition %in% conditions, ]

prepare_fields <- function(frame) {
  frame$dataset_label <- unname(dataset_labels[frame$dataset])
  frame$dataset_label <- factor(frame$dataset_label, levels = rev(dataset_order))
  frame$condition <- factor(frame$condition, levels = conditions)
  frame
}

bootstrap <- prepare_fields(bootstrap)
summary <- prepare_fields(summary)

plot <- ggplot(
  bootstrap,
  aes(
    x = distance,
    y = dataset_label,
    fill = dataset_label,
    color = dataset_label
  )
) +
  geom_violin(
    orientation = "y",
    scale = "width",
    trim = FALSE,
    width = 0.82,
    linewidth = 0.5,
    alpha = 0.64
  ) +
  geom_errorbarh(
    data = summary,
    aes(y = dataset_label, xmin = ci_low, xmax = ci_high),
    inherit.aes = FALSE,
    height = 0,
    linewidth = 0.65,
    color = "#2B2B2B"
  ) +
  geom_point(
    data = summary,
    aes(x = distance, y = dataset_label),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    color = "#2B2B2B",
    stroke = 0.75,
    size = 2.4
  ) +
  facet_wrap(~condition, ncol = 2, scales = "free_x") +
  scale_fill_manual(values = dataset_colors, drop = FALSE) +
  scale_color_manual(values = dataset_colors, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.08))) +
  labs(
    title = "Original sccGAN critic scores for four natural E4.5 conditions",
    subtitle = paste0(
      "Original critic.pth with condition-specific natural calibration; ",
      "each blastoid is shown in a distinct color"
    ),
    x = "Natural-calibrated critic-output W1 distance (lower is closer)",
    y = NULL,
    fill = NULL,
    color = NULL,
    caption = paste0(
      "Violin width shows 2,000 cell-level bootstrap replicates; white dots are point estimates and bars are 95% CIs. ",
      "PE-like cells are evaluated separately against Pre_VE and PaE; TE-like cells separately against mTE and pTE.\n",
      "All datasets use identical fixed-depth preprocessing and the original 2,000-gene order."
    )
  ) +
  theme_classic(base_size = 11.5) +
  theme(
    strip.background = element_rect(fill = "#EFEFEF", color = NA),
    strip.text = element_text(face = "bold", size = 12.5),
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 10.6, color = "#4D4D4D"),
    plot.caption = element_text(hjust = 0, size = 8.7, color = "#555555"),
    axis.text.y = element_text(color = "#222222"),
    legend.position = "bottom",
    legend.key.width = grid::unit(0.7, "cm")
  ) +
  guides(
    fill = guide_legend(nrow = 2, byrow = TRUE),
    color = guide_legend(nrow = 2, byrow = TRUE)
  )

png_path <- file.path(input_dir, "figure_original_critic_four_conditions_violin.png")
pdf_path <- file.path(input_dir, "figure_original_critic_four_conditions_violin.pdf")

ggsave(png_path, plot, width = 13.5, height = 10, dpi = 320, bg = "white")
ggsave(pdf_path, plot, width = 13.5, height = 10, device = cairo_pdf, bg = "white")

cat("Wrote:", png_path, "\n")
cat("Wrote:", pdf_path, "\n")
