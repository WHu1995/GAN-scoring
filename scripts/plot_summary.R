#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[[1]] else file.path("results", "calibrated")

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

assign_group <- function(dataset) {
  ifelse(
    dataset == "Natural E4.5 holdout",
    "Natural E4.5",
    ifelse(
      dataset == "Jorgensen_iGata4_blastoids",
      "Jorgensen iGata4",
      "Other blastoid"
    )
  )
}

group_values <- c(
  "Natural E4.5" = "#0072B2",
  "Jorgensen iGata4" = "#009E73",
  "Other blastoid" = "#8A8A8A"
)

lineage_boot <- read.csv(
  file.path(input_dir, "lineage_bootstrap.csv"),
  stringsAsFactors = FALSE
)
lineage_summary <- read.csv(
  file.path(input_dir, "condition_and_lineage_distances.csv"),
  stringsAsFactors = FALSE
)
lineage_summary <- lineage_summary[lineage_summary$level == "lineage", ]

overall_boot <- read.csv(
  file.path(input_dir, "overall_bootstrap.csv"),
  stringsAsFactors = FALSE
)
overall_summary <- read.csv(
  file.path(input_dir, "overall_distances.csv"),
  stringsAsFactors = FALSE
)

lineage_boot$panel <- lineage_boot$lineage
lineage_summary$panel <- lineage_summary$lineage
overall_boot$panel <- "Overall"
overall_summary$panel <- "Overall"

bootstrap <- rbind(
  lineage_boot[, c("dataset", "panel", "iteration", "distance")],
  overall_boot[, c("dataset", "panel", "iteration", "distance")]
)
summary <- rbind(
  lineage_summary[, c("dataset", "panel", "distance", "ci_low", "ci_high")],
  overall_summary[, c("dataset", "panel", "distance", "ci_low", "ci_high")]
)

prepare_fields <- function(frame) {
  frame$dataset_label <- unname(dataset_labels[frame$dataset])
  frame$dataset_label <- factor(frame$dataset_label, levels = rev(dataset_order))
  frame$panel <- factor(frame$panel, levels = c("Overall", "EPI", "PE", "TE"))
  frame$group <- assign_group(frame$dataset)
  frame
}

bootstrap <- prepare_fields(bootstrap)
summary <- prepare_fields(summary)

plot <- ggplot(
  bootstrap,
  aes(x = distance, y = dataset_label, fill = group, color = group)
) +
  geom_violin(
    orientation = "y",
    scale = "width",
    trim = FALSE,
    width = 0.84,
    linewidth = 0.45,
    alpha = 0.62
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
    size = 2.3
  ) +
  facet_wrap(~panel, ncol = 2, scales = "free_x") +
  scale_fill_manual(values = group_values) +
  scale_color_manual(values = group_values) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.08))) +
  labs(
    title = "Original sccGAN critic, calibrated to natural E4.5",
    subtitle = paste0(
      "The original network and critic.pth are unchanged; raw critic values are converted to ",
      "condition-specific distribution distances"
    ),
    x = "Natural-calibrated critic-output W1 distance (lower is closer)",
    y = NULL,
    fill = NULL,
    color = NULL,
    caption = paste0(
      "Violin width shows 2,000 cell-level bootstrap replicates; white dots are point estimates and bars are 95% CIs. ",
      "All datasets use the same fixed-depth preprocessing and original 2,000-gene order.\n",
      "PE/TE scores use the closest corresponding natural stage. The natural split is a scoring-calibration holdout, not an external model-training holdout."
    )
  ) +
  theme_classic(base_size = 11.5) +
  theme(
    strip.background = element_rect(fill = "#EFEFEF", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 10.6, color = "#4D4D4D"),
    plot.caption = element_text(hjust = 0, size = 8.7, color = "#555555"),
    axis.text.y = element_text(color = "#222222"),
    legend.position = "bottom"
  )

png_path <- file.path(input_dir, "figure_original_critic_natural_calibrated_violin.png")
pdf_path <- file.path(input_dir, "figure_original_critic_natural_calibrated_violin.pdf")

ggsave(png_path, plot, width = 13.5, height = 10, dpi = 320, bg = "white")
ggsave(pdf_path, plot, width = 13.5, height = 10, device = cairo_pdf, bg = "white")

cat("Wrote:", png_path, "\n")
cat("Wrote:", pdf_path, "\n")
