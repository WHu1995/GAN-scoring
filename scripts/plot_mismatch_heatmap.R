#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[[1]] else file.path("results", "calibrated")
input_file <- file.path(input_dir, "blastoid_celltype_mismatch_matrix.csv")

condition_order <- c("Pre_EPI", "Pre_VE", "PaE", "mTE", "pTE")
lineage_order <- c("EPI", "PE", "TE")
dataset_order <- c(
  "Jorgensen_iGata4_blastoids",
  "Xu_TPS_blastoids",
  "Li_EPSC_blastoids",
  "Bing_TBLC_blastoids",
  "Sozen_EPSC_TSC",
  "Rivron_ESC_TSC_blastoids"
)
dataset_labels <- c(
  "Jorgensen_iGata4_blastoids" = "Jorgensen iGata4",
  "Xu_TPS_blastoids" = "Xu TPS",
  "Li_EPSC_blastoids" = "Li EPSC",
  "Bing_TBLC_blastoids" = "Bing TBLC",
  "Sozen_EPSC_TSC" = "Sozen EPSC + TSC",
  "Rivron_ESC_TSC_blastoids" = "Rivron ESC + TSC"
)
lineage_labels <- c(EPI = "EPI-like", PE = "PE-like", TE = "TE-like")

matrix <- read.csv(input_file, stringsAsFactors = FALSE)
matrix$dataset_label <- unname(dataset_labels[matrix$dataset])
matrix$dataset_label <- factor(
  matrix$dataset_label,
  levels = unname(dataset_labels[dataset_order])
)
matrix$annotated_lineage <- factor(
  unname(lineage_labels[matrix$annotated_lineage]),
  levels = rev(unname(lineage_labels[lineage_order]))
)
matrix$supplied_condition <- factor(
  matrix$supplied_condition,
  levels = condition_order
)
matrix$within_annotated_lineage <- tolower(
  as.character(matrix$within_annotated_lineage)
) == "true"
matrix$label <- sprintf("%.2f", matrix$distance)

fill_limit <- max(matrix$ci_high, matrix$distance, na.rm = TRUE)
text_switch <- 0.58 * fill_limit

figure <- ggplot(
  matrix,
  aes(x = supplied_condition, y = annotated_lineage, fill = distance)
) +
  geom_tile(color = "white", linewidth = 0.9) +
  geom_tile(
    data = matrix[matrix$within_annotated_lineage, ],
    fill = NA,
    color = "#111111",
    linewidth = 1.15
  ) +
  geom_text(
    aes(label = label, color = ifelse(distance > text_switch, "light", "dark")),
    fontface = "bold",
    size = 3.45
  ) +
  facet_wrap(~dataset_label, ncol = 3, drop = FALSE) +
  scale_fill_gradientn(
    colors = c("#2166AC", "#67A9CF", "#F7F7F7", "#F4A582", "#B2182B"),
    limits = c(0, fill_limit),
    name = "Calibrated W1\ndistance"
  ) +
  scale_color_manual(
    values = c(dark = "#111111", light = "white"),
    guide = "none"
  ) +
  labs(
    title = "Cell-type mismatch challenge using the original sccGAN critic",
    subtitle = paste0(
      "Each annotated blastoid lineage is scored against all five natural E4.5 conditions; ",
      "lower values indicate greater similarity"
    ),
    x = "Natural E4.5 condition supplied to critic",
    y = "Annotated blastoid cell type",
    caption = paste0(
      "Values are robustly standardized W1 distances between blastoid and natural-calibration critic-output distributions. ",
      "Black borders mark labels belonging to the annotated lineage; unbordered cells are mismatch tests.\n",
      "A well condition-specific critic should show lower distances inside black borders and higher distances for cross-lineage labels. ",
      "Rivron has no annotated PE-like population."
    )
  ) +
  coord_fixed(ratio = 0.82) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "#EFEFEF", color = NA),
    strip.text = element_text(face = "bold", size = 10.5),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 10, color = "#4D4D4D"),
    plot.caption = element_text(size = 8.4, color = "#555555", hjust = 0),
    axis.text.x = element_text(face = "bold", angle = 35, hjust = 1),
    axis.text.y = element_text(color = "#222222"),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.spacing = grid::unit(0.85, "lines"),
    legend.position = "right"
  )

png_path <- file.path(input_dir, "figure_original_critic_celltype_mismatch_heatmap.png")
pdf_path <- file.path(input_dir, "figure_original_critic_celltype_mismatch_heatmap.pdf")

ggsave(png_path, figure, width = 14, height = 8.7, dpi = 320, bg = "white")
ggsave(pdf_path, figure, width = 14, height = 8.7, device = cairo_pdf, bg = "white")

cat("Wrote:", png_path, "\n")
cat("Wrote:", pdf_path, "\n")
