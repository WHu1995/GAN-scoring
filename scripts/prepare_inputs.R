#!/usr/bin/env Rscript

# Prepare natural E4.5 and blastoid cells on one comparable input scale.
#
# Every eligible cell is multinomially downsampled to the same total count,
# normalized to 10,000 counts, log1p transformed, and written in the same
# 2,000-gene order.  Repeated downsampling gives a simple measure of score
# stability.  Cells below target_depth are retained but marked ineligible for
# cross-study aggregate scores.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[[1]] else "data"
output_dir <- if (length(args) >= 2) args[[2]] else "prepared"
target_depth <- if (length(args) >= 3) as.integer(args[[3]]) else 1000L
natural_repeats <- if (length(args) >= 4) as.integer(args[[4]]) else 4L
blastoid_repeats <- if (length(args) >= 5) as.integer(args[[5]]) else 2L
seed <- if (length(args) >= 6) as.integer(args[[6]]) else 20260930L

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

reference <- readRDS(file.path(data_dir, "E4.5.rds"))

# Preserve the original model feature set, but undo names that R's
# make.names() or spreadsheet date conversion changed in the old CSV.
csv_header <- names(read.csv(
  file.path(data_dir, "E450_expression.csv"), nrows = 1L, check.names = FALSE
))
model_gene_names <- csv_header[!(csv_header %in% c("", "cell_type"))]
reference_names <- rownames(reference)
reference_rows <- match(model_gene_names, reference_names)
missing <- is.na(reference_rows)
reference_rows[missing] <- match(
  model_gene_names[missing], make.names(reference_names, unique = FALSE)
)
reference_rows[model_gene_names == "1-Sep"] <- match("Sept1", reference_names)
reference_rows[model_gene_names == "9-Sep"] <- match("Sept9", reference_names)
if (length(model_gene_names) != 2000L || anyNA(reference_rows)) {
  stop(
    "Could not recover all 2,000 original model genes from E4.5.rds; missing: ",
    paste(model_gene_names[is.na(reference_rows)], collapse = ", ")
  )
}
genes <- reference_names[reference_rows]
if (anyDuplicated(genes)) {
  stop("Recovered canonical gene names are not unique")
}
writeLines(genes, file.path(output_dir, "genes.txt"), useBytes = TRUE)
write.csv(
  data.frame(model_gene_name = model_gene_names, canonical_gene = genes),
  file.path(output_dir, "gene_name_map.csv"),
  row.names = FALSE,
  quote = TRUE
)

natural_condition <- function(cell_type) {
  result <- as.character(cell_type)
  result[result == "PE"] <- "Pre_VE"
  result
}

to_lineage <- function(cell_type) {
  value <- as.character(cell_type)
  out <- rep(NA_character_, length(value))
  out[value %in% c("Pre_EPI", "EPI-like")] <- "EPI"
  out[value %in% c("PaE", "PE", "Pre_VE", "PE-like")] <- "PE"
  out[value %in% c("mTE", "pTE", "TE-like")] <- "TE"
  out
}

write_fixed_depth <- function(object, dataset, repeats, is_reference = FALSE) {
  counts <- GetAssayData(object, assay = "RNA", layer = "counts")
  if (!inherits(counts, "sparseMatrix")) {
    counts <- as(counts, "dgCMatrix")
  }

  gene_rows <- match(genes, rownames(counts))
  source_to_feature <- integer(nrow(counts))
  available <- !is.na(gene_rows)
  source_to_feature[gene_rows[available]] <- which(available)

  n_cells <- ncol(counts)
  n_samples <- n_cells * repeats
  output <- matrix(0, nrow = n_samples, ncol = length(genes))
  cell_ids <- rep(colnames(counts), each = repeats)
  replicate_id <- rep(seq_len(repeats), times = n_cells)
  original_depth <- rep(as.numeric(Matrix::colSums(counts)), each = repeats)
  eligible <- original_depth >= target_depth

  row_out <- 0L
  for (j in seq_len(n_cells)) {
    column <- counts[, j, drop = FALSE]
    source_rows <- column@i + 1L
    source_values <- column@x
    total <- sum(source_values)
    destinations <- source_to_feature[source_rows]
    keep <- destinations > 0L

    for (r in seq_len(repeats)) {
      row_out <- row_out + 1L
      if (total >= target_depth) {
        sampled <- as.numeric(rmultinom(1L, target_depth, prob = source_values))
        denominator <- target_depth
      } else {
        sampled <- source_values
        denominator <- total
      }
      positive <- keep & sampled > 0
      if (denominator > 0 && any(positive)) {
        output[row_out, destinations[positive]] <- log1p(
          sampled[positive] / denominator * 10000
        )
      }
    }
  }

  raw_type <- object[["cell_type"]][, 1]
  if (is_reference) {
    condition <- natural_condition(raw_type)
  } else {
    condition <- rep(NA_character_, length(raw_type))
  }
  meta <- data.frame(
    dataset = dataset,
    cell_id = cell_ids,
    replicate = replicate_id,
    cell_type_raw = rep(as.character(raw_type), each = repeats),
    condition = rep(condition, each = repeats),
    lineage = rep(to_lineage(raw_type), each = repeats),
    original_depth = original_depth,
    eligible_depth = eligible,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  binary_path <- file.path(output_dir, paste0(dataset, ".f32"))
  connection <- file(binary_path, open = "wb")
  on.exit(close(connection), add = TRUE)
  # R matrices are column-major.  Transpose before vectorizing so Python can
  # reshape the file directly as (sample, gene) in C order.
  writeBin(as.vector(t(output)), connection, size = 4L, endian = "little")
  close(connection)
  on.exit(NULL, add = FALSE)
  write.csv(
    meta,
    file.path(output_dir, paste0(dataset, ".meta.csv")),
    row.names = FALSE,
    quote = TRUE
  )
  cat(
    dataset,
    ": cells=", n_cells,
    ", repeats=", repeats,
    ", mapped_genes=", sum(available),
    ", eligible=", sum(meta$eligible_depth), "/", nrow(meta),
    "\n",
    sep = ""
  )
}

write_fixed_depth(reference, "E4_5_reference", natural_repeats, is_reference = TRUE)
rm(reference)
gc()

blastoid_files <- setdiff(
  list.files(data_dir, pattern = "\\.rds$", full.names = TRUE),
  file.path(data_dir, "E4.5.rds")
)
for (path in sort(blastoid_files)) {
  object <- readRDS(path)
  dataset <- sub("\\.rds$", "", basename(path))
  write_fixed_depth(object, dataset, blastoid_repeats, is_reference = FALSE)
  rm(object)
  gc()
}

manifest <- data.frame(
  key = c(
    "target_depth", "scale_factor", "genes", "natural_repeats",
    "blastoid_repeats", "seed", "reference_file"
  ),
  value = c(
    target_depth, 10000, length(genes), natural_repeats,
    blastoid_repeats, seed, "E4.5.rds"
  ),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(output_dir, "manifest.csv"), row.names = FALSE, quote = TRUE)

cat("Prepared inputs written to ", normalizePath(output_dir), "\n", sep = "")
