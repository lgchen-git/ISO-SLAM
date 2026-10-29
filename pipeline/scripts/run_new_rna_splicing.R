#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag) {
  pos <- match(flag, args)
  if (is.na(pos) || pos == length(args)) stop("Missing argument: ", flag)
  args[[pos + 1L]]
}
samples <- fread(value_after("--samples"))
psi <- fread(value_after("--psi"), na.strings = c("NA", ""))
contrasts <- fromJSON(value_after("--contrasts-json"), simplifyVector = FALSE)
output <- value_after("--out")
required <- c("sample_id", "group")
if (!all(required %in% names(samples))) stop("Sample sheet lacks sample_id/group")
if (!all(c("event_id", "gene_id", "sample_id", "inclusion_new_read_equivalents", "exclusion_new_read_equivalents", "new_RNA_PSI") %in% names(psi))) stop("PSI table schema is incomplete")
sample_ids <- samples$sample_id
psi_wide <- dcast(psi, event_id + gene_id ~ sample_id, value.var = "new_RNA_PSI")
support_wide <- dcast(psi, event_id + gene_id ~ sample_id, value.var = "total_new_read_equivalents")
rownames(psi_wide) <- psi_wide$event_id
rownames(support_wide) <- support_wide$event_id
psi_mat <- as.matrix(psi_wide[, ..sample_ids])
support_mat <- as.matrix(support_wide[, ..sample_ids])
rownames(psi_mat) <- psi_wide$event_id
rownames(support_mat) <- support_wide$event_id
storage.mode(psi_mat) <- "numeric"
storage.mode(support_mat) <- "numeric"
annotation <- psi_wide[, .(event_id, gene_id)]
results <- list()
for (item in contrasts) {
  g1 <- as.character(item$treatment)
  g2 <- as.character(item$control)
  ids1 <- samples[group == g1, sample_id]
  ids2 <- samples[group == g2, sample_id]
  if (length(ids1) < 2L || length(ids2) < 2L) stop("At least two replicates per group are required")
  eligible <- rowSums(support_mat[, ids1, drop = FALSE] >= 20, na.rm = TRUE) >= 2L &
    rowSums(support_mat[, ids2, drop = FALSE] >= 20, na.rm = TRUE) >= 2L
  if (!any(eligible)) next
  selected <- c(ids1, ids2)
  p <- psi_mat[eligible, selected, drop = FALSE]
  n <- support_mat[eligible, selected, drop = FALSE]
  transformed <- qlogis((p * n + 0.5) / (n + 1))
  weights <- pmax(n * p * (1 - p), 1)
  weights[!is.finite(transformed)] <- 0
  transformed[!is.finite(transformed)] <- 0
  local_groups <- factor(c(rep(g1, length(ids1)), rep(g2, length(ids2))), levels = c(g1, g2))
  design <- model.matrix(~0 + local_groups)
  colnames(design) <- c("treatment", "control")
  fit <- lmFit(transformed, design, weights = weights)
  contrast_matrix <- matrix(c(1, -1), ncol = 1, dimnames = list(colnames(design), "treatment_minus_control"))
  fit <- eBayes(contrasts.fit(fit, contrast_matrix), trend = nrow(transformed) > 1L)
  tab <- as.data.table(topTable(fit, number = Inf, sort.by = "none"), keep.rownames = "event_id")
  row_idx <- match(tab$event_id, rownames(psi_mat))
  tab[, `:=`(group_1 = g1, group_2 = g2,
    mean_PSI_group_1 = rowMeans(psi_mat[row_idx, ids1, drop = FALSE], na.rm = TRUE),
    mean_PSI_group_2 = rowMeans(psi_mat[row_idx, ids2, drop = FALSE], na.rm = TRUE))]
  tab[, `:=`(contrast = paste0(g1, "_vs_", g2), delta_PSI = mean_PSI_group_1 - mean_PSI_group_2,
    FDR = p.adjust(P.Value, method = "BH"),
    gene_id = annotation$gene_id[match(event_id, annotation$event_id)])]
  tab[, significant := FDR < 0.05 & abs(delta_PSI) >= 0.10]
  tab[, support_group_1 := apply(support_mat[row_idx, ids1, drop = FALSE], 1, paste, collapse = ";")]
  tab[, support_group_2 := apply(support_mat[row_idx, ids2, drop = FALSE], 1, paste, collapse = ";")]
  results[[length(results) + 1L]] <- tab[, .(contrast, event_id, gene_id, group_1, group_2,
    mean_PSI_group_1, mean_PSI_group_2, delta_PSI, P.Value, FDR, significant, support_group_1, support_group_2)]
}
if (length(results)) out <- rbindlist(results, use.names = TRUE, fill = TRUE) else out <- data.table(
  contrast=character(), event_id=character(), gene_id=character(), group_1=character(), group_2=character(),
  mean_PSI_group_1=numeric(), mean_PSI_group_2=numeric(), delta_PSI=numeric(), P.Value=numeric(), FDR=numeric(),
  significant=logical(), support_group_1=character(), support_group_2=character())
setnames(out, "P.Value", "p_value")
fwrite(out, output, sep = "\t", na = "NA")
