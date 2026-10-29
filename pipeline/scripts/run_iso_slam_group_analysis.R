#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag) {
  pos <- match(flag, args)
  if (is.na(pos) || pos == length(args)) stop("Missing argument: ", flag)
  args[[pos + 1L]]
}
optional_value_after <- function(flag, default = NULL) {
  if (flag %in% args) value_after(flag) else default
}
samples_path <- value_after("--samples")
outdir <- value_after("--outdir")
gtf_path <- value_after("--gtf")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
tables_dir <- file.path(outdir, "tables")
figures_dir <- file.path(outdir, "figures")
metadata_dir <- file.path(outdir, "metadata")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

samples <- fread(samples_path)
required <- c("sample_id", "group", "source_results_root", "caller_results_root")
if (!all(required %in% names(samples))) stop("Sample table is missing required columns")
label_control_group <- optional_value_after("--label-control-group", optional_value_after("--control-group", ""))
if ("--contrasts-json" %in% args) {
  contrast_specs <- jsonlite::fromJSON(value_after("--contrasts-json"), simplifyVector = FALSE)
  if (!length(contrast_specs)) stop("At least one configured contrast is required")
  if (any(!vapply(contrast_specs, function(x) all(c("treatment", "control") %in% names(x)), logical(1)))) {
    stop("Each configured contrast must include treatment and control")
  }
  contrasts <- lapply(contrast_specs, function(x) c(as.character(x$treatment), as.character(x$control)))
} else {
  control_group <- value_after("--control-group")
  treatment_group <- value_after("--treatment-group")
  contrasts <- list(c(treatment_group, control_group))
  if (!nzchar(label_control_group)) label_control_group <- control_group
}
groups <- unique(as.character(samples$group))
if (length(groups) < 2L || anyNA(groups) || any(!nzchar(groups))) stop("Sample sheet must contain at least two named groups")
if (any(vapply(contrasts, function(x) length(x) != 2L || anyNA(x) || any(!nzchar(x)) || x[[1]] == x[[2]], logical(1)))) {
  stop("Each contrast must name two distinct, non-empty groups")
}
contrast_names <- vapply(contrasts, function(x) paste0(x[[1]], "_vs_", x[[2]]), character(1))
if (anyDuplicated(contrast_names)) stop("Configured contrasts must be unique")
contrast_groups <- unique(unlist(contrasts, use.names = FALSE))
if (!all(contrast_groups %in% groups)) stop("Configured contrast refers to a group missing from the sample sheet")
if (nzchar(label_control_group) && !label_control_group %in% groups) stop("Label control group is missing from the sample sheet")
replicates <- samples[, .N, by = group]
if (any(replicates$N < 2L)) stop("At least two biological replicates are required")
samples[, group := factor(group, levels = groups)]
sample_ids <- samples$sample_id
group_palette <- c("#3C5488", "#E64B35", "#00A087", "#4DBBD5", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000")
if (length(groups) > length(group_palette)) stop("Group palette supports at most eight groups")
group_colors <- setNames(group_palette[seq_along(groups)], groups)
layer_titles <- c(gene_abundance = "Gene abundance", gene_labeling = "Gene labeling fraction", isoform_abundance = "Transcript abundance", total_isoform_usage = "Total RNA transcript usage", new_isoform_usage = "New RNA transcript usage")

read_gene_annotation <- function(path) {
  if (!file.exists(path)) stop("Gene annotation GTF does not exist: ", path)
  gtf <- fread(cmd = paste("zgrep -v '^#'", shQuote(path)), sep = "\t", header = FALSE, select = c(3, 9), quote = "", fill = TRUE, showProgress = FALSE)
  setnames(gtf, c("feature_type", "attributes"))
  gtf <- gtf[feature_type == "gene"]
  extract_attr <- function(x, key) {
    pattern <- paste0("(?:^|;[[:space:]]*)", key, " \\\"([^\\\"]+)\\\"")
    found <- regexec(pattern, x, perl = TRUE)
    matches <- regmatches(x, found)
    vapply(matches, function(hit) if (length(hit) >= 2L) hit[[2L]] else NA_character_, character(1))
  }
  annotation <- gtf[, .(
    gene_id = extract_attr(attributes, "gene_id"),
    gene_name = extract_attr(attributes, "gene_name"),
    gene_type = extract_attr(attributes, "gene_type")
  )]
  annotation[is.na(gene_name) | gene_name == "", gene_name := gene_id]
  annotation[is.na(gene_type), gene_type := ""]
  unique(annotation, by = "gene_id")
}

reference_transcripts <- fread(cmd = paste("zgrep -v '^#'", shQuote(gtf_path)), sep="\t", header=FALSE, select=c(3,9), quote="", fill=TRUE, showProgress=FALSE)
reference_transcripts <- reference_transcripts[V3=="transcript", sub('.*transcript_id "([^";]+)".*', '\\1', V9)]
gene_annotation <- read_gene_annotation(gtf_path)
fwrite(gene_annotation, file.path(tables_dir, "gene_annotation.tsv"), sep = "\t")

read_sample <- function(row) {
  sid <- row$sample_id
  root <- file.path(row$source_results_root, sid)
  caller <- file.path(row$caller_results_root, sid, "05_slam_model/04_tool_output/caller")
  gene_ntr <- fread(file.path(root, "05_slam_model/01_tables", paste0(sid, ".grand3_gene_ntr.tsv")))
  gene_dyn <- fread(file.path(root, "06_isoform_dynamics/01_tables", paste0(sid, ".gene_dynamics.tsv")))
  iso_dyn <- fread(file.path(root, "06_isoform_dynamics/01_tables", paste0(sid, ".isoform_dynamics.tsv")))
  catalog <- fread(file.path(root, "04_isoform_catalog/01_tables", paste0(sid, ".isoform_catalog.tsv")))
  catalog <- catalog[filter_status == "PASS" & structural_category == "full-splice_match" & !is.na(associated_transcript) & associated_transcript != "."]
  reference_ids <- unique(catalog[, .(isoform_id, transcript_id = associated_transcript)])
  reference_ids <- reference_ids[transcript_id %in% reference_transcripts]
  iso_dyn <- merge(iso_dyn, reference_ids, by = "isoform_id")
  iso_dyn <- iso_dyn[, .(total_reads = sum(as.numeric(total_reads)), expected_new_reads = sum(as.numeric(expected_new_reads))), by = .(transcript_id, gene_id)]
  iso_dyn[, posterior_new_fraction := fifelse(total_reads > 0, expected_new_reads / total_reads, NA_real_)]
  iso_dyn[, total_isoform_usage := total_reads / sum(total_reads), by = gene_id]
  iso_dyn[, new_isoform_usage := if (sum(expected_new_reads) > 0) expected_new_reads / sum(expected_new_reads) else rep(NA_real_, .N), by = gene_id]
  substitution_path <- file.path(caller, paste0(sid, ".substitution_summary.tsv"))
  if (file.exists(substitution_path)) {
    substitutions <- fread(substitution_path)
    tc <- substitutions[substitution == "T>C"]
    bg <- substitutions[substitution == "diagnostic_background"]
    if (nrow(tc) != 1L || nrow(bg) != 1L || bg$source != "max_ta_tg_rate_background_site_filtered") stop(sid, ": invalid T>C or diagnostic background summary")
    tc_rate <- as.numeric(tc$rate)
    background_rate <- as.numeric(bg$rate)
  } else {
    qc <- fromJSON(file.path(caller, paste0(sid, ".bulk_slam_qc.json")), simplifyVector = FALSE)
    if (!identical(qc$bg_rate_source, "max_ta_tg_rate_background_site_filtered") || !isTRUE(qc$calling_parameters$bg_from_background_sites)) stop(sid, ": legacy QC lacks auditable diagnostic background")
    tc_rate <- if (is.null(qc$global_tc_rate)) 0 else as.numeric(qc$global_tc_rate)
    if (!is.finite(tc_rate) || tc_rate == 0) tc_rate <- as.numeric(qc$total_tc_count) / as.numeric(qc$t_covered)
    background_rate <- as.numeric(qc$bg_tc_rate)
  }
  list(
    gene_ntr = gene_ntr[, .(sample_id = sid, gene_id, support = as.numeric(read_count), value = as.numeric(grand3_binom_ntr_map))],
    gene_abundance = gene_dyn[, .(sample_id = sid, feature_id = gene_id, gene_id, count = as.numeric(assigned_reads))],
    isoform = iso_dyn[, .(
      sample_id = sid, feature_id = transcript_id, gene_id,
      count = as.numeric(total_reads), support = as.numeric(total_reads),
      ntr = as.numeric(posterior_new_fraction), total_usage = as.numeric(total_isoform_usage),
      new_usage = as.numeric(new_isoform_usage), new_support = as.numeric(expected_new_reads)
    )],
    signal = data.table(sample_id = sid, group = as.character(row$group), tc_rate = tc_rate, diagnostic_background_rate = background_rate)
  )
}

loaded <- lapply(seq_len(nrow(samples)), function(i) read_sample(samples[i]))
gene_ntr_long <- rbindlist(lapply(loaded, `[[`, "gene_ntr"))
gene_abundance_long <- rbindlist(lapply(loaded, `[[`, "gene_abundance"))
iso_long <- rbindlist(lapply(loaded, `[[`, "isoform"))
signal <- rbindlist(lapply(loaded, `[[`, "signal"))

matrix_from_long <- function(dt, value_col) {
  if (!nrow(dt)) return(list(matrix=matrix(numeric(), nrow=0, ncol=length(sample_ids), dimnames=list(character(), sample_ids)), annotation=data.table(feature_id=character(),gene_id=character())))
  wide <- dcast(dt, feature_id + gene_id ~ sample_id, value.var = value_col, fill = 0)
  annotation <- wide[, .(feature_id, gene_id)]
  mat <- as.matrix(wide[, ..sample_ids])
  rownames(mat) <- wide$feature_id
  storage.mode(mat) <- "numeric"
  list(matrix = mat, annotation = annotation)
}

gene_counts <- matrix_from_long(gene_abundance_long, "count")
iso_counts <- matrix_from_long(iso_long, "count")
iso_ntr <- matrix_from_long(iso_long, "ntr")
iso_support <- matrix_from_long(iso_long, "support")
iso_total_usage <- matrix_from_long(iso_long, "total_usage")
iso_new_usage <- matrix_from_long(iso_long, "new_usage")
iso_new_support <- matrix_from_long(iso_long, "new_support")

gene_ntr_wide <- dcast(gene_ntr_long, gene_id ~ sample_id, value.var = "value", fill = NA_real_)
gene_support_wide <- dcast(gene_ntr_long, gene_id ~ sample_id, value.var = "support", fill = 0)
gene_ntr_mat <- as.matrix(gene_ntr_wide[, ..sample_ids]); rownames(gene_ntr_mat) <- gene_ntr_wide$gene_id
gene_support_mat <- as.matrix(gene_support_wide[, ..sample_ids]); rownames(gene_support_mat) <- gene_support_wide$gene_id

design <- model.matrix(~0 + group, data = samples)
colnames(design) <- groups

dge <- DGEList(round(gene_counts$matrix), genes = gene_counts$annotation)
keep <- filterByExpr(dge, design = design)
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
dge <- estimateDisp(dge, design, robust = TRUE)
gene_fit <- glmQLFit(dge, design, robust = TRUE)
gene_logcpm <- cpm(dge, log = TRUE, prior.count = 2)

iso_dge <- NULL; iso_fit <- NULL
iso_logcpm <- matrix(numeric(),nrow=0,ncol=length(sample_ids),dimnames=list(character(),sample_ids))
if (nrow(iso_counts$matrix) && all(colSums(iso_counts$matrix)>0)) {
  candidate <- DGEList(round(iso_counts$matrix), genes=iso_counts$annotation)
  iso_keep <- filterByExpr(candidate, design=design)
  if (any(iso_keep)) {
    iso_dge <- calcNormFactors(candidate[iso_keep,,keep.lib.sizes=FALSE])
    iso_dge <- estimateDisp(iso_dge,design,robust=TRUE)
    iso_fit <- glmQLFit(iso_dge,design,robust=TRUE)
    iso_logcpm <- cpm(iso_dge,log=TRUE,prior.count=2)
  }
}

pca <- prcomp(t(gene_logcpm), scale. = TRUE)
pca_scores <- data.table(sample_id = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2])
pca_scores <- merge(pca_scores, samples[, .(sample_id, group)], by = "sample_id", sort = FALSE)
variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
fwrite(pca_scores, file.path(tables_dir, "sample_pca_scores.tsv"), sep = "\t")

place_pca_labels <- function(dt) {
  xr <- diff(range(dt$PC1)); yr <- diff(range(dt$PC2))
  if (!is.finite(xr) || xr == 0) xr <- 1
  if (!is.finite(yr) || yr == 0) yr <- 1
  candidates <- expand.grid(radius = c(0.075, 0.12, 0.17, 0.23), angle = seq(0, 330, by = 30))
  placed <- matrix(numeric(0), ncol = 2)
  dt[, `:=`(label_x = PC1, label_y = PC2)]
  for (i in seq_len(nrow(dt))) {
    best <- NULL
    for (j in seq_len(nrow(candidates))) {
      x <- dt$PC1[i] + cos(candidates$angle[j] * pi / 180) * candidates$radius[j] * xr
      y <- dt$PC2[i] + sin(candidates$angle[j] * pi / 180) * candidates$radius[j] * yr
      label_clear <- !nrow(placed) || all(abs((placed[, 1] - x) / xr) > 0.13 | abs((placed[, 2] - y) / yr) > 0.07)
      point_clear <- all(abs((dt$PC1 - x) / xr) > 0.045 | abs((dt$PC2 - y) / yr) > 0.045)
      if (label_clear && point_clear) { best <- c(x, y); break }
    }
    if (is.null(best)) best <- c(x, y)
    dt[i, `:=`(label_x = best[1], label_y = best[2])]
    placed <- rbind(placed, best)
  }
  dt[]
}
pca_scores <- place_pca_labels(pca_scores)
fwrite(pca_scores, file.path(tables_dir, "sample_pca_scores.tsv"), sep = "\t")

cor_mat <- cor(gene_logcpm, method = "spearman")
cor_long <- as.data.table(as.table(cor_mat)); setnames(cor_long, c("sample_1", "sample_2", "correlation"))
fwrite(cor_long, file.path(tables_dir, "sample_gene_abundance_correlations.tsv"), sep = "\t")

signal <- merge(signal, samples[, .(sample_id, group)], by = c("sample_id", "group"), sort = FALSE)
fwrite(signal, file.path(tables_dir, "sample_labeling_signal.tsv"), sep = "\t")
within_signal <- signal[, .(
  n = .N,
  tc_rate_mean = mean(tc_rate), tc_rate_sd = sd(tc_rate), tc_rate_cv = sd(tc_rate) / mean(tc_rate),
  diagnostic_background_mean = mean(diagnostic_background_rate), diagnostic_background_sd = sd(diagnostic_background_rate),
  diagnostic_background_cv = sd(diagnostic_background_rate) / mean(diagnostic_background_rate)
), by = group]
fwrite(within_signal, file.path(tables_dir, "within_group_signal_summary.tsv"), sep = "\t")

theme_report <- theme_classic(base_size = 11) + theme(legend.position = "right", plot.title = element_text(face = "bold", size = 12))
save_plot <- function(plot, stem, width = 7.4, height = 6.2) {
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 450, bg = "white")
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, bg = "white")
}

p_pca <- ggplot(pca_scores, aes(PC1, PC2, color = group)) +
  geom_segment(aes(xend = label_x, yend = label_y), color = "#7B8794", linewidth = 0.3, show.legend = FALSE) +
  geom_point(size = 3.2, alpha = 0.88) +
  geom_label(aes(label = sample_id, x = label_x, y = label_y), size = 2.8, fill = "white", linewidth = 0.18, label.padding = unit(0.12, "lines"), show.legend = FALSE) +
  scale_color_manual(values = group_colors) + coord_equal() +
  labs(title = "Gene abundance PCA", x = sprintf("PC1 (%.1f%%)", variance[1]), y = sprintf("PC2 (%.1f%%)", variance[2])) + theme_report
save_plot(p_pca, file.path(figures_dir, "gene_abundance_pca"))

p_cor <- ggplot(cor_long, aes(sample_1, sample_2, fill = correlation)) + geom_tile() +
  scale_fill_gradientn(colors = c("#E5E7EB", "#4DBBD5", "#3C5488", "#E64B35"), limits = c(min(cor_long$correlation), 1)) +
  coord_equal() + labs(title = "Sample correlation", x = NULL, y = NULL, fill = "Spearman") +
  theme_report + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_cor, file.path(figures_dir, "gene_abundance_correlation"), 7.4, 6.7)

p_signal <- melt(signal, id.vars = c("sample_id", "group"), measure.vars = c("tc_rate", "diagnostic_background_rate"), variable.name = "metric", value.name = "rate")
p_signal[, metric := factor(metric, levels = c("tc_rate", "diagnostic_background_rate"), labels = c("T>C conversion", "Background"))]
p_sig <- ggplot(p_signal, aes(group, rate * 100, color = group)) + geom_point(size = 2.7, alpha = 0.9) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.45, linewidth = 0.4) + facet_wrap(~metric, scales = "free_y") +
  scale_color_manual(values = group_colors, guide = "none") + labs(title = "Within-group labeling signal", x = NULL, y = "Rate (%)") + theme_report
save_plot(p_sig, file.path(figures_dir, "within_group_labeling_signal"), 8.2, 4.8)

contrast_vector <- function(g1, g2) {
  value <- setNames(rep(0, length(groups)), groups)
  value[[g1]] <- 1
  value[[g2]] <- -1
  value
}

empty_result <- function(g1, g2, layer) {
  data.table(feature_id=character(), gene_id=character(), group_1=character(), group_2=character(), layer=character(), effect=numeric(), effect_type=character(), mean_group_1=numeric(), mean_group_2=numeric(), FDR=numeric(), significant=logical())
}

edge_result <- function(fit, dge_object, g1, g2, layer) {
  if (is.null(fit)) return(empty_result(g1,g2,layer))
  qlf <- glmQLFTest(fit, contrast = contrast_vector(g1, g2))
  tab <- as.data.table(topTags(qlf, n = Inf, sort.by = "none")$table, keep.rownames = "row_id")
  if ("feature_id" %in% names(tab)) {
    tab[, row_id := NULL]
  } else {
    setnames(tab, "row_id", "feature_id")
  }
  tab[, `:=`(group_1 = g1, group_2 = g2, layer = layer, effect = logFC, effect_type = "log2FC")]
  ids1 <- samples[group == g1, sample_id]; ids2 <- samples[group == g2, sample_id]
  abundance <- cpm(dge_object, log = TRUE, prior.count = 2)
  tab[, `:=`(mean_group_1 = rowMeans(abundance[feature_id, ids1, drop = FALSE]), mean_group_2 = rowMeans(abundance[feature_id, ids2, drop = FALSE]))]
  tab[, significant := FDR < 0.05 & abs(effect) >= 1]
  tab[]
}

proportion_result <- function(value_mat, support_mat, annotation, g1, g2, layer, delta_threshold, min_support = 20) {
  ids1 <- samples[group == g1, sample_id]; ids2 <- samples[group == g2, sample_id]
  common <- rownames(value_mat)
  support_mat <- support_mat[common, , drop = FALSE]
  eligible <- rowSums(support_mat[, ids1, drop = FALSE] >= min_support) >= 2 & rowSums(support_mat[, ids2, drop = FALSE] >= min_support) >= 2
  if (!any(eligible)) return(empty_result(g1, g2, layer))
  p <- value_mat[eligible, c(ids1, ids2), drop = FALSE]
  n <- support_mat[eligible, c(ids1, ids2), drop = FALSE]
  p[!is.finite(p)] <- NA_real_
  if (!any(is.finite(p))) return(empty_result(g1, g2, layer))
  transformed <- qlogis((p * n + 0.5) / (n + 1))
  weights <- pmax(n * p * (1 - p), 1)
  weights[!is.finite(transformed)] <- 0
  transformed[!is.finite(transformed)] <- 0
  local_groups <- factor(c(rep(g1, length(ids1)), rep(g2, length(ids2))), levels = c(g1, g2))
  local_design <- model.matrix(~0 + local_groups); colnames(local_design) <- c("treatment", "control")
  fit <- lmFit(transformed, local_design, weights = weights)
  fit <- contrasts.fit(fit, contrasts = matrix(c(1, -1), ncol = 1, dimnames = list(colnames(local_design), "treatment_minus_control")))
  fit <- eBayes(fit, trend = nrow(transformed)>1)
  tab <- as.data.table(topTable(fit, number = Inf, sort.by = "none"), keep.rownames = "feature_id")
  raw <- value_mat[tab$feature_id, , drop = FALSE]
  tab[, `:=`(
    group_1 = g1, group_2 = g2, layer = layer,
    mean_group_1 = rowMeans(raw[, ids1, drop = FALSE], na.rm = TRUE),
    mean_group_2 = rowMeans(raw[, ids2, drop = FALSE], na.rm = TRUE)
  )]
  tab[, `:=`(effect = mean_group_1 - mean_group_2, effect_type = "delta_proportion", FDR = p.adjust(P.Value, method = "BH"))]
  tab[, significant := FDR < 0.05 & abs(effect) >= delta_threshold]
  if (!is.null(annotation)) tab <- merge(tab, annotation, by = "feature_id", all.x = TRUE, sort = FALSE)
  tab[]
}

plot_differential <- function(tab, matrix_values, sample_order, stem, title) {
  tab[, log10_fdr := -log10(pmax(FDR, 1e-300))]
  p <- ggplot(tab, aes(effect, log10_fdr, color = significant)) + geom_point(size = 1.1, alpha = 0.68) +
    scale_color_manual(values = c(`FALSE` = "#9AA4B2", `TRUE` = "#E64B35")) +
    labs(title = title, x = unique(tab$effect_type), y = "-log10(FDR)", color = "Significant") + theme_report
  save_plot(p, paste0(stem, "_volcano"), 7.4, 5.8)
  selected <- head(tab[order(FDR, -abs(effect)), feature_id], 20L)
  selected <- selected[selected %in% rownames(matrix_values)]
  if (!length(selected)) {
    fwrite(data.table(feature_id=character(), sample_id=character(), z_score=numeric(), gene_id=character(), gene_name=character(), group=character()), paste0(stem, "_heatmap_data.tsv"), sep="\t")
    return(invisible(NULL))
  }
  z <- matrix_values[selected, sample_order, drop = FALSE]
  z <- t(scale(t(z))); z[!is.finite(z)] <- 0
  heat <- as.data.table(as.table(z)); setnames(heat, c("feature_id", "sample_id", "z_score"))
  heat <- merge(heat, unique(tab[, .(feature_id, gene_id, gene_name)]), by = "feature_id", all.x = TRUE, sort = FALSE)
  heat[, group := as.character(samples$group[match(sample_id, samples$sample_id)])]
  fwrite(heat, paste0(stem, "_heatmap_data.tsv"), sep = "\t")
  heat[, display_label := fifelse(!is.na(gene_name) & gene_name != gene_id, paste0(gene_name, " (", feature_id, ")"), feature_id)]
  display_levels <- unique(heat$display_label[match(rev(selected), heat$feature_id)])
  heat[, display_label := factor(display_label, levels = display_levels)]
  ph <- ggplot(heat, aes(sample_id, display_label, fill = z_score)) + geom_tile() +
    scale_fill_gradient2(low = "#3C5488", mid = "#F7F7F7", high = "#E64B35", midpoint = 0) +
    labs(title = paste(title, "Top features"), x = NULL, y = NULL, fill = "Row z-score") +
    theme_report + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 7))
  save_plot(ph, paste0(stem, "_heatmap"), 8.4, 6.8)
}

contrast_summary <- list()
for (pair in contrasts) {
  g1 <- pair[[1]]; g2 <- pair[[2]]; cname <- paste0(g1, "_vs_", g2)
  cdir <- file.path(tables_dir, "differential", cname); fdir <- file.path(figures_dir, "differential", cname)
  dir.create(cdir, recursive = TRUE, showWarnings = FALSE); dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
  results <- list(
    gene_abundance = edge_result(gene_fit, dge, g1, g2, "gene_abundance"),
    gene_labeling = proportion_result(gene_ntr_mat, gene_support_mat, NULL, g1, g2, "gene_labeling", 0.05),
    isoform_abundance = edge_result(iso_fit, iso_dge, g1, g2, "isoform_abundance"),
    total_isoform_usage = proportion_result(iso_total_usage$matrix, iso_support$matrix, iso_total_usage$annotation, g1, g2, "total_isoform_usage", 0.10),
    new_isoform_usage = proportion_result(iso_new_usage$matrix, iso_new_support$matrix, iso_new_usage$annotation, g1, g2, "new_isoform_usage", 0.10)
  )
  matrices <- list(gene_abundance = gene_logcpm, gene_labeling = gene_ntr_mat, isoform_abundance = iso_logcpm, total_isoform_usage = iso_total_usage$matrix, new_isoform_usage = iso_new_usage$matrix)
  for (layer in names(results)) {
    tab <- results[[layer]]
    if (!"gene_id" %in% names(tab)) tab[, gene_id := feature_id]
    tab <- merge(tab, gene_annotation, by = "gene_id", all.x = TRUE, sort = FALSE)
    tab[is.na(gene_name) | gene_name == "", gene_name := gene_id]
    fwrite(tab, file.path(cdir, paste0(layer, ".tsv")), sep = "\t")
    contrast_summary[[length(contrast_summary) + 1L]] <- data.table(contrast = cname, group_1 = g1, group_2 = g2, layer = layer, tested_features = nrow(tab), significant_features = sum(tab$significant, na.rm = TRUE))
    plot_differential(tab, matrices[[layer]], sample_ids, file.path(fdir, layer), paste(gsub("_vs_", " vs ", cname), layer_titles[[layer]]))
  }
}
contrast_summary <- rbindlist(contrast_summary)
fwrite(contrast_summary, file.path(tables_dir, "contrast_summary.tsv"), sep = "\t")

contrast_summary[, layer_display := factor(layer_titles[layer], levels = unname(layer_titles))]
p_counts <- ggplot(contrast_summary, aes(reorder(gsub("_vs_", " vs ", contrast), significant_features), significant_features, fill = layer_display)) +
  geom_col(position = "dodge") + coord_flip() + labs(title = "Significant differential features", x = NULL, y = "Feature count", fill = "Layer") + theme_report
save_plot(p_counts, file.path(figures_dir, "differential_feature_counts"), 9.2, 6.4)

tested_gene_ids <- unique(unlist(lapply(list(gene_counts$annotation$gene_id, iso_counts$annotation$gene_id), as.character)))
annotation_coverage <- mean(tested_gene_ids %in% gene_annotation$gene_id)
if (!is.finite(annotation_coverage) || annotation_coverage < 0.95) stop(sprintf("Gene annotation coverage %.2f%% is below 95%%", 100 * annotation_coverage))
write_json(list(
  status = "PASS", sample_count = nrow(samples), groups = groups,
  biological_replicates_per_group = setNames(as.list(replicates$N), as.character(replicates$group)), label_control_group = label_control_group, transcript_comparison_scope = "unique_reference_full_splice_matches", usage_denominator = "comparable_reference_transcripts_within_gene", batch_covariate = FALSE,
  contrast_details = lapply(contrasts, function(x) list(treatment = x[[1]], control = x[[2]])),
  contrasts = as.list(contrast_names),
  layers = unique(contrast_summary$layer),
  gene_features_union = nrow(gene_ntr_mat), isoform_features_union = nrow(iso_counts$matrix),
  gene_features_present_in_all_samples = sum(rowSums(gene_support_mat > 0) == ncol(gene_support_mat)),
  isoform_features_present_in_all_samples = sum(rowSums(iso_counts$matrix > 0) == ncol(iso_counts$matrix)),
  gene_annotation_coverage = annotation_coverage,
  thresholds = list(abundance = "FDR<0.05 and abs(log2FC)>=1", gene_labeling = "FDR<0.05 and abs(delta_NTR)>=0.05", isoform_usage = "FDR<0.05 and abs(delta_usage)>=0.10")
), file.path(metadata_dir, "group_analysis_validation.json"), auto_unbox = TRUE, pretty = TRUE)
