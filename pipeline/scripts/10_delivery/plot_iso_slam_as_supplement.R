#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: plot_iso_slam_as_supplement.R <differential.tsv> <sample_psi.tsv> <output_dir>")
differential_path <- args[[1]]
sample_psi_path <- args[[2]]
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

group_palette <- c("#3C5488", "#E64B35", "#00A087", "#4DBBD5", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000")
event_palette <- c(SE = "#3C5488", A3 = "#E64B35", A5 = "#00A087", MX = "#4DBBD5", RI = "#F39B7F")
theme_report <- theme_bw(base_size = 10) + theme(
  panel.grid.minor = element_blank(),
  panel.grid.major.y = element_blank(),
  strip.background = element_rect(fill = "#F2F4F7", colour = NA),
  strip.text = element_text(face = "bold", colour = "#243247"),
  axis.text = element_text(colour = "#263445"),
  axis.title = element_text(colour = "#263445"),
  plot.title = element_text(face = "bold", size = 13, colour = "#243247"),
  legend.position = "bottom"
)

diff <- fread(differential_path, sep = "\t", na.strings = c("NA", ""), showProgress = FALSE)
required_diff <- c("contrast", "event_id", "gene_name", "delta_PSI", "FDR", "significant")
if (!all(required_diff %in% names(diff))) stop("Differential table is missing required columns")
diff[, `:=`(delta_PSI = as.numeric(delta_PSI), FDR = as.numeric(FDR))]
diff[, significant := tolower(as.character(significant)) == "true"]
contrast_order <- unique(diff$contrast)
tested <- diff[, .(count = .N), by = "contrast"][match(contrast_order, contrast)]
tested[, metric := "Events tested"]
sig_counts <- diff[significant == TRUE, .(count = .N), by = "contrast"]
sig_counts <- merge(data.table(contrast = contrast_order), sig_counts, by = "contrast", all.x = TRUE, sort = FALSE)
sig_counts[is.na(count), count := 0L]
sig_counts[, metric := "Significant events"]
overview <- rbindlist(list(tested, sig_counts), use.names = TRUE)
overview[, contrast := factor(contrast, levels = contrast_order)]
overview[, metric := factor(metric, levels = c("Events tested", "Significant events"))]
overview_plot <- ggplot(overview, aes(x = contrast, y = count, fill = metric)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = count), vjust = -0.35, size = 3.2, colour = "#263445") +
  facet_wrap(~metric, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("Events tested" = "#3C5488", "Significant events" = "#E64B35"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "New-RNA alternative-splicing comparison overview", x = NULL, y = "Number of events") +
  theme_report + theme(axis.text.x = element_text(angle = 28, hjust = 1))
ggsave(file.path(output_dir, "AS_contrast_overview.png"), overview_plot, width = 10, height = 6, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "AS_contrast_overview.pdf"), overview_plot, width = 10, height = 6, bg = "white")

sig <- copy(diff[significant == TRUE])
sig[, event_type := sub(".*;([A-Z0-9]+):.*", "\\1", event_id)]
if (nrow(sig)) {
  setorder(sig, contrast, delta_PSI)
  sig[, row_id := seq_len(.N), by = "contrast"]
  sig[, event_label := paste0(row_id, "  ", gene_name, " [", event_type, "]")]
  sig[, event_key := paste(contrast, row_id, sep = "::")]
  sig[, event_key := factor(event_key, levels = unique(event_key))]
  delta_plot <- ggplot(sig, aes(x = delta_PSI, y = event_key, colour = event_type)) +
    geom_vline(xintercept = 0, colour = "#AEB7C2", linewidth = 0.4) +
    geom_segment(aes(x = 0, xend = delta_PSI, yend = event_key), linewidth = 0.65, alpha = 0.8) +
    geom_point(size = 2.2) +
    facet_wrap(~contrast, scales = "free_y", ncol = 1) +
    scale_y_discrete(labels = setNames(sig$event_label, sig$event_key)) +
    scale_colour_manual(values = event_palette, drop = FALSE) +
    labs(title = "Significant new-RNA splicing events", x = "Delta PSI (group 1 - group 2)", y = NULL, colour = "Event type") +
    theme_report + theme(axis.text.y = element_text(size = 7.5), legend.position = "bottom")
  delta_height <- max(8, 2 + nrow(sig) * 0.22)
  ggsave(file.path(output_dir, "AS_significant_delta_psi.png"), delta_plot, width = 11, height = delta_height, dpi = 300, limitsize = FALSE, bg = "white")
  ggsave(file.path(output_dir, "AS_significant_delta_psi.pdf"), delta_plot, width = 11, height = delta_height, limitsize = FALSE, bg = "white")
} else {
  empty_plot <- ggplot() + annotate("text", x = 0, y = 0, label = "No significant events") +
    labs(title = "Significant new-RNA splicing events", x = NULL, y = NULL) + theme_void() + theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(output_dir, "AS_significant_delta_psi.png"), empty_plot, width = 9, height = 3, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, "AS_significant_delta_psi.pdf"), empty_plot, width = 9, height = 3, bg = "white")
}

if (nrow(sig)) {
  candidate_events <- unique(sig[order(FDR), .(event_id, gene_name, event_type)], by = "event_id")
  candidate_events[, event_order := seq_len(.N)]
  candidate_events[, event_label := sprintf("%02d  %s [%s]", event_order, gene_name, event_type)]
  psi <- fread(
    sample_psi_path, sep = "\t",
    select = c("event_id", "gene_name", "sample_id", "group", "new_RNA_PSI"),
    na.strings = c("NA", ""), showProgress = FALSE
  )
  psi <- psi[event_id %chin% candidate_events$event_id]
  psi <- merge(psi, candidate_events[, .(event_id, event_order, event_label)], by = "event_id", all = FALSE, sort = FALSE)
  psi[, `:=`(new_RNA_PSI = as.numeric(new_RNA_PSI), group = factor(group, levels = unique(group)))]
  psi[, replicate := factor(seq_len(.N), levels = seq_len(max(.N))), by = .(event_id, group)]
  psi[, event_label := factor(event_label, levels = candidate_events$event_label)]
  groups <- levels(psi$group)
  if (length(groups) > length(group_palette)) stop("Group palette supports at most eight groups")
  group_colors <- setNames(group_palette[seq_along(groups)], groups)
  page_size <- 7L
  page_count <- ceiling(nrow(candidate_events) / page_size)
  pdf(file.path(output_dir, "AS_candidate_sample_psi.pdf"), width = 11, height = 9, onefile = TRUE)
  for (page in seq_len(page_count)) {
    first <- (page - 1L) * page_size + 1L
    last <- min(page * page_size, nrow(candidate_events))
    labels <- candidate_events$event_label[first:last]
    page_data <- psi[event_label %in% labels]
    psi_plot <- ggplot(page_data, aes(x = group, y = new_RNA_PSI, colour = group, shape = replicate)) +
      geom_point(position = position_jitter(width = 0.08, height = 0, seed = 17), size = 2.1, na.rm = TRUE) +
      facet_wrap(~event_label, ncol = 1) +
      scale_colour_manual(values = group_colors, drop = FALSE) +
      scale_shape_manual(values = c("1" = 16, "2" = 17, "3" = 15, "4" = 18)) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0.02, 0.06))) +
      labs(title = "Per-sample new-RNA PSI for significant candidate events", x = NULL, y = "New-RNA PSI", colour = "Group", shape = "Replicate") +
      theme_report + theme(axis.text.x = element_text(angle = 25, hjust = 1), strip.text = element_text(size = 9))
    page_name <- sprintf("AS_candidate_sample_psi_page%d.png", page)
    ggsave(file.path(output_dir, page_name), psi_plot, width = 11, height = 9, dpi = 300, bg = "white")
    print(psi_plot)
  }
  dev.off()
}
