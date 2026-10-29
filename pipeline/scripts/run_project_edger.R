#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(edgeR))

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) stop("Usage: run_project_edger.R <results_root> <samples.tsv> <mode>")
results_root <- args[[1]]
samples_path <- args[[2]]
mode <- args[[3]]
if (mode != "full") quit(status=0)

cohort <- file.path(results_root, "_project", "06_cohort_qc")
de_stage <- file.path(results_root, "_project", "07_differential_expression")
for (d in c("01_tables", "02_figures", "03_logs", "04_tool_output", "05_metadata")) dir.create(file.path(cohort,d), recursive=TRUE, showWarnings=FALSE)
for (d in c("01_tables", "02_figures", "03_logs", "04_tool_output", "05_metadata")) dir.create(file.path(de_stage,d), recursive=TRUE, showWarnings=FALSE)
samples <- read.delim(samples_path, check.names=FALSE, stringsAsFactors=FALSE)
samples$group <- factor(samples$group, levels=c("A","B","AB","NC"))
rownames(samples) <- samples$sample_id
group_colors <- c(A="#E64B35", B="#4DBBD5", AB="#00A087", NC="#3C5488")

read_counts <- function(path) {
  x <- read.delim(path, check.names=FALSE, stringsAsFactors=FALSE)
  gene_name <- x$gene_name
  mat <- as.matrix(x[, samples$sample_id, drop=FALSE])
  storage.mode(mat) <- "integer"
  rownames(mat) <- x$gene_id
  list(counts=mat, gene_name=setNames(gene_name, x$gene_id))
}

plot_qc <- function(dge, prefix, title) {
  logcpm <- cpm(dge, log=TRUE, prior.count=2)
  pca <- prcomp(t(logcpm), scale.=FALSE)
  pdf(file.path(cohort,"02_figures",paste0(prefix,"_pca.pdf")), 7.4, 6.2)
  plot(pca$x[,1], pca$x[,2], pch=19, cex=1.4, col=group_colors[as.character(samples[colnames(logcpm),"group"])],
       xlab=sprintf("PC1 (%.1f%%)",100*summary(pca)$importance[2,1]), ylab=sprintf("PC2 (%.1f%%)",100*summary(pca)$importance[2,2]), main=title)
  text(pca$x[,1], pca$x[,2], labels=rownames(pca$x), pos=3, cex=.8, col="#333333")
  legend("topright", legend=names(group_colors), col=group_colors, pch=19, bty="n")
  dev.off()
  png(file.path(cohort,"02_figures",paste0(prefix,"_pca.png")), width=3330, height=2790, res=450, bg="white")
  plot(pca$x[,1], pca$x[,2], pch=19, cex=1.4, col=group_colors[as.character(samples[colnames(logcpm),"group"])],
       xlab=sprintf("PC1 (%.1f%%)",100*summary(pca)$importance[2,1]), ylab=sprintf("PC2 (%.1f%%)",100*summary(pca)$importance[2,2]), main=title)
  text(pca$x[,1], pca$x[,2], labels=rownames(pca$x), pos=3, cex=.8, col="#333333")
  legend("topright", legend=names(group_colors), col=group_colors, pch=19, bty="n"); dev.off()
  write.table(cor(logcpm, method="spearman"), file.path(cohort,"01_tables",paste0(prefix,"_spearman.tsv")), sep="\t", quote=FALSE, col.names=NA)
}

run_one <- function(kind, path) {
  obj <- read_counts(path)
  y <- DGEList(obj$counts, group=samples[colnames(obj$counts),"group"])
  keep <- filterByExpr(y, group=y$samples$group)
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y, method="TMM")
  design <- model.matrix(~0+group, data=samples[colnames(y),,drop=FALSE])
  colnames(design) <- sub("group", "", colnames(design))
  y <- estimateDisp(y, design, robust=TRUE)
  fit <- glmQLFit(y, design, robust=TRUE)
  plot_qc(y, kind, paste(kind, "RNA PCA"))
  contrasts <- c("A-B","A-AB","A-NC","B-AB","B-NC","AB-NC")
  index <- data.frame()
  for (label in contrasts) {
    con <- makeContrasts(contrasts=label, levels=design)
    qlf <- glmQLFTest(fit, contrast=con)
    tab <- topTags(qlf, n=Inf, sort.by="PValue")$table
    tab$gene_id <- rownames(tab)
    tab$gene_name <- unname(obj$gene_name[tab$gene_id])
    tab$significant_fdr_0_05 <- tab$FDR < 0.05
    tab <- tab[,c("gene_id","gene_name","logFC","logCPM","F","PValue","FDR","significant_fdr_0_05")]
    safe <- gsub("-","_vs_",label)
    out <- file.path(de_stage,"01_tables",paste0(kind,"_",safe,".tsv"))
    write.table(tab, out, sep="\t", quote=FALSE, row.names=FALSE)
    index <- rbind(index, data.frame(rna_type=kind, comparison=label, tested=nrow(tab), significant=sum(tab$FDR<0.05), path=out))
  }
  write.table(index, file.path(de_stage,"01_tables",paste0(kind,"_contrast_index.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
}

run_one("total", file.path(cohort,"01_tables","total_rna_counts.tsv"))
run_one("new", file.path(cohort,"01_tables","new_rna_counts.tsv"))
