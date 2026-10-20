#!/usr/bin/env Rscript
# Heatmaps of the protein-coding RNA-seq response (IFNB vs Mock), across six conditions, columns left to right:
#   Mock (1-3) | IFNB (10-12) | IFNB+EBOV (16-18) | EBOV (4-6) | IFNB+MARV (22-24) | MARV (7-9)
#
#   1. heatmap_top50_ISGs_...      the 50 most upregulated IFNB-vs-Mock genes (padj < 0.05, log2FC > 1), ranked by
#                                  log2FC: first row = most upregulated. Gene names (and the IFNB log2FC) on the rows.
#   2. heatmap_top50_IDGs_...      the 50 most downregulated genes (padj < 0.05, log2FC < -1), most negative first.
#   3. heatmap_all_protein_coding  every protein-coding gene, no gene labels. Rows are grouped by k-means cluster
#                                  (K_CLUSTERS), clusters ordered from most to least IFNB-induced, genes ordered
#                                  by similarity within a cluster.
# Colour = log2 change of each sample's VST expression relative to the MEAN of the three Mock samples (blue = lower than
# Mock, red = higher), capped at the 98th percentile of |change| (automatic; override with CAP_*). Set
# VALUE <- "zscore" to show per-gene z-scores across the 18 samples instead.
# Gene rankings come from DE_v2_IFNB_vs_Mock_protein_coding.csv (shrunk log2FC, from 04/07/08); expression from vsd_v2.rds.
# A check table compares the sample numbers with the group recorded for each in your data.
# Needs pheatmap. Usage: Rscript 47_heatmaps_isg_idg_transcriptome.R   (PROJ from 'source setup_session.sh')
# Writes to <project>/requant_v2_outputs/figures_heatmaps/ (PNG + PDF) and the gene lists as CSV.

suppressPackageStartupMessages({ library(DESeq2); library(readr); library(dplyr); library(tidyr) })
if (!requireNamespace("pheatmap", quietly = TRUE)) stop("Install pheatmap (conda install -c conda-forge r-pheatmap)")

N_TOP <- 50; LFC_CUT <- 1; PADJ_CUT <- 0.05; K_CLUSTERS <- 12
VALUE <- "log2_vs_mock"          # or "zscore"
CAP_ISG <- NULL; CAP_IDG <- NULL; CAP_ALL <- NULL

SETS <- list(Mock = c("1", "2", "3"), IFNB = c("10", "11", "12"), `IFNB+EBOV` = c("16", "17", "18"),
             EBOV = c("4", "5", "6"), `IFNB+MARV` = c("22", "23", "24"), MARV = c("7", "8", "9"))   # column order
EXPECTED <- c(Mock = "Mock_Mock", IFNB = "IFNB_Mock", `IFNB+EBOV` = "IFNB_EBOV", EBOV = "Mock_EBOV",
              `IFNB+MARV` = "IFNB_MARV", MARV = "Mock_MARV")
GCOL <- c(Mock = "#7F7F7F", IFNB = "#2B6CB0", `IFNB+EBOV` = "#9BD08E", EBOV = "#1B7837", `IFNB+MARV` = "#C9A6DB", MARV = "#762A83")

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_heatmaps"); dir.create(OUT, showWarnings = FALSE)

vsd <- readRDS(file.path(RNA, "vsd_v2.rds")); vm <- assay(vsd)
sel <- tibble(sample = unlist(SETS), group = rep(names(SETS), lengths(SETS)))
sel$data_group <- as.character(vsd$group)[match(sel$sample, colnames(vm))]
if (any(is.na(sel$data_group))) stop("Samples not in the VST object: ", paste(sel$sample[is.na(sel$data_group)], collapse = ", "))
message("Check - chosen group vs group recorded in the data:")
print(as.data.frame(sel |> transmute(sample, chosen_as = group, group_in_data = data_group,
                                     match = ifelse(data_group == EXPECTED[group], "ok", "DIFFERENT"))), row.names = FALSE)
if (any(sel$data_group != EXPECTED[sel$group])) message("WARNING: some samples are recorded under a different group. Using YOUR sample numbers.")

de <- read_csv(file.path(RNA, "DE_v2_IFNB_vs_Mock_protein_coding.csv"), show_col_types = FALSE)
de <- de |> filter(gene_id %in% rownames(vm), !is.na(padj))
x <- vm[, sel$sample, drop = FALSE]
mock_cols <- sel$sample[sel$group == "Mock"]

value_matrix <- function(ids) {
  m <- x[ids, , drop = FALSE]
  if (VALUE == "zscore") { z <- t(scale(t(m))); z[stats::complete.cases(z), , drop = FALSE] }
  else m - rowMeans(m[, mock_cols, drop = FALSE])
}
auto_cap <- function(m, override) if (!is.null(override)) override else min(8, max(2, ceiling(quantile(abs(m), 0.98, na.rm = TRUE))))
ann_col <- data.frame(Group = factor(sel$group, levels = names(SETS)), row.names = sel$sample)
gaps <- cumsum(lengths(SETS))[-length(SETS)]
pal <- colorRampPalette(c("#2C5F9E", "#F7F7F7", "#B2382C"))(101)
lab_unit <- if (VALUE == "zscore") "z-score of VST" else "log2 change vs Mock mean (VST)"

draw <- function(m, rows_lab, name, title, cap, w, h, extra = list()) {
  brk <- seq(-cap, cap, length.out = 102)
  args <- c(list(m, color = pal, breaks = brk, cluster_rows = FALSE, cluster_cols = FALSE, annotation_col = ann_col,
                 annotation_colors = list(Group = GCOL), gaps_col = unname(gaps), show_colnames = TRUE, border_color = NA,
                 main = title, silent = TRUE, fontsize = 9), extra)
  if (!is.null(rows_lab)) { args$show_rownames <- TRUE; args$labels_row <- rows_lab; args$fontsize_row <- 8 } else args$show_rownames <- FALSE
  p <- do.call(pheatmap::pheatmap, args)
  for (ext in c("png", "pdf")) {
    f <- file.path(OUT, paste0(name, ".", ext))
    if (ext == "png") png(f, width = w, height = h, units = "in", res = 300) else pdf(f, width = w, height = h)
    grid::grid.newpage(); grid::grid.draw(p$gtable); dev.off()
  }
}

# ---- 1 and 2: top ISGs / IDGs ------------------------------------------------------------------------------------------
top_set <- function(up) {
  g <- if (up) de |> filter(padj < PADJ_CUT, log2FoldChange > LFC_CUT) |> arrange(desc(log2FoldChange))
       else de |> filter(padj < PADJ_CUT, log2FoldChange < -LFC_CUT) |> arrange(log2FoldChange)
  message(sprintf("%s: %d genes pass padj<%g & log2FC %s %d; plotting the top %d", if (up) "ISGs" else "IDGs", nrow(g), PADJ_CUT,
                  if (up) ">" else "<", if (up) LFC_CUT else -LFC_CUT, min(N_TOP, nrow(g))))
  head(g, N_TOP)
}
for (up in c(TRUE, FALSE)) {
  g <- top_set(up); tag <- if (up) "ISGs" else "IDGs"
  m <- value_matrix(g$gene_id); g <- g[match(rownames(m), g$gene_id), ]
  labs <- make.unique(sprintf("%s  (%+.1f)", g$gene_name, g$log2FoldChange))
  write_csv(g |> transmute(rank = row_number(), gene_id, gene_name, log2FoldChange, padj),
            file.path(OUT, paste0("top", N_TOP, "_", tag, "_IFNB_vs_Mock.csv")))
  draw(m, labs, paste0("heatmap_top", N_TOP, "_", tag, "_IFNB_vs_Mock"),
       sprintf("Top %d %s by IFNB vs Mock log2FC (%s; label = gene (IFNB log2FC))", nrow(m), if (up) "upregulated genes (ISGs)" else "downregulated genes (IDGs)", lab_unit),
       auto_cap(m, if (up) CAP_ISG else CAP_IDG), 8, 0.2 * nrow(m) + 3.2)
}

# ---- 3: whole protein-coding transcriptome -----------------------------------------------------------------------------
allm <- value_matrix(de$gene_id)
message("Whole-transcriptome heatmap: ", nrow(allm), " protein-coding genes")
set.seed(1)
km <- kmeans(allm, centers = K_CLUSTERS, nstart = 10, iter.max = 100)
ifnb_cols <- sel$sample[sel$group == "IFNB"]
score <- vapply(seq_len(K_CLUSTERS), function(k) mean(allm[km$cluster == k, ifnb_cols, drop = FALSE]), numeric(1))
new_id <- match(seq_len(K_CLUSTERS), order(score, decreasing = TRUE))
cl <- new_id[km$cluster]; names(cl) <- rownames(allm)
ord <- unlist(lapply(seq_len(K_CLUSTERS), function(k) {
  g <- names(cl)[cl == k]
  if (length(g) > 2) g[hclust(dist(allm[g, , drop = FALSE]), method = "ward.D2")$order] else g
}))
allm <- allm[ord, , drop = FALSE]
cl_f <- factor(paste0("C", cl[ord]), levels = paste0("C", seq_len(K_CLUSTERS)))
sym <- setNames(de$gene_name, de$gene_id)
write_csv(tibble(gene_id = rownames(allm), gene_name = unname(sym[rownames(allm)]), cluster = as.character(cl_f)),
          file.path(OUT, "heatmap_all_protein_coding_clusters.csv"))
ann_row <- data.frame(Cluster = cl_f, row.names = rownames(allm))
cl_cols <- setNames(grDevices::hcl.colors(K_CLUSTERS, "Dark 3"), levels(cl_f))
ann_colors <- list(Group = GCOL, Cluster = cl_cols)
brk <- seq(-auto_cap(allm, CAP_ALL), auto_cap(allm, CAP_ALL), length.out = 102)
p <- pheatmap::pheatmap(allm, color = pal, breaks = brk, cluster_rows = FALSE, cluster_cols = FALSE,
                        annotation_col = ann_col, annotation_row = ann_row, annotation_colors = ann_colors,
                        gaps_col = unname(gaps), gaps_row = unname(cumsum(table(cl_f))[-K_CLUSTERS]),
                        show_rownames = FALSE, show_colnames = TRUE, border_color = NA, fontsize = 9, silent = TRUE,
                        main = sprintf("All protein-coding genes (%d): %s", nrow(allm), lab_unit))
for (ext in c("png", "pdf")) {
  f <- file.path(OUT, paste0("heatmap_all_protein_coding.", ext))
  if (ext == "png") png(f, width = 8, height = 11, units = "in", res = 300) else pdf(f, width = 8, height = 11)
  grid::grid.newpage(); grid::grid.draw(p$gtable); dev.off()
}
message("\nHeatmaps and gene lists in: ", OUT)
