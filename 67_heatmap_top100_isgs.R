#!/usr/bin/env Rscript
# Heatmap of the top TOP_N IFNB ISGs (IFNB vs Mock: padj < PADJ_CUT, shrunk log2FC > LFC_CUT, protein-coding, ranked by log2FC) across
# Untreated (Mock) | IFNB | IFNB+EBOV. Values = VST expression minus the gene's mean in the Untreated samples (log2 change vs untreated), so the
# Untreated column sits near 0; colour scale is symmetric blue-white-red and capped at the 98th percentile. Rows are in order of IFNB log2FC, not clustered.
#   heatmap_top100_ISGs_groups.png / .pdf          three columns (group means), numbers in cells
#   heatmap_top100_ISGs_replicates.png / .pdf      the same genes, every replicate (3 + 3 + 3 columns)
#   heatmap_top100_ISGs.csv                        gene, log2FC / padj IFNB vs Mock, and the group means shown
# Writes to requant_v2_outputs/figures_heatmaps/   Usage: Rscript 67_heatmap_top100_isgs.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(pheatmap); library(DESeq2); library(grid) })
select <- dplyr::select; filter <- dplyr::filter

TOP_N <- 100; PADJ_CUT <- 0.05; LFC_CUT <- 1; CAP_Q <- 0.98
G <- c(Untreated = "Mock_Mock", IFNB = "IFNB_Mock", `IFNB+EBOV` = "IFNB_EBOV")      # column order, left to right
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_heatmaps"); dir.create(OUT, showWarnings = FALSE)

pcg <- read_csv(file.path(RNA, "DE_v2_IFNB_vs_Mock_protein_coding.csv"), show_col_types = FALSE) |> filter(!is.na(padj))
isg <- pcg |> filter(padj < PADJ_CUT, log2FoldChange > LFC_CUT) |> arrange(desc(log2FoldChange))
message(sprintf("ISGs: %d; using the top %d by log2FC", nrow(isg), min(TOP_N, nrow(isg))))
top <- isg |> slice_head(n = TOP_N) |> mutate(label = make.unique(ifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name)))

vsd <- readRDS(file.path(RNA, "vsd_v2.rds")); vm <- SummarizedExperiment::assay(vsd)
grp <- as.character(vsd$group)
cols <- lapply(G, function(g) colnames(vm)[grp == g])
if (any(lengths(cols) == 0)) stop("Missing groups in the VST object: ", paste(names(G)[lengths(cols) == 0], collapse = ", "))
message("Samples used: ", paste(sprintf("%s = %s", names(cols), vapply(cols, paste, "", collapse = ",")), collapse = " | "))
X <- vm[top$gene_id, unlist(cols), drop = FALSE]
X <- X - rowMeans(X[, cols[[1]], drop = FALSE])                       # minus the Untreated mean
rownames(X) <- top$label
Mg <- sapply(cols, function(cn) rowMeans(X[, cn, drop = FALSE]))      # group means
write_csv(bind_cols(top |> select(gene_id, gene_name, log2FC_IFNB_vs_Mock = log2FoldChange, padj_IFNB_vs_Mock = padj), as.data.frame(Mg) |> setNames(paste0("mean_vs_untreated_", sub("\\+", "_", names(cols))))),
          file.path(OUT, "heatmap_top100_ISGs.csv"))

lim <- as.numeric(quantile(abs(X), CAP_Q)); brk <- seq(-lim, lim, length.out = 101)
pal <- colorRampPalette(c("#2B6CB0", "white", "#C53030"))(100)
rowfs <- if (nrow(X) > 60) 6.5 else 9
save_ph <- function(mat, name, w, h, ...) {
  ph <- pheatmap(mat, color = pal, breaks = brk, cluster_rows = FALSE, cluster_cols = FALSE, border_color = NA, fontsize_row = rowfs, fontsize_col = 10, fontsize = 9,
                 angle_col = 0, silent = TRUE, legend_breaks = c(-lim, 0, lim), legend_labels = c(sprintf("%.1f", -lim), "0", sprintf("%.1f", lim)),
                 main = sprintf("Top %d IFNB ISGs: log2 change vs Untreated", nrow(mat)), ...)
  png(file.path(OUT, paste0(name, ".png")), width = w, height = h, units = "in", res = 300); grid.newpage(); grid.draw(ph$gtable); invisible(dev.off())
  pdf(file.path(OUT, paste0(name, ".pdf")), width = w, height = h); grid.newpage(); grid.draw(ph$gtable); invisible(dev.off())
}
save_ph(Mg, "heatmap_top100_ISGs_groups", 4.5, 0.17 * nrow(Mg) + 2, cellwidth = 55, cellheight = 9.5, display_numbers = matrix(sprintf("%.1f", Mg), nrow(Mg)), number_color = "grey20", fontsize_number = 4.5)
nrep <- lengths(cols)
save_ph(X, "heatmap_top100_ISGs_replicates", 0.55 * ncol(X) + 2.5, 0.17 * nrow(X) + 2, cellwidth = 24, cellheight = 9.5, gaps_col = cumsum(nrep)[-length(nrep)],
        labels_col = paste0(rep(names(cols), nrep), " ", unlist(lapply(nrep, seq_len))), angle_col = 45)
message("Colour scale: +/- ", round(lim, 2), " log2 units (", CAP_Q * 100, "th percentile of |change|)\nSaved to: ", OUT)
