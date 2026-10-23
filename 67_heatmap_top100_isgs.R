#!/usr/bin/env Rscript
# Heatmap of the top TOP_N IFNB ISGs (IFNB vs Mock: padj < PADJ_CUT, shrunk log2FC > LFC_CUT, protein-coding, ranked by log2FC) showing EXPRESSION in every replicate:
# columns (left to right) Untreated 1-3 | IFNB 1-3 | IFNB+EBOV 1-3. Rows are in order of IFNB log2FC, not clustered. Colour = VST expression (log2-scale, library-size normalised), not a fold-change.
#   heatmap_top100_ISGs_expression.png / .pdf        colour = VST expression (absolute level; white = low, dark red = high; scale from the 2nd to 98th percentile)
#   heatmap_top100_ISGs_expression_rowZ.png / .pdf   same values scaled per gene (row z-score: blue = below that gene's mean, red = above), which makes each gene's
#                                                     pattern across the 9 samples comparable even when genes differ a lot in absolute level
#   heatmap_top100_ISGs_expression_mockgrey.png / .pdf   every replicate coloured against that gene's Untreated (Mock) mean: Mock = light grey/white, above Mock = red, below = blue
#   heatmap_top100_ISGs_expression.csv               gene, IFNB-vs-Mock log2FC / padj, and the VST value of every sample
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
nrep <- lengths(cols)
X <- vm[top$gene_id, unlist(cols), drop = FALSE]; rownames(X) <- top$label
colnames(X) <- paste0(rep(names(cols), nrep), " ", unlist(lapply(nrep, seq_len)))
write_csv(bind_cols(top |> select(gene_id, gene_name, log2FC_IFNB_vs_Mock = log2FoldChange, padj_IFNB_vs_Mock = padj), as.data.frame(X)), file.path(OUT, "heatmap_top100_ISGs_expression.csv"))

rowfs <- if (nrow(X) > 60) 6.5 else 9
W <- 0.5 * ncol(X) + 2.5; H <- 0.17 * nrow(X) + 2.2
save_ph <- function(mat, name, pal, brk, lab, title) {
  ph <- pheatmap(mat, color = pal, breaks = brk, cluster_rows = FALSE, cluster_cols = FALSE, border_color = NA, fontsize_row = rowfs, fontsize_col = 9, fontsize = 9,
                 angle_col = 45, cellwidth = 24, cellheight = 9.5, gaps_col = cumsum(nrep)[-length(nrep)], silent = TRUE,
                 legend_breaks = lab$at, legend_labels = lab$text, main = title)
  png(file.path(OUT, paste0(name, ".png")), width = W, height = H, units = "in", res = 300); grid.newpage(); grid.draw(ph$gtable); invisible(dev.off())
  pdf(file.path(OUT, paste0(name, ".pdf")), width = W, height = H); grid.newpage(); grid.draw(ph$gtable); invisible(dev.off())
}
# 1. absolute VST expression
lo <- as.numeric(quantile(X, 1 - CAP_Q)); hi <- as.numeric(quantile(X, CAP_Q))
save_ph(X, "heatmap_top100_ISGs_expression", colorRampPalette(c("#FFFFFF", "#FEE0D2", "#FC9272", "#DE2D26", "#67000D"))(100), seq(lo, hi, length.out = 101),
        list(at = round(c(lo, (lo + hi) / 2, hi), 1), text = sprintf("%.1f", c(lo, (lo + hi) / 2, hi))),
        sprintf("Top %d IFNB ISGs: expression (VST, log2 scale)", nrow(X)))
# 2. row z-score
Z <- t(scale(t(X))); lim <- as.numeric(quantile(abs(Z), CAP_Q))
save_ph(Z, "heatmap_top100_ISGs_expression_rowZ", colorRampPalette(c("#2B6CB0", "white", "#C53030"))(100), seq(-lim, lim, length.out = 101),
        list(at = c(-lim, 0, lim), text = c(sprintf("%.1f", -lim), "0", sprintf("%.1f", lim))),
        sprintf("Top %d IFNB ISGs: expression, row z-score", nrow(Z)))
# 3. Mock-anchored: each replicate's VST minus the gene's Untreated mean, Mock level = light grey / white
Xm <- X - rowMeans(X[, seq_len(nrep[[1]]), drop = FALSE]); limm <- as.numeric(quantile(abs(Xm), CAP_Q))
save_ph(Xm, "heatmap_top100_ISGs_expression_mockgrey", colorRampPalette(c("#2B6CB0", "#F2F2F2", "#C53030"))(100), seq(-limm, limm, length.out = 101),
        list(at = c(-limm, 0, limm), text = c(sprintf("%.1f", -limm), "0 (Mock mean)", sprintf("%.1f", limm))),
        sprintf("Top %d IFNB ISGs: each replicate vs Untreated mean (Mock = grey)", nrow(Xm)))
message("VST colour range ", round(lo, 1), " to ", round(hi, 1), "; z-score range +/- ", round(lim, 2), "\nSaved to: ", OUT)
