#!/usr/bin/env Rscript
# Clean version of the whole protein-coding transcriptome heatmap for Mock | IFNB | IFNB+EBOV | EBOV.
# Same gene ORDER as the organised heatmap from 51_mock_ifnb_ebov_analysis.R (IFNB-up genes first, then IFNB-down, then genes
# only EBOV changes, then unchanged genes; inside the IFNB groups EBOV-suppressed genes come before partly- and non-suppressed
# ones), but with NO row annotation bars and NO row labels, so the patterns can be read from the layout alone.
# Only the sample-group colour bar above the columns is kept. Thin white gaps separate the blocks (GAPS <- FALSE removes them).
# Reads figures_mock_ifnb_ebov/transcriptome_blocks.csv (written by 51), so run 51 first. Nothing is recomputed.
# Colour = log2 change of each sample's VST vs the mean of the Mock samples (blue lower, red higher).
# Usage: Rscript 52_heatmap_all_clean.R   (PROJ from 'source setup_session.sh')
# Writes figures_mock_ifnb_ebov/heatmap_all_protein_coding_clean.png / .pdf

suppressPackageStartupMessages({ library(DESeq2); library(readr); library(dplyr) })
if (!requireNamespace("pheatmap", quietly = TRUE)) stop("Install pheatmap (conda install -c conda-forge r-pheatmap)")
select <- dplyr::select; filter <- dplyr::filter

GAPS <- TRUE; CAP <- NULL           # CAP: fix the colour limit (log2 units); NULL = automatic
SETS <- list(Mock = c("1", "2", "3"), IFNB = c("10", "11", "12"), `IFNB+EBOV` = c("16", "17", "18"), EBOV = c("4", "5", "6"))
GCOL <- c(Mock = "#7F7F7F", IFNB = "#2B6CB0", `IFNB+EBOV` = "#9BD08E", EBOV = "#1B7837")
BLOCKS <- tibble::tribble(~resp, ~eff, ~title,
  "Up", "Suppressed", "IFNB-up, suppressed by EBOV", "Up", "Partly suppressed", "IFNB-up, partly suppressed", "Up", "Not suppressed", "IFNB-up, NOT suppressed",
  "Down", "Suppressed", "IFNB-down, reversed by EBOV", "Down", "Partly suppressed", "IFNB-down, partly reversed", "Down", "Not suppressed", "IFNB-down, NOT reversed",
  "None", "EBOV-altered (no IFNB response)", "EBOV-altered only", "None", "Unchanged", "Unchanged")

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs"); OUT <- file.path(RNA, "figures_mock_ifnb_ebov")
bf <- file.path(OUT, "transcriptome_blocks.csv")
if (!file.exists(bf)) stop("Run 51_mock_ifnb_ebov_analysis.R first (missing ", bf, ")")
B <- read_csv(bf, show_col_types = FALSE)

vsd <- readRDS(file.path(RNA, "vsd_v2.rds")); vm <- assay(vsd)
cols <- unlist(SETS)
if (!all(cols %in% colnames(vm))) stop("Samples missing from the VST object: ", paste(setdiff(cols, colnames(vm)), collapse = ", "))
mock_cols <- SETS$Mock

# same ordering as 51: block by block; IFNB-down by most negative first, EBOV-altered by EBOV log2FC, the rest by IFNB log2FC
ordered <- bind_rows(lapply(seq_len(nrow(BLOCKS)), function(i) {
  b <- B |> filter(block == BLOCKS$title[i])
  b[order(if (BLOCKS$resp[i] == "Down") b$lfc_B else if (BLOCKS$resp[i] == "None" && BLOCKS$eff[i] != "Unchanged") -b$lfc_E else -b$lfc_B), ]
})) |> filter(gene_id %in% rownames(vm))
blk <- factor(ordered$block, levels = BLOCKS$title)
message("Genes: ", nrow(ordered)); print(table(blk))

m <- vm[ordered$gene_id, cols, drop = FALSE]
m <- m - rowMeans(m[, mock_cols, drop = FALSE])
cap <- if (!is.null(CAP)) CAP else min(8, max(2, ceiling(quantile(abs(m), 0.98, na.rm = TRUE))))
grp <- rep(names(SETS), lengths(SETS))
ann_col <- data.frame(Group = factor(grp, levels = names(SETS)), row.names = cols)
cnt <- as.integer(table(blk))
gaps_row <- if (GAPS) cumsum(cnt)[-length(cnt)][cnt[-length(cnt)] > 0] else NULL

h <- pheatmap::pheatmap(m, color = colorRampPalette(c("#2C5F9E", "#F7F7F7", "#B2382C"))(101), breaks = seq(-cap, cap, length.out = 102),
                        cluster_rows = FALSE, cluster_cols = FALSE, annotation_col = ann_col, annotation_colors = list(Group = GCOL),
                        gaps_col = unname(cumsum(lengths(SETS))[-length(SETS)]), gaps_row = gaps_row,
                        show_rownames = FALSE, show_colnames = TRUE, border_color = NA, fontsize = 9, silent = TRUE,
                        main = sprintf("All protein-coding genes (%d): log2 change vs Mock mean (VST)", nrow(m)))
for (ext in c("png", "pdf")) {
  f <- file.path(OUT, paste0("heatmap_all_protein_coding_clean.", ext))
  if (ext == "png") png(f, width = 7, height = 11, units = "in", res = 300) else pdf(f, width = 7, height = 11)
  grid::grid.newpage(); grid::grid.draw(h$gtable); dev.off()
}
message("Saved heatmap_all_protein_coding_clean.png / .pdf in ", OUT)
