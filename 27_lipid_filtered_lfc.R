#!/usr/bin/env Rscript
# Lipid-related genes that respond to IFNB: padj < 0.05 AND |log2FC| > 0.7
# (about 1.6-fold) in the shrunk IFNB vs Mock results. Applied to
#   (a) all lipid-metabolism genes (script 13)  and
#   (b) the endocytic-vesicle / lysosome lipid enzymes (script 16).
# Change LFC_CUT / PADJ_CUT below to use other thresholds.

suppressPackageStartupMessages({ library(readr); library(dplyr) })
LFC_CUT <- 0.7; PADJ_CUT <- 0.05

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs")

jobs <- list(
  list(inp = "DE_v2_IFNB_vs_Mock_protein_coding_lipid_metabolism.csv",
       out = "DE_v2_IFNB_vs_Mock_lipid_genes_log2FC_gt0.7_padj.csv", label = "all lipid-metabolism genes"),
  list(inp = "DE_v2_IFNB_vs_Mock_lipid_enzymes_endocytic_vesicle.csv",
       out = "DE_v2_IFNB_vs_Mock_lipid_endocytic_enzymes_log2FC_gt0.7_padj.csv", label = "endocytic-vesicle lipid enzymes"))

for (j in jobs) {
  f <- file.path(OUT, j$inp)
  if (!file.exists(f)) { message("Skipping ", j$label, " (", j$inp, " not found - run the earlier script)"); next }
  d <- read_csv(f, show_col_types = FALSE)
  keep <- d |> filter(!is.na(padj), padj < PADJ_CUT, abs(log2FoldChange) > LFC_CUT) |>
    mutate(direction = ifelse(log2FoldChange > 0, "up", "down")) |>
    relocate(direction, .after = padj) |>
    arrange(desc(log2FoldChange))
  write_csv(keep, file.path(OUT, j$out))
  message(sprintf("%-34s %4d genes in, %3d pass (up %d, down %d)  ->  %s", j$label, nrow(d), nrow(keep),
                  sum(keep$direction == "up"), sum(keep$direction == "down"), j$out))
}
