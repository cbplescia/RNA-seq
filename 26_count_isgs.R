#!/usr/bin/env Rscript
# How many ISGs (and IDGs) does IFNB induce? Counts genes in the shrunk
# IFNB vs Mock results at a few log2FC cutoffs (padj < 0.05), for all genes
# and for protein-coding genes only. Default definition used in the figures:
# padj < 0.05 and log2FC > 1.
suppressPackageStartupMessages({ library(readr); library(dplyr) })
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs")
all_g <- read_csv(file.path(OUT, "DE_v2_IFNB_vs_Mock_annotated.csv"), show_col_types = FALSE)
pc_g  <- all_g |> filter(gene_biotype == "protein_coding")
count <- function(d, label) bind_rows(lapply(c(0.585, 1, 2), function(l) tibble(
  gene_set = label, log2FC_cutoff = l, fold_change = round(2^l, 2),
  ISGs_up = sum(d$padj < 0.05 & d$log2FoldChange >  l, na.rm = TRUE),
  IDGs_down = sum(d$padj < 0.05 & d$log2FoldChange < -l, na.rm = TRUE))))
res <- bind_rows(count(all_g, "all genes tested"), count(pc_g, "protein-coding only"))
print(as.data.frame(res), row.names = FALSE)
message("\nDefault definition (padj<0.05, log2FC>1): ",
        sum(all_g$padj < 0.05 & all_g$log2FoldChange > 1, na.rm = TRUE), " ISGs among ", nrow(all_g),
        " genes tested; ", sum(pc_g$padj < 0.05 & pc_g$log2FoldChange > 1, na.rm = TRUE),
        " among ", nrow(pc_g), " protein-coding genes.")
