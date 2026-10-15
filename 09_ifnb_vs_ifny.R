#!/usr/bin/env Rscript
# Which protein-coding genes respond differently to IFNB than to IFNY?
#
# Tests IFNB_Mock vs IFNY_Mock directly in the saved DESeq2 model (Mock cancels
# out), restricted to protein-coding genes, with ashr-shrunk log2FC. Positive
# log2FC = higher in IFNB than IFNY. Then lines each gene up against its
# IFNB-vs-Mock and IFNY-vs-Mock results and labels the response pattern.
# Needs 07_add_gene_symbols.R and 08_split_protein_coding.R to have been run.

suppressPackageStartupMessages({
  library(DESeq2); library(readr); library(dplyr); library(ggplot2)
})

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs")

LFC_CUT  <- 1      # same thresholds as 05_compare_old_vs_new.R
PADJ_CUT <- 0.05

rd <- function(nm) read_csv(file.path(OUT, paste0("DE_v2_", nm, "_protein_coding.csv")),
                            show_col_types = FALSE)
b <- rd("IFNB_vs_Mock"); y <- rd("IFNY_vs_Mock")
pc_ids <- intersect(b$gene_id, y$gene_id)

dds <- readRDS(file.path(OUT, "dds_v2_decoyAware_hostNormalised.rds"))
dds <- dds[rownames(dds) %in% pc_ids, ]        # protein-coding only, so padj is
message("Protein-coding genes tested: ", nrow(dds))   # corrected over this set

ct  <- c("group", "IFNB_Mock", "IFNY_Mock")
res <- results(dds, contrast = ct)
shr <- if (requireNamespace("ashr", quietly = TRUE)) {
  lfcShrink(dds, contrast = ct, res = res, type = "ashr")
} else {
  lfcShrink(dds, contrast = ct, type = "normal")
}

diff <- as.data.frame(shr) |> tibble::rownames_to_column("gene_id") |>
  transmute(gene_id, baseMean,
            log2FC_IFNB_vs_IFNY = log2FoldChange,
            pvalue_diff = res[gene_id, "pvalue"],
            padj_diff   = res[gene_id, "padj"])

sig <- function(l, p) !is.na(p) & p < PADJ_CUT & abs(l) > LFC_CUT
dir <- function(l) ifelse(l > 0, "up", "down")

tab <- diff |>
  left_join(b |> select(gene_id, gene_name, gene_biotype,
                        log2FC_IFNB_vs_Mock = log2FoldChange, padj_IFNB_vs_Mock = padj),
            by = "gene_id") |>
  left_join(y |> select(gene_id, log2FC_IFNY_vs_Mock = log2FoldChange,
                        padj_IFNY_vs_Mock = padj), by = "gene_id") |>
  mutate(
    sigB = sig(log2FC_IFNB_vs_Mock, padj_IFNB_vs_Mock),
    sigY = sig(log2FC_IFNY_vs_Mock, padj_IFNY_vs_Mock),
    response = case_when(
      sigB & sigY & sign(log2FC_IFNB_vs_Mock) == sign(log2FC_IFNY_vs_Mock) ~
        paste0("both_", dir(log2FC_IFNB_vs_Mock)),
      sigB & sigY ~ "opposite_direction",
      sigB ~ paste0("IFNB_only_", dir(log2FC_IFNB_vs_Mock)),
      sigY ~ paste0("IFNY_only_", dir(log2FC_IFNY_vs_Mock)),
      TRUE ~ "neither"),
    differs_sig = sig(log2FC_IFNB_vs_IFNY, padj_diff),
    higher_in = ifelse(differs_sig, ifelse(log2FC_IFNB_vs_IFNY > 0, "IFNB", "IFNY"), NA)
  ) |>
  select(gene_id, gene_name, gene_biotype, baseMean,
         log2FC_IFNB_vs_IFNY, padj_diff, differs_sig, higher_in,
         log2FC_IFNB_vs_Mock, padj_IFNB_vs_Mock,
         log2FC_IFNY_vs_Mock, padj_IFNY_vs_Mock, response) |>
  arrange(padj_diff, desc(abs(log2FC_IFNB_vs_IFNY)))

write_csv(tab, file.path(OUT, "DE_v2_IFNB_vs_IFNY_protein_coding.csv"))

message("\n--- genes that differ between IFNB and IFNY (padj<", PADJ_CUT,
        ", |log2FC|>", LFC_CUT, ") ---")
print(as.data.frame(tab |> filter(differs_sig) |> count(higher_in)), row.names = FALSE)
message("\n--- response pattern vs Mock ---")
print(as.data.frame(tab |> count(response, sort = TRUE)), row.names = FALSE)

p <- ggplot(tab, aes(log2FC_IFNY_vs_Mock, log2FC_IFNB_vs_Mock)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey70", linewidth = 0.4) +
  geom_point(aes(colour = differs_sig), alpha = 0.5, size = 0.7) +
  scale_colour_manual(values = c(`FALSE` = "grey75", `TRUE` = "#C53030"),
                      labels = c("not different", "IFNB vs IFNY significant")) +
  labs(title = "Protein-coding genes: IFNB vs IFNY response",
       x = "log2FC, IFNY vs Mock", y = "log2FC, IFNB vs Mock", colour = NULL) +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
ggsave(file.path(OUT, "IFNB_vs_IFNY_scatter.png"), p, width = 6, height = 5.5, dpi = 300)
message("\nOutputs in: ", OUT)
