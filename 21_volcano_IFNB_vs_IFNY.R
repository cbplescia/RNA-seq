#!/usr/bin/env Rscript
# Volcano plot of the direct IFNB vs IFNY comparison (protein-coding genes), from
# 09_ifnb_vs_ifny.R. log2FC > 0 = higher in IFNB; log2FC < 0 = higher in IFNY.
# Coloured when padj < 0.05 and |log2FC| > 1. The top 10 genes higher in IFNB and
# the top 10 higher in IFNY (largest |log2FC| among the significant ones) are labelled.
# Needs 09_ifnb_vs_ifny.R to have been run, and the ggrepel package.

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })
if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Install ggrepel (conda install -c conda-forge r-ggrepel)")

LFC_CUT <- 1; PADJ_CUT <- 0.05; N_LABEL <- 10
INK <- "#1F2933"
COL <- c(`Higher in IFNB` = "#2B6CB0", `Higher in IFNγ` = "#D69E2E", `Not significant` = "#B8B8B8")

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs")
FIG <- file.path(OUT, "figures_protein_coding")
dir.create(FIG, showWarnings = FALSE)

d <- read_csv(file.path(OUT, "DE_v2_IFNB_vs_IFNY_protein_coding.csv"), show_col_types = FALSE) |>
  transmute(gene_name, lfc = log2FC_IFNB_vs_IFNY, padj = padj_diff) |> filter(!is.na(padj)) |>
  mutate(y = -log10(pmax(padj, 1e-300)),
         status = factor(case_when(padj < PADJ_CUT & lfc >  LFC_CUT ~ "Higher in IFNB",
                                   padj < PADJ_CUT & lfc < -LFC_CUT ~ "Higher in IFNγ",
                                   TRUE ~ "Not significant"), levels = names(COL)))
lab <- bind_rows(d |> filter(status == "Higher in IFNB") |> slice_max(lfc, n = N_LABEL),
                 d |> filter(status == "Higher in IFNγ") |> slice_min(lfc, n = N_LABEL))
message("Genes tested: ", nrow(d))
print(table(d$status))
message("Top ", N_LABEL, " higher in IFNB: ", paste(lab$gene_name[lab$status == "Higher in IFNB"], collapse = ", "))
message("Top ", N_LABEL, " higher in IFNγ: ", paste(lab$gene_name[lab$status == "Higher in IFNγ"], collapse = ", "))

xr <- max(abs(d$lfc)) * 1.05
p <- ggplot(d, aes(lfc, y)) +
  geom_hline(yintercept = -log10(PADJ_CUT), linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_point(aes(colour = status), size = 0.9, alpha = 0.6) +
  geom_point(data = lab, aes(colour = status), size = 1.6) +
  ggrepel::geom_text_repel(data = lab, aes(label = gene_name), colour = INK, size = 3,
                           max.overlaps = Inf, min.segment.length = 0, segment.colour = "grey60",
                           segment.size = 0.25, box.padding = 0.4, seed = 1) +
  scale_colour_manual(values = COL, name = NULL) +
  scale_x_continuous(limits = c(-xr, xr)) +
  labs(title = "IFNB vs IFNγ: differential expression of protein-coding genes",
       subtitle = sprintf("shrunk log2FC; coloured when padj<%g & |log2FC|>%g; top %d genes higher in each interferon labelled",
                          PADJ_CUT, LFC_CUT, N_LABEL),
       x = "log2 fold-change (IFNB vs IFNγ)   ← higher in IFNγ  |  higher in IFNB →",
       y = expression(-log[10](padj))) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.border = element_rect(colour = "grey70"),
        axis.text = element_text(colour = INK), axis.title = element_text(colour = INK),
        plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35"),
        legend.position = "bottom")
ggsave(file.path(FIG, "volcano_IFNB_vs_IFNY.png"), p, width = 8, height = 7, dpi = 300)
ggsave(file.path(FIG, "volcano_IFNB_vs_IFNY.pdf"), p, width = 8, height = 7)
message("\nSaved to: ", FIG)
