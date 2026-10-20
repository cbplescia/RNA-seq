#!/usr/bin/env Rscript
# One pie of all tested protein-coding genes (IFNB vs Mock). The "significantly altered" wedge (padj < PADJ_CUT) is drawn to scale as thin
# slices, ONE PER GENE, ordered by log2FC (largest decrease first, largest increase last) and coloured by that gene's log2FC:
#   navy -> blue = IDGs (log2FC < -LFC_CUT) | pale middle band = changes that did not meet the FC threshold | red -> dark red = ISGs (log2FC > LFC_CUT)
# So the width of each colour band is the number of genes and the colour is their fold-change. ISG and IDG ends are labelled with bracket
# ticks; the grey wedge is the genes that are not significant. Shrunk log2FC from DE_v2_IFNB_vs_Mock_protein_coding.csv.
# Writes pie_IFNB_significant_log2FC_gradient.png / .pdf to requant_v2_outputs/figures_protein_coding/
# Usage: Rscript 62_pie_significant_log2fc_gradient.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter; arrange <- dplyr::arrange

PADJ_CUT <- 0.05; LFC_CUT <- 1
INK <- "#1F2933"; C_NS <- "#D9D9D9"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs", "figures_protein_coding"); dir.create(OUT, showWarnings = FALSE)

d <- read_csv(file.path(PROJ, "requant_v2_outputs", "DE_v2_IFNB_vs_Mock_protein_coding.csv"), show_col_types = FALSE) |> filter(!is.na(padj))
N <- nrow(d)
S <- d |> filter(padj < PADJ_CUT) |> arrange(log2FoldChange) |> mutate(i = row_number())
nSig <- nrow(S); nNS <- N - nSig
nISG <- sum(S$log2FoldChange > LFC_CUT); nIDG <- sum(S$log2FoldChange < -LFC_CUT); nOth <- nSig - nISG - nIDG
message(sprintf("tested %d | significant %d | ISGs %d | IDGs %d | significant but below FC threshold %d | not significant %d", N, nSig, nISG, nIDG, nOth, nNS))

lo <- min(S$log2FoldChange); hi <- max(S$log2FoldChange)
cols <- c("#08306B", "#2C7FB8", "#9ECAE1", "#FF1A2B", "#8B0000")                       # navy, blue (IDG) | pale blue (middle) | red, dark red (ISG)
vals <- scales::rescale(c(lo, -LFC_CUT, 0, LFC_CUT, hi), from = c(lo, hi))              # pale blue sits at 0, blue / red at the FC thresholds
mid <- function(a, b) (a + b) / 2
p <- ggplot() +
  geom_rect(data = S, aes(xmin = 0, xmax = 1, ymin = i - 1, ymax = i, fill = log2FoldChange, colour = log2FoldChange), linewidth = 0.1) +   # colour = fill hides seams
  geom_rect(aes(xmin = 0, xmax = 1, ymin = nSig, ymax = N), fill = C_NS, colour = "white", linewidth = 0.8) +
  annotate("text", x = 0.55, y = nSig + nNS / 2, label = sprintf("Not significant\n%d (%.1f%%)", nNS, 100 * nNS / N), colour = INK, fontface = "bold", size = 4.4, lineheight = 0.95) +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = nSig - nISG, ymax = nSig), fill = NA, colour = "black", linewidth = 0.5) +                    # ISG sector border
  geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = nIDG), fill = NA, colour = "black", linewidth = 0.5) +        # IDG sector border
  # bracket around the significantly altered wedge + labels
  geom_path(data = tibble::tibble(x = c(0.97, 1.04, 1.04, 0.97), y = c(0, 0, nSig, nSig)), aes(x, y), colour = "black", linewidth = 0.5, lineend = "round") +   # thin outline round the periphery, ends reach back toward the pie
  annotate("segment", x = 1.04, xend = 1.3, y = nSig / 2, yend = nSig / 2, colour = INK, linewidth = 0.4) +
  annotate("text", x = 1.58, y = nSig / 2, label = sprintf("Significantly altered\n%d (%.1f%%)", nSig, 100 * nSig / N), colour = INK, fontface = "bold", size = 4.2, lineheight = 0.95) +
  annotate("segment", x = 1.0, xend = 1.22, y = nSig - nISG / 2, yend = nSig - nISG / 2, colour = "#8B0000", linewidth = 0.4) +
  annotate("text", x = 1.38, y = nSig - nISG / 2, label = sprintf("ISGs (log2FC > %g)\n%d", LFC_CUT, nISG), colour = "#8B0000", fontface = "bold", size = 3.8, lineheight = 0.95) +
  annotate("segment", x = 1.0, xend = 1.22, y = nIDG / 2, yend = nIDG / 2, colour = "#08306B", linewidth = 0.4) +
  annotate("text", x = 1.42, y = nIDG / 2, label = sprintf("IDGs (log2FC < -%g)\n%d", LFC_CUT, nIDG), colour = "#08306B", fontface = "bold", size = 3.8, lineheight = 0.95) +
  scale_fill_gradientn(colours = cols, values = vals, limits = c(lo, hi), name = "log2FC\n(IFNB vs Mock)", breaks = c(-LFC_CUT, 0, LFC_CUT, floor(hi))) +
  scale_colour_gradientn(colours = cols, values = vals, limits = c(lo, hi), guide = "none") +
  coord_polar(theta = "y") + scale_x_continuous(limits = c(0, 1.85), expand = c(0, 0)) + scale_y_continuous(limits = c(0, N), expand = c(0, 0)) +
  guides(fill = guide_colourbar(barwidth = 14, barheight = 0.9, title.position = "left", title.vjust = 0.8)) +
  labs(title = sprintf("IFNB vs Mock: all protein-coding genes tested (n = %d)", N),
       subtitle = paste(strwrap(sprintf("significantly altered = padj < %g. Inside that wedge each gene is one thin slice ordered by log2FC and coloured by it; pale middle band (|log2FC| <= %g, %d genes) = change did not meet the FC threshold",
                                   PADJ_CUT, LFC_CUT, nOth), width = 110), collapse = "\n")) +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(colour = INK, face = "bold", hjust = 0.5), plot.subtitle = element_text(colour = "grey35", size = 8, hjust = 0.5),
        legend.position = "bottom", legend.title = element_text(colour = INK, size = 9), legend.text = element_text(colour = INK, size = 9), plot.margin = margin(10, 10, 10, 10))
ggsave(file.path(OUT, "pie_IFNB_significant_log2FC_gradient.png"), p, width = 8.5, height = 7.5, dpi = 300)
ggsave(file.path(OUT, "pie_IFNB_significant_log2FC_gradient.pdf"), p, width = 8.5, height = 7.5)
message("Saved to: ", OUT)
