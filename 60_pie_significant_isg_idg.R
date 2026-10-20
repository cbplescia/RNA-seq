#!/usr/bin/env Rscript
# Pie / donut charts for IFNB vs Mock, all tested protein-coding genes:
#   which are significantly altered (padj < PADJ_CUT), and of those, which are altered enough to be ISGs / IDGs (shrunk |log2FC| > LFC_CUT).
#   pie_IFNB_significant_genes            significantly altered vs not significant
#   donut_IFNB_significant_ISG_IDG        inner ring: not significant / significantly altered; outer ring splits the significant genes into
#                                         ISGs, IDGs and smaller changes up / down (|log2FC| <= LFC_CUT)
# Uses requant_v2_outputs/DE_v2_IFNB_vs_Mock_protein_coding.csv (genes with padj = NA are not counted, as in the other figures).
# Writes to requant_v2_outputs/figures_protein_coding/ (PNG + PDF) and IFNB_significance_pie_counts.csv
# Usage: Rscript 60_pie_significant_isg_idg.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter

PADJ_CUT <- 0.05; LFC_CUT <- 1
INK <- "#1F2933"; C_ISG <- "#C53030"; C_IDG <- "#2B6CB0"; C_NS <- "#D9D9D9"; C_SIG <- "#6B7280"; C_UP <- "#F0A8A8"; C_DN <- "#A9C5E6"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs", "figures_protein_coding"); dir.create(OUT, showWarnings = FALSE)

d <- read_csv(file.path(PROJ, "requant_v2_outputs", "DE_v2_IFNB_vs_Mock_protein_coding.csv"), show_col_types = FALSE) |> filter(!is.na(padj))
N <- nrow(d)
sig <- d$padj < PADJ_CUT
nNS <- sum(!sig); nSig <- sum(sig)
nISG <- sum(sig & d$log2FoldChange >  LFC_CUT); nIDG <- sum(sig & d$log2FoldChange < -LFC_CUT)
nUp <- sum(sig & d$log2FoldChange > 0 & d$log2FoldChange <=  LFC_CUT); nDn <- sum(sig & d$log2FoldChange <= 0 & d$log2FoldChange >= -LFC_CUT)
message(sprintf("Protein-coding genes tested: %d | significant (padj < %g): %d (%.1f%%) | ISGs %d | IDGs %d | smaller up %d | smaller down %d | not significant %d",
                N, PADJ_CUT, nSig, 100 * nSig / N, nISG, nIDG, nUp, nDn, nNS))
stopifnot(nISG + nIDG + nUp + nDn == nSig)
write_csv(tibble(group = c("All tested protein-coding genes", "Not significant", "Significantly altered", "ISGs", "IDGs", "Significant, smaller increase", "Significant, smaller decrease"),
                 genes = c(N, nNS, nSig, nISG, nIDG, nUp, nDn),
                 percent_of_all = round(100 * c(N, nNS, nSig, nISG, nIDG, nUp, nDn) / N, 1),
                 percent_of_significant = c(NA, NA, 100, round(100 * c(nISG, nIDG, nUp, nDn) / nSig, 1))),
          file.path(OUT, "IFNB_significance_pie_counts.csv"))

ring <- function(label, n, fill, txt, xmin, xmax) tibble(label = label, n = n, fill = fill, txt = txt, xmin = xmin, xmax = xmax) |> mutate(ymax = cumsum(n), ymin = ymax - n, pct = 100 * n / N)
draw <- function(rings, title, subtitle, min_label, file, hole = 0) {
  dd <- bind_rows(rings)
  p <- ggplot(dd) + geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), colour = "white", linewidth = 1) +
    geom_text(data = dd |> filter(pct >= min_label), aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = sprintf("%s\n%d (%.1f%%)", label, n, pct), colour = txt), size = 3.2, fontface = "bold", lineheight = 0.9) +
    scale_fill_identity() + scale_colour_identity() + coord_polar(theta = "y") + scale_x_continuous(limits = c(hole, max(dd$xmax))) +
    labs(title = title, subtitle = subtitle) + theme_void(base_size = 12) +
    theme(plot.title = element_text(colour = INK, face = "bold", hjust = 0.5), plot.subtitle = element_text(colour = "grey35", size = 8.5, hjust = 0.5), plot.margin = margin(10, 10, 10, 10))
  ggsave(file.path(OUT, paste0(file, ".png")), p, width = 7.5, height = 7.5, dpi = 300)
  ggsave(file.path(OUT, paste0(file, ".pdf")), p, width = 7.5, height = 7.5)
}
# 1. simple pie
draw(list(ring(c("Significantly\naltered", "Not\nsignificant"), c(nSig, nNS), c(C_SIG, C_NS), c("white", INK), 0, 1)),
     sprintf("IFNB vs Mock: all protein-coding genes tested (n = %d)", N),
     sprintf("significantly altered = padj < %g", PADJ_CUT), 0, "pie_IFNB_significant_genes")
# 2. nested donut: significant genes split into ISG / IDG / smaller changes
inner <- ring(c("Significantly\naltered", "Not\nsignificant"), c(nSig, nNS), c(C_SIG, C_NS), c("white", INK), 1, 2)
outer <- ring(c("ISGs", "IDGs", "Smaller\nincrease", "Smaller\ndecrease", "Not significant"), c(nISG, nIDG, nUp, nDn, nNS), c(C_ISG, C_IDG, C_UP, C_DN, "#EFEFEF"), c("white", "white", INK, INK, INK), 2.05, 3.05)
draw(list(inner, outer), sprintf("IFNB vs Mock: significantly altered genes and ISGs / IDGs (n = %d protein-coding genes)", N),
     sprintf("significant: padj < %g; ISG: log2FC > %g; IDG: log2FC < -%g; smaller increase / decrease: significant with |log2FC| <= %g", PADJ_CUT, LFC_CUT, LFC_CUT, LFC_CUT),
     1, "donut_IFNB_significant_ISG_IDG")
message("Saved to: ", OUT)
