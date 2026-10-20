#!/usr/bin/env Rscript
# The "significantly altered vs not significant" pie for IFNB vs Mock (all tested protein-coding genes) with a zoom-in inset on the
# significantly altered wedge: the altered genes split into ISGs (red), IDGs (blue) and smaller increases / decreases (|log2FC| <= LFC_CUT).
# Dashed lines connect the inset to the wedge it magnifies. Same numbers as 60_pie_significant_isg_idg.R.
# Writes pie_IFNB_significant_with_ISG_IDG_inset.png / .pdf (11 x 7 in) to requant_v2_outputs/figures_protein_coding/
# Usage: Rscript 61_pie_significant_with_isg_idg_inset.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2); library(grid) })
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
N <- nrow(d); sig <- d$padj < PADJ_CUT; nSig <- sum(sig); nNS <- N - nSig
nISG <- sum(sig & d$log2FoldChange > LFC_CUT); nIDG <- sum(sig & d$log2FoldChange < -LFC_CUT)
nUp <- sum(sig & d$log2FoldChange > 0 & d$log2FoldChange <= LFC_CUT); nDn <- nSig - nISG - nIDG - nUp
message(sprintf("tested %d | significant %d | ISGs %d | IDGs %d | smaller up %d | smaller down %d | not significant %d", N, nSig, nISG, nIDG, nUp, nDn, nNS))

pie_theme <- theme_void() + theme(plot.margin = margin(0, 0, 0, 0), legend.position = "none")
build <- function(df, total, xmax, text_size) {
  df <- df |> mutate(ymax = cumsum(n), ymin = ymax - n)
  ggplot(df) + geom_rect(aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax, fill = fill), colour = "white", linewidth = 1) +
    geom_text(aes(x = lx, y = (ymin + ymax) / 2, label = lab, colour = txt), size = text_size, fontface = "bold", lineheight = 0.9) +
    scale_fill_identity() + scale_colour_identity() + coord_polar(theta = "y") +
    scale_x_continuous(limits = c(0, xmax), expand = c(0, 0)) + scale_y_continuous(limits = c(0, total), expand = c(0, 0)) + pie_theme
}
p_main <- build(tibble(n = c(nSig, nNS), fill = c(C_SIG, C_NS), txt = c("white", INK), lx = 0.6,
                       lab = c(sprintf("Significantly\naltered\n%d (%.1f%%)", nSig, 100 * nSig / N), sprintf("Not\nsignificant\n%d (%.1f%%)", nNS, 100 * nNS / N))), N, 1, 4.6)
pct <- function(x) sprintf("%d (%.1f%%)", x, 100 * x / nSig)
p_in <- build(tibble(n = c(nISG, nIDG, nUp, nDn), fill = c(C_ISG, C_IDG, C_UP, C_DN), txt = c(C_ISG, C_IDG, INK, INK), lx = c(1.3, 1.3, 0.55, 0.55),
                     lab = c(paste0("ISGs\n", pct(nISG)), paste0("IDGs\n", pct(nIDG)), paste0("Smaller\nincrease\n", pct(nUp)), paste0("Smaller\ndecrease\n", pct(nDn)))), nSig, 1.55, 3.4)

# layout in inches on an 11 x 7 canvas (a coord_polar pie fills 0.4 of its square panel, which fixes where the wedge edges are)
W <- 11; H <- 7
cx <- 3.3; cy <- 3.4; S <- 6.2; r <- 0.4 * S
cx2 <- 8.5; cy2 <- 3.4; S2 <- 4.4; r2 <- 0.4 / 1.55 * S2
th <- 2 * pi * nSig / N                                  # the altered wedge runs clockwise from 12 o'clock to angle th
u <- function(x) unit(x, "in")
draw <- function() {
  grid.newpage()
  grid.text("IFNB vs Mock: all protein-coding genes tested", x = u(0.4), y = u(H - 0.35), just = "left", gp = gpar(fontface = "bold", fontsize = 17, col = INK))
  grid.text(sprintf("n = %d; significantly altered = padj < %g; inset: ISG = log2FC > %g, IDG = log2FC < -%g (shrunk log2FC)", N, PADJ_CUT, LFC_CUT, LFC_CUT),
            x = u(0.4), y = u(H - 0.75), just = "left", gp = gpar(fontsize = 9.5, col = "grey35"))
  print(p_main, vp = viewport(x = u(cx), y = u(cy), width = u(S), height = u(S)))
  grid.lines(x = u(c(cx, cx2)), y = u(c(cy + r, cy2 + r2)), gp = gpar(col = "grey55", lty = 2, lwd = 1))                                   # top edge of the wedge -> top of inset
  grid.lines(x = u(c(cx + r * sin(th), cx2)), y = u(c(cy + r * cos(th), cy2 - r2)), gp = gpar(col = "grey55", lty = 2, lwd = 1))             # other edge -> bottom of inset
  print(p_in, vp = viewport(x = u(cx2), y = u(cy2), width = u(S2), height = u(S2)))
  grid.text(sprintf("Significantly altered genes (n = %d)", nSig), x = u(cx2), y = u(cy2 + r2 + 1.05), gp = gpar(fontface = "bold", fontsize = 11, col = C_SIG))
}
png(file.path(OUT, "pie_IFNB_significant_with_ISG_IDG_inset.png"), width = W, height = H, units = "in", res = 300); draw(); invisible(dev.off())
pdf(file.path(OUT, "pie_IFNB_significant_with_ISG_IDG_inset.pdf"), width = W, height = H); draw(); invisible(dev.off())
message("Saved to: ", OUT)
