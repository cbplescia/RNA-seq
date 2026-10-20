#!/usr/bin/env Rscript
# Volcano plots of the IFNB response, combined from the earlier 56 and 57 scripts. Makes two figures:
#
#  1. volcano_IFNB_vs_Mock_ISG_IDG        IFNB vs Mock, protein-coding genes: red = ISGs (padj < PADJ_CUT, log2FC > LFC_CUT),
#                                         blue = IDGs (padj < PADJ_CUT, log2FC < -LFC_CUT), black = everything else (genes beyond
#                                         +/-LFC_CUT whose padj is not significant stay black).
#  2. volcano_four_panel                  the same genes, coloured by their IFNB result in EVERY panel, across four comparisons:
#                                         IFNB vs Mock | EBOV vs Mock | IFNB+EBOV vs Mock | IFNB+EBOV vs IFNB (EBOV on top of IFNB).
#                                         Each panel counts how many ISGs / IDGs are still significantly changed (bottom right:
#                                         pulled back toward Mock). Needs figures_mock_ifnb_ebov/transcriptome_blocks.csv from
#                                         51_mock_ifnb_ebov_analysis.R; without it only figure 1 is made.
# Figure sizes (inches) are set just below: the heights were reduced from 7 -> 5.5 (single) and 10 -> 7.5 (four-panel).
# Writes figure 1 to requant_v2_outputs/figures_protein_coding/ and figure 2 + counts to requant_v2_outputs/figures_mock_ifnb_ebov/ (PNG + PDF).
# Usage: Rscript 58_volcano_ifnb_and_four_panel.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter

LFC_CUT <- 1; PADJ_CUT <- 0.05
SINGLE_W <- 8;  SINGLE_H <- 5.5      # figure 1
FOUR_W   <- 12; FOUR_H   <- 7.5      # figure 2
INK <- "#1F2933"; RED <- "#D7191C"; BLUE <- "#2B6CB0"

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT1 <- file.path(RNA, "figures_protein_coding");  dir.create(OUT1, showWarnings = FALSE)
OUT2 <- file.path(RNA, "figures_mock_ifnb_ebov");  dir.create(OUT2, showWarnings = FALSE)
save2 <- function(p, dir, name, w, h) {
  ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h)
}
base_theme <- function(size = 12) theme_bw(base_size = size) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey94", linewidth = 0.3), panel.border = element_rect(colour = "grey70"),
        axis.text = element_text(colour = INK), axis.title = element_text(colour = INK), plot.title = element_text(colour = INK, face = "bold"),
        plot.subtitle = element_text(colour = "grey35", size = 8.5), strip.background = element_rect(fill = "grey95", colour = "grey70"),
        strip.text = element_text(colour = INK, face = "bold"))

# ======================= 1. single IFNB vs Mock volcano ============================================================================
d <- read_csv(file.path(RNA, "DE_v2_IFNB_vs_Mock_protein_coding.csv"), show_col_types = FALSE) |>
  filter(!is.na(padj)) |>
  mutate(y = -log10(pmax(padj, 1e-300)),
         class = case_when(padj < PADJ_CUT & log2FoldChange >  LFC_CUT ~ "ISG",
                           padj < PADJ_CUT & log2FoldChange < -LFC_CUT ~ "IDG", TRUE ~ "other"))
n_cls <- table(factor(d$class, levels = c("ISG", "IDG", "other")))
message(sprintf("Protein-coding genes: %d | ISGs (red): %d | IDGs (blue): %d | other (black): %d", nrow(d), n_cls[["ISG"]], n_cls[["IDG"]], n_cls[["other"]]))
message(sprintf("Genes beyond +/-%g log2FC but not significant (kept black): %d", LFC_CUT, sum(abs(d$log2FoldChange) > LFC_CUT & d$padj >= PADJ_CUT)))
write_csv(d |> filter(class != "other") |> arrange(class, desc(abs(log2FoldChange))) |> select(gene_id, gene_name, class, log2FoldChange, padj),
          file.path(OUT1, "volcano_IFNB_vs_Mock_ISGs_IDGs.csv"))
xr <- max(abs(d$log2FoldChange)) * 1.05
p1 <- ggplot(d, aes(log2FoldChange, y)) +
  geom_point(data = d |> filter(class == "other"), colour = "black", size = 0.7, alpha = 0.45) +
  geom_point(data = d |> filter(class != "other"), aes(colour = class), size = 1.1, alpha = 0.85) +
  geom_hline(yintercept = -log10(PADJ_CUT), linetype = "dashed", colour = "grey50", linewidth = 0.35) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", colour = "grey50", linewidth = 0.35) +
  scale_colour_manual(values = c(ISG = RED, IDG = BLUE), guide = "none") +
  scale_x_continuous(limits = c(-xr, xr)) +
  annotate("text", x = xr * 0.97, y = Inf, label = sprintf("ISGs: %d", n_cls[["ISG"]]), colour = RED, hjust = 1, vjust = 1.6, fontface = "bold", size = 4.2) +
  annotate("text", x = -xr * 0.97, y = Inf, label = sprintf("IDGs: %d", n_cls[["IDG"]]), colour = BLUE, hjust = 0, vjust = 1.6, fontface = "bold", size = 4.2) +
  labs(title = "IFNB vs Mock: protein-coding transcripts",
       subtitle = sprintf("red = ISGs (log2FC > %g, padj < %g); blue = IDGs (log2FC < -%g, padj < %g); black = all other genes; %d genes", LFC_CUT, PADJ_CUT, LFC_CUT, PADJ_CUT, nrow(d)),
       x = "log2 fold-change (IFNB vs Mock)", y = expression(-log[10](padj))) + base_theme(12) +
  theme(plot.subtitle = element_text(colour = "grey35", size = 9))
save2(p1, OUT1, "volcano_IFNB_vs_Mock_ISG_IDG", SINGLE_W, SINGLE_H)
message("Figure 1 saved in ", OUT1)

# ======================= 2. four-panel volcano =====================================================================================
bf <- file.path(OUT2, "transcriptome_blocks.csv")
if (!file.exists(bf)) {
  message("transcriptome_blocks.csv not found - run 51_mock_ifnb_ebov_analysis.R to make the four-panel figure. Done with figure 1 only.")
} else {
  B <- read_csv(bf, show_col_types = FALSE) |> mutate(class = case_when(resp == "Up" ~ "ISG", resp == "Down" ~ "IDG", TRUE ~ "other"))
  nI <- sum(B$class == "ISG"); nD <- sum(B$class == "IDG")
  PANELS <- c("IFNB vs Mock", "EBOV vs Mock", "IFNB+EBOV vs Mock", "IFNB+EBOV vs IFNB  (EBOV on top of IFNB)")
  long <- bind_rows(
    B |> transmute(gene_id, gene_name, class, lfc = lfc_B,    padj = padj_B,    panel = PANELS[1], lfc_ref = lfc_B),
    B |> transmute(gene_id, gene_name, class, lfc = lfc_E,    padj = padj_E,    panel = PANELS[2], lfc_ref = lfc_B),
    B |> transmute(gene_id, gene_name, class, lfc = lfc_BE,   padj = padj_BE,   panel = PANELS[3], lfc_ref = lfc_B),
    B |> transmute(gene_id, gene_name, class, lfc = lfc_BE_B, padj = padj_BE_B, panel = PANELS[4], lfc_ref = lfc_B)) |>
    filter(!is.na(lfc), !is.na(padj)) |>
    mutate(panel = factor(panel, levels = PANELS), y = -log10(pmax(padj, 1e-300)))
  cnt <- long |> filter(class != "other") |> group_by(panel, class) |>
    summarise(n = if (dplyr::first(as.character(panel)) == PANELS[4]) sum(padj < PADJ_CUT & sign(lfc) == -sign(lfc_ref))              # pulled back toward Mock
                  else sum(padj < PADJ_CUT & abs(lfc) > LFC_CUT & sign(lfc) == sign(lfc_ref)),                                           # still changed in the IFNB direction
              .groups = "drop")
  get_n <- function(pn, cl) { v <- cnt$n[cnt$panel == pn & cnt$class == cl]; if (length(v)) v else 0L }
  ann <- tibble(panel = factor(PANELS, levels = PANELS)) |>
    mutate(txt_isg = vapply(as.character(panel), function(pn) if (pn == PANELS[4]) sprintf("ISGs pulled down: %d / %d", get_n(pn, "ISG"), nI) else sprintf("ISGs still changed: %d / %d", get_n(pn, "ISG"), nI), character(1)),
           txt_idg = vapply(as.character(panel), function(pn) if (pn == PANELS[4]) sprintf("IDGs pulled up: %d / %d", get_n(pn, "IDG"), nD) else sprintf("IDGs still changed: %d / %d", get_n(pn, "IDG"), nD), character(1)))
  write_csv(cnt |> mutate(of = ifelse(class == "ISG", nI, nD)), file.path(OUT2, "volcano_four_panel_counts.csv"))
  print(as.data.frame(cnt |> pivot_wider(names_from = class, values_from = n)), row.names = FALSE)

  xr4 <- max(abs(long$lfc[long$panel == PANELS[1]])) * 1.05
  ymax <- as.numeric(quantile(long$y, 0.999))
  p2 <- ggplot(long, aes(lfc, y)) +
    geom_point(data = long |> filter(class == "other"), colour = "black", size = 0.5, alpha = 0.35) +
    geom_point(data = long |> filter(class != "other"), aes(colour = class), size = 0.9, alpha = 0.8) +
    geom_hline(yintercept = -log10(PADJ_CUT), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    geom_text(data = ann, aes(x = xr4 * 0.97, y = Inf, label = txt_isg), inherit.aes = FALSE, colour = RED, hjust = 1, vjust = 1.5, size = 3.1, fontface = "bold") +
    geom_text(data = ann, aes(x = -xr4 * 0.97, y = Inf, label = txt_idg), inherit.aes = FALSE, colour = BLUE, hjust = 0, vjust = 1.5, size = 3.1, fontface = "bold") +
    facet_wrap(~panel, ncol = 2) +
    scale_colour_manual(values = c(ISG = RED, IDG = BLUE), guide = "none") +
    scale_x_continuous(limits = c(-xr4, xr4), oob = scales::squish) +
    scale_y_continuous(limits = c(0, ymax), oob = scales::squish) +
    labs(title = "How EBOV modulates the IFNB response: four-panel volcano",
         subtitle = "colour is fixed by the IFNB-vs-Mock result in every panel: red = IFNB ISGs, blue = IFNB IDGs, black = all other protein-coding genes; dashed = log2FC +/-1 and padj 0.05",
         x = "log2 fold-change (shrunk)", y = expression(-log[10](padj))) + base_theme(11)
  save2(p2, OUT2, "volcano_four_panel", FOUR_W, FOUR_H)
  message("Figure 2 saved in ", OUT2)
}
