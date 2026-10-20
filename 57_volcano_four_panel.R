#!/usr/bin/env Rscript
# Four-panel volcano: how the IFNB response is (or is not) modulated by EBOV alone and by IFNB+EBOV.
# Every gene keeps its colour from the IFNB analysis in ALL panels, so you can follow the IFNB-defined genes:
#   red  = IFNB ISGs (IFNB vs Mock: padj < PADJ_CUT, log2FC >  LFC_CUT)      blue = IFNB IDGs (log2FC < -LFC_CUT)      black = all other genes
# Panels (shrunk log2FC from DESeq2; protein-coding genes):
#   top left      IFNB vs Mock                 the response itself (this is where red / blue are defined)
#   top right     EBOV vs Mock                 does EBOV alone move the IFNB-defined genes?
#   bottom left   IFNB+EBOV vs Mock            how much of the IFNB response is left with EBOV present
#   bottom right  IFNB+EBOV vs IFNB            what EBOV does on top of IFNB: red genes shifting left / blue genes shifting right = suppression
# Each panel reports how many of the ISGs / IDGs are still significant (padj < PADJ_CUT, same direction, |log2FC| > LFC_CUT); for the
# bottom-right panel it reports how many are significantly pulled back toward Mock. All panels share the same x range.
# Uses figures_mock_ifnb_ebov/transcriptome_blocks.csv (written by 51_mock_ifnb_ebov_analysis.R; run 51 first).
# Writes to figures_mock_ifnb_ebov/volcano_four_panel.png / .pdf
# Usage: Rscript 57_volcano_four_panel.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter

LFC_CUT <- 1; PADJ_CUT <- 0.05
INK <- "#1F2933"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs", "figures_mock_ifnb_ebov")
bf <- file.path(OUT, "transcriptome_blocks.csv")
if (!file.exists(bf)) stop("Run 51_mock_ifnb_ebov_analysis.R first (missing ", bf, ")")
B <- read_csv(bf, show_col_types = FALSE) |>
  mutate(class = case_when(resp == "Up" ~ "ISG", resp == "Down" ~ "IDG", TRUE ~ "other"))
nI <- sum(B$class == "ISG"); nD <- sum(B$class == "IDG")
message(sprintf("Protein-coding genes: %d | IFNB ISGs %d | IFNB IDGs %d", nrow(B), nI, nD))

PANELS <- c("IFNB vs Mock", "EBOV vs Mock", "IFNB+EBOV vs Mock", "IFNB+EBOV vs IFNB  (EBOV on top of IFNB)")
long <- bind_rows(
  B |> transmute(gene_id, gene_name, class, lfc = lfc_B,    padj = padj_B,    panel = PANELS[1], lfc_ref = lfc_B),
  B |> transmute(gene_id, gene_name, class, lfc = lfc_E,    padj = padj_E,    panel = PANELS[2], lfc_ref = lfc_B),
  B |> transmute(gene_id, gene_name, class, lfc = lfc_BE,   padj = padj_BE,   panel = PANELS[3], lfc_ref = lfc_B),
  B |> transmute(gene_id, gene_name, class, lfc = lfc_BE_B, padj = padj_BE_B, panel = PANELS[4], lfc_ref = lfc_B)) |>
  filter(!is.na(lfc), !is.na(padj)) |>
  mutate(panel = factor(panel, levels = PANELS), y = -log10(pmax(padj, 1e-300)))

# per-panel counts for the annotations
cnt <- long |> filter(class != "other") |> group_by(panel, class) |>
  summarise(n = if (first(as.character(panel)) == PANELS[4])
                  sum(padj < PADJ_CUT & sign(lfc) == -sign(lfc_ref))                       # significantly pulled back toward Mock
                else sum(padj < PADJ_CUT & abs(lfc) > LFC_CUT & sign(lfc) == sign(lfc_ref)),  # still significantly changed in the IFNB direction
            .groups = "drop")
get_n <- function(pn, cl) { v <- cnt$n[cnt$panel == pn & cnt$class == cl]; if (length(v)) v else 0L }
ann <- tibble(panel = factor(PANELS, levels = PANELS)) |>
  mutate(txt_isg = vapply(as.character(panel), function(pn) if (pn == PANELS[4]) sprintf("ISGs pulled down: %d / %d", get_n(pn, "ISG"), nI) else sprintf("ISGs still changed: %d / %d", get_n(pn, "ISG"), nI), character(1)),
         txt_idg = vapply(as.character(panel), function(pn) if (pn == PANELS[4]) sprintf("IDGs pulled up: %d / %d", get_n(pn, "IDG"), nD) else sprintf("IDGs still changed: %d / %d", get_n(pn, "IDG"), nD), character(1)))
write_csv(cnt |> mutate(of = ifelse(class == "ISG", nI, nD)), file.path(OUT, "volcano_four_panel_counts.csv"))
print(as.data.frame(cnt |> pivot_wider(names_from = class, values_from = n)), row.names = FALSE)

xr <- max(abs(long$lfc[long$panel == PANELS[1]])) * 1.05
ymax <- as.numeric(quantile(long$y, 0.999))
p <- ggplot(long, aes(lfc, y)) +
  geom_point(data = long |> filter(class == "other"), colour = "black", size = 0.5, alpha = 0.35) +
  geom_point(data = long |> filter(class != "other"), aes(colour = class), size = 0.9, alpha = 0.8) +
  geom_hline(yintercept = -log10(PADJ_CUT), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_text(data = ann, aes(x = xr * 0.97, y = Inf, label = txt_isg), inherit.aes = FALSE, colour = "#D7191C", hjust = 1, vjust = 1.5, size = 3.1, fontface = "bold") +
  geom_text(data = ann, aes(x = -xr * 0.97, y = Inf, label = txt_idg), inherit.aes = FALSE, colour = "#2B6CB0", hjust = 0, vjust = 1.5, size = 3.1, fontface = "bold") +
  facet_wrap(~panel, ncol = 2) +
  scale_colour_manual(values = c(ISG = "#D7191C", IDG = "#2B6CB0"), guide = "none") +
  scale_x_continuous(limits = c(-xr, xr), oob = scales::squish) +
  scale_y_continuous(limits = c(0, ymax), oob = scales::squish) +
  labs(title = "How EBOV modulates the IFNB response: four-panel volcano",
       subtitle = "colour is fixed by the IFNB-vs-Mock result in every panel: red = IFNB ISGs, blue = IFNB IDGs, black = all other protein-coding genes; dashed = log2FC +/-1 and padj 0.05",
       x = "log2 fold-change (shrunk)", y = expression(-log[10](padj))) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey94", linewidth = 0.3), panel.border = element_rect(colour = "grey70"),
        axis.text = element_text(colour = INK), axis.title = element_text(colour = INK), plot.title = element_text(colour = INK, face = "bold"),
        plot.subtitle = element_text(colour = "grey35", size = 8.5), strip.background = element_rect(fill = "grey95", colour = "grey70"),
        strip.text = element_text(colour = INK, face = "bold"))
ggsave(file.path(OUT, "volcano_four_panel.png"), p, width = 12, height = 10, dpi = 300)
ggsave(file.path(OUT, "volcano_four_panel.pdf"), p, width = 12, height = 10)
message("\nSaved volcano_four_panel.png / .pdf in ", OUT)
