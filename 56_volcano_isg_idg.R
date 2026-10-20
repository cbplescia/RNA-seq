#!/usr/bin/env Rscript
# Volcano plot of IFNB vs Mock (protein-coding genes): ISGs red, IDGs blue, everything else black.
#   red   = ISG:  padj < PADJ_CUT and log2FC >  LFC_CUT
#   blue  = IDG:  padj < PADJ_CUT and log2FC < -LFC_CUT
#   black = all other genes - including genes beyond +/-LFC_CUT whose padj is NOT significant (colour is only given when the padj supports it)
# Dashed lines: log2FC = +/-LFC_CUT and the padj threshold. Shrunk log2FC from DE_v2_IFNB_vs_Mock_protein_coding.csv.
# Writes to <project>/requant_v2_outputs/figures_protein_coding/volcano_IFNB_vs_Mock_ISG_IDG.png / .pdf
# Usage: Rscript 56_volcano_isg_idg.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter

LFC_CUT <- 1; PADJ_CUT <- 0.05
INK <- "#1F2933"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs", "figures_protein_coding"); dir.create(OUT, showWarnings = FALSE)

d <- read_csv(file.path(PROJ, "requant_v2_outputs", "DE_v2_IFNB_vs_Mock_protein_coding.csv"), show_col_types = FALSE) |>
  filter(!is.na(padj)) |>
  mutate(y = -log10(pmax(padj, 1e-300)),
         class = case_when(padj < PADJ_CUT & log2FoldChange >  LFC_CUT ~ "ISG",
                           padj < PADJ_CUT & log2FoldChange < -LFC_CUT ~ "IDG", TRUE ~ "other"))
n <- table(factor(d$class, levels = c("ISG", "IDG", "other")))
message(sprintf("Protein-coding genes: %d | ISGs (red): %d | IDGs (blue): %d | other (black): %d", nrow(d), n[["ISG"]], n[["IDG"]], n[["other"]]))
message(sprintf("Genes beyond +/-%g log2FC but not significant (kept black): %d", LFC_CUT, sum(abs(d$log2FoldChange) > LFC_CUT & d$padj >= PADJ_CUT)))
write_csv(d |> filter(class != "other") |> arrange(class, desc(abs(log2FoldChange))) |> select(gene_id, gene_name, class, log2FoldChange, padj),
          file.path(OUT, "volcano_IFNB_vs_Mock_ISGs_IDGs.csv"))

xr <- max(abs(d$log2FoldChange)) * 1.05
p <- ggplot(d, aes(log2FoldChange, y)) +
  geom_point(data = d |> filter(class == "other"), colour = "black", size = 0.7, alpha = 0.45) +
  geom_point(data = d |> filter(class != "other"), aes(colour = class), size = 1.1, alpha = 0.85) +
  geom_hline(yintercept = -log10(PADJ_CUT), linetype = "dashed", colour = "grey50", linewidth = 0.35) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", colour = "grey50", linewidth = 0.35) +
  scale_colour_manual(values = c(ISG = "#D7191C", IDG = "#2B6CB0"), guide = "none") +
  scale_x_continuous(limits = c(-xr, xr)) +
  annotate("text", x = xr * 0.97, y = Inf, label = sprintf("ISGs: %d", n[["ISG"]]), colour = "#D7191C", hjust = 1, vjust = 1.6, fontface = "bold", size = 4.2) +
  annotate("text", x = -xr * 0.97, y = Inf, label = sprintf("IDGs: %d", n[["IDG"]]), colour = "#2B6CB0", hjust = 0, vjust = 1.6, fontface = "bold", size = 4.2) +
  labs(title = "IFNB vs Mock: protein-coding transcripts",
       subtitle = sprintf("red = ISGs (log2FC > %g, padj < %g); blue = IDGs (log2FC < -%g, padj < %g); black = all other genes; %d genes", LFC_CUT, PADJ_CUT, LFC_CUT, PADJ_CUT, nrow(d)),
       x = "log2 fold-change (IFNB vs Mock)", y = expression(-log[10](padj))) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey94", linewidth = 0.3),
        panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK), axis.title = element_text(colour = INK),
        plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 9))
ggsave(file.path(OUT, "volcano_IFNB_vs_Mock_ISG_IDG.png"), p, width = 8, height = 7, dpi = 300)
ggsave(file.path(OUT, "volcano_IFNB_vs_Mock_ISG_IDG.pdf"), p, width = 8, height = 7)
message("Saved to: ", OUT)
