#!/usr/bin/env Rscript
# Volcano plots of all protein-coding genes for IFNB vs Mock and IFNY vs Mock, with
# the top 10 most upregulated genes of EACH contrast labelled (highest log2FC among
# genes with padj < 0.05 and log2FC > 1). Two versions of the same data:
#   *_panels   side by side, shared axes; red = significantly up, blue = down
#   *_overlay  both interferons in one plot; blue = IFNB, orange = IFNY
# Needs 07 and 08 to have been run, and the ggrepel package.

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })
if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Install ggrepel (conda install -c conda-forge r-ggrepel)")

LFC_CUT <- 1; PADJ_CUT <- 0.05; N_LABEL <- 10
INK <- "#1F2933"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs")
FIG <- file.path(OUT, "figures_protein_coding")
dir.create(FIG, showWarnings = FALSE)

rd <- function(nm, lab) read_csv(file.path(OUT, paste0("DE_v2_", nm, "_protein_coding.csv")),
                                 show_col_types = FALSE) |>
  transmute(gene_name, lfc = log2FoldChange, padj, contrast = lab) |> filter(!is.na(padj))
d <- bind_rows(rd("IFNB_vs_Mock", "IFNB vs Mock"), rd("IFNY_vs_Mock", "IFNγ vs Mock")) |>
  mutate(y = -log10(pmax(padj, 1e-300)),
         status = case_when(padj < PADJ_CUT & lfc >  LFC_CUT ~ "Up",
                            padj < PADJ_CUT & lfc < -LFC_CUT ~ "Down", TRUE ~ "Not significant"),
         contrast = factor(contrast, levels = c("IFNB vs Mock", "IFNγ vs Mock")))
lab <- d |> filter(status == "Up") |> group_by(contrast) |> slice_max(lfc, n = N_LABEL) |> ungroup()
message("Top ", N_LABEL, " upregulated per contrast:")
for (cn in levels(d$contrast)) message("  ", cn, ": ", paste(lab$gene_name[lab$contrast == cn], collapse = ", "))
message("Significant up/down: ", paste(capture.output(print(table(d$contrast, d$status))), collapse = "\n"))

xr <- max(abs(d$lfc)) * 1.05
base <- list(
  geom_hline(yintercept = -log10(PADJ_CUT), linetype = "dashed", colour = "grey60", linewidth = 0.3),
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", colour = "grey60", linewidth = 0.3),
  scale_x_continuous(limits = c(-xr, xr)),
  labs(x = "log2 fold-change vs Mock (shrunk)", y = expression(-log[10](padj))),
  theme_bw(base_size = 12),
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.border = element_rect(colour = "grey70"),
        axis.text = element_text(colour = INK), axis.title = element_text(colour = INK),
        plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35"),
        strip.background = element_rect(fill = "grey95", colour = "grey70"),
        strip.text = element_text(colour = INK, face = "bold"), legend.position = "bottom"))
repel <- function(data) ggrepel::geom_text_repel(data = data, aes(label = gene_name), colour = INK, size = 3,
  max.overlaps = Inf, min.segment.length = 0, segment.colour = "grey60", segment.size = 0.25,
  box.padding = 0.4, seed = 1)
save2 <- function(p, base_name, w, h) {
  ggsave(file.path(FIG, paste0(base_name, ".png")), p, width = w, height = h, dpi = 300)
  ggsave(file.path(FIG, paste0(base_name, ".pdf")), p, width = w, height = h)
}
sub <- sprintf("protein-coding genes; coloured when padj<%g & |log2FC|>%g; top %d upregulated genes per interferon labelled",
               PADJ_CUT, LFC_CUT, N_LABEL)

# ---- side-by-side panels ---------------------------------------------------------------
p1 <- ggplot(d, aes(lfc, y)) + base +
  geom_point(aes(colour = status), size = 0.8, alpha = 0.55) +
  geom_point(data = lab, aes(colour = status), size = 1.4) + repel(lab) +
  facet_wrap(~contrast, nrow = 1) +
  scale_colour_manual(values = c(Up = "#C53030", Down = "#2B6CB0", `Not significant` = "#B8B8B8"), name = NULL) +
  labs(title = "Interferon-induced changes in protein-coding genes", subtitle = sub)
save2(p1, "volcano_IFNB_IFNY_panels", 12, 6.5)

# ---- single overlaid plot ----------------------------------------------------------------
ov_col <- c(`IFNB vs Mock` = "#2B6CB0", `IFNγ vs Mock` = "#D69E2E")
sig <- d |> filter(status != "Not significant")
p2 <- ggplot(d, aes(lfc, y)) + base +
  geom_point(data = d |> filter(status == "Not significant"), colour = "#D5D5D5", size = 0.7, alpha = 0.5) +
  geom_point(data = sig, aes(colour = contrast), size = 0.9, alpha = 0.6) +
  geom_point(data = lab, aes(colour = contrast), size = 1.6) + repel(lab) +
  scale_colour_manual(values = ov_col, name = NULL) +
  labs(title = "IFNB and IFNγ vs Mock: protein-coding genes", subtitle = sub)
save2(p2, "volcano_IFNB_IFNY_overlay", 8.5, 7)
message("\nSaved to: ", FIG)
