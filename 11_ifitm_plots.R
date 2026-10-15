#!/usr/bin/env Rscript
# IFITM genes, Mock vs IFNB:
#   (1) boxplot of VST expression with every sample as a point, all IFITMs in
#       one panel (Mock and IFNB side by side for each gene)
#   (2) one bar graph of the IFNB vs Mock log2 fold-change for each IFITM
#       (shrunk values from DE_v2_IFNB_vs_Mock; bars show +/- lfcSE; stars show padj)
# IFITMs are found by symbol (Ifitm1, Ifitm2, Ifitm3, ...). Any that are not in
# the data (e.g. dropped by the low-expression filter in 04) are reported.
# Needs 07_add_gene_symbols.R to have been run.

suppressPackageStartupMessages({
  library(DESeq2); library(readr); library(dplyr); library(ggplot2)
})

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs")
FIG <- file.path(OUT, "figures_ifitm")
dir.create(FIG, showWarnings = FALSE)

INK <- "#1F2933"
COL <- c(Mock = "#8A8A8A", IFNB = "#2B6CB0")

de <- read_csv(file.path(OUT, "DE_v2_IFNB_vs_Mock_annotated.csv"), show_col_types = FALSE)
ifitm <- de |> filter(grepl("^Ifitm[0-9]+$", gene_name))
ifitm <- ifitm |> arrange(as.numeric(sub("Ifitm", "", gene_name)))
if (!nrow(ifitm)) stop("No genes named Ifitm* found in the DE table.")
ord <- ifitm$gene_name
message("IFITM genes in the DE table: ", paste(ord, collapse = ", "))
message("(Ifitm genes absent here were removed by the expression filter in 04 and cannot be plotted.)")

theme_fig <- function() {
  theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
          panel.border = element_rect(colour = "grey70"),
          axis.text = element_text(colour = INK), axis.title = element_text(colour = INK),
          plot.title = element_text(colour = INK, face = "bold"),
          plot.subtitle = element_text(colour = "grey35"),
          axis.text.x = element_text(face = "italic"),
          legend.position = "bottom")
}
save_gg <- function(p, base, w, h) {
  ggsave(paste0(base, ".png"), p, width = w, height = h, dpi = 300)
  ggsave(paste0(base, ".pdf"), p, width = w, height = h)
}

# ---- (1) VST boxplot -----------------------------------------------------------
vsd  <- readRDS(file.path(OUT, "vsd_v2.rds"))
vmat <- assay(vsd)
samp <- tibble(sample_id = colnames(vmat), group = as.character(vsd$group)) |>
  filter(group %in% c("Mock_Mock", "IFNB_Mock")) |>
  mutate(treatment = factor(ifelse(group == "Mock_Mock", "Mock", "IFNB"), levels = c("Mock", "IFNB")))

ids <- ifitm$gene_id[ifitm$gene_id %in% rownames(vmat)]
long <- as.data.frame(t(vmat[ids, samp$sample_id, drop = FALSE])) |>
  tibble::rownames_to_column("sample_id") |>
  tidyr::pivot_longer(-sample_id, names_to = "gene_id", values_to = "vst") |>
  left_join(samp, by = "sample_id") |>
  left_join(ifitm |> select(gene_id, gene_name), by = "gene_id") |>
  mutate(gene_name = factor(gene_name, levels = ord[ord %in% gene_name]))
write_csv(long, file.path(FIG, "ifitm_vst_values.csv"))

p1 <- ggplot(long, aes(gene_name, vst, fill = treatment, colour = treatment)) +
  geom_boxplot(alpha = 0.25, outlier.shape = NA, width = 0.65, linewidth = 0.5,
               position = position_dodge(width = 0.75)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75, seed = 1),
             size = 1.8, alpha = 0.9) +
  scale_fill_manual(values = COL, name = NULL) +
  scale_colour_manual(values = COL, name = NULL) +
  labs(title = "IFITM expression: Mock vs IFNB",
       subtitle = "VST-normalised expression; each point is one sample (n = 3 per group)",
       x = NULL, y = "VST expression (log2 scale)") +
  theme_fig()
save_gg(p1, file.path(FIG, "ifitm_vst_boxplot_Mock_vs_IFNB"), 7, 5.5)

# ---- (2) log2FC bar graph ---------------------------------------------------------
bar <- ifitm |>
  mutate(gene_name = factor(gene_name, levels = ord),
         se = if ("lfcSE" %in% names(ifitm)) lfcSE else 0,
         stars = case_when(is.na(padj) ~ "n/a", padj < 0.001 ~ "***", padj < 0.01 ~ "**",
                           padj < 0.05 ~ "*", TRUE ~ "ns"),
         lab_y = ifelse(log2FoldChange >= 0, log2FoldChange + se, log2FoldChange - se))
write_csv(bar |> select(gene_id, gene_name, log2FoldChange, any_of("lfcSE"), padj, stars),
          file.path(FIG, "ifitm_log2FC_IFNB_vs_Mock.csv"))

p2 <- ggplot(bar, aes(gene_name, log2FoldChange)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_col(fill = COL[["IFNB"]], width = 0.6) +
  { if (any(bar$se > 0)) geom_errorbar(aes(ymin = log2FoldChange - se, ymax = log2FoldChange + se),
                                       width = 0.2, colour = INK, linewidth = 0.4) } +
  geom_text(aes(y = lab_y, label = stars, vjust = ifelse(log2FoldChange >= 0, -0.4, 1.3)),
            colour = INK, size = 3.6) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
  labs(title = "IFITM induction by IFNB",
       subtitle = "log2 fold-change, IFNB vs Mock (shrunk); error bars = lfcSE; * padj<0.05, ** <0.01, *** <0.001",
       x = NULL, y = "log2 fold-change (IFNB vs Mock)") +
  theme_fig()
save_gg(p2, file.path(FIG, "ifitm_log2FC_bar_IFNB_vs_Mock"), 6.5, 5.5)

message("\nPlotted: ", paste(levels(long$gene_name), collapse = ", "))
message("Figures and the plotted values (csv) are in: ", FIG)
