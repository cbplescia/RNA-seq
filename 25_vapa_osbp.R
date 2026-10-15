#!/usr/bin/env Rscript
# Vapa and Osbp, Mock vs IFNB:
#   (1) VST expression boxplots, every sample as a point
#   (2) one bar graph of shrunk log2FC, IFNB vs Mock
# Needs 07_add_gene_symbols.R to have been run. A gene missing from the data
# (removed by the low-expression filter in 04) is reported and left off.

suppressPackageStartupMessages({
  library(DESeq2); library(readr); library(dplyr); library(ggplot2)
})

GENES <- c("Vapa", "Osbp")
LFC_CUT <- 1; PADJ_CUT <- 0.05
INK <- "#1F2933"
COL <- c(Mock = "#8A8A8A", IFNB = "#2B6CB0")
DIR_COL <- c(Up = "#C53030", Down = "#2B6CB0", `Not significant` = "#B8B8B8")

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs")
FIG <- file.path(OUT, "figures_genes")
dir.create(FIG, showWarnings = FALSE)

de <- read_csv(file.path(OUT, "DE_v2_IFNB_vs_Mock_annotated.csv"), show_col_types = FALSE)
if (!"lfcSE" %in% names(de)) de$lfcSE <- NA_real_
hit <- tibble(symbol = GENES) |> mutate(key = toupper(symbol)) |>
  left_join(de |> mutate(key = toupper(gene_name)) |> distinct(key, .keep_all = TRUE) |>
              select(key, gene_id, log2FoldChange, lfcSE, padj), by = "key")
if (any(is.na(hit$gene_id))) message("Not in the data (left off): ", paste(hit$symbol[is.na(hit$gene_id)], collapse = ", "))
hit <- hit |> filter(!is.na(gene_id))
if (!nrow(hit)) stop("Neither gene was found in the data.")

fmt_p <- function(p) ifelse(is.na(p), "padj n/a", ifelse(p < 0.001, "padj < 0.001", sprintf("padj = %.3g", p)))
hit <- hit |> mutate(panel = sprintf("%s\nlog2FC %+.2f, %s", symbol, log2FoldChange, fmt_p(padj)))
hit$panel <- factor(hit$panel, levels = hit$panel)

# ---- (1) VST boxplots -----------------------------------------------------------------
vsd  <- readRDS(file.path(OUT, "vsd_v2.rds"))
vmat <- assay(vsd)
samp <- tibble(sample_id = colnames(vmat), group = as.character(vsd$group)) |>
  filter(group %in% c("Mock_Mock", "IFNB_Mock")) |>
  mutate(treatment = factor(ifelse(group == "Mock_Mock", "Mock", "IFNB"), levels = c("Mock", "IFNB")))
long <- as.data.frame(t(vmat[hit$gene_id, samp$sample_id, drop = FALSE])) |>
  tibble::rownames_to_column("sample_id") |>
  tidyr::pivot_longer(-sample_id, names_to = "gene_id", values_to = "vst") |>
  left_join(samp, by = "sample_id") |>
  left_join(hit |> select(gene_id, symbol, panel), by = "gene_id")
write_csv(long |> select(-panel), file.path(FIG, "vapa_osbp_vst_values.csv"))

p1 <- ggplot(long, aes(treatment, vst, fill = treatment, colour = treatment)) +
  geom_boxplot(alpha = 0.25, outlier.shape = NA, width = 0.55, linewidth = 0.5) +
  geom_jitter(width = 0.1, height = 0, size = 2, alpha = 0.9) +
  facet_wrap(~panel, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = COL, guide = "none") + scale_colour_manual(values = COL, guide = "none") +
  labs(title = "Vapa and Osbp: Mock vs IFNB",
       subtitle = "VST-normalised expression; each point is one sample (n = 3 per group)",
       x = NULL, y = "VST expression (log2 scale)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        panel.border = element_rect(colour = "grey70"),
        strip.background = element_rect(fill = "grey95", colour = "grey70"),
        strip.text = element_text(colour = INK, face = "bold"),
        axis.text = element_text(colour = INK), axis.title = element_text(colour = INK),
        plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35"))
w <- 2.8 * nrow(hit) + 1
ggsave(file.path(FIG, "vapa_osbp_vst_boxplot.png"), p1, width = w, height = 5, dpi = 300)
ggsave(file.path(FIG, "vapa_osbp_vst_boxplot.pdf"), p1, width = w, height = 5)

# ---- (2) log2FC bar graph ---------------------------------------------------------------
bar <- hit |> mutate(
  symbol = factor(symbol, levels = GENES), se = coalesce(lfcSE, 0),
  direction = factor(case_when(!is.na(padj) & padj < PADJ_CUT & log2FoldChange >  LFC_CUT ~ "Up",
                               !is.na(padj) & padj < PADJ_CUT & log2FoldChange < -LFC_CUT ~ "Down",
                               TRUE ~ "Not significant"), levels = names(DIR_COL)),
  stars = case_when(is.na(padj) ~ "", padj < 0.001 ~ "***", padj < 0.01 ~ "**", padj < 0.05 ~ "*", TRUE ~ "ns"),
  top = ifelse(log2FoldChange >= 0, log2FoldChange + se, log2FoldChange - se))
write_csv(bar |> select(symbol, gene_id, log2FoldChange, lfcSE, padj, stars, direction),
          file.path(FIG, "vapa_osbp_log2FC_IFNB_vs_Mock.csv"))
p2 <- ggplot(bar, aes(symbol, log2FoldChange, fill = direction)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = log2FoldChange - se, ymax = log2FoldChange + se), width = 0.2, colour = INK, linewidth = 0.4) +
  geom_text(aes(y = top, label = stars, vjust = ifelse(log2FoldChange >= 0, -0.4, 1.3)), colour = INK, size = 3.6) +
  scale_fill_manual(values = DIR_COL, drop = FALSE, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.15))) +
  labs(title = "Vapa and Osbp: IFNB vs Mock",
       subtitle = sprintf("shrunk log2 fold-change ± lfcSE; coloured when padj<%g & |log2FC|>%g; * padj<0.05, ** <0.01, *** <0.001", PADJ_CUT, LFC_CUT),
       x = NULL, y = "log2 fold-change (IFNB vs Mock)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        panel.border = element_rect(colour = "grey70"),
        axis.text = element_text(colour = INK), axis.text.x = element_text(face = "italic"),
        axis.title = element_text(colour = INK),
        plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 8.5),
        legend.position = "bottom")
ggsave(file.path(FIG, "vapa_osbp_log2FC_bar_IFNB_vs_Mock.png"), p2, width = 4.5, height = 5.5, dpi = 300)
ggsave(file.path(FIG, "vapa_osbp_log2FC_bar_IFNB_vs_Mock.pdf"), p2, width = 4.5, height = 5.5)
print(as.data.frame(bar |> select(symbol, log2FoldChange, padj, stars)), row.names = FALSE)
message("\nSaved to: ", FIG)
