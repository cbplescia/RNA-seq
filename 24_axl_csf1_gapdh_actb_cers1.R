#!/usr/bin/env Rscript
# Axl, Csf1 (M-CSF), Gapdh, Actb (beta-actin) and Cers1 (ceramide synthase 1).
#   (1) VST expression boxplots, Mock vs IFNB vs IFNy, every sample as a point
#   (2) one grouped bar graph of shrunk log2FC: IFNB vs Mock and IFNy vs Mock
# Needs 07_add_gene_symbols.R to have been run. A gene missing from the data
# (removed by the low-expression filter in 04) is reported and left off.

suppressPackageStartupMessages({
  library(DESeq2); library(readr); library(dplyr); library(ggplot2)
})

GENES <- tibble::tribble(
  ~symbol, ~title,
  "Axl",   "Axl",
  "Csf1",  "Csf1\n(M-CSF)",
  "Gapdh", "Gapdh",
  "Actb",  "Actb\n(\u03b2-actin)",
  "Cers1", "Cers1\n(ceramide synthase 1)")
NCOL <- 3      # boxplot panels per row; genes are drawn in the order above

INK <- "#1F2933"
G_COL <- c(Mock = "#8A8A8A", IFNB = "#2B6CB0", `IFNγ` = "#D69E2E")
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs")
FIG <- file.path(OUT, "figures_genes")
dir.create(FIG, showWarnings = FALSE)

rd <- function(nm, lab) read_csv(file.path(OUT, paste0("DE_v2_", nm, "_annotated.csv")), show_col_types = FALSE) |>
  mutate(key = toupper(gene_name), contrast = lab) |>
  select(key, gene_id, log2FoldChange, any_of("lfcSE"), padj, contrast)
de <- bind_rows(rd("IFNB_vs_Mock", "IFNB vs Mock"), rd("IFNY_vs_Mock", "IFNγ vs Mock"))
if (!"lfcSE" %in% names(de)) de$lfcSE <- NA_real_

genes <- GENES |> mutate(key = toupper(symbol)) |>
  left_join(de |> distinct(key, gene_id), by = "key")
miss <- genes$symbol[is.na(genes$gene_id)]
if (length(miss)) message("Not in the data (left off): ", paste(miss, collapse = ", "))
genes <- genes |> filter(!is.na(gene_id))
if (!nrow(genes)) stop("Neither gene was found in the data.")

# ---- (1) VST boxplots --------------------------------------------------------------------
vsd  <- readRDS(file.path(OUT, "vsd_v2.rds"))
vmat <- assay(vsd)
samp <- tibble(sample_id = colnames(vmat), group = as.character(vsd$group)) |>
  filter(group %in% c("Mock_Mock", "IFNB_Mock", "IFNY_Mock")) |>
  mutate(treatment = factor(recode(group, Mock_Mock = "Mock", IFNB_Mock = "IFNB", IFNY_Mock = "IFNγ"),
                            levels = names(G_COL)))
long <- as.data.frame(t(vmat[genes$gene_id, samp$sample_id, drop = FALSE])) |>
  tibble::rownames_to_column("sample_id") |>
  tidyr::pivot_longer(-sample_id, names_to = "gene_id", values_to = "vst") |>
  left_join(samp, by = "sample_id") |>
  left_join(genes |> select(gene_id, symbol, title), by = "gene_id")
long$title <- factor(long$title, levels = genes$title)
write_csv(long |> select(-title), file.path(FIG, "axl_csf1_gapdh_actb_cers1_vst_values.csv"))

p1 <- ggplot(long, aes(treatment, vst, fill = treatment, colour = treatment)) +
  geom_boxplot(alpha = 0.25, outlier.shape = NA, width = 0.55, linewidth = 0.5) +
  geom_jitter(width = 0.1, height = 0, size = 2, alpha = 0.9) +
  facet_wrap(~title, ncol = NCOL, scales = "free_y") +
  scale_fill_manual(values = G_COL, guide = "none") + scale_colour_manual(values = G_COL, guide = "none") +
  labs(title = "Selected genes: Mock vs IFNB vs IFNγ",
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
n <- nrow(genes)
w1 <- 3.2 * min(n, NCOL) + 1
h1 <- 4.2 * ceiling(n / NCOL) + 1.2
ggsave(file.path(FIG, "axl_csf1_gapdh_actb_cers1_vst_boxplot.png"), p1, width = w1, height = h1, dpi = 300)
ggsave(file.path(FIG, "axl_csf1_gapdh_actb_cers1_vst_boxplot.pdf"), p1, width = w1, height = h1)

# ---- (2) grouped log2FC bar graph ------------------------------------------------------
bar <- de |> filter(gene_id %in% genes$gene_id) |>
  left_join(genes |> select(gene_id, symbol), by = "gene_id") |>
  mutate(symbol = factor(symbol, levels = genes$symbol),
         contrast = factor(contrast, levels = c("IFNB vs Mock", "IFNγ vs Mock")),
         se = coalesce(lfcSE, 0),
         stars = case_when(is.na(padj) ~ "", padj < 0.001 ~ "***", padj < 0.01 ~ "**", padj < 0.05 ~ "*", TRUE ~ "ns"),
         top = ifelse(log2FoldChange >= 0, log2FoldChange + se, log2FoldChange - se))
write_csv(bar |> select(symbol, contrast, gene_id, log2FoldChange, any_of("lfcSE"), padj, stars),
          file.path(FIG, "axl_csf1_gapdh_actb_cers1_log2FC.csv"))
dodge <- position_dodge(width = 0.75)
p2 <- ggplot(bar, aes(symbol, log2FoldChange, fill = contrast)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_col(width = 0.7, position = dodge) +
  geom_errorbar(aes(ymin = log2FoldChange - se, ymax = log2FoldChange + se), width = 0.2,
                position = dodge, colour = INK, linewidth = 0.4) +
  geom_text(aes(y = top, label = stars, vjust = ifelse(log2FoldChange >= 0, -0.4, 1.3)),
            position = dodge, colour = INK, size = 3.4) +
  scale_fill_manual(values = c(`IFNB vs Mock` = G_COL[["IFNB"]], `IFNγ vs Mock` = G_COL[["IFNγ"]]), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.15))) +
  labs(title = "Selected genes: response to IFNB and IFNγ",
       subtitle = "shrunk log2 fold-change ± lfcSE vs Mock; * padj<0.05, ** <0.01, *** <0.001, ns = not significant",
       x = NULL, y = "log2 fold-change vs Mock") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        panel.border = element_rect(colour = "grey70"),
        axis.text = element_text(colour = INK), axis.text.x = element_text(face = "italic"),
        axis.title = element_text(colour = INK),
        plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 9),
        legend.position = "bottom")
ggsave(file.path(FIG, "axl_csf1_gapdh_actb_cers1_log2FC_bar.png"), p2, width = 1.6 * n + 2.5, height = 5.5, dpi = 300)
ggsave(file.path(FIG, "axl_csf1_gapdh_actb_cers1_log2FC_bar.pdf"), p2, width = 1.6 * n + 2.5, height = 5.5)
print(as.data.frame(bar |> select(symbol, contrast, log2FoldChange, padj, stars)), row.names = FALSE)
message("\nSaved to: ", FIG)
