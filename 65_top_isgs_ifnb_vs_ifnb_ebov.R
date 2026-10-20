#!/usr/bin/env Rscript
# The top N_TOP IFNB ISGs (IFNB vs Mock: padj < PADJ_CUT, shrunk log2FC > LFC_CUT, ranked by log2FC; protein-coding genes, the same ranking as the volcano labels).
# For each gene two bars: log2FC of IFNB vs Mock (blue) and IFNB+EBOV vs Mock (light green), error bars +/- 1 lfcSE. A bracket joins the two bars and gives the
# stat of IFNB+EBOV vs IFNB: DESeq2 Wald test on the difference between the two log2FCs (contrast IFNB_EBOV vs IFNB_Mock, BH padj across all genes),
# i.e. whether EBOV significantly changes the IFNB induction of that gene (* < 0.05, ** < 0.01, *** < 0.001, ns). Stars above each bar (small) = padj vs Mock.
# Writes top_ISGs_IFNB_vs_IFNB_EBOV.png / .pdf and top_ISGs_IFNB_vs_IFNB_EBOV.csv to requant_v2_outputs/figures_gene_bars/
# Usage: Rscript 65_top_isgs_ifnb_vs_ifnb_ebov.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2); library(DESeq2) })
select <- dplyr::select; filter <- dplyr::filter

N_TOP <- 10; PADJ_CUT <- 0.05; LFC_CUT <- 1
INK <- "#1F2933"; C_B <- "#2B6CB0"; C_BE <- "#9BD08E"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_gene_bars"); dir.create(OUT, showWarnings = FALSE)
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

pcg <- read_csv(file.path(RNA, "DE_v2_IFNB_vs_Mock_protein_coding.csv"), show_col_types = FALSE) |> filter(!is.na(padj))
top <- pcg |> filter(padj < PADJ_CUT, log2FoldChange > LFC_CUT) |> arrange(desc(log2FoldChange)) |> slice_head(n = N_TOP) |>
  transmute(gene_id, gene_name = ifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name), lfc_B = log2FoldChange, se_B = lfcSE, padj_B = padj)
be <- read_csv(file.path(RNA, "DE_v2_IFNB_EBOV_vs_Mock_annotated.csv"), show_col_types = FALSE) |> transmute(gene_id, lfc_BE = log2FoldChange, se_BE = lfcSE, padj_BE = padj)
dds <- readRDS(file.path(RNA, "dds_v2_decoyAware_hostNormalised.rds"))
df <- as.data.frame(results(dds, contrast = c("group", "IFNB_EBOV", "IFNB_Mock"))) |> tibble::rownames_to_column("gene_id") |> transmute(gene_id, diff_lfc = log2FoldChange, diff_padj = padj)
T <- top |> left_join(be, by = "gene_id") |> left_join(df, by = "gene_id") |> mutate(gene_name = factor(gene_name, levels = gene_name))
write_csv(T |> transmute(gene_name, gene_id, log2FC_IFNB_vs_Mock = lfc_B, padj_IFNB_vs_Mock = padj_B, log2FC_IFNB_EBOV_vs_Mock = lfc_BE, padj_IFNB_EBOV_vs_Mock = padj_BE,
                         log2FC_diff_IFNB_EBOV_vs_IFNB_unshrunk = diff_lfc, padj_IFNB_EBOV_vs_IFNB = diff_padj), file.path(OUT, "top_ISGs_IFNB_vs_IFNB_EBOV.csv"))
print(as.data.frame(T |> transmute(gene = gene_name, IFNB = round(lfc_B, 2), IFNB_EBOV = round(lfc_BE, 2), padj_diff = signif(diff_padj, 3), stars = stars(diff_padj))), row.names = FALSE)

L <- bind_rows(T |> transmute(gene_name, cond = "IFNB vs Mock", lfc = lfc_B, se = se_B, padj = padj_B),
               T |> transmute(gene_name, cond = "IFNB+EBOV vs Mock", lfc = lfc_BE, se = se_BE, padj = padj_BE)) |>
  mutate(cond = factor(cond, levels = c("IFNB vs Mock", "IFNB+EBOV vs Mock")), top = lfc + se)
pad <- 0.03 * max(L$top, na.rm = TRUE)
bk <- L |> group_by(gene_name) |> summarise(y = max(top, na.rm = TRUE) + 2.2 * pad, .groups = "drop") |> left_join(T |> select(gene_name, diff_padj), by = "gene_name") |>
  mutate(x = as.numeric(gene_name), lab = stars(diff_padj))
dg <- position_dodge(width = 0.8)
p <- ggplot(L, aes(gene_name, lfc, fill = cond)) + geom_col(position = dg, width = 0.8) +
  geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), position = dg, width = 0.22, colour = INK, linewidth = 0.35) +
  geom_text(aes(y = top + 0.4 * pad, label = stars(padj)), position = dg, colour = "grey40", size = 2.4, vjust = 0) +
  geom_segment(data = bk, aes(x = x - 0.2, xend = x + 0.2, y = y, yend = y), inherit.aes = FALSE, colour = INK, linewidth = 0.4) +
  geom_segment(data = bk, aes(x = x - 0.2, xend = x - 0.2, y = y, yend = y - pad), inherit.aes = FALSE, colour = INK, linewidth = 0.4) +
  geom_segment(data = bk, aes(x = x + 0.2, xend = x + 0.2, y = y, yend = y - pad), inherit.aes = FALSE, colour = INK, linewidth = 0.4) +
  geom_text(data = bk, aes(x = x, y = y + 0.4 * pad, label = lab), inherit.aes = FALSE, colour = INK, size = 4, vjust = 0, fontface = "bold") +
  geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) + scale_fill_manual(values = c(`IFNB vs Mock` = C_B, `IFNB+EBOV vs Mock` = C_BE), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(title = sprintf("Top %d IFNB ISGs: IFNB vs IFNB+EBOV", N_TOP),
       subtitle = "bars = shrunk log2FC vs Mock +/- 1 SE (small grey stars = padj vs Mock); bracket = IFNB+EBOV vs IFNB (does EBOV change the induction?)", x = NULL, y = "log2 fold-change vs Mock") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK), axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
        axis.title = element_text(colour = INK), plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 8.5), legend.position = "top")
ggsave(file.path(OUT, "top_ISGs_IFNB_vs_IFNB_EBOV.png"), p, width = 9, height = 5.5, dpi = 300)
ggsave(file.path(OUT, "top_ISGs_IFNB_vs_IFNB_EBOV.pdf"), p, width = 9, height = 5.5)
message("Saved to: ", OUT)
