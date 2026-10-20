#!/usr/bin/env Rscript
# Bar graphs of log2FC for selected genes (edit GENES). Shrunk log2FC from the DESeq2 results (DE_v2_*_annotated.csv), error bars = +/- 1 lfcSE,
# stars = padj of that comparison (* < 0.05, ** < 0.01, *** < 0.001).
#   gene_bars_Gm43302_Gpsm1_Dnmt3b_IFNB_vs_Mock   one bar per gene, IFNB vs Mock (red = up, blue = down)
#   gene_bars_Gm43302_IFNB_EBOV_sets             one gene (GENE2): IFNB vs Mock | IFNB+EBOV vs Mock | EBOV vs Mock
# In the second graph a bracket joins the IFNB and IFNB+EBOV bars and gives the test of IFNB+EBOV vs IFNB: the difference between the
# two log2FCs (DESeq2 Wald test, contrast IFNB_EBOV vs IFNB_Mock, BH padj across all genes), i.e. whether EBOV changes the IFNB response of that gene.
# Writes PNG + PDF and gene_bars_Gm43302_Gpsm1_Dnmt3b_log2fc.csv to requant_v2_outputs/figures_gene_bars/
# Usage: Rscript 63_gene_bars_log2fc.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2); library(DESeq2) })
select <- dplyr::select; filter <- dplyr::filter

GENES <- c("Gm43302", "Gpsm1", "Dnmt3b")      # matched ignoring case (GPSM1 -> Gpsm1)
INK <- "#1F2933"; RED <- "#C53030"; BLUE <- "#2B6CB0"
COMP <- c(IFNB_vs_Mock = "IFNB vs Mock", IFNB_EBOV_vs_Mock = "IFNB+EBOV vs Mock", EBOV_vs_Mock = "EBOV vs Mock")
CCOL <- c(`IFNB vs Mock` = "#2B6CB0", `IFNB+EBOV vs Mock` = "#9BD08E", `EBOV vs Mock` = "#1B7837")
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_gene_bars"); dir.create(OUT, showWarnings = FALSE)
save2 <- function(p, name, w, h) { ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300); ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h) }
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
base_theme <- theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
  panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK), axis.title = element_text(colour = INK), plot.title = element_text(colour = INK, face = "bold"),
  plot.subtitle = element_text(colour = "grey35", size = 9), legend.position = "top")

D <- bind_rows(lapply(names(COMP), function(nm) {
  read_csv(file.path(RNA, paste0("DE_v2_", nm, "_annotated.csv")), show_col_types = FALSE) |>
    mutate(key = tolower(gene_name)) |> filter(key %in% tolower(GENES)) |>
    transmute(key, gene_id, gene_name, comparison = COMP[[nm]], lfc = log2FoldChange, se = lfcSE, padj)
}))
missing <- setdiff(tolower(GENES), unique(D$key))
if (length(missing)) message("Not found in the DE tables (check the symbol): ", paste(GENES[tolower(GENES) %in% missing], collapse = ", "))
disp <- setNames(GENES, tolower(GENES))
D <- D |> mutate(gene = factor(disp[key], levels = GENES), comparison = factor(comparison, levels = COMP)) |> filter(!is.na(gene)) |>
  mutate(lab = stars(padj), ypos = ifelse(lfc >= 0, lfc + se, lfc - se))
write_csv(D |> arrange(gene, comparison) |> select(gene, gene_id, comparison, log2FC = lfc, lfcSE = se, padj), file.path(OUT, "gene_bars_Gm43302_Gpsm1_Dnmt3b_log2fc.csv"))
print(as.data.frame(D |> arrange(gene, comparison) |> transmute(gene, comparison, log2FC = round(lfc, 2), padj = signif(padj, 3))), row.names = FALSE)

# test of IFNB+EBOV vs IFNB (difference between the IFNB vs Mock and IFNB+EBOV vs Mock effects)
dds <- readRDS(file.path(RNA, "dds_v2_decoyAware_hostNormalised.rds"))
rb <- as.data.frame(results(dds, contrast = c("group", "IFNB_EBOV", "IFNB_Mock"))) |> tibble::rownames_to_column("gene_id")
BR <- D |> distinct(gene, gene_id) |> left_join(rb |> select(gene_id, diff_lfc = log2FoldChange, diff_padj = padj), by = "gene_id") |> mutate(blab = stars(diff_padj))
write_csv(BR |> select(gene, gene_id, `IFNB+EBOV vs IFNB log2FC (unshrunk)` = diff_lfc, `padj` = diff_padj), file.path(OUT, "gene_bars_Gm43302_Gpsm1_Dnmt3b_IFNB_EBOV_vs_IFNB.csv"))
message("\nIFNB+EBOV vs IFNB (does EBOV change the IFNB response?):"); print(as.data.frame(BR |> transmute(gene, log2FC = round(diff_lfc, 2), padj = signif(diff_padj, 3), stars = blab)), row.names = FALSE)

# 1. IFNB vs Mock
d1 <- D |> filter(comparison == COMP[["IFNB_vs_Mock"]]) |> mutate(dir = ifelse(lfc >= 0, "up", "down"))
pad <- 0.08 * diff(range(c(d1$lfc - d1$se, d1$lfc + d1$se, 0)))
p1 <- ggplot(d1, aes(gene, lfc, fill = dir)) + geom_col(width = 0.65) +
  geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), width = 0.18, colour = INK, linewidth = 0.4) +
  geom_text(aes(y = ypos + ifelse(lfc >= 0, pad, -pad), label = lab), colour = INK, size = 4, vjust = 0.5) +
  geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) + scale_fill_manual(values = c(up = RED, down = BLUE), guide = "none") +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  labs(title = "IFNB vs Mock: log2 fold-change of selected genes", subtitle = "shrunk log2FC +/- 1 SE; stars = padj (ns = not significant)", x = NULL, y = "log2 fold-change (IFNB vs Mock)") + base_theme
save2(p1, "gene_bars_Gm43302_Gpsm1_Dnmt3b_IFNB_vs_Mock", 5, 4.5)

# 2. one gene (GENE2): IFNB, IFNB+EBOV and EBOV vs Mock; stars + padj above each bar (vs Mock), bracket = IFNB+EBOV vs IFNB
GENE2 <- "Gm43302"
d2 <- D |> filter(gene == GENE2)
if (!nrow(d2)) stop(GENE2, " not found in the DE tables")
fp <- function(p) ifelse(is.na(p), "", sprintf("%s\npadj = %s", stars(p), formatC(p, digits = 2, format = "g")))
d2 <- d2 |> mutate(x = as.numeric(comparison), blab = fp(padj))
pad2 <- 0.05 * diff(range(c(d2$lfc - d2$se, d2$lfc + d2$se, 0)))
top <- max(d2$ypos[d2$x %in% 1:2]) + 7 * pad2
br <- BR |> filter(gene == GENE2)
p2 <- ggplot(d2, aes(x, lfc, fill = comparison)) + geom_col(width = 0.65) +
  geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), width = 0.15, colour = INK, linewidth = 0.4) +
  geom_text(aes(y = ypos + ifelse(lfc >= 0, pad2, -pad2), label = blab), colour = INK, size = 3.4, lineheight = 0.9, vjust = ifelse(d2$lfc >= 0, 0, 1)) +
  annotate("segment", x = 1, xend = 2, y = top, yend = top, colour = INK, linewidth = 0.4) +
  annotate("segment", x = c(1, 2), xend = c(1, 2), y = top, yend = top - 1.5 * pad2, colour = INK, linewidth = 0.4) +
  annotate("text", x = 1.5, y = top + 0.6 * pad2, label = sprintf("IFNB+EBOV vs IFNB\n%s", fp(br$diff_padj)), colour = INK, size = 3.4, lineheight = 0.9, vjust = 0) +
  geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) + scale_fill_manual(values = CCOL, guide = "none") +
  scale_x_continuous(breaks = d2$x, labels = as.character(d2$comparison), expand = expansion(add = 0.6)) + scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  labs(title = sprintf("%s: IFNB, IFNB+EBOV and EBOV vs Mock", GENE2),
       subtitle = "shrunk log2FC +/- 1 SE; text above each bar = padj vs Mock; bracket = test of IFNB+EBOV vs IFNB (does EBOV change the IFNB response?)", x = NULL, y = "log2 fold-change vs Mock") + base_theme +
  theme(plot.subtitle = element_text(colour = "grey35", size = 8))
save2(p2, sprintf("gene_bars_%s_IFNB_EBOV_sets", GENE2), 5.5, 5.2)
message("Saved to: ", OUT)
