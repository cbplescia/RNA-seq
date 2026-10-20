#!/usr/bin/env Rscript
# Bar graphs of log2FC for selected genes (edit GENES). Shrunk log2FC from the DESeq2 results (DE_v2_*_annotated.csv), error bars = +/- 1 lfcSE,
# stars = padj of that comparison (* < 0.05, ** < 0.01, *** < 0.001).
#   gene_bars_Gm43302_Gpsm1_Dnmt3b_IFNB_vs_Mock   one bar per gene, IFNB vs Mock (red = up, blue = down)
#   gene_bars_Gm43302_Gpsm1_Dnmt3b_grouped        per gene, three bars: IFNB vs Mock | IFNB+EBOV vs Mock | EBOV vs Mock
# Writes PNG + PDF and gene_bars_Gm43302_Gpsm1_Dnmt3b_log2fc.csv to requant_v2_outputs/figures_gene_bars/
# Usage: Rscript 63_gene_bars_log2fc.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
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

# 2. grouped: IFNB, IFNB+EBOV, EBOV vs Mock
dg <- position_dodge(width = 0.8)
pad2 <- 0.05 * diff(range(c(D$lfc - D$se, D$lfc + D$se, 0)))
p2 <- ggplot(D, aes(gene, lfc, fill = comparison)) + geom_col(position = dg, width = 0.72) +
  geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), position = dg, width = 0.2, colour = INK, linewidth = 0.35) +
  geom_text(aes(y = ypos + ifelse(lfc >= 0, pad2, -pad2), label = lab), position = dg, colour = INK, size = 3.2, vjust = 0.5) +
  geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) + scale_fill_manual(values = CCOL, name = NULL) +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  labs(title = "Selected genes: IFNB, IFNB+EBOV and EBOV vs Mock", subtitle = "shrunk log2FC +/- 1 SE; stars = padj of each comparison against Mock (ns = not significant)", x = NULL, y = "log2 fold-change vs Mock") + base_theme
save2(p2, "gene_bars_Gm43302_Gpsm1_Dnmt3b_grouped", 6.5, 4.8)
message("Saved to: ", OUT)
