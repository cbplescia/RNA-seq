#!/usr/bin/env Rscript
# Bar graphs of log2FC for genes that report on TNFAIP3 (A20) regulation: is Tnfaip3 transcription altered, and is it consistent with
# (a) TLR4 -> MyD88/TRAF6 -> NF-kB signalling, or (b) the Galphai3 / cAMP / PKA / CREB axis?  Genes are grouped in panels (edit GENE_SETS):
#   Genes of interest           Gm43302, Gpsm1 (AGS3, a Galphai regulator), Dnmt3b
#   TNFAIP3 (A20)               the gene in question; it is both an NF-kB target and a CREB-responsive gene
#   TLR4 / NF-kB                Tlr4, Myd88, Traf6 (pathway components), Nfkb1, Rela (NF-kB subunits), Nfkbia (IkBa) and Tnf, Il6 (classic NF-kB target readouts)
#   Galphai3 / cAMP / PKA / CREB  Gnai3, Prkaca (PKA catalytic subunit), Creb1, and the CREB target / cAMP-feedback genes Crem, Nr4a1, Pde4b
# Read this as a TRANSCRIPT readout: induction of NF-kB targets (Nfkbia, Tnf, Il6, Tnfaip3) argues for NF-kB output; induction of Crem / Nr4a1 / Pde4b
# argues for CREB output; Creb1 / Prkaca / Gnai3 / Tlr4 mRNA changes alone do not show pathway activity (that is post-translational: phospho-CREB, cAMP, etc.).
# Shrunk log2FC from DE_v2_*_annotated.csv, error bars +/- 1 lfcSE, stars = padj (* < 0.05, ** < 0.01, *** < 0.001, ns).
#   gene_bars_tnfaip3_IFNB_vs_Mock          one bar per gene, IFNB vs Mock
#   gene_bars_tnfaip3_grouped               three bars per gene: IFNB vs Mock | IFNB+EBOV vs Mock | EBOV vs Mock
# Writes PNG + PDF and gene_bars_tnfaip3_log2fc.csv to requant_v2_outputs/figures_gene_bars/
# Usage: Rscript 64_gene_bars_tnfaip3_axis.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter

GENE_SETS <- list(`Genes of interest`              = c("Gm43302", "Gpsm1", "Dnmt3b"),
                  `TNFAIP3 (A20)`                  = c("Tnfaip3"),
                  `TLR4 / NF-kB`                   = c("Tlr4", "Myd88", "Traf6", "Nfkb1", "Rela", "Nfkbia", "Tnf", "Il6"),
                  `Galphai3 / cAMP / PKA / CREB`   = c("Gnai3", "Prkaca", "Creb1", "Crem", "Nr4a1", "Pde4b"))
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
base_theme <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
  panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK), axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_text(colour = INK),
  plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 9), legend.position = "top",
  strip.background = element_rect(fill = "grey95", colour = "grey70"), strip.text = element_text(colour = INK, face = "bold", size = 8.5), panel.spacing.x = unit(0.2, "lines"))

GT <- bind_rows(lapply(names(GENE_SETS), function(s) tibble(set = s, gene = GENE_SETS[[s]]))) |> distinct(gene, .keep_all = TRUE) |>
  mutate(key = tolower(gene), set = factor(set, levels = names(GENE_SETS)))
D <- bind_rows(lapply(names(COMP), function(nm) {
  read_csv(file.path(RNA, paste0("DE_v2_", nm, "_annotated.csv")), show_col_types = FALSE) |>
    mutate(key = tolower(gene_name)) |> filter(key %in% GT$key) |> distinct(key, .keep_all = TRUE) |>
    transmute(key, gene_id, comparison = COMP[[nm]], lfc = log2FoldChange, se = lfcSE, padj)
}))
miss <- setdiff(GT$key, unique(D$key))
if (length(miss)) message("Not found in the DE tables (check the symbol): ", paste(GT$gene[GT$key %in% miss], collapse = ", "))
D <- D |> left_join(GT |> select(key, gene, set), by = "key") |>
  mutate(gene = factor(gene, levels = GT$gene), comparison = factor(comparison, levels = COMP), lab = stars(padj), ypos = ifelse(lfc >= 0, lfc + se, lfc - se))
write_csv(D |> arrange(set, gene, comparison) |> select(set, gene, gene_id, comparison, log2FC = lfc, lfcSE = se, padj), file.path(OUT, "gene_bars_tnfaip3_log2fc.csv"))
print(as.data.frame(D |> arrange(set, gene, comparison) |> transmute(set, gene, comparison, log2FC = round(lfc, 2), padj = signif(padj, 3))), row.names = FALSE)
nG <- n_distinct(D$gene)
facet <- facet_grid(~set, scales = "free_x", space = "free_x", labeller = label_wrap_gen(14))

# 1. IFNB vs Mock
d1 <- D |> filter(comparison == COMP[["IFNB_vs_Mock"]]) |> mutate(dir = ifelse(lfc >= 0, "up", "down"))
pad <- 0.05 * diff(range(c(d1$lfc - d1$se, d1$lfc + d1$se, 0)))
p1 <- ggplot(d1, aes(gene, lfc, fill = dir)) + geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), width = 0.2, colour = INK, linewidth = 0.4) +
  geom_text(aes(y = ypos + ifelse(lfc >= 0, pad, -pad), label = lab), colour = INK, size = 3.4, vjust = 0.5) +
  geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) + facet + scale_fill_manual(values = c(up = RED, down = BLUE), guide = "none") +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  labs(title = "IFNB vs Mock: TNFAIP3, TLR4 / NF-kB and Galphai3 / cAMP / PKA / CREB readouts", subtitle = "shrunk log2FC +/- 1 SE; stars = padj (ns = not significant)", x = NULL, y = "log2 fold-change (IFNB vs Mock)") + base_theme
save2(p1, "gene_bars_tnfaip3_IFNB_vs_Mock", max(8, 0.55 * nG + 2.5), 5.2)

# 2. grouped
dg <- position_dodge(width = 0.85)
pad2 <- 0.04 * diff(range(c(D$lfc - D$se, D$lfc + D$se, 0)))
p2 <- ggplot(D, aes(gene, lfc, fill = comparison)) + geom_col(position = dg, width = 0.78) +
  geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), position = dg, width = 0.25, colour = INK, linewidth = 0.3) +
  geom_text(aes(y = ypos + ifelse(lfc >= 0, pad2, -pad2), label = lab), position = dg, colour = INK, size = 2.4, vjust = 0.5) +
  geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) + facet + scale_fill_manual(values = CCOL, name = NULL) +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  labs(title = "IFNB, IFNB+EBOV and EBOV vs Mock: TNFAIP3, TLR4 / NF-kB and Galphai3 / cAMP / PKA / CREB readouts", subtitle = "shrunk log2FC +/- 1 SE; stars = padj of each comparison against Mock (ns = not significant)", x = NULL, y = "log2 fold-change vs Mock") + base_theme
save2(p2, "gene_bars_tnfaip3_grouped", max(11, 1.0 * nG + 2.5), 5.6)
message("Saved to: ", OUT)
