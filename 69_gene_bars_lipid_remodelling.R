#!/usr/bin/env Rscript
# Bar graphs of differential expression for host genes that remodel lipid bilayers, especially endocytic / endolysosomal membranes (edit GENE_SETS):
#   Flippases (P4-ATPases)      Atp8a1/2, Atp8b1-4, Atp9a/b, Atp10a/b/d, Atp11a/b/c  - move PS / PE to the cytosolic leaflet
#   Flippase subunits           Tmem30a/b/c (CDC50A/B/C), required for P4-ATPase trafficking and activity
#   Scramblases                 Plscr1-4, Ano6 (TMEM16F) and other anoctamins, Xkr4-9, Tmem41b, Vmp1, Atg9a  - dissipate leaflet asymmetry
#   Floppases (ABC)             Abca1/3/7, Abcb1a/1b/4, Abcc1, Abcg1/2  - move lipids to the exoplasmic / luminal leaflet
#   Endolysosomal lipid handling Npc1, Npc2, Lipa, Smpd1 (acid sphingomyelinase), Asah1, Pikfyve, Pla2g15, Cln5, Pld3  (cholesterol, sphingolipid, PI(3,5)P2, BMP)
#   PS receptors / entry        Axl, Mertk, Tyro3, Gas6, Pros1, Havcr1 (TIM-1), Timd4, Havcr2  - read PS exposed on virions (apoptotic mimicry)
# Genes removed by the pre-filter in 04 (no counts) cannot be plotted; they are listed as "Not found". Genes with baseMean < MIN_BASEMEAN are left off the figures
# (too few counts for a reliable log2FC) but stay in the CSV. Shrunk log2FC from DE_v2_*_annotated.csv, error bars +/- 1 lfcSE, stars = padj (* < 0.05, ** < 0.01, *** < 0.001, ns).
#   gene_bars_lipid_remodelling_IFNB_vs_Mock          one bar per gene, IFNB vs Mock
#   gene_bars_lipid_remodelling_grouped               three bars per gene: IFNB vs Mock | IFNB+EBOV vs Mock | EBOV vs Mock
# Writes PNG + PDF and gene_bars_lipid_remodelling_log2fc.csv to requant_v2_outputs/figures_gene_bars/
# Usage: Rscript 69_gene_bars_lipid_remodelling.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter

GENE_SETS <- list(`Flippases (P4-ATPases)`        = c("Atp8a1", "Atp8a2", "Atp8b1", "Atp8b2", "Atp8b3", "Atp8b4", "Atp9a", "Atp9b", "Atp10a", "Atp10b", "Atp10d", "Atp11a", "Atp11b", "Atp11c"),
                  `Flippase subunits`             = c("Tmem30a", "Tmem30b", "Tmem30c"),
                  `Scramblases`                   = c("Plscr1", "Plscr2", "Plscr3", "Plscr4", "Ano6", "Ano10", "Ano3", "Ano4", "Ano5", "Ano7", "Ano8", "Ano9", "Xkr4", "Xkr5", "Xkr6", "Xkr7", "Xkr8", "Xkr9", "Tmem41b", "Vmp1", "Atg9a"),
                  `Floppases (ABC transporters)`  = c("Abca1", "Abca3", "Abca7", "Abcb1a", "Abcb1b", "Abcb4", "Abcc1", "Abcg1", "Abcg2"),
                  `Endolysosomal lipid handling`  = c("Npc1", "Npc2", "Lipa", "Smpd1", "Asah1", "Pikfyve", "Pla2g15", "Cln5", "Pld3"),
                  `PS receptors / viral entry`    = c("Axl", "Mertk", "Tyro3", "Gas6", "Pros1", "Havcr1", "Timd4", "Havcr2"))
LOW_EXPR <- 50                                          # baseMean below this: dagger on the label (log2FC less reliable)
MIN_BASEMEAN <- 10                                      # genes below this are not drawn (still in the CSV)
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
    transmute(key, gene_id, baseMean, comparison = COMP[[nm]], lfc = log2FoldChange, se = lfcSE, padj)
}))
miss <- setdiff(GT$key, unique(D$key))
if (length(miss)) message("Not found in the DE tables (check the symbol): ", paste(GT$gene[GT$key %in% miss], collapse = ", "))
D <- D |> left_join(GT |> select(key, gene, set), by = "key") |>
  mutate(gene = factor(gene, levels = GT$gene), comparison = factor(comparison, levels = COMP), lab = stars(padj), ypos = ifelse(lfc >= 0, lfc + se, lfc - se))
write_csv(D |> arrange(set, gene, comparison) |> select(set, gene, gene_id, baseMean, comparison, log2FC = lfc, lfcSE = se, padj), file.path(OUT, "gene_bars_lipid_remodelling_log2fc.csv"))
print(as.data.frame(D |> arrange(set, gene, comparison) |> transmute(set, gene, comparison, log2FC = round(lfc, 2), padj = signif(padj, 3))), row.names = FALSE)
nOff <- D |> group_by(gene) |> summarise(bm = max(baseMean, na.rm = TRUE), .groups = "drop") |> filter(bm < MIN_BASEMEAN)
if (nrow(nOff)) message("Left off the figures (baseMean < ", MIN_BASEMEAN, "; still in the CSV): ", paste(nOff$gene, collapse = ", "))
D <- D |> filter(!gene %in% nOff$gene) |> mutate(gene = droplevels(gene), set = droplevels(set))
if (!nrow(D)) stop("No genes left to plot")
nG <- n_distinct(D$gene)
bm <- D |> group_by(gene) |> summarise(bm = max(baseMean, na.rm = TRUE), .groups = "drop")
xlab <- setNames(paste0(bm$gene, ifelse(bm$bm < LOW_EXPR, "\u2020", "")), bm$gene)
if (any(bm$bm < LOW_EXPR)) message("Low expression (baseMean < ", LOW_EXPR, ", marked with a dagger; treat log2FC cautiously): ", paste(bm$gene[bm$bm < LOW_EXPR], collapse = ", "))
facet <- list(facet_grid(~set, scales = "free_x", space = "free_x", labeller = label_wrap_gen(16)), scale_x_discrete(labels = xlab))

# 1. IFNB vs Mock
d1 <- D |> filter(comparison == COMP[["IFNB_vs_Mock"]]) |> mutate(dir = ifelse(lfc >= 0, "up", "down"))
pad <- 0.05 * diff(range(c(d1$lfc - d1$se, d1$lfc + d1$se, 0)))
p1 <- ggplot(d1, aes(gene, lfc, fill = dir)) + geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), width = 0.2, colour = INK, linewidth = 0.4) +
  geom_text(aes(y = ypos + ifelse(lfc >= 0, pad, -pad), label = lab), colour = INK, size = 3.4, vjust = 0.5) +
  geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) + facet + scale_fill_manual(values = c(up = RED, down = BLUE), guide = "none") +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  labs(title = "IFNB vs Mock: genes that remodel host lipid bilayers (flippases, scramblases, floppases, endolysosomal lipid handling)", subtitle = "shrunk log2FC +/- 1 SE; stars = padj (ns = not significant); dagger = low expression (baseMean < 50)", x = NULL, y = "log2 fold-change (IFNB vs Mock)") + base_theme
save2(p1, "gene_bars_lipid_remodelling_IFNB_vs_Mock", min(26, max(8, 0.5 * nG + 2.5)), 5.4)

# 2. grouped
dg <- position_dodge(width = 0.85)
pad2 <- 0.04 * diff(range(c(D$lfc - D$se, D$lfc + D$se, 0)))
p2 <- ggplot(D, aes(gene, lfc, fill = comparison)) + geom_col(position = dg, width = 0.78) +
  geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), position = dg, width = 0.25, colour = INK, linewidth = 0.3) +
  geom_text(aes(y = ypos + ifelse(lfc >= 0, pad2, -pad2), label = lab), position = dg, colour = INK, size = 2.4, vjust = 0.5) +
  geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) + facet + scale_fill_manual(values = CCOL, name = NULL) +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  labs(title = "IFNB, IFNB+EBOV and EBOV vs Mock: genes that remodel host lipid bilayers", subtitle = "shrunk log2FC +/- 1 SE; stars = padj of each comparison against Mock (ns = not significant); dagger = low expression (baseMean < 50)", x = NULL, y = "log2 fold-change vs Mock") + base_theme
save2(p2, "gene_bars_lipid_remodelling_grouped", min(30, max(11, 0.75 * nG + 2.5)), 5.8)
message("Saved to: ", OUT)
