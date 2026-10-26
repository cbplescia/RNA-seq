#!/usr/bin/env Rscript
# Bar graphs of the differential expression of the host endocytic / endolysosomal PROTEASES in IFNB vs Mock, EBOV vs Mock and IFNB+EBOV vs Mock (three bars per gene).
# Groups (edit GENE_SETS; my literature-based list):
#   cysteine cathepsins      Ctsb, Ctsc, Ctsf, Ctsh, Ctsk, Ctsl, Ctso, Ctss, Ctsw, Ctsz (Ctsb / Ctsl prime the EBOV glycoprotein in late endosomes)
#   aspartyl proteases       Ctsd, Ctse, Napsa, Bace1/2 and the gamma-secretase subunits Psen1/2, Psenen, Aph1a, Ncstn
#   serine proteases         Ctsa, Tpp1, Prcp, Dpp7, Scpep1 (lysosomal carboxy- / aminopeptidases)
#   asparaginyl endopeptidase Lgmn (legumain)
#   proprotein convertases and glycoprotein-processing proteases   Furin, Pcsk5/6/7, Tmprss2 (endosomal / TGN processing of viral and host precursors)
#   metallopeptidases and sheddases   Ece1, Ide, Nrd1, Thop1, Adam10, Adam17
# Shrunk log2FC from DE_v2_*_annotated.csv, error bars +/- 1 lfcSE, stars = padj of each comparison against Mock (* < 0.05, ** < 0.01, *** < 0.001, ns). Genes removed by the pre-filter
# in 04 are listed as "Not found"; genes with baseMean < MIN_BASEMEAN are left off the figure (still in the CSV); a dagger marks baseMean < LOW_EXPR (log2FC less reliable).
# Writes PNG + PDF and endocytic_proteases_log2fc.csv to requant_v2_outputs/figures_gene_bars/
# Usage: Rscript 85_gene_bars_endocytic_proteases.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter

GENE_SETS <- list(`Cysteine cathepsins`              = c("Ctsb", "Ctsc", "Ctsf", "Ctsh", "Ctsk", "Ctsl", "Ctso", "Ctss", "Ctsw", "Ctsz"),
                  `Aspartyl proteases / gamma-secretase` = c("Ctsd", "Ctse", "Napsa", "Bace1", "Bace2", "Psen1", "Psen2", "Psenen", "Aph1a", "Ncstn"),
                  `Serine proteases`                  = c("Ctsa", "Tpp1", "Prcp", "Dpp7", "Scpep1"),
                  `Asparaginyl endopeptidase`         = c("Lgmn"),
                  `Proprotein convertases / TMPRSS`   = c("Furin", "Pcsk5", "Pcsk6", "Pcsk7", "Tmprss2"),
                  `Metallopeptidases and sheddases`   = c("Ece1", "Ide", "Nrd1", "Thop1", "Adam10", "Adam17"))
LOW_EXPR <- 50                                          # baseMean below this: dagger on the label (log2FC less reliable)
MIN_BASEMEAN <- 10                                      # genes below this are not drawn (still in the CSV)
INK <- "#1F2933"; RED <- "#C53030"; BLUE <- "#2B6CB0"
COMP <- c(IFNB_vs_Mock = "IFNB vs Mock", EBOV_vs_Mock = "EBOV vs Mock", IFNB_EBOV_vs_Mock = "IFNB+EBOV vs Mock")
CCOL <- c(`IFNB vs Mock` = "#2B6CB0", `EBOV vs Mock` = "#1B7837", `IFNB+EBOV vs Mock` = "#9BD08E")
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
write_csv(D |> arrange(set, gene, comparison) |> select(set, gene, gene_id, baseMean, comparison, log2FC = lfc, lfcSE = se, padj), file.path(OUT, "endocytic_proteases_log2fc.csv"))
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

# grouped: IFNB, EBOV and IFNB+EBOV vs Mock
dg <- position_dodge(width = 0.85)
pad2 <- 0.04 * diff(range(c(D$lfc - D$se, D$lfc + D$se, 0)))
p2 <- ggplot(D, aes(gene, lfc, fill = comparison)) + geom_col(position = dg, width = 0.78) +
  geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), position = dg, width = 0.25, colour = INK, linewidth = 0.3) +
  geom_text(aes(y = ypos + ifelse(lfc >= 0, pad2, -pad2), label = lab), position = dg, colour = INK, size = 2.4, vjust = 0.5) +
  geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) + facet + scale_fill_manual(values = CCOL, name = NULL) +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  labs(title = "IFNB, EBOV and IFNB+EBOV vs Mock: endocytic / endolysosomal proteases", subtitle = "shrunk log2FC +/- 1 SE; stars = padj of each comparison against Mock (ns = not significant); dagger = low expression (baseMean < 50)", x = NULL, y = "log2 fold-change vs Mock") + base_theme
save2(p2, "gene_bars_endocytic_proteases", min(30, max(11, 0.75 * nG + 2.5)), 5.8)
message("Saved to: ", OUT)
