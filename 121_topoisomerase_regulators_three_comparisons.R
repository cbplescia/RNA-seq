#!/usr/bin/env Rscript
# DNA topoisomerases and the genes that can modify their activity, function or stability, three bars per gene, left to right IFNB vs Mock (blue), EBOV vs Mock (dark green), IFNB+EBOV vs Mock (light green).
# Genes are grouped in panels by how they act on topoisomerases (CLASSES below; edit freely - this is a curated list from the literature, NOT an exhaustive or database-derived one, and many genes act in several ways):
#   topoisomerases | TOP3 partners (BLM-TOP3A-RMI, TOP3B-TDRD3-FMR1) | repair of TOP1 / TOP2 cleavage complexes (TDP1, TDP2, PARP1/XRCC1, MRE11-RAD50-NBN-CtIP) | DNA-protein-crosslink proteolysis and extraction (SPRTN, TEX264, GCNA, VCP) |
#   SUMO / ubiquitin control (UBE2I, SUMO, PIAS, TOPORS, RNF4, ZNF451, SENP1) | kinases (CK2, CDK1, PLK1, AURKB) | chromatin / condensin partners (SMC2, SMC4, KIF4A) | other TOP1/TOP2 interactors (NCL, NPM1, WRN, TP53, BRCA1, HMGB1)
# Shrunk log2FC from DE_v2_*_annotated.csv, error bars +/- 1 lfcSE, stars = padj of each comparison against Mock (* < 0.05, ** < 0.01, *** < 0.001, ns), dagger on the gene name = baseMean < 50. Genes missing from the tables (removed by the pre-filter in 04, or a
# different symbol) are listed as "Not found". Symbols are matched ignoring case; old / alternative names are in ALIAS.
# Env GENES="Label=Gene1,Gene2;Label2=Gene3" (any gene sets), OUT_NAME, TITLE.   Writes bars_<OUT_NAME>_pageN.png/.pdf and topoisomerase_regulators_three_comparisons.csv to requant_v2_outputs/figures_gene_bars/.   Env: PER_PAGE (default 26 genes per page), ONLY_SIG=1 (genes with padj < 0.05 in at least one comparison)
# Usage: Rscript 121_topoisomerase_regulators_three_comparisons.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter

CLASSES <- list(
  "Topoisomerases"                              = c("Top1", "Top1mt", "Top2a", "Top2b", "Top3a", "Top3b"),
  "TOP3 partners"                               = c("Blm", "Rmi1", "Rmi2", "Tdrd3", "Fmr1"),
  "Repair of TOP1 / TOP2 cleavage complexes"    = c("Tdp1", "Tdp2", "Parp1", "Xrcc1", "Lig3", "Pnkp", "Aptx", "Mre11a", "Rad50", "Nbn", "Rbbp8"),
  "DNA-protein crosslink proteolysis / extraction" = c("Sprtn", "Tex264", "Gcna", "Vcp"),
  "SUMO / ubiquitin control"                    = c("Ube2i", "Sumo1", "Sumo2", "Sumo3", "Pias1", "Pias4", "Topors", "Rnf4", "Zfp451", "Senp1"),
  "Kinases"                                     = c("Csnk2a1", "Csnk2a2", "Csnk2b", "Cdk1", "Plk1", "Aurkb"),
  "Chromatin / condensin partners"              = c("Smc2", "Smc4", "Kif4a"),
  "Other TOP1 / TOP2 interactors"               = c("Ncl", "Npm1", "Wrn", "Trp53", "Brca1", "Hmgb1"))
if (nzchar(Sys.getenv("GENES"))) {   # GENES="Label=Gene1,Gene2;Label2=Gene3" replaces the list above with any gene sets (one panel per label)
  parts <- strsplit(strsplit(Sys.getenv("GENES"), ";")[[1]], "=")
  CLASSES <- setNames(lapply(parts, function(z) trimws(strsplit(z[2], ",")[[1]])), trimws(vapply(parts, `[`, "", 1)))
}
OUT_NAME <- Sys.getenv("OUT_NAME", "topoisomerase_regulators"); TITLE <- Sys.getenv("TITLE", "Topoisomerases and genes that modify their activity")
ALIAS <- list(Rbl1 = "p107", Rbl2 = "p130", Mre11a = "Mre11", Zfp451 = "Znf451", Kif4a = "Kif4", Trp53 = "Tp53", Csnk2b = "Ck2b")
LOW_EXPR <- 50
COMP <- c(IFNB_vs_Mock = "IFNB vs Mock", EBOV_vs_Mock = "EBOV vs Mock", IFNB_EBOV_vs_Mock = "IFNB+EBOV vs Mock"); CCOL <- c(`IFNB vs Mock` = "#2B6CB0", `EBOV vs Mock` = "#1B7837", `IFNB+EBOV vs Mock` = "#9BD08E")
INK <- "#1F2933"; PER_PAGE <- as.integer(Sys.getenv("PER_PAGE", "26")); ONLY_SIG <- Sys.getenv("ONLY_SIG", "0") == "1"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs"); OUT <- file.path(RNA, "figures_gene_bars"); dir.create(OUT, showWarnings = FALSE)
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

tabs <- lapply(names(COMP), function(nm) read_csv(file.path(RNA, paste0("DE_v2_", nm, "_annotated.csv")), show_col_types = FALSE) |> mutate(key = tolower(gene_name))); names(tabs) <- names(COMP)
find_gene <- function(q) {
  syms <- tolower(if (!is.null(ALIAS[[q]])) c(q, ALIAS[[q]]) else q)
  h <- tabs[[1]] |> filter(key %in% syms | tolower(sub("\\..*$", "", gene_id)) %in% syms) |> distinct(gene_id, .keep_all = TRUE)
  if (nrow(h)) h$gene_id[1] else NA_character_
}
want <- bind_rows(lapply(names(CLASSES), function(g) tibble(class = g, query = CLASSES[[g]]))) |> mutate(gene_id = vapply(query, find_gene, ""))
miss <- want$query[is.na(want$gene_id)]; if (length(miss)) message("Not found in the DE tables (removed by the pre-filter, or a different symbol): ", paste(miss, collapse = ", "))
want <- want |> filter(!is.na(gene_id)) |> distinct(gene_id, .keep_all = TRUE) |> mutate(class = factor(class, levels = names(CLASSES))) |> arrange(class) |> mutate(pos = row_number())
D <- bind_rows(lapply(names(COMP), function(nm) tabs[[nm]] |> inner_join(want, by = "gene_id") |> transmute(class, query, gene = gene_name, gene_id, pos, baseMean, comparison = COMP[[nm]], lfc = log2FoldChange, se = lfcSE, padj)))
if (ONLY_SIG) { keep <- D |> group_by(gene_id) |> summarise(s = any(!is.na(padj) & padj < 0.05), .groups = "drop") |> filter(s) |> pull(gene_id); D <- D |> filter(gene_id %in% keep) }
D <- D |> mutate(comparison = factor(comparison, levels = COMP), lab = stars(padj), vj = ifelse(lfc >= 0, 0, 1), ys = lfc + ifelse(lfc >= 0, 1, -1) * (se + 0.03),
                 shown = ifelse(is.na(gene) | gene == "", query, gene), xlab = paste0(shown, ifelse(baseMean < LOW_EXPR, "†", "")))
if (!nrow(D)) stop("No genes to plot")
# pages: whole classes until PER_PAGE genes would be exceeded
gsz <- D |> distinct(class, gene_id) |> dplyr::count(class) |> arrange(class); pg <- 1; used <- 0; gpage <- c()
for (i in seq_len(nrow(gsz))) { if (used > 0 && used + gsz$n[i] > PER_PAGE) { pg <- pg + 1; used <- 0 }; gpage[as.character(gsz$class[i])] <- pg; used <- used + gsz$n[i] }
D$page <- unname(gpage[as.character(D$class)])
print(as.data.frame(D |> transmute(class, gene = shown, baseMean = round(baseMean, 1), comparison, log2FC = round(lfc, 2), padj = signif(padj, 2)) |> tidyr::pivot_wider(names_from = comparison, values_from = c(log2FC, padj))), row.names = FALSE)
write_csv(D |> select(class, gene = shown, gene_id, baseMean, comparison, log2FC = lfc, lfcSE = se, padj, page), file.path(OUT, paste0(OUT_NAME, "_three_comparisons.csv")))
npg <- max(D$page); dg <- position_dodge(width = 0.85)
for (pgi in seq_len(npg)) {
  Dp <- D |> filter(page == pgi); Dp$xlab <- factor(Dp$xlab, levels = unique(Dp$xlab[order(Dp$pos)]))
  p <- ggplot(Dp, aes(xlab, lfc, fill = comparison)) + geom_col(position = dg, width = 0.8) +
    geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), position = dg, width = 0.25, colour = INK, linewidth = 0.35) +
    geom_text(aes(y = ys, label = lab, vjust = vj), position = dg, colour = INK, size = 2.6) + geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) +
    facet_grid(. ~ class, scales = "free_x", space = "free_x", labeller = labeller(class = label_wrap_gen(20))) +
    scale_fill_manual(values = CCOL, name = NULL) + scale_y_continuous(expand = expansion(mult = 0.12)) +
    labs(title = paste0(TITLE, ": IFNB, EBOV and IFNB+EBOV vs Mock", if (npg > 1) paste0("  [page ", pgi, " of ", npg, "]") else ""),
         subtitle = "left to right: IFNB, EBOV, IFNB+EBOV vs Mock; shrunk log2FC +/- 1 SE; stars = padj of each comparison (ns = not significant); dagger = baseMean < 50", x = NULL, y = "log2 fold-change vs Mock") +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3), panel.border = element_rect(colour = "grey70"),
          axis.text = element_text(colour = INK), axis.text.x = element_text(face = "italic", angle = 45, hjust = 1), plot.title = element_text(colour = INK, face = "bold", size = 11), plot.subtitle = element_text(colour = "grey35", size = 8.5),
          legend.position = "top", strip.background = element_rect(fill = "grey93", colour = "grey70"), strip.text = element_text(colour = INK, face = "bold", size = 8))
  w <- max(7, 0.85 * nlevels(Dp$xlab) + 3)
  ggsave(file.path(OUT, paste0("bars_", OUT_NAME, "_page", pgi, ".png")), p, width = w, height = 5.4, dpi = 300); ggsave(file.path(OUT, paste0("bars_", OUT_NAME, "_page", pgi, ".pdf")), p, width = w, height = 5.4)
}
message("\nSaved ", npg, " page(s) to: ", OUT)
