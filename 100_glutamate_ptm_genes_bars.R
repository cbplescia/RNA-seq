#!/usr/bin/env Rscript
# Bar graphs (IFNB vs Mock next to EBOV vs Mock, vertical) for genes that modify glutamic acid residues in proteins. Panels (my classification; edit PANELS):
#   gamma-carboxylation (Gla)           Ggcx (the carboxylase) and the vitamin K cycle that feeds it, Vkorc1, Vkorc1l1
#   polyglutamylation (glutamate chains) Ttll1/2/4/5/6/7/9/11/13 (tubulin tyrosine ligase-like glutamylases)
#   glycylation (glycine on Glu)         Ttll3, Ttll8
#   deglutamylation / deglycylation      Agtpbp1 (CCP1), Agbl1-5 (cytosolic carboxypeptidases)
#   ADP-ribosylation of Glu / Asp        writers Parp1, Parp2 (+ the switch factor Hpf1), erasers Parg, Targ1 (Oarsd1), Macrod1, Macrod2, Adprh, Adprhl2
#   N-terminal pyroglutamate             Qpct, Qpctl (glutaminyl cyclases; they also convert N-terminal glutamate)
#   other, from GO                       every other gene annotated to "peptidyl-glutamic acid modification" (GO:0018200, children included) in org.Mm.eg.db, if clusterProfiler / org.Mm.eg.db are installed
# Shrunk log2FC from DE_v2_IFNB_vs_Mock / DE_v2_EBOV_vs_Mock _annotated.csv, +/- 1 lfcSE, stars = padj against Mock (* < 0.05, ** < 0.01, *** < 0.001, ns). A dagger marks baseMean < LOW_EXPR.
# Genes missing from the tables (removed by the pre-filter in 04, or a different symbol) are listed as "Not found"; they cannot be plotted. All plotted genes are in glutamate_ptm_genes_log2fc.csv.
# Writes bars_glutamate_ptm_genes.png / .pdf and the csv to requant_v2_outputs/figures_gene_bars/.   Usage: Rscript 100_glutamate_ptm_genes_bars.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter

PANELS <- list(`gamma-carboxylation` = c("Ggcx", "Vkorc1", "Vkorc1l1"),
               `polyglutamylation` = c("Ttll1", "Ttll2", "Ttll4", "Ttll5", "Ttll6", "Ttll7", "Ttll9", "Ttll11", "Ttll13"),
               `glycylation (Ttll3/8)` = c("Ttll3", "Ttll8"),
               `deglutamylation (CCPs)` = c("Agtpbp1", "Agbl1", "Agbl2", "Agbl3", "Agbl4", "Agbl5"),
               `ADP-ribosylation of Glu / Asp` = c("Parp1", "Parp2", "Hpf1", "Parg", "Targ1", "Macrod1", "Macrod2", "Adprh", "Adprhl2"),
               `N-terminal pyroglutamate` = c("Qpct", "Qpctl"))
ALIAS <- list(TARG1 = c("Targ1", "Oarsd1", "C6orf130"), HPF1 = c("Hpf1", "C4orf27", "1700123O20Rik"), AGTPBP1 = c("Agtpbp1", "Nna1", "Ccp1"))
LOW_EXPR <- 50
COMP <- c(IFNB_vs_Mock = "IFNB vs Mock", EBOV_vs_Mock = "EBOV vs Mock"); CCOL <- c(`IFNB vs Mock` = "#2B6CB0", `EBOV vs Mock` = "#1B7837"); INK <- "#1F2933"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_gene_bars"); dir.create(OUT, showWarnings = FALSE)
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
tabs <- lapply(names(COMP), function(nm) read_csv(file.path(RNA, paste0("DE_v2_", nm, "_annotated.csv")), show_col_types = FALSE) |> mutate(key = tolower(gene_name), gid = sub("\\..*$", "", gene_id))); names(tabs) <- names(COMP)

# ---- extra genes from GO "peptidyl-glutamic acid modification" ----
if (requireNamespace("clusterProfiler", quietly = TRUE) && requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
  g2a <- as.list(org.Mm.eg.db::org.Mm.egGO2ALLEGS)
  e <- unique(as.character(unlist(g2a[c("GO:0018200", "GO:0018095", "GO:0018094", "GO:0017187")]))); e <- e[!is.na(e)]
  go_ens <- if (length(e)) unique(suppressMessages(suppressWarnings(clusterProfiler::bitr(e, "ENTREZID", "ENSEMBL", org.Mm.eg.db::org.Mm.eg.db)))$ENSEMBL) else character()
  curated <- tolower(unlist(c(PANELS, ALIAS)))
  extra <- tabs[[1]] |> filter(gid %in% go_ens, !key %in% curated) |> distinct(gene_id, .keep_all = TRUE) |> pull(gene_name)
  message(sprintf("GO peptidyl-glutamic acid modification: %d genes in the DE tables, %d not in my curated panels: %s", sum(tabs[[1]]$gid %in% go_ens), length(extra), paste(extra, collapse = ", ")))
  if (length(extra)) PANELS[["other (GO: glutamic acid modification)"]] <- extra
} else message("clusterProfiler / org.Mm.eg.db not installed - only the curated panels are drawn")

find_gene <- function(q) {
  syms <- tolower(if (!is.null(ALIAS[[toupper(q)]])) c(q, ALIAS[[toupper(q)]]) else q)
  h <- tabs[[1]] |> filter(key %in% syms | tolower(gid) %in% syms) |> distinct(gene_id, .keep_all = TRUE)
  if (nrow(h)) h$gene_id[1] else NA_character_
}
GT <- bind_rows(lapply(names(PANELS), function(p) tibble(panel = p, query = PANELS[[p]]))) |> distinct(query, .keep_all = TRUE); GT$gene_id <- vapply(GT$query, find_gene, "")
miss <- GT$query[is.na(GT$gene_id)]; if (length(miss)) message("Not found in the DE tables (removed by the pre-filter, or a different symbol): ", paste(miss, collapse = ", "))
GT <- GT |> filter(!is.na(gene_id)) |> mutate(pos = row_number())
if (!nrow(GT)) stop("None of the genes were found")
D <- bind_rows(lapply(names(COMP), function(nm) tabs[[nm]] |> inner_join(GT, by = "gene_id") |> transmute(panel, query, gene = gene_name, gene_id, pos, baseMean, comparison = COMP[[nm]], lfc = log2FoldChange, se = lfcSE, padj)))
D <- D |> mutate(panel = factor(panel, levels = names(PANELS)), comparison = factor(comparison, levels = COMP), lab = stars(padj), vj = ifelse(lfc >= 0, 0, 1), ys = lfc + ifelse(lfc >= 0, 1, -1) * (se + 0.03),
                 xl = paste0(gene, ifelse(baseMean < LOW_EXPR, "†", "")), key = paste(panel, xl, sep = "||"))
D$key <- factor(D$key, levels = unique(D$key[order(D$pos)])); D$panel <- droplevels(D$panel)
write_csv(D |> arrange(pos, comparison) |> select(panel, gene, gene_id, baseMean, comparison, log2FC = lfc, lfcSE = se, padj), file.path(OUT, "glutamate_ptm_genes_log2fc.csv"))
print(as.data.frame(D |> transmute(gene, baseMean = round(baseMean, 1), comparison, log2FC = round(lfc, 2), padj = signif(padj, 2)) |> pivot_wider(names_from = comparison, values_from = c(log2FC, padj))), row.names = FALSE)

dg <- position_dodge(width = 0.85)
p <- ggplot(D, aes(key, lfc, fill = comparison)) + geom_col(position = dg, width = 0.8) + geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), position = dg, width = 0.25, colour = INK, linewidth = 0.3) +
  geom_text(aes(y = ys, label = lab, vjust = vj), position = dg, colour = INK, size = 2.3) + geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) +
  facet_grid(. ~ panel, scales = "free_x", space = "free_x", labeller = labeller(panel = label_wrap_gen(14))) + scale_fill_manual(values = CCOL, name = NULL) +
  scale_x_discrete(labels = function(x) sub("^.*\\|\\|", "", x)) + scale_y_continuous(expand = expansion(mult = 0.12)) +
  labs(title = "Genes that modify glutamic acid residues: IFNB vs Mock and EBOV vs Mock", subtitle = sprintf("shrunk log2FC +/- 1 SE; stars = padj of each comparison (ns = not significant); dagger = baseMean < %d; %d genes found, %d not in the data", LOW_EXPR, nrow(GT), length(miss)),
       x = NULL, y = "log2 fold-change vs Mock") + theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3), panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK),
        axis.text.x = element_text(face = "italic", angle = 50, hjust = 1, size = 8), plot.title = element_text(colour = INK, face = "bold", size = 11), plot.subtitle = element_text(colour = "grey35", size = 8),
        strip.background = element_rect(fill = "grey95", colour = "grey70"), strip.text = element_text(colour = INK, face = "bold", size = 7.5), legend.position = "top", panel.spacing.x = unit(0.2, "lines"))
w <- max(7, 0.5 * nlevels(D$key) + 1.3 * nlevels(D$panel) + 2)
ggsave(file.path(OUT, "bars_glutamate_ptm_genes.png"), p, width = w, height = 5.4, dpi = 300, limitsize = FALSE); ggsave(file.path(OUT, "bars_glutamate_ptm_genes.pdf"), p, width = w, height = 5.4, limitsize = FALSE)
message("Saved to: ", OUT)
