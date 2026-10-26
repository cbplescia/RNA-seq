#!/usr/bin/env Rscript
# Vertical bar graphs: for each gene THREE bars side by side, left to right IFNB vs Mock (blue), EBOV vs Mock (dark green), IFNB+EBOV vs Mock (light green), genes in the order you list them (edit SETS):
#   set 1  App, Ctsb, Ctsd, Ctsl, Psap, Sort1, Plau, Cst3, Klc4, Lgmn, Apoa5, Apoc2, Klc3, Aph1b, Serpinf2, Nts, Itm2b, Fah, Ctsf, Mmel1, Fpr2
#   set 2  Suox, Sptlc1, Tst, Smpd2, Aox1, Aldh2, Cbs, Asah2, Aldh1b1, Cerk, Akr1b7, Sqor, Gm2a, Hexa, Pdxp, Glb1, Mocs3, Limk1
# Shrunk log2FC from DE_v2_*_annotated.csv, error bars +/- 1 lfcSE, stars = padj of each comparison against Mock (* < 0.05, ** < 0.01, *** < 0.001, ns). A dagger on a gene name marks baseMean < LOW_EXPR
# (few reads, so the log2FC is less reliable); genes missing from the tables (removed by the pre-filter in 04, or a different symbol) are listed as "Not found". Mouse symbols are matched ignoring case.
# Writes bars_<set>.png / .pdf and gene_bars_three_comparisons.csv to requant_v2_outputs/figures_gene_bars/.   Usage: Rscript 94_gene_bars_three_comparisons.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter

SETS <- list(gene_set_1 = list(title = "Gene set 1: IFNB, EBOV and IFNB+EBOV vs Mock", genes = c("App", "Ctsb", "Ctsd", "Ctsl", "Psap", "Sort1", "Plau", "Cst3", "Klc4", "Lgmn", "Apoa5", "Apoc2", "Klc3", "Aph1b", "Serpinf2", "Nts", "Itm2b", "Fah", "Ctsf", "Mmel1", "Fpr2")),
             gene_set_2 = list(title = "Gene set 2: IFNB, EBOV and IFNB+EBOV vs Mock", genes = c("Suox", "Sptlc1", "Tst", "Smpd2", "Aox1", "Aldh2", "Cbs", "Asah2", "Aldh1b1", "Cerk", "Akr1b7", "Sqor", "Gm2a", "Hexa", "Pdxp", "Glb1", "Mocs3", "Limk1")))
ALIAS <- list()
LOW_EXPR <- 50
COMP <- c(IFNB_vs_Mock = "IFNB vs Mock", EBOV_vs_Mock = "EBOV vs Mock", IFNB_EBOV_vs_Mock = "IFNB+EBOV vs Mock"); CCOL <- c(`IFNB vs Mock` = "#2B6CB0", `EBOV vs Mock` = "#1B7837", `IFNB+EBOV vs Mock` = "#9BD08E")
INK <- "#1F2933"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_gene_bars"); dir.create(OUT, showWarnings = FALSE)
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

tabs <- lapply(names(COMP), function(nm) read_csv(file.path(RNA, paste0("DE_v2_", nm, "_annotated.csv")), show_col_types = FALSE) |> mutate(key = tolower(gene_name))); names(tabs) <- names(COMP)
find_gene <- function(q) {
  syms <- tolower(if (!is.null(ALIAS[[q]])) c(q, ALIAS[[q]]) else q)
  h <- tabs[[1]] |> filter(key %in% syms | tolower(sub("\\..*$", "", gene_id)) %in% syms) |> distinct(gene_id, .keep_all = TRUE)
  if (nrow(h)) h$gene_id[1] else NA_character_
}
allrows <- list()
for (s in names(SETS)) {
  S <- SETS[[s]]; gids <- vapply(S$genes, find_gene, ""); miss <- S$genes[is.na(gids)]
  message("\n== ", S$title); if (length(miss)) message("   Not found in the DE tables (removed by the pre-filter, or a different symbol): ", paste(miss, collapse = ", "))
  keep <- !is.na(gids); if (!any(keep)) next
  ord <- tibble(query = S$genes[keep], gene_id = gids[keep], pos = seq_len(sum(keep)))
  D <- bind_rows(lapply(names(COMP), function(nm) tabs[[nm]] |> inner_join(ord, by = "gene_id") |> transmute(query, gene = gene_name, gene_id, pos, baseMean, comparison = COMP[[nm]], lfc = log2FoldChange, se = lfcSE, padj)))
  D <- D |> mutate(comparison = factor(comparison, levels = COMP), lab = stars(padj), vj = ifelse(lfc >= 0, 0, 1), ys = lfc + ifelse(lfc >= 0, 1, -1) * (se + 0.03),
                   shown = ifelse(is.na(gene) | gene == "", query, gene), xlab = paste0(shown, ifelse(baseMean < LOW_EXPR, "†", "")))
  D$xlab <- factor(D$xlab, levels = unique(D$xlab[order(D$pos)]))
  allrows[[s]] <- D |> mutate(set = s)
  print(as.data.frame(D |> transmute(gene = shown, baseMean = round(baseMean, 1), comparison, log2FC = round(lfc, 2), padj = signif(padj, 2)) |> tidyr::pivot_wider(names_from = comparison, values_from = c(log2FC, padj))), row.names = FALSE)
  dg <- position_dodge(width = 0.85)
  p <- ggplot(D, aes(xlab, lfc, fill = comparison)) + geom_col(position = dg, width = 0.8) +
    geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), position = dg, width = 0.25, colour = INK, linewidth = 0.35) +
    geom_text(aes(y = ys, label = lab, vjust = vj), position = dg, colour = INK, size = 2.8) + geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) +
    scale_fill_manual(values = CCOL, name = NULL) + scale_y_continuous(expand = expansion(mult = 0.12)) +
    labs(title = S$title, subtitle = "left to right: IFNB, EBOV, IFNB+EBOV vs Mock; shrunk log2FC +/- 1 SE; stars = padj of each comparison (ns = not significant); dagger = baseMean < 50", x = NULL, y = "log2 fold-change vs Mock") +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3), panel.border = element_rect(colour = "grey70"),
          axis.text = element_text(colour = INK), axis.text.x = element_text(face = "italic", angle = 45, hjust = 1), plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 8.5), legend.position = "top")
  w <- max(6, 0.9 * nlevels(D$xlab) + 2)
  ggsave(file.path(OUT, paste0("bars_", s, ".png")), p, width = w, height = 5, dpi = 300); ggsave(file.path(OUT, paste0("bars_", s, ".pdf")), p, width = w, height = 5)
}
if (length(allrows)) write_csv(bind_rows(allrows) |> select(set, query, gene, gene_id, baseMean, comparison, log2FC = lfc, lfcSE = se, padj), file.path(OUT, "gene_bars_three_comparisons.csv"))
message("\nSaved to: ", OUT)
