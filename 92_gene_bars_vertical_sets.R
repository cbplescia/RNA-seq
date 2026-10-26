#!/usr/bin/env Rscript
# Vertical bar graphs: for each gene two bars side by side, IFNB vs Mock (blue) and EBOV vs Mock (green), in the order you list them (edit SETS):
#   glutathione / detoxification set   GSTP1, GSTM1, GSTA4, GSTT1, GSTM2, EPHX1, GSTA3, GSTP2, ALDH3A1, GGT7
#   secreted / receptor / lipoprotein  CLU, LRP2, AGRN, RAPSN, TCN2, APOD, CHRNA1, APOA5, C1S1
#   iron / heme / stress set           HMOX1, FTL1, FTH1, SLC40A1, TFR2, BLVRB, UGT1A6B, EPOR, FTL1-PS1, GCLM, LIF, HSD11B1, CTF1, TRF
# Shrunk log2FC from DE_v2_*_annotated.csv, error bars +/- 1 lfcSE, stars = padj of each comparison against Mock (* < 0.05, ** < 0.01, *** < 0.001, ns). A dagger on a gene name marks baseMean < LOW_EXPR
# (few reads, so the log2FC is less reliable); genes missing from the tables (removed by the pre-filter in 04, or a different symbol) are listed as "Not found". Mouse symbols are matched ignoring case;
# human-style names are mapped to the mouse gene in ALIAS (e.g. C1S1 -> C1s1).
# Writes bars_<set>.png / .pdf and gene_bars_vertical_sets.csv to requant_v2_outputs/figures_gene_bars/.   Usage: Rscript 92_gene_bars_vertical_sets.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter

SETS <- list(glutathione_detox = list(title = "Glutathione / detoxification genes", genes = c("GSTP1", "GSTM1", "GSTA4", "GSTT1", "GSTM2", "EPHX1", "GSTA3", "GSTP2", "ALDH3A1", "GGT7")),
             secreted_receptor = list(title = "CLU, LRP2, AGRN, RAPSN, TCN2, APOD, CHRNA1, APOA5, C1S1", genes = c("CLU", "LRP2", "AGRN", "RAPSN", "TCN2", "APOD", "CHRNA1", "APOA5", "C1S1")),
             iron_heme_stress = list(title = "Iron / heme / stress-response genes", genes = c("HMOX1", "FTL1", "FTH1", "SLC40A1", "TFR2", "BLVRB", "UGT1A6B", "EPOR", "FTL1-PS1", "GCLM", "LIF", "HSD11B1", "CTF1", "TRF")))
ALIAS <- list(C1S1 = c("C1s1", "C1s"), UGT1A6B = c("Ugt1a6b", "Ugt1a6"), `FTL1-PS1` = c("Ftl1-ps1"), TRF = c("Trf", "Tf"), FTL1 = c("Ftl1", "Ftl"), FTH1 = c("Fth1", "Fth"), GGT7 = c("Ggt7"))
LOW_EXPR <- 50
COMP <- c(IFNB_vs_Mock = "IFNB vs Mock", EBOV_vs_Mock = "EBOV vs Mock"); CCOL <- c(`IFNB vs Mock` = "#2B6CB0", `EBOV vs Mock` = "#1B7837")
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
  print(as.data.frame(D |> filter(comparison == COMP[["IFNB_vs_Mock"]]) |> transmute(gene = shown, baseMean = round(baseMean, 1), IFNB_log2FC = round(lfc, 2), IFNB_padj = signif(padj, 2)) |>
                        left_join(D |> filter(comparison == COMP[["EBOV_vs_Mock"]]) |> transmute(gene = shown, EBOV_log2FC = round(lfc, 2), EBOV_padj = signif(padj, 2)), by = "gene")), row.names = FALSE)
  dg <- position_dodge(width = 0.8)
  p <- ggplot(D, aes(xlab, lfc, fill = comparison)) + geom_col(position = dg, width = 0.78) +
    geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), position = dg, width = 0.25, colour = INK, linewidth = 0.35) +
    geom_text(aes(y = ys, label = lab, vjust = vj), position = dg, colour = INK, size = 2.8) + geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) +
    scale_fill_manual(values = CCOL, name = NULL) + scale_y_continuous(expand = expansion(mult = 0.12)) +
    labs(title = S$title, subtitle = "IFNB vs Mock and EBOV vs Mock side by side; shrunk log2FC +/- 1 SE; stars = padj of each comparison (ns = not significant); dagger = baseMean < 50", x = NULL, y = "log2 fold-change vs Mock") +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3), panel.border = element_rect(colour = "grey70"),
          axis.text = element_text(colour = INK), axis.text.x = element_text(face = "italic", angle = 45, hjust = 1), plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 8.5), legend.position = "top")
  w <- max(5.5, 0.62 * nlevels(D$xlab) + 2)
  ggsave(file.path(OUT, paste0("bars_", s, ".png")), p, width = w, height = 5, dpi = 300); ggsave(file.path(OUT, paste0("bars_", s, ".pdf")), p, width = w, height = 5)
}
if (length(allrows)) write_csv(bind_rows(allrows) |> select(set, query, gene, gene_id, baseMean, comparison, log2FC = lfc, lfcSE = se, padj), file.path(OUT, "gene_bars_vertical_sets.csv"))
message("\nSaved to: ", OUT)
