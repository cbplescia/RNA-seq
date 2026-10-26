#!/usr/bin/env Rscript
# Same three-bar graph as 105 (left to right IFNB, EBOV, IFNB+EBOV vs Mock; shrunk log2FC +/- 1 lfcSE; stars = padj vs Mock; dagger = baseMean < 50) but for ALL cytochrome P450 genes (symbols Cyp<number>..., e.g. Cyp1a1, Cyp27a1, Cyp2e1)
# present in the DE tables (genes removed by the pre-filter in 04 are not shown). Genes are ordered by family / subfamily number and split into pages of PER_PAGE genes (env PER_PAGE, default 20).
# Optional: PATTERN=regex env to change the symbol match (default "^Cyp[0-9]"), ONLY_SIG=1 to keep only Cyp genes significant (padj < 0.05) in at least one comparison.
# Writes bars_cyp_pageN.png/.pdf and cyp_genes_three_comparisons.csv to requant_v2_outputs/figures_gene_bars/.   Usage: Rscript 110_cyp_genes_three_comparisons.R   (PROJ from 'source setup_session.sh')
suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter

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
PER_PAGE <- as.integer(Sys.getenv("PER_PAGE", "20")); PATTERN <- Sys.getenv("PATTERN", "^Cyp[0-9]"); ONLY_SIG <- Sys.getenv("ONLY_SIG", "0") == "1"
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

tabs <- lapply(names(COMP), function(nm) read_csv(file.path(RNA, paste0("DE_v2_", nm, "_annotated.csv")), show_col_types = FALSE)); names(tabs) <- names(COMP)
cyp <- tabs[[1]] |> filter(!is.na(gene_name), grepl(PATTERN, gene_name, ignore.case = TRUE)) |> distinct(gene_id, .keep_all = TRUE)
if (ONLY_SIG) { sg <- unique(unlist(lapply(tabs, function(t) t$gene_id[!is.na(t$padj) & t$padj < 0.05]))); cyp <- cyp |> filter(gene_id %in% sg) }
if (!nrow(cyp)) stop("No genes matching ", PATTERN, " in the DE tables")
# natural order: family number, subfamily letter(s), member number
num <- function(x, i) suppressWarnings(as.integer(regmatches(x, regexec("^[Cc]yp([0-9]+)([a-zA-Z]*)([0-9]*)", x))[[1]][i + 1]))
cyp <- cyp |> rowwise() |> mutate(f1 = num(gene_name, 1), f2 = num(gene_name, 2), f3 = { v <- regmatches(gene_name, regexec("^[Cc]yp([0-9]+)([a-zA-Z]*)([0-9]*)", gene_name))[[1]][3]; ifelse(is.na(v), "", v) }, f4 = { v <- num(gene_name, 3); ifelse(is.na(v), 0L, v) }) |> ungroup() |>
  arrange(f1, f3, f4, gene_name) |> mutate(pos = row_number())
message("\n== ", nrow(cyp), " Cyp genes in the DE tables: ", paste(cyp$gene_name, collapse = ", "))
D <- bind_rows(lapply(names(COMP), function(nm) tabs[[nm]] |> inner_join(cyp |> select(gene_id, pos), by = "gene_id") |> transmute(gene = gene_name, gene_id, pos, baseMean, comparison = COMP[[nm]], lfc = log2FoldChange, se = lfcSE, padj)))
D <- D |> mutate(comparison = factor(comparison, levels = COMP), lab = stars(padj), vj = ifelse(lfc >= 0, 0, 1), ys = lfc + ifelse(lfc >= 0, 1, -1) * (se + 0.03), xlab = paste0(gene, ifelse(baseMean < LOW_EXPR, "†", "")),
                 page = (pos - 1) %/% PER_PAGE + 1)
print(as.data.frame(D |> transmute(gene, baseMean = round(baseMean, 1), comparison, log2FC = round(lfc, 2), padj = signif(padj, 2)) |> tidyr::pivot_wider(names_from = comparison, values_from = c(log2FC, padj))), row.names = FALSE)
write_csv(D |> select(gene, gene_id, baseMean, comparison, log2FC = lfc, lfcSE = se, padj, page), file.path(OUT, "cyp_genes_three_comparisons.csv"))
npg <- max(D$page); dg <- position_dodge(width = 0.85)
for (pg in seq_len(npg)) {
  Dp <- D |> filter(page == pg); Dp$xlab <- factor(Dp$xlab, levels = unique(Dp$xlab[order(Dp$pos)]))
  p <- ggplot(Dp, aes(xlab, lfc, fill = comparison)) + geom_col(position = dg, width = 0.8) +
    geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), position = dg, width = 0.25, colour = INK, linewidth = 0.35) +
    geom_text(aes(y = ys, label = lab, vjust = vj), position = dg, colour = INK, size = 2.8) + geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) +
    scale_fill_manual(values = CCOL, name = NULL) + scale_y_continuous(expand = expansion(mult = 0.12)) +
    labs(title = paste0("Cytochrome P450 genes (Cyp): IFNB, EBOV and IFNB+EBOV vs Mock", if (npg > 1) paste0("  [page ", pg, " of ", npg, "]") else ""),
         subtitle = "left to right: IFNB, EBOV, IFNB+EBOV vs Mock; shrunk log2FC +/- 1 SE; stars = padj of each comparison (ns = not significant); dagger = baseMean < 50", x = NULL, y = "log2 fold-change vs Mock") +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3), panel.border = element_rect(colour = "grey70"),
          axis.text = element_text(colour = INK), axis.text.x = element_text(face = "italic", angle = 45, hjust = 1), plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 8.5), legend.position = "top")
  w <- max(6, 0.9 * nlevels(Dp$xlab) + 2)
  ggsave(file.path(OUT, paste0("bars_cyp_page", pg, ".png")), p, width = w, height = 5, dpi = 300); ggsave(file.path(OUT, paste0("bars_cyp_page", pg, ".pdf")), p, width = w, height = 5)
}
message("\nSaved ", npg, " page(s) to: ", OUT)
