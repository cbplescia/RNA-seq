#!/usr/bin/env Rscript
# Same three-bar graph as 105 (left to right IFNB, EBOV, IFNB+EBOV vs Mock; shrunk log2FC +/- 1 lfcSE; stars = padj vs Mock; dagger = baseMean < 50) but for ALL cytochrome P450 genes (symbols Cyp<number>..., e.g. Cyp1a1, Cyp27a1, Cyp2e1)
# present in the DE tables (genes removed by the pre-filter in 04 are not shown). Genes are grouped by what they mainly metabolize (GROUPS below; edit to taste), ordered by family number inside each group, shown as panels per group, and split into pages of up to PER_PAGE genes (env PER_PAGE, default 25) without splitting a group across pages unless it is larger than a page.
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
PER_PAGE <- as.integer(Sys.getenv("PER_PAGE", "25")); PATTERN <- Sys.getenv("PATTERN", "^Cyp[0-9]"); ONLY_SIG <- Sys.getenv("ONLY_SIG", "0") == "1"
# substrate class by gene-symbol pattern (first match wins; mainly after the Nebert/Guengerich P450 classification, mouse paralogues share the class of their family). Many P450s act on several substrates; the group is the dominant/best-known one.
GROUPS <- list(
  "Steroid hormone synthesis"             = "^Cyp(11a1|11b[12]|17a1|19a1|21a1)",
  "Cholesterol, oxysterols & bile acids"  = "^Cyp(7a1|7b1|8b1|27a1|39a1|46a1|51(a1)?)$",
  "Vitamin D & retinoic acid"             = "^Cyp(2r1|24a1|27b1|26[abc]1)$",
  "Fatty acids & eicosanoids"             = "^Cyp(4a|4f|4v|4x|2j|2u1|5a1|8a1)",
  "Drugs & xenobiotics (liver-type)"      = "^Cyp(1a|1b|2a|2b|2c|2d|2e|2f|2g|2s|2t|3a|4b)",
  "Other / orphan"                        = ".")
grp_of <- function(x) { for (g in names(GROUPS)) if (grepl(GROUPS[[g]], x, ignore.case = TRUE)) return(g); "Other / orphan" }
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

tabs <- lapply(names(COMP), function(nm) read_csv(file.path(RNA, paste0("DE_v2_", nm, "_annotated.csv")), show_col_types = FALSE)); names(tabs) <- names(COMP)
cyp <- tabs[[1]] |> filter(!is.na(gene_name), grepl(PATTERN, gene_name, ignore.case = TRUE)) |> distinct(gene_id, .keep_all = TRUE)
if (ONLY_SIG) { sg <- unique(unlist(lapply(tabs, function(t) t$gene_id[!is.na(t$padj) & t$padj < 0.05]))); cyp <- cyp |> filter(gene_id %in% sg) }
if (!nrow(cyp)) stop("No genes matching ", PATTERN, " in the DE tables")
# natural order: family number, subfamily letter(s), member number
num <- function(x, i) suppressWarnings(as.integer(regmatches(x, regexec("^[Cc]yp([0-9]+)([a-zA-Z]*)([0-9]*)", x))[[1]][i + 1]))
cyp <- cyp |> rowwise() |> mutate(f1 = num(gene_name, 1), f2 = num(gene_name, 2), f3 = { v <- regmatches(gene_name, regexec("^[Cc]yp([0-9]+)([a-zA-Z]*)([0-9]*)", gene_name))[[1]][3]; ifelse(is.na(v), "", v) }, f4 = { v <- num(gene_name, 3); ifelse(is.na(v), 0L, v) }) |> ungroup() |>
  mutate(group = factor(vapply(gene_name, grp_of, ""), levels = names(GROUPS))) |> arrange(group, f1, f3, f4, gene_name) |> mutate(pos = row_number())
message("\n== ", nrow(cyp), " Cyp genes in the DE tables, by group:"); for (g in levels(cyp$group)) if (any(cyp$group == g)) message("  ", g, " (", sum(cyp$group == g), "): ", paste(cyp$gene_name[cyp$group == g], collapse = ", "))
D <- bind_rows(lapply(names(COMP), function(nm) tabs[[nm]] |> inner_join(cyp |> select(gene_id, pos, group), by = "gene_id") |> transmute(group, gene = gene_name, gene_id, pos, baseMean, comparison = COMP[[nm]], lfc = log2FoldChange, se = lfcSE, padj)))
D <- D |> mutate(comparison = factor(comparison, levels = COMP), lab = stars(padj), vj = ifelse(lfc >= 0, 0, 1), ys = lfc + ifelse(lfc >= 0, 1, -1) * (se + 0.03), xlab = paste0(gene, ifelse(baseMean < LOW_EXPR, "†", "")))
# pages: fill whole groups until PER_PAGE is exceeded
gsz <- D |> distinct(group, gene) |> dplyr::count(group) |> arrange(group); pg <- 1; used <- 0; gpage <- c()
for (i in seq_len(nrow(gsz))) { if (used > 0 && used + gsz$n[i] > PER_PAGE) { pg <- pg + 1; used <- 0 }; gpage[as.character(gsz$group[i])] <- pg; used <- used + gsz$n[i] }
D$page <- unname(gpage[as.character(D$group)])
print(as.data.frame(D |> transmute(group, gene, baseMean = round(baseMean, 1), comparison, log2FC = round(lfc, 2), padj = signif(padj, 2)) |> tidyr::pivot_wider(names_from = comparison, values_from = c(log2FC, padj))), row.names = FALSE)
write_csv(D |> select(group, gene, gene_id, baseMean, comparison, log2FC = lfc, lfcSE = se, padj, page), file.path(OUT, "cyp_genes_three_comparisons.csv"))
npg <- max(D$page); dg <- position_dodge(width = 0.85)
for (pg in seq_len(npg)) {
  Dp <- D |> filter(page == pg); Dp$xlab <- factor(Dp$xlab, levels = unique(Dp$xlab[order(Dp$pos)]))
  p <- ggplot(Dp, aes(xlab, lfc, fill = comparison)) + geom_col(position = dg, width = 0.8) +
    geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), position = dg, width = 0.25, colour = INK, linewidth = 0.35) +
    geom_text(aes(y = ys, label = lab, vjust = vj), position = dg, colour = INK, size = 2.8) + geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) +
    facet_grid(. ~ group, scales = "free_x", space = "free_x", labeller = labeller(group = label_wrap_gen(22))) +
    scale_fill_manual(values = CCOL, name = NULL) + scale_y_continuous(expand = expansion(mult = 0.12)) +
    labs(title = paste0("Cytochrome P450 genes (Cyp): IFNB, EBOV and IFNB+EBOV vs Mock", if (npg > 1) paste0("  [page ", pg, " of ", npg, "]") else ""),
         subtitle = "left to right: IFNB, EBOV, IFNB+EBOV vs Mock; shrunk log2FC +/- 1 SE; stars = padj of each comparison (ns = not significant); dagger = baseMean < 50", x = NULL, y = "log2 fold-change vs Mock") +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3), panel.border = element_rect(colour = "grey70"),
          axis.text = element_text(colour = INK), axis.text.x = element_text(face = "italic", angle = 45, hjust = 1), plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 8.5), legend.position = "top", strip.background = element_rect(fill = "grey93", colour = "grey70"), strip.text = element_text(colour = INK, face = "bold", size = 8.5))
  w <- max(6, 0.9 * nlevels(Dp$xlab) + 2)
  ggsave(file.path(OUT, paste0("bars_cyp_page", pg, ".png")), p, width = w, height = 5, dpi = 300); ggsave(file.path(OUT, paste0("bars_cyp_page", pg, ".pdf")), p, width = w, height = 5)
}
message("\nSaved ", npg, " page(s) to: ", OUT)
