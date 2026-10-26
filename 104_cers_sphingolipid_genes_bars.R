#!/usr/bin/env Rscript
# Vertical bar graphs: for each gene THREE bars side by side, left to right IFNB vs Mock (blue), EBOV vs Mock (dark green), IFNB+EBOV vs Mock (light green).
#  (1) all ceramide synthase genes Cers1-Cers6 (old names Lass1-6 are matched too)                                  -> bars_Cers_genes
#  (2) all sphingolipid-affiliated genes, in panels (my curated list; edit SPH):                                      -> bars_sphingolipid_genes_p<n> (paged)
#      de novo synthesis (Sptlc1/2/3, Sptssa/b, Kdsr), ceramide synthases and regulators (Cers1-6, Ormdl1-3, Degs1/2, Fa2h, Elovl1/4), sphingomyelin synthesis and breakdown (Sgms1/2, Samd8, Smpd1-4, Smpdl3a/b, Enpp7),
#      ceramidases (Asah1/2, Acer1-3), ceramide transport and phosphorylation (Cert1, Cptp, Cerk, Cerkl), sphingosine and S1P metabolism and transport (Sphk1/2, Sgpl1, Sgpp1/2, Plpp1/3, Spns2, Mfsd2b, Abcc1),
#      S1P receptors (S1pr1-5), glycosphingolipid synthesis (Ugcg, B4galt5/6, A4galt, B3galt4, B4galnt1, St3gal5, St3gal2, St8sia1, Ugt8, Gal3st1) and lysosomal glycosphingolipid breakdown
#      (Gba1, Gba2, Galc, Arsa, Gla, Glb1, Hexa, Hexb, Gm2a, Psap, Neu1-4)
# Shrunk log2FC from DE_v2_*_annotated.csv, +/- 1 lfcSE, stars = padj of each comparison against Mock (* < 0.05, ** < 0.01, *** < 0.001, ns), dagger = baseMean < LOW_EXPR (log2FC less reliable).
# Genes removed by the pre-filter in 04 (no counts) cannot be plotted and are listed as "Not found". Writes the figures (PNG + PDF) and cers_genes_log2fc.csv / sphingolipid_genes_log2fc.csv to requant_v2_outputs/figures_gene_bars/.
# Usage: Rscript 104_cers_sphingolipid_genes_bars.R   (PROJ from 'source setup_session.sh';  PER_PAGE=30 optional)

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter

CERS <- list(`ceramide synthases` = c("Cers1", "Cers2", "Cers3", "Cers4", "Cers5", "Cers6"))
SPH <- list(`de novo synthesis` = c("Sptlc1", "Sptlc2", "Sptlc3", "Sptssa", "Sptssb", "Kdsr"),
            `ceramide synthases + regulators` = c("Cers1", "Cers2", "Cers3", "Cers4", "Cers5", "Cers6", "Ormdl1", "Ormdl2", "Ormdl3", "Degs1", "Degs2", "Fa2h", "Elovl1", "Elovl4"),
            `sphingomyelin synthesis / breakdown` = c("Sgms1", "Sgms2", "Samd8", "Smpd1", "Smpd2", "Smpd3", "Smpd4", "Smpdl3a", "Smpdl3b", "Enpp7"),
            `ceramidases` = c("Asah1", "Asah2", "Acer1", "Acer2", "Acer3"),
            `ceramide transport / phosphorylation` = c("Cert1", "Cptp", "Cerk", "Cerkl"),
            `sphingosine / S1P metabolism + transport` = c("Sphk1", "Sphk2", "Sgpl1", "Sgpp1", "Sgpp2", "Plpp1", "Plpp3", "Spns2", "Mfsd2b", "Abcc1"),
            `S1P receptors` = c("S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5"),
            `glycosphingolipid synthesis` = c("Ugcg", "B4galt5", "B4galt6", "A4galt", "B3galt4", "B4galnt1", "St3gal5", "St3gal2", "St8sia1", "Ugt8", "Gal3st1"),
            `lysosomal glycosphingolipid breakdown` = c("Gba1", "Gba2", "Galc", "Arsa", "Gla", "Glb1", "Hexa", "Hexb", "Gm2a", "Psap", "Neu1", "Neu2", "Neu3", "Neu4"))
ALIAS <- list(Cers1 = "Lass1", Cers2 = "Lass2", Cers3 = "Lass3", Cers4 = "Lass4", Cers5 = "Lass5", Cers6 = "Lass6", Gba1 = "Gba", Cert1 = "Col4a3bp", Cptp = "Gltpd1", Samd8 = "Sgms-r", Plpp1 = "Ppap2a", Plpp3 = "Ppap2b", Spns2 = "Spns2")
PER_PAGE <- as.integer(Sys.getenv("PER_PAGE", unset = "40")); LOW_EXPR <- 50
COMP <- c(IFNB_vs_Mock = "IFNB vs Mock", EBOV_vs_Mock = "EBOV vs Mock", IFNB_EBOV_vs_Mock = "IFNB+EBOV vs Mock"); CCOL <- c(`IFNB vs Mock` = "#2B6CB0", `EBOV vs Mock` = "#1B7837", `IFNB+EBOV vs Mock` = "#9BD08E"); INK <- "#1F2933"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_gene_bars"); dir.create(OUT, showWarnings = FALSE)
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
tabs <- lapply(names(COMP), function(nm) read_csv(file.path(RNA, paste0("DE_v2_", nm, "_annotated.csv")), show_col_types = FALSE) |> mutate(key = tolower(gene_name))); names(tabs) <- names(COMP)
find_gene <- function(q) { syms <- tolower(c(q, ALIAS[[q]])); h <- tabs[[1]] |> filter(key %in% syms) |> distinct(gene_id, .keep_all = TRUE); if (nrow(h)) h$gene_id[1] else NA_character_ }
gene_table <- function(L) { GT <- bind_rows(lapply(names(L), function(p) tibble(panel = p, query = L[[p]]))) |> distinct(query, .keep_all = TRUE); GT$gene_id <- vapply(GT$query, find_gene, "")
  miss <- GT$query[is.na(GT$gene_id)]; if (length(miss)) message("Not found in the DE tables (removed by the pre-filter, or a different symbol): ", paste(miss, collapse = ", "))
  list(GT = GT |> filter(!is.na(gene_id)) |> mutate(pos = row_number()), miss = miss) }
collect <- function(GT) bind_rows(lapply(names(COMP), function(nm) tabs[[nm]] |> inner_join(GT, by = "gene_id") |>
  transmute(panel, gene = gene_name, gene_id, pos, baseMean, comparison = COMP[[nm]], lfc = log2FoldChange, se = lfcSE, padj))) |>
  mutate(comparison = factor(comparison, levels = COMP), lab = stars(padj), vj = ifelse(lfc >= 0, 0, 1), ys = lfc + ifelse(lfc >= 0, 1, -1) * (se + 0.03), xl = paste0(gene, ifelse(baseMean < LOW_EXPR, "†", "")))
draw <- function(D, title, fname, extra = "") {
  D$panel <- factor(D$panel, levels = unique(D$panel[order(D$pos)])); D$key <- paste(D$panel, D$xl, sep = "||"); D$key <- factor(D$key, levels = unique(D$key[order(D$pos)]))
  dg <- position_dodge(width = 0.85)
  p <- ggplot(D, aes(key, lfc, fill = comparison)) + geom_col(position = dg, width = 0.8) + geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), position = dg, width = 0.25, colour = INK, linewidth = 0.3) +
    geom_text(aes(y = ys, label = lab, vjust = vj), position = dg, colour = INK, size = 2.2) + geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) +
    facet_grid(. ~ panel, scales = "free_x", space = "free_x", labeller = labeller(panel = label_wrap_gen(15))) + scale_fill_manual(values = CCOL, name = NULL) +
    scale_x_discrete(labels = function(x) sub("^.*\\|\\|", "", x)) + scale_y_continuous(expand = expansion(mult = 0.12)) +
    labs(title = title, subtitle = sprintf("left to right: IFNB, EBOV, IFNB+EBOV vs Mock; shrunk log2FC +/- 1 SE; stars = padj of each comparison (ns = not significant); dagger = baseMean < %d%s", LOW_EXPR, extra), x = NULL, y = "log2 fold-change vs Mock") +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3), panel.border = element_rect(colour = "grey70"),
          axis.text = element_text(colour = INK), axis.text.x = element_text(face = "italic", angle = 50, hjust = 1, size = 7.5), plot.title = element_text(colour = INK, face = "bold", size = 11), plot.subtitle = element_text(colour = "grey35", size = 7.5),
          strip.background = element_rect(fill = "grey95", colour = "grey70"), strip.text = element_text(colour = INK, face = "bold", size = 7), legend.position = "top", panel.spacing.x = unit(0.15, "lines"))
  w <- max(6.5, 0.62 * nlevels(D$key) + 1.2 * nlevels(D$panel) + 2)
  ggsave(file.path(OUT, paste0(fname, ".png")), p, width = w, height = 5.5, dpi = 300, limitsize = FALSE); ggsave(file.path(OUT, paste0(fname, ".pdf")), p, width = w, height = 5.5, limitsize = FALSE)
}
tsv <- function(D) D |> arrange(pos, comparison) |> select(panel, gene, gene_id, baseMean, comparison, log2FC = lfc, lfcSE = se, padj)

# ---- (1) Cers ----
g1 <- gene_table(CERS)
if (nrow(g1$GT)) { D1 <- collect(g1$GT); write_csv(tsv(D1), file.path(OUT, "cers_genes_log2fc.csv"))
  print(as.data.frame(D1 |> transmute(gene, baseMean = round(baseMean, 1), comparison, log2FC = round(lfc, 2), padj = signif(padj, 2)) |> pivot_wider(names_from = comparison, values_from = c(log2FC, padj))), row.names = FALSE)
  draw(D1, "Ceramide synthase genes (Cers1-6)", "bars_Cers_genes", sprintf("; %d of 6 found", nrow(g1$GT))) } else message("No Cers gene found in the DE tables")

# ---- (2) all sphingolipid genes ----
g2 <- gene_table(SPH); GT <- g2$GT; if (!nrow(GT)) stop("No sphingolipid gene found")
D2 <- collect(GT); write_csv(tsv(D2), file.path(OUT, "sphingolipid_genes_log2fc.csv"))
message(sprintf("sphingolipid genes: %d found of %d listed", nrow(GT), nrow(GT) + length(g2$miss)))
pg <- GT |> mutate(page = ceiling(row_number() / PER_PAGE)); np <- max(pg$page)
for (k in seq_len(np)) {
  ids <- pg$gene_id[pg$page == k]; d <- D2 |> filter(gene_id %in% ids)
  draw(d, "Sphingolipid-affiliated genes", if (np > 1) sprintf("bars_sphingolipid_genes_p%d", k) else "bars_sphingolipid_genes",
       sprintf("; %d genes found, %d not in the data%s", nrow(GT), length(g2$miss), if (np > 1) sprintf("; page %d of %d", k, np) else ""))
}
message("Saved to: ", OUT)
