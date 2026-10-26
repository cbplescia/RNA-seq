#!/usr/bin/env Rscript
# Differential expression bar graphs, IFNB vs Mock and IFNG vs Mock, for selected sphingolipid genes (edit SPEC to add or move genes):
#   sphingosine kinases Sphk1, Sphk2 (make S1P) | S1P degradation Sgpl1 and phosphatases Sgpp1, Sgpp2 | ceramide synthesis Cers1, Cers3, Degs1 and the SPT / ceramide-synthesis
#   inhibitors Ormdl1-3 | S1P exporter Spns2 | ceramide / GlcCer hydrolysis Asah1 (acid ceramidase), Gba2 (non-lysosomal glucocerebrosidase)
# One PNG with two panels side by side (IFNB vs Mock | IFNG vs Mock), genes in the same order in both, grouped by role and sorted by IFNB log2FC.
# Bars = shrunk log2FC (red up, blue down), lines = +/- 1 lfcSE, stars = padj (* < 0.05, ** < 0.01, *** < 0.001, ns). Genes removed by the pre-filter in 04 are listed as "Not found"
# and genes with baseMean < MIN_BASEMEAN are left off the figure (still in the CSV).
# Writes PNG + PDF to requant_v2_outputs/figures_sphingolipid_genes/ plus sphingolipid_genes_log2fc.csv and genes_not_found.csv
# Usage: Rscript 73_gene_bars_sphingolipid_genes.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter

MIN_BASEMEAN <- 10
INK <- "#1F2933"; RED <- "#C53030"; BLUE <- "#2B6CB0"
COMP <- c(IFNB_vs_Mock = "IFNB vs Mock", IFNY_vs_Mock = "IFNG vs Mock")
ALIAS <- list(GBA2 = c("Gba2"))   # trailing * = every symbol starting with that; add human -> mouse names here if a symbol is not found
SPEC <- list(
  A = list(file = "sphingolipid_genes", title = "Sphingolipid genes: S1P, ceramide and ceramide / GlcCer hydrolysis", sets = list(
    `Sphingosine kinases (S1P synthesis)` = c("SPHK1", "SPHK2"),
    `S1P degradation and phosphatases` = c("SGPL1", "SGPP1", "SGPP2"),
    `Ceramide synthesis and its regulation` = c("CERS1", "CERS3", "ORMDL1", "ORMDL2", "ORMDL3", "DEGS1"),
    `S1P export` = c("SPNS2"),
    `Ceramide / glucosylceramide hydrolysis` = c("ASAH1", "GBA2")))
)

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_sphingolipid_genes"); dir.create(OUT, showWarnings = FALSE)
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

tabs <- lapply(names(COMP), function(nm) read_csv(file.path(RNA, paste0("DE_v2_", nm, "_annotated.csv")), show_col_types = FALSE)); names(tabs) <- names(COMP)
allg <- tabs[[1]] |> mutate(key = tolower(gene_name))

# ---- map every queried symbol to gene ids ----
Q <- bind_rows(lapply(names(SPEC), function(k) bind_rows(lapply(names(SPEC[[k]]$sets), function(s) tibble(group = k, set = s, query = SPEC[[k]]$sets[[s]])))))
find <- function(q) {
  nm <- if (!is.null(ALIAS[[q]])) ALIAS[[q]] else q
  ex <- tolower(nm[!grepl("\\*$", nm)]); wl <- tolower(sub("\\*$", "", nm[grepl("\\*$", nm)]))
  hit <- allg$key %in% ex | (if (length(wl)) Reduce(`|`, lapply(wl, function(p) startsWith(allg$key, p))) else FALSE)
  allg |> filter(hit) |> distinct(gene_id, .keep_all = TRUE) |> select(gene_id, gene_name)
}
H <- bind_rows(lapply(seq_len(nrow(Q)), function(i) { h <- find(Q$query[i]); if (nrow(h)) bind_cols(Q[rep(i, nrow(h)), ], h) else NULL }))
nf <- Q |> filter(!query %in% H$query)
write_csv(nf, file.path(OUT, "genes_not_found.csv"))
message(sprintf("%d of %d queried genes found; not found (removed by the pre-filter, not expressed, or a different symbol): %s", n_distinct(H$query), nrow(Q), paste(nf$query, collapse = ", ")))

# ---- results for both contrasts ----
D <- bind_rows(lapply(names(COMP), function(nm) tabs[[nm]] |> filter(gene_id %in% H$gene_id) |>
  transmute(gene_id, comparison = COMP[[nm]], baseMean, lfc = log2FoldChange, se = lfcSE, padj)))
D <- H |> inner_join(D, by = "gene_id", relationship = "many-to-many") |> mutate(comparison = factor(comparison, levels = COMP), lab = stars(padj))
write_csv(D |> arrange(group, set, gene_name, comparison) |> select(group, set, gene_name, gene_id, baseMean, comparison, log2FC = lfc, lfcSE = se, padj), file.path(OUT, "sphingolipid_genes_log2fc.csv"))
low <- D |> group_by(gene_id, gene_name) |> summarise(bm = max(baseMean, na.rm = TRUE), .groups = "drop") |> filter(bm < MIN_BASEMEAN)
if (nrow(low)) message("Left off the figures (baseMean < ", MIN_BASEMEAN, "): ", paste(low$gene_name, collapse = ", "))
D <- D |> filter(!gene_id %in% low$gene_id)

# ---- one PNG per group ----
for (k in names(SPEC)) {
  d <- D |> filter(group == k)
  if (!nrow(d)) { message("Group ", k, ": no genes to plot"); next }
  ord <- d |> filter(comparison == COMP[["IFNB_vs_Mock"]]) |> group_by(set, gene_id) |> summarise(o = max(lfc, na.rm = TRUE), .groups = "drop")
  d <- d |> left_join(ord, by = c("set", "gene_id")) |> mutate(set = factor(set, levels = names(SPEC[[k]]$sets)), key = paste0(set, "||", gene_name), dir = ifelse(lfc >= 0, "up", "down"))
  lv <- d |> distinct(set, key, o) |> arrange(set, o) |> pull(key)                                 # lowest at the bottom, highest at the top of each block
  d <- d |> mutate(key = factor(key, levels = lv), hj = ifelse(lfc >= 0, 0, 1), xs = lfc + ifelse(lfc >= 0, 1, -1) * (se + 0.05))
  n <- n_distinct(d$key)
  p <- ggplot(d, aes(lfc, key, fill = dir)) + geom_col(width = 0.75) +
    geom_segment(aes(x = lfc - se, xend = lfc + se, yend = key), colour = INK, linewidth = 0.3) +
    geom_text(aes(x = xs, label = lab, hjust = hj), colour = INK, size = 2.3, vjust = 0.45) +
    geom_vline(xintercept = 0, colour = INK, linewidth = 0.4) +
    facet_grid(set ~ comparison, scales = "free_y", space = "free_y", labeller = labeller(set = label_wrap_gen(24))) +
    scale_fill_manual(values = c(up = RED, down = BLUE), guide = "none") +
    scale_y_discrete(labels = function(x) sub("^.*\\|\\|", "", x)) + scale_x_continuous(expand = expansion(mult = 0.1)) +
    labs(title = SPEC[[k]]$title, subtitle = "shrunk log2FC +/- 1 SE vs Mock; stars = padj (ns = not significant); genes grouped by function, sorted by IFNB log2FC",
         x = "log2 fold-change vs Mock", y = NULL) +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
      panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK), axis.text.y = element_text(face = "italic", size = if (n > 80) 6.5 else 8),
      plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 8.5), strip.background = element_rect(fill = "grey95", colour = "grey70"),
      strip.text = element_text(colour = INK, face = "bold", size = 8), strip.text.y = element_text(angle = 0, hjust = 0), panel.spacing.y = unit(0.15, "lines"))
  h <- min(48, max(4, 0.15 * n + 0.45 * n_distinct(d$set) + 1.8))
  ggsave(file.path(OUT, paste0("bars_", SPEC[[k]]$file, ".png")), p, width = 11, height = h, dpi = 300, limitsize = FALSE)
  ggsave(file.path(OUT, paste0("bars_", SPEC[[k]]$file, ".pdf")), p, width = 11, height = h, limitsize = FALSE)
  message(sprintf("Group %s: %d genes plotted", k, n))
}
message("Saved to: ", OUT)
