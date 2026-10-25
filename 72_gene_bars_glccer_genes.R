#!/usr/bin/env Rscript
# Differential expression bar graphs, IFNB vs Mock and IFNG vs Mock, for glucosylceramide (GlcCer) genes, grouped by role (edit SPEC to add or move genes):
#   ceramide supply (Sptlc, Cers, Degs, Kdsr, Cert1) and ceramide from sphingomyelin / recycling (Smpd, Sgms, Asah, Acer)
#   GlcCer synthesis (Ugcg, with the UDP-glucose supplier Ugp2)
#   GlcCer degradation and lysosomal handling (Gba1 = lysosomal glucocerebrosidase, Gba2 = non-lysosomal, Psap/saposin C, Scarb2/LIMP-2 that traffics GBA1, Npc1/2)
#   GlcCer transport and flopping (Abcb1 / MDR1, Abca12, the glycolipid transfer protein Gltp, Plekha8 = FAPP2)
#   downstream glycosylation (LacCer synthases B4galt5/6, A4galt, gangliosides) and glycolipid catabolism (sialidases Neu1-4, Hexa/b, Glb1, Gla)
#   galactosylceramide / sulfatide genes (Ugt8, Gal3st1, Galc, Arsa) as related galactolipid metabolism
# One PNG with two panels side by side (IFNB vs Mock | IFNG vs Mock), genes in the same order in both, grouped by role and sorted by IFNB log2FC.
# Bars = shrunk log2FC (red up, blue down), lines = +/- 1 lfcSE, stars = padj (* < 0.05, ** < 0.01, *** < 0.001, ns). Mouse symbols are matched ignoring case; human names are
# mapped to the mouse gene (ALIAS). Genes removed by the pre-filter in 04 are listed as "Not found" and genes with baseMean < MIN_BASEMEAN are left off the figure (still in the CSV).
# Writes PNG + PDF to requant_v2_outputs/figures_glccer_genes/ plus glccer_genes_log2fc.csv and genes_not_found.csv
# Usage: Rscript 72_gene_bars_glccer_genes.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter

MIN_BASEMEAN <- 10
INK <- "#1F2933"; RED <- "#C53030"; BLUE <- "#2B6CB0"
COMP <- c(IFNB_vs_Mock = "IFNB vs Mock", IFNY_vs_Mock = "IFNG vs Mock")
ALIAS <- list(GBA1 = c("Gba1", "Gba"), CERT1 = c("Cert1", "Col4a3bp"), ABCB1 = "Abcb1*", PLEKHA8 = c("Plekha8", "Fapp2"), GBA3 = c("Gba3"))   # trailing * = every symbol starting with that
SPEC <- list(
  A = list(file = "GlcCer_genes", title = "Glucosylceramide (GlcCer): synthesis, metabolism, transport and modification", sets = list(
    `Ceramide supply (precursor)` = c("SPTLC1", "SPTLC2", "SPTLC3", "SPTSSA", "SPTSSB", "KDSR", "CERS1", "CERS2", "CERS3", "CERS4", "CERS5", "CERS6", "DEGS1", "DEGS2", "CERT1"),
    `Ceramide from sphingomyelin and recycling` = c("SMPD1", "SMPD2", "SMPD3", "SMPD4", "SGMS1", "SGMS2", "ASAH1", "ASAH2", "ACER1", "ACER2", "ACER3"),
    `GlcCer synthesis` = c("UGCG", "UGP2"),
    `GlcCer degradation and lysosomal handling` = c("GBA1", "GBA2", "GBA3", "PSAP", "SCARB2", "NPC1", "NPC2"),
    `GlcCer transport and flopping` = c("ABCB1", "ABCA12", "GLTP", "PLEKHA8"),
    `Downstream glycosylation (LacCer, gangliosides, globosides)` = c("B4GALT5", "B4GALT6", "A4GALT", "B3GALT4", "B3GNT5", "B4GALNT1", "ST3GAL5", "ST3GAL2", "ST8SIA1"),
    `Glycolipid catabolism (sialidases, hexosaminidases, galactosidases)` = c("NEU1", "NEU2", "NEU3", "NEU4", "HEXA", "HEXB", "GLB1", "GLA"),
    `Galactosylceramide and sulfatide (related)` = c("UGT8", "GAL3ST1", "GALC", "ARSA")))
)

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_glccer_genes"); dir.create(OUT, showWarnings = FALSE)
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
write_csv(D |> arrange(group, set, gene_name, comparison) |> select(group, set, gene_name, gene_id, baseMean, comparison, log2FC = lfc, lfcSE = se, padj), file.path(OUT, "glccer_genes_log2fc.csv"))
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
