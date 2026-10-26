#!/usr/bin/env Rscript
# Differential expression bar graphs, IFNB vs Mock and EBOV vs Mock, for host genes that synthesise, alter, modify or transport lipids in the ENDOCYTIC system
# (my classification, merged from the endocytic / lipid gene lists used earlier - check it and edit SPEC):
#   A  endolysosomal, in or facing the LUMEN: lysosomal hydrolases (cathepsins, Gba1, Hexa/b, Lipa, Smpd1, Asah1 ...), membrane proteins (Lamp1/2, Scarb2, Npc1, tetraspanins, Ifitm1-3, Cybb ...),
#     TRPML / TPC channels and V-ATPase subunits
#   B  lipid TRANSLOCASES (flippases + CDC50, ABC floppases, scramblases) - change leaflet composition incl. the lumen-facing leaflet, whichever side the protein faces
#   C  endocytic lipid machinery on the CYTOPLASMIC face: phosphoinositide kinases / phosphatases, lipid transfer and contact-site proteins, LC3 lipidation (Atg3/7), cytoplasmic-leaflet enzymes
# One PNG per group with two panels side by side (IFNB vs Mock | EBOV vs Mock), genes in the same order in both, grouped by function and sorted by IFNB log2FC.
# Bars = shrunk log2FC (red up, blue down), lines = +/- 1 lfcSE, stars = padj (* < 0.05, ** < 0.01, *** < 0.001, ns). Mouse symbols are matched ignoring case; human names are mapped to the mouse gene
# (ALIAS). Genes removed by the pre-filter in 04 are listed as "Not found"; genes with baseMean < MIN_BASEMEAN are left off the figures (still in the CSV).
# Writes PNG + PDF to requant_v2_outputs/figures_endocytic_lipid_ifnb_ebov/ plus endocytic_lipid_genes_log2fc.csv and genes_not_found.csv
# Usage: Rscript 86_gene_bars_endocytic_lipid_ifnb_ebov.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter

MIN_BASEMEAN <- 10
INK <- "#1F2933"; RED <- "#C53030"; BLUE <- "#2B6CB0"
COMP <- c(IFNB_vs_Mock = "IFNB vs Mock", EBOV_vs_Mock = "EBOV vs Mock")
ALIAS <- list(GBA1 = c("Gba1", "Gba"), CERT1 = c("Cert1", "Col4a3bp"), ABCB1 = "Abcb1*", COQ8A = c("Coq8a", "Adck3"), COQ8B = c("Coq8b", "Adck4"), PRELID3B = c("Prelid3b", "Slmo2"),
              PEDS1 = c("Peds1", "Tmem189"), TAFAZZIN = c("Tafazzin", "Taz"), PLPP1 = c("Plpp1", "Ppap2a"), PLPP2 = c("Plpp2", "Ppap2c"), PLPP3 = c("Plpp3", "Ppap2b"),
              SELENOI = c("Selenoi", "Ept1"), TMEM24 = c("Tmem24", "C2cd2l"))   # trailing * = every symbol starting with that
SPEC <- list(
  A = list(file = "A_endolysosomal_lumen_facing", title = "Endolysosomal: lumen / lumen-facing lipid genes", sets = list(
    `Lysosomal hydrolases and lumenal enzymes` = c("ARSA", "ASAH1", "CTSA", "CTSB", "CTSD", "CTSL", "GALC", "GBA1", "GLB1", "GM2A", "GRN", "HEXA", "HEXB", "LIPA", "NEU1", "PLA2G15", "PLD3", "PLD4", "PSAP", "SMPD1", "NPC2"),
    `Endosomal / lysosomal membrane and lumen-facing proteins` = c("CLN3", "CLN5", "LAMP1", "LAMP2", "SCARB2", "NPC1", "TMEM106B", "SPNS1", "PLPP1", "PLPP2", "PLPP3", "CD63", "CD81", "CD9", "TSPAN6", "IFITM1", "IFITM2", "IFITM3", "CYBB", "MCOLN1", "TPCN1", "TPCN2", "SLC38A9", "TLR4"),
    `V-ATPase (lumen acidification)` = c("ATP6AP1", "ATP6V0A1", "ATP6V0D1", "ATP6V1A", "ATP6V1B2"))),
  B = list(file = "B_lipid_translocases", title = "Lipid translocases (flippases, floppases, scramblases): change leaflet asymmetry incl. the lumen-facing leaflet", sets = list(
    `Flippases (P4-ATPases) and CDC50 subunits` = c("ATP8A1", "ATP8A2", "ATP8B1", "ATP8B2", "ATP8B3", "ATP8B4", "ATP9A", "ATP9B", "ATP10A", "ATP10B", "ATP10D", "ATP11A", "ATP11B", "ATP11C", "TMEM30A", "TMEM30B", "TMEM30C"),
    `Floppases (ABC transporters)` = c("ABCA1", "ABCA3", "ABCA4", "ABCA7", "ABCA12", "ABCB1", "ABCB4", "ABCC1", "ABCG1", "ABCG2", "ABCG4", "ABCG5", "ABCG8"),
    `Scramblases and lipid-transporting membrane proteins` = c("ANO3", "ANO4", "ANO5", "ANO6", "ANO7", "ANO9", "ANO10", "PLSCR1", "PLSCR2", "PLSCR3", "PLSCR4", "XKR4", "XKR8", "XKR9", "TMEM63A", "TMEM63B", "ATG9A", "ATG9B", "TMEM41B", "VMP1", "SPNS2", "MFSD2A", "MFSD2B"))),
  C = list(file = "C_endocytic_cytoplasmic_side", title = "Endocytic lipid machinery acting on the cytoplasmic face", sets = list(
    `Phosphoinositide kinases and phosphatases` = c("FIG4", "INPP4A", "INPP4B", "INPP5B", "INPP5D", "INPP5E", "INPP5F", "INPPL1", "MTM1", "MTMR1", "MTMR2", "MTMR3", "MTMR4", "MTMR6", "OCRL", "PI4K2A", "PI4K2B", "PI4KA", "PI4KB", "PIK3C2A", "PIK3C2B", "PIK3C2G", "PIK3C3", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIKFYVE", "PIP4K2A", "PIP4K2B", "PIP4K2C", "PIP5K1A", "PIP5K1B", "PIP5K1C", "PTEN", "SACM1L", "SYNJ1", "SYNJ2", "VAC14"),
    `Lipid transfer and contact-site proteins` = c("ATG2A", "ATG2B", "CERT1", "CPTP", "ESYT1", "ESYT2", "ESYT3", "GLTP", "GRAMD1A", "GRAMD1B", "GRAMD1C", "OSBP", "OSBPL1A", "OSBPL2", "OSBPL3", "OSBPL5", "OSBPL8", "OSBPL9", "OSBPL10", "OSBPL11", "PITPNA", "PITPNB", "PITPNM1", "PITPNM2", "PITPNM3", "STARD2", "STARD3", "STARD3NL", "STARD4", "STARD5", "STARD10", "TMEM24", "VAPA", "VAPB", "VPS13A", "VPS13B", "VPS13C", "VPS13D"),
    `LC3 lipidation (autophagy conjugation)` = c("ATG3", "ATG7"),
    `Cytoplasmic-leaflet lipid enzymes and signalling` = c("PLD1", "PLD2", "PLCB1", "PLCB2", "PLCB3", "PLCD1", "PLCE1", "PLCG1", "PLCG2", "DAGLA", "DAGLB", "MGLL", "ABHD6", "DGKA", "DGKB", "DGKD", "DGKE", "DGKG", "DGKH", "DGKI", "DGKK", "DGKQ", "DGKZ", "SPHK1", "SPHK2", "SMPD3", "PLA2G4A", "PLA2G6", "DDHD1", "DDHD2", "LRRK2")))
)

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_endocytic_lipid_ifnb_ebov"); dir.create(OUT, showWarnings = FALSE)
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
write_csv(D |> arrange(group, set, gene_name, comparison) |> select(group, set, gene_name, gene_id, baseMean, comparison, log2FC = lfc, lfcSE = se, padj), file.path(OUT, "endocytic_lipid_genes_log2fc.csv"))
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
