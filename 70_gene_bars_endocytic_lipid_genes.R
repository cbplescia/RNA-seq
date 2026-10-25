#!/usr/bin/env Rscript
# Differential expression bar graphs, IFNB vs Mock and IFNG vs Mock, for 242 host genes involved in endocytic / lysosomal membranes and lipid handling,
# split into three PNGs by where they act (my classification - check it and move genes between groups in SPEC below):
#   A  endolysosomal, in or facing the LUMEN: lumenal hydrolases, membrane proteins with luminal domains (LAMPs, tetraspanins, NPC1, IFITMs, CYBB, TLR4 ...),
#      V-ATPase subunits, and lipid floppases / exporters (ABCA1, ABCG1) - these act on the lumen-facing leaflet even though the protein sits on the cytoplasmic side
#   B  endocytic machinery acting on the CYTOPLASMIC face only: ESCRT, Rab/Arf/SNARE/tethering, phosphoinositide kinases, mTORC1/Ragulator, autophagy,
#      lipid transfer proteins and cytoplasmic-leaflet lipid enzymes, kinases / phosphatases
#   C  not endocytic: transcription factors, ER sterol handling, plasma-membrane receptors, cytosolic IFN effectors and chaperones
# Each PNG has two panels side by side (IFNB vs Mock | IFNG vs Mock) with the genes in the same order in both, grouped by function, sorted by IFNB log2FC.
# Bars = shrunk log2FC (red up, blue down), lines = +/- 1 lfcSE, stars = padj (* < 0.05, ** < 0.01, *** < 0.001, ns). Mouse symbols are matched ignoring case; a few human
# names are mapped to the mouse gene (ALIAS). Genes removed by the pre-filter in 04 are listed as "Not found" and genes with baseMean < MIN_BASEMEAN are left off the figures
# (still in the CSV). Writes PNG + PDF to requant_v2_outputs/figures_lipid_endocytic_genes/ plus lipid_endocytic_genes_log2fc.csv and genes_not_found.csv
# Usage: Rscript 70_gene_bars_endocytic_lipid_genes.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter

MIN_BASEMEAN <- 10
INK <- "#1F2933"; RED <- "#C53030"; BLUE <- "#2B6CB0"
COMP <- c(IFNB_vs_Mock = "IFNB vs Mock", IFNY_vs_Mock = "IFNG vs Mock")
ALIAS <- list(GBA1 = c("Gba1", "Gba"), CERT1 = c("Cert1", "Col4a3bp"), OAS1 = "Oas1*", IRGM = "Irgm*", RUBCN = c("Rubcn", "Kiaa0226"), LAMTOR1 = c("Lamtor1", "C11orf59"),
              LAMTOR2 = c("Lamtor2", "Rbsn"), LAMTOR3 = c("Lamtor3", "Mapkap1"), SDCBP = c("Sdcbp", "Sdcbp1"), VPS37A = c("Vps37a", "Hcrp1"))   # trailing * = every symbol starting with that
SPEC <- list(
  A = list(file = "A_endolysosomal_lumen_facing", title = "Endocytic, in or facing the lumen", sets = list(
    `Lumenal hydrolases and lipid-degrading enzymes` = c("ASAH1", "CTSB", "CTSD", "CTSL", "GBA1", "GLB1", "GM2A", "HEXA", "HEXB", "LIPA", "PSAP", "SMPD1", "NPC2"),
    `Endosomal / lysosomal membrane proteins facing the lumen` = c("LAMP1", "LAMP2", "SCARB2", "CD63", "CD81", "CD9", "TSPAN6", "NPC1", "IFITM1", "IFITM2", "IFITM3", "CYBB", "MCOLN1", "TPCN1", "TPCN2", "SLC38A9", "TLR4"),
    `V-ATPase (lumen acidification)` = c("ATP6AP1", "ATP6V0A1", "ATP6V0D1", "ATP6V1A", "ATP6V1B2"),
    `Lipid floppases / exporters (lumen-facing leaflet)` = c("ABCA1", "ABCG1"))),
  B = list(file = "B_endocytic_cytoplasmic_side", title = "Endocytic machinery, cytoplasmic face only", sets = list(
    `ESCRT / MVB sorting` = c("CHMP1A", "CHMP1B", "CHMP2A", "CHMP2B", "CHMP3", "CHMP4A", "CHMP4B", "CHMP4C", "CHMP5", "CHMP6", "HGS", "STAM", "STAM2", "TSG101", "VPS25", "VPS28", "SNF8", "VPS36", "MVB12A", "VPS37A", "IST1", "VTA1", "VPS4A", "VPS4B", "PDCD6IP", "SDCBP"),
    `Rab / Arf GTPases, tethering and SNAREs` = c("RAB4A", "RAB5A", "RAB5B", "RAB5C", "RAB7A", "RAB9A", "RAB11A", "ARF6", "ARL8B", "EEA1", "RILP", "WDFY2", "ZFYVE1", "ZFYVE16", "VPS11", "VPS16", "VPS18", "VPS33A", "VPS39", "VPS41", "STX7", "STX8", "VAMP7", "VAMP8", "VTI1B", "SYT7", "UNC13A", "UNC13B", "SNX1", "SNX2", "SNX3", "SNX5", "SNX6", "SNX27"),
    `Phosphoinositide kinases / phosphatases` = c("PIK3C3", "PIK3R4", "UVRAG", "PIKFYVE", "VAC14", "FIG4", "PI4K2A", "PI4K2B", "PI4KA", "PI4KB", "PIP5K1C"),
    `mTORC1 / Ragulator (lysosomal nutrient sensing)` = c("MTOR", "RPTOR", "MLST8", "RHEB", "RRAGA", "RRAGB", "RRAGC", "RRAGD", "LAMTOR1", "LAMTOR2", "LAMTOR3", "LAMTOR4", "LAMTOR5", "FLCN", "FNIP1", "FNIP2", "TSC1", "TSC2", "GPR155"),
    `Autophagy machinery` = c("ATG3", "ATG5", "ATG7", "ATG10", "ATG12", "ATG14", "ATG16L1", "BECN1", "GABARAP", "GABARAPL1", "GABARAPL2", "MAP1LC3A", "MAP1LC3B", "MAP1LC3C", "RUBCN", "TECPR1"),
    `Lipid transfer proteins and cytoplasmic-leaflet lipid enzymes` = c("OSBP", "OSBPL1A", "VAPA", "VAPB", "CERT1", "CPTP", "STARD3", "STARD3NL", "GRAMD1A", "GRAMD1B", "GRAMD1C", "PLD1", "PLD2", "SPHK1", "SPHK2", "SMPD3", "DAGLA", "DAGLB", "MGLL"),
    `Kinases, phosphatases and regulators` = c("PRKCA", "PRKCB", "PRKCD", "PRKCE", "PRKCQ", "PRKCZ", "PRKD1", "PRKD2", "PRKD3", "PPP1CA", "PPP2CA", "PPP2CB", "PPP3CA", "RAF1", "KSR1", "RASGRP1", "RASGRP3", "CHN1", "CHN2", "CAMKK2", "LRRK2", "NCF1", "NCF2", "NCF4"))),
  C = list(file = "C_not_endocytic", title = "Not endocytic (nuclear, ER, plasma membrane, cytosolic IFN effectors)", sets = list(
    `Transcription factors and IFN signalling` = c("STAT1", "STAT2", "IRF1", "IRF9", "TFEB", "TFE3", "MITF", "NR1H2", "NR1H3", "SREBF2", "SET"),
    `ER sterol / lipid handling` = c("SCAP", "INSIG1", "INSIG2", "SOAT1", "SOAT2", "CH25H", "CYP7B1", "SGPL1"),
    `Plasma-membrane receptors and lipid-signalling GPCRs` = c("S1PR1", "S1PR2", "S1PR3", "S1PR4", "S1PR5", "LPAR1", "LPAR2", "LPAR3", "LPAR4", "LPAR5", "LPAR6", "GPR4", "GPR55", "GPR132", "NTRK1", "CD40", "FAS", "CAV1", "FLOT1", "FLOT2", "SCARB1", "LDLR", "ENPP2"),
    `IFN-stimulated effectors and chaperones` = c("ISG15", "MX1", "OAS1", "GBP1", "GBP2", "GBP3", "GBP4", "GBP5", "GBP6", "GBP7", "IRGM", "HSPA1A", "HSPA1B", "HSPA8")))
)

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_lipid_endocytic_genes"); dir.create(OUT, showWarnings = FALSE)
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
write_csv(D |> arrange(group, set, gene_name, comparison) |> select(group, set, gene_name, gene_id, baseMean, comparison, log2FC = lfc, lfcSE = se, padj), file.path(OUT, "lipid_endocytic_genes_log2fc.csv"))
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
