#!/usr/bin/env Rscript
# Bar graphs of the differential expression of host lipid-metabolism, lipid-modification and lipid-transport genes in EBOV vs Mock (426 genes in five groups; my classification -
# check it and move genes between groups in SPEC below):
#   A  endolysosomal, in or facing the LUMEN (lysosomal hydrolases and membrane proteins)      B  lipid TRANSLOCASES (flippases, floppases, scramblases)
#   C  endocytic / contact-site machinery acting on the CYTOPLASMIC face (PI kinases / phosphatases, lipid transfer proteins, cytoplasmic-leaflet enzymes)
#   D  not endocytic: ER / Golgi lipid synthesis (cholesterol, phospholipids, sphingolipids, fatty acids) and SREBP / LXR / PPAR regulators
#   E  not endocytic: mitochondria, peroxisome, lipid droplets, plasma-membrane receptors, secreted lipases, eicosanoid and redox enzymes
# One PNG per group, genes grouped by function and sorted by EBOV log2FC. Bars = shrunk log2FC (red up, blue down), lines = +/- 1 lfcSE, stars = padj (* < 0.05, ** < 0.01, *** < 0.001, ns).
# Genes removed by the pre-filter in 04 are listed as "Not found"; genes with baseMean < MIN_BASEMEAN are left off the figures (still in the CSV).
# SIG_ONLY=1 draws only the genes with padj < 0.05 (a shorter, focused set of figures); the CSV always has every gene.
# Writes PNG + PDF to requant_v2_outputs/figures_ebov_lipid_genes/ plus ebov_lipid_genes_log2fc.csv and genes_not_found.csv
# Usage: Rscript 80_gene_bars_lipid_genes_ebov_vs_mock.R   (PROJ from 'source setup_session.sh';  SIG_ONLY=1 optional)

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter

MIN_BASEMEAN <- 10
INK <- "#1F2933"; RED <- "#C53030"; BLUE <- "#2B6CB0"
COMP <- c(EBOV_vs_Mock = "EBOV vs Mock"); SIG_ONLY <- identical(Sys.getenv("SIG_ONLY"), "1")
ALIAS <- list(GBA1 = c("Gba1", "Gba"), CERT1 = c("Cert1", "Col4a3bp"), ABCB1 = "Abcb1*", COQ8A = c("Coq8a", "Adck3"), COQ8B = c("Coq8b", "Adck4"), PRELID3B = c("Prelid3b", "Slmo2"),
              PEDS1 = c("Peds1", "Tmem189"), TAFAZZIN = c("Tafazzin", "Taz"), PLPP1 = c("Plpp1", "Ppap2a"), PLPP2 = c("Plpp2", "Ppap2c"), PLPP3 = c("Plpp3", "Ppap2b"),
              SELENOI = c("Selenoi", "Ept1"), TMEM24 = c("Tmem24", "C2cd2l"))   # trailing * = every symbol starting with that
SPEC <- list(
  A = list(file = "A_endolysosomal_lumen_facing", title = "Endolysosomal: in or facing the lumen", sets = list(
    `Lysosomal hydrolases and lumenal enzymes` = c("ARSA", "ASAH1", "CTSA", "GALC", "GBA1", "GLB1", "GM2A", "GRN", "HEXA", "HEXB", "LIPA", "NEU1", "PLA2G15", "PLD3", "PLD4", "PSAP", "SMPD1", "NPC2"),
    `Lysosomal / endosomal membrane and lumen-facing proteins` = c("CLN3", "CLN5", "LAMP2", "SCARB2", "NPC1", "TMEM106B", "SPNS1", "PLPP1", "PLPP2", "PLPP3"))),
  B = list(file = "B_lipid_translocases", title = "Lipid translocases (flippases, floppases, scramblases): change leaflet asymmetry incl. the lumen-facing leaflet", sets = list(
    `Flippases (P4-ATPases) and CDC50 subunits` = c("ATP8A1", "ATP8A2", "ATP8B1", "ATP8B2", "ATP8B3", "ATP8B4", "ATP9A", "ATP9B", "ATP10A", "ATP10B", "ATP10D", "ATP11A", "ATP11B", "ATP11C", "TMEM30A", "TMEM30B", "TMEM30C"),
    `Floppases (ABC transporters)` = c("ABCA1", "ABCA3", "ABCA4", "ABCA7", "ABCA12", "ABCB1", "ABCB4", "ABCC1", "ABCG1", "ABCG2", "ABCG4", "ABCG5", "ABCG8"),
    `Scramblases and lipid-transporting membrane proteins` = c("ANO3", "ANO4", "ANO5", "ANO6", "ANO7", "ANO9", "ANO10", "PLSCR1", "PLSCR2", "PLSCR3", "PLSCR4", "XKR4", "XKR8", "XKR9", "TMEM63A", "TMEM63B", "ATG9A", "ATG9B", "TMEM41B", "VMP1", "SPNS2", "MFSD2A", "MFSD2B"))),
  C = list(file = "C_endocytic_cytoplasmic_side", title = "Endocytic / membrane-contact machinery acting on the cytoplasmic face", sets = list(
    `Phosphoinositide kinases and phosphatases` = c("FIG4", "INPP4A", "INPP4B", "INPP5B", "INPP5D", "INPP5E", "INPP5F", "INPPL1", "MTM1", "MTMR1", "MTMR2", "MTMR3", "MTMR4", "MTMR6", "OCRL", "PI4K2A", "PI4K2B", "PI4KA", "PI4KB", "PIK3C2A", "PIK3C2B", "PIK3C2G", "PIK3C3", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIKFYVE", "PIP4K2A", "PIP4K2B", "PIP4K2C", "PIP5K1A", "PIP5K1B", "PIP5K1C", "PTEN", "SACM1L", "SYNJ1", "SYNJ2"),
    `Lipid transfer and contact-site proteins` = c("ATG2A", "ATG2B", "CERT1", "CPTP", "ESYT1", "ESYT2", "ESYT3", "GLTP", "GRAMD1A", "GRAMD1B", "GRAMD1C", "OSBP", "OSBPL1A", "OSBPL2", "OSBPL3", "OSBPL5", "OSBPL8", "OSBPL9", "OSBPL10", "OSBPL11", "PITPNA", "PITPNB", "PITPNM1", "PITPNM2", "PITPNM3", "STARD2", "STARD3", "STARD4", "STARD5", "STARD10", "TMEM24", "VPS13A", "VPS13B", "VPS13C", "VPS13D"),
    `Autophagy conjugation` = c("ATG3", "ATG7"),
    `Cytoplasmic-leaflet lipid enzymes and signalling` = c("PLD1", "PLD2", "PLCB1", "PLCB2", "PLCB3", "PLCD1", "PLCE1", "PLCG1", "PLCG2", "DAGLA", "DAGLB", "MGLL", "ABHD6", "DGKA", "DGKB", "DGKD", "DGKE", "DGKG", "DGKH", "DGKI", "DGKK", "DGKQ", "DGKZ", "SPHK1", "SPHK2", "PLA2G4A", "PLA2G6", "DDHD1", "DDHD2", "LRRK2"))),
  D = list(file = "D_ER_Golgi_lipid_synthesis", title = "Not endocytic: ER / Golgi lipid synthesis and sterol regulation", sets = list(
    `Cholesterol and isoprenoid synthesis, sterol handling` = c("CYP51A1", "CYP46A1", "DHCR24", "DHCR7", "EBP", "FDFT1", "FDPS", "HMGCR", "HMGCS1", "HSD17B7", "LSS", "MSMO1", "MVD", "MVK", "NSDHL", "PMVK", "SC5D", "SQLE", "TM7SF2", "LBR", "ARV1", "SOAT1", "SOAT2", "NCEH1", "CES1", "CH25H", "CYP7B1"),
    `SREBP / LXR / PPAR regulators` = c("SCAP", "INSIG1", "INSIG2", "MBTPS1", "MBTPS2", "SREBF1", "SREBF2", "NR1H2", "NR1H3", "PPARA", "PPARG", "PPARGC1A"),
    `Phospholipid and neutral-lipid synthesis (ER / Golgi)` = c("AGPAT1", "AGPAT2", "AGPAT3", "AGPAT4", "AGPAT5", "CDIPT", "CDS1", "CDS2", "CEPT1", "CHKA", "CHKB", "CHPT1", "PCYT1A", "PCYT1B", "PCYT2", "PEMT", "PTDSS1", "PTDSS2", "SELENOI", "LPCAT1", "LPCAT2", "LPCAT3", "LPCAT4", "LPGAT1", "LCLAT1", "MBOAT1", "MBOAT2", "MBOAT7", "GPAT2", "GPAT3", "GPAT4", "LPIN1", "LPIN2", "LPIN3", "DGAT1", "DGAT2", "MOGAT1", "MOGAT2", "PEDS1", "GPCPD1", "GDPD1", "GDPD3", "ABHD12", "PNPLA6", "PNPLA7"),
    `Sphingolipid synthesis and metabolism` = c("ACER1", "ACER2", "ACER3", "B3GALT4", "B4GALNT1", "B4GALT5", "B4GALT6", "CERK", "CERS1", "CERS2", "CERS3", "CERS4", "CERS5", "CERS6", "DEGS1", "DEGS2", "FA2H", "GAL3ST1", "GBA2", "KDSR", "ORMDL1", "ORMDL2", "ORMDL3", "SAMD8", "SGMS1", "SGMS2", "SGPL1", "SGPP1", "SGPP2", "SMPD2", "SMPD3", "SMPD4", "SPTLC1", "SPTLC2", "SPTLC3", "SPTSSA", "SPTSSB", "ST3GAL2", "ST3GAL5", "ST8SIA1", "UGCG", "UGT8"),
    `Fatty-acid activation, elongation and desaturation` = c("ACSL1", "ACSL3", "ACSL4", "ACSL5", "ACSL6", "ELOVL1", "ELOVL2", "ELOVL3", "ELOVL4", "ELOVL5", "ELOVL6", "ELOVL7", "FADS1", "FADS2", "FADS3", "SCD", "SCD5"))),
  E = list(file = "E_mito_peroxisome_LD_PM_other", title = "Not endocytic: mitochondria, peroxisome, lipid droplets, plasma membrane and secreted", sets = list(
    `Mitochondrial lipid metabolism` = c("ACADM", "ACADVL", "AGK", "COQ10A", "COQ10B", "COQ2", "COQ3", "COQ4", "COQ5", "COQ6", "COQ7", "COQ8A", "COQ8B", "COQ9", "PDSS1", "PDSS2", "CPT1A", "CPT1B", "CPT1C", "CPT2", "CRAT", "CROT", "HADHA", "HADHB", "ECI1", "SLC25A20", "SLC22A5", "CRLS1", "PGS1", "PISD", "PTPMT1", "TAMM41", "TAFAZZIN", "GPAT1", "PRELID1", "PRELID3B", "TRIAP1", "PPTC7", "CYP27A1", "PNPLA8"),
    `Peroxisome` = c("ABCD1", "ABCD2", "ABCD3", "PEX1", "PEX12", "PEX14", "PEX19", "PEX2", "PEX3", "PEX5", "PEX7", "GNPAT", "AGPS", "FAR1", "FAR2"),
    `Lipid droplets and lipolysis` = c("BSCL2", "CIDEA", "CIDEB", "CIDEC", "FITM1", "FITM2", "G0S2", "PLIN1", "PLIN2", "PLIN3", "PLIN4", "PLIN5", "PNPLA2", "LIPE", "ABHD5", "FABP4", "FABP5"),
    `Plasma-membrane receptors, secreted lipases, eicosanoid and redox enzymes` = c("CD36", "BSG", "NPTN", "LDLR", "SCARB1", "GPR55", "ENPP2", "ENPP7", "ASAH2", "LCAT", "LPL", "PLA1A", "PLA2G10", "PLA2G2A", "PLA2G5", "PLA2G7", "PTGS1", "PTGS2", "ALOX5", "ALOX12", "ALOX15", "AIFM2", "GPX4", "NEU2", "NEU3", "NEU4")))
)

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_ebov_lipid_genes"); dir.create(OUT, showWarnings = FALSE)
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
write_csv(D |> arrange(group, set, gene_name, comparison) |> select(group, set, gene_name, gene_id, baseMean, comparison, log2FC = lfc, lfcSE = se, padj), file.path(OUT, "ebov_lipid_genes_log2fc.csv"))
low <- D |> group_by(gene_id, gene_name) |> summarise(bm = max(baseMean, na.rm = TRUE), .groups = "drop") |> filter(bm < MIN_BASEMEAN)
if (nrow(low)) message("Left off the figures (baseMean < ", MIN_BASEMEAN, "): ", paste(low$gene_name, collapse = ", "))
D <- D |> filter(!gene_id %in% low$gene_id)
if (SIG_ONLY) { D <- D |> filter(!is.na(padj), padj < 0.05); message("SIG_ONLY: ", n_distinct(D$gene_id), " genes with padj < 0.05 drawn") }

# ---- one PNG per group ----
for (k in names(SPEC)) {
  d <- D |> filter(group == k)
  if (!nrow(d)) { message("Group ", k, ": no genes to plot"); next }
  ord <- d |> filter(comparison == COMP[["EBOV_vs_Mock"]]) |> group_by(set, gene_id) |> summarise(o = max(lfc, na.rm = TRUE), .groups = "drop")
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
    labs(title = SPEC[[k]]$title, subtitle = "shrunk log2FC +/- 1 SE vs Mock; stars = padj (ns = not significant); genes grouped by function, sorted by EBOV log2FC",
         x = "log2 fold-change vs Mock", y = NULL) +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
      panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK), axis.text.y = element_text(face = "italic", size = if (n > 80) 6.5 else 8),
      plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 8.5), strip.background = element_rect(fill = "grey95", colour = "grey70"),
      strip.text = element_text(colour = INK, face = "bold", size = 8), strip.text.y = element_text(angle = 0, hjust = 0), panel.spacing.y = unit(0.15, "lines"))
  h <- min(48, max(4, 0.15 * n + 0.45 * n_distinct(d$set) + 1.8))
  ggsave(file.path(OUT, paste0("bars_EBOV_vs_Mock_", SPEC[[k]]$file, ".png")), p, width = 7.5, height = h, dpi = 300, limitsize = FALSE)
  ggsave(file.path(OUT, paste0("bars_EBOV_vs_Mock_", SPEC[[k]]$file, ".pdf")), p, width = 7.5, height = h, limitsize = FALSE)
  message(sprintf("Group %s: %d genes plotted", k, n))
}
message("Saved to: ", OUT)
