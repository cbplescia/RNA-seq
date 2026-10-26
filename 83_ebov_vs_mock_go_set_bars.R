#!/usr/bin/env Rscript
# Bar plots of the EBOV vs Mock differential expression of three GO-defined gene sets (protein-coding genes; GO annotation from org.Mm.eg.db, which includes every child term):
#   (1) calcium_binding     all genes annotated to "calcium ion binding"                       GO:0005509
#   (2) ion_channels        all genes annotated to "ion channel activity" (monoatomic ion channel activity)  GO:0005216, split by channel type
#   (3) endolysosomal       all genes annotated to endosome / lysosome compartments (endosome GO:0005768, lysosome GO:0005764, late endosome GO:0005770, multivesicular body GO:0005771,
#                           endolysosome GO:0036019), split by sub-compartment (lysosomal lumen / membrane, endosome membrane, late / early endosome, MVB, other)
# These sets are large (hundreds of genes), so a figure shows the genes that are significantly changed (padj < PADJ_PLOT, default 0.05; |log2FC| > LFC_PLOT, default 0 = any fold-change),
# sorted within each group by log2FC and paged at PER_PAGE genes per PNG. The CSV for each set has EVERY gene in the set with its log2FC and padj.
# Bars = shrunk log2FC (red up, blue down), lines = +/- 1 lfcSE, stars = padj (* < 0.05, ** < 0.01, *** < 0.001).
# Writes to requant_v2_outputs/figures_ebov_go_sets/: bars_EBOV_vs_Mock_<set>_p<n>.png/.pdf and EBOV_vs_Mock_<set>_all_genes.csv
# Usage: Rscript 83_ebov_vs_mock_go_set_bars.R   (PROJ from 'source setup_session.sh';  LFC_PLOT=0.58 PER_PAGE=60 optional)

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
for (pkg in c("clusterProfiler", "org.Mm.eg.db")) if (!requireNamespace(pkg, quietly = TRUE)) stop("Install ", pkg)
select <- dplyr::select; filter <- dplyr::filter; slice <- dplyr::slice
OrgDb <- org.Mm.eg.db::org.Mm.eg.db

PADJ_PLOT <- 0.05; LFC_PLOT <- as.numeric(Sys.getenv("LFC_PLOT", unset = "0")); PER_PAGE <- as.integer(Sys.getenv("PER_PAGE", unset = "70"))
INK <- "#1F2933"; RED <- "#C53030"; BLUE <- "#2B6CB0"
SETS <- list(
  calcium_binding = list(title = "Calcium ion binding genes", cats = list(`Calcium ion binding` = "GO:0005509")),
  ion_channels = list(title = "Ion channel genes", cats = list(`Calcium channels` = "GO:0005262", `Potassium channels` = "GO:0005267", `Sodium channels` = "GO:0005272", `Chloride channels` = "GO:0005254",
                                                                `Ligand-gated channels` = "GO:0015276", `Voltage-gated channels` = "GO:0022832", `Other ion channels` = "GO:0005216")),
  endolysosomal = list(title = "Endolysosomal-associated genes", cats = list(`Lysosomal lumen` = "GO:0043202", `Lysosomal membrane` = "GO:0005765", `Late endosome` = "GO:0005770",
                                                                         `Early endosome` = "GO:0005769", `Multivesicular body` = "GO:0005771", `Endosome membrane` = "GO:0010008", `Lysosome (other)` = "GO:0005764",
                                                                         `Endosome (other)` = "GO:0005768", `Endolysosome` = "GO:0036019")))
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_ebov_go_sets"); dir.create(OUT, showWarnings = FALSE)
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

pcf <- file.path(RNA, "DE_v2_EBOV_vs_Mock_protein_coding.csv")
d <- if (file.exists(pcf)) read_csv(pcf, show_col_types = FALSE) else read_csv(file.path(RNA, "DE_v2_EBOV_vs_Mock_annotated.csv"), show_col_types = FALSE) |> filter(gene_biotype == "protein_coding")
d <- d |> filter(!is.na(padj), grepl("^ENSMUSG", gene_id)) |> mutate(ens = sub("\\..*$", "", gene_id), label = ifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name))
message(sprintf("%d protein-coding genes tested in EBOV vs Mock", nrow(d)))

g2a <- as.list(org.Mm.eg.db::org.Mm.egGO2ALLEGS)
ens_of <- function(go) {
  e <- unique(as.character(unlist(g2a[go]))); e <- e[!is.na(e)]; if (!length(e)) return(character())
  m <- suppressMessages(suppressWarnings(clusterProfiler::bitr(e, "ENTREZID", "ENSEMBL", OrgDb))); unique(m$ENSEMBL)
}
for (s in names(SETS)) {
  S <- SETS[[s]]; message("\n== ", S$title)
  members <- lapply(S$cats, ens_of)
  if (!any(lengths(members) > 0)) { message("   no genes returned for these GO terms - skipped"); next }
  assigned <- character(); cat_of <- setNames(rep(NA_character_, 0), character())
  for (cn in names(members)) { g <- setdiff(members[[cn]], assigned); cat_of[g] <- cn; assigned <- c(assigned, g) }     # first matching category in the order above wins
  D <- d |> filter(ens %in% names(cat_of)) |> mutate(cat = factor(cat_of[ens], levels = names(S$cats)))
  message(sprintf("   %d genes in the set among the %d tested; %d with padj < %g", nrow(D), nrow(d), sum(D$padj < PADJ_PLOT), PADJ_PLOT))
  write_csv(D |> arrange(cat, desc(log2FoldChange)) |> select(category = cat, gene_id, gene_name, baseMean, log2FC = log2FoldChange, lfcSE, padj), file.path(OUT, paste0("EBOV_vs_Mock_", s, "_all_genes.csv")))
  P <- D |> filter(padj < PADJ_PLOT, abs(log2FoldChange) > LFC_PLOT)
  message(sprintf("   drawing %d genes (padj < %g, |log2FC| > %g)", nrow(P), PADJ_PLOT, LFC_PLOT))
  if (!nrow(P)) next
  P <- P |> arrange(cat, desc(log2FoldChange)) |> mutate(page = ceiling(row_number() / PER_PAGE), dir = ifelse(log2FoldChange >= 0, "up", "down"),
                                                        key = paste0(cat, "||", label), hj = ifelse(log2FoldChange >= 0, 0, 1), xs = log2FoldChange + ifelse(log2FoldChange >= 0, 1, -1) * (lfcSE + 0.04))
  np <- max(P$page)
  for (pg in seq_len(np)) {
    Q <- P |> filter(page == pg); Q$key <- factor(Q$key, levels = rev(Q$key)); Q$cat <- droplevels(Q$cat)
    p <- ggplot(Q, aes(log2FoldChange, key, fill = dir)) + geom_col(width = 0.75) +
      geom_segment(aes(x = log2FoldChange - lfcSE, xend = log2FoldChange + lfcSE, yend = key), colour = INK, linewidth = 0.3) +
      geom_text(aes(x = xs, label = stars(padj), hjust = hj), colour = INK, size = 2.3, vjust = 0.45) + geom_vline(xintercept = 0, colour = INK, linewidth = 0.4) +
      facet_grid(cat ~ ., scales = "free_y", space = "free_y", labeller = labeller(cat = label_wrap_gen(16))) + scale_fill_manual(values = c(up = RED, down = BLUE), guide = "none") +
      scale_y_discrete(labels = function(x) sub("^.*\\|\\|", "", x)) + scale_x_continuous(expand = expansion(mult = 0.1)) +
      labs(title = sprintf("EBOV vs Mock: %s", S$title), subtitle = sprintf("%d of %d tested genes in the set are shown (padj < %g%s); shrunk log2FC +/- 1 SE; stars = padj; grouped by %s, sorted by log2FC%s",
                                                                        nrow(P), nrow(D), PADJ_PLOT, if (LFC_PLOT > 0) sprintf(", |log2FC| > %g", LFC_PLOT) else "", if (length(S$cats) > 1) "sub-category" else "set",
                                                                        if (np > 1) sprintf("; page %d of %d", pg, np) else ""), x = "log2 fold-change (EBOV vs Mock)", y = NULL) +
      theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3), panel.border = element_rect(colour = "grey70"),
            axis.text = element_text(colour = INK), axis.text.y = element_text(face = "italic", size = 7.5), plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 8),
            strip.background = element_rect(fill = "grey95", colour = "grey70"), strip.text = element_text(colour = INK, face = "bold", size = 8), strip.text.y = element_text(angle = 0, hjust = 0), panel.spacing.y = unit(0.15, "lines"))
    h <- max(3.5, 0.16 * nrow(Q) + 0.4 * n_distinct(Q$cat) + 1.8); nm <- sprintf("bars_EBOV_vs_Mock_%s_p%d", s, pg)
    ggsave(file.path(OUT, paste0(nm, ".png")), p, width = 8, height = h, dpi = 300, limitsize = FALSE); ggsave(file.path(OUT, paste0(nm, ".pdf")), p, width = 8, height = h, limitsize = FALSE)
  }
  message(sprintf("   %d page(s) written", np))
}
message("\nSaved to: ", OUT)
