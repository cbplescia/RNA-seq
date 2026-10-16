#!/usr/bin/env Rscript
# Genes of the KEGG "Calcium signaling pathway" (mmu04020) among the IFNB responders.
#
# 1. <TAG>_calcium_signaling_up_down_genes.csv   the up/down genes (from 28's annotation table)
#    that KEGG places in the pathway - the list you asked for
# 2. <TAG>_calcium_signaling_all_pathway_genes.csv   EVERY tested gene in the pathway, with
#    log2FC, padj and a significant/not flag (needs internet to fetch KEGG's member list)
# 3. <TAG>_calcium_signaling_volcano.png/.pdf   volcano of the pathway genes: red = significantly
#    up, blue = down (padj < 0.05, |log2FC| > 1), grey = not significant; significant genes
#    labelled; all other tested genes shown faintly behind for context
#
# Usage:  Rscript 31_calcium_signaling.R [TAG]
#   TAG = non_ifn_antiviral (default), non_canonical_isgs, or all_enzymes  - i.e. the gene set
#   whose 28_lipid_function_analysis.R run you are following up. Run 28 for that TAG first.

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })
for (pkg in c("ggrepel", "clusterProfiler", "org.Mm.eg.db"))
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install ", pkg)
select <- dplyr::select; first <- dplyr::first; rename <- dplyr::rename
slice <- dplyr::slice; filter <- dplyr::filter; desc <- dplyr::desc; count <- dplyr::count

PATHWAY_NAME <- "Calcium signaling"     # matched inside the KEGG pathway names
PATHWAY_ID   <- "mmu04020"
LFC_CUT <- 1; PADJ_CUT <- 0.05; MAX_LABELS <- 40
INK <- "#1F2933"
COL <- c(Up = "#C53030", Down = "#2B6CB0", `Not significant` = "#B8B8B8")

TAG <- "non_ifn_antiviral"
a <- commandArgs(trailingOnly = TRUE); if (length(a)) TAG <- a[1]
background_file <- c(
  non_ifn_antiviral  = "DE_v2_IFNB_vs_Mock_protein_coding_without_ISG_and_antiviral_genes.csv",
  non_canonical_isgs = "DE_v2_IFNB_vs_Mock_protein_coding_without_canonical_type_I_ISGs.csv",
  all_enzymes        = "DE_v2_IFNB_vs_Mock_protein_coding_enzymes.csv")[[TAG]]
if (is.null(background_file)) stop("Unknown TAG '", TAG, "'")

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs")
RES <- file.path(OUT, "lipid_function_analysis")

# ---- 1. responders in the pathway (from the 28 annotation table) -------------------------------
annf <- file.path(RES, paste0(TAG, "_gene_annotations.csv"))
if (!file.exists(annf)) stop("Missing ", annf, " - run  Rscript 28_lipid_function_analysis.R ", TAG, "  first")
ann <- read_csv(annf, show_col_types = FALSE)
ca <- ann |> filter(!is.na(KEGG_pathways), grepl(PATHWAY_NAME, KEGG_pathways, ignore.case = TRUE)) |>
  arrange(desc(log2FoldChange))
write_csv(ca, file.path(RES, paste0(TAG, "_calcium_signaling_up_down_genes.csv")))
message(sprintf("%s: %d responding genes are in the KEGG '%s' pathway (%d up, %d down)",
                TAG, nrow(ca), PATHWAY_NAME, sum(ca$direction == "up"), sum(ca$direction == "down")))
message("Up:   ", paste(ca$gene_name[ca$direction == "up"], collapse = ", "))
message("Down: ", paste(ca$gene_name[ca$direction == "down"], collapse = ", "))

# ---- 2. all tested genes in the pathway ------------------------------------------------------
bg <- read_csv(file.path(OUT, background_file), show_col_types = FALSE)
bg$.key <- sub("\\.[0-9]+$", "", bg$gene_id)
members <- tryCatch({
  kk <- clusterProfiler::download_KEGG("mmu")
  p2g <- kk$KEGGPATHID2EXTID
  unique(p2g$to[p2g$from %in% c(PATHWAY_ID, sub("^mmu", "", PATHWAY_ID))])
}, error = function(e) { message("Could not fetch KEGG membership (", conditionMessage(e), ")"); NULL })

if (!is.null(members) && length(members)) {
  map <- suppressMessages(suppressWarnings(clusterProfiler::bitr(unique(bg$.key), fromType = "ENSEMBL",
                                                                 toType = "ENTREZID", OrgDb = org.Mm.eg.db::org.Mm.eg.db))) |>
    distinct(ENSEMBL, .keep_all = TRUE)
  bg$in_pathway <- map$ENTREZID[match(bg$.key, map$ENSEMBL)] %in% members
  message(sprintf("KEGG lists %d genes in the pathway; %d of them are in the tested gene set", length(members), sum(bg$in_pathway)))
} else {
  message("Falling back to the responding genes only for the volcano")
  bg$in_pathway <- bg$gene_id %in% ca$gene_id
}
bg <- bg |> mutate(
  status = case_when(!is.na(padj) & padj < PADJ_CUT & log2FoldChange >  LFC_CUT ~ "Up",
                     !is.na(padj) & padj < PADJ_CUT & log2FoldChange < -LFC_CUT ~ "Down",
                     TRUE ~ "Not significant"),
  y = -log10(pmax(padj, 1e-300)))
pw <- bg |> filter(in_pathway, !is.na(padj))
write_csv(pw |> select(-.key, -y, -in_pathway) |> arrange(desc(log2FoldChange)),
          file.path(RES, paste0(TAG, "_calcium_signaling_all_pathway_genes.csv")))
message(sprintf("Pathway genes tested: %d  (significant up %d, down %d)", nrow(pw), sum(pw$status == "Up"), sum(pw$status == "Down")))

# ---- 3. volcano ---------------------------------------------------------------------------------
pw$status <- factor(pw$status, levels = names(COL))
lab <- pw |> filter(status != "Not significant") |> slice_max(abs(log2FoldChange), n = MAX_LABELS, with_ties = FALSE)
ctx <- bg |> filter(!in_pathway, !is.na(padj))
xr <- max(abs(bg$log2FoldChange), na.rm = TRUE) * 1.05
p <- ggplot(pw, aes(log2FoldChange, y)) +
  geom_hline(yintercept = -log10(PADJ_CUT), linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_point(data = ctx, colour = "#E6E6E6", size = 0.5, alpha = 0.6) +
  geom_point(aes(colour = status), size = 2, alpha = 0.85) +
  ggrepel::geom_text_repel(data = lab, aes(label = gene_name), colour = INK, size = 3, max.overlaps = Inf,
                           min.segment.length = 0, segment.colour = "grey60", segment.size = 0.25,
                           box.padding = 0.4, seed = 1) +
  scale_colour_manual(values = COL, name = NULL, drop = FALSE) +
  scale_x_continuous(limits = c(-xr, xr)) +
  labs(title = sprintf("KEGG %s pathway genes: IFNB vs Mock", PATHWAY_NAME),
       subtitle = sprintf("%s gene set; %d pathway genes tested (faint grey = all other tested genes); coloured when padj<%g & |log2FC|>%g%s",
                          TAG, nrow(pw), PADJ_CUT, LFC_CUT,
                          if (nrow(lab) == MAX_LABELS) sprintf("; top %d by |log2FC| labelled", MAX_LABELS) else ""),
       x = "log2 fold-change (IFNB vs Mock, shrunk)", y = expression(-log[10](padj))) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey94", linewidth = 0.3),
        panel.border = element_rect(colour = "grey70"),
        axis.text = element_text(colour = INK), axis.title = element_text(colour = INK),
        plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 8.5),
        legend.position = "bottom")
base <- file.path(RES, paste0(TAG, "_calcium_signaling_volcano"))
ggsave(paste0(base, ".png"), p, width = 8, height = 7, dpi = 300)
ggsave(paste0(base, ".pdf"), p, width = 8, height = 7)
message("\nSaved to: ", RES)
