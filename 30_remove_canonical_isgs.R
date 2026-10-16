#!/usr/bin/env Rscript
# IFNB vs Mock, protein-coding genes, with the KNOWN canonical type I ISGs removed.
#
# "Canonical type I ISG" = a gene in ANY of these curated sets (mouse orthologs via msigdbr):
#   hallmark  MSigDB Hallmark HALLMARK_INTERFERON_ALPHA_RESPONSE
#   reactome  Reactome "Interferon alpha/beta signaling"
#             (REACTOME_INTERFERON_ALPHA_BETA_SIGNALING)
#   extra     EXTRA_ISGS below: classic mouse ISGs added by hand (edit freely)
# A column records which set(s) flagged each removed gene so you can audit the call.
#
# Caveat: the Hallmark and Reactome sets are human-derived, so mouse-specific ISGs (e.g.
# members of expanded mouse families) may be missing and stay in the "remaining" list.
# "Not on the list" means not a KNOWN canonical ISG, not that a gene is unresponsive.
#
# Input : DE_v2_IFNB_vs_Mock_protein_coding.csv   (shrunk log2FC; from 07 + 08)
# Output: DE_v2_IFNB_vs_Mock_protein_coding_without_canonical_type_I_ISGs.csv   (what you asked for)
#         DE_v2_IFNB_vs_Mock_protein_coding_canonical_type_I_ISGs_removed.csv   (the removed genes + source)
# Needs the msigdbr package:  conda install -c conda-forge r-msigdbr   (or install.packages("msigdbr"))

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr) })

EXTRA_ISGS <- c("Mx1", "Mx2", "Isg15", "Isg20", "Ifit1", "Ifit2", "Ifit3", "Oas1a", "Oas2", "Oas3",
                "Oasl1", "Oasl2", "Rsad2", "Usp18", "Irf7", "Irf9", "Stat1", "Stat2", "Ifi44",
                "Ifi27l2a", "Ddx58", "Ifih1", "Bst2", "Ifitm3", "Cxcl10", "Zbp1", "Eif2ak2",
                "Herc6", "Trim30a", "Samhd1", "Ifi35", "Cmpk2", "Parp9", "Dtx3l")
LFC_CUT <- 1; PADJ_CUT <- 0.05

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs")

d <- read_csv(file.path(OUT, "DE_v2_IFNB_vs_Mock_protein_coding.csv"), show_col_types = FALSE)
d$.key <- sub("\\.[0-9]+$", "", d$gene_id)
d$.sym <- toupper(d$gene_name)
message("Protein-coding genes in IFNB vs Mock: ", nrow(d))

# ---- fetch the curated sets (handles old and new msigdbr interfaces) ------------------------------
get_set <- function(collection, subcollection, set_name) {
  m <- tryCatch(msigdbr::msigdbr(species = "Mus musculus", collection = collection, subcollection = subcollection),
                error = function(e) tryCatch(msigdbr::msigdbr(species = "Mus musculus", category = collection,
                                                              subcategory = subcollection),
                                             error = function(e2) NULL))
  if (is.null(m)) return(NULL)
  m <- m |> filter(gs_name == set_name)
  if (!nrow(m)) return(NULL)
  list(sym = toupper(unique(m$gene_symbol)),
       ens = if ("ensembl_gene" %in% names(m)) unique(na.omit(m$ensembl_gene)) else character(0))
}
sets <- list()
if (requireNamespace("msigdbr", quietly = TRUE)) {
  sets$hallmark <- get_set("H", NULL, "HALLMARK_INTERFERON_ALPHA_RESPONSE")
  sets$reactome <- get_set("C2", "CP:REACTOME", "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING")
  for (n in c("hallmark", "reactome"))
    if (is.null(sets[[n]])) message("WARNING: could not load the ", n, " set from msigdbr - it is skipped")
    else message(sprintf("%-9s %3d mouse genes", n, length(sets[[n]]$sym)))
} else message("WARNING: msigdbr is not installed - using ONLY the EXTRA_ISGS list (install r-msigdbr for the curated sets)")
sets$extra <- list(sym = toupper(EXTRA_ISGS), ens = character(0))
sets <- Filter(Negate(is.null), sets)

flag <- function(s) d$.sym %in% s$sym | d$.key %in% s$ens
hits <- sapply(sets, flag)
if (is.null(dim(hits))) hits <- matrix(hits, ncol = length(sets), dimnames = list(NULL, names(sets)))
d$canonical_ISG_source <- apply(hits, 1, function(r) paste(colnames(hits)[r], collapse = "+"))
is_isg <- d$canonical_ISG_source != ""

removed   <- d[is_isg, ]  |> select(-.key, -.sym) |> relocate(canonical_ISG_source, .after = gene_name) |> arrange(desc(log2FoldChange))
remaining <- d[!is_isg, ] |> select(-.key, -.sym, -canonical_ISG_source) |> arrange(desc(log2FoldChange))
write_csv(remaining, file.path(OUT, "DE_v2_IFNB_vs_Mock_protein_coding_without_canonical_type_I_ISGs.csv"))
write_csv(removed,   file.path(OUT, "DE_v2_IFNB_vs_Mock_protein_coding_canonical_type_I_ISGs_removed.csv"))

up <- function(x) sum(!is.na(x$padj) & x$padj < PADJ_CUT & x$log2FoldChange >  LFC_CUT)
dn <- function(x) sum(!is.na(x$padj) & x$padj < PADJ_CUT & x$log2FoldChange < -LFC_CUT)
message(sprintf("\nRemoved as canonical type I ISGs: %d genes (%d of them significantly up, %d down)", nrow(removed), up(removed), dn(removed)))
message(sprintf("Remaining: %d genes (significant: %d up, %d down at padj<%g & |log2FC|>%g)",
                nrow(remaining), up(remaining), dn(remaining), PADJ_CUT, LFC_CUT))
message("Top 10 remaining upregulated: ", paste(head(remaining$gene_name[!is.na(remaining$padj) & remaining$padj < PADJ_CUT], 10), collapse = ", "))
message("\nWrote:\n  ", file.path(OUT, "DE_v2_IFNB_vs_Mock_protein_coding_without_canonical_type_I_ISGs.csv"),
        "\n  ", file.path(OUT, "DE_v2_IFNB_vs_Mock_protein_coding_canonical_type_I_ISGs_removed.csv"))
