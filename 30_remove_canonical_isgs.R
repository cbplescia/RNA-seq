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
#         DE_v2_IFNB_vs_Mock_non_canonical_up_down.csv   (significant up/down genes among the remaining; for 28)
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

# responders among the non-canonical genes, for 28_lipid_function_analysis.R non_canonical_isgs
resp <- remaining |> filter(!is.na(padj), padj < PADJ_CUT, abs(log2FoldChange) > LFC_CUT) |>
  mutate(direction = ifelse(log2FoldChange > 0, "up", "down")) |>
  relocate(direction, .after = padj)
write_csv(resp, file.path(OUT, "DE_v2_IFNB_vs_Mock_non_canonical_up_down.csv"))

up <- function(x) sum(!is.na(x$padj) & x$padj < PADJ_CUT & x$log2FoldChange >  LFC_CUT)
dn <- function(x) sum(!is.na(x$padj) & x$padj < PADJ_CUT & x$log2FoldChange < -LFC_CUT)
message(sprintf("\nRemoved as canonical type I ISGs: %d genes (%d of them significantly up, %d down)", nrow(removed), up(removed), dn(removed)))
message(sprintf("Remaining: %d genes (significant: %d up, %d down at padj<%g & |log2FC|>%g)",
                nrow(remaining), up(remaining), dn(remaining), PADJ_CUT, LFC_CUT))
message("Top 10 remaining upregulated: ", paste(head(remaining$gene_name[!is.na(remaining$padj) & remaining$padj < PADJ_CUT], 10), collapse = ", "))
message("\nNext:  Rscript 28_lipid_function_analysis.R non_canonical_isgs")
message("\nWrote:\n  ", file.path(OUT, "DE_v2_IFNB_vs_Mock_protein_coding_without_canonical_type_I_ISGs.csv"),
        "\n  ", file.path(OUT, "DE_v2_IFNB_vs_Mock_protein_coding_canonical_type_I_ISGs_removed.csv"))

# =====================================================================================================
# Second, stricter removal: ALSO drop every gene annotated in GO (Biological Process, including all
# child terms) to an interferon-response or antiviral term, so the remaining list is free of the
# "response to IFN / defense against virus" genes and any pathway analysis on it shows what changes
# outside those programmes. Edit GO_REMOVE_TERMS to change what counts. Produces a second set of files:
#   DE_v2_IFNB_vs_Mock_protein_coding_without_ISG_and_antiviral_genes.csv
#   DE_v2_IFNB_vs_Mock_protein_coding_ISG_and_antiviral_genes_removed.csv   (with removed_by)
#   DE_v2_IFNB_vs_Mock_non_ISG_antiviral_up_down.csv                         (for 28: non_ifn_antiviral)
# =====================================================================================================
BROAD_IMMUNE <- FALSE    # TRUE also removes everything under "innate immune response" and "defense response"
GO_REMOVE_TERMS <- c(
  "GO:0034340",  # response to type I interferon
  "GO:0071357",  # cellular response to type I interferon
  "GO:0060337",  # type I interferon-mediated signaling pathway
  "GO:0035455",  # response to interferon-alpha
  "GO:0035457",  # cellular response to interferon-alpha
  "GO:0035456",  # response to interferon-beta
  "GO:0035458",  # cellular response to interferon-beta
  "GO:0032606",  # type I interferon production
  "GO:0009615",  # response to virus
  "GO:0051607",  # defense response to virus
  "GO:0140374",  # antiviral innate immune response
  "GO:0050688",  # regulation of defense response to virus
  "GO:0050792",  # regulation of viral process
  "GO:0045069")  # regulation of viral genome replication
if (BROAD_IMMUNE) GO_REMOVE_TERMS <- c(GO_REMOVE_TERMS, "GO:0045087", "GO:0006952")  # innate immune response, defense response

for (pkg in c("org.Mm.eg.db", "GO.db", "AnnotationDbi"))
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install ", pkg, " (BiocManager::install) for the GO-based removal")
select <- dplyr::select; first <- dplyr::first; rename <- dplyr::rename
slice <- dplyr::slice; filter <- dplyr::filter; desc <- dplyr::desc; count <- dplyr::count

go_hit <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys = unique(d$.key), keytype = "ENSEMBL",
                                                 columns = c("GOALL", "ONTOLOGYALL"))) |>
  filter(ONTOLOGYALL == "BP", GOALL %in% GO_REMOVE_TERMS) |> distinct(ENSEMBL, GOALL)
nm <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = GO_REMOVE_TERMS, keytype = "GOID", columns = "TERM"))
go_hit <- go_hit |> left_join(nm, by = c("GOALL" = "GOID"))
by_gene <- go_hit |> group_by(ENSEMBL) |> summarise(GO_terms = paste(sort(unique(TERM)), collapse = "; "), .groups = "drop")
d2 <- d |> left_join(by_gene, by = c(".key" = "ENSEMBL"))
d2$GO_terms[is.na(d2$GO_terms)] <- ""
is_go <- d2$GO_terms != ""
d2$removed_by <- trimws(paste(ifelse(d2$canonical_ISG_source != "", paste0("curated set: ", d2$canonical_ISG_source), ""),
                              ifelse(is_go, paste0("GO: ", d2$GO_terms), ""), sep = " | "), whitespace = "[ |]")
drop2 <- is_isg | is_go

removed2   <- d2[drop2, ]  |> select(-.key, -.sym, -canonical_ISG_source, -GO_terms) |> relocate(removed_by, .after = gene_name) |> arrange(desc(log2FoldChange))
remaining2 <- d2[!drop2, ] |> select(-.key, -.sym, -canonical_ISG_source, -GO_terms, -removed_by) |> arrange(desc(log2FoldChange))
write_csv(remaining2, file.path(OUT, "DE_v2_IFNB_vs_Mock_protein_coding_without_ISG_and_antiviral_genes.csv"))
write_csv(removed2,   file.path(OUT, "DE_v2_IFNB_vs_Mock_protein_coding_ISG_and_antiviral_genes_removed.csv"))
resp2 <- remaining2 |> filter(!is.na(padj), padj < PADJ_CUT, abs(log2FoldChange) > LFC_CUT) |>
  mutate(direction = ifelse(log2FoldChange > 0, "up", "down")) |> relocate(direction, .after = padj)
write_csv(resp2, file.path(OUT, "DE_v2_IFNB_vs_Mock_non_ISG_antiviral_up_down.csv"))

message(sprintf("\n--- Stricter removal (curated ISG sets + GO interferon/antiviral terms%s) ---", if (BROAD_IMMUNE) " + innate immune/defense" else ""))
message(sprintf("Removed: %d genes in total (%d only via the GO terms; %d significantly up, %d down among all removed)",
                nrow(removed2), sum(is_go & !is_isg), up(removed2), dn(removed2)))
message(sprintf("Remaining: %d genes; significant: %d up, %d down (padj<%g & |log2FC|>%g)",
                nrow(remaining2), up(remaining2), dn(remaining2), PADJ_CUT, LFC_CUT))
message("Genes removed per GO term:")
print(as.data.frame(go_hit |> count(TERM, sort = TRUE)), row.names = FALSE)
message("\nNext:  Rscript 28_lipid_function_analysis.R non_ifn_antiviral")
