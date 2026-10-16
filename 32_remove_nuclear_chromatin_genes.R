#!/usr/bin/env Rscript
# Further refine the IFNB vs Mock list (interferon / antiviral genes already removed by
# 30_remove_canonical_isgs.R) by dropping genes that are BOTH
#   (a) fully nuclear   : annotated to the nucleus (GO:0005634) and to NO non-nuclear compartment
#                         (cytoplasm GO:0005737, membrane GO:0016020, extracellular region GO:0005576;
#                         these include all their child terms, e.g. mitochondrion, ER, plasma membrane), and
#   (b) chromatin/DNA-only in function : after ignoring generic terms (protein binding, metal/zinc ion
#                         binding, ATP binding, ...), EVERY direct GO molecular-function annotation falls
#                         under DNA binding, DNA-binding transcription factor / transcription regulator
#                         activity, chromatin / histone / nucleosome binding, histone-modifying activity
#                         or catalytic activity acting on DNA (polymerases, helicases, methyltransferases,
#                         nucleases, ...).
# A nuclear gene that also binds RNA, phosphorylates proteins, is cytoplasmic in any annotation, or has no
# molecular-function annotation at all is KEPT. All GO evidence (including inferred, IEA) is used, so any
# cytoplasmic or non-DNA annotation keeps a gene in the list - the removal is deliberately conservative.
#
# Input : DE_v2_IFNB_vs_Mock_protein_coding_without_ISG_and_antiviral_genes.csv   (from 30)
# Output: DE_v2_IFNB_vs_Mock_protein_coding_without_ISG_antiviral_nuclear_genes.csv  (the refined list)
#         DE_v2_IFNB_vs_Mock_protein_coding_nuclear_chromatin_genes_removed.csv      (with the reasons)
#         DE_v2_IFNB_vs_Mock_refined_up_down.csv    (significant up/down; for 28: non_ifn_nuclear_filtered)
# Edit CORE_MF / GENERIC_MF below to change what counts.

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr) })
for (pkg in c("org.Mm.eg.db", "GO.db", "AnnotationDbi"))
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install ", pkg, " (BiocManager::install)")
OrgDb <- org.Mm.eg.db::org.Mm.eg.db
select <- dplyr::select; first <- dplyr::first; rename <- dplyr::rename
slice <- dplyr::slice; filter <- dplyr::filter; desc <- dplyr::desc; count <- dplyr::count

LFC_CUT <- 1; PADJ_CUT <- 0.05

NUCLEUS    <- "GO:0005634"
NONNUCLEAR <- c("GO:0005737", "GO:0016020", "GO:0005576")  # cytoplasm, membrane, extracellular region
CORE_MF <- c("GO:0003677",  # DNA binding
             "GO:0003700",  # DNA-binding transcription factor activity
             "GO:0140110",  # transcription regulator activity
             "GO:0003682",  # chromatin binding
             "GO:0042393",  # histone binding
             "GO:0140993",  # histone modifying activity
             "GO:0031491",  # nucleosome binding
             "GO:0140097",  # catalytic activity, acting on DNA
             "GO:0000976",  # transcription cis-regulatory region binding
             "GO:0043565")  # sequence-specific DNA binding
GENERIC_MF <- c("GO:0005488", "GO:0005515", "GO:0042802", "GO:0042803", "GO:0046982", "GO:0044877",  # (protein) binding
                "GO:0046872", "GO:0043169", "GO:0043167", "GO:0046914", "GO:0008270",                # metal / zinc ion binding
                "GO:0003676", "GO:0097159", "GO:1901363",                                           # nucleic acid / compound binding
                "GO:0005524", "GO:0000166", "GO:1901265", "GO:0036094", "GO:0043168")                # ATP / nucleotide binding

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs")

d <- read_csv(file.path(OUT, "DE_v2_IFNB_vs_Mock_protein_coding_without_ISG_and_antiviral_genes.csv"), show_col_types = FALSE)
d$.key <- sub("\\.[0-9]+$", "", d$gene_id)
message("Genes in the starting list: ", nrow(d))

# ---- (a) fully nuclear ------------------------------------------------------------------------------
cc <- suppressMessages(AnnotationDbi::select(OrgDb, keys = unique(d$.key), keytype = "ENSEMBL",
                                             columns = c("GOALL", "ONTOLOGYALL"))) |> filter(ONTOLOGYALL == "CC")
in_nuc  <- unique(cc$ENSEMBL[cc$GOALL == NUCLEUS])
in_non  <- unique(cc$ENSEMBL[cc$GOALL %in% NONNUCLEAR])
nuc_only <- setdiff(in_nuc, in_non)

# ---- (b) chromatin/DNA-only function ------------------------------------------------------------------
offspring <- AnnotationDbi::as.list(GO.db::GOMFOFFSPRING)
allowed <- unique(c(CORE_MF, na.omit(unlist(offspring[intersect(CORE_MF, names(offspring))]))))
mf <- suppressMessages(AnnotationDbi::select(OrgDb, keys = unique(d$.key), keytype = "ENSEMBL",
                                             columns = c("GO", "ONTOLOGY"))) |>
  filter(ONTOLOGY == "MF", !is.na(GO), !GO %in% GENERIC_MF) |> distinct(ENSEMBL, GO)
mf_sum <- mf |> group_by(ENSEMBL) |> summarise(n = n(), n_ok = sum(GO %in% allowed), .groups = "drop")
chrom_only <- mf_sum$ENSEMBL[mf_sum$n > 0 & mf_sum$n == mf_sum$n_ok]

names_tab <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = unique(mf$GO), keytype = "GOID", columns = "TERM"))
mf_terms <- mf |> left_join(names_tab, by = c("GO" = "GOID")) |> group_by(ENSEMBL) |>
  summarise(direct_MF_terms = paste(sort(unique(TERM)), collapse = "; "), .groups = "drop")

d <- d |> mutate(fully_nuclear = .key %in% nuc_only, chromatin_DNA_only_function = .key %in% chrom_only) |>
  left_join(mf_terms, by = c(".key" = "ENSEMBL"))
drop <- d$fully_nuclear & d$chromatin_DNA_only_function

removed   <- d[drop, ]  |> select(-.key) |> relocate(fully_nuclear, chromatin_DNA_only_function, direct_MF_terms, .after = gene_name) |>
  arrange(desc(log2FoldChange))
remaining <- d[!drop, ] |> select(-.key, -fully_nuclear, -chromatin_DNA_only_function, -direct_MF_terms) |> arrange(desc(log2FoldChange))
write_csv(remaining, file.path(OUT, "DE_v2_IFNB_vs_Mock_protein_coding_without_ISG_antiviral_nuclear_genes.csv"))
write_csv(removed,   file.path(OUT, "DE_v2_IFNB_vs_Mock_protein_coding_nuclear_chromatin_genes_removed.csv"))
resp <- remaining |> filter(!is.na(padj), padj < PADJ_CUT, abs(log2FoldChange) > LFC_CUT) |>
  mutate(direction = ifelse(log2FoldChange > 0, "up", "down")) |> relocate(direction, .after = padj)
write_csv(resp, file.path(OUT, "DE_v2_IFNB_vs_Mock_refined_up_down.csv"))

up <- function(x) sum(!is.na(x$padj) & x$padj < PADJ_CUT & x$log2FoldChange >  LFC_CUT)
dn <- function(x) sum(!is.na(x$padj) & x$padj < PADJ_CUT & x$log2FoldChange < -LFC_CUT)
message(sprintf("\nFully nuclear (nucleus, no non-nuclear annotation): %d genes", sum(d$fully_nuclear)))
message(sprintf("Chromatin/DNA-only molecular function:               %d genes", sum(d$chromatin_DNA_only_function)))
message(sprintf("BOTH (removed):                                       %d genes (%d significantly up, %d down)", nrow(removed), up(removed), dn(removed)))
message(sprintf("Remaining: %d genes; significant: %d up, %d down (padj<%g & |log2FC|>%g)",
                nrow(remaining), up(remaining), dn(remaining), PADJ_CUT, LFC_CUT))
sig_rm <- removed |> filter(!is.na(padj), padj < PADJ_CUT, abs(log2FoldChange) > LFC_CUT)
if (nrow(sig_rm)) message("Significant genes removed: ", paste(head(sig_rm$gene_name, 30), collapse = ", "),
                          if (nrow(sig_rm) > 30) paste0(" ... (+", nrow(sig_rm) - 30, " more)") else "")
message("\nNext:  Rscript 28_lipid_function_analysis.R non_ifn_nuclear_filtered")
