#!/usr/bin/env Rscript
# Why is a gene (not) in the lipid / endocytic-vesicle enzyme list?
# Traces one gene through every filter used by scripts 07, 08, 13 and 16 and says
# which step dropped it.   Usage:  Rscript 17_diagnose_gene.R Naglu
suppressPackageStartupMessages({ library(readr); library(dplyr) })
for (pkg in c("org.Mm.eg.db", "GO.db", "AnnotationDbi"))
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install ", pkg, " (BiocManager::install)")

sym <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(sym)) sym <- "Naglu"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs")
rd <- function(f) { p <- file.path(OUT, f); if (file.exists(p)) read_csv(p, show_col_types = FALSE) else NULL }
find <- function(d) if (is.null(d)) NULL else d |> filter(toupper(gene_name) == toupper(sym))
yn <- function(x) if (x) "YES" else "NO"

message("=== ", sym, " ===")
a <- find(rd("DE_v2_IFNB_vs_Mock_annotated.csv"))
message("1. In the IFNB vs Mock DE table (survived the expression filter in 04): ", yn(!is.null(a) && nrow(a) > 0))
if (is.null(a) || !nrow(a)) { message("   -> dropped at step 04 (too lowly expressed) or no gene with this symbol. Nothing downstream can include it."); quit(status = 0) }
message(sprintf("   gene_id %s | biotype %s | log2FC %.2f | padj %s", a$gene_id[1], a$gene_biotype[1], a$log2FoldChange[1], format(a$padj[1], digits = 3)))
message("2. protein_coding biotype: ", yn(identical(a$gene_biotype[1], "protein_coding")))
l <- find(rd("DE_v2_IFNB_vs_Mock_protein_coding_lipid_metabolism.csv"))
in_lip <- !is.null(l) && nrow(l) > 0
message("3. In the GO lipid-metabolism list (script 13): ", yn(in_lip))

key <- sub("\\.[0-9]+$", "", a$gene_id[1])
go <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys = key, keytype = "ENSEMBL",
        columns = c("GOALL", "ONTOLOGYALL", "EVIDENCEALL"))) |> filter(!is.na(GOALL))
if (!nrow(go)) { message("   org.Mm.eg.db has NO GO annotations for this Ensembl ID."); quit(status = 0) }
nm <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = unique(go$GOALL), keytype = "GOID", columns = "TERM"))
go <- left_join(go, nm, by = c("GOALL" = "GOID"))

lipid_under <- c("GO:0006629", unlist(as.list(GO.db::GOBPOFFSPRING)["GO:0006629"]))
bp_lip <- go |> filter(ONTOLOGYALL == "BP", GOALL %in% lipid_under) |> distinct(GOALL, TERM)
message("   GO lipid terms for this gene: ", if (nrow(bp_lip)) paste(bp_lip$TERM, collapse = "; ") else "none")
if (!in_lip && !nrow(bp_lip)) message("   -> NOT annotated to any lipid metabolic process term in GO, so script 13 never selects it.\n      (script 16 can still bring it in via EXTRA_GENES)")

cc <- go |> filter(ONTOLOGYALL == "CC") |> distinct(GOALL, TERM, EVIDENCEALL)
vesicle_ids <- c("GO:0030139", "GO:0005768", "GO:0005771", "GO:0045335")
comp <- cc |> filter(grepl("endo|lysos|vesicle|vacuol|phagos|membrane|lumen", TERM, ignore.case = TRUE))
message("4. Relevant cellular-component annotations (term | evidence codes):")
print(as.data.frame(comp |> group_by(GOALL, TERM) |> summarise(evidence = paste(sort(unique(EVIDENCEALL)), collapse = ","), .groups = "drop")), row.names = FALSE)
v <- any(cc$GOALL %in% vesicle_ids); ly <- any(cc$GOALL == "GO:0005764")
message("   endocytic vesicle / endosome / phagosome / MVB: ", yn(v), "   lysosome: ", yn(ly))
enz <- any(go$ONTOLOGYALL == "MF" & go$GOALL == "GO:0003824")
message("5. Enzyme (catalytic activity, GO:0003824): ", yn(enz))
message("\nVerdict:")
if (!in_lip) message(" - not in the lipid list (step 13)")
if (!v && !ly) message(" - no endosome/vesicle/lysosome annotation (step 16 compartment filter)")
if (!v && ly)  message(" - lysosome-only: dropped by old script 16 default (INCLUDE_LYSOSOME = FALSE); kept now")
if (!enz) message(" - no catalytic-activity annotation (step 16 enzyme filter)")
if (in_lip && (v || ly) && enz) message(" - passes every filter")
