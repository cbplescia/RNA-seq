#!/usr/bin/env Rscript
# Subset the IFNB vs Mock results (protein-coding genes) to genes involved in
# lipid metabolism, using Gene Ontology (GO) Biological Process annotations.
#
# How the list is compiled
#   - GO term GO:0006629 "lipid metabolic process" is the umbrella term. Using
#     org.Mm.eg.db's GOALL mapping, a gene counts if it is annotated to that term
#     OR to any more specific child term (e.g. fatty acid biosynthesis, cholesterol
#     catabolism, phospholipid metabolism), because GOALL includes all ancestors.
#   - That umbrella already covers both synthesis (GO:0008610 "lipid biosynthetic
#     process") and breakdown (GO:0016042 "lipid catabolic process"); both get their
#     own TRUE/FALSE column so you can filter to either.
#   - Extra TRUE/FALSE columns flag common sub-areas (fatty acid, sterol/cholesterol,
#     steroid, phospholipid, sphingolipid, glycerolipid, triglyceride).
#   - Set EXCLUDE_IEA <- TRUE to drop annotations that are purely computational
#     (evidence code IEA) and keep only curated/experimental evidence.
#
# Input : DE_v2_IFNB_vs_Mock_protein_coding.csv   (from 07 + 08; shrunk log2FC)
# Output: DE_v2_IFNB_vs_Mock_protein_coding_lipid_metabolism.csv
#         lipid_metabolism_gene_GO_terms.csv      (which GO terms put each gene on the list)
# Needs the Bioconductor packages org.Mm.eg.db and GO.db:
#   BiocManager::install(c("org.Mm.eg.db", "GO.db"))   or  conda install -c bioconda bioconductor-org.mm.eg.db bioconductor-go.db

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr) })
for (pkg in c("org.Mm.eg.db", "GO.db", "AnnotationDbi"))
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Package '", pkg, "' is needed. Install with:  BiocManager::install(c('org.Mm.eg.db','GO.db'))")

EXCLUDE_IEA <- FALSE

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs")

de <- read_csv(file.path(OUT, "DE_v2_IFNB_vs_Mock_protein_coding.csv"), show_col_types = FALSE)
de$.key <- sub("\\.[0-9]+$", "", de$gene_id)

# sub-area flags: column name -> GO ID
flags <- c(
  lipid_metabolic     = "GO:0006629",   # umbrella term
  lipid_biosynthetic  = "GO:0008610",
  lipid_catabolic     = "GO:0016042",
  fatty_acid          = "GO:0006631",
  sterol_cholesterol  = "GO:0016125",   # sterol metabolic process
  steroid             = "GO:0008202",
  phospholipid        = "GO:0006644",
  sphingolipid        = "GO:0006665",
  glycerolipid        = "GO:0046486",
  triglyceride        = "GO:0006641")

go <- suppressMessages(AnnotationDbi::select(
  org.Mm.eg.db::org.Mm.eg.db, keys = unique(de$.key), keytype = "ENSEMBL",
  columns = c("GOALL", "EVIDENCEALL", "ONTOLOGYALL"))) |>
  filter(ONTOLOGYALL == "BP", !is.na(GOALL))
if (EXCLUDE_IEA) go <- go |> filter(EVIDENCEALL != "IEA")

hits <- go |> filter(GOALL %in% flags) |> distinct(ENSEMBL, GOALL)
flag_tab <- lapply(names(flags), function(n)
  tibble(.key = unique(hits$ENSEMBL[hits$GOALL == flags[[n]]]), !!n := TRUE)) |>
  Reduce(function(a, b) full_join(a, b, by = ".key"), x = _)

lip <- de |> inner_join(flag_tab, by = ".key") |>
  mutate(across(all_of(names(flags)), ~ tidyr::replace_na(.x, FALSE))) |>
  select(-.key) |> arrange(padj, desc(abs(log2FoldChange)))
write_csv(lip, file.path(OUT, "DE_v2_IFNB_vs_Mock_protein_coding_lipid_metabolism.csv"))

# which specific GO terms (children included) put each gene on the list
all_lip <- go |> filter(GOALL %in% flags) |> distinct(ENSEMBL) 
child_terms <- go |> filter(ENSEMBL %in% all_lip$ENSEMBL) |> distinct(ENSEMBL, GOALL)
# Keep only terms that sit under the umbrella term (so unrelated GO terms are not listed)
umbrella_desc <- c("GO:0006629", unlist(as.list(GO.db::GOBPOFFSPRING)["GO:0006629"]))
child_terms <- child_terms |> filter(GOALL %in% umbrella_desc)
names_tab <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = unique(child_terms$GOALL),
                                                    keytype = "GOID", columns = "TERM"))
terms_out <- child_terms |> left_join(names_tab, by = c("GOALL" = "GOID")) |>
  left_join(de |> select(.key, gene_id, gene_name), by = c("ENSEMBL" = ".key")) |>
  transmute(gene_id, gene_name, go_id = GOALL, go_term = TERM) |>
  arrange(gene_name, go_id)
write_csv(terms_out, file.path(OUT, "lipid_metabolism_gene_GO_terms.csv"))

nsig <- function(d) sum(!is.na(d$padj) & d$padj < 0.05 & abs(d$log2FoldChange) > 1)
message("Protein-coding genes in IFNB vs Mock table: ", nrow(de))
message("Annotated to lipid metabolic process (incl. all child terms): ", nrow(lip),
        if (EXCLUDE_IEA) "  [IEA annotations excluded]" else "")
message("  biosynthetic: ", sum(lip$lipid_biosynthetic), "   catabolic: ", sum(lip$lipid_catabolic))
message("  of these, significant (padj<0.05 & |log2FC|>1): ", nsig(lip),
        "  (up: ", sum(lip$padj < 0.05 & lip$log2FoldChange > 1, na.rm = TRUE),
        ", down: ", sum(lip$padj < 0.05 & lip$log2FoldChange < -1, na.rm = TRUE), ")")
message("\nWrote:\n  ", file.path(OUT, "DE_v2_IFNB_vs_Mock_protein_coding_lipid_metabolism.csv"),
        "\n  ", file.path(OUT, "lipid_metabolism_gene_GO_terms.csv"))
