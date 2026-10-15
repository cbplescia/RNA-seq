#!/usr/bin/env Rscript
# Write annotated copies of the differential expression CSVs with a gene_name
# (symbol) column next to gene_id. Symbols come from the same Ensembl GRCm39
# release 109 GTF used to build the reference, so they match the gene IDs.
# Originals are untouched; copies are written as <name>_annotated.csv.

suppressPackageStartupMessages({ library(readr); library(dplyr) })

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs")
GTF <- file.path(PROJ, "ref", "Mus_musculus.GRCm39.109.gtf")
stopifnot(file.exists(GTF))

# ---- gene_id -> gene_name from the GTF "gene" lines --------------------------
ln <- readLines(GTF)
ln <- ln[!startsWith(ln, "#")]
ln <- ln[grepl("\tgene\t", ln, fixed = TRUE)]
attr <- sub("^([^\t]*\t){8}", "", ln)
grab <- function(key) {
  m <- regmatches(attr, regexpr(paste0(key, ' "[^"]*"'), attr))
  out <- rep(NA_character_, length(attr))
  out[grepl(paste0(key, ' "'), attr)] <- sub(paste0(key, ' "([^"]*)"'), "\\1", m)
  out
}
map <- tibble(gene_id = grab("gene_id"), gene_name = grab("gene_name"),
              gene_biotype = grab("gene_biotype")) |> distinct(gene_id, .keep_all = TRUE)
message("Genes in GTF: ", nrow(map), "  (with a symbol: ", sum(!is.na(map$gene_name)), ")")

strip_ver <- function(x) sub("\\.[0-9]+$", "", x)

# ---- annotate every DE table --------------------------------------------------
files <- list.files(OUT, pattern = "^(DE_v2_.*|changed_calls_.*)\\.csv$", full.names = TRUE)
files <- files[!grepl("_annotated\\.csv$", files)]
if (!length(files)) stop("No DE csv files found in ", OUT)

for (f in files) {
  d <- read_csv(f, show_col_types = FALSE)
  d$gene_id <- as.character(d$gene_id)
  a <- d |>
    mutate(.key = strip_ver(gene_id)) |>
    left_join(map |> mutate(.key = strip_ver(gene_id)) |> select(.key, gene_name, gene_biotype),
              by = ".key") |>
    # construct features and genes with no symbol fall back to their ID
    mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name)) |>
    select(gene_id, gene_name, gene_biotype, everything(), -.key)
  outf <- sub("\\.csv$", "_annotated.csv", f)
  write_csv(a, outf)
  message(sprintf("%-55s %6d rows, %d without a symbol",
                  basename(outf), nrow(a), sum(a$gene_name == a$gene_id)))
}
message("\nWrote annotated copies to: ", OUT)
