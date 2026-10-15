#!/usr/bin/env Rscript
# Split the shrunk IFNB/IFNY vs Mock DE tables (annotated by 07) into
# protein-coding genes and everything else. Run 07_add_gene_symbols.R first.
# Rows with no biotype (e.g. recombinant construct features) go in the
# non-protein-coding file.

suppressPackageStartupMessages({ library(readr); library(dplyr) })

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs")

for (nm in c("IFNB_vs_Mock", "IFNY_vs_Mock")) {
  f <- file.path(OUT, paste0("DE_v2_", nm, "_annotated.csv"))
  if (!file.exists(f)) stop("Missing ", f, " - run 07_add_gene_symbols.R first")
  d  <- read_csv(f, show_col_types = FALSE)
  pc <- !is.na(d$gene_biotype) & d$gene_biotype == "protein_coding"

  write_csv(d[pc, ],  file.path(OUT, paste0("DE_v2_", nm, "_protein_coding.csv")))
  write_csv(d[!pc, ], file.path(OUT, paste0("DE_v2_", nm, "_non_protein_coding.csv")))

  message("\n=== ", nm, " ===")
  message("protein coding: ", sum(pc), "   other: ", sum(!pc), "   total: ", nrow(d))
  print(as.data.frame(d[!pc, ] |> count(gene_biotype, sort = TRUE) |> head(8)), row.names = FALSE)
}
message("\nWrote files to: ", OUT)
