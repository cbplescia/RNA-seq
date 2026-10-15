#!/usr/bin/env Rscript
# Is the old-vs-new difference just lfcShrink? Re-computes the IFNB/IFNY vs Mock
# contrasts UNSHRUNK from the saved DESeq2 object, compares them to the original
# results, and prints size factors by group. Does not re-run the pipeline.

suppressPackageStartupMessages({
  library(DESeq2); library(readr); library(dplyr)
})

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs")

dds <- readRDS(file.path(OUT, "dds_v2_decoyAware_hostNormalised.rds"))

message("\n--- size factors by group (should be similar across groups) ---")
sf <- data.frame(group = dds$group, sizeFactor = sizeFactors(dds))
print(as.data.frame(sf |> group_by(group) |>
  summarise(n = n(), mean_sf = mean(sizeFactor), min_sf = min(sizeFactor),
            max_sf = max(sizeFactor))), row.names = FALSE)

read_old <- function(path) {
  d <- suppressWarnings(read_csv(path, show_col_types = FALSE))
  id <- intersect(c("gene_id", "gene", "...1", "GENEID", "row"), names(d))[1]
  d |> transmute(gene_id = .data[[id]], lfc_old = log2FoldChange, padj_old = padj)
}

pairs <- list(
  IFNB_vs_Mock = list(old = "DE_IFNB_vs_Mock_fresh.csv", ct = c("group", "IFNB_Mock", "Mock_Mock")),
  IFNY_vs_Mock = list(old = "DE_IFNY_vs_Mock_fresh.csv", ct = c("group", "IFNY_Mock", "Mock_Mock"))
)

rows <- list()
for (nm in names(pairs)) {
  res <- results(dds, contrast = pairs[[nm]]$ct)
  new <- as.data.frame(res) |> tibble::rownames_to_column("gene_id") |>
    transmute(gene_id, lfc_new = log2FoldChange, padj_new = padj)
  write_csv(new, file.path(OUT, paste0("DE_v2_unshrunk_", nm, ".csv")))

  m <- inner_join(read_old(file.path(PROJ, pairs[[nm]]$old)), new, by = "gene_id") |>
    mutate(sig_old = !is.na(padj_old) & padj_old < 0.05 & abs(lfc_old) > 1,
           sig_new = !is.na(padj_new) & padj_new < 0.05 & abs(lfc_new) > 1)
  rows[[nm]] <- data.frame(
    contrast = nm, n = nrow(m),
    pearson_r = round(cor(m$lfc_old, m$lfc_new, use = "complete.obs"), 4),
    kept = sum(m$sig_old & m$sig_new),
    lost = sum(m$sig_old & !m$sig_new),
    gained = sum(!m$sig_old & m$sig_new))
}
message("\n--- OLD vs NEW (unshrunk) ---")
print(bind_rows(rows), row.names = FALSE)
message("\nCompare with the shrunk run: r = 0.793 (IFNB), 0.846 (IFNY).")
