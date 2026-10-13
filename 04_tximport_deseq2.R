#!/usr/bin/env Rscript
# Import the re-quantified data, separate construct-derived transcripts from
# host genes, and build a DESeq2 object normalised on host genes only.
#
# Run from the repo directory with the project path available, e.g.
#   Rscript 04_tximport_deseq2.R
# or set PROJ first:
#   PROJ=/path/to/Maury_25161 Rscript 04_tximport_deseq2.R

suppressPackageStartupMessages({
  library(tximport); library(DESeq2); library(readr); library(dplyr)
  library(tibble); library(tidyr)
})

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
stopifnot(dir.exists(PROJ))

QUANT_DIR <- file.path(PROJ, "salmon_quant_v2")
TX2GENE   <- file.path(PROJ, "tx2gene_GRCm39.109_plus_model.tsv")
VIRAL_FA  <- Sys.getenv("VIRAL_FA", unset = file.path(PROJ, "ref", "recombinant_model.fa"))
OUT       <- file.path(PROJ, "requant_v2_outputs")
dir.create(OUT, showWarnings = FALSE)

# ---- samples ---------------------------------------------------------------
# sample_info.csv is tab-separated despite the extension; detect either way.
si_path <- file.path(PROJ, "sample_info.csv")
first   <- readLines(si_path, n = 1)
si <- if (grepl("\t", first)) read_tsv(si_path, show_col_types = FALSE) else
                              read_csv(si_path, show_col_types = FALSE)
si <- si |>
  mutate(
    sample_id  = as.character(sample_id),
    ifn        = ifelse(condition1 %in% c("None", "none", NA), "Mock", condition1),
    virus      = ifelse(condition2 %in% c("None", "none", NA), "Mock", condition2),
    group      = factor(paste(ifn, virus, sep = "_")),
    ifn        = relevel(factor(ifn), ref = "Mock"),
    virus      = relevel(factor(virus), ref = "Mock"),
    replicate  = factor(replicate)
  ) |>
  mutate(group = relevel(group, ref = "Mock_Mock"))

files <- file.path(QUANT_DIR, si$sample_id, "quant.sf")
names(files) <- si$sample_id
missing <- files[!file.exists(files)]
if (length(missing)) stop("Missing quant.sf for: ", paste(names(missing), collapse = ", "))

# ---- which targets are construct-derived -----------------------------------
viral_ids <- character(0)
if (file.exists(VIRAL_FA)) {
  hdr <- grep("^>", readLines(VIRAL_FA), value = TRUE)
  viral_ids <- sub("^>", "", sub("[[:space:]].*$", "", hdr))
  message("Construct features in reference: ", paste(viral_ids, collapse = ", "))
} else {
  warning("No construct FASTA at ", VIRAL_FA, " - proceeding host-only.")
}

# ---- import ----------------------------------------------------------------
tx2gene <- read_tsv(TX2GENE, show_col_types = FALSE)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

# ---- the construct deliverable --------------------------------------------
# Transcript-level numbers for the construct, straight from quant.sf, so the
# per-feature detail is not lost to gene-level summarisation.
if (length(viral_ids)) {
  per_sample <- lapply(si$sample_id, function(s) {
    q <- read_tsv(file.path(QUANT_DIR, s, "quant.sf"), show_col_types = FALSE)
    q |> filter(Name %in% viral_ids) |>
      transmute(feature = Name, length = Length, eff_length = EffectiveLength,
                TPM, reads = NumReads, sample_id = s)
  }) |> bind_rows()

  construct_tbl <- per_sample |>
    left_join(si |> select(sample_id, ifn, virus, group, replicate), by = "sample_id") |>
    arrange(feature, group, replicate)

  write_csv(construct_tbl, file.path(OUT, "construct_features_per_sample.csv"))

  # Wide reads matrix + a per-group summary, which is what you actually read.
  wide <- construct_tbl |>
    select(feature, sample_id, reads) |>
    pivot_wider(names_from = sample_id, values_from = reads)
  write_csv(wide, file.path(OUT, "construct_features_reads_matrix.csv"))

  summary_tbl <- construct_tbl |>
    group_by(feature, group) |>
    summarise(mean_reads = mean(reads), mean_TPM = mean(TPM),
              n_detected = sum(reads > 0), n = n(), .groups = "drop") |>
    arrange(feature, desc(mean_TPM))
  write_csv(summary_tbl, file.path(OUT, "construct_features_by_group.csv"))

  message("\n--- Construct detection by group (mean TPM) ---")
  print(as.data.frame(summary_tbl), row.names = FALSE)
  message("\nMock_Mock rows should be ~0. If they are not, the construct FASTA ",
          "shares sequence with the mouse transcriptome - see README.")
}

# ---- DESeq2, normalised on host genes only ---------------------------------
# Construct RNA can be a large and wildly uneven share of the library in
# infected samples. If it is allowed into size-factor estimation it drags the
# normalisation of those samples and manufactures apparent host-gene changes.
# So: keep the construct in the object (so it can be inspected and tested) but
# estimate size factors from host genes alone via controlGenes.
cd <- as.data.frame(si)
rownames(cd) <- cd$sample_id
cd <- cd[colnames(txi$counts), , drop = FALSE]   # keep colData aligned to txi
stopifnot(identical(rownames(cd), colnames(txi$counts)))

dds <- DESeqDataSetFromTximport(txi, colData = cd, design = ~ group)

host_idx <- which(!rownames(dds) %in% viral_ids)
stopifnot(length(host_idx) > 1000)
dds <- estimateSizeFactors(dds, controlGenes = host_idx)

# Light pre-filter: drop genes with no real signal anywhere.
keep <- rowSums(counts(dds, normalized = TRUE) >= 10) >= 3
message("\nGenes kept after filtering: ", sum(keep), " of ", length(keep))
dds <- dds[keep | rownames(dds) %in% viral_ids, ]

dds <- DESeq(dds)
saveRDS(dds, file.path(OUT, "dds_v2_decoyAware_hostNormalised.rds"))

vsd <- vst(dds, blind = FALSE)
saveRDS(vsd, file.path(OUT, "vsd_v2.rds"))

# ---- contrasts -------------------------------------------------------------
# The interferon comparisons this dataset is for, plus the virus arms.
contrasts <- list(
  IFNB_vs_Mock      = c("group", "IFNB_Mock",  "Mock_Mock"),
  IFNY_vs_Mock      = c("group", "IFNY_Mock",  "Mock_Mock"),
  EBOV_vs_Mock      = c("group", "Mock_EBOV",  "Mock_Mock"),
  MARV_vs_Mock      = c("group", "Mock_MARV",  "Mock_Mock"),
  EVLP_vs_Mock      = c("group", "Mock_EVLP",  "Mock_Mock"),
  MVLP_vs_Mock      = c("group", "Mock_MVLP",  "Mock_Mock"),
  IFNB_EBOV_vs_Mock = c("group", "IFNB_EBOV",  "Mock_Mock"),
  IFNB_MARV_vs_Mock = c("group", "IFNB_MARV",  "Mock_Mock"),
  IFNY_EBOV_vs_Mock = c("group", "IFNY_EBOV",  "Mock_Mock"),
  IFNY_MARV_vs_Mock = c("group", "IFNY_MARV",  "Mock_Mock")
)

avail <- levels(dds$group)
for (nm in names(contrasts)) {
  ct <- contrasts[[nm]]
  if (!all(ct[2:3] %in% avail)) { message("skip ", nm, " (group absent)"); next }
  res <- if (requireNamespace("ashr", quietly = TRUE)) {
    lfcShrink(dds, contrast = ct, type = "ashr")
  } else {
    lfcShrink(dds, contrast = ct, type = "normal")
  }
  out <- as.data.frame(res) |> rownames_to_column("gene_id") |>
    arrange(padj, desc(abs(log2FoldChange)))
  write_csv(out, file.path(OUT, paste0("DE_v2_", nm, ".csv")))
  message(sprintf("%-20s padj<0.05 & |LFC|>1: %d", nm,
                  sum(out$padj < 0.05 & abs(out$log2FoldChange) > 1, na.rm = TRUE)))
}

message("\nOutputs in: ", OUT)
