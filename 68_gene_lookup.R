#!/usr/bin/env Rscript
# Find out why a gene has no counts: wrong symbol? removed by the pre-filter? or really not expressed?   (default gene: Ch25h, the 25-hydroxycholesterol / 25HC enzyme)
# For each requested gene it reports
#   1. SYMBOL    the matching genes in the Ensembl GTF used for the reference (exact symbol; if none, similar symbols), with gene_id and biotype
#   2. PIPELINE  whether the gene is still in the DESeq2 object / DE tables. 04_tximport_deseq2.R drops every gene that does not have normalised counts >= 10 in at
#                least 3 samples, and those genes are then missing from ALL DE_v2_* tables, heatmaps and bar graphs - this is the usual reason for "no counts"
#   3. SALMON    the unfiltered Salmon quantification (quant.sf of every sample): NumReads and TPM per sample and group mean, summed over the gene's transcripts
# so you can see whether the gene is detected at all, in which groups, and whether the filter would have removed it.
# Usage:   GENE=Ch25h Rscript 68_gene_lookup.R        (comma-separate several: GENE=Ch25h,Cyp7b1)       PROJ from 'source setup_session.sh'
# Writes requant_v2_outputs/gene_lookup_<GENE>.csv (per sample) and prints the tables.

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr) })
select <- dplyr::select; filter <- dplyr::filter

GENES <- trimws(strsplit(Sys.getenv("GENE", unset = "Ch25h"), ",")[[1]])
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
strip <- function(x) sub("\\..*$", "", x)
options(width = 200, tibble.width = Inf)

# ---- 1. symbol lookup in the GTF ----
GTF <- file.path(PROJ, "ref", "Mus_musculus.GRCm39.109.gtf")
ln <- readLines(GTF); ln <- ln[!startsWith(ln, "#")]
gl <- ln[grepl("\tgene\t", ln, fixed = TRUE)]; rm(ln)
grab <- function(k) ifelse(grepl(paste0(k, ' "'), gl, fixed = TRUE), sub(paste0('.*', k, ' "([^"]+)".*'), "\\1", gl), NA_character_)
gmap <- tibble(gene_id = grab("gene_id"), gene_name = grab("gene_name"), biotype = grab("gene_biotype"), chr = sub("\t.*$", "", gl))
message(sprintf("GTF: %d genes (%d with a symbol)", nrow(gmap), sum(!is.na(gmap$gene_name))))

# ---- sample sheet, filtered objects ----
si_path <- file.path(PROJ, "sample_info.csv"); first <- readLines(si_path, n = 1)
si <- if (grepl("\t", first)) read_tsv(si_path, show_col_types = FALSE) else read_csv(si_path, show_col_types = FALSE)
si <- si |> mutate(sample_id = as.character(sample_id), ifn = ifelse(condition1 %in% c("None", "none", NA), "Mock", condition1),
                   virus = ifelse(condition2 %in% c("None", "none", NA), "Mock", condition2), group = paste(ifn, virus, sep = "_"))
dds_ids <- NULL
dds_f <- file.path(RNA, "dds_v2_decoyAware_hostNormalised.rds")
if (file.exists(dds_f)) { suppressPackageStartupMessages(library(DESeq2)); dds <- readRDS(dds_f); dds_ids <- strip(rownames(dds)) }
de_f <- file.path(RNA, "DE_v2_IFNB_vs_Mock_annotated.csv")
de <- if (file.exists(de_f)) read_csv(de_f, show_col_types = FALSE) |> mutate(gid = strip(gene_id)) else NULL
tx2gene <- read_tsv(file.path(PROJ, "tx2gene_GRCm39.109_plus_model.tsv"), show_col_types = FALSE); names(tx2gene)[1:2] <- c("tx", "gene_id")
tx2gene <- tx2gene |> mutate(tx = strip(tx), gene_id = strip(gene_id))
qfiles <- file.path(PROJ, "salmon_quant_v2", si$sample_id, "quant.sf"); names(qfiles) <- si$sample_id

for (G in GENES) {
  message("\n==================== ", G, " ====================")
  q <- tolower(G)
  hit <- gmap |> filter(tolower(gene_name) == q | tolower(strip(gene_id)) == q)
  if (!nrow(hit)) {
    sim <- gmap |> filter(grepl(q, tolower(gene_name), fixed = TRUE) | startsWith(tolower(gene_name), substr(q, 1, min(3, nchar(q))))) |> arrange(gene_name)
    message("1. SYMBOL: no gene with the exact symbol '", G, "' in the GTF. Similar symbols (", nrow(sim), "):")
    print(as.data.frame(head(sim, 25)), row.names = FALSE)
    message("   -> rerun with one of these (or its ENSMUSG id): GENE=<symbol> Rscript 68_gene_lookup.R")
    next
  }
  message("1. SYMBOL: found in the GTF:"); print(as.data.frame(hit), row.names = FALSE)
  for (k in seq_len(nrow(hit))) {
    gid <- strip(hit$gene_id[k]); nm <- hit$gene_name[k]
    message(sprintf("\n--- %s (%s, %s) ---", nm, hit$gene_id[k], hit$biotype[k]))
    # 2. pipeline
    in_dds <- !is.null(dds_ids) && gid %in% dds_ids
    message("2. PIPELINE: in the DESeq2 object / DE tables: ", if (in_dds) "YES" else "NO - removed by the pre-filter in 04 (normalised counts >= 10 in >= 3 samples) or absent from the reference")
    if (!is.null(de) && gid %in% de$gid) print(as.data.frame(de |> filter(gid == !!gid) |> select(gene_name, baseMean, log2FoldChange, lfcSE, padj)), row.names = FALSE)
    # 3. salmon
    txs <- tx2gene$tx[tx2gene$gene_id == gid]
    if (!length(txs)) { message("3. SALMON: no transcripts for this gene in tx2gene"); next }
    Q <- bind_rows(lapply(names(qfiles), function(s) {
      qq <- read_tsv(qfiles[[s]], col_types = cols_only(Name = col_character(), NumReads = col_double(), TPM = col_double()), progress = FALSE) |> mutate(tx = strip(Name)) |> filter(tx %in% txs)
      tibble(sample_id = s, n_transcripts = nrow(qq), NumReads = sum(qq$NumReads), TPM = sum(qq$TPM))
    })) |> left_join(si |> select(sample_id, group), by = "sample_id") |> mutate(gene = nm, gene_id = hit$gene_id[k])
    message(sprintf("3. SALMON (unfiltered quant.sf, %d transcript(s)): reads per sample", length(txs)))
    print(as.data.frame(Q |> select(sample_id, group, NumReads, TPM) |> mutate(NumReads = round(NumReads, 1), TPM = round(TPM, 2))), row.names = FALSE)
    message("   group means:")
    print(as.data.frame(Q |> group_by(group) |> summarise(n = n(), mean_NumReads = round(mean(NumReads), 1), mean_TPM = round(mean(TPM), 2), .groups = "drop")), row.names = FALSE)
    nOK <- sum(Q$NumReads >= 10)
    message(sprintf("   samples with >= 10 reads: %d of %d (the 04 pre-filter keeps a gene only if >= 3 samples have >= 10 normalised counts) -> %s", nOK, nrow(Q),
                    if (nOK >= 3) "passes the filter" else "WOULD BE REMOVED by the filter (too few reads)"))
    write_csv(Q, file.path(RNA, paste0("gene_lookup_", nm, ".csv")))
  }
}
