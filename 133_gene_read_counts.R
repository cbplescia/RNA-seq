#!/usr/bin/env Rscript
# What are the actual read counts of a list of genes, sample by sample - including genes the low-count pre-filter removed?   Reads the UNFILTERED Salmon output, so every gene in the GTF can be looked up.
# Three units are reported for every gene and sample:
#   raw reads     NumReads summed over the gene's transcripts (what was actually sequenced; best for judging whether a gene is on at all)
#   TPM           transcripts per million (corrects for transcript length; best for comparing how much different genes are expressed within a sample)
#   normalized    raw reads / DESeq2 size factor (the scale of the pre-filter and of DESeq2; best for comparing the SAME gene between conditions)
# The pre-filter in 04 kept a gene only if it had >= 10 normalized reads in >= 3 of the 33 samples; genes that fail are not in the DE tables or the VST object, which is why they cannot be drawn on the VST plot. VST is a good scale for genes that passed (clustering, plots) but is not defined for removed genes.
# Outputs (requant_v2_outputs/figures_gene_counts/): gene_read_counts_long.tsv (gene x sample: reads, TPM, normalized), gene_read_counts_by_group.tsv (mean / min / max per group), reads_per_sample (raw reads, log scale), normalized_counts_vs_filter (with the pre-filter cut-off drawn), tpm_by_group (.png/.pdf).
# Env: OUT_NAME (suffix for the output file names), GENE (comma list of mouse symbols, case as in the GTF), GROUPS (default Mock_Mock,IFNB_Mock,IFNY_Mock; use GROUPS=all for every group), NCOL (panels per row, 6), TITLE, PROJ.   Usage: GENE=Alb,Ins2,Axl Rscript 133_gene_read_counts.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2); library(DESeq2) })
select <- dplyr::select; filter <- dplyr::filter; count <- dplyr::count; summarise <- dplyr::summarise; mutate <- dplyr::mutate; arrange <- dplyr::arrange   # attached Bioconductor / other packages can mask these verbs
GENES <- trimws(strsplit(Sys.getenv("GENE", "Alb,Ins2,Krt14,Snap25,Axl"), ",")[[1]]); GENES <- GENES[nzchar(GENES)]
ORD <- c("Mock_Mock", "IFNB_Mock", "IFNY_Mock", "Mock_EBOV", "IFNB_EBOV", "IFNY_EBOV", "Mock_MARV", "IFNB_MARV", "IFNY_MARV", "Mock_EVLP", "Mock_MVLP")
LAB <- c(Mock_Mock = "Mock", IFNB_Mock = "IFNB", IFNY_Mock = "IFNG", Mock_EBOV = "EBOV", IFNB_EBOV = "IFNB+EBOV", IFNY_EBOV = "IFNG+EBOV", Mock_MARV = "MARV", IFNB_MARV = "IFNB+MARV", IFNY_MARV = "IFNG+MARV", Mock_EVLP = "EVLP", Mock_MVLP = "MVLP")
COL <- c(Mock = "#7F7F7F", EBOV = "#1B7837", IFNB = "#2B6CB0", `IFNB+EBOV` = "#9BD08E", IFNG = "#D69E2E", `IFNG+EBOV` = "#E8C77A", MARV = "#8C510A", `IFNB+MARV` = "#C7A27C", `IFNG+MARV` = "#DCC9A8", EVLP = "#6A51A3", MVLP = "#CBA4D8")
NCOL <- as.integer(Sys.getenv("NCOL", "6")); TITLE <- Sys.getenv("TITLE", ""); INK <- "#1F2933"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs"); QUANT <- file.path(PROJ, "salmon_quant_v2"); GTF <- file.path(PROJ, "ref", "Mus_musculus.GRCm39.109.gtf"); OUT <- file.path(RNA, "figures_gene_counts"); dir.create(OUT, showWarnings = FALSE)
TAG <- Sys.getenv("OUT_NAME", ""); suf <- function(x) if (nzchar(TAG)) paste0(x, "_", TAG) else x   # OUT_NAME=zeb2 -> reads_per_sample_zeb2.png etc. (does not overwrite earlier runs)
save2 <- function(p, name, w, h) { name <- suf(name); ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300); ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h) }
nov <- function(x) sub("\\..*$", "", x); gl <- function(g) ifelse(g %in% names(LAB), LAB[g], g)
theme_t <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK), plot.title = element_text(colour = INK, face = "bold", size = 11), plot.subtitle = element_text(colour = "grey35", size = 8.5), strip.background = element_rect(fill = "grey95", colour = "grey70"), strip.text = element_text(face = "bold.italic", colour = INK, size = 8))

dds <- readRDS(file.path(RNA, "dds_v2_decoyAware_hostNormalised.rds")); cd <- as.data.frame(colData(dds)); cd$sample <- rownames(cd); cd$group <- as.character(cd$group)
sf <- sizeFactors(dds); if (is.null(sf)) { nf <- normalizationFactors(dds); sf <- exp(colMeans(log(nf))); sf <- sf / exp(mean(log(sf))); message("size factors derived from the normalizationFactors (tximport object)") }; names(sf) <- colnames(dds); kept <- nov(rownames(dds))
GRP <- if (identical(tolower(Sys.getenv("GROUPS", "Mock_Mock,IFNB_Mock,IFNY_Mock")), "all")) ORD[ORD %in% cd$group] else trimws(strsplit(Sys.getenv("GROUPS", "Mock_Mock,IFNB_Mock,IFNY_Mock"), ",")[[1]]); cd <- cd |> filter(group %in% GRP)
if (!nrow(cd)) stop("no samples in the requested groups")

# ---- gene -> transcripts, one pass over the GTF ----
pat <- tempfile(); writeLines(sprintf('gene_name "%s"', GENES), pat); message("looking the genes up in the GTF (one pass) ...")
ln <- tryCatch(readLines(pipe(sprintf("grep -F -f '%s' '%s'", pat, GTF))), error = function(e) character(0), warning = function(w) character(0)); f <- strsplit(ln, "\t"); f <- f[vapply(f, length, 1L) >= 9]; f <- f[vapply(f, function(x) x[3] == "transcript", TRUE)]
gv <- function(x, k) sub(paste0('.*', k, ' "([^"]+)".*'), "\\1", x)
tab <- tibble(gene = vapply(f, function(x) gv(x[9], "gene_name"), ""), gene_id = nov(vapply(f, function(x) gv(x[9], "gene_id"), "")), tx = nov(vapply(f, function(x) gv(x[9], "transcript_id"), ""))) |> filter(tolower(gene) %in% tolower(GENES))
miss <- GENES[!tolower(GENES) %in% tolower(tab$gene)]; if (length(miss)) message("not found in the GTF under that symbol: ", paste(miss, collapse = ", "))
if (!nrow(tab)) stop("none of the genes were found"); gene_ids <- tab |> distinct(gene, gene_id); gene_ids <- gene_ids |> arrange(match(tolower(gene), tolower(GENES))); gene_ids$gene <- factor(gene_ids$gene, levels = unique(gene_ids$gene))
# ---- counts from every sample's quant.sf ----
rows <- list(); message("reading ", nrow(cd), " quant.sf files ...")
for (s in cd$sample) { q <- read_tsv(file.path(QUANT, s, "quant.sf"), col_types = "cdddd", progress = FALSE); q$tx <- nov(q$Name)
  m <- q |> inner_join(tab |> select(gene, gene_id, tx), by = "tx") |> group_by(gene, gene_id) |> summarise(reads = sum(NumReads), TPM = sum(TPM), n_transcripts = n(), .groups = "drop"); m$sample <- s; rows[[s]] <- m }
L <- bind_rows(rows) |> right_join(expand_grid(gene_ids, sample = cd$sample), by = c("gene", "gene_id", "sample")) |> mutate(reads = coalesce(reads, 0), TPM = coalesce(TPM, 0)) |>
  left_join(cd |> select(sample, group), by = "sample") |> mutate(norm = reads / sf[sample], passes_prefilter = gene_id %in% kept, group = factor(gl(group), levels = gl(GRP)), gene = factor(gene, levels = levels(gene_ids$gene)))
write_tsv(L |> arrange(gene, group, sample), file.path(OUT, paste0(suf("gene_read_counts_long"), ".tsv")))
S <- L |> group_by(gene, gene_id, passes_prefilter, group) |> summarise(n = n(), mean_reads = mean(reads), min_reads = min(reads), max_reads = max(reads), mean_TPM = mean(TPM), mean_norm = mean(norm), samples_ge10_norm = sum(norm >= 10), .groups = "drop")
write_tsv(S, file.path(OUT, paste0(suf("gene_read_counts_by_group"), ".tsv")))
cat("\nmean raw reads / mean TPM per group (samples with >= 10 normalized reads in brackets); 'kept' = passed the pre-filter and is in the VST object\n")
print(as.data.frame(S |> mutate(cell = sprintf("%s reads, %.2f TPM (%d/%d)", format(round(mean_reads, 1), big.mark = ","), mean_TPM, samples_ge10_norm, n)) |> select(gene, kept = passes_prefilter, group, cell) |> pivot_wider(names_from = group, values_from = cell)), row.names = FALSE)

# ---- figures ----
lab <- L |> distinct(gene, passes_prefilter) |> mutate(panel = paste0(gene, ifelse(passes_prefilter, "", "  [removed by pre-filter]"))); L <- L |> left_join(lab |> select(gene, panel), by = "gene"); L$panel <- factor(L$panel, levels = lab$panel[match(levels(L$gene), lab$gene)])
ng <- n_distinct(L$gene); ncl <- min(ng, NCOL); nrw <- ceiling(ng / ncl); w <- 1.2 + ncl * (0.6 * length(GRP) + 0.9); h <- 1.8 + 2.7 * nrw; ttl <- if (nzchar(TITLE)) TITLE else paste0(ng, " genes")
brk <- c(0, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000, 300000)
base <- function(y, ylab, sub) ggplot(L, aes(group, .data[[y]], fill = group)) + geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.8, colour = INK, linewidth = 0.35) + geom_point(position = position_jitter(width = 0.08, height = 0, seed = 1), size = 1.8, shape = 21, colour = INK, fill = "white", stroke = 0.5) +
  facet_wrap(~panel, ncol = ncl, scales = "free_y") + scale_fill_manual(values = COL, guide = "none") + scale_y_continuous(trans = scales::pseudo_log_trans(base = 10), breaks = brk, labels = scales::label_comma(), expand = expansion(mult = c(0.05, 0.12))) + labs(subtitle = sub, x = NULL, y = ylab) + theme_t + theme(axis.text.x = element_text(angle = 40, hjust = 1))
p1 <- base("reads", "raw reads (log scale)", "every dot = one sample; number of reads actually assigned to the gene (Salmon NumReads, summed over its transcripts)") + labs(title = paste0(ttl, ": raw reads per sample"))
save2(p1, "reads_per_sample", w, h)
p2 <- base("norm", "normalized reads (reads / size factor, log scale)", "dashed line = the pre-filter cut-off (10 normalized reads; a gene had to reach it in >= 3 of the 33 samples to be kept)") + geom_hline(yintercept = 10, linetype = "dashed", colour = "#C53030") + labs(title = paste0(ttl, ": normalized reads vs the pre-filter"))
save2(p2, "normalized_counts_vs_filter", w, h)
p3 <- base("TPM", "TPM (log scale)", "transcripts per million: length-corrected, comparable between genes within a sample") + labs(title = paste0(ttl, ": TPM"))
save2(p3, "tpm_by_group", w, h)
# one combined figure: raw reads (top row) and normalized reads (bottom row), same groups, for a handful of genes
if (ng <= NCOL) {
  Lc <- L |> pivot_longer(c(reads, norm), names_to = "measure", values_to = "value") |> mutate(measure = factor(measure, levels = c("reads", "norm"), labels = c("raw reads", "normalized reads (reads / size factor)")))
  mm <- Lc |> group_by(panel, group, measure) |> summarise(top = max(value), m = mean(value), .groups = "drop") |> mutate(txt = format(round(m), big.mark = ","))
  pc <- ggplot(Lc, aes(group, value, fill = group)) + geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.8, colour = INK, linewidth = 0.35) + geom_point(position = position_jitter(width = 0.08, height = 0, seed = 1), size = 1.8, shape = 21, colour = INK, fill = "white", stroke = 0.5) +
    geom_text(data = mm, mapping = aes(x = group, y = top, label = txt), inherit.aes = FALSE, vjust = -0.9, size = 2.7, colour = INK) + facet_grid(measure ~ panel, scales = "free_y") + scale_fill_manual(values = COL, guide = "none") +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10), breaks = brk, labels = scales::label_comma(), expand = expansion(mult = c(0.05, 0.18))) +
    labs(title = paste0(ttl, ": raw and normalized reads"), subtitle = "every dot = one sample; number above each box = mean; log scale; normalized = raw reads / DESeq2 size factor (host-normalized)", x = NULL, y = NULL) + theme_t + theme(axis.text.x = element_text(angle = 40, hjust = 1))
  save2(pc, "reads_and_normalized", 1.4 + ncl * (0.62 * length(GRP) + 0.9), 7.4) }
message("\nSaved to: ", OUT)
