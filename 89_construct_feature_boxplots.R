#!/usr/bin/env Rscript
# Box plots of every feature (CDS / transcribed unit) of your recombinant construct, for each condition. The construct is the set of records in ref/recombinant_model.fa, quantified by Salmon next
# to the mouse transcriptome and kept in the DESeq2 object (size factors come from host genes only), so its features are the rows of the VST / count tables that are not ENSMUSG gene IDs.
#   construct_boxplot_VST.png / .pdf        VST expression (log2 scale, host-normalised) per feature, one panel per feature, one box per condition, replicates as points
#   construct_boxplot_TPM.png / .pdf        log10(TPM + 1) from the Salmon quantification (construct_features_per_sample.csv) - the same plot on a scale that does not depend on the VST;
#                                           reads per million of the whole library, so it is comparable between conditions but TPM sums to 1e6 within each sample
#   construct_feature_values.csv            VST, TPM and reads for every feature in every sample
# Mock (no construct) samples should sit near zero; if they do not, the construct shares sequence with mouse transcripts (see viral_refs/README.md) and that feature is not interpretable alone.
# Conditions are the 11 groups in the data (Mock, IFNB, IFNG, EBOV, IFNB+EBOV, IFNG+EBOV, MARV, IFNB+MARV, IFNG+MARV, EVLP, MVLP); restrict them with GROUPS, e.g.
#   GROUPS=Mock_Mock,Mock_EBOV,IFNB_EBOV Rscript 89_construct_feature_boxplots.R
# Writes to requant_v2_outputs/figures_construct/.   Usage: Rscript 89_construct_feature_boxplots.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2); library(DESeq2) })
select <- dplyr::select; filter <- dplyr::filter

ORDER <- c("Mock_Mock", "IFNB_Mock", "IFNY_Mock", "Mock_EBOV", "IFNB_EBOV", "IFNY_EBOV", "Mock_MARV", "IFNB_MARV", "IFNY_MARV", "Mock_EVLP", "Mock_MVLP")
LAB <- c(Mock_Mock = "Mock", IFNB_Mock = "IFNB", IFNY_Mock = "IFNG", Mock_EBOV = "EBOV", IFNB_EBOV = "IFNB+EBOV", IFNY_EBOV = "IFNG+EBOV", Mock_MARV = "MARV", IFNB_MARV = "IFNB+MARV", IFNY_MARV = "IFNG+MARV", Mock_EVLP = "EVLP", Mock_MVLP = "MVLP")
COL <- c(Mock = "#7F7F7F", IFNB = "#2B6CB0", IFNG = "#D69E2E", EBOV = "#1B7837", `IFNB+EBOV` = "#9BD08E", `IFNG+EBOV` = "#E8C77A", MARV = "#8C510A", `IFNB+MARV` = "#C7A27C", `IFNG+MARV` = "#DCC9A8", EVLP = "#762A83", MVLP = "#C2A5CF")
INK <- "#1F2933"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_construct"); dir.create(OUT, showWarnings = FALSE)
save2 <- function(p, name, w, h) { ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300); ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h) }

vsd <- readRDS(file.path(RNA, "vsd_v2.rds")); vm <- SummarizedExperiment::assay(vsd); grp <- setNames(as.character(vsd$group), colnames(vm))
feat <- rownames(vm)[!grepl("^ENSMUSG", rownames(vm))]
if (!length(feat)) stop("No construct features in the VST object (rows that are not ENSMUSG gene IDs). Check ", file.path(RNA, "construct_features_per_sample.csv"))
message("Construct features found (", length(feat), "): ", paste(feat, collapse = ", "))
G <- if (nzchar(Sys.getenv("GROUPS"))) trimws(strsplit(Sys.getenv("GROUPS"), ",")[[1]]) else ORDER
G <- G[G %in% grp]; if (!length(G)) stop("None of the requested groups are in the data: ", paste(unique(grp), collapse = ", "))
lab_of <- function(g) ifelse(g %in% names(LAB), LAB[g], g)

V <- as_tibble(vm[feat, , drop = FALSE], rownames = "feature") |> pivot_longer(-feature, names_to = "sample_id", values_to = "vst") |> mutate(group = grp[sample_id])
cf <- file.path(RNA, "construct_features_per_sample.csv")
if (file.exists(cf)) {
  T <- read_csv(cf, show_col_types = FALSE) |> transmute(feature, sample_id = as.character(sample_id), TPM, reads, length)
  V <- V |> left_join(T, by = c("feature", "sample_id"))
} else { message("construct_features_per_sample.csv not found - TPM plot skipped"); V$TPM <- NA_real_; V$reads <- NA_real_; V$length <- NA_real_ }
V <- V |> filter(group %in% G) |> mutate(cond = factor(lab_of(group), levels = lab_of(G)), feature = factor(feature, levels = feat), l10tpm = log10(TPM + 1))
write_csv(V |> arrange(feature, cond, sample_id) |> select(feature, sample_id, group, vst, TPM, reads, length), file.path(OUT, "construct_feature_values.csv"))
print(as.data.frame(V |> group_by(feature, cond) |> summarise(n = n(), mean_vst = round(mean(vst), 2), mean_TPM = round(mean(TPM), 2), mean_reads = round(mean(reads), 1), .groups = "drop")), row.names = FALSE)

cols <- COL[levels(V$cond)]; nf <- n_distinct(V$feature); ncol <- min(3, nf); nrow <- ceiling(nf / ncol)
mk <- function(y, ylab, title) ggplot(V, aes(cond, .data[[y]], fill = cond)) + geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8, colour = INK, linewidth = 0.35) +
  geom_point(position = position_jitter(width = 0.1, height = 0, seed = 1), size = 1.4, shape = 21, colour = INK, fill = "white", stroke = 0.4) +
  facet_wrap(~feature, scales = "free_y", ncol = ncol) + scale_fill_manual(values = cols, guide = "none") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(title = title, subtitle = "box = median and quartiles of the replicates; points = individual samples", x = NULL, y = ylab) +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3), panel.border = element_rect(colour = "grey70"),
        axis.text = element_text(colour = INK), axis.text.x = element_text(angle = 50, hjust = 1, size = 8), plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 8.5),
        strip.background = element_rect(fill = "grey95", colour = "grey70"), strip.text = element_text(colour = INK, face = "bold"))
w <- max(6, 0.45 * length(G) * ncol + 1.5); h <- 3.2 * nrow + 1.2
save2(mk("vst", "VST expression (log2 scale)", "Construct features: VST expression by condition"), "construct_boxplot_VST", w, h)
if (any(!is.na(V$TPM))) save2(mk("l10tpm", "log10(TPM + 1)", "Construct features: Salmon TPM by condition"), "construct_boxplot_TPM", w, h)
message("Saved to: ", OUT)
