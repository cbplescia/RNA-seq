#!/usr/bin/env Rscript
# Abundance of every construct feature (each transcribed unit / CDS of ref/recombinant_model.fa) in each condition, ranked from most to least abundant WITHIN each condition.
# One panel per condition; bar = mean TPM of the replicates (Salmon, construct_features_per_sample.csv), error bar = +/- 1 SD, points = individual replicates; every feature keeps the same colour
# in every panel so you can see how the ranking changes between conditions.
#   construct_abundance_TPM.png / .pdf         linear TPM scale
#   construct_abundance_log10TPM.png / .pdf    log10(TPM + 1), which shows the low-abundance features when one feature dominates
#   construct_abundance_ranked.csv             for every condition: feature, rank, mean TPM, SD, mean reads, % of the construct TPM in that condition
# TPM compares features WITHIN a sample (it corrects for transcript length and sums to 1e6 per sample); the percentage column is each feature's share of the construct signal in that condition.
# Mock-type conditions without the construct should be ~0; if not, the construct shares sequence with mouse transcripts (viral_refs/README.md).
# Restrict conditions with GROUPS, e.g.  GROUPS=Mock_EBOV,IFNB_EBOV Rscript 90_construct_abundance_by_condition.R
# Writes to requant_v2_outputs/figures_construct/.   Usage: Rscript 90_construct_abundance_by_condition.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter

ORDER <- c("Mock_Mock", "IFNB_Mock", "IFNY_Mock", "Mock_EBOV", "IFNB_EBOV", "IFNY_EBOV", "Mock_MARV", "IFNB_MARV", "IFNY_MARV", "Mock_EVLP", "Mock_MVLP")
LAB <- c(Mock_Mock = "Mock", IFNB_Mock = "IFNB", IFNY_Mock = "IFNG", Mock_EBOV = "EBOV", IFNB_EBOV = "IFNB+EBOV", IFNY_EBOV = "IFNG+EBOV", Mock_MARV = "MARV", IFNB_MARV = "IFNB+MARV", IFNY_MARV = "IFNG+MARV", Mock_EVLP = "EVLP", Mock_MVLP = "MVLP")
INK <- "#1F2933"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_construct"); dir.create(OUT, showWarnings = FALSE)
save2 <- function(p, name, w, h) { ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300); ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h) }

cf <- file.path(RNA, "construct_features_per_sample.csv")
if (!file.exists(cf)) stop("Missing ", cf, " - it is written by 04_tximport_deseq2.R when the construct FASTA is present")
T <- read_csv(cf, show_col_types = FALSE) |> mutate(group = as.character(group))
G <- if (nzchar(Sys.getenv("GROUPS"))) trimws(strsplit(Sys.getenv("GROUPS"), ",")[[1]]) else ORDER
G <- G[G %in% T$group]; if (!length(G)) stop("None of the requested groups are in the data: ", paste(unique(T$group), collapse = ", "))
lab_of <- function(g) ifelse(g %in% names(LAB), LAB[g], g)
T <- T |> filter(group %in% G) |> mutate(cond = factor(lab_of(group), levels = lab_of(G)))
feat <- unique(T$feature); message("Construct features (", length(feat), "): ", paste(feat, collapse = ", "))

S <- T |> group_by(cond, feature) |> summarise(mean_TPM = mean(TPM), sd_TPM = sd(TPM), mean_reads = mean(reads), n = n(), .groups = "drop") |>
  group_by(cond) |> arrange(desc(mean_TPM), .by_group = TRUE) |> mutate(rank = row_number(), pct_of_construct = ifelse(sum(mean_TPM) > 0, 100 * mean_TPM / sum(mean_TPM), NA_real_)) |> ungroup()
write_csv(S |> select(condition = cond, rank, feature, mean_TPM, sd_TPM, mean_reads, pct_of_construct, n_replicates = n), file.path(OUT, "construct_abundance_ranked.csv"))
print(as.data.frame(S |> group_by(cond) |> summarise(construct_TPM = round(sum(mean_TPM), 1), top_feature = feature[1], top_TPM = round(mean_TPM[1], 1), .groups = "drop")), row.names = FALSE)

fcol <- setNames(scales::hue_pal(l = 55, c = 90)(length(feat)), feat)
S <- S |> mutate(key = paste(cond, feature, sep = "||")); T <- T |> mutate(key = paste(cond, feature, sep = "||"))
S$key <- factor(S$key, levels = S$key[order(S$cond, S$rank)]); T$key <- factor(T$key, levels = levels(S$key))
mk <- function(log) {
  y <- if (log) "l10" else "TPM"; yb <- if (log) "log10(mean TPM + 1)" else "mean TPM"
  Sx <- S |> mutate(l10 = log10(mean_TPM + 1), lo = if (log) log10(pmax(mean_TPM - sd_TPM, 0) + 1) else pmax(mean_TPM - sd_TPM, 0), hi = if (log) log10(mean_TPM + sd_TPM + 1) else mean_TPM + sd_TPM)
  Tx <- T |> mutate(l10 = log10(TPM + 1))
  ggplot() + geom_col(data = Sx, aes(key, if (log) l10 else mean_TPM, fill = feature), width = 0.75) +
    geom_errorbar(data = Sx, aes(key, ymin = lo, ymax = hi), width = 0.25, colour = INK, linewidth = 0.3) +
    geom_point(data = Tx, aes(key, if (log) l10 else TPM), position = position_jitter(width = 0.08, height = 0, seed = 1), size = 1, shape = 21, colour = INK, fill = "white", stroke = 0.3) +
    facet_wrap(~cond, scales = "free_x", ncol = 3) + scale_x_discrete(labels = function(x) sub("^.*\\|\\|", "", x)) + scale_fill_manual(values = fcol, name = "construct feature") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + guides(fill = guide_legend(nrow = 2)) +
    labs(title = "Construct feature abundance in each condition, ranked high to low", subtitle = sprintf("bar = mean of the replicates, line = +/- 1 SD, points = replicates; Salmon TPM (%s scale)", if (log) "log10(TPM + 1)" else "linear"),
         x = NULL, y = yb) +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3), panel.border = element_rect(colour = "grey70"),
          axis.text = element_text(colour = INK), axis.text.x = element_text(angle = 50, hjust = 1, size = 7.5), plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 8.5),
          strip.background = element_rect(fill = "grey95", colour = "grey70"), strip.text = element_text(colour = INK, face = "bold"), legend.position = "bottom")
}
nr <- ceiling(length(G) / 3); w <- max(8, 0.32 * length(feat) * 3 + 2.5); h <- 2.9 * nr + 1.8
save2(mk(FALSE), "construct_abundance_TPM", w, h)
save2(mk(TRUE), "construct_abundance_log10TPM", w, h)
message("Saved to: ", OUT)
