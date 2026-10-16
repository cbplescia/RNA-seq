#!/usr/bin/env Rscript
# Lysosomal-membrane lipids in IFNB-treated cells: are SHORTER chains relatively more abundant?
#
# Input : 01_lipidomics_ranked_IFNB_vs_Untreated.csv  (Lipid_ID, Class, log2FC_IFNB_vs_Untreated, ...)
# Usage : Rscript 34_lipidomics_chain_length.R [path/to/file.csv]
# log2FC > 0 = higher in IFNB. This file has no per-lipid p-values, so the tests below ask whether the
# PATTERN across species in a class is chain-length dependent; they do not test single lipids.
#
# Method
#  1. Drop internal standards (names containing +D7:(s), -D7 ...), oxidised species (COOH / OOH in the
#     name) and rows with missing values in untreated or IFNB (their log2FC is inflated by zeros).
#     Each can be switched back on below.
#  2. Chain measure = total acyl carbons summed over all chains in the name (PC(16:0_18:1) and PC(34:1)
#     both give 34, so differently formatted names agree). For pooling across classes, carbons per chain =
#     total / number of chains for that class.
#  3. Lysosomal-membrane classes (edit the lists below):
#       core    : LBPA (= BMP), SM, Cer, Hex1Cer / Hex2Cer, gangliosides GM1-GM3, sulfatide (ST), Ch
#       general : PC, PE, PS, PI, PG, LPC, LPE - in the lysosomal limiting membrane but not specific to it
#  4. Per class: Spearman rho and regression slope (log2FC per extra carbon, 95% CI) of log2FC on chain
#     length, plus mean log2FC of shorter vs longer half. BH-adjusted across classes. NEGATIVE slope =
#     shorter chains relatively higher in IFNB.
#  5. Pooled sets use class-adjusted regression, and "core vs other" tests whether the slopes differ.

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })

EXCLUDE_STANDARDS <- TRUE; EXCLUDE_OXIDIZED <- TRUE; EXCLUDE_MISSING <- TRUE
MIN_SPECIES <- 5; ALPHA <- 0.05
CORE    <- c("LBPA", "BMP", "SM", "Cer", "Hex1Cer", "Hex2Cer", "GM1", "GM2", "GM3", "ST", "Ch")
GENERAL <- c("PC", "PE", "PS", "PI", "PG", "LPC", "LPE")
NCH <- c(TG = 3, CL = 4, LPC = 1, LPE = 1, LPI = 1, LPG = 1, LPA = 1, ChE = 1, FA = 1, AcCa = 1, ZyE = 1, SPH = 1)

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "lipidomics_chain_length")
dir.create(OUT, showWarnings = FALSE)

a <- commandArgs(trailingOnly = TRUE)
fname <- if (length(a)) a[1] else "01_lipidomics_ranked_IFNB_vs_Untreated.csv"
cand <- c(fname, file.path(PROJ, fname), file.path(getwd(), fname))
f <- cand[file.exists(cand)][1]
if (is.na(f)) stop("Cannot find ", fname, " - give its full path as the first argument")
raw <- read_csv(f, show_col_types = FALSE)
need <- c("Lipid_ID", "Class", "log2FC_IFNB_vs_Untreated")
if (!all(need %in% names(raw))) stop("Expected columns not found: ", paste(setdiff(need, names(raw)), collapse = ", "))
message("Read ", nrow(raw), " lipids from ", f)

d <- raw |> transmute(lipid = Lipid_ID, class = Class, lfc = log2FC_IFNB_vs_Untreated,
                      miss = (if ("any_missing_Untreated" %in% names(raw)) any_missing_Untreated else FALSE) |
                             (if ("any_missing_IFNB" %in% names(raw)) any_missing_IFNB else FALSE))
d$standard <- grepl("\\+D[0-9]|-D[0-9]|\\(s\\)", d$lipid)
d$oxidized <- grepl("COOH|OOH", d$lipid)

chain_info <- function(s) {
  m <- regmatches(s, gregexpr("(?<![0-9])[0-9]{1,2}:[0-9]{1,2}", s, perl = TRUE))
  t(vapply(m, function(z) {
    if (!length(z)) return(c(NA_real_, NA_real_))
    p <- do.call(rbind, strsplit(z, ":", fixed = TRUE)); c(sum(as.numeric(p[, 1])), sum(as.numeric(p[, 2])))
  }, numeric(2)))
}
ci <- chain_info(d$lipid)
d$total_C <- ci[, 1]; d$total_DB <- ci[, 2]
d$n_chains <- ifelse(d$class %in% names(NCH), NCH[d$class], 2)
d$carbons_per_chain <- d$total_C / d$n_chains
d$tier <- case_when(d$class %in% CORE ~ "Core lysosomal", d$class %in% GENERAL ~ "General membrane", TRUE ~ "Other")
write_csv(d, file.path(OUT, "lipids_parsed.csv"))

message(sprintf("Excluded: %d standards, %d oxidised, %d with missing values; %d have no chain info (e.g. free cholesterol)",
                sum(d$standard), sum(d$oxidized & !d$standard), sum(d$miss), sum(is.na(d$total_C))))
keep <- !is.na(d$total_C) & !is.na(d$lfc)
if (EXCLUDE_STANDARDS) keep <- keep & !d$standard
if (EXCLUDE_OXIDIZED) keep <- keep & !d$oxidized
if (EXCLUDE_MISSING) keep <- keep & !d$miss
x <- d[keep, ]
message("Lipids analysed: ", nrow(x))
chol <- d |> filter(class == "Ch", !standard)
if (nrow(chol)) message("Free cholesterol (no chain length): log2FC = ", paste(round(chol$lfc, 2), collapse = ", "))

# ---- per-class statistics --------------------------------------------------------------------------------------
class_stats <- function(g, label, tier, xvar = "total_C", adjust = FALSE) {
  g <- g[!is.na(g[[xvar]]), ]
  n <- nrow(g)
  if (n < MIN_SPECIES || n_distinct(g[[xvar]]) < 3) return(tibble(group = label, tier = tier, n = n))
  xv <- g[[xvar]]
  ct <- suppressWarnings(cor.test(xv, g$lfc, method = "spearman", exact = FALSE))
  fit <- if (adjust) lm(lfc ~ xv + class, data = g) else lm(lfc ~ xv, data = g)
  ci <- suppressMessages(confint(fit))["xv", ]
  med <- median(xv); sh <- g$lfc[xv <= med]; lo <- g$lfc[xv > med]
  tibble(group = label, tier = tier, n = n, chain_min = min(xv), chain_max = max(xv),
         spearman_rho = unname(ct$estimate), p_spearman = ct$p.value,
         slope_per_carbon = unname(coef(fit)["xv"]), ci_low = ci[[1]], ci_high = ci[[2]],
         p_slope = summary(fit)$coefficients["xv", 4],
         mean_lfc_shorter_half = mean(sh), mean_lfc_longer_half = mean(lo), short_minus_long = mean(sh) - mean(lo))
}
lyso <- x |> filter(tier != "Other")
by_class <- bind_rows(lapply(split(lyso, lyso$class), function(g) class_stats(g, g$class[1], g$tier[1])))
tested <- !is.na(by_class$slope_per_carbon)
by_class$padj_BH <- NA_real_
by_class$padj_BH[tested] <- p.adjust(by_class$p_spearman[tested], "BH")

pooled <- bind_rows(
  class_stats(x |> filter(tier == "Core lysosomal"), "Core lysosomal (class-adjusted)", "Pooled", "carbons_per_chain", TRUE),
  class_stats(lyso, "All lysosomal-membrane (class-adjusted)", "Pooled", "carbons_per_chain", TRUE),
  class_stats(x |> filter(tier == "Other", class != "Co"), "Non-lysosomal lipids (class-adjusted)", "Pooled", "carbons_per_chain", TRUE))

# do core lysosomal lipids differ from the rest? (class-centred, interaction of chain length x group)
cc <- x |> filter(tier %in% c("Core lysosomal", "Other"), class != "Co") |> group_by(class) |>
  filter(n() >= MIN_SPECIES) |> mutate(xc = carbons_per_chain - mean(carbons_per_chain), yc = lfc - mean(lfc)) |> ungroup()
cc$grp <- factor(ifelse(cc$tier == "Core lysosomal", "core", "other"), levels = c("other", "core"))
inter <- summary(lm(yc ~ xc * grp, data = cc))$coefficients
slope_diff <- tibble(comparison = "slope(core) - slope(non-lysosomal)", estimate = inter["xc:grpcore", 1],
                     p = inter["xc:grpcore", 4])

res <- bind_rows(pooled, by_class)
p_use <- ifelse(is.na(res$padj_BH), res$p_slope, res$padj_BH)
res$verdict <- ifelse(is.na(res$slope_per_carbon), "too few species / chain lengths",
               ifelse(p_use < ALPHA & res$slope_per_carbon < 0, "SHORTER chains higher in IFNB",
               ifelse(p_use < ALPHA & res$slope_per_carbon > 0, "LONGER chains higher in IFNB", "no clear trend")))
write_csv(res, file.path(OUT, "chain_length_stats.csv"))
write_csv(slope_diff, file.path(OUT, "core_vs_other_slope_difference.csv"))

message("\n=== chain length vs log2FC  (negative slope = shorter chains higher in IFNB) ===")
print(as.data.frame(res |> transmute(group, n, rho = round(spearman_rho, 2), slope = round(slope_per_carbon, 3),
                                     p_adj = signif(p_use, 2), verdict)), row.names = FALSE)
message(sprintf("\nCore lysosomal vs non-lysosomal slope difference: %.3f log2FC per carbon, p = %.3g",
                slope_diff$estimate, slope_diff$p))

# ---- figures ----------------------------------------------------------------------------------------------------------
INK <- "#1F2933"
theme_f <- function() theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(),
  panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK),
  strip.background = element_rect(fill = "grey95", colour = "grey70"), strip.text = element_text(colour = INK, face = "bold"),
  plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 9),
  legend.position = "bottom")
ok_classes <- by_class$group[tested]
pl <- lyso |> filter(class %in% ok_classes)
if (nrow(pl)) {
  ann <- by_class |> filter(group %in% ok_classes) |>
    transmute(class = group, lab = sprintf("rho %.2f, slope %.3f/C\nadj. p %.2g", spearman_rho, slope_per_carbon, padj_BH))
  p1 <- ggplot(pl, aes(total_C, lfc)) +
    geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "#2B6CB0", fill = "#2B6CB0", alpha = 0.15, linewidth = 0.6) +
    geom_point(aes(colour = tier), size = 1.8, alpha = 0.85) +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab), hjust = -0.05, vjust = 1.2, size = 2.8,
              colour = INK, inherit.aes = FALSE) +
    facet_wrap(~class, scales = "free") +
    scale_colour_manual(values = c(`Core lysosomal` = "#C53030", `General membrane` = "#6B7280"), name = NULL) +
    labs(title = "Lysosomal-membrane lipids: log2FC (IFNB vs untreated) by chain length",
         subtitle = "each dot is one lipid species; downward slope = shorter chains relatively higher in IFNB",
         x = "total acyl carbons", y = "log2 fold-change (IFNB vs untreated)") + theme_f()
  nf <- length(ok_classes)
  ggsave(file.path(OUT, "lysosomal_lipids_logFC_vs_chain_length.png"), p1,
         width = 3.2 * min(nf, 4) + 1, height = 3 * ceiling(nf / 4) + 1.6, dpi = 300, limitsize = FALSE)
}
fp <- res |> filter(!is.na(slope_per_carbon)) |>
  mutate(tier = factor(tier, levels = c("Pooled", "Core lysosomal", "General membrane")))
if (nrow(fp)) {
  fp$group <- factor(fp$group, levels = rev(fp$group[order(fp$tier, fp$slope_per_carbon)]))
  p2 <- ggplot(fp, aes(slope_per_carbon, group, colour = tier)) +
    geom_vline(xintercept = 0, colour = "grey50") +
    geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0.25, linewidth = 0.5, orientation = "y") +
    geom_point(size = 2.5) +
    scale_colour_manual(values = c(Pooled = "#1F2933", `Core lysosomal` = "#C53030", `General membrane` = "#6B7280"),
                        name = NULL, drop = FALSE) +
    labs(title = "Does chain length predict the IFNB response?",
         subtitle = "change in log2FC per extra carbon (95% CI); left of zero = shorter chains higher in IFNB",
         x = "slope (log2FC per carbon)", y = NULL) + theme_f()
  ggsave(file.path(OUT, "chain_length_slopes.png"), p2, width = 8, height = 0.38 * nrow(fp) + 2.5, dpi = 300)
}
message("\nResults in: ", OUT)
