#!/usr/bin/env Rscript
# Does the IFNB response depend on how thick a membrane each lipid species would make?
#
# Per-species proxy for bilayer thickness = EFFECTIVE acyl chain length:
#     effective carbons per chain = (total carbons - DB_PENALTY x double bonds) / number of chains
# Longer, more saturated chains -> thicker; each cis double bond kinks a chain and shortens it, thinning the bilayer.
# (Sphingolipids: the sphingoid base's one trans double bond is not counted - trans bonds do not kink.)
# An approximate RELATIVE thickness in Angstrom is also given: 2 x (1.5 + 1.265 x effective carbons), from Tanford's
# extended-chain length. Real bilayers are thinner and also depend on cholesterol, headgroup and temperature, so
# use it only to compare species, never as a measured thickness.
# DB_PENALTY is an assumption (default 1 carbon per double bond), so a sensitivity table over 0-2 is produced.
# Species are binned thin / intermediate / thick by THIN_BELOW and THICK_ABOVE (effective carbons per chain).
#
# Needs 34_lipidomics_chain_length.R to have been run (lipids_parsed.csv) and the original CSV (for abundances).
# Writes to <project>/lipidomics_chain_length/ (PNG + PDF + CSV):
#   J_log2FC_by_thickness_category   log2FC of thin / intermediate / thick species, per class
#   K_thickness_distribution         abundance-weighted thickness distribution, untreated vs IFNB (descriptive, n = 1)
#   L_penalty_sensitivity            Spearman rho (effective length vs log2FC) for penalties 0 - 2 per class
#   bilayer_thickness_stats.csv      per class: rho, p, BH, mean log2FC per category, weighted thickness shift
#   bilayer_thickness_species.csv    every species with its effective length, relative thickness and category

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })

DB_PENALTY <- 1.0; PENALTIES <- c(0, 0.5, 1, 1.5, 2)
THIN_BELOW <- 15.5; THICK_ABOVE <- 17.5      # effective carbons per chain
MIN_SPECIES <- 5
SPHINGO <- c("SM", "Cer", "Hex1Cer", "Hex2Cer", "GM1", "GM2", "GM3", "ST", "SPH")

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "lipidomics_chain_length")
pf <- file.path(OUT, "lipids_parsed.csv")
if (!file.exists(pf)) stop("Run 34_lipidomics_chain_length.R first (missing ", pf, ")")
a <- commandArgs(trailingOnly = TRUE)
fname <- if (length(a)) a[1] else "01_lipidomics_ranked_IFNB_vs_Untreated.csv"
cand <- c(fname, file.path(PROJ, fname), file.path(getwd(), fname))
rf <- cand[file.exists(cand)][1]
if (is.na(rf)) stop("Cannot find ", fname, " - give its full path as the first argument")
raw <- read_csv(rf, show_col_types = FALSE) |> select(lipid = Lipid_ID, untreated = Untreated_Norm, ifnb = B_Norm)

x <- read_csv(pf, show_col_types = FALSE) |>
  filter(!standard, !oxidized, !miss, !is.na(total_C), !is.na(lfc), !(class %in% c("Co", "Ch"))) |>
  left_join(raw, by = "lipid")
x$db_eff <- ifelse(x$class %in% SPHINGO, pmax(x$total_DB - 1, 0), x$total_DB)
eff_len <- function(pen) (x$total_C - pen * x$db_eff) / x$n_chains
x$eff_len <- eff_len(DB_PENALTY)
x$thickness_A <- 2 * (1.5 + 1.265 * x$eff_len)
x$category <- cut(x$eff_len, c(-Inf, THIN_BELOW, THICK_ABOVE, Inf), labels = c("thin", "intermediate", "thick"))
write_csv(x |> select(lipid, class, tier, lfc, total_C, total_DB, n_chains, eff_len, thickness_A, category, untreated, ifnb),
          file.path(OUT, "bilayer_thickness_species.csv"))
message(sprintf("Species: %d.  Effective chain length %.1f-%.1f carbons (median %.1f).  thin <= %.1f < intermediate <= %.1f < thick",
                nrow(x), min(x$eff_len), max(x$eff_len), median(x$eff_len), THIN_BELOW, THICK_ABOVE))

lyso <- x |> filter(tier != "Other")
big <- lyso |> count(class) |> filter(n >= MIN_SPECIES) |> pull(class)
lyso <- lyso |> filter(class %in% big)

# ---- per-class statistics ---------------------------------------------------------------------------------------
stats_row <- function(g) {
  ct <- suppressWarnings(cor.test(g$eff_len, g$lfc, method = "spearman", exact = FALSE))
  cat_mean <- tapply(g$lfc, g$category, mean); cat_n <- table(g$category)
  wu <- weighted.mean(g$eff_len, g$untreated); wi <- weighted.mean(g$eff_len, g$ifnb)
  tibble(class = g$class[1], tier = g$tier[1], n = nrow(g), spearman_rho = unname(ct$estimate), p = ct$p.value,
         mean_lfc_thin = cat_mean[["thin"]], n_thin = as.integer(cat_n[["thin"]]),
         mean_lfc_intermediate = cat_mean[["intermediate"]], n_intermediate = as.integer(cat_n[["intermediate"]]),
         mean_lfc_thick = cat_mean[["thick"]], n_thick = as.integer(cat_n[["thick"]]),
         weighted_len_untreated = wu, weighted_len_IFNB = wi, weighted_shift_carbons = wi - wu,
         weighted_shift_A_approx = 2 * 1.265 * (wi - wu))
}
st <- bind_rows(lapply(split(lyso, lyso$class), stats_row))
st$padj_BH <- p.adjust(st$p, "BH")
st$reading <- ifelse(st$padj_BH < 0.05 & st$spearman_rho < 0, "thinner-bilayer species relatively higher in IFNB",
              ifelse(st$padj_BH < 0.05 & st$spearman_rho > 0, "thicker-bilayer species relatively higher in IFNB", "no clear trend"))
write_csv(st, file.path(OUT, "bilayer_thickness_stats.csv"))
message("\nEffective chain length vs log2FC (negative rho = thinner-bilayer species higher in IFNB), penalty = ", DB_PENALTY)
print(as.data.frame(st |> transmute(class, n, rho = round(spearman_rho, 2), p = signif(p, 2), p_BH = signif(padj_BH, 2),
                                    wt_shift_C = round(weighted_shift_carbons, 3), reading)), row.names = FALSE)

sens <- bind_rows(lapply(PENALTIES, function(pen) {
  lyso$len <- (lyso$total_C - pen * lyso$db_eff) / lyso$n_chains
  bind_rows(lapply(split(lyso, lyso$class), function(g) {
    ct <- suppressWarnings(cor.test(g$len, g$lfc, method = "spearman", exact = FALSE))
    tibble(penalty = pen, class = g$class[1], rho = unname(ct$estimate), p = ct$p.value)
  }))
}))
write_csv(sens, file.path(OUT, "bilayer_thickness_penalty_sensitivity.csv"))

# ---- figures ---------------------------------------------------------------------------------------------------------
INK <- "#1F2933"
theme_f <- function() theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(),
  panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK),
  strip.background = element_rect(fill = "grey95", colour = "grey70"), strip.text = element_text(colour = INK, face = "bold"),
  plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 9),
  legend.position = "bottom")
save2 <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300, limitsize = FALSE)
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, limitsize = FALSE)
}
cat_cols <- c(thin = "#9ECAE1", intermediate = "#4292C6", thick = "#08519C")
nf <- length(big); ncol_f <- min(nf, 4); nrow_f <- ceiling(nf / 4)

# J
cn <- lyso |> count(class, category)
pJ <- ggplot(lyso, aes(category, lfc)) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_boxplot(aes(fill = category), outlier.shape = NA, alpha = 0.6, width = 0.6) +
  geom_jitter(width = 0.12, height = 0, size = 1.2, alpha = 0.6, colour = INK) +
  geom_text(data = cn, aes(y = -Inf, label = paste0("n=", n)), vjust = -0.6, size = 2.8, colour = "grey35") +
  facet_wrap(~class, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = cat_cols, guide = "none") +
  labs(title = "IFNB response by how thick a bilayer each species would make",
       subtitle = sprintf("thin <= %.1f, thick > %.1f effective carbons per chain (double-bond penalty %.1f per bond)", THIN_BELOW, THICK_ABOVE, DB_PENALTY),
       x = NULL, y = "log2 fold-change (IFNB vs untreated)") + theme_f()
save2(pJ, "J_log2FC_by_thickness_category", 3.2 * ncol_f + 1, 3.2 * nrow_f + 1.6)

# K
kd <- lyso |> select(class, eff_len, untreated, ifnb) |>
  pivot_longer(c(untreated, ifnb), names_to = "condition", values_to = "abund") |>
  group_by(class, condition) |> mutate(w = abund / sum(abund)) |> ungroup() |>
  mutate(condition = factor(ifelse(condition == "untreated", "untreated", "IFNB"), levels = c("untreated", "IFNB")))
km <- kd |> group_by(class, condition) |> summarise(m = weighted.mean(eff_len, abund), .groups = "drop")
pK <- ggplot(kd, aes(eff_len, weight = w, colour = condition, fill = condition)) +
  geom_density(alpha = 0.18, linewidth = 0.7, adjust = 1.2) +
  geom_vline(data = km, aes(xintercept = m, colour = condition), linetype = "dashed", linewidth = 0.5) +
  facet_wrap(~class, scales = "free", ncol = 4) +
  scale_colour_manual(values = c(untreated = "#8A8A8A", IFNB = "#2B6CB0"), name = NULL) +
  scale_fill_manual(values = c(untreated = "#8A8A8A", IFNB = "#2B6CB0"), name = NULL) +
  scale_x_continuous(sec.axis = sec_axis(~ 2 * (1.5 + 1.265 * .), name = "approx. relative thickness (Angstrom)")) +
  labs(title = "Abundance-weighted thickness distribution of each class",
       subtitle = "dashed = abundance-weighted mean; descriptive only (one measurement per condition, intensities are semi-quantitative)",
       x = "effective carbons per chain (lower = thinner bilayer)", y = "share of class signal") + theme_f()
save2(pK, "K_thickness_distribution", 3.6 * ncol_f + 1, 3.2 * nrow_f + 1.8)

# L
sens$class <- factor(sens$class, levels = rev(big)); sens$star <- ifelse(sens$p < 0.05, "*", "")
pL <- ggplot(sens, aes(factor(penalty), class, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = paste0(sprintf("%.2f", rho), star)), size = 3, colour = INK) +
  scale_fill_gradient2(low = "#2C5F9E", mid = "#F7F7F7", high = "#B2382C", midpoint = 0, limits = c(-0.6, 0.6),
                       oob = scales::squish, name = "Spearman rho") +
  labs(title = "How much does the result depend on the double-bond assumption?",
       subtitle = "rho of effective chain length vs log2FC; * = nominal p < 0.05. Negative (blue) = thinner-bilayer species higher in IFNB",
       x = "carbons removed per double bond (0 = ignore unsaturation)", y = NULL) + theme_f()
save2(pL, "L_penalty_sensitivity", 7.5, 0.4 * length(big) + 2.5)
message("\nFigures and tables in: ", OUT)
