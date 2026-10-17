#!/usr/bin/env Rscript
# Delta-delta of the pilot lipidomics:  DD = (IFNB vs untreated) - (IFNy vs untreated)   (= IFNB vs IFNy).
#   DD > 0  the lipid went up more (or down less) with IFNB   -> IFNB-skewed
#   DD < 0  the lipid went up more (or down less) with IFNy   -> IFNy-skewed
#
# IMPORTANT - global offset: total lipid signal is ~4-fold lower in IFNy than in untreated (and ~1.3-fold lower in
# IFNB), so almost every species has DD > 0. To see genuine differences, species DD is also shown CENTRED on the
# median DD (CENTER = TRUE), and class / group changes are shown as a share of all lipid signal as well as totals.
#
# Figures (PNG + PDF + CSV in <project>/lipidomics_chain_length/):
#   DD_A_species_scatter    x = IFNB vs untreated, y = IFNy vs untreated, one dot per species, coloured by centred DD;
#                           grey line = DD of 0, dashed line = the median offset; top 10 each side labelled
#   DD_B_species_ranked     the 20 most IFNB-skewed and 20 most IFNy-skewed species by centred DD
#   DD_C_class_dumbbell     each class: IFNB vs untreated and IFNy vs untreated (all species summed), joined, DD printed
#   DD_D_group_dumbbell     the same for charge groups, packing-shape groups and the core lysosomal group
# Filters: standards removed; rows with missing values removed; species must reach MIN_GROUP_INTENSITY in untreated,
# IFNB and IFNy (otherwise a near-zero value sets the fold change) - those are saved to DD_species_excluded.csv.
# One measurement per condition and no p-values: descriptive only. Intensities are semi-quantitative.
# Usage: Rscript 44_lipidomics_delta_delta.R [path/to/01_lipidomics_ranked_IFNB_vs_Untreated.csv]

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Install ggrepel (conda install -c conda-forge r-ggrepel)")

CENTER <- TRUE; MIN_GROUP_INTENSITY <- 10000; N_LABEL <- 10; N_RANK <- 20; MIN_CLASS_SPECIES <- 3
CHARGE <- list(Anionic = c("PS", "PI", "PG", "PA", "cPA", "LBPA", "BMP", "CL", "ST", "GM1", "GM2", "GM3", "LPI", "LPG", "LPA", "LPS"),
               Zwitterionic = c("PC", "PE", "SM", "LPC", "LPE", "AcCa"),
               Neutral = c("TG", "DG", "Cer", "Hex1Cer", "Hex2Cer", "Ch", "ChE", "CE", "ZyE"))
SHAPE <- list(`Cone-shaped` = c("PE", "PA", "CL", "DG", "Cer", "Ch", "LBPA", "BMP"),
              `Cylindrical` = c("PC", "PS", "PI", "PG", "SM", "Hex1Cer", "Hex2Cer", "ST"),
              `Inverted-cone` = c("LPC", "LPE", "LPI", "LPG", "LPA", "LPS", "cPA", "AcCa", "GM1", "GM2", "GM3"),
              `Neutral storage (TG/CE)` = c("TG", "ChE", "CE", "ZyE"))
LYSO <- list(`Core lysosomal-membrane lipids` = c("LBPA", "BMP", "SM", "Cer", "Hex1Cer", "Hex2Cer", "GM1", "GM2", "GM3", "ST", "Ch"))
FAMILIES <- list(`Charge` = CHARGE, `Packing shape` = SHAPE, `Lysosome` = LYSO)
INK <- "#1F2933"; B_COL <- "#2B6CB0"; Y_COL <- "#D69E2E"

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "lipidomics_chain_length"); dir.create(OUT, showWarnings = FALSE)
save2 <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300, limitsize = FALSE)
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, limitsize = FALSE)
}
theme_f <- function() theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.border = element_rect(colour = "grey70"),
  axis.text = element_text(colour = INK), strip.background = element_rect(fill = "grey95", colour = "grey70"),
  strip.text = element_text(colour = INK, face = "bold"), plot.title = element_text(colour = INK, face = "bold"),
  plot.subtitle = element_text(colour = "grey35", size = 8.5), legend.position = "bottom")

a <- commandArgs(trailingOnly = TRUE)
fname <- if (length(a)) a[1] else "01_lipidomics_ranked_IFNB_vs_Untreated.csv"
cand <- c(fname, file.path(PROJ, fname), file.path(getwd(), fname)); f <- cand[file.exists(cand)][1]
if (is.na(f)) stop("Cannot find ", fname, " - give its full path as the first argument")
raw <- read_csv(f, show_col_types = FALSE)
need <- c("Lipid_ID", "Class", "Untreated_Norm", "B_Norm", "Y_Norm", "log2FC_IFNB_vs_Untreated", "log2FC_IFNY_vs_Untreated")
if (!all(need %in% names(raw))) stop("Missing columns: ", paste(setdiff(need, names(raw)), collapse = ", "))
miss <- (if ("any_missing_Untreated" %in% names(raw)) raw$any_missing_Untreated else FALSE) |
        (if ("any_missing_IFNB" %in% names(raw)) raw$any_missing_IFNB else FALSE) |
        (if ("any_missing_IFNY" %in% names(raw)) raw$any_missing_IFNY else FALSE)
d <- raw |> transmute(lipid = Lipid_ID, class = Class, U = Untreated_Norm, B = B_Norm, Y = Y_Norm,
                      b = log2FC_IFNB_vs_Untreated, y = log2FC_IFNY_vs_Untreated, missing = miss)
d$standard <- grepl("\\+D[0-9]|-D[0-9]|\\(s\\)", d$lipid) | grepl("-D[0-9]$", d$class)
d$dd <- d$b - d$y

# ---- species level ------------------------------------------------------------------------------------------------
sp <- d |> filter(!standard, !missing, is.finite(dd))
det <- sp$U >= MIN_GROUP_INTENSITY & sp$B >= MIN_GROUP_INTENSITY & sp$Y >= MIN_GROUP_INTENSITY
write_csv(sp[!det, ] |> select(lipid, class, U, B, Y, b, y, dd), file.path(OUT, "DD_species_excluded.csv"))
sp <- sp[det, ]
offset <- if (CENTER) median(sp$dd) else 0
sp$dd_centred <- sp$dd - offset
message(sprintf("Species plotted: %d (excluded %d below the detection floor %d). Median DD = %+.2f", nrow(sp), sum(!det), MIN_GROUP_INTENSITY, median(sp$dd)))
write_csv(sp |> arrange(desc(dd_centred)) |> select(lipid, class, b_IFNB_vs_untreated = b, y_IFNy_vs_untreated = y, DD_raw = dd, DD_centred = dd_centred),
          file.path(OUT, "DD_species.csv"))
lab <- bind_rows(sp |> slice_max(dd_centred, n = N_LABEL, with_ties = FALSE), sp |> slice_min(dd_centred, n = N_LABEL, with_ties = FALSE))
lim <- ceiling(max(abs(c(sp$b, sp$y))) * 1.05)
pA <- ggplot(sp, aes(b, y)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) + geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, colour = "grey55", linewidth = 0.4) +
  { if (CENTER) geom_abline(slope = 1, intercept = -offset, colour = "black", linetype = "dashed", linewidth = 0.4) } +
  geom_point(aes(colour = dd_centred), size = 1.8, alpha = 0.85) +
  ggrepel::geom_text_repel(data = lab, aes(label = lipid), size = 2.6, colour = INK, max.overlaps = Inf, min.segment.length = 0,
                           segment.colour = "grey60", segment.size = 0.25, box.padding = 0.3, seed = 1) +
  scale_colour_gradient2(low = Y_COL, mid = "#BDBDBD", high = B_COL, midpoint = 0, limits = c(-3, 3), oob = scales::squish,
                         name = ifelse(CENTER, "DD (centred)", "DD")) +
  coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(title = "Delta-delta: IFNB vs untreated against IFNγ vs untreated, per lipid species",
       subtitle = sprintf("grey line: equal change in both; dashed: the median offset (%+.2f), i.e. the global shift; blue = IFNB-skewed, orange = IFNγ-skewed; top %d each side labelled", offset, N_LABEL),
       x = "log2FC, IFNB vs untreated", y = "log2FC, IFNγ vs untreated") + theme_f()
save2(pA, "DD_A_species_scatter", 9, 9)

rk <- bind_rows(sp |> slice_max(dd_centred, n = N_RANK, with_ties = FALSE), sp |> slice_min(dd_centred, n = N_RANK, with_ties = FALSE)) |>
  arrange(dd_centred) |> mutate(skew = ifelse(dd_centred >= 0, "IFNB-skewed", "IFNγ-skewed"))
rk$lab <- factor(make.unique(rk$lipid), levels = make.unique(rk$lipid))
pB <- ggplot(rk, aes(dd_centred, lab, fill = skew)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) + geom_col(width = 0.75) +
  scale_fill_manual(values = c(`IFNB-skewed` = B_COL, `IFNγ-skewed` = Y_COL), name = NULL) +
  labs(title = "Most IFNB-skewed and IFNγ-skewed lipid species", x = ifelse(CENTER, "delta-delta, centred on the median offset (log2)", "delta-delta (log2)"), y = NULL,
       subtitle = "DD = (IFNB vs untreated) - (IFNγ vs untreated)") + theme_f()
save2(pB, "DD_B_species_ranked", 7.5, 0.24 * nrow(rk) + 2.2)

# ---- class and group level -----------------------------------------------------------------------------------------
dd_tab <- function(tab, tot) {
  tab |> mutate(b_total = log2(B / U), y_total = log2(Y / U),
                b_share = log2((B / tot[["B"]]) / (U / tot[["U"]])), y_share = log2((Y / tot[["Y"]]) / (U / tot[["U"]])),
                dd_total = b_total - y_total, dd_share = b_share - y_share)
}
all_d <- d |> filter(!standard)
ct <- all_d |> group_by(class) |> summarise(n = n(), U = sum(U), B = sum(B), Y = sum(Y), .groups = "drop")
tot <- c(U = sum(ct$U), B = sum(ct$B), Y = sum(ct$Y))
cl <- dd_tab(ct, tot) |> filter(n >= MIN_CLASS_SPECIES)
gl <- bind_rows(lapply(names(FAMILIES), function(fam) bind_rows(lapply(names(FAMILIES[[fam]]), function(g) {
  s <- ct |> filter(class %in% FAMILIES[[fam]][[g]]); tibble(family = fam, class = g, n = sum(s$n), U = sum(s$U), B = sum(s$B), Y = sum(s$Y))
})))) |> filter(n > 0) |> dd_tab(tot)
write_csv(cl |> select(class, n_species = n, b_total, y_total, dd_total, b_share, y_share, dd_share), file.path(OUT, "DD_classes.csv"))
write_csv(gl |> select(family, group = class, n_species = n, b_total, y_total, dd_total, b_share, y_share, dd_share), file.path(OUT, "DD_groups.csv"))
message(sprintf("Global shift: IFNB %+.2f, IFNγ %+.2f log2 (total lipid signal vs untreated)", log2(tot[["B"]] / tot[["U"]]), log2(tot[["Y"]] / tot[["U"]])))

dumbbell <- function(t, title, ncol_label) {
  lg <- t |> pivot_longer(c(b_total, y_total, b_share, y_share), names_to = c("who", "basis"), names_sep = "_", values_to = "lfc") |>
    mutate(basis = factor(ifelse(basis == "total", "total signal", "share of total lipid signal"), levels = c("total signal", "share of total lipid signal")),
           who = factor(ifelse(who == "b", "IFNB vs untreated", "IFNγ vs untreated"), levels = c("IFNB vs untreated", "IFNγ vs untreated")))
  ddl <- t |> pivot_longer(c(dd_total, dd_share), names_to = "basis", values_to = "dd") |>
    mutate(basis = factor(ifelse(basis == "dd_total", "total signal", "share of total lipid signal"), levels = levels(lg$basis)))
  ddl <- left_join(ddl, lg |> group_by(class, basis) |> summarise(xmax = max(lfc), .groups = "drop"), by = c("class", "basis"))
  ord <- t |> arrange(dd_share) |> pull(class)
  lg$class <- factor(lg$class, levels = ord); ddl$class <- factor(ddl$class, levels = ord)
  ggplot(lg, aes(lfc, class)) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_line(aes(group = class), colour = "grey55", linewidth = 0.8) +
    geom_point(aes(colour = who), size = 3) +
    geom_text(data = ddl, aes(x = xmax, label = sprintf("  DD %+.2f", dd)), hjust = 0, size = 2.9, colour = INK) +
    facet_wrap(~basis) +
    scale_colour_manual(values = c(`IFNB vs untreated` = B_COL, `IFNγ vs untreated` = Y_COL), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.3))) +
    labs(title = title, subtitle = "all species of a class summed; DD = IFNB change minus IFNγ change; classes ordered by DD of the share",
         x = "log2 fold-change vs untreated", y = NULL) + theme_f()
}
save2(dumbbell(cl, "Class-level delta-delta", 1), "DD_C_class_dumbbell", 11, 0.3 * nrow(cl) + 2.5)
gl2 <- gl |> mutate(class = paste0(class, " (", n, ")"))
save2(dumbbell(gl2, "Delta-delta for biophysical lipid groups (number of species in brackets)", 1), "DD_D_group_dumbbell", 11, 0.35 * nrow(gl2) + 2.5)

message("\nClasses by delta-delta of share (negative = IFNγ-skewed):")
print(as.data.frame(cl |> arrange(dd_share) |> transmute(class, species = n, IFNB = round(b_share, 2), IFNy = round(y_share, 2), DD = round(dd_share, 2))), row.names = FALSE)
message("\nFigures and tables in: ", OUT)
