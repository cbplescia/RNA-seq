#!/usr/bin/env Rscript
# Lipidomics at the CLASS and BIOPHYSICAL-GROUP level: every species of a class is summed into one value (all LBPA
# species = one point, etc.), and classes are grouped by charge (anionic / zwitterionic / neutral) and by packing
# shape (cone-shaped, cylindrical, inverted-cone, neutral storage), plus a "core lysosomal-membrane lipids" group.
#
# Each change is shown two ways (panels / rows):
#   total signal              log2(sum in condition A / sum in condition B)
#   share of total lipid      the same after dividing by the sum of ALL lipids in each condition. This removes an
#                             overall shift (e.g. IFNy total lipid signal is ~4-fold lower than untreated), so it
#                             asks whether a class gained or lost ground relative to the whole lipidome.
# Three contrasts: IFNB vs untreated, IFNy vs untreated, IFNB vs IFNy (positive = higher in the first-named).
# Intensities are semi-quantitative (ionisation differs between classes), so compare a class between conditions,
# not between classes. One measurement per condition: descriptive only, no statistics.
# The charge / shape assignments are the usual textbook ones and are approximate (curvature depends on the acyl
# chains and environment); edit CHARGE, SHAPE and LYSO_CORE below to change them.
#
# Usage: Rscript 43_lipidomics_class_and_property_groups.R [path/to/01_lipidomics_ranked_IFNB_vs_Untreated.csv]
# Writes to <project>/lipidomics_chain_length/ (PNG + PDF + CSV):
#   class_level_<contrast>            one bubble per class (size = number of species; colour = charge)
#   property_groups_<contrast>        bars for charge groups, shape groups and the core lysosomal group
#   property_group_composition        share of total signal of each group in untreated / IFNB / IFNy
#   class_level_changes.csv, property_group_changes.csv

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Install ggrepel (conda install -c conda-forge r-ggrepel)")

EXCLUDE_STANDARDS <- TRUE
CHARGE <- list(
  Anionic      = c("PS", "PI", "PG", "PA", "cPA", "LBPA", "BMP", "CL", "ST", "GM1", "GM2", "GM3", "LPI", "LPG", "LPA", "LPS"),
  Zwitterionic = c("PC", "PE", "SM", "LPC", "LPE", "AcCa"),
  Neutral      = c("TG", "DG", "Cer", "Hex1Cer", "Hex2Cer", "Ch", "ChE", "CE", "ZyE"))
SHAPE <- list(
  `Cone-shaped (negative curvature)`   = c("PE", "PA", "CL", "DG", "Cer", "Ch", "LBPA", "BMP"),
  `Cylindrical (bilayer-forming)`      = c("PC", "PS", "PI", "PG", "SM", "Hex1Cer", "Hex2Cer", "ST"),
  `Inverted-cone (positive curvature)` = c("LPC", "LPE", "LPI", "LPG", "LPA", "LPS", "cPA", "AcCa", "GM1", "GM2", "GM3"),
  `Neutral storage (TG / CE)`          = c("TG", "ChE", "CE", "ZyE"))
LYSO_CORE <- list(`Core lysosomal-membrane lipids` = c("LBPA", "BMP", "SM", "Cer", "Hex1Cer", "Hex2Cer", "GM1", "GM2", "GM3", "ST", "Ch"))
FAMILIES <- list(`Charge` = CHARGE, `Packing shape` = SHAPE, `Lysosome` = LYSO_CORE)

CONTRASTS <- list(
  list(id = "IFNB_vs_untreated", label = "IFNB vs untreated", num = "B", den = "U"),
  list(id = "IFNY_vs_untreated", label = "IFNγ vs untreated", num = "Y", den = "U"),
  list(id = "IFNB_vs_IFNY",      label = "IFNB vs IFNγ (positive = higher in IFNB)", num = "B", den = "Y"))

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "lipidomics_chain_length"); dir.create(OUT, showWarnings = FALSE)
INK <- "#1F2933"; UP <- "#C53030"; DOWN <- "#2B6CB0"
BASIS <- c(lfc_total = "total signal", lfc_relative = "share of total lipid signal")

a <- commandArgs(trailingOnly = TRUE)
fname <- if (length(a)) a[1] else "01_lipidomics_ranked_IFNB_vs_Untreated.csv"
cand <- c(fname, file.path(PROJ, fname), file.path(getwd(), fname)); f <- cand[file.exists(cand)][1]
if (is.na(f)) stop("Cannot find ", fname, " - give its full path as the first argument")
raw <- read_csv(f, show_col_types = FALSE)
need <- c("Lipid_ID", "Class", "Untreated_Norm", "B_Norm", "Y_Norm")
if (!all(need %in% names(raw))) stop("Missing columns: ", paste(setdiff(need, names(raw)), collapse = ", "))
d <- raw |> transmute(lipid = Lipid_ID, class = Class, U = Untreated_Norm, B = B_Norm, Y = Y_Norm)
std <- grepl("\\+D[0-9]|-D[0-9]|\\(s\\)", d$lipid) | grepl("-D[0-9]$", d$class)
if (EXCLUDE_STANDARDS) d <- d[!std, ]
message("Lipid species used: ", nrow(d), " in ", n_distinct(d$class), " classes (", sum(std), " standards ", if (EXCLUDE_STANDARDS) "excluded" else "kept", ")")

ct <- d |> group_by(class) |> summarise(n = n(), U = sum(U), B = sum(B), Y = sum(Y), .groups = "drop")
tot <- c(U = sum(ct$U), B = sum(ct$B), Y = sum(ct$Y))
message(sprintf("Total lipid signal vs untreated: IFNB %+.2f log2, IFNγ %+.2f log2", log2(tot[["B"]] / tot[["U"]]), log2(tot[["Y"]] / tot[["U"]])))

charge_of <- function(cl) { for (g in names(CHARGE)) if (cl %in% CHARGE[[g]]) return(g); "Other / unassigned" }
unassigned <- ct$class[vapply(ct$class, charge_of, character(1)) == "Other / unassigned"]
if (length(unassigned)) message("Classes not assigned to a charge group: ", paste(unassigned, collapse = ", "))
ct$charge <- vapply(ct$class, charge_of, character(1))

change <- function(tab, num, den) {
  tab |> mutate(lfc_total = log2(.data[[num]] / .data[[den]]),
                lfc_relative = log2((.data[[num]] / tot[[num]]) / (.data[[den]] / tot[[den]])),
                abundance = log10((.data[[num]] + .data[[den]]) / 2))
}

# ---- group-level sums -----------------------------------------------------------------------------------------------
gl <- bind_rows(lapply(names(FAMILIES), function(fam) {
  bind_rows(lapply(names(FAMILIES[[fam]]), function(g) {
    s <- ct |> filter(class %in% FAMILIES[[fam]][[g]])
    tibble(family = fam, group = g, n_classes = nrow(s), n = sum(s$n), U = sum(s$U), B = sum(s$B), Y = sum(s$Y))
  }))
})) |> filter(n > 0)

ch_all <- list(); gr_all <- list()
for (cn in CONTRASTS) {
  cc <- change(ct, cn$num, cn$den) |> mutate(contrast = cn$label)
  gg <- change(gl, cn$num, cn$den) |> mutate(contrast = cn$label)
  ch_all[[cn$id]] <- cc; gr_all[[cn$id]] <- gg

  # class bubbles
  cl <- cc |> pivot_longer(c(lfc_total, lfc_relative), names_to = "basis", values_to = "lfc") |>
    mutate(basis = factor(BASIS[basis], levels = BASIS), charge = factor(charge, levels = c(names(CHARGE), "Other / unassigned")))
  xr <- max(abs(cl$lfc[is.finite(cl$lfc)])) * 1.1
  p <- ggplot(cl, aes(lfc, abundance)) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
    geom_point(aes(size = n, colour = charge), alpha = 0.8) +
    ggrepel::geom_text_repel(aes(label = class), size = 3, colour = INK, max.overlaps = Inf, min.segment.length = 0,
                             segment.colour = "grey70", segment.size = 0.25, box.padding = 0.3, seed = 1) +
    facet_wrap(~basis) +
    scale_colour_manual(values = c(Anionic = "#7B3294", Zwitterionic = "#4D9DE0", Neutral = "#6B7280", `Other / unassigned` = "#C8C8C8"), name = "charge", drop = FALSE) +
    scale_size_area(max_size = 11, breaks = c(1, 10, 50, 250), name = "species in class") +
    scale_x_continuous(limits = c(-xr, xr)) +
    labs(title = paste0("Lipid classes (all species summed): ", cn$label),
         subtitle = "left panel = change in the class total; right panel = change in the class's share of all lipid signal (corrects for an overall shift)",
         x = paste0("log2 fold-change (", cn$label, ")"), y = "log10 mean class signal") +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.border = element_rect(colour = "grey70"),
      strip.background = element_rect(fill = "grey95", colour = "grey70"), strip.text = element_text(colour = INK, face = "bold"),
      axis.text = element_text(colour = INK), plot.title = element_text(colour = INK, face = "bold"),
      plot.subtitle = element_text(colour = "grey35", size = 8.5), legend.position = "bottom")
  ggsave(file.path(OUT, paste0("class_level_", cn$id, ".png")), p, width = 12, height = 6.5, dpi = 300)
  ggsave(file.path(OUT, paste0("class_level_", cn$id, ".pdf")), p, width = 12, height = 6.5)

  # property-group bars
  gb <- gg |> pivot_longer(c(lfc_total, lfc_relative), names_to = "basis", values_to = "lfc") |>
    mutate(basis = factor(BASIS[basis], levels = BASIS), family = factor(family, levels = names(FAMILIES)),
           dir = ifelse(lfc >= 0, "higher", "lower"), lab = paste0(group, "\n(", n, " species)"))
  gb$lab <- factor(gb$lab, levels = unique(gb$lab))
  pb <- ggplot(gb, aes(lab, lfc, fill = dir)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) + geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%+.2f", lfc), vjust = ifelse(lfc >= 0, -0.4, 1.3)), size = 3, colour = INK) +
    facet_grid(basis ~ family, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = c(higher = UP, lower = DOWN), guide = "none") +
    scale_x_discrete(labels = scales::label_wrap(16)) + scale_y_continuous(expand = expansion(mult = 0.18)) +
    labs(title = paste0("Biophysical lipid groups: ", cn$label),
         subtitle = "all species of the member classes summed; top row = group total, bottom row = group share of all lipid signal",
         x = NULL, y = "log2 fold-change") +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
      panel.border = element_rect(colour = "grey70"), strip.background = element_rect(fill = "grey95", colour = "grey70"),
      strip.text = element_text(colour = INK, face = "bold"), axis.text = element_text(colour = INK, size = 8.5),
      plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 8.5))
  ggsave(file.path(OUT, paste0("property_groups_", cn$id, ".png")), pb, width = 12, height = 7, dpi = 300)
  ggsave(file.path(OUT, paste0("property_groups_", cn$id, ".pdf")), pb, width = 12, height = 7)
}
write_csv(bind_rows(ch_all) |> select(contrast, class, charge, n_species = n, U, B, Y, lfc_total, lfc_relative), file.path(OUT, "class_level_changes.csv"))
write_csv(bind_rows(gr_all) |> select(contrast, family, group, n_species = n, U, B, Y, lfc_total, lfc_relative), file.path(OUT, "property_group_changes.csv"))

# ---- composition: share of total signal per group, per condition -----------------------------------------------------------
comp <- gl |> pivot_longer(c(U, B, Y), names_to = "cond", values_to = "signal") |>
  mutate(share = 100 * signal / tot[cond], cond = factor(c(U = "untreated", B = "IFNB", Y = "IFNγ")[cond], levels = c("untreated", "IFNB", "IFNγ")),
         group = factor(group, levels = unique(gl$group)))
pc <- ggplot(comp, aes(cond, share, fill = cond)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", share)), vjust = -0.3, size = 2.8, colour = INK) +
  facet_wrap(~group, scales = "free_y", labeller = labeller(group = scales::label_wrap_gen(28))) +
  scale_fill_manual(values = c(untreated = "#8A8A8A", IFNB = DOWN, `IFNγ` = "#D69E2E"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Share of total lipid signal carried by each group", x = NULL, y = "% of all lipid signal in that condition",
       subtitle = "each panel has its own y scale; a group can gain share while its own total falls") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    panel.border = element_rect(colour = "grey70"), strip.background = element_rect(fill = "grey95", colour = "grey70"),
    strip.text = element_text(colour = INK, face = "bold", size = 8.5), axis.text = element_text(colour = INK),
    plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 8.5))
ng <- n_distinct(comp$group)
ggsave(file.path(OUT, "property_group_composition.png"), pc, width = 12, height = 3.2 * ceiling(ng / 4) + 1.2, dpi = 300)
ggsave(file.path(OUT, "property_group_composition.pdf"), pc, width = 12, height = 3.2 * ceiling(ng / 4) + 1.2)

message("\nIFNB vs untreated, property groups (log2FC; total | share of all lipid):")
print(as.data.frame(gr_all[["IFNB_vs_untreated"]] |> transmute(family, group, species = n, total = round(lfc_total, 2), share = round(lfc_relative, 2))), row.names = FALSE)
message("\nFigures and tables in: ", OUT)
