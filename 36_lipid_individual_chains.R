#!/usr/bin/env Rscript
# Individual-chain view of the IFNB lipid response: which acyl CHAINS (16:0, 18:1, 22:6 ...) and which
# chain PAIRINGS (e.g. 18:1/18:1) go up or down - instead of the summed "total carbons".
# Run 34_lipidomics_chain_length.R first (reads its lipids_parsed.csv; same filters).
#
# Only species whose name lists every chain can be used (e.g. LBPA(18:1_22:6), PS(18:0_20:3)). Names that give
# only the sum (PC(34:1), SM(d43:5)) say nothing about the individual chains and are left out here. Sphingolipids
# (SM, Cer, ...) are skipped: their "chains" are a sphingoid base plus an N-acyl chain, which are different things.
#
# Writes to <project>/lipidomics_chain_length/ (PNG + PDF):
#   F_chain_effects_<class>      mean log2FC of all species containing each individual chain
#   G_chain_pair_heatmap         chain A x chain B for each diacyl class; tile = log2FC (number printed if <= 50 species)
#   H_ranked_species_<class>     every species of small classes ranked by log2FC, labelled with its chains
#   I_shorter_vs_longer_chain    log2FC against the SHORTER and the LONGER chain length of each species
#   chain_level_stats.csv        Spearman rho of log2FC vs shorter / longer / mean chain length per class
# log2FC > 0 (red) = higher in IFNB.

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })

MIN_SPECIES <- 8     # explicit-chain species needed for a class to be plotted/tested
MAX_BAR <- 40        # classes up to this size also get a ranked bar chart
SPHINGO <- c("SM", "Cer", "Hex1Cer", "Hex2Cer", "GM1", "GM2", "GM3", "ST", "SPH")
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "lipidomics_chain_length")
pf <- file.path(OUT, "lipids_parsed.csv")
if (!file.exists(pf)) stop("Run 34_lipidomics_chain_length.R first (missing ", pf, ")")
x <- read_csv(pf, show_col_types = FALSE) |> filter(!standard, !oxidized, !miss, !is.na(lfc), tier != "Other")

INK <- "#1F2933"; UP <- "#C53030"; DOWN <- "#2B6CB0"
theme_f <- function() theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(),
  panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK),
  strip.background = element_rect(fill = "grey95", colour = "grey70"), strip.text = element_text(colour = INK, face = "bold"),
  plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 9),
  legend.position = "bottom")
save2 <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300, limitsize = FALSE)
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, limitsize = FALSE)
}

# ---- species whose name lists every chain, as long format ---------------------------------------------------------
found <- regmatches(x$lipid, gregexpr("(?<![0-9])[0-9]{1,2}:[0-9]{1,2}", x$lipid, perl = TRUE))
x$n_found <- lengths(found)
x$explicit <- x$n_found >= x$n_chains & x$n_found >= 1 & !(x$class %in% SPHINGO)
ex <- x[x$explicit, ]
found <- found[x$explicit]
ex$chains <- lapply(found, function(z) {
  p <- do.call(rbind, strsplit(z, ":", fixed = TRUE)); o <- order(as.numeric(p[, 1]), as.numeric(p[, 2])); z[o]
})
ex$pair <- vapply(ex$chains, paste, character(1), collapse = "/")
cnt <- ex |> count(class) |> filter(n >= MIN_SPECIES)
message("Classes with at least ", MIN_SPECIES, " explicit-chain species: ", paste0(cnt$class, " (", cnt$n, ")", collapse = ", "))
ex <- ex[ex$class %in% cnt$class, ]
if (!nrow(ex)) stop("No class has enough species with explicit chains.")
chain_levels <- unique(unlist(ex$chains))
cp <- do.call(rbind, strsplit(chain_levels, ":", fixed = TRUE))
chain_levels <- chain_levels[order(as.numeric(cp[, 1]), as.numeric(cp[, 2]))]

# ---- F: effect of each individual chain --------------------------------------------------------------------------
long <- ex |> select(lipid, class, lfc, chains) |> tidyr::unnest(chains) |> rename(chain = chains) |>
  distinct(lipid, chain, .keep_all = TRUE)
eff <- long |> group_by(class, chain) |>
  summarise(mean_lfc = mean(lfc), se = ifelse(n() > 1, sd(lfc) / sqrt(n()), NA_real_), n = n(), .groups = "drop") |>
  filter(n >= 2) |> mutate(chain = factor(chain, levels = chain_levels), dir = ifelse(mean_lfc >= 0, "higher in IFNB", "lower in IFNB"))
for (cl in unique(eff$class)) {
  e <- eff[eff$class == cl, ]
  p <- ggplot(e, aes(chain, mean_lfc, colour = dir)) +
    geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
    geom_linerange(aes(ymin = mean_lfc - se, ymax = mean_lfc + se), linewidth = 0.6, na.rm = TRUE) +
    geom_point(aes(size = n)) +
    scale_colour_manual(values = c(`higher in IFNB` = UP, `lower in IFNB` = DOWN), name = NULL) +
    scale_size_continuous(range = c(1.5, 5), name = "species with\nthis chain") +
    labs(title = paste0(cl, ": log2FC of species containing each chain"),
         subtitle = "chains ordered by length then double bonds; bars = SE; only chains found in 2+ species",
         x = "acyl chain (carbons:double bonds)", y = "mean log2FC (IFNB vs untreated)") + theme_f() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save2(p, paste0("F_chain_effects_", cl), max(6, 0.45 * nrow(e) + 2.5), 4.8)
}

# ---- G: chain A x chain B heatmap -------------------------------------------------------------------------------
di <- ex[ex$n_chains == 2 & lengths(ex$chains) == 2, ]
if (nrow(di)) {
  di$A <- vapply(di$chains, `[`, character(1), 1); di$B <- vapply(di$chains, `[`, character(1), 2)
  gt <- di |> group_by(class, A, B) |> summarise(lfc = mean(lfc), n = n(), .groups = "drop") |>
    mutate(A = factor(A, levels = chain_levels), B = factor(B, levels = chain_levels))
  sz <- gt |> group_by(class) |> summarise(tot = sum(n))
  gt <- gt |> left_join(sz, by = "class") |> mutate(lab = ifelse(tot <= 50, sprintf("%.1f", lfc), NA_character_))
  pG <- ggplot(gt, aes(A, B, fill = lfc)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = lab), size = 2.4, colour = INK, na.rm = TRUE) +
    facet_wrap(~class, scales = "free") +
    scale_fill_gradient2(low = "#2C5F9E", mid = "#F7F7F7", high = "#B2382C", midpoint = 0, limits = c(-1.5, 1.5),
                         oob = scales::squish, name = "log2FC") +
    labs(title = "Which chain pairings go up or down with IFNB?",
         subtitle = "each tile is one species, e.g. x = 16:0, y = 18:1 is 16:0/18:1; blank = not measured; shorter-chain partner on x",
         x = "first (shorter) chain", y = "second (longer) chain") + theme_f() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.background = element_rect(fill = "grey95"),
          panel.grid = element_blank())
  save2(pG, "G_chain_pair_heatmap", 5 * min(n_distinct(gt$class), 3) + 1, 4.8 * ceiling(n_distinct(gt$class) / 3) + 1)
}

# ---- H: ranked species for small classes ---------------------------------------------------------------------
for (cl in cnt$class[cnt$n <= MAX_BAR]) {
  s <- ex[ex$class == cl, ] |> arrange(lfc) |> mutate(dir = ifelse(lfc >= 0, "higher in IFNB", "lower in IFNB"))
  s$label <- make.unique(s$pair); s$label <- factor(s$label, levels = s$label)
  p <- ggplot(s, aes(lfc, label, fill = dir)) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) + geom_col(width = 0.75) +
    scale_fill_manual(values = c(`higher in IFNB` = UP, `lower in IFNB` = DOWN), name = NULL) +
    labs(title = paste0(cl, " species ranked by log2FC"), subtitle = "label = the two acyl chains (shorter first)",
         x = "log2 fold-change (IFNB vs untreated)", y = NULL) + theme_f()
  save2(p, paste0("H_ranked_species_", cl), 6.5, 0.24 * nrow(s) + 2)
}

# ---- I: shorter vs longer chain, and statistics ------------------------------------------------------------------
sl <- di |> mutate(shorter = as.numeric(sub(":.*$", "", A)), longer = as.numeric(sub(":.*$", "", B)),
                   mean_chain = (shorter + longer) / 2)
st <- bind_rows(lapply(split(sl, sl$class), function(g) {
  f <- function(v) suppressWarnings(cor.test(v, g$lfc, method = "spearman", exact = FALSE))
  a <- f(g$shorter); b <- f(g$longer); m <- f(g$mean_chain)
  tibble(class = g$class[1], n = nrow(g),
         rho_shorter_chain = unname(a$estimate), p_shorter = a$p.value,
         rho_longer_chain = unname(b$estimate), p_longer = b$p.value,
         rho_mean_chain = unname(m$estimate), p_mean = m$p.value)
}))
for (cn in c("p_shorter", "p_longer", "p_mean")) st[[sub("^p", "padj", cn)]] <- p.adjust(st[[cn]], "BH")
write_csv(st, file.path(OUT, "chain_level_stats.csv"))
message("\nSpearman rho of log2FC vs chain length (explicit-chain species; negative = shorter chains higher in IFNB):")
print(as.data.frame(st |> transmute(class, n, shorter_chain = round(rho_shorter_chain, 2), p1 = signif(p_shorter, 2),
                                    longer_chain = round(rho_longer_chain, 2), p2 = signif(p_longer, 2),
                                    mean_chain = round(rho_mean_chain, 2), p3 = signif(p_mean, 2))), row.names = FALSE)

pl <- sl |> select(class, lfc, shorter, longer) |> pivot_longer(c(shorter, longer), names_to = "which", values_to = "carbons") |>
  mutate(which = factor(ifelse(which == "shorter", "shorter chain of the pair", "longer chain of the pair"),
                        levels = c("shorter chain of the pair", "longer chain of the pair")))
pI <- ggplot(pl, aes(carbons, lfc)) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "#2B6CB0", fill = "#2B6CB0", alpha = 0.15, linewidth = 0.6) +
  geom_jitter(width = 0.15, height = 0, size = 1.5, alpha = 0.7, colour = INK) +
  facet_grid(which ~ class, scales = "free_x") +
  labs(title = "Log2FC against the length of each individual chain",
       subtitle = "each dot is one species; carbons of its shorter (top) and longer (bottom) acyl chain",
       x = "carbons in that chain", y = "log2 fold-change (IFNB vs untreated)") + theme_f()
save2(pI, "I_shorter_vs_longer_chain", 2.8 * n_distinct(pl$class) + 1.5, 6.5)
message("\nFigures written to: ", OUT)
