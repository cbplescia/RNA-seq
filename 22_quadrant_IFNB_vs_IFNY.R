#!/usr/bin/env Rscript
# Four-quadrant scatter of protein-coding genes: x = log2FC IFNB vs Mock, y = log2FC
# IFNY vs Mock, one dot per gene.
#   top right    up in both          bottom left   down in both
#   bottom right up with IFNB, down/flat with IFNY     top left  the reverse
# Dots above the diagonal respond more to IFNY, below it more to IFNB.
# Colour = response class (significance vs Mock, padj<0.05 & |log2FC|>1, from
# 09_ifnb_vs_ifny.R). The top N_LABEL genes per class, farthest from the origin, are
# labelled. Needs 09 to have been run, and ggrepel.

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })
if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Install ggrepel (conda install -c conda-forge r-ggrepel)")

N_LABEL <- 5; LFC_CUT <- 1
INK <- "#1F2933"
COL <- c(`Up in both` = "#C53030", `Down in both` = "#6B46C1", `IFNB only` = "#2B6CB0",
         `IFNγ only` = "#D69E2E", `Opposite directions` = "#2F855A", `Neither` = "#C8C8C8")

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs")
FIG <- file.path(OUT, "figures_protein_coding")
dir.create(FIG, showWarnings = FALSE)

d <- read_csv(file.path(OUT, "DE_v2_IFNB_vs_IFNY_protein_coding.csv"), show_col_types = FALSE) |>
  transmute(gene_name, x = log2FC_IFNB_vs_Mock, y = log2FC_IFNY_vs_Mock, response) |>
  filter(!is.na(x), !is.na(y)) |>
  mutate(class = factor(case_when(response == "both_up" ~ "Up in both",
                                  response == "both_down" ~ "Down in both",
                                  grepl("^IFNB_only", response) ~ "IFNB only",
                                  grepl("^IFNY_only", response) ~ "IFNγ only",
                                  response == "opposite_direction" ~ "Opposite directions",
                                  TRUE ~ "Neither"), levels = names(COL)),
         dist = sqrt(x^2 + y^2))
lab <- d |> filter(class != "Neither") |> group_by(class) |> slice_max(dist, n = N_LABEL) |> ungroup()
message("Genes plotted: ", nrow(d)); print(table(d$class))
for (cl in levels(lab$class)) message(cl, ": ", paste(lab$gene_name[lab$class == cl], collapse = ", "))

lim <- ceiling(max(abs(c(d$x, d$y))) * 1.05)
corner <- data.frame(x = c(lim, -lim, -lim, lim) * 0.97, y = c(lim, lim, -lim, -lim) * 0.97,
                     hj = c(1, 0, 0, 1), vj = c(1, 1, 0, 0),
                     t = c("Up in both", "Up with IFNγ, not IFNB", "Down in both", "Up with IFNB, not IFNγ"))
p <- ggplot(d, aes(x, y)) +
  geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey55") +
  geom_hline(yintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", colour = "grey80", linewidth = 0.3) +
  geom_point(data = d |> filter(class == "Neither"), colour = COL[["Neither"]], size = 0.7, alpha = 0.5) +
  geom_point(data = d |> filter(class != "Neither"), aes(colour = class), size = 1, alpha = 0.7) +
  geom_point(data = lab, aes(colour = class), size = 1.8) +
  ggrepel::geom_text_repel(data = lab, aes(label = gene_name), colour = INK, size = 3,
                           max.overlaps = Inf, min.segment.length = 0, segment.colour = "grey60",
                           segment.size = 0.25, box.padding = 0.4, seed = 1) +
  geom_text(data = corner, aes(x, y, label = t, hjust = hj, vjust = vj), colour = "grey45",
            size = 3, fontface = "italic", inherit.aes = FALSE) +
  scale_colour_manual(values = COL, breaks = setdiff(names(COL), "Neither"), name = NULL) +
  coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(title = "Response to IFNB vs IFNγ: protein-coding genes",
       subtitle = sprintf("each dot is one gene (shrunk log2FC vs Mock); colours = significant (padj<0.05 & |log2FC|>%g) vs Mock; dotted line = equal response; top %d per class labelled",
                          LFC_CUT, N_LABEL),
       x = "log2 fold-change, IFNB vs Mock", y = "log2 fold-change, IFNγ vs Mock") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey94", linewidth = 0.3),
        panel.border = element_rect(colour = "grey70"),
        axis.text = element_text(colour = INK), axis.title = element_text(colour = INK),
        plot.title = element_text(colour = INK, face = "bold"),
        plot.subtitle = element_text(colour = "grey35", size = 8.5),
        legend.position = "bottom")
ggsave(file.path(FIG, "quadrant_IFNB_vs_IFNY.png"), p, width = 8, height = 8.5, dpi = 300)
ggsave(file.path(FIG, "quadrant_IFNB_vs_IFNY.pdf"), p, width = 8, height = 8.5)
message("\nSaved to: ", FIG)
