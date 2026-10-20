#!/usr/bin/env Rscript
# Pie / donut charts of the whole protein-coding transcriptome:  (i) ISGs, (ii) IDGs, (iii) not IFNB-responsive, and
# which of the ISGs / IDGs have their IFNB change ANTAGONISED by EBOV.
#   ISG = IFNB vs Mock padj < 0.05 and log2FC > 1;  IDG = padj < 0.05 and log2FC < -1.
#   ANTAGONISED by EBOV = retention <= ANTAG_RETENTION (default 0.5: EBOV removes at least half of the change) AND the
#     IFNB+EBOV-vs-IFNB contrast moves the gene back toward Mock with padj < 0.05  (ISGs: induction suppressed; IDGs: repression reversed)
#     retention = log2FC(IFNB+EBOV vs Mock) / log2FC(IFNB vs Mock).  Everything else = little to not antagonised.
# Uses figures_mock_ifnb_ebov/transcriptome_blocks.csv written by 51_mock_ifnb_ebov_analysis.R (run 51 first).
# Writes to figures_mock_ifnb_ebov/ (PNG + PDF):
#   pie_protein_coding_response           all protein-coding genes: ISGs / IDGs / not IFNB-responsive
#   donut_IFNB_responsive_antagonised     ISGs and IDGs only: inner ring ISG vs IDG, outer ring antagonised vs little/not  (the readable one)
#   donut_all_protein_coding_antagonised  everything: inner ring ISG / IDG / not responsive, outer ring adds the EBOV effect (small slices unlabelled)
# plus transcriptome_response_summary.csv and genes_IDG_antagonised_by_EBOV.csv / genes_ISG_antagonised_by_EBOV.csv
# Usage: Rscript 55_pie_transcriptome_response.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter

ANTAG_RETENTION <- 0.5; PADJ_CUT <- 0.05
INK <- "#1F2933"
C_ISG <- "#C53030"; C_IDG <- "#2B6CB0"; C_NONE <- "#D9D9D9"; C_ANT <- "#762A83"; C_LIT <- "#F1A340"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs", "figures_mock_ifnb_ebov")
bf <- file.path(OUT, "transcriptome_blocks.csv")
if (!file.exists(bf)) stop("Run 51_mock_ifnb_ebov_analysis.R first (missing ", bf, ")")
B <- read_csv(bf, show_col_types = FALSE) |>
  mutate(toward_mock = !is.na(lfc_BE_B) & sign(lfc_BE_B) == -sign(lfc_B) & !is.na(padj_BE_B) & padj_BE_B < PADJ_CUT,
         antag = resp != "None" & !is.na(retain) & retain <= ANTAG_RETENTION & toward_mock)
N <- nrow(B)
nI <- sum(B$resp == "Up"); nD <- sum(B$resp == "Down"); nN <- sum(B$resp == "None")
aI <- sum(B$resp == "Up" & B$antag); aD <- sum(B$resp == "Down" & B$antag)
message(sprintf("Protein-coding genes: %d | ISGs %d | IDGs %d | not IFNB-responsive %d", N, nI, nD, nN))
message(sprintf("Antagonised by EBOV (retention <= %.2f, significant): ISGs %d of %d (%.0f%%), IDGs %d of %d (%.0f%%)",
                ANTAG_RETENTION, aI, nI, 100 * aI / nI, aD, nD, 100 * aD / nD))
summ <- tibble(group = c("ISGs", "ISGs: suppressed by EBOV", "ISGs: little to not suppressed", "IDGs", "IDGs: reversed by EBOV", "IDGs: little to not reversed",
                         "Not IFNB-responsive", "All protein-coding"),
               genes = c(nI, aI, nI - aI, nD, aD, nD - aD, nN, N))
write_csv(summ |> mutate(percent_of_parent = c(100 * nI / N, 100 * aI / nI, 100 * (nI - aI) / nI, 100 * nD / N, 100 * aD / nD, 100 * (nD - aD) / nD, 100 * nN / N, 100) |> round(1)),
          file.path(OUT, "transcriptome_response_summary.csv"))
for (r in c("Up", "Down")) write_csv(B |> filter(resp == r, antag) |> arrange(if (r == "Up") desc(lfc_B) else lfc_B) |>
  select(gene_id, gene_name, lfc_B, padj_B, lfc_BE, retain, lfc_BE_B, padj_BE_B),
  file.path(OUT, paste0("genes_", if (r == "Up") "ISG" else "IDG", "_antagonised_by_EBOV.csv")))

# ring builder: segments in order, cumulative angle
ring <- function(label, n, fill, txt, xmin, xmax, total = sum(n)) {
  tibble(label = label, n = n, fill = fill, txt = txt, xmin = xmin, xmax = xmax) |>
    mutate(ymax = cumsum(n), ymin = ymax - n, pct = 100 * n / total)
}
draw <- function(rings, title, subtitle, min_label = 2.5, hole = 0.4, file) {
  d <- bind_rows(rings)
  p <- ggplot(d) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), colour = "white", linewidth = 1) +
    geom_text(data = d |> filter(pct >= min_label), aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = sprintf("%s\n%d (%.0f%%)", label, n, pct), colour = txt),
              size = 3.1, fontface = "bold", lineheight = 0.9) +
    scale_fill_identity() + scale_colour_identity() +
    coord_polar(theta = "y") + scale_x_continuous(limits = c(0, max(d$xmax))) +
    labs(title = title, subtitle = subtitle) + theme_void(base_size = 12) +
    theme(plot.title = element_text(colour = INK, face = "bold", hjust = 0.5), plot.subtitle = element_text(colour = "grey35", size = 8.5, hjust = 0.5),
          plot.margin = margin(10, 10, 10, 10))
  ggsave(file.path(OUT, paste0(file, ".png")), p, width = 8, height = 8, dpi = 300)
  ggsave(file.path(OUT, paste0(file, ".pdf")), p, width = 8, height = 8)
}

# 1. simple pie of everything
draw(list(ring(c("ISGs", "IDGs", "Not IFNB-\nresponsive"), c(nI, nD, nN), c(C_ISG, C_IDG, C_NONE), c("white", "white", INK), 0, 1)),
     sprintf("All protein-coding transcripts (n = %d)", N),
     sprintf("ISG: padj < %g, log2FC > 1;  IDG: padj < %g, log2FC < -1 (IFNB vs Mock). ISGs %d, IDGs %d, not responsive %d", PADJ_CUT, PADJ_CUT, nI, nD, nN),
     min_label = 0, file = "pie_protein_coding_response")

# 2. responsive genes only, with the EBOV antagonism ring (the readable one)
inner <- ring(c("ISGs", "IDGs"), c(nI, nD), c(C_ISG, C_IDG), c("white", "white"), 1, 2)
outer <- ring(c("ISGs\nsuppressed", "ISGs\nlittle/not", "IDGs\nreversed", "IDGs\nlittle/not"), c(aI, nI - aI, aD, nD - aD),
              c(C_ANT, C_LIT, C_ANT, C_LIT), c("white", INK, "white", INK), 2.05, 3.05)
draw(list(inner, outer), sprintf("IFNB-responsive genes (n = %d): which does EBOV antagonise?", nI + nD),
     sprintf("inner ring: ISGs vs IDGs; outer ring: purple = change antagonised by EBOV (>= %.0f%% removed, significant), orange = little to not antagonised",
             100 * (1 - ANTAG_RETENTION)), min_label = 1, file = "donut_IFNB_responsive_antagonised")

# 3. everything, nested
inner3 <- ring(c("ISGs", "IDGs", "Not IFNB-\nresponsive"), c(nI, nD, nN), c(C_ISG, C_IDG, C_NONE), c("white", "white", INK), 1, 2)
outer3 <- ring(c("ISGs suppressed", "ISGs little/not", "IDGs reversed", "IDGs little/not", "Not responsive"), c(aI, nI - aI, aD, nD - aD, nN),
               c(C_ANT, C_LIT, C_ANT, C_LIT, "#EFEFEF"), c("white", INK, "white", INK, INK), 2.05, 3.05)
draw(list(inner3, outer3), sprintf("All protein-coding transcripts (n = %d) and the EBOV effect on the IFNB response", N),
     sprintf("outer ring: purple = antagonised by EBOV, orange = little to not; slices under 2.5%% are unlabelled (ISGs %d: %d suppressed; IDGs %d: %d reversed)", nI, aI, nD, aD),
     min_label = 2.5, file = "donut_all_protein_coding_antagonised")
message("\nSaved to: ", OUT)
