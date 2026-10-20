#!/usr/bin/env Rscript
# Pie chart: of all IFNB ISGs (IFNB vs Mock: padj < 0.05, log2FC > 1), how many are
#   (i)  little to not suppressed by EBOV   vs   (ii) suppressed by EBOV
# Uses figures_mock_ifnb_ebov/transcriptome_blocks.csv written by 51_mock_ifnb_ebov_analysis.R (run 51 first).
# retention = log2FC(IFNB+EBOV vs Mock) / log2FC(IFNB vs Mock)   (1 = full induction kept, 0 = gone)
#   SUPPRESSED (ii)               retention <= SUPP_RETENTION (default 0.5, i.e. EBOV removes at least half of the induction) AND the
#                                 IFNB+EBOV-vs-IFNB contrast moves the gene back toward Mock with padj < 0.05
#   LITTLE TO NOT SUPPRESSED (i)  every other ISG (retention above the cutoff, or no significant effect of EBOV on the IFNB response)
# A second pie splits the ISGs three ways with the stricter cut-offs used in the heatmaps (suppressed: retention <= 0.25 and significant;
# not suppressed: retention >= 0.75; partly suppressed: in between). Change the cut-offs below.
# Writes to figures_mock_ifnb_ebov/: pie_ISGs_suppressed_vs_not (.png/.pdf), pie_ISGs_three_way (.png/.pdf),
# ISG_suppression_summary.csv and gene lists genes_ISG_<group>.csv
# Usage: Rscript 54_pie_isg_suppression.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter

SUPP_RETENTION <- 0.5; R_STRONG <- 0.25; R_FAIL <- 0.75; PADJ_CUT <- 0.05
INK <- "#1F2933"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs", "figures_mock_ifnb_ebov")
bf <- file.path(OUT, "transcriptome_blocks.csv")
if (!file.exists(bf)) stop("Run 51_mock_ifnb_ebov_analysis.R first (missing ", bf, ")")
B <- read_csv(bf, show_col_types = FALSE)
isg <- B |> filter(resp == "Up", !is.na(retain), !is.na(lfc_BE_B)) |>
  mutate(toward_mock = sign(lfc_BE_B) == -sign(lfc_B) & !is.na(padj_BE_B) & padj_BE_B < PADJ_CUT,
         two_way = ifelse(retain <= SUPP_RETENTION & toward_mock, "Suppressed by EBOV", "Little to not suppressed"),
         three_way = case_when(retain <= R_STRONG & toward_mock ~ "Suppressed", retain >= R_FAIL ~ "Not suppressed", TRUE ~ "Partly suppressed"))
n <- nrow(isg)
message(sprintf("IFNB ISGs (padj<%g, log2FC>1): %d", PADJ_CUT, n))

tab2 <- isg |> count(two_way) |> mutate(pct = 100 * n / sum(n))
tab3 <- isg |> count(three_way) |> mutate(pct = 100 * n / sum(n))
write_csv(bind_rows(tab2 |> transmute(split = sprintf("two-way (suppressed: retention <= %.2f and significant)", SUPP_RETENTION), group = two_way, genes = n, percent = round(pct, 1)),
                    tab3 |> transmute(split = sprintf("three-way (suppressed <= %.2f, not suppressed >= %.2f)", R_STRONG, R_FAIL), group = three_way, genes = n, percent = round(pct, 1))),
          file.path(OUT, "ISG_suppression_summary.csv"))
print(as.data.frame(tab2), row.names = FALSE); print(as.data.frame(tab3), row.names = FALSE)
for (g in unique(isg$two_way)) write_csv(isg |> filter(two_way == g) |> arrange(desc(lfc_B)) |> select(gene_id, gene_name, lfc_B, padj_B, lfc_BE, retain, lfc_BE_B, padj_BE_B),
                                         file.path(OUT, paste0("genes_ISG_", gsub("[^A-Za-z]+", "_", g), ".csv")))

pie <- function(tab, grp, cols, title, subtitle) {
  d <- tab |> mutate(group = factor(.data[[grp]], levels = names(cols))) |> arrange(group) |>
    mutate(lab = sprintf("%s\n%d genes (%.0f%%)", group, n, pct))
  ggplot(d, aes(x = "", y = n, fill = group)) +
    geom_col(width = 1, colour = "white", linewidth = 1.2) + coord_polar(theta = "y", start = 0) +
    geom_text(aes(label = lab), position = position_stack(vjust = 0.5), size = 4, colour = ifelse(d$group %in% c("Suppressed by EBOV", "Suppressed"), "white", INK), fontface = "bold") +
    scale_fill_manual(values = cols, guide = "none") +
    labs(title = title, subtitle = subtitle) + theme_void(base_size = 12) +
    theme(plot.title = element_text(colour = INK, face = "bold", hjust = 0.5), plot.subtitle = element_text(colour = "grey35", size = 9, hjust = 0.5),
          plot.margin = margin(10, 10, 10, 10))
}
p2 <- pie(tab2, "two_way", c(`Little to not suppressed` = "#F1A340", `Suppressed by EBOV` = "#762A83"),
          sprintf("IFNB ISGs (n = %d): how many does EBOV suppress?", n),
          sprintf("suppressed = EBOV removes at least %.0f%% of the IFNB induction (significant); everything else = little to not suppressed", 100 * (1 - SUPP_RETENTION)))
ggsave(file.path(OUT, "pie_ISGs_suppressed_vs_not.png"), p2, width = 6.5, height = 6.5, dpi = 300)
ggsave(file.path(OUT, "pie_ISGs_suppressed_vs_not.pdf"), p2, width = 6.5, height = 6.5)
p3 <- pie(tab3, "three_way", c(`Not suppressed` = "#F1A340", `Partly suppressed` = "#D9C3E8", Suppressed = "#762A83"),
          sprintf("IFNB ISGs (n = %d): three-way split", n),
          sprintf("suppressed: >= %.0f%% of induction lost (significant); not suppressed: >= %.0f%% kept; partly: in between", 100 * (1 - R_STRONG), 100 * R_FAIL))
ggsave(file.path(OUT, "pie_ISGs_three_way.png"), p3, width = 6.5, height = 6.5, dpi = 300)
ggsave(file.path(OUT, "pie_ISGs_three_way.pdf"), p3, width = 6.5, height = 6.5)
message("\nSaved to: ", OUT)
