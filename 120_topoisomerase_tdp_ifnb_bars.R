#!/usr/bin/env Rscript
# Bar graph of IFNB vs Mock differential expression for the DNA topoisomerases (Top1, Top1mt, Top2a, Top2b, Top3a, Top3b) and the tyrosyl-DNA phosphodiesterases (Tdp1, Tdp2), one bar per gene in two panels (edit GENES to add more).
# Shrunk log2FC from DE_v2_IFNB_vs_Mock_annotated.csv, error bars +/- 1 lfcSE, stars = padj (* < 0.05, ** < 0.01, *** < 0.001, ns), dagger on the gene name = baseMean < LOW_EXPR (few reads, log2FC less reliable). Genes missing from the table (removed by the
# pre-filter in 04, or a different symbol) are listed as "Not found". Symbols are matched ignoring case.
# Env: COMP (default IFNB_vs_Mock; e.g. COMP=EBOV_vs_Mock for another contrast in the same layout).   Writes bars_topoisomerase_tdp.png / .pdf and topoisomerase_tdp_ifnb_vs_mock.csv to requant_v2_outputs/figures_gene_bars/.
# Usage: Rscript 120_topoisomerase_tdp_ifnb_bars.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter

GENES <- list("DNA topoisomerases" = c("Top1", "Top1mt", "Top2a", "Top2b", "Top3a", "Top3b"), "Tyrosyl-DNA phosphodiesterases" = c("Tdp1", "Tdp2"))
COMP <- Sys.getenv("COMP", "IFNB_vs_Mock"); COMP_LAB <- gsub("_vs_", " vs ", gsub("_", "+", COMP)); COMP_LAB <- gsub("\\+vs\\+", " vs ", COMP_LAB)
BARCOL <- c(IFNB_vs_Mock = "#2B6CB0", EBOV_vs_Mock = "#1B7837", IFNB_EBOV_vs_Mock = "#9BD08E"); FILL <- if (COMP %in% names(BARCOL)) BARCOL[[COMP]] else "#2B6CB0"
LOW_EXPR <- 50; INK <- "#1F2933"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs"); OUT <- file.path(RNA, "figures_gene_bars"); dir.create(OUT, showWarnings = FALSE)
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

tab <- read_csv(file.path(RNA, paste0("DE_v2_", COMP, "_annotated.csv")), show_col_types = FALSE) |> mutate(key = tolower(gene_name))
want <- bind_rows(lapply(names(GENES), function(g) tibble(class = g, query = GENES[[g]], pos = seq_along(GENES[[g]]))))
D <- want |> mutate(key = tolower(query)) |> left_join(tab |> filter(!is.na(gene_name)) |> distinct(key, .keep_all = TRUE) |> select(key, gene_id, gene = gene_name, baseMean, lfc = log2FoldChange, se = lfcSE, padj), by = "key")
miss <- D$query[is.na(D$gene_id)]; if (length(miss)) message("Not found in the DE table (removed by the pre-filter, or a different symbol): ", paste(miss, collapse = ", "))
D <- D |> filter(!is.na(gene_id)); if (!nrow(D)) stop("None of the genes were found in ", COMP)
D <- D |> mutate(class = factor(class, levels = names(GENES)), lab = stars(padj), vj = ifelse(lfc >= 0, 0, 1), ys = lfc + ifelse(lfc >= 0, 1, -1) * (se + 0.03), xlab = paste0(gene, ifelse(baseMean < LOW_EXPR, "†", "")))
D <- D |> arrange(class, pos); D$xlab <- factor(D$xlab, levels = unique(D$xlab))
print(as.data.frame(D |> transmute(class, gene, baseMean = round(baseMean, 1), log2FC = round(lfc, 2), lfcSE = round(se, 2), padj = signif(padj, 2))), row.names = FALSE)
write_csv(D |> select(class, gene, gene_id, baseMean, log2FC = lfc, lfcSE = se, padj), file.path(OUT, "topoisomerase_tdp_ifnb_vs_mock.csv"))

p <- ggplot(D, aes(xlab, lfc)) + geom_col(fill = FILL, width = 0.7) + geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), width = 0.25, colour = INK, linewidth = 0.35) +
  geom_text(aes(y = ys, label = lab, vjust = vj), colour = INK, size = 3) + geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) +
  facet_grid(. ~ class, scales = "free_x", space = "free_x") + scale_y_continuous(expand = expansion(mult = 0.12)) +
  labs(title = paste0("Topoisomerases and tyrosyl-DNA phosphodiesterases: ", COMP_LAB), subtitle = "shrunk log2FC +/- 1 SE; stars = padj (ns = not significant); dagger = baseMean < 50", x = NULL, y = paste0("log2 fold-change (", COMP_LAB, ")")) +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3), panel.border = element_rect(colour = "grey70"),
        axis.text = element_text(colour = INK), axis.text.x = element_text(face = "italic", angle = 45, hjust = 1), plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 8.5),
        strip.background = element_rect(fill = "grey93", colour = "grey70"), strip.text = element_text(colour = INK, face = "bold"))
w <- max(5, 0.8 * nrow(D) + 2.5)
ggsave(file.path(OUT, "bars_topoisomerase_tdp.png"), p, width = w, height = 4.8, dpi = 300); ggsave(file.path(OUT, "bars_topoisomerase_tdp.pdf"), p, width = w, height = 4.8)
message("\nSaved to: ", OUT)
