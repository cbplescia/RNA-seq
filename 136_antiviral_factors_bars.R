#!/usr/bin/env Rscript
# Bar graphs of differential expression (shrunk log2FC +/- 1 SE, padj stars, dagger = baseMean < 50) for antiviral / egress-restriction genes: Marchf8 (MARCH8), Gbp5, Selplg (PSGL-1), Serinc3, Serinc5, Apobec3, Bst2 (tetherin), Cgas (cGAS; old symbol Mb21d1)
# and ALL protein-coding TRIM genes found in the DE tables (symbols Trim<number>..., ordered by number). One bar per gene, one figure per contrast: IFNB vs Mock and, separately, IFNG (IFNY) vs Mock.
# Notes: cGAMP is a molecule, not a gene (cGAS makes it); SERINC3 and SERINC5 are two genes; mouse has a single Apobec3; MARCH8 is Marchf8 in the current mouse annotation. Genes that are not in the DE tables (removed by the low-count pre-filter, or another symbol) are listed in the console.
# Writes bars_antiviral_IFNB_vs_Mock.png/.pdf, bars_antiviral_IFNY_vs_Mock.png/.pdf and antiviral_factors_log2FC.csv to requant_v2_outputs/figures_gene_bars/.   Env: TRIM_PATTERN (default ^Trim[0-9]), PROJ.   Usage: Rscript 136_antiviral_factors_bars.R

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter; count <- dplyr::count; summarise <- dplyr::summarise; mutate <- dplyr::mutate; arrange <- dplyr::arrange
FACT <- c("Marchf8", "Gbp5", "Selplg", "Serinc3", "Serinc5", "Apobec3", "Bst2", "Cgas"); ALIAS <- list(Marchf8 = "March8", Cgas = "Mb21d1", Selplg = "Psgl1")
TRIM_PATTERN <- Sys.getenv("TRIM_PATTERN", "^Trim[0-9]"); LOW_EXPR <- 50; INK <- "#1F2933"
CON <- list(IFNB_vs_Mock = list(lab = "IFNB vs Mock", col = "#2B6CB0"), IFNY_vs_Mock = list(lab = "IFNG vs Mock", col = "#D97706"))
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs"); OUT <- file.path(RNA, "figures_gene_bars"); dir.create(OUT, showWarnings = FALSE)
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
tn <- function(x) suppressWarnings(as.numeric(sub("^[A-Za-z]+([0-9]+).*$", "\\1", x)))
allr <- list()
for (nm in names(CON)) {
  d <- read_csv(file.path(RNA, paste0("DE_v2_", nm, "_annotated.csv")), show_col_types = FALSE) |> filter(!is.na(gene_name), gene_name != "") |> distinct(gene_id, .keep_all = TRUE); key <- tolower(d$gene_name)
  fac <- bind_rows(lapply(FACT, function(g) { syms <- tolower(c(g, ALIAS[[g]])); h <- d[key %in% syms, ]; if (!nrow(h)) { message("[", CON[[nm]]$lab, "] not in the DE table: ", g); return(NULL) }; h[1, ] |> mutate(query = g) })) |> mutate(class = "Antiviral / restriction factors", ord = match(query, FACT))
  tr <- d |> filter(grepl(TRIM_PATTERN, gene_name), tolower(gene_biotype) == "protein_coding" | is.na(gene_biotype)) |> mutate(query = gene_name, class = "TRIM family", ord = tn(gene_name) * 10 + 1000 + rank(gene_name) / 1000)
  D <- bind_rows(fac, tr) |> transmute(class, gene = gene_name, baseMean, lfc = log2FoldChange, se = lfcSE, padj, ord, contrast = CON[[nm]]$lab)
  D <- D |> mutate(class = factor(class, levels = c("Antiviral / restriction factors", "TRIM family")), lab = stars(padj), vj = ifelse(lfc >= 0, 0, 1), ys = lfc + ifelse(lfc >= 0, 1, -1) * (se + 0.03), xlab = paste0(gene, ifelse(baseMean < LOW_EXPR, "†", ""))) |> arrange(class, ord)
  D$xlab <- factor(D$xlab, levels = unique(D$xlab)); allr[[nm]] <- D; message("[", CON[[nm]]$lab, "] ", sum(D$class == "Antiviral / restriction factors"), " factors and ", sum(D$class == "TRIM family"), " TRIM genes in the DE table")
  p <- ggplot(D, aes(xlab, lfc)) + geom_col(fill = CON[[nm]]$col, width = 0.75) + geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), width = 0.3, colour = INK, linewidth = 0.3) + geom_text(aes(y = ys, label = lab, vjust = vj), colour = INK, size = 2.4) + geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) +
    facet_grid(. ~ class, scales = "free_x", space = "free_x") + scale_y_continuous(expand = expansion(mult = 0.12)) +
    labs(title = paste0("Antiviral / restriction factors and TRIM genes: ", CON[[nm]]$lab), subtitle = "one bar per gene; shrunk log2FC +/- 1 SE; stars = padj (ns = not significant); dagger = baseMean < 50; TRIM genes ordered by number", x = NULL, y = paste0("log2 fold-change (", CON[[nm]]$lab, ")")) +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3), panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK), axis.text.x = element_text(face = "italic", angle = 55, hjust = 1, size = 7.5),
          plot.title = element_text(colour = INK, face = "bold", size = 11), plot.subtitle = element_text(colour = "grey35", size = 8.5), strip.background = element_rect(fill = "grey93", colour = "grey70"), strip.text = element_text(colour = INK, face = "bold", size = 9))
  w <- max(8, 0.3 * nrow(D) + 3); fn <- paste0("bars_antiviral_", nm); ggsave(file.path(OUT, paste0(fn, ".png")), p, width = w, height = 5.6, dpi = 300); ggsave(file.path(OUT, paste0(fn, ".pdf")), p, width = w, height = 5.6)
}
write_csv(bind_rows(allr) |> select(contrast, class, gene, baseMean, log2FC = lfc, lfcSE = se, padj), file.path(OUT, "antiviral_factors_log2FC.csv")); message("\nSaved to: ", OUT)
