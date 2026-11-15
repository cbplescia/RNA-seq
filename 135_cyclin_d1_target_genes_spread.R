#!/usr/bin/env Rscript
# Genes under cyclin D1 control, spread across the genome, and their differential expression (default IFNB, EBOV and IFNB+EBOV vs Mock; bars, shrunk log2FC +/- 1 SE, padj stars, dagger = baseMean < 50).
# Interpretation: cyclin D1 (Ccnd1) + CDK4/6 phosphorylate RB, which releases E2F; so the canonical genes "under cyclin D1 control" are the E2F targets. The script takes the MSigDB mouse Hallmark E2F_TARGETS set (SET=... for another set, e.g.
#   HALLMARK_G2M_CHECKPOINT), keeps the genes that are expressed in your cells (baseMean >= MIN_BASEMEAN, default 100), looks up each gene's chromosome in the GTF and picks the best-expressed gene on each chromosome
#   (N_PER_CHR per chromosome, default 1), so the genes come from all over the genome and not from one cluster. Ccnd1 itself is plotted first as a reference. Bars are ordered by chromosome (chr1 ... chr19, chrX), x labels show the chromosome.
# If instead you meant genes driven by the Ccnd1 PROMOTER (a promoter-reporter, or the TFs that act on it), say so; that is a different list (see 132 / the earlier TF bar graph).
# Needs msigdbr (you installed it for the GSEA).  Env: SET (HALLMARK_E2F_TARGETS), N_PER_CHR (1), MIN_BASEMEAN (100), COMPS (default IFNB_vs_Mock,EBOV_vs_Mock,IFNB_EBOV_vs_Mock), PROJ.
# Writes bars_cyclin_d1_targets_spread.png / .pdf and cyclin_d1_targets_spread.csv (every candidate with chromosome, baseMean and whether it was picked) to requant_v2_outputs/figures_gene_bars/.   Usage: Rscript 135_cyclin_d1_target_genes_spread.R

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter; count <- dplyr::count; summarise <- dplyr::summarise; mutate <- dplyr::mutate; arrange <- dplyr::arrange
if (!requireNamespace("msigdbr", quietly = TRUE)) stop("Install msigdbr:  Rscript -e 'install.packages(\"msigdbr\", repos = \"https://cloud.r-project.org\")'")
SET <- Sys.getenv("SET", "HALLMARK_E2F_TARGETS"); N_PER_CHR <- as.integer(Sys.getenv("N_PER_CHR", "1")); MIN_BASEMEAN <- as.numeric(Sys.getenv("MIN_BASEMEAN", "100")); LOW_EXPR <- 50; INK <- "#1F2933"
ALLC <- c(IFNB_vs_Mock = "IFNB vs Mock", IFNY_vs_Mock = "IFNG vs Mock", EBOV_vs_Mock = "EBOV vs Mock", IFNB_EBOV_vs_Mock = "IFNB+EBOV vs Mock", IFNY_EBOV_vs_Mock = "IFNG+EBOV vs Mock")
ALLCOL <- c(`IFNB vs Mock` = "#2B6CB0", `IFNG vs Mock` = "#D97706", `EBOV vs Mock` = "#1B7837", `IFNB+EBOV vs Mock` = "#9BD08E", `IFNG+EBOV vs Mock` = "#F6C177")
SEL <- trimws(strsplit(Sys.getenv("COMPS", "IFNB_vs_Mock,EBOV_vs_Mock,IFNB_EBOV_vs_Mock"), ",")[[1]]); COMP <- ALLC[SEL]; CCOL <- ALLCOL[COMP]
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs"); OUT <- file.path(RNA, "figures_gene_bars"); dir.create(OUT, showWarnings = FALSE); GTF <- file.path(PROJ, "ref", "Mus_musculus.GRCm39.109.gtf")
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

get_set <- function(nm) { for (arg in c("collection", "category")) { r <- tryCatch(suppressWarnings(do.call(msigdbr::msigdbr, setNames(list("Mus musculus", "H"), c("species", arg)))), error = function(e) NULL); if (!is.null(r) && nrow(r)) return(r |> filter(gs_name == nm) |> pull(gene_symbol) |> unique()) }
  stop("could not load the msigdbr sets") }
g <- get_set(SET); if (!length(g)) stop("set '", SET, "' not found among the Hallmark sets (names like HALLMARK_E2F_TARGETS)"); message(SET, ": ", length(g), " genes")
tabs <- lapply(names(COMP), function(nm) read_csv(file.path(RNA, paste0("DE_v2_", nm, "_annotated.csv")), show_col_types = FALSE) |> mutate(key = tolower(gene_name))); names(tabs) <- names(COMP)
base <- tabs[[1]] |> filter(!is.na(gene_name), tolower(gene_name) %in% tolower(g), baseMean >= MIN_BASEMEAN) |> distinct(gene_id, .keep_all = TRUE); message(nrow(base), " of them expressed (baseMean >= ", MIN_BASEMEAN, ")")
# chromosome from the GTF (one grep pass over the gene lines of the candidates)
pat <- tempfile(); writeLines(c(sprintf('gene_name "%s"', base$gene_name), 'gene_name "Ccnd1"'), pat)
ln <- tryCatch(readLines(pipe(sprintf("grep -F -f '%s' '%s'", pat, GTF))), error = function(e) character(0), warning = function(w) character(0)); f <- strsplit(ln, "\t"); f <- f[vapply(f, length, 1L) >= 9]; f <- f[vapply(f, function(x) x[3] == "gene", TRUE)]
gv <- function(x, k) sub(paste0('.*', k, ' "([^"]+)".*'), "\\1", x)
loc <- tibble(gene_name = vapply(f, function(x) gv(x[9], "gene_name"), ""), chr = vapply(f, function(x) x[1], ""), start = as.integer(vapply(f, function(x) x[4], ""))) |> distinct(gene_name, .keep_all = TRUE)
cand <- base |> inner_join(loc, by = "gene_name") |> filter(chr %in% c(as.character(1:19), "X")) |> mutate(chr_o = ifelse(chr == "X", 20L, as.integer(chr)))
pick <- cand |> group_by(chr) |> slice_max(baseMean, n = N_PER_CHR, with_ties = FALSE) |> ungroup() |> arrange(chr_o, start); message(nrow(pick), " genes picked on ", n_distinct(pick$chr), " chromosomes")
write_csv(cand |> mutate(picked = gene_id %in% pick$gene_id) |> arrange(chr_o, desc(baseMean)) |> select(gene_name, gene_id, chr, start, baseMean, picked), file.path(OUT, "cyclin_d1_targets_spread.csv"))
want <- bind_rows(tibble(class = "Cyclin D1", gene_name = "Ccnd1", chr = loc$chr[tolower(loc$gene_name) == "ccnd1"][1], ord = 0L), pick |> transmute(class = paste0(SET, " (one per chromosome)"), gene_name, chr, ord = chr_o * 1000L + row_number()))
cls <- unique(want$class); want <- want |> mutate(class = factor(class, levels = cls), pos = row_number())
D <- bind_rows(lapply(names(COMP), function(nm) tabs[[nm]] |> inner_join(want |> mutate(key = tolower(gene_name)) |> select(key, class, chr, pos), by = "key") |> transmute(class, gene = gene_name, chr, pos, baseMean, comparison = COMP[[nm]], lfc = log2FoldChange, se = lfcSE, padj)))
D <- D |> mutate(comparison = factor(comparison, levels = COMP), lab = stars(padj), vj = ifelse(lfc >= 0, 0, 1), ys = lfc + ifelse(lfc >= 0, 1, -1) * (se + 0.03), xlab = paste0(gene, ifelse(baseMean < LOW_EXPR, "†", ""), "\nchr", chr))
D$xlab <- factor(D$xlab, levels = unique(D$xlab[order(D$pos)])); print(as.data.frame(D |> transmute(class, gene, chr, baseMean = round(baseMean, 1), comparison, log2FC = round(lfc, 2), padj = signif(padj, 2)) |> pivot_wider(names_from = comparison, values_from = c(log2FC, padj))), row.names = FALSE)
dg <- position_dodge(width = 0.85)
p <- ggplot(D, aes(xlab, lfc, fill = comparison)) + geom_col(position = dg, width = 0.8) + geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), position = dg, width = 0.25, colour = INK, linewidth = 0.35) + geom_text(aes(y = ys, label = lab, vjust = vj), position = dg, colour = INK, size = 2.5) +
  geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) + facet_grid(. ~ class, scales = "free_x", space = "free_x", labeller = labeller(class = label_wrap_gen(22))) + scale_fill_manual(values = CCOL, name = NULL) + scale_y_continuous(expand = expansion(mult = 0.12)) +
  labs(title = "Cyclin D1 and E2F target genes spread across the genome", subtitle = paste0(SET, ": best-expressed gene on each chromosome; shrunk log2FC +/- 1 SE; stars = padj (ns = not significant); dagger = baseMean < 50; label shows the chromosome"), x = NULL, y = "log2 fold-change vs Mock") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3), panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK), axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.title = element_text(colour = INK, face = "bold", size = 11), plot.subtitle = element_text(colour = "grey35", size = 8), legend.position = "top", strip.background = element_rect(fill = "grey93", colour = "grey70"), strip.text = element_text(colour = INK, face = "bold", size = 8))
w <- max(8, 0.34 * nrow(D) + 3); ggsave(file.path(OUT, "bars_cyclin_d1_targets_spread.png"), p, width = w, height = 5.6, dpi = 300); ggsave(file.path(OUT, "bars_cyclin_d1_targets_spread.pdf"), p, width = w, height = 5.6)
message("\nSaved to: ", OUT)
