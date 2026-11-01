#!/usr/bin/env Rscript
# Box-and-whisker plot of VST expression of one or more genes in any number of groups (default Mock, EBOV, MARV), every replicate drawn as a point over the box.
# Above each non-reference group: the DESeq2 result for that group vs Mock (shrunk log2FC and BH-adjusted padj stars) read from DE_v2_<group>_vs_Mock_annotated.csv; the first group (Mock_Mock) is the reference.
# Env: GENE (comma list, matched by symbol ignoring case; default Top1), GROUPS (comma list of group names in the data, default Mock_Mock,Mock_EBOV,Mock_MARV; e.g. add IFNB_Mock,IFNB_EBOV), PROJ.
# Group names: Mock_Mock, IFNB_Mock, IFNY_Mock, Mock_EBOV, IFNB_EBOV, IFNY_EBOV, Mock_MARV, IFNB_MARV, IFNY_MARV, Mock_EVLP, Mock_MVLP. Box = median and quartiles, whiskers = 1.5 x IQR (with 3 replicates the box is simply the spread of the three values).
# Writes boxplot_vst_<genes>_<groups>.png / .pdf and a .csv of the VST values to requant_v2_outputs/figures_gene_bars/.   Usage: GENE=Top1 Rscript 123_vst_boxplot_groups.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2); library(DESeq2) })
select <- dplyr::select; filter <- dplyr::filter

GENES <- trimws(strsplit(Sys.getenv("GENE", unset = "Top1"), ",")[[1]])
GRP <- trimws(strsplit(Sys.getenv("GROUPS", unset = "Mock_Mock,Mock_EBOV,Mock_MARV"), ",")[[1]])
LAB <- c(Mock_Mock = "Mock", Mock_EBOV = "EBOV", IFNB_Mock = "IFNB", IFNB_EBOV = "IFNB+EBOV", IFNY_Mock = "IFNG", IFNY_EBOV = "IFNG+EBOV", Mock_MARV = "MARV", IFNB_MARV = "IFNB+MARV", IFNY_MARV = "IFNG+MARV", Mock_EVLP = "EVLP", Mock_MVLP = "MVLP")
COL <- c(Mock = "#7F7F7F", EBOV = "#1B7837", IFNB = "#2B6CB0", `IFNB+EBOV` = "#9BD08E", IFNG = "#D69E2E", `IFNG+EBOV` = "#E8C77A", MARV = "#8C510A", `IFNB+MARV` = "#C7A27C", `IFNG+MARV` = "#DCC9A8", EVLP = "#6A51A3", MVLP = "#CBA4D8")
INK <- "#1F2933"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs"); OUT <- file.path(RNA, "figures_gene_bars"); dir.create(OUT, showWarnings = FALSE)
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
gl <- function(g) ifelse(g %in% names(LAB), LAB[g], g)

vsd <- readRDS(file.path(RNA, "vsd_v2.rds")); vm <- SummarizedExperiment::assay(vsd); grp <- as.character(vsd$group)
if (!all(GRP %in% grp)) stop("Groups not in the data: ", paste(setdiff(GRP, grp), collapse = ", "), ". Available: ", paste(sort(unique(grp)), collapse = ", "))
cn <- function(g) { p <- strsplit(g, "_")[[1]]; if (p[1] == "Mock") paste0(p[2], "_vs_Mock") else if (p[2] == "Mock") paste0(p[1], "_vs_Mock") else paste0(p[1], "_", p[2], "_vs_Mock") }
stat_tab <- if (GRP[1] == "Mock_Mock") bind_rows(lapply(GRP[-1], function(g) { f <- file.path(RNA, paste0("DE_v2_", cn(g), "_annotated.csv")); if (!file.exists(f)) { message("No DE table ", basename(f), " - no stat for ", g); return(NULL) }
  read_csv(f, show_col_types = FALSE) |> transmute(gene_id, code = g, baseMean, lfc = log2FoldChange, padj) })) else NULL
allg <- read_csv(file.path(RNA, "DE_v2_EBOV_vs_Mock_annotated.csv"), show_col_types = FALSE) |> mutate(key = tolower(gene_name), gid = sub("\\..*$", "", gene_id))

rows <- list()
for (g in GENES) {
  h <- allg |> filter(key == tolower(g) | tolower(gid) == tolower(g)) |> distinct(gene_id, .keep_all = TRUE); h <- h[h$gene_id %in% rownames(vm), , drop = FALSE]
  if (!nrow(h)) { message("'", g, "' is not in the VST object (not found, or removed by the pre-filter in 04)"); next }
  for (k in seq_len(nrow(h))) { smp <- colnames(vm)[grp %in% GRP]; cd <- grp[match(smp, colnames(vm))]
    rows[[length(rows) + 1]] <- tibble(gene = h$gene_name[k], gene_id = h$gene_id[k], sample = smp, code = cd, group = factor(gl(cd), levels = gl(GRP)), vst = as.numeric(vm[h$gene_id[k], smp])) }
}
if (!length(rows)) stop("None of the genes were found")
D <- bind_rows(rows)
ann <- D |> group_by(gene, gene_id, code, group) |> summarise(top = max(vst), .groups = "drop")
if (!is.null(stat_tab)) ann <- ann |> left_join(stat_tab, by = c("gene_id", "code")) |> filter(code != GRP[1]) |> mutate(txt = sprintf("log2FC %+.2f\n%s", lfc, stars(padj))) else ann <- ann |> mutate(txt = "")
write_csv(D |> arrange(gene, group), file.path(OUT, paste0("boxplot_vst_", paste(GENES, collapse = "_"), "_", paste(GRP, collapse = "_"), ".csv")))
print(as.data.frame(D |> group_by(gene, group) |> summarise(n = n(), mean_vst = round(mean(vst), 2), min = round(min(vst), 2), max = round(max(vst), 2), .groups = "drop")), row.names = FALSE)
if (!is.null(stat_tab)) print(as.data.frame(ann |> transmute(gene, group, baseMean = round(baseMean, 1), log2FC = round(lfc, 2), padj = signif(padj, 2))), row.names = FALSE)

cols <- COL[gl(GRP)]; names(cols) <- gl(GRP); yr <- range(D$vst); off <- 0.04 * max(diff(yr), 0.5)
p <- ggplot(D, aes(group, vst, fill = group)) + geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.8, colour = INK, linewidth = 0.4) +
  geom_point(position = position_jitter(width = 0.08, height = 0, seed = 1), size = 2.2, shape = 21, colour = INK, fill = "white", stroke = 0.6) +
  geom_text(data = ann, aes(x = group, y = top + off, label = txt), inherit.aes = FALSE, colour = INK, size = 3.2, lineheight = 0.95, vjust = 0) +
  facet_wrap(~gene, scales = "free_y", nrow = 1) + scale_fill_manual(values = cols, guide = "none") + scale_y_continuous(expand = expansion(mult = c(0.08, 0.3))) +
  labs(title = sprintf("%s: VST expression, %s", paste(unique(D$gene), collapse = ", "), paste(gl(GRP), collapse = ", ")), subtitle = "box = median and quartiles of the replicates; points = individual samples; stat = DESeq2 vs Mock (shrunk log2FC, BH padj; ns = not significant)",
       x = NULL, y = "VST expression (log2 scale)") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK), axis.title = element_text(colour = INK), plot.title = element_text(colour = INK, face = "bold"),
        plot.subtitle = element_text(colour = "grey35", size = 8.5), strip.background = element_rect(fill = "grey95", colour = "grey70"), strip.text = element_text(colour = INK, face = "bold.italic"))
nm <- paste0("boxplot_vst_", paste(unique(D$gene), collapse = "_"), "_", paste(gl(GRP), collapse = "_"))
w <- 1.8 + 0.9 * length(GRP) * n_distinct(D$gene)
ggsave(file.path(OUT, paste0(nm, ".png")), p, width = w, height = 5, dpi = 300); ggsave(file.path(OUT, paste0(nm, ".pdf")), p, width = w, height = 5)
message("Saved to: ", OUT)
