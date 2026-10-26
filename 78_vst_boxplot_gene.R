#!/usr/bin/env Rscript
# Box plot of VST expression of one or more genes in Mock (untreated) and EBOV samples; every replicate is drawn as a point over the box.
# The stat shown above the boxes is the DESeq2 result for that gene from DE_v2_EBOV_vs_Mock_annotated.csv (shrunk log2FC and padj, Wald test BH-adjusted over all genes).
# Genes are matched by symbol ignoring case (default Npc1); comma-separate several to get one panel each.
# Change the groups with GROUPS (names as in the data, e.g. Mock_Mock, Mock_EBOV, IFNB_Mock, IFNB_EBOV).
# Writes boxplot_vst_<gene>_<groupA>_vs_<groupB>.png / .pdf and a .csv of the VST values to requant_v2_outputs/figures_gene_bars/
# Usage: GENE=Npc1 Rscript 78_vst_boxplot_gene.R     (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2); library(DESeq2) })
select <- dplyr::select; filter <- dplyr::filter

GENES <- trimws(strsplit(Sys.getenv("GENE", unset = "Npc1"), ",")[[1]])
GRP <- trimws(strsplit(Sys.getenv("GROUPS", unset = "Mock_Mock,Mock_EBOV"), ",")[[1]])
LAB <- c(Mock_Mock = "Mock", Mock_EBOV = "EBOV", IFNB_Mock = "IFNB", IFNB_EBOV = "IFNB+EBOV", IFNY_Mock = "IFNG", IFNY_EBOV = "IFNG+EBOV", Mock_MARV = "MARV", IFNB_MARV = "IFNB+MARV", IFNY_MARV = "IFNG+MARV")
COL <- c(Mock = "#7F7F7F", EBOV = "#1B7837", IFNB = "#2B6CB0", `IFNB+EBOV` = "#9BD08E", IFNG = "#D69E2E", `IFNG+EBOV` = "#E8C77A", MARV = "#8C510A", `IFNB+MARV` = "#C7A27C", `IFNG+MARV` = "#DCC9A8")
INK <- "#1F2933"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_gene_bars"); dir.create(OUT, showWarnings = FALSE)
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

vsd <- readRDS(file.path(RNA, "vsd_v2.rds")); vm <- SummarizedExperiment::assay(vsd); grp <- as.character(vsd$group)
if (!all(GRP %in% grp)) stop("Groups not in the data: ", paste(setdiff(GRP, grp), collapse = ", "), ". Available: ", paste(sort(unique(grp)), collapse = ", "))
gl <- function(g) ifelse(g %in% names(LAB), LAB[g], g)
# results table for the pair (the file name follows the DE_v2_* naming; EBOV vs Mock is the default pair)
pair_file <- if (identical(GRP, c("Mock_Mock", "Mock_EBOV"))) "DE_v2_EBOV_vs_Mock_annotated.csv" else NULL
de <- if (!is.null(pair_file) && file.exists(file.path(RNA, pair_file))) read_csv(file.path(RNA, pair_file), show_col_types = FALSE) |> mutate(key = tolower(gene_name)) else NULL
if (is.null(de)) message("No DE table is attached to this group pair - plotting without a stat (the default Mock vs EBOV pair uses DE_v2_EBOV_vs_Mock_annotated.csv)")
allg <- read_csv(file.path(RNA, "DE_v2_EBOV_vs_Mock_annotated.csv"), show_col_types = FALSE) |> mutate(key = tolower(gene_name), gid = sub("\\..*$", "", gene_id))

rows <- list()
for (g in GENES) {
  h <- allg |> filter(key == tolower(g) | tolower(gid) == tolower(g)) |> distinct(gene_id, .keep_all = TRUE)
  h <- h[h$gene_id %in% rownames(vm), , drop = FALSE]
  if (!nrow(h)) { message("'", g, "' is not in the VST object (not found, or removed by the pre-filter in 04)"); next }
  for (k in seq_len(nrow(h))) {
    smp <- colnames(vm)[grp %in% GRP]
    rows[[length(rows) + 1]] <- tibble(gene = h$gene_name[k], gene_id = h$gene_id[k], sample = smp, group = factor(gl(grp[match(smp, colnames(vm))]), levels = gl(GRP)), vst = as.numeric(vm[h$gene_id[k], smp]))
  }
}
if (!length(rows)) stop("None of the genes were found")
D <- bind_rows(rows)
ann <- D |> distinct(gene, gene_id) |> left_join(if (!is.null(de)) de |> select(gene_id, baseMean, log2FoldChange, padj) else tibble(gene_id = character()), by = "gene_id") |>
  left_join(D |> group_by(gene) |> summarise(top = max(vst), .groups = "drop"), by = "gene") |>
  mutate(txt = if (!is.null(de)) sprintf("log2FC = %.2f\npadj = %s  %s", log2FoldChange, formatC(padj, digits = 2, format = "g"), stars(padj)) else "")
write_csv(D |> arrange(gene, group), file.path(OUT, paste0("boxplot_vst_", paste(GENES, collapse = "_"), "_", paste(GRP, collapse = "_vs_"), ".csv")))
print(as.data.frame(D |> group_by(gene, group) |> summarise(n = n(), mean_vst = round(mean(vst), 2), min = round(min(vst), 2), max = round(max(vst), 2), .groups = "drop")), row.names = FALSE)
if (!is.null(de)) print(as.data.frame(ann |> select(gene, baseMean, log2FoldChange, padj) |> mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)

cols <- COL[gl(GRP)]; names(cols) <- gl(GRP)
yr <- range(D$vst); pad <- 0.35 * diff(yr) + 0.15
p <- ggplot(D, aes(group, vst, fill = group)) + geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.8, colour = INK, linewidth = 0.4) +
  geom_point(position = position_jitter(width = 0.08, height = 0, seed = 1), size = 2.2, shape = 21, colour = INK, fill = "white", stroke = 0.6) +
  geom_text(data = ann, aes(x = 1.5, y = top + 0.12 * pad, label = txt), inherit.aes = FALSE, colour = INK, size = 3.4, lineheight = 0.95, vjust = 0) +
  facet_wrap(~gene, scales = "free_y", nrow = 1) + scale_fill_manual(values = cols, guide = "none") + scale_y_continuous(expand = expansion(mult = c(0.08, 0.3))) +
  labs(title = sprintf("%s: VST expression, %s vs %s", paste(unique(D$gene), collapse = ", "), gl(GRP[1]), gl(GRP[2])), subtitle = "box = median and quartiles of the replicates; points = individual samples; stat = DESeq2 (shrunk log2FC, BH padj)",
       x = NULL, y = "VST expression (log2 scale)") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK), axis.title = element_text(colour = INK), plot.title = element_text(colour = INK, face = "bold"),
        plot.subtitle = element_text(colour = "grey35", size = 8.5), strip.background = element_rect(fill = "grey95", colour = "grey70"), strip.text = element_text(colour = INK, face = "bold"))
nm <- paste0("boxplot_vst_", paste(unique(D$gene), collapse = "_"), "_", paste(gl(GRP), collapse = "_vs_"))
ggsave(file.path(OUT, paste0(nm, ".png")), p, width = 2.6 + 2.4 * n_distinct(D$gene), height = 5, dpi = 300)
ggsave(file.path(OUT, paste0(nm, ".pdf")), p, width = 2.6 + 2.4 * n_distinct(D$gene), height = 5)
message("Saved to: ", OUT)
