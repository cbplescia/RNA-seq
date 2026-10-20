#!/usr/bin/env Rscript
# PCA of six conditions, chosen by sample number:
#   Mock 1-3 | EBOV 4-6 | MARV 7-9 | IFNB 10-12 | IFNB+EBOV 16-18 | IFNB+MARV 22-24
# Dots are colour-coded by group and labelled with the REPLICATE number (1, 2, 3) inside the dot, not the sample ID.
# IFNB-treated groups are squares, the others circles; a light shaded triangle joins the three replicates of each group.
# PCA = the 500 most variable protein-coding genes (VST, vsd_v2.rds) across these 18 samples, centred, via prcomp.
# A check table is printed: it compares the groups you chose with the group each sample has in sample_info.csv.
#
# Writes to <project>/requant_v2_outputs/figures_pca/ (PNG + PDF):
#   PCA_six_groups_PC1_PC2, PCA_six_groups_PC1_PC3, PCA_six_groups_coordinates.csv, PCA_six_groups_top_loadings.csv
# Usage: Rscript 46_pca_six_groups.R   (PROJ comes from 'source setup_session.sh')

suppressPackageStartupMessages({ library(DESeq2); library(readr); library(dplyr); library(tidyr); library(ggplot2) })

NTOP <- 500
SETS <- list(Mock = c("1", "2", "3"), EBOV = c("4", "5", "6"), MARV = c("7", "8", "9"),
             IFNB = c("10", "11", "12"), `IFNB+EBOV` = c("16", "17", "18"), `IFNB+MARV` = c("22", "23", "24"))
EXPECTED <- c(Mock = "Mock_Mock", EBOV = "Mock_EBOV", MARV = "Mock_MARV", IFNB = "IFNB_Mock",
              `IFNB+EBOV` = "IFNB_EBOV", `IFNB+MARV` = "IFNB_MARV")
COLS <- c(Mock = "#7F7F7F", EBOV = "#E69F00", MARV = "#009E73", IFNB = "#0072B2", `IFNB+EBOV` = "#D55E00", `IFNB+MARV` = "#CC79A7")
INK <- "#1F2933"

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_pca"); dir.create(OUT, showWarnings = FALSE)

vsd <- readRDS(file.path(RNA, "vsd_v2.rds"))
vm <- assay(vsd)
sel <- tibble(sample = unlist(SETS), treatment = rep(names(SETS), lengths(SETS)),
              replicate = unlist(lapply(SETS, seq_along)))
sel$data_group <- as.character(vsd$group)[match(sel$sample, colnames(vm))]
sel$expected <- EXPECTED[sel$treatment]
if (any(is.na(sel$data_group))) stop("Samples not found in the VST object: ", paste(sel$sample[is.na(sel$data_group)], collapse = ", "),
                                    "\nSample IDs available: ", paste(colnames(vm), collapse = ", "))
message("Check - chosen group vs the group recorded for that sample in the data:")
print(as.data.frame(sel |> transmute(sample, chosen_as = treatment, replicate, group_in_data = data_group,
                                     match = ifelse(data_group == expected, "ok", "DIFFERENT"))), row.names = FALSE)
if (any(sel$data_group != sel$expected)) message("WARNING: some samples are recorded under a different group than you specified. Using YOUR selection.")

pcf <- file.path(RNA, "DE_v2_IFNB_vs_Mock_protein_coding.csv")
sym <- NULL
if (file.exists(pcf)) {
  pc <- read_csv(pcf, show_col_types = FALSE)
  vm <- vm[rownames(vm) %in% pc$gene_id, , drop = FALSE]; sym <- setNames(pc$gene_name, pc$gene_id)
} else message("Protein-coding table not found - using every gene in the VST object")
x <- vm[, sel$sample, drop = FALSE]
v <- apply(x, 1, var); x <- x[order(v, decreasing = TRUE)[seq_len(min(NTOP, sum(v > 0)))], , drop = FALSE]
p <- prcomp(t(x), center = TRUE, scale. = FALSE)
pct <- 100 * p$sdev^2 / sum(p$sdev^2)
message(sprintf("\n%d samples, %d genes. Variance explained: PC1 %.1f%%, PC2 %.1f%%, PC3 %.1f%%, PC4 %.1f%%", nrow(sel), nrow(x), pct[1], pct[2], pct[3], pct[4]))

co <- as.data.frame(p$x[, 1:5]) |> tibble::rownames_to_column("sample") |> left_join(sel |> select(sample, treatment, replicate), by = "sample") |>
  mutate(treatment = factor(treatment, levels = names(SETS)), ifn = factor(ifelse(grepl("IFNB", treatment), "IFNB-treated", "no IFN"), levels = c("no IFN", "IFNB-treated")))
write_csv(co, file.path(OUT, "PCA_six_groups_coordinates.csv"))
ld <- as.data.frame(p$rotation[, 1:3]) |> tibble::rownames_to_column("gene_id") |> mutate(gene_name = if (is.null(sym)) gene_id else unname(sym[gene_id]))
write_csv(bind_rows(lapply(c("PC1", "PC2", "PC3"), function(pcn) ld |> slice_max(abs(.data[[pcn]]), n = 20) |> mutate(driver_of = pcn))),
          file.path(OUT, "PCA_six_groups_top_loadings.csv"))

plot_pc <- function(px, py) {
  hull <- co |> group_by(treatment) |> slice(chull(.data[[px]], .data[[py]])) |> ungroup()
  ix <- as.integer(sub("PC", "", px)); iy <- as.integer(sub("PC", "", py))
  ggplot(co, aes(.data[[px]], .data[[py]])) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) + geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_polygon(data = hull, aes(fill = treatment, group = treatment), alpha = 0.13, colour = NA) +
    geom_point(aes(fill = treatment, shape = ifn), size = 8, colour = "white", stroke = 0.9) +
    geom_text(aes(label = replicate), colour = "white", size = 3.6, fontface = "bold") +
    scale_fill_manual(values = COLS, name = "group") +
    scale_shape_manual(values = c(`no IFN` = 21, `IFNB-treated` = 22), name = NULL) +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 5, alpha = 1), order = 1),
           shape = guide_legend(override.aes = list(fill = "grey60", size = 5), order = 2)) +
    scale_x_continuous(expand = expansion(mult = 0.12)) + scale_y_continuous(expand = expansion(mult = 0.12)) +
    labs(title = "PCA: Mock, EBOV, MARV, IFNB, IFNB+EBOV and IFNB+MARV",
         subtitle = sprintf("top %d most variable protein-coding genes (VST), 18 samples; number inside each dot = replicate", nrow(x)),
         x = sprintf("%s (%.1f%% of variance)", px, pct[ix]), y = sprintf("%s (%.1f%% of variance)", py, pct[iy])) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK),
          plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 9),
          legend.position = "right")
}
for (pp in list(c("PC1", "PC2"), c("PC1", "PC3"))) {
  g <- plot_pc(pp[1], pp[2]); stub <- paste0("PCA_six_groups_", pp[1], "_", pp[2])
  ggsave(file.path(OUT, paste0(stub, ".png")), g, width = 9, height = 6.8, dpi = 300)
  ggsave(file.path(OUT, paste0(stub, ".pdf")), g, width = 9, height = 6.8)
}
message("\nFigures and tables in: ", OUT)
