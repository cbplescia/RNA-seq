#!/usr/bin/env Rscript
# PCA of the RNA-seq samples.
#   1. Mock, IFNB and IFNy samples (9 samples: the interferon-only groups, no virus)
#   2. Mock and IFNB samples only (6 samples)
# Each PCA is recomputed on its own samples: the 500 most variable protein-coding genes (variance of the VST values
# across THAT set of samples) are centred and rotated with prcomp. Uses vsd_v2.rds from 04_tximport_deseq2.R.
# Change NTOP below to use more or fewer genes (NTOP <- Inf = all protein-coding genes).
#
# Writes to <project>/requant_v2_outputs/figures_pca/ (PNG + PDF + CSV):
#   PCA_Mock_IFNB_IFNY, PCA_Mock_IFNB            the plots (PC1 vs PC2, % variance on the axes, samples labelled)
#   PCA_<set>_coordinates.csv                    PC1..PC5 of every sample
#   PCA_<set>_top_loadings.csv                   the genes that drive PC1 and PC2 most
# Usage: Rscript 45_pca_mock_ifnb_ifny.R   (PROJ is read from the environment; run 'source setup_session.sh' first)

suppressPackageStartupMessages({ library(DESeq2); library(readr); library(dplyr); library(tidyr); library(ggplot2) })
if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Install ggrepel (conda install -c conda-forge r-ggrepel)")

NTOP <- 500
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_pca"); dir.create(OUT, showWarnings = FALSE)
INK <- "#1F2933"
COLS <- c(Mock = "#8A8A8A", IFNB = "#2B6CB0", `IFNγ` = "#D69E2E")

vsd <- readRDS(file.path(RNA, "vsd_v2.rds"))
vm <- assay(vsd)
meta <- tibble(sample = colnames(vm), group = as.character(vsd$group))
pcf <- file.path(RNA, "DE_v2_IFNB_vs_Mock_protein_coding.csv")
if (file.exists(pcf)) {
  pc <- read_csv(pcf, show_col_types = FALSE)
  vm <- vm[rownames(vm) %in% pc$gene_id, , drop = FALSE]; sym <- setNames(pc$gene_name, pc$gene_id)
} else { message("Protein-coding table not found - using every gene in the VST object"); sym <- setNames(rownames(vm), rownames(vm)) }
message("Genes available for PCA: ", nrow(vm))

label_of <- c(Mock_Mock = "Mock", IFNB_Mock = "IFNB", IFNY_Mock = "IFNγ")
run_pca <- function(groups, file_stub, title) {
  m <- meta |> filter(group %in% groups)
  if (!nrow(m)) stop("No samples found for ", paste(groups, collapse = ", "))
  x <- vm[, m$sample, drop = FALSE]
  v <- apply(x, 1, var); x <- x[order(v, decreasing = TRUE)[seq_len(min(NTOP, sum(v > 0)))], , drop = FALSE]
  p <- prcomp(t(x), center = TRUE, scale. = FALSE)
  pct <- 100 * p$sdev^2 / sum(p$sdev^2)
  co <- as.data.frame(p$x[, seq_len(min(5, ncol(p$x))), drop = FALSE]) |> tibble::rownames_to_column("sample") |>
    left_join(m, by = "sample") |> mutate(treatment = factor(label_of[group], levels = c("Mock", "IFNB", "IFNγ")))
  write_csv(co, file.path(OUT, paste0("PCA_", file_stub, "_coordinates.csv")))
  ld <- as.data.frame(p$rotation[, 1:2, drop = FALSE]) |> tibble::rownames_to_column("gene_id") |>
    mutate(gene_name = unname(sym[gene_id]))
  write_csv(bind_rows(ld |> slice_max(abs(PC1), n = 20) |> mutate(driver_of = "PC1"),
                      ld |> slice_max(abs(PC2), n = 20) |> mutate(driver_of = "PC2")),
            file.path(OUT, paste0("PCA_", file_stub, "_top_loadings.csv")))
  message(sprintf("\n%s: %d samples, %d genes; variance explained PC1 %.1f%%, PC2 %.1f%%, PC3 %.1f%%", title, nrow(m), nrow(x), pct[1], pct[2], pct[3]))
  print(as.data.frame(co |> select(sample, treatment, PC1, PC2) |> mutate(across(c(PC1, PC2), ~ round(.x, 1)))), row.names = FALSE)
  g <- ggplot(co, aes(PC1, PC2, colour = treatment)) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) + geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_point(size = 4.2, alpha = 0.9) +
    ggrepel::geom_text_repel(aes(label = sample), size = 3, colour = INK, show.legend = FALSE, box.padding = 0.4, seed = 1) +
    scale_colour_manual(values = COLS, breaks = levels(droplevels(co$treatment)), name = NULL) +
    scale_x_continuous(expand = expansion(mult = 0.15)) + scale_y_continuous(expand = expansion(mult = 0.15)) +
    labs(title = title, subtitle = sprintf("top %d most variable protein-coding genes (VST), %d samples; each dot is one sample", nrow(x), nrow(m)),
         x = sprintf("PC1 (%.1f%% of variance)", pct[1]), y = sprintf("PC2 (%.1f%% of variance)", pct[2])) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK),
          plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 9),
          legend.position = "bottom")
  ggsave(file.path(OUT, paste0("PCA_", file_stub, ".png")), g, width = 7.5, height = 6.5, dpi = 300)
  ggsave(file.path(OUT, paste0("PCA_", file_stub, ".pdf")), g, width = 7.5, height = 6.5)
}
run_pca(c("Mock_Mock", "IFNB_Mock", "IFNY_Mock"), "Mock_IFNB_IFNY", "PCA: Mock, IFNB and IFNγ")
run_pca(c("Mock_Mock", "IFNB_Mock"), "Mock_IFNB", "PCA: Mock and IFNB")
message("\nFigures and tables in: ", OUT)
