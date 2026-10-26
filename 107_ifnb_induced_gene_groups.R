#!/usr/bin/env Rscript
# Which strongly IFNB-induced genes does EBOV bring back to baseline, restrict partly, or leave alone?   (all conditions vs Mock; protein-coding genes)
# Uses figures_mock_ifnb_ebov/transcriptome_blocks.csv written by 51_mock_ifnb_ebov_analysis.R (run 51 first): log2FC and padj of IFNB vs Mock, IFNB+EBOV vs Mock and IFNB+EBOV vs IFNB.
#   highly induced = IFNB vs Mock padj < PADJ_CUT and log2FC >= HIGH (default 1)
#   retention      = log2FC(IFNB+EBOV vs Mock) / log2FC(IFNB vs Mock): 1 = induction fully kept, 0 = back to baseline
#   FULLY SUPPRESSED   retention <= R_FULL (0.25) and IFNB+EBOV is significantly lower than IFNB (padj < PADJ_CUT)         -> "back to baseline"
#   PARTLY RESTRICTED  R_FULL < retention < R_NOT (0.75)  and IFNB+EBOV lower than IFNB                                  -> "partly restricted"
#   NOT SUPPRESSED     retention >= R_NOT (includes genes EBOV raises further)                                              -> "not suppressed"
#   UNCLEAR            retention <= R_FULL / < R_NOT but the IFNB+EBOV vs IFNB difference is not significant (noisy)          -> listed, not drawn
# Writes to requant_v2_outputs/figures_ifnb_ebov_retention/ (PNG + PDF):
#   scatter_IFNB_vs_IFNB_EBOV          every highly induced gene: x = log2FC IFNB, y = log2FC IFNB+EBOV, coloured by group, guide lines at retention 0, 0.25, 0.75, 1; top genes labelled
#   profiles_by_group                  per group: each gene's log2FC at Mock (0) -> IFNB -> IFNB+EBOV (thin lines) and the group mean (thick line)
#   heatmap_groups                     all grouped genes x every replicate of Untreated | IFNB | IFNB+EBOV (each replicate vs the Untreated mean, Mock = grey), rows in group blocks sorted by IFNB log2FC
#   group_counts                       how many genes per group
#   ifnb_induced_groups.csv            every gene with its group, log2FC / padj values and retention (the input of 108_promoter_tf_groups.py)
# Usage: Rscript 107_ifnb_induced_gene_groups.R   (PROJ from 'source setup_session.sh';  HIGH=1.5  R_FULL=0.25  R_NOT=0.75 optional)

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2); library(ggrepel); library(pheatmap); library(grid) })
select <- dplyr::select; filter <- dplyr::filter

PADJ_CUT <- 0.05; HIGH <- as.numeric(Sys.getenv("HIGH", unset = "1")); R_FULL <- as.numeric(Sys.getenv("R_FULL", unset = "0.25")); R_NOT <- as.numeric(Sys.getenv("R_NOT", unset = "0.75")); N_LAB <- 8; CAP_Q <- 0.98
INK <- "#1F2933"; GC <- c(`fully suppressed` = "#762A83", `partly restricted` = "#E08214", `not suppressed` = "#1B7837")
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_ifnb_ebov_retention"); dir.create(OUT, showWarnings = FALSE)
save2 <- function(p, name, w, h) { ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300); ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h) }
bf <- file.path(RNA, "figures_mock_ifnb_ebov", "transcriptome_blocks.csv"); if (!file.exists(bf)) stop("Run 51_mock_ifnb_ebov_analysis.R first (missing ", bf, ")")
B <- read_csv(bf, show_col_types = FALSE)

H <- B |> filter(!is.na(lfc_B), !is.na(padj_B), padj_B < PADJ_CUT, lfc_B >= HIGH, !is.na(lfc_BE)) |>
  mutate(retention = lfc_BE / lfc_B, lower = !is.na(padj_BE_B) & padj_BE_B < PADJ_CUT & !is.na(lfc_BE_B) & lfc_BE_B < 0,
         group = case_when(retention >= R_NOT ~ "not suppressed", retention <= R_FULL & lower ~ "fully suppressed", retention < R_NOT & retention > R_FULL & lower ~ "partly restricted", TRUE ~ "unclear"))
tab <- table(factor(H$group, levels = c(names(GC), "unclear")))
message(sprintf("highly induced genes (IFNB padj < %g, log2FC >= %g): %d | fully suppressed %d | partly restricted %d | not suppressed %d | unclear %d", PADJ_CUT, HIGH, nrow(H), tab[[1]], tab[[2]], tab[[3]], tab[[4]]))
write_csv(H |> arrange(factor(group, levels = c(names(GC), "unclear")), desc(lfc_B)) |> select(gene_id, gene_name, group, retention, lfc_B, padj_B, lfc_BE, padj_BE, lfc_BE_B, padj_BE_B), file.path(OUT, "ifnb_induced_groups.csv"))
G <- H |> filter(group %in% names(GC)) |> mutate(group = factor(group, levels = names(GC)), label = ifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name))
if (!nrow(G)) stop("No gene falls in a group at these cut-offs")
base_theme <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey94", linewidth = 0.3), panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK),
  axis.title = element_text(colour = INK), plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 8.5), strip.background = element_rect(fill = "grey95", colour = "grey70"), strip.text = element_text(colour = INK, face = "bold"))

# ---- counts ----
pc <- ggplot(G |> count(group), aes(group, n, fill = group)) + geom_col(width = 0.65) + geom_text(aes(label = n), vjust = -0.3) + scale_fill_manual(values = GC, guide = "none") + scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Strongly IFNB-induced genes by how EBOV affects them", subtitle = sprintf("IFNB vs Mock log2FC >= %g, padj < %g; %d genes unclear (not drawn)", HIGH, PADJ_CUT, tab[["unclear"]]), x = NULL, y = "genes") + base_theme
save2(pc, "group_counts", 5.5, 4)

# ---- scatter ----
lab <- G |> group_by(group) |> slice_max(lfc_B, n = N_LAB, with_ties = FALSE) |> ungroup(); lim <- max(G$lfc_B, na.rm = TRUE) * 1.05
ps <- ggplot(G, aes(lfc_B, lfc_BE, colour = group)) + geom_abline(slope = c(0, R_FULL, R_NOT, 1), intercept = 0, linetype = c("solid", "dashed", "dashed", "dotted"), colour = "grey55", linewidth = 0.4) + geom_point(size = 1.4, alpha = 0.85) +
  geom_text_repel(data = lab, aes(label = label), size = 2.8, fontface = "italic", max.overlaps = Inf, min.segment.length = 0, segment.size = 0.2, segment.colour = "grey60", box.padding = 0.3, seed = 1, show.legend = FALSE) +
  scale_colour_manual(values = GC, name = NULL) + coord_cartesian(xlim = c(HIGH * 0.9, lim), ylim = c(min(0, min(G$lfc_BE, na.rm = TRUE)) - 0.2, lim)) +
  labs(title = "Strongly IFNB-induced genes: is the induction kept when EBOV is present?", subtitle = sprintf("lines = retention 0 (back to Mock), %g, %g and 1 (kept in full)", R_FULL, R_NOT), x = "log2FC IFNB vs Mock", y = "log2FC IFNB+EBOV vs Mock") + base_theme + theme(legend.position = "bottom")
save2(ps, "scatter_IFNB_vs_IFNB_EBOV", 7.5, 7)

# ---- profiles ----
Lp <- bind_rows(G |> transmute(label, group, condition = "Mock", lfc = 0), G |> transmute(label, group, condition = "IFNB", lfc = lfc_B), G |> transmute(label, group, condition = "IFNB+EBOV", lfc = lfc_BE)) |> mutate(condition = factor(condition, levels = c("Mock", "IFNB", "IFNB+EBOV")))
Mn <- Lp |> group_by(group, condition) |> summarise(lfc = mean(lfc), .groups = "drop")
pp <- ggplot() + geom_line(data = Lp, aes(condition, lfc, group = label, colour = group), alpha = 0.25, linewidth = 0.3) + geom_line(data = Mn, aes(condition, lfc, group = group), colour = INK, linewidth = 1.4) + geom_point(data = Mn, aes(condition, lfc), colour = INK, size = 2.5) +
  facet_wrap(~group, nrow = 1) + scale_colour_manual(values = GC, guide = "none") + labs(title = "Expression profile of each group (log2FC vs Mock)", subtitle = "thin lines = single genes, thick line = group mean", x = NULL, y = "log2FC vs Mock") + base_theme
save2(pp, "profiles_by_group", 9, 4.5)

# ---- heatmap ----
suppressPackageStartupMessages(library(DESeq2)); vsd <- readRDS(file.path(RNA, "vsd_v2.rds")); vm <- SummarizedExperiment::assay(vsd); grp <- as.character(vsd$group)
Gs <- c(Untreated = "Mock_Mock", IFNB = "IFNB_Mock", `IFNB+EBOV` = "IFNB_EBOV"); cols <- lapply(Gs, function(g) colnames(vm)[grp == g]); nrep <- lengths(cols)
Gh <- G |> arrange(group, desc(lfc_B)) |> filter(gene_id %in% rownames(vm)); X <- vm[Gh$gene_id, unlist(cols), drop = FALSE]; X <- X - rowMeans(X[, cols[[1]], drop = FALSE])
rn <- make.unique(Gh$label); rownames(X) <- rn; colnames(X) <- paste0(rep(names(cols), nrep), " ", unlist(lapply(nrep, seq_len)))
show <- Gh |> group_by(group) |> mutate(r = row_number()) |> ungroup() |> mutate(lbl = ifelse(r <= N_LAB, label, "")); ann <- data.frame(group = Gh$group, row.names = rn); lim_h <- as.numeric(quantile(abs(X), CAP_Q))
ph <- pheatmap(X, color = colorRampPalette(c("#2B6CB0", "#F2F2F2", "#C53030"))(100), breaks = seq(-lim_h, lim_h, length.out = 101), cluster_rows = FALSE, cluster_cols = FALSE, border_color = NA, annotation_row = ann, annotation_colors = list(group = GC),
               labels_row = show$lbl, fontsize_row = 7, fontsize_col = 9, angle_col = 45, cellwidth = 22, cellheight = 2.6, gaps_row = cumsum(table(Gh$group))[-length(table(Gh$group))], gaps_col = cumsum(nrep)[-length(nrep)], silent = TRUE,
               main = sprintf("Strongly IFNB-induced genes by EBOV effect (%d genes): each replicate vs the Untreated mean", nrow(X)))
W <- 0.5 * ncol(X) + 3.2; H2 <- 0.0365 * nrow(X) + 2.6
png(file.path(OUT, "heatmap_groups.png"), width = W, height = H2, units = "in", res = 300); grid.newpage(); grid.draw(ph$gtable); invisible(dev.off())
pdf(file.path(OUT, "heatmap_groups.pdf"), width = W, height = H2); grid.newpage(); grid.draw(ph$gtable); invisible(dev.off())
message("Saved to: ", OUT, "\nNext: python 108_promoter_tf_groups.py reads ifnb_induced_groups.csv for the promoter / transcription factor comparison")
