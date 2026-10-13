#!/usr/bin/env Rscript
# How much did the fixes actually move the results? Compares the original
# IFNB vs Mock and IFNY vs Mock calls against the re-quantified ones.
#
# This is the file to look at before redoing any downstream figure: if the
# interferon log2 fold-changes are essentially unchanged, your existing ISG and
# IDG conclusions stand and only the marginal genes need revisiting.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(tibble)
})

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs")
dir.create(OUT, showWarnings = FALSE)

pairs <- list(
  IFNB_vs_Mock = list(old = "DE_IFNB_vs_Mock_fresh.csv", new = "DE_v2_IFNB_vs_Mock.csv"),
  IFNY_vs_Mock = list(old = "DE_IFNY_vs_Mock_fresh.csv", new = "DE_v2_IFNY_vs_Mock.csv")
)

read_de <- function(path) {
  if (!file.exists(path)) return(NULL)
  d <- suppressWarnings(read_csv(path, show_col_types = FALSE))
  id <- intersect(c("gene_id", "gene", "...1", "GENEID", "row"), names(d))[1]
  if (is.na(id)) stop("Cannot find a gene id column in ", path)
  d |> transmute(gene_id = .data[[id]], log2FoldChange, padj)
}

summary_rows <- list()

for (nm in names(pairs)) {
  old <- read_de(file.path(PROJ, pairs[[nm]]$old))
  new <- read_de(file.path(OUT,  pairs[[nm]]$new))
  if (is.null(old) || is.null(new)) { message("skip ", nm, " (file missing)"); next }

  m <- inner_join(old, new, by = "gene_id", suffix = c("_old", "_new")) |>
    mutate(
      sig_old = !is.na(padj_old) & padj_old < 0.05 & abs(log2FoldChange_old) > 1,
      sig_new = !is.na(padj_new) & padj_new < 0.05 & abs(log2FoldChange_new) > 1,
      status  = case_when(sig_old &  sig_new ~ "kept",
                          sig_old & !sig_new ~ "lost",
                         !sig_old &  sig_new ~ "gained",
                          TRUE               ~ "ns")
    )

  r <- cor(m$log2FoldChange_old, m$log2FoldChange_new,
           use = "complete.obs", method = "pearson")
  tally <- m |> count(status)
  message("\n=== ", nm, " ===")
  message("genes compared: ", nrow(m), "   log2FC Pearson r: ", round(r, 4))
  print(as.data.frame(tally), row.names = FALSE)

  summary_rows[[nm]] <- tibble(
    contrast = nm, n_compared = nrow(m), pearson_r = r,
    kept   = sum(m$status == "kept"),
    lost   = sum(m$status == "lost"),
    gained = sum(m$status == "gained")
  )

  write_csv(m |> filter(status %in% c("lost", "gained")) |> arrange(status, padj_new),
            file.path(OUT, paste0("changed_calls_", nm, ".csv")))

  p <- ggplot(m, aes(log2FoldChange_old, log2FoldChange_new)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey70", linewidth = 0.4) +
    geom_point(aes(colour = status), alpha = 0.5, size = 0.7) +
    scale_colour_manual(values = c(kept = "#2B6CB0", lost = "#C53030",
                                   gained = "#2F855A", ns = "grey80")) +
    labs(title = paste0(nm, "  (r = ", round(r, 3), ")"),
         subtitle = "original index vs decoy-aware + bias-corrected",
         x = "log2FC, original", y = "log2FC, re-quantified", colour = NULL) +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
  ggsave(file.path(OUT, paste0("compare_", nm, ".png")), p,
         width = 6, height = 5.5, dpi = 300)
}

if (length(summary_rows)) {
  s <- bind_rows(summary_rows)
  write_csv(s, file.path(OUT, "old_vs_new_summary.csv"))
  message("\n--- summary ---")
  print(as.data.frame(s), row.names = FALSE)
}
message("\nOutputs in: ", OUT)
