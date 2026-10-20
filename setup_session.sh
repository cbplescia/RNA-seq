#!/usr/bin/env bash
# Re-establish the working session after a restart. Run from INSIDE the Maury_25161 folder with:
#     source setup_session.sh
# (source it, do not execute it, so the conda environment and PROJ stay active in your terminal)
#   1. activates the rnaseq_env conda environment
#   2. sets PROJ to the current folder, which every analysis script reads
#   3. checks that the key input files are present
#   4. checks which R packages are installed and prints the install command for any that are missing

conda activate rnaseq_env || { echo "Could not activate rnaseq_env - run 'conda init zsh' once, or open a new terminal"; return 1 2>/dev/null; }
export PROJ="$(pwd)"
echo "Environment: $CONDA_DEFAULT_ENV    PROJ = $PROJ"

echo; echo "---- key files ----"
for f in requant_v2_outputs/vsd_v2.rds requant_v2_outputs/dds_v2_decoyAware_hostNormalised.rds \
         requant_v2_outputs/DE_v2_IFNB_vs_Mock_protein_coding.csv requant_v2_outputs/DE_v2_IFNB_vs_Mock.csv \
         sample_info.csv 01_lipidomics_ranked_IFNB_vs_Untreated.csv; do
  if [ -e "$f" ]; then echo "  OK       $f"; else echo "  MISSING  $f"; fi
done

echo; echo "---- R packages ----"
Rscript -e '
core <- c("DESeq2", "readr", "dplyr", "tidyr", "ggplot2", "ggrepel")
opt  <- c("pheatmap", "patchwork", "clusterProfiler", "enrichplot", "org.Mm.eg.db", "GO.db", "msigdbr", "ReactomePA", "STRINGdb")
have <- function(p) requireNamespace(p, quietly = TRUE)
mc <- core[!vapply(core, have, logical(1))]; mo <- opt[!vapply(opt, have, logical(1))]
cat("  core packages missing:    ", if (length(mc)) paste(mc, collapse = ", ") else "none", "\n")
cat("  optional packages missing:", if (length(mo)) paste(mo, collapse = ", ") else "none", "\n")
if (length(mc) || length(mo)) cat("  install with e.g.:  conda install -c conda-forge -c bioconda",
  paste0(c("r-", "bioconductor-")[1], tolower(c(mc, mo)[1]), " ..."), "\n   (r-<name> for CRAN packages, bioconductor-<name> for Bioconductor ones)\n")
' 2>/dev/null
echo; echo "Ready. Scripts read PROJ, so run them from anywhere, e.g.:  Rscript 45_pca_mock_ifnb_ifny.R"
