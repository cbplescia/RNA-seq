#!/usr/bin/env bash
# Shared settings for the Maury_25161 re-quantification.
# Edit PROJ if the project moves. Everything else is derived from it.

# Root of the results directory that holds raw_fastq/, trimmed_fastq/, ref/, logs/
PROJ="${PROJ:-$HOME/Documents/aviti_rnaseq/2025/aviti/20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz/results/Maury_25161}"

# Existing Ensembl GRCm39 release 109 reference files (already present)
GENOME_FA="$PROJ/ref/Mus_musculus.GRCm39.dna.primary_assembly.fa"
TX_FA="$PROJ/ref/Mus_musculus.GRCm39.transcripts.fa"
GTF="$PROJ/ref/Mus_musculus.GRCm39.109.gtf"

# Recombinant construct / viral sequences. One FASTA, one record per feature.
# See viral_refs/README.md in the git repo for how to build this.
VIRAL_FA="${VIRAL_FA:-$PROJ/ref/recombinant_model.fa}"

# Outputs of the rebuild (kept separate from the originals so nothing is overwritten)
GENTROME="$PROJ/ref/gentrome_GRCm39.109_plus_model.fa"
DECOYS="$PROJ/ref/decoys.txt"
IDX="$PROJ/salmon_index_decoy_model"
QUANT_DIR="$PROJ/salmon_quant_v2"
TX2GENE="$PROJ/tx2gene_GRCm39.109_plus_model.tsv"

# Library type: detected unambiguously as ISR on sample 1 (dUTP-style stranded,
# read 1 on the reverse strand). Set explicitly so it cannot vary per sample.
LIBTYPE="ISR"

THREADS="${THREADS:-8}"

# Sample IDs, in the order used by sample_info.csv
SAMPLES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21_D 22 23 24 \
         25_A 26_B 27_C 28 29 30 31 32 33)
