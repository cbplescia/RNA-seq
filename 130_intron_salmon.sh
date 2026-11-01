#!/usr/bin/env bash
# Step 2 of the intron / exon check: build a Salmon index of (transcripts + construct + the intron subset from 129 + genome decoys) and quantify the Mock and EBOV samples listed in intron_samples.txt against it.
# Same settings as 03_salmon_quant.sh (library type ISR, --seqBias --gcBias, k = 31). Only the 6 samples of the two groups are quantified. Resumable (skips a finished index / sample).
# Needs: python 129_build_intron_reference.py already run (ref/introns_subset.fa, intron_samples.txt).   Index build ~30-60 min and about as much RAM as your first index; the gentrome file (~3.5 GB) can be deleted afterwards.
# Usage: bash 130_intron_salmon.sh   (after 'source setup_session.sh'; config.sh must sit next to this script)
set -euo pipefail
cd "$(dirname "$0")" && source ./config.sh
INTR="$PROJ/ref/introns_subset.fa"; SAMP="$PROJ/intron_samples.txt"; IDX_I="$PROJ/salmon_index_introns"; QD="$PROJ/salmon_quant_introns"; GENTI="$PROJ/ref/gentrome_with_introns.fa"
[[ -s "$INTR" && -s "$SAMP" ]] || { echo "Run  python 129_build_intron_reference.py  first (missing $INTR or $SAMP)"; exit 1; }
for f in "$TX_FA" "$GENOME_FA"; do [[ -s "$f" ]] || { echo "Missing $f"; exit 1; }; done
[[ -s "$DECOYS" ]] || grep '^>' "$GENOME_FA" | cut -d ' ' -f1 | sed 's/^>//' > "$DECOYS"
if [[ ! -d "$IDX_I" ]]; then
    echo "Building the gentrome (transcripts + construct + introns + genome decoys) ..."
    if [[ -s "$VIRAL_FA" ]]; then cat "$TX_FA" "$VIRAL_FA" "$INTR" "$GENOME_FA" > "$GENTI"; else cat "$TX_FA" "$INTR" "$GENOME_FA" > "$GENTI"; fi
    echo "Building the Salmon index (slow) ..."
    salmon index -t "$GENTI" -d "$DECOYS" -i "$IDX_I" -k 31 -p "$THREADS"
fi
mkdir -p "$QD"
while IFS=$'\t' read -r s grp; do
    [[ -n "$s" ]] || continue
    out="$QD/$s"
    if [[ -s "$out/quant.sf" ]]; then echo "[$s] done, skipping"; continue; fi
    r1="$PROJ/trimmed_fastq/${s}_R1.trim.fastq.gz"; r2="$PROJ/trimmed_fastq/${s}_R2.trim.fastq.gz"
    [[ -s "$r1" && -s "$r2" ]] || { echo "[$s] MISSING FASTQ, skipping"; continue; }
    echo "[$s] ($grp) quantifying ..."
    salmon quant -i "$IDX_I" -l "$LIBTYPE" -1 "$r1" -2 "$r2" -p "$THREADS" --seqBias --gcBias -o "$out"
done < "$SAMP"
echo; echo "Done. Next: Rscript 131_intron_exon_analysis.R"
