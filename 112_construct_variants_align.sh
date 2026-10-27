#!/usr/bin/env bash
# Base-level look at the recombinant construct: align every sample's trimmed reads end-to-end to ref/recombinant_model.fa ONLY (bowtie2), keep confidently placed reads (MAPQ >= MIN_MAPQ), and write a samtools pileup per sample.
# Why a new alignment: Salmon (quasi-mapping) does not write alignments, so it cannot tell you which bases differ from the reference. Reads from the mouse transcriptome are not removed here, but reads that are not
# construct-like simply do not align end-to-end; any host read that does is a sign that the construct shares sequence with mouse (check Mock samples: they should have ~0 depth).
# Output (in $PROJ/construct_variants/): bam/<sample>.bam(.bai), pileup/<sample>.pileup.txt, alignment_summary.tsv.     Next: python 113_construct_variants.py
# Needs bowtie2 + samtools:  conda install -c conda-forge -c bioconda bowtie2 samtools      Resumable (skips samples with a pileup).   Usage: bash 112_construct_variants_align.sh   (after 'source setup_session.sh')
set -euo pipefail
cd "$(dirname "$0")" && source ./config.sh
MIN_MAPQ="${MIN_MAPQ:-10}"; MIN_BQ="${MIN_BQ:-20}"
command -v bowtie2 >/dev/null && command -v samtools >/dev/null || { echo "Need bowtie2 and samtools: conda install -c conda-forge -c bioconda bowtie2 samtools"; exit 1; }
[[ -s "$VIRAL_FA" ]] || { echo "Missing construct FASTA: $VIRAL_FA"; exit 1; }
W="$PROJ/construct_variants"; mkdir -p "$W/bam" "$W/pileup" "$W/index"
cp "$VIRAL_FA" "$W/index/model.fa"; samtools faidx "$W/index/model.fa"
[[ -s "$W/index/model.1.bt2" ]] || bowtie2-build --quiet "$W/index/model.fa" "$W/index/model"
printf "sample\tpairs_in\tpairs_aligned_construct\tpct\treads_MAPQ_ok\n" > "$W/alignment_summary.tsv.tmp"
for s in "${SAMPLES[@]}"; do
    pu="$W/pileup/$s.pileup.txt"
    r1="$PROJ/trimmed_fastq/${s}_R1.trim.fastq.gz"; r2="$PROJ/trimmed_fastq/${s}_R2.trim.fastq.gz"
    if [[ -s "$pu" ]]; then echo "[$s] done, skipping"; continue; fi
    [[ -s "$r1" && -s "$r2" ]] || { echo "[$s] MISSING FASTQ, skipping"; continue; }
    echo "[$s] aligning..."
    bowtie2 --end-to-end --very-sensitive --no-unal --no-mixed --no-discordant -x "$W/index/model" -1 "$r1" -2 "$r2" -p "$THREADS" 2> "$W/bam/$s.bowtie2.log" \
      | samtools view -b -q "$MIN_MAPQ" - | samtools sort -@ 4 -o "$W/bam/$s.bam" -
    samtools index "$W/bam/$s.bam"
    # -aa: report every reference position (also depth 0); -B: no BAQ (raw counts); overlapping mates are counted once
    samtools mpileup -aa -B -Q "$MIN_BQ" -d 10000000 -f "$W/index/model.fa" "$W/bam/$s.bam" > "$pu" 2>/dev/null
    pin=$(grep -m1 "reads; of these" "$W/bam/$s.bowtie2.log" | awk '{print $1}'); pal=$(grep -E "aligned concordantly (exactly 1 time|>1 times)" "$W/bam/$s.bowtie2.log" | awk 'NR>0{n+=$1} END{print n+0}')
    ok=$(samtools view -c -f 64 "$W/bam/$s.bam")
    printf "%s\t%s\t%s\t%.4f\t%s\n" "$s" "$pin" "$pal" "$(awk -v a="$pal" -v b="$pin" 'BEGIN{print (b>0)?100*a/b:0}')" "$ok" >> "$W/alignment_summary.tsv.tmp"
done
mv "$W/alignment_summary.tsv.tmp" "$W/alignment_summary.tsv"
echo "Done. Pileups in $W/pileup. Next: python 113_construct_variants.py"
