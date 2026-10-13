#!/usr/bin/env bash
# Does ANY read in these libraries come from the GFP reporter?
#
# Primary evidence: exact k-mer matching. Every read is checked for an exact
# 31-mer tiled from your GFP sequence, on both strands. A 31-mer of jellyfish-
# derived GFP does not occur in the mouse genome, so a read carrying one is
# diagnostic - there is no alignment score to argue about and no index that can
# force-assign a host read to the only sequence on offer. Zero hits is a real
# zero (for the variant you supplied); non-zero hits are real reads you can
# look at by eye.
#
# This answers presence/absence only. Abundance comes later, from the full
# rebuild where the whole transcriptome competes for every read.
set -euo pipefail
cd "$(dirname "$0")" && source ../config.sh

: "${GFP_FA:=$PROJ/ref/gfp.fa}"
: "${MAX_READS:=4000000}"   # reads screened per mate; 0 = every read
: "${KEEP_HITS:=200}"       # matching reads saved per sample for inspection
: "${K:=31}"
: "${STEP:=10}"

OUT="$PROJ/gfp_screen"
TSV="$OUT/gfp_kmer_hits.tsv"
PAT="$OUT/gfp_kmers.txt"

if [[ ! -s "$GFP_FA" ]]; then
    cat >&2 <<EOF
No GFP sequence found at:
  $GFP_FA

Put your reporter's nucleotide sequence there (FASTA, one record), or set
GFP_FA to point at it. See gfp/README.md - it matters that this is the variant
actually in your construct, since the screen is exact-match.
EOF
    exit 1
fi

mkdir -p "$OUT"

echo "===== validating $GFP_FA ====="
python3 -I gfp_tools.py validate "$GFP_FA"
echo

python3 -I gfp_tools.py kmers "$GFP_FA" "$K" "$STEP" > "$PAT"
echo "Search set: $(wc -l < "$PAT" | tr -d ' ') ${K}-mers (both strands, every ${STEP} nt)"

if [[ "$MAX_READS" -eq 0 ]]; then
    echo "Depth: every read in every sample (slow - hours; full sensitivity)"
else
    printf "Depth: first %s reads per mate per sample (fast pass)\n" "$MAX_READS"
    echo "       re-run with MAX_READS=0 for full depth once you know where to look"
fi
echo

printf "sample\treads_screened\tR1_hits\tR2_hits\ttotal_hits\thits_per_10M\n" > "$TSV"

screen_one_mate() {  # $1 = fastq.gz   $2 = file to collect matching reads into
    local fq="$1" sink="$2"
    local ctr; ctr="$(mktemp)"
    set +o pipefail
    gzip -cd "$fq" 2>/dev/null \
        | awk -v max="$MAX_READS" -v cf="$ctr" \
              'NR%4==2 { n++; print; if (max > 0 && n >= max) exit }
               END     { print n+0 > cf }' \
        | grep -F -f "$PAT" >> "$sink" || true
    set -o pipefail
    MATE_READS="$(cat "$ctr" 2>/dev/null || echo 0)"
    rm -f "$ctr"
}

for s in "${SAMPLES[@]}"; do
    r1="$PROJ/trimmed_fastq/${s}_R1.trim.fastq.gz"
    r2="$PROJ/trimmed_fastq/${s}_R2.trim.fastq.gz"
    if [[ ! -s "$r1" || ! -s "$r2" ]]; then
        echo "  sample $s: trimmed FASTQs missing - skipped" >&2
        continue
    fi

    hits_file="$OUT/hits_${s}.txt"
    : > "$hits_file"

    printf "  sample %-5s " "$s"

    screen_one_mate "$r1" "$hits_file"; reads_r1="$MATE_READS"
    h1=$(wc -l < "$hits_file" | tr -d ' ')

    screen_one_mate "$r2" "$hits_file"; reads_r2="$MATE_READS"
    h_all=$(wc -l < "$hits_file" | tr -d ' ')
    h2=$(( h_all - h1 ))

    screened=$(( reads_r1 + reads_r2 ))
    if [[ "$screened" -gt 0 ]]; then
        per10m=$(python3 -c "print(f'{$h_all / $screened * 1e7:.2f}')")
    else
        per10m="NA"
    fi

    printf "screened %12s   hits %6s  (%s per 10M)\n" "$screened" "$h_all" "$per10m"
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$s" "$screened" "$h1" "$h2" "$h_all" "$per10m" >> "$TSV"

    # Keep a bounded sample of the matching reads; drop the file if empty.
    if [[ "$h_all" -eq 0 ]]; then
        rm -f "$hits_file"
    else
        head -n "$KEEP_HITS" "$hits_file" > "$hits_file.keep" && mv "$hits_file.keep" "$hits_file"
    fi
done

echo
echo "===== SUMMARY ====="
column -t -s "$(printf '\t')" "$TSV" 2>/dev/null || cat "$TSV"

echo
echo "Table:          $TSV"
echo "Matching reads: $OUT/hits_<sample>.txt  (sequence lines, for eyeballing)"
cat <<'EOF'

How to read this
----------------
Samples 1-3 and 10-15 never saw the construct, so they are your negative
controls. Whatever number they show is the floor: index hopping, carryover
between libraries on the same run, or a mis-specified reference. Judge every
other sample against that floor, not against zero.

  all samples ~0        no GFP-derived RNA detectable. Either the reporter is
                        not transcribed in this system, the samples that should
                        carry it do not, or the sequence you supplied is not the
                        variant in the construct. The last one is the cheapest
                        to rule out - re-read the validate output above.

  only expected
  samples positive      real signal. Look at hits_<sample>.txt, then move on to
                        the replicon rebuild for actual numbers.

  negatives positive
  at similar levels     contamination or hopping rather than biology. Compare
                        against the run's index-assignment metrics before
                        drawing any conclusion from the positives.
EOF
