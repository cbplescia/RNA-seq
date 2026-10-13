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
: "${MODE:=nt}"             # nt = exact nucleotide; aa = translated (see below)
: "${K:=31}"                # nucleotide k-mer length
: "${STEP:=10}"             # nucleotide tiling stride
: "${PEPK:=9}"              # peptide window length for MODE=aa

OUT="$PROJ/gfp_screen"
TSV="$OUT/gfp_kmer_hits.tsv"
PAT="$OUT/gfp_kmers.txt"
PEP="$OUT/gfp_peptides.txt"

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

SCRIPT_DIR="$PWD"

if [[ "$MODE" == "aa" ]]; then
    python3 -I gfp_tools.py peptides "$GFP_FA" "$PEPK" > "$PEP"
    echo "Search set: $(wc -l < "$PEP" | tr -d ' ') ${PEPK}-aa peptide windows"
    echo "Mode: TRANSLATED - six-frame translation of every read."
    echo "      Finds a GFP whose codon usage differs from the reference"
    echo "      (EGFP vs wild-type avGFP), at maybe 10x the runtime of nt mode."
    filter_stage() { python3 -I "$SCRIPT_DIR/aa_filter.py" "$PEP"; }
else
    python3 -I gfp_tools.py kmers "$GFP_FA" "$K" "$STEP" > "$PAT"
    echo "Search set: $(wc -l < "$PAT" | tr -d ' ') ${K}-mers (both strands, every ${STEP} nt)"
    echo "Mode: EXACT nucleotide - only finds the variant you supplied."
    filter_stage() { grep -F -f "$PAT"; }
fi

# Positive control on the search machinery itself. A stretch of the reference
# must come back through the filter. This cannot tell you whether the reference
# is the right GFP variant - it catches the other way of getting a zero, where
# the filter is silently broken and every sample looks empty.
CONTROL_SEQ="$(python3 -I - "$GFP_FA" <<'PY'
import sys
sys.path.insert(0, ".")
from gfp_tools import longest_orf, read_fasta, translate
seq = read_fasta(sys.argv[1])[0][1]
if translate(seq)[:-1].count("*"):
    found = longest_orf(seq)
    if found:
        seq = found[0]
print(seq[:72])
PY
)"
if [[ -z "$(printf '%s\n' "$CONTROL_SEQ" | filter_stage || true)" ]]; then
    echo >&2
    echo "FAILED self-check: a sequence taken straight from the reference did" >&2
    echo "not come back through the filter. The search is broken, so every" >&2
    echo "sample would read as zero. Not screening." >&2
    exit 1
fi
echo "Self-check: reference sequence is detected by the filter."

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
        | filter_stage >> "$sink" || true
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
    if [[ "$screened" -eq 0 ]]; then
        echo >&2
        echo "ABORT: sample $s counted zero reads screened, but its FASTQs are" >&2
        echo "non-empty. The read stream broke rather than finding nothing." >&2
        exit 1
    fi
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
