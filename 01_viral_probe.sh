#!/usr/bin/env bash
# Cheap presence/absence screen: is there ANY construct-derived RNA in the
# libraries? Builds a tiny index from the construct alone and quantifies a few
# samples. Minutes, not hours - run this before the full rebuild.
#
# IMPORTANT: a construct-only index has no host sequence to compete, so host
# reads that partially resemble the construct have nowhere else to go and will
# be force-assigned to it. Read the output as "is there a signal worth chasing",
# never as a count. The real numbers come from 03_salmon_quant.sh, where the
# full transcriptome and genome decoys compete for every read.
set -euo pipefail
cd "$(dirname "$0")" && source ./config.sh

[[ -s "$VIRAL_FA" ]] || { echo "No construct FASTA at $VIRAL_FA - see viral_refs/README.md"; exit 1; }

PROBE_IDX="$PROJ/ref/probe_index_model_only"
PROBE_OUT="$PROJ/salmon_probe_model_only"

# Expect signal in: an infected sample, a VLP sample, and a mock as the negative.
PROBE_SAMPLES=(4 16 28 1)

echo "Building construct-only probe index..."
salmon index -t "$VIRAL_FA" -i "$PROBE_IDX" -k 31 -p "$THREADS"

mkdir -p "$PROBE_OUT"
for s in "${PROBE_SAMPLES[@]}"; do
    echo "Probing sample $s ..."
    salmon quant -i "$PROBE_IDX" -l "$LIBTYPE" \
        -1 "$PROJ/trimmed_fastq/${s}_R1.trim.fastq.gz" \
        -2 "$PROJ/trimmed_fastq/${s}_R2.trim.fastq.gz" \
        -p "$THREADS" -o "$PROBE_OUT/$s" >/dev/null 2>&1
done

echo
echo "===== READS ASSIGNED TO THE CONSTRUCT (upper bound) ====="
printf "%-8s %12s %12s\n" sample processed mapped
for s in "${PROBE_SAMPLES[@]}"; do
    python3 - "$PROBE_OUT/$s/aux_info/meta_info.json" "$s" <<'PY'
import json, sys
m = json.load(open(sys.argv[1]))
print(f"{sys.argv[2]:<8} {m['num_processed']:>12,} {m['num_mapped']:>12,}")
PY
done

echo
echo "===== PER-FEATURE TPM ====="
for s in "${PROBE_SAMPLES[@]}"; do
    echo "--- sample $s"
    awk 'NR==1 || $4>0 {printf "    %-30s len=%-8s TPM=%-12s reads=%s\n",$1,$2,$4,$5}' \
        "$PROBE_OUT/$s/quant.sf"
done

cat <<'EOF'

Interpretation:
  Mock (sample 1) should be near zero. If mock carries as many reads as the
  infected samples, the construct FASTA shares sequence with the mouse
  transcriptome - most likely a mouse-derived cDNA, promoter or UTR in the
  backbone. Find which record it is and split or mask it before the full run.
EOF
