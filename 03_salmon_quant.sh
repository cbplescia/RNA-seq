#!/usr/bin/env bash
# Re-quantify all 33 samples against the decoy-aware index.
#
# Changes from the original run:
#   -l ISR        explicit, instead of -l A. ISR was detected unambiguously on
#                 sample 1 (compatible_fragment_ratio 1.0, strand bias 1.4e-5).
#                 Fixing it removes the chance that one sample gets typed
#                 differently and has its counts silently flipped.
#   --seqBias     corrects for random-hexamer priming bias at fragment starts.
#   --gcBias      corrects fragment-level GC bias; recommended for any DE run.
#   --writeUnmappedNames  so the unmapped fraction can be inspected rather than
#                 assumed. Needed to confirm the construct reads are now being
#                 captured.
#   --validateMappings is the default in 1.10 and is no longer passed explicitly.
#
# Resumable: skips any sample that already has a complete quant.sf.
set -euo pipefail
cd "$(dirname "$0")" && source ./config.sh

[[ -d "$IDX" ]] || { echo "No index at $IDX - run 02_build_reference.sh first"; exit 1; }
mkdir -p "$QUANT_DIR"

for s in "${SAMPLES[@]}"; do
    out="$QUANT_DIR/$s"
    if [[ -s "$out/quant.sf" && -s "$out/aux_info/meta_info.json" ]]; then
        echo "[$s] already done, skipping"
        continue
    fi
    r1="$PROJ/trimmed_fastq/${s}_R1.trim.fastq.gz"
    r2="$PROJ/trimmed_fastq/${s}_R2.trim.fastq.gz"
    [[ -s "$r1" && -s "$r2" ]] || { echo "[$s] MISSING FASTQ, skipping"; continue; }

    echo "[$s] quantifying..."
    salmon quant \
        -i "$IDX" \
        -l "$LIBTYPE" \
        -1 "$r1" -2 "$r2" \
        -p "$THREADS" \
        --seqBias \
        --gcBias \
        --writeUnmappedNames \
        -o "$out"
done

echo
echo "===== MAPPING SUMMARY (new run) ====="
printf "%-8s %12s %12s %8s %10s %12s\n" sample processed mapped pct libtype decoy_frags
for s in "${SAMPLES[@]}"; do
    f="$QUANT_DIR/$s/aux_info/meta_info.json"
    [[ -s "$f" ]] || continue
    python3 - "$f" "$s" <<'PY'
import json, sys
m = json.load(open(sys.argv[1]))
print(f"{sys.argv[2]:<8} {m['num_processed']:>12,} {m['num_mapped']:>12,} "
      f"{m['percent_mapped']:>7.2f}% {','.join(m['library_types']):>10} "
      f"{m.get('num_decoy_fragments',0):>12,}")
PY
done

cat <<'EOF'

What to expect:
  - Mapping rate DROPS by a few points versus the original 85-91%. That is the
    fix working: reads that were previously forced onto transcripts are now
    being absorbed by the genome decoys. num_decoy_fragments shows how many.
  - library_types must read ISR for all 33.
EOF
