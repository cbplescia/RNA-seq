#!/usr/bin/env bash
# Check that the inputs are all present and that there is room to do the rebuild.
# Read-only: changes nothing.
set -uo pipefail
cd "$(dirname "$0")" && source ./config.sh

echo "PROJ = $PROJ"
echo
echo "===== REQUIRED INPUTS ====="
for f in "$GENOME_FA" "$TX_FA" "$GTF"; do
    if [[ -s "$f" ]]; then echo "  OK      $f"; else echo "  MISSING $f"; fi
done
if [[ -s "$VIRAL_FA" ]]; then
    echo "  OK      $VIRAL_FA"
    echo "          $(grep -c '^>' "$VIRAL_FA") records:"
    grep '^>' "$VIRAL_FA" | sed 's/^/            /'
else
    echo "  MISSING $VIRAL_FA   <-- build this first, see viral_refs/README.md"
fi

echo
echo "===== TRIMMED FASTQ PAIRS ====="
missing=0
for s in "${SAMPLES[@]}"; do
    for r in R1 R2; do
        f="$PROJ/trimmed_fastq/${s}_${r}.trim.fastq.gz"
        [[ -s "$f" ]] || { echo "  MISSING $f"; missing=$((missing+1)); }
    done
done
echo "  ${#SAMPLES[@]} samples expected, $missing files missing"

echo
echo "===== INTEGRITY OF TRIMMED FASTQ (gzip streams) ====="
echo "  (slow - a few minutes. Ctrl-C to skip; rerun this section alone later.)"
for s in "${SAMPLES[@]}"; do
    for r in R1 R2; do
        f="$PROJ/trimmed_fastq/${s}_${r}.trim.fastq.gz"
        [[ -s "$f" ]] || continue
        if gzip -t "$f" 2>/dev/null; then :; else echo "  CORRUPT $f"; fi
    done
done
echo "  done"

echo
echo "===== DISK AND MEMORY ====="
du -sh "$PROJ/raw_fastq" "$PROJ/trimmed_fastq" "$PROJ/ref" 2>/dev/null
df -h "$PROJ" | tail -2
if command -v sysctl >/dev/null; then
    echo "  physical RAM: $(( $(sysctl -n hw.memsize) / 1073741824 )) GB"
fi
cat <<'EOF'

  A decoy-aware mouse index needs roughly 25-30 GB of free disk for the
  gentrome plus index, and peaks around 16 GB of RAM while building.
  If the build is killed for memory, add --sparse to the salmon index call
  in 02_build_reference.sh (smaller index, slightly slower quant, same results).
EOF

echo
echo "===== TOOLS ====="
echo "  salmon: $(salmon --version 2>&1 || echo NOT FOUND)"
echo "  (build the index from the env that has salmon 1.10.3: conda activate salmon_env)"
