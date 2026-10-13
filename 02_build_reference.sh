#!/usr/bin/env bash
# Build a decoy-aware Salmon index over the Ensembl GRCm39 release 109
# transcriptome PLUS the recombinant construct.
#
# Two fixes over the original index:
#   1. Genome decoys. The original had num_decoys: 0, so reads from intronic,
#      intergenic and unprocessed-transcript sequence had nowhere to go and were
#      force-assigned to whichever transcript they partially resembled. Adding
#      the primary assembly as decoys gives those reads a home and they drop out.
#   2. The construct is added as real quantification targets, so construct-derived
#      RNA is counted instead of being silently distributed over host transcripts.
#
# Writes to new paths - the original salmon_index/ is left untouched.
set -euo pipefail
cd "$(dirname "$0")" && source ./config.sh

[[ -s "$GENOME_FA" ]] || { echo "Missing genome: $GENOME_FA"; exit 1; }
[[ -s "$TX_FA"     ]] || { echo "Missing transcriptome: $TX_FA"; exit 1; }
[[ -s "$VIRAL_FA"  ]] || { echo "Missing construct FASTA: $VIRAL_FA - see viral_refs/README.md"; exit 1; }

# --- sanity: no duplicate sequence names across the three inputs -------------
echo "Checking for duplicate sequence names..."
dups=$( { grep '^>' "$TX_FA" "$VIRAL_FA" "$GENOME_FA" | sed 's/.*>//; s/[[:space:]].*//'; } \
        | sort | uniq -d )
if [[ -n "$dups" ]]; then
    echo "Duplicate sequence names found - Salmon will refuse or silently drop these:"
    echo "$dups" | sed 's/^/  /'
    exit 1
fi

# --- gentrome: targets first, decoys last -----------------------------------
echo "Building gentrome (transcripts + construct + genome decoys)..."
cat "$TX_FA" "$VIRAL_FA" "$GENOME_FA" > "$GENTROME"

# Decoy list = genome records only. The construct records are NOT decoys.
grep '^>' "$GENOME_FA" | cut -d ' ' -f1 | sed 's/^>//' > "$DECOYS"
echo "  $(wc -l < "$DECOYS") decoy sequences"

# --- index -------------------------------------------------------------------
# k=31 is right for 75 bp reads. Add --sparse here if the build runs out of RAM.
echo "Building index (this is the slow step - 30-60 min, ~16 GB RAM)..."
salmon index \
    -t "$GENTROME" \
    -d "$DECOYS" \
    -i "$IDX" \
    -k 31 \
    -p "$THREADS"

# --- tx2gene -----------------------------------------------------------------
# Host entries come from the existing map. Construct features are their own
# "gene" so they survive gene-level summarisation individually.
echo "Building tx2gene..."
{
    # existing host map, minus its header
    tail -n +2 "$PROJ/tx2gene_mouse_from_fasta.tsv"
    # construct: feature maps to itself
    grep '^>' "$VIRAL_FA" | cut -d ' ' -f1 | sed 's/^>//' | awk '{print $1"\t"$1}'
} | awk 'BEGIN{print "TXNAME\tGENEID"} {print}' > "$TX2GENE"
echo "  $(( $(wc -l < "$TX2GENE") - 1 )) transcript->gene pairs"

# Record what went in, for the methods section.
cat > "$PROJ/ref/REBUILD_PROVENANCE.txt" <<EOF
Built: $(date)
salmon: $(salmon --version 2>&1)
transcriptome: $TX_FA
  records: $(grep -c '^>' "$TX_FA")
construct:    $VIRAL_FA
  records: $(grep -c '^>' "$VIRAL_FA")
$(grep '^>' "$VIRAL_FA" | sed 's/^/    /')
decoys:       $GENOME_FA ($(wc -l < "$DECOYS") sequences)
index:        $IDX (k=31, decoy-aware)
tx2gene:      $TX2GENE
EOF

echo
echo "Done. Index: $IDX"
echo "Provenance: $PROJ/ref/REBUILD_PROVENANCE.txt"
echo "The gentrome ($(du -h "$GENTROME" | cut -f1)) can be deleted once the index is built."
