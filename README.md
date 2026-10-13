# Maury_25161 — re-quantification from trimmed reads

Rebuild of the Salmon quantification for the 33-sample murine interferon /
filovirus RNA-seq dataset, fixing two problems in the original run and adding
the recombinant construct as a quantification target.

Nothing here overwrites the original outputs. The new index, quants and results
go to `salmon_index_decoy_model/`, `salmon_quant_v2/` and
`requant_v2_outputs/`, so the old and new results can be compared directly and
the old ones can be fallen back on.

## How far back the rebuild goes

| Stage | Reused or redone | Why |
|---|---|---|
| FASTQ generation (bases2fastq) | reused | Produced by the core from the AVITI run; nothing to improve. |
| cutadapt trimming | **reused** | The parameters were correct: Illumina-compatible adapter on both mates, `-q 20,20`, minimum length 50 on 75 bp reads. Re-trimming would consume hours and produce the same files. |
| Transcriptome FASTA | reused | Ensembl GRCm39 release 109, built with gffread. Staying on 109 keeps the new results comparable to the existing analyses. |
| Salmon index | **redone** | Adds genome decoys and the construct. |
| Salmon quant | **redone** | Explicit `-l ISR`, plus sequence and GC bias correction. |
| tximport / DESeq2 | **redone** | Host-only size-factor estimation. |

Going further back than the trimmed FASTQs buys nothing. The two real problems
were both downstream of trimming.

## What is being fixed

**1. The index had no decoys** (`num_decoys: 0`). Reads originating from
intronic, intergenic and unprocessed-transcript sequence had nowhere else to
go, so selective alignment assigned them to whichever transcript they partially
resembled. Adding the primary assembly as decoys gives those reads a correct
destination and they drop out of the counts. Expect the mapping rate to fall a
few points from the original 85–91%; that decrease is the fix working, and
`num_decoy_fragments` quantifies it.

**2. Bias correction was off.** `--seqBias` and `--gcBias` are now on.
Random-hexamer priming bias and fragment GC bias are both sample-specific, so
leaving them uncorrected lets technical variation enter the fold-changes.

**3. Library type was re-detected per sample.** `-l A` inferred the type
independently for each run. It landed on ISR for sample 1 with no ambiguity
(`compatible_fragment_ratio` 1.0, strand bias 1.4e-5), but a single sample
typed differently would have had its counts silently reversed. Now fixed as
`-l ISR` for all 33.

**4. The construct was invisible.** Construct-derived RNA had no target to map
to, so those reads were either discarded or scattered across host transcripts.
They are now counted as named features.

**5. Size factors included construct RNA.** Construct transcripts can be a
large and very uneven share of the library in infected samples. Letting them
into median-of-ratios normalisation drags those samples' size factors and
manufactures apparent host-gene changes. Size factors are now estimated from
host genes only, via `estimateSizeFactors(controlGenes = ...)`, while the
construct stays in the object so it can still be inspected.

## Just checking for GFP reporter reads? Start in `gfp/`

The replicon rebuild below is staged and ready, but it needs the construct
sequences and an hour of index building. If the immediate question is only
*does anything map to the GFP reporter*, `gfp/` answers it on its own — exact
k-mer matching against the trimmed FASTQs, no index, no reindexing, nothing
touched that the interferon analyses depend on.

```bash
cd gfp
./00_make_gfp_fasta.sh ~/path/to/construct.gb   # or set GFP_FA to a FASTA
./01_gfp_kmer_screen.sh                         # all 33 samples, fast pass
python3 -I 02_gfp_hit_profile.py "$PROJ/ref/gfp.fa" "$PROJ"/gfp_screen/hits_*.txt
```

See `gfp/README.md`. The rebuild below is unaffected either way, and the GFP
record validated there drops straight into `recombinant_model.fa` when you come
back to the replicon.

## Running the full rebuild

```bash
conda activate salmon_env          # salmon 1.10.3

# 1. Put the construct FASTA in place first - see viral_refs/README.md
# 2. Check inputs, disk and memory (read-only)
./00_preflight.sh

# 3. Cheap screen: is there construct RNA in the libraries at all? (minutes)
./01_viral_probe.sh

# 4. Build the decoy-aware index (30-60 min, ~16 GB RAM, ~30 GB disk)
./02_build_reference.sh

# 5. Re-quantify all 33 samples (~5 min each)
./03_salmon_quant.sh

# 6. Import, normalise on host genes, run the contrasts
conda activate rnaseq_env
Rscript 04_tximport_deseq2.R

# 7. See how much actually changed
Rscript 05_compare_old_vs_new.R
```

`config.sh` holds the paths; edit `PROJ` there if the project moves.
`03_salmon_quant.sh` is resumable and skips samples already finished.

## The construct deliverable

`04_tximport_deseq2.R` writes three files to `requant_v2_outputs/`:

- `construct_features_per_sample.csv` — every construct feature in every
  sample, with length, effective length, TPM and read count, annotated with
  the interferon and virus condition
- `construct_features_reads_matrix.csv` — features × samples read counts
- `construct_features_by_group.csv` — mean TPM and detection count per
  condition group, which is the one to read first

Mock samples should be at or near zero. If they are not, the construct FASTA
shares sequence with the mouse transcriptome — almost always a mouse-derived
cDNA or UTR in the backbone. `viral_refs/README.md` covers how to find and fix
that.

## Reading `05_compare_old_vs_new.R` output

`old_vs_new_summary.csv` gives, per contrast, the Pearson correlation of old
against new log2 fold-changes and the number of significant calls kept, lost
and gained. If the interferon contrasts correlate at r > 0.98, the existing ISG
and IDG conclusions hold and only the genes in `changed_calls_*.csv` need
revisiting. Those files are the short list worth actually checking by hand.

## Note on normalisation and the interferon question

This dataset is for murine interferon-stimulated and interferon-downregulated
genes, so one detail is worth keeping in view. Type I and II interferon induce
a large, coordinated transcriptional programme. When a treatment shifts a
substantial share of the transcriptome in one direction, median-of-ratios
normalisation can absorb part of that shift into the size factors and
compress the very fold-changes you are measuring. It is usually modest and
DESeq2's assumption of a mostly-unchanged majority typically survives, but it
is worth confirming rather than assuming: compare the size factors across
groups, and if the interferon-treated samples sit systematically off, consider
estimating size factors from a set of genes known not to respond to interferon.
Worth checking once, before building figures on top of it.
