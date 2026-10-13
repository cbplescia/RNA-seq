# Is there any GFP reporter RNA in these libraries?

A standalone presence/absence screen. It does not touch the index, the
quantification, or anything the interferon analyses depend on, and it does not
require the replicon work to be set up first. Run it, get an answer, decide
later whether a full replicon-aware rebuild is worth the hours.

## Why this is a k-mer screen and not an alignment

Every read is checked for an **exact 31-mer** tiled from your GFP sequence, on
both strands. The logic is that a 31-mer of jellyfish-derived GFP does not
occur in the mouse genome, so a read carrying one came from the reporter. There
is no alignment score to tune, no index that can force a host read onto the one
sequence on offer, and nothing to argue about: the hits are reads you can open
and look at.

That is the opposite trade-off from `01_viral_probe.sh` one level up, which
builds a construct-only Salmon index. That approach gives you numbers, but a
single-sequence index has no host to compete, so host reads that merely
resemble the construct get assigned to it. For the question "is there anything
there at all", exact matching is both stricter and simpler.

The cost is that exact matching is **variant-specific**, and more brutally so
than it sounds. EGFP is codon-optimised relative to wild-type avGFP: the
proteins are nearly identical but a fifth to a quarter of the nucleotides
differ, scattered across third codon positions. A shared 31-mer needs about ten
consecutive untouched codons, which essentially never happens. In a synthetic
test at 66% nucleotide identity with an identical protein, the two variants
shared **zero** exact 31-mers, and nucleotide mode returned zero hits on data
containing 120 planted reads.

So a zero in nucleotide mode means *"no reads from this exact variant"*, not
*"no GFP"*. Two ways out:

- supply the GFP actually in your construct (best: most specific, and fastest)
- run `MODE=aa`, which searches translated reads and survives codon
  differences — see **Variant-tolerant mode** below

## Variant-tolerant mode (`MODE=aa`)

```bash
MODE=aa ./01_gfp_kmer_screen.sh
```

Each read is translated in all six frames and checked for a 9-residue window of
the reference protein. GFP variants differ at a handful of residues even when
their codon usage differs everywhere, so a protein-level search finds one
variant using another as the reference. A 9-mer peptide match is still
specific — the mouse proteome contains no GFP.

Validated: with a reference sharing **zero** 31-mers with the planted reads, nt
mode found 0 and aa mode found 240 of 240, with the negative control still at
zero and full coverage across the record.

The trade-off is speed. It pushes every read through Python at roughly 30,000
reads/s, against `grep`'s streaming speed. A 500k-read-per-mate pass over 33
samples is around 15–20 minutes; full depth would be a day or more, so aa mode
is for the fast pass. If you need full depth, get the construct's own sequence
and use nucleotide mode.

Both modes run a **self-check** before screening: a stretch of the reference is
pushed through the filter and must come back. That catches a broken search —
which otherwise reports every sample as zero, indistinguishable from a real
negative. It does not, and cannot, tell you whether the reference is the right
variant.

## If your reference is a GenBank mRNA record

An NCBI mRNA entry (`...complete cds`) includes the UTRs, so the record is not
the bare CDS and frame 1 is not the coding frame. The validator detects this,
finds the real ORF, confirms the chromophore inside it, and tells you to
extract it:

```bash
python3 -I gfp_tools.py orf GFP.FASTA > ref/gfp.fa
python3 -I gfp_tools.py validate ref/gfp.fa
```

`orf` takes the longest ATG-to-stop ORF on either strand, so it also fixes a
reverse-complemented record.

## Step 1 — get the reporter sequence in place

Expected at `ref/gfp.fa` inside the project directory, or point `GFP_FA`
elsewhere. From a SnapGene, Benchling or Addgene GenBank export:

```bash
cd gfp
./00_make_gfp_fasta.sh ~/path/to/my_construct.gb          # finds a GFP-like CDS
./00_make_gfp_fasta.sh ~/path/to/my_construct.gb eGFP     # or name the feature
```

It reads the `FEATURES` table, handles `complement(...)` and `join(...)`
locations, takes the CDS over a bare annotation, and prints the other
candidates it ignored. If no GFP-like feature is found it lists every feature
name in the file so you can pick one.

If you already have a plain FASTA of the reporter, skip this and set
`GFP_FA=/path/to/it`.

**No construct map at hand?** Fetch the reference EGFP from NCBI on your Mac:

```bash
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=U55762.1&rettype=fasta&retmode=text" \
  > ~/Downloads/pEGFP-N1.fa
```

That is the pEGFP-N1 vector, so you then need to cut out the EGFP CDS rather
than screen with the whole plasmid. Your own map is better — use it if you have
it. A deliberately missing sequence beats a guessed one, which is why there is
no GFP sequence committed to this repo: a wrong reference returns zero hits that
are indistinguishable from a true negative.

## Step 2 — validate the record

Automatic inside the screen, and runnable on its own:

```bash
python3 -I gfp_tools.py validate /path/to/gfp.fa
```

It reports length, translation and the position of the GFP chromophore motif
(Thr-Tyr-Gly around residue 65, Ser-Tyr-Gly in wild-type), which is the
variant-independent signature that the record really is a GFP. It fails on
non-ACGTN characters and on in-frame internal stops, which is what a
reverse-complemented or out-of-frame record looks like.

Orientation and frame warnings do not invalidate the screen — it searches both
strands and ignores frame. They tell you whether the record is the clean CDS
you meant to use.

## Step 3 — run the screen

```bash
cd gfp
./01_gfp_kmer_screen.sh                 # fast pass: first 4M reads per mate
MAX_READS=0 ./01_gfp_kmer_screen.sh     # full depth, every read
```

All 33 samples, R1 and R2 separately. Missing samples are skipped, not fatal.
Output goes to `gfp_screen/` in the project directory:

- `gfp_kmer_hits.tsv` — reads screened, hits per mate, hits per 10M
- `hits_<sample>.txt` — up to 200 matching read sequences per positive sample

The fast pass is for presence/absence. A transcribed reporter in an infected
sample is not subtle and will show up in the first few million reads. Use full
depth when you need to say something about a sample that looked empty, since
absence is the claim that needs the depth behind it.

Runtime is dominated by gzip decompression, so expect a few minutes per sample
at full depth on a laptop.

## Step 4 — check where the hits land

```bash
python3 -I 02_gfp_hit_profile.py /path/to/gfp.fa /path/to/gfp_screen/hits_*.txt
```

Prints a coverage sparkline across the GFP record per sample. Reads from
genuinely transcribed GFP spread across the whole coding sequence. Hits piled
into one narrow window usually mean that window is not GFP-specific after all
— a linker, a restriction site, or a low-complexity stretch — so the script
says so when coverage is under 25%.

This is the step that turns "there are 40 hits" into "there is reporter RNA".

## Recovering the construct's own sequence (`03_recover_construct_gfp.py`)

If `MODE=aa` finds reads but `MODE=nt` finds none, the construct carries a
different codon variant than your reference — and those reads are the
construct's sequence. This anchors each read to the reference through its
translation, piles them up and calls a consensus:

```bash
KEEP_HITS=0 MODE=aa MAX_READS=2000000 ./01_gfp_kmer_screen.sh
python3 -I 03_recover_construct_gfp.py "$PROJ/ref/gfp.fa" "$PROJ"/gfp_screen/hits_*.txt
```

`KEEP_HITS=0` keeps every matching read rather than the first 200, which is
what gives the pileup its depth. Pool all positive samples — they share the
same construct, so they contribute to one consensus.

The report gives positions called, mean depth, nucleotide and amino-acid
identity to the reference, the residue differences, whether the chromophore
survived, and any coverage gaps. Output goes to
`ref/construct_gfp_recovered.fa`.

Validated against ground truth on synthetic data: with a reference at 67%
nucleotide identity and 0.2% per-base sequencing error, it recovered 711 of 717
positions with **zero** errors at called positions. Uncalled positions are the
extreme termini, where reads need 27 nt inside the ORF to be anchored; more
depth fills them.

It assumes no indels relative to the reference protein. Each read is placed
independently and must have its peptide anchors agree, so an indel costs
individual reads rather than corrupting the consensus — it shows up as a run of
disagreement in the per-position report.

This sequence is what goes into `recombinant_model.fa` as the reporter record,
and it is empirically derived from your own libraries rather than from a map
that may not match what is in the freezer.

## Interpreting the table

Samples 1–3 and 10–15 never saw the construct. They are the negative controls,
and whatever they show is the floor: index hopping, carryover on a shared run,
or a reference problem. Judge everything else against that floor.

| Pattern | Reading |
|---|---|
| Everything ~0 | No detectable reporter RNA. Could be real (protein-only VLPs carry no reporter transcript), or the wrong GFP variant. Re-check the validate output first — it is the cheapest explanation to rule out. |
| Only expected samples positive | Real signal. Worth the replicon rebuild for actual numbers. |
| Negatives positive at similar levels | Contamination or index hopping, not biology. Check `Maury_25161_IndexAssignment.csv` before concluding anything from the positives. |
| Positives, but coverage concentrated | Suspect. Read the hit sequences before believing them. |

One biological point: if the EVLP and MVLP samples received protein-only
virus-like particles, zero reporter RNA is the **correct** result, not a
failure. If they were transcription- and replication-competent trVLPs, the
reporter should be detectable. Knowing which system you used is what makes this
screen a test rather than a number.

## What this deliberately does not do

No abundance estimates, no normalisation, no differential expression. Hits per
10M reads is a screening statistic, not a quantity — a transcript's
abundance depends on effective length and on competing targets that a k-mer
count knows nothing about. Real numbers come from the replicon rebuild one
level up, where the full transcriptome and genome decoys compete for every
read.

## When you come back to the replicon

Nothing here needs undoing. The pipeline in the parent directory is already
staged for it:

1. Build `ref/recombinant_model.fa` with one record per transcribed feature —
   the replicon/minigenome, the viral genes, and the reporter. See
   `../viral_refs/README.md`.
2. The GFP record you validated here goes straight into that file.
3. Then `../00_preflight.sh` → `../02_build_reference.sh` →
   `../03_salmon_quant.sh` → `../04_tximport_deseq2.R`.

`config.sh` already carries `GFP_FA` next to `VIRAL_FA`, so the screen and the
rebuild stay independent and neither blocks the other.
