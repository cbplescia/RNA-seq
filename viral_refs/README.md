# Building the recombinant construct FASTA

This is the one input I cannot generate for you. Everything else in the
rebuild is derived from files you already have.

Put the finished file at `ref/recombinant_model.fa` inside the project
directory (or point `VIRAL_FA` somewhere else).

## The rule that matters: one record per transcribed feature, not one per plasmid

It is tempting to drop in the whole plasmid as a single sequence. Don't. A
single record gives you one number that mixes everything the construct
transcribes, and Salmon's effective-length correction will be computed against
a length that no real transcript has. Split it into the units that are actually
transcribed and that you would want to report separately, for example:

```
>EBOV_GP
ATG...
>EBOV_VP40
ATG...
>EBOV_NP
ATG...
>reporter_eGFP
ATG...
>selection_PuroR
ATG...
```

Use short names with no spaces, no `|`, and no characters that need quoting.
The name is what appears in your count matrix.

## What to include

Anything the construct can transcribe in the transfected or infected cell:

- the viral genes themselves (GP, VP40, NP, VP35, VP30, VP24, L, as applicable)
- any reporter (eGFP, mCherry, luciferase)
- any selection marker (PuroR, NeoR, BlastR)
- for a minigenome or trVLP system, the minigenome transcript as its own record
- the vector backbone elements only if they end up in mature RNA (WPRE, a
  retained intron, a non-mammalian polyA signal)

## What to leave out, and why it matters

**Any mouse-derived sequence in the backbone.** If the construct carries a
mouse cDNA, a mouse promoter with transcribed 5' UTR, or a mouse polyA signal,
those regions are identical to host transcripts. Reads from the endogenous
mouse gene will be split between the host transcript and your construct
record, which quietly deflates the host gene and inflates the construct — in
*every* sample, mock included. Either trim those regions out of the record or
keep them and remember that the feature is not interpretable on its own.

Vector regions that are never transcribed (bacterial ori, AmpR, the plasmid
backbone outside the expression cassette) are harmless but add noise; leaving
them out keeps the output cleaner.

**Promoter sequence itself** (CMV, EF1a) is upstream of the transcription start
and should not be in the record.

## Getting the sequences out of a GenBank or SnapGene file

If you have `.gb` / `.gbk` files, this extracts each CDS as its own record:

```bash
python3 -I - <<'PY'
from Bio import SeqIO
import sys
recs = []
for path in ["construct1.gb", "construct2.gb"]:          # edit
    gb = SeqIO.read(path, "genbank")
    for f in gb.features:
        if f.type != "CDS":
            continue
        name = (f.qualifiers.get("gene") or f.qualifiers.get("label")
                or f.qualifiers.get("product") or ["unnamed"])[0]
        name = name.replace(" ", "_").replace("|", "_")
        print(f">{name}")
        print(str(f.extract(gb.seq)))
PY
```

Redirect that to `ref/recombinant_model.fa`, then read the headers back and
check they are the features you expect.

SnapGene can also export directly: File → Export → FASTA, with features as
separate sequences.

## Checking it before the full rebuild

```bash
./00_preflight.sh      # confirms the file is found and lists its records
./01_viral_probe.sh    # minutes: is there construct RNA in the libraries at all?
```

The probe quantifies a few samples against the construct alone. Mock should be
near zero. If mock looks like the infected samples, you have host-homologous
sequence in a record — find it and fix it before building the real index.

## One biological check worth doing first

If your "VLP" samples received protein-only virus-like particles, there may be
no construct RNA in them at all, and a count of zero is the correct answer
rather than a bug. If they are transcription- and replication-competent VLPs
(trVLPs), the minigenome is present and should be detectable. Knowing which
system you used tells you what the probe output *should* look like, which is
what makes it a useful test.
