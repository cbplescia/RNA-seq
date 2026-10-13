#!/usr/bin/env python3
"""Quick screen: do any reads map to your replicon reference, and to what?

Exact 31-mer matching against model_reference/genome.fa, with every hit
attributed to a feature from model_reference/annotation.gtf. Presence/absence
and per-feature read counts - not abundance.

  python3 -I model_ref/01_quick_screen.py

  REF_DIR    directory holding genome.fa and annotation.gtf
             default: $PROJ/model_reference
  SAMPLES    comma-separated sample ids, default 1,2,3,4,5,6,7,8,9
  MAX_READS  reads per mate per sample, default 1000000; 0 = every read
  K, STEP    k-mer length and tiling stride, default 31 and 10
  KEEP_HITS  matching reads saved per sample, default 200; 0 = all
  OUT_DIR    results directory, default $PROJ/model_ref_screen

Why exact matching rather than an aligner: a replicon-only index has no host
sequence to compete, so host reads that merely resemble the construct get
force-assigned to it. An exact 31-mer of construct sequence is diagnostic on
its own, and the mock samples measure the floor empirically.
"""

import gzip
import os
import sys
from collections import Counter, defaultdict

PROJ = os.environ.get("PROJ") or os.path.join(
    os.path.expanduser("~"),
    "Documents/aviti_rnaseq/2025/aviti",
    "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
    "results/Maury_25161")

REF_DIR = os.environ.get("REF_DIR", os.path.join(PROJ, "model_reference"))
SAMPLES = os.environ.get("SAMPLES", "1,2,3,4,5,6,7,8,9").split(",")
MAX_READS = int(os.environ.get("MAX_READS", 1_000_000))
K = int(os.environ.get("K", 31))
STEP = int(os.environ.get("STEP", 10))
KEEP_HITS = int(os.environ.get("KEEP_HITS", 200))
OUT_DIR = os.environ.get("OUT_DIR", os.path.join(PROJ, "model_ref_screen"))

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    recs, name, seq = [], None, []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(seq)))
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.upper())
    if name is not None:
        recs.append((name, "".join(seq)))
    return recs


def parse_gtf(path):
    """Features as (seqname, start0, end0, name, type), 0-based half-open."""
    feats = []
    wanted = {"gene", "transcript", "exon", "CDS", "five_prime_utr",
              "three_prime_utr", "mRNA"}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            ftype = f[2]
            if wanted and ftype not in wanted:
                continue
            try:
                start, end = int(f[3]) - 1, int(f[4])
            except ValueError:
                continue
            attrs = f[8]
            name = None
            for key in ("gene_name", "gene_id", "transcript_id", "ID", "Name",
                        "locus_tag", "product"):
                # GTF style: key "value";   GFF style: key=value;
                for sep in ('"', "="):
                    tag = f"{key} {sep}" if sep == '"' else f"{key}{sep}"
                    i = attrs.find(tag)
                    if i >= 0:
                        rest = attrs[i + len(tag):]
                        end_ch = '"' if sep == '"' else ";"
                        name = rest.split(end_ch)[0].strip()
                        break
                if name:
                    break
            feats.append((f[0], start, end, name or f"{ftype}@{f[3]}", ftype))
    return feats


def main():
    gfa = os.path.join(REF_DIR, "genome.fa")
    gtf = os.path.join(REF_DIR, "annotation.gtf")
    for p in (gfa, gtf):
        if not os.path.exists(p):
            sys.exit(f"not found: {p}\nSet REF_DIR if your reference is elsewhere.")

    recs = read_fasta(gfa)
    if not recs:
        sys.exit(f"no FASTA records in {gfa}")

    print("===== REFERENCE =====")
    total_nt = 0
    for name, seq in recs:
        amb = len(seq) - sum(seq.count(b) for b in "ACGT")
        total_nt += len(seq)
        print(f"  {name:<28} {len(seq):>9,} nt"
              + (f"   ({amb} non-ACGT)" if amb else ""))

    feats = parse_gtf(gtf)
    by_type = Counter(f[4] for f in feats)
    print(f"\n  annotation.gtf: {len(feats)} features  "
          + ", ".join(f"{t}={n}" for t, n in sorted(by_type.items())))
    seqs_in_gtf = {f[0] for f in feats}
    seqs_in_fa = {n for n, _ in recs}
    orphan = seqs_in_gtf - seqs_in_fa
    if orphan:
        print(f"  WARNING: GTF references sequence names absent from genome.fa: "
              f"{sorted(orphan)[:5]}")
        print(f"           genome.fa has: {sorted(seqs_in_fa)[:5]}")
        print("           Feature attribution will be empty until these match.")

    # Does this reference carry the reporter variant that is actually in the
    # libraries? If the replicon was drafted from published sequence, its
    # reporter may be a different codon variant and will not match by k-mer.
    rec_gfp = os.path.join(PROJ, "ref", "construct_gfp_recovered.fa")
    if os.path.exists(rec_gfp):
        gfp = read_fasta(rec_gfp)[0][1]
        probe = gfp[60:150]
        joined = "\n".join(s for _, s in recs)
        present = probe in joined or revcomp(probe) in joined
        print(f"\n  recovered reporter sequence present in this reference: "
              f"{'YES' if present else 'NO'}")
        if not present:
            print("    The reporter in this reference differs from the one in the")
            print("    libraries (or is absent). Reporter-derived reads will not")
            print("    match it by k-mer - that is the variant trap again.")

    # ---- k-mer index over the reference, both strands ----------------------
    pos_of = {}
    for name, seq in recs:
        starts = list(range(0, max(0, len(seq) - K + 1), STEP))
        last = len(seq) - K
        if starts and starts[-1] != last:
            starts.append(last)
        for i in starts:
            km = seq[i:i + K]
            if "N" in km or len(km) < K:
                continue
            pos_of.setdefault(km, (name, i))
            pos_of.setdefault(revcomp(km), (name, i))
    if not pos_of:
        sys.exit("no usable k-mers from the reference")
    print(f"\n  search set: {len(pos_of):,} {K}-mers (both strands, every {STEP} nt)")

    # Self-check: reference sequence must be detected.
    ctrl = recs[0][1][:200]
    if not any(ctrl[i:i + K] in pos_of for i in range(len(ctrl) - K + 1)):
        sys.exit("FAILED self-check: reference sequence not detected by its own "
                 "k-mer set. Not screening.")
    print("  self-check: reference is detected by its own search set")

    # ---- feature lookup ----------------------------------------------------
    feats_by_seq = defaultdict(list)
    for sname, s, e, fname, ftype in feats:
        feats_by_seq[sname].append((s, e, fname, ftype))

    def features_at(sname, pos, length):
        out = set()
        for s, e, fname, ftype in feats_by_seq.get(sname, ()):
            if pos < e and (pos + length) > s:
                out.add(fname)
        return out

    os.makedirs(OUT_DIR, exist_ok=True)
    rows = []
    per_feature = defaultdict(Counter)   # feature -> sample -> reads

    print(f"\n===== SCREEN (samples {','.join(SAMPLES)}) =====")
    if MAX_READS:
        print(f"depth: first {MAX_READS:,} reads per mate")
    else:
        print("depth: every read")

    for s in SAMPLES:
        hits_path = os.path.join(OUT_DIR, f"hits_{s}.txt")
        kept, screened, nhits = [], 0, 0
        placed_rec = Counter()

        for mate in ("R1", "R2"):
            fq = os.path.join(PROJ, "trimmed_fastq", f"{s}_{mate}.trim.fastq.gz")
            if not os.path.exists(fq):
                print(f"  sample {s}: missing {os.path.basename(fq)} - skipped")
                screened = -1
                break
            with gzip.open(fq, "rt") as fh:
                for i, line in enumerate(fh):
                    if i % 4 != 1:
                        continue
                    screened += 1
                    r = line.strip().upper()
                    for j in range(0, len(r) - K + 1):
                        got = pos_of.get(r[j:j + K])
                        if got is None:
                            continue
                        sname, p = got
                        nhits += 1
                        placed_rec[sname] += 1
                        start = max(0, p - j)
                        for fname in features_at(sname, start, len(r)) or {"(unannotated)"}:
                            per_feature[fname][s] += 1
                        if KEEP_HITS == 0 or len(kept) < KEEP_HITS:
                            kept.append(r)
                        break
                    if MAX_READS and screened >= MAX_READS * (1 if mate == "R1" else 2):
                        break

        if screened < 0:
            continue
        per10m = nhits / screened * 1e7 if screened else 0
        print(f"  sample {s:<5} screened {screened:>12,}   hits {nhits:>6}"
              f"  ({per10m:,.2f} per 10M)")
        rows.append((s, screened, nhits, per10m))
        if kept:
            with open(hits_path, "w") as fh:
                fh.write("\n".join(kept) + "\n")

    # ---- output ------------------------------------------------------------
    tsv = os.path.join(OUT_DIR, "model_ref_hits.tsv")
    with open(tsv, "w") as fh:
        fh.write("sample\treads_screened\thits\thits_per_10M\n")
        for s, sc, h, p in rows:
            fh.write(f"{s}\t{sc}\t{h}\t{p:.2f}\n")

    if per_feature:
        print("\n===== READS PER FEATURE =====")
        order = [s for s, _, _, _ in rows]
        width = max(len(f) for f in per_feature) + 2
        print("  " + "feature".ljust(width) + "".join(f"{s:>8}" for s in order))
        for fname, counts in sorted(per_feature.items(),
                                    key=lambda kv: -sum(kv[1].values())):
            print("  " + fname.ljust(width)
                  + "".join(f"{counts.get(s, 0):>8}" for s in order))
        fpath = os.path.join(OUT_DIR, "model_ref_hits_by_feature.tsv")
        with open(fpath, "w") as fh:
            fh.write("feature\t" + "\t".join(order) + "\n")
            for fname, counts in sorted(per_feature.items(),
                                        key=lambda kv: -sum(kv[1].values())):
                fh.write(fname + "\t" + "\t".join(str(counts.get(s, 0))
                                                  for s in order) + "\n")
        print(f"\n  {fpath}")

    print(f"  {tsv}")
    print(f"  {OUT_DIR}/hits_<sample>.txt")

    print("""
How to read this
  Samples 1-3 are mock and never saw the construct. Their count is the floor -
  index hopping, carryover, or sequence shared with the mouse transcriptome.
  Judge 4-9 against that floor rather than against zero.

  A feature row that is high in mock as well as infected samples is almost
  certainly host-homologous sequence in the reference, not construct RNA.
  "(unannotated)" counts reads landing outside any GTF feature - backbone,
  intergenic, or regions the annotation does not cover.""")


if __name__ == "__main__":
    main()
