#!/usr/bin/env python3
"""Turn a replicon genome + GTF into Salmon quantification targets.

Salmon quantifies transcripts, not genomes. Handing it an 18 kb genomic contig
as a single record collapses every gene into one number and computes effective
length against a molecule that does not exist. This extracts one record per
transcript, exon-joined and strand-aware, so each viral gene is counted
separately.

  python3 -I model_ref/02_make_targets.py

  REF_DIR    directory holding the replicon FASTA and annotation
  GENOME_FA  the FASTA, if the directory holds more than one
  GTF        the annotation, if the directory holds more than one
  PREFIX     prepended to every record name, default "model"
  NAME_MAP   TSV of <gtf_name><tab><new_name> to rename features
  EXCLUDE    comma-separated gtf names to leave out
  REPORTER   reporter FASTA to append, default ref/construct_gfp_recovered.fa
             (set REPORTER=none to omit)
  OUT        output FASTA, default ref/recombinant_model.fa
  FEATURE    GTF feature type to join, default exon (try CDS if exon is absent)

Overlapping transcripts are reported rather than silently emitted: near-
identical records compete for the same reads, so Salmon splits them and
neither number means much on its own.
"""

import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from importlib.machinery import SourceFileLoader  # noqa: E402

_screen = SourceFileLoader(
    "screen", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "01_quick_screen.py")).load_module()

PROJ = _screen.PROJ
REF_DIR = _screen.REF_DIR
read_fasta = _screen.read_fasta
revcomp = _screen.revcomp
find_one = _screen.find_one

PREFIX = os.environ.get("PREFIX", "model")
FEATURE = os.environ.get("FEATURE", "exon")
EXCLUDE = {x.strip() for x in os.environ.get("EXCLUDE", "").split(",") if x.strip()}
OUT = os.environ.get("OUT", os.path.join(PROJ, "ref", "recombinant_model.fa"))
REPORTER = os.environ.get(
    "REPORTER", os.path.join(PROJ, "ref", "construct_gfp_recovered.fa"))


def attr(attrs, key):
    for sep, end in (('"', '"'), ("=", ";")):
        tag = f"{key} {sep}" if sep == '"' else f"{key}{sep}"
        i = attrs.find(tag)
        if i >= 0:
            return attrs[i + len(tag):].split(end)[0].strip()
    return None


def main():
    gfa = find_one(REF_DIR, _screen.FASTA_EXT, "GENOME_FA", "FASTA file")
    gtf = find_one(REF_DIR, _screen.GTF_EXT, "GTF", "GTF/GFF annotation")
    print(f"FASTA {gfa}\nGTF   {gtf}\n")

    seqs = dict(read_fasta(gfa))

    # Group the chosen feature type by transcript.
    tx = defaultdict(lambda: {"exons": [], "strand": "+", "seqname": None,
                              "gene": None})
    seen_types = set()
    with open(gtf) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            seen_types.add(f[2])
            if f[2] != FEATURE:
                continue
            tid = (attr(f[8], "transcript_id") or attr(f[8], "ID")
                   or attr(f[8], "gene_id") or f"{f[0]}:{f[3]}")
            gene = (attr(f[8], "gene_name") or attr(f[8], "gene_id") or tid)
            rec = tx[tid]
            rec["exons"].append((int(f[3]) - 1, int(f[4])))
            rec["strand"] = f[6] if f[6] in "+-" else "+"
            rec["seqname"] = f[0]
            rec["gene"] = gene

    if not tx:
        sys.exit(f"no '{FEATURE}' features in the GTF.\n"
                 f"  types present: {', '.join(sorted(seen_types))}\n"
                 f"  Set FEATURE to one of those.")

    name_map = {}
    nm = os.environ.get("NAME_MAP")
    if nm:
        with open(nm) as fh:
            for line in fh:
                p = line.split()
                if len(p) >= 2:
                    name_map[p[0]] = p[1]

    records, skipped = [], []
    for tid, rec in tx.items():
        gene = rec["gene"]
        key = name_map.get(tid) or name_map.get(gene) or gene
        if gene in EXCLUDE or tid in EXCLUDE or key in EXCLUDE:
            skipped.append(gene)
            continue
        if rec["seqname"] not in seqs:
            sys.exit(f"GTF sequence name {rec['seqname']!r} is not in the FASTA "
                     f"(which has {list(seqs)[:4]}). Fix the names first.")
        g = seqs[rec["seqname"]]
        exons = sorted(rec["exons"])
        seq = "".join(g[s:e] for s, e in exons)
        if rec["strand"] == "-":
            seq = revcomp(seq)
        safe = "".join(c if c.isalnum() or c in "._-" else "_" for c in str(key))
        records.append({
            "name": f"{PREFIX}_{safe}", "seq": seq, "gene": gene,
            "span": (exons[0][0] + 1, exons[-1][1]), "n_exons": len(exons),
            "strand": rec["strand"], "seqname": rec["seqname"],
        })

    records.sort(key=lambda r: r["span"])

    print(f"{len(records)} target(s) from '{FEATURE}' features"
          + (f"; excluded {', '.join(skipped)}" if skipped else ""))
    print(f"\n  {'record':<28}{'length':>8}{'exons':>7}{'strand':>8}   span")
    for r in records:
        print(f"  {r['name']:<28}{len(r['seq']):>8}{r['n_exons']:>7}"
              f"{r['strand']:>8}   {r['span'][0]}-{r['span'][1]}")

    # Duplicate names would be silently dropped by salmon index.
    names = [r["name"] for r in records]
    dups = {n for n in names if names.count(n) > 1}
    if dups:
        sys.exit(f"\nduplicate record names: {sorted(dups)}\n"
                 "  Salmon will refuse or drop these. Use NAME_MAP to "
                 "disambiguate, or FEATURE=CDS / a different id attribute.")

    # Overlap report: near-identical targets split reads between themselves.
    overlaps = []
    for i in range(len(records)):
        for j in range(i + 1, len(records)):
            a, b = records[i], records[j]
            if a["seqname"] != b["seqname"]:
                continue
            lo = max(a["span"][0], b["span"][0])
            hi = min(a["span"][1], b["span"][1])
            if hi > lo:
                frac = (hi - lo) / min(a["span"][1] - a["span"][0],
                                       b["span"][1] - b["span"][0])
                overlaps.append((a["name"], b["name"], hi - lo, frac))
    if overlaps:
        print("\nOVERLAPPING TARGETS - these compete for the same reads, so "
              "Salmon\nsplits them and neither count stands alone:")
        for n1, n2, nt, frac in sorted(overlaps, key=lambda o: -o[3]):
            print(f"  {n1} / {n2}: {nt} nt ({frac:.0%} of the shorter)")
        print("  Keep them only if you specifically want to resolve those "
              "products;\n  otherwise EXCLUDE all but one and quantify the "
              "locus once.")

    if REPORTER and REPORTER.lower() != "none" and os.path.exists(REPORTER):
        for name, seq in read_fasta(REPORTER):
            rn = "reporter_EmeraldGFP" if "GFP" in name.upper() else f"reporter_{name}"
            records.append({"name": rn, "seq": seq, "gene": rn,
                            "span": (0, 0), "n_exons": 1, "strand": "+",
                            "seqname": "-"})
            print(f"\nappended reporter: {rn} ({len(seq)} nt) from "
                  f"{os.path.basename(REPORTER)}")
    elif REPORTER and REPORTER.lower() != "none":
        print(f"\nNOTE: no reporter FASTA at {REPORTER} - nothing appended.")

    os.makedirs(os.path.dirname(OUT) or ".", exist_ok=True)
    with open(OUT, "w") as fh:
        for r in records:
            fh.write(f">{r['name']}\n")
            for i in range(0, len(r["seq"]), 60):
                fh.write(r["seq"][i:i + 60] + "\n")

    total = sum(len(r["seq"]) for r in records)
    print(f"\nwrote {OUT}")
    print(f"  {len(records)} records, {total:,} nt total")
    print("\nThese become the construct targets in the rebuild. Check the names "
          "above:\nthey are what appear in your count matrix.")


if __name__ == "__main__":
    main()
