#!/usr/bin/env python3
"""Helpers for the GFP screen: validate a GFP FASTA, and tile it into k-mers.

Two subcommands:

  validate <gfp.fa>          sanity-check the records before you trust a result
  kmers <gfp.fa> [k] [step]  write search patterns (both strands) to stdout

The validator exists because the failure mode of this whole screen is silent:
a truncated, reverse-complemented or simply wrong sequence produces zero hits
that look exactly like a real negative. It checks each record for the things a
functional GFP coding sequence must have, so a bad input fails loudly here
instead of quietly later.
"""

import sys

CODONS = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    recs, name, seq = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(seq)))
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.upper().replace(" ", ""))
    if name is not None:
        recs.append((name, "".join(seq)))
    return recs


def translate(seq):
    return "".join(CODONS.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


def validate(path):
    recs = read_fasta(path)
    if not recs:
        print(f"FAIL: no FASTA records found in {path}")
        return 1

    problems, warnings = [], []
    print(f"{len(recs)} record(s) in {path}\n")

    for name, seq in recs:
        print(f"--- {name}")
        print(f"    length      {len(seq)} nt")

        bad = set(seq) - set("ACGTN")
        if bad:
            problems.append(f"{name}: non-ACGTN characters {sorted(bad)}")
        n_frac = seq.count("N") / max(len(seq), 1)
        if n_frac > 0.01:
            warnings.append(f"{name}: {n_frac:.1%} Ns")

        # A GFP-family coding sequence is ~717 nt (239 aa) plus a stop codon.
        if not 600 <= len(seq) <= 900:
            warnings.append(
                f"{name}: length {len(seq)} is outside the ~700-750 nt expected "
                "for a GFP CDS - fine if this is a fusion or a fragment, "
                "a problem if you meant the bare CDS")

        aa = translate(seq)
        internal_stops = aa[:-1].count("*")
        print(f"    translates  {len(aa)} aa, starts {aa[:12]!r}")

        if not seq.startswith("ATG"):
            warnings.append(f"{name}: does not start with ATG")
        if internal_stops:
            problems.append(
                f"{name}: {internal_stops} internal stop codon(s) in frame 1 - "
                "the record is probably reverse-complemented, out of frame, "
                "or not a clean CDS")

        # The GFP chromophore is a Thr-Tyr-Gly tripeptide around residue 65
        # (Ser-Tyr-Gly in wild-type avGFP). Its presence is a strong, sequence-
        # variant-independent signal that this really is a GFP.
        window = aa[55:80]
        motif = next((m for m in ("TYG", "SYG", "AYG", "CYG") if m in window), None)
        if motif:
            pos = aa.index(motif, 55) + 1
            print(f"    chromophore {motif} at residue {pos}  (looks like GFP)")
        else:
            warnings.append(
                f"{name}: no GFP chromophore motif (TYG/SYG) near residue 65 - "
                "check this is really the fluorophore CDS and is in frame")

    print()
    for w in warnings:
        print(f"WARN: {w}")
    for p in problems:
        print(f"FAIL: {p}")

    if problems:
        print("\nFix the FAIL lines before screening - a bad reference gives "
              "zero hits that look just like a true negative.")
        return 1
    print("\nOK to screen." + ("  Read the WARN lines first." if warnings else ""))
    if warnings:
        print("Orientation and frame do not affect detection - the screen "
              "searches both strands and ignores frame. They only tell you "
              "whether the record is the clean CDS you meant to use.")
    return 0


def kmers(path, k=31, step=10):
    recs = read_fasta(path)
    out = []
    seen = set()
    for _, seq in recs:
        if len(seq) < k:
            continue
        starts = list(range(0, len(seq) - k + 1, step))
        last = len(seq) - k
        if starts and starts[-1] != last:
            starts.append(last)
        for i in starts:
            kmer = seq[i:i + k]
            if "N" in kmer:
                continue
            for s in (kmer, revcomp(kmer)):
                if s not in seen:
                    seen.add(s)
                    out.append(s)
    if not out:
        sys.exit(f"no usable {k}-mers from {path}")
    print("\n".join(out))


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    cmd, fa = sys.argv[1], sys.argv[2]
    if cmd == "validate":
        sys.exit(validate(fa))
    elif cmd == "kmers":
        k = int(sys.argv[3]) if len(sys.argv) > 3 else 31
        step = int(sys.argv[4]) if len(sys.argv) > 4 else 10
        kmers(fa, k, step)
    else:
        sys.exit(__doc__)
