#!/usr/bin/env python3
"""Helpers for the GFP screen: validate a reference, extract its ORF, and tile
it into nucleotide or peptide search patterns.

  validate <fa>              sanity-check a record before trusting a result
  orf <fa>                   extract the longest ORF as a clean CDS FASTA
  kmers <fa> [k] [step]      nucleotide patterns, both strands (exact search)
  peptides <fa> [k]          peptide patterns (variant-tolerant search)

The validator exists because the failure mode of this screen is silent: a
truncated, reverse-complemented, UTR-padded or simply wrong-variant reference
produces zero hits that look exactly like a real negative.
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

# The GFP chromophore is a Thr-Tyr-Gly tripeptide around residue 65
# (Ser-Tyr-Gly in wild-type avGFP). Finding it is variant-independent evidence
# that a record really is a GFP, whatever the codon usage.
CHROMO = ("TYG", "SYG", "AYG", "CYG", "GYG")


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


def find_chromophore(aa, lo=40, hi=95):
    """Position (1-based) of the chromophore motif, or None."""
    window = aa[lo:hi]
    for m in CHROMO:
        if m in window:
            return aa.index(m, lo) + 1, m
    return None


def longest_orf(seq):
    """Longest ATG->stop ORF over both strands.

    Returns (cds_including_stop, strand, start_0based) or None. An mRNA record
    from GenBank carries UTRs, so the CDS rarely begins at position 1.
    """
    best = None
    for strand, s in (("+", seq), ("-", revcomp(seq))):
        for frame in range(3):
            i = frame
            while i < len(s) - 2:
                if s[i:i + 3] != "ATG":
                    i += 3
                    continue
                j = i
                while j < len(s) - 2:
                    if CODONS.get(s[j:j + 3]) == "*":
                        cand = s[i:j + 3]
                        if best is None or len(cand) > len(best[0]):
                            best = (cand, strand, i)
                        break
                    j += 3
                i += 3
    return best


def validate(path):
    recs = read_fasta(path)
    if not recs:
        print(f"FAIL: no FASTA records found in {path}")
        return 1

    problems, warnings, hints = [], [], []
    print(f"{len(recs)} record(s) in {path}\n")

    for name, seq in recs:
        print(f"--- {name}")
        print(f"    length      {len(seq)} nt")

        bad = set(seq) - set("ACGTN")
        if bad:
            problems.append(f"{name}: non-ACGTN characters {sorted(bad)}")
            continue
        if seq.count("N") / max(len(seq), 1) > 0.01:
            warnings.append(f"{name}: over 1% Ns")

        aa = translate(seq)
        frame1_clean = (seq.startswith("ATG")
                        and aa[:-1].count("*") == 0
                        and find_chromophore(aa) is not None)

        if frame1_clean:
            pos, motif = find_chromophore(aa)
            print(f"    translates  {len(aa)} aa, starts {aa[:12]!r}")
            print(f"    chromophore {motif} at residue {pos}  (looks like GFP)")
            if not 600 <= len(seq) <= 900:
                warnings.append(
                    f"{name}: length {len(seq)} is outside the ~700-750 nt of a "
                    "bare GFP CDS - fine for a fusion, not for the CDS alone")
            continue

        # Frame 1 is not a clean CDS. Before complaining, look for a real ORF -
        # a GenBank mRNA record includes UTRs, so this is the normal case.
        print(f"    frame 1     {len(aa)} aa, {aa[:-1].count('*')} internal stops"
              f" -> not the coding frame")
        orf = longest_orf(seq)
        if orf:
            cds, strand, start = orf
            caa = translate(cds)
            chrom = find_chromophore(caa)
            print(f"    longest ORF {len(cds)} nt ({len(caa) - 1} aa) on the "
                  f"{strand} strand at position {start + 1}")
            if chrom:
                pos, motif = chrom
                print(f"    chromophore {motif} at residue {pos} of that ORF"
                      f"  (it is a GFP - just not trimmed to the CDS)")
                problems.append(
                    f"{name}: the record carries UTR or flanking sequence, so "
                    "it is not the bare CDS")
                hints.append(
                    "This looks like a GenBank mRNA record. Extract the CDS:\n"
                    f"    python3 -I gfp_tools.py orf {path} > ref/gfp.fa")
            else:
                problems.append(
                    f"{name}: longest ORF carries no GFP chromophore motif - "
                    "this may not be a GFP at all")
        else:
            problems.append(f"{name}: no ATG-to-stop ORF found in either strand")

    print()
    for w in warnings:
        print(f"WARN: {w}")
    for p in problems:
        print(f"FAIL: {p}")
    for h in hints:
        print(f"\n{h}")

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


def orf(path):
    recs = read_fasta(path)
    if not recs:
        sys.exit(f"no records in {path}")
    name, seq = recs[0]
    found = longest_orf(seq)
    if not found:
        sys.exit(f"{name}: no ATG-to-stop ORF found")
    cds, strand, start = found
    aa = translate(cds)
    chrom = find_chromophore(aa)
    label = f"{name}_CDS"
    print(f">{label}")
    for i in range(0, len(cds), 60):
        print(cds[i:i + 60])
    msg = (f"extracted {len(cds)} nt ({len(aa) - 1} aa) from {strand} strand "
           f"position {start + 1}-{start + len(cds)}")
    print(msg, file=sys.stderr)
    if chrom:
        print(f"chromophore {chrom[1]} at residue {chrom[0]} - confirmed GFP",
              file=sys.stderr)
    else:
        print("WARNING: no chromophore motif in the extracted ORF", file=sys.stderr)


def kmers(path, k=31, step=10):
    recs = read_fasta(path)
    out, seen = [], set()
    for _, seq in recs:
        if len(seq) < k:
            continue
        starts = list(range(0, len(seq) - k + 1, step))
        if starts and starts[-1] != len(seq) - k:
            starts.append(len(seq) - k)
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


def peptides(path, k=9):
    """Every k-aa window of the reference protein.

    Peptide windows survive codon optimisation: EGFP and wild-type avGFP differ
    at hundreds of nucleotides but only a handful of residues, so a protein-level
    search finds one using the other as reference. Nucleotide k-mers cannot.
    """
    recs = read_fasta(path)
    out, seen = [], set()
    for _, seq in recs:
        aa = translate(seq)
        # If the record is not a clean CDS, use its longest ORF instead.
        if aa[:-1].count("*"):
            found = longest_orf(seq)
            if found:
                aa = translate(found[0])
        aa = aa.rstrip("*")
        for i in range(len(aa) - k + 1):
            pep = aa[i:i + k]
            if "X" in pep or "*" in pep:
                continue
            if pep not in seen:
                seen.add(pep)
                out.append(pep)
    if not out:
        sys.exit(f"no usable {k}-aa peptides from {path}")
    print("\n".join(out))


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    cmd, fa = sys.argv[1], sys.argv[2]
    if cmd == "validate":
        sys.exit(validate(fa))
    elif cmd == "orf":
        orf(fa)
    elif cmd == "kmers":
        kmers(fa,
              int(sys.argv[3]) if len(sys.argv) > 3 else 31,
              int(sys.argv[4]) if len(sys.argv) > 4 else 10)
    elif cmd == "peptides":
        peptides(fa, int(sys.argv[3]) if len(sys.argv) > 3 else 9)
    else:
        sys.exit(__doc__)
