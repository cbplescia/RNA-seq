#!/usr/bin/env python3
"""Where along GFP do the matching reads land?

Run after 01_gfp_kmer_screen.sh, on the hits_<sample>.txt files it leaves
behind. Takes the GFP reference and one or more hit files, and prints a
coverage profile across the record.

This is the check that separates real reporter RNA from an artifact. Reads from
a transcribed GFP spread across the whole coding sequence. Hits that pile onto
one narrow window usually mean that window is not GFP-specific after all - a
linker, a restriction site, a low-complexity stretch, or a region shared with
something else in the library.

  python3 -I 02_gfp_hit_profile.py ref/gfp.fa gfp_screen/hits_*.txt
"""

import sys
from collections import Counter

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from gfp_tools import CODONS, longest_orf, read_fasta, revcomp, translate  # noqa: E402

BARS = " ▁▂▃▄▅▆▇█"


def place_by_peptide(read, aa_pos_of, pepk):
    """Nucleotide offset of a read placed through its translation.

    Used for reads found in translated mode: they encode the reference protein
    but, being a different codon variant, share no exact nucleotide k-mer with
    it. Returns an offset into the reference CDS, or None.
    """
    for strand in (read, revcomp(read)):
        for frame in range(3):
            s = strand[frame:]
            n = (len(s) // 3) * 3
            if n < pepk * 3:
                continue
            aa = "".join(CODONS.get(s[i:i + 3], "X") for i in range(0, n, 3))
            for i in range(len(aa) - pepk + 1):
                p = aa_pos_of.get(aa[i:i + pepk])
                if p is not None:
                    # Residue index -> nucleotide offset, minus the read's own
                    # lead-in so the footprint starts where the read starts.
                    return max(0, (p - i) * 3)
    return None


def profile(ref_path, hit_paths, k=31):
    recs = read_fasta(ref_path)
    if not recs:
        sys.exit(f"no records in {ref_path}")
    name, ref = recs[0][0], recs[0][1]

    # Every k-mer of the reference, forward and reverse, mapped to its position.
    pos_of = {}
    for i in range(len(ref) - k + 1):
        kmer = ref[i:i + k]
        pos_of.setdefault(kmer, i)
        pos_of.setdefault(revcomp(kmer), i)

    # Peptide fallback, for hits found in translated mode against a different
    # GFP variant: those reads share no exact nucleotide k-mer with the
    # reference, so they can only be placed through the protein.
    ref_cds = ref
    if translate(ref)[:-1].count("*"):
        found = longest_orf(ref)
        if found:
            ref_cds = found[0]
    ref_aa = translate(ref_cds).rstrip("*")
    PEPK = 9
    aa_pos_of = {ref_aa[i:i + PEPK]: i for i in range(len(ref_aa) - PEPK + 1)}

    print(f"reference: {name}  ({len(ref)} nt)\n")

    for path in hit_paths:
        try:
            reads = [l.strip().upper() for l in open(path) if l.strip()]
        except OSError as e:
            print(f"{path}: {e}")
            continue

        cov = Counter()
        placed = unplaced = 0
        placed_by_peptide = 0
        for r in reads:
            found = None
            for i in range(len(r) - k + 1):
                p = pos_of.get(r[i:i + k])
                if p is not None:
                    found = p
                    break
            if found is None:
                found = place_by_peptide(r, aa_pos_of, PEPK)
                if found is not None:
                    placed_by_peptide += 1
            if found is None:
                unplaced += 1
                continue
            placed += 1
            # Credit the read's whole footprint, clipped to the reference.
            start = max(0, found - 0)
            for j in range(start, min(len(ref), start + len(r))):
                cov[j] += 1

        label = path.rsplit("/", 1)[-1]
        extra = f", {placed_by_peptide} via protein" if placed_by_peptide else ""
        print(f"--- {label}: {len(reads)} reads, {placed} placed{extra}, "
              f"{unplaced} unplaced")
        if not placed:
            print("    nothing to plot\n")
            continue

        # 60 bins across the record
        nbins, peak = 60, 0
        bins = [0] * nbins
        width = len(ref) / nbins
        for p, c in cov.items():
            bins[min(nbins - 1, int(p / width))] += c
        peak = max(bins)
        line = "".join(BARS[min(8, int(b / peak * 8 + 0.5))] if peak else " " for b in bins)
        covered = sum(1 for b in bins if b) / nbins

        print(f"    5'|{line}|3'")
        print(f"    {covered:.0%} of the record covered, peak depth {peak}")
        if covered < 0.25:
            print("    NOTE: hits are concentrated in a small part of the record - "
                  "inspect those reads before believing the signal")
        print()


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    profile(sys.argv[1], sys.argv[2:])
