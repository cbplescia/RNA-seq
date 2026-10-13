#!/usr/bin/env python3
"""Protein-level read filter: six-frame translate each read, keep the ones
carrying a peptide window from the reference.

Reads FASTQ sequence lines on stdin, writes matching lines to stdout - a drop-in
replacement for `grep -F -f <nucleotide k-mers>` in 01_gfp_kmer_screen.sh.

Why this exists: nucleotide k-mer matching is variant-specific. EGFP is codon-
optimised relative to wild-type avGFP, so the two share almost no exact 31-mers
while encoding nearly the same protein. Searching translated reads finds a GFP
using a different GFP variant as the reference.

  usage: aa_filter.py peptides.txt < seqlines > matching_seqlines
"""

import os
import sys

# `python3 -I` deliberately leaves the script's own directory off sys.path, so
# the sibling import has to be made explicit.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from gfp_tools import CODONS, revcomp  # noqa: E402

STOP = "*"


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)

    with open(sys.argv[1]) as fh:
        peps = {l.strip() for l in fh if l.strip()}
    if not peps:
        sys.exit("no peptides to search for")
    k = len(next(iter(peps)))
    if any(len(p) != k for p in peps):
        sys.exit("peptide patterns must all be the same length")

    get = CODONS.get
    out = sys.stdout.write

    for line in sys.stdin:
        seq = line.strip().upper()
        if len(seq) < k * 3:
            continue
        hit = False
        for strand in (seq, revcomp(seq)):
            for frame in range(3):
                s = strand[frame:]
                n = (len(s) // 3) * 3
                if n < k * 3:
                    continue
                aa = "".join(get(s[i:i + 3], "X") for i in range(0, n, 3))
                # Stop codons split the read into translatable stretches; a real
                # coding read has at least one stretch of k residues.
                for piece in aa.replace("X", STOP).split(STOP):
                    if len(piece) < k:
                        continue
                    for i in range(len(piece) - k + 1):
                        if piece[i:i + k] in peps:
                            hit = True
                            break
                    if hit:
                        break
                if hit:
                    break
            if hit:
                break
        if hit:
            out(seq + "\n")


if __name__ == "__main__":
    main()
