#!/usr/bin/env python3
"""Predict lysine modifications on one or more protein sequences and compare two proteins (e.g. eGP vs mGP).

IMPORTANT - what this is and is not.  MultiLyGAN (https://github.com/Lab-Xu/MultiLyGAN) is a research pipeline (MATLAB 2016a + TensorFlow 1.6 feature / GAN code) that publishes its TRAINING DATA
(Data/*.mat: 17-residue windows centred on lysines that are / are not acetylated, glycated, malonylated, methylated, succinylated, SUMOylated or ubiquitinated) but no pretrained model and no
command-line predictor. This script therefore does NOT run MultiLyGAN itself. It downloads MultiLyGAN's published training data and trains a Random Forest per modification (the classifier family
MultiLyGAN uses) on sequence-window features (one-hot + k-spaced amino-acid pair composition, k = 0 and 1), then scores every lysine of your proteins. No CGAN augmentation, no structure /
AAindex features, so expect it to be weaker than the published pipeline. Treat the output as a ranked hypothesis list and confirm with an independent predictor (MusiteDeep, GPS) and, ideally, MS / K->R mutagenesis.

What it does
  1. downloads the seven training files into ./multilygan_data/ (once), trains one Random Forest per modification on a balanced sample; the evaluation split is by PROTEIN ID (no protein is in
     both training and test), so the reported AUC is not inflated by homologous windows; the call threshold is the score that gives FPR (default 5%) on held-out non-modified lysines
  2. scores every lysine (17-mer window, X padded at the termini) in each FASTA
  3. with TWO FASTAs: aligns the proteins (global alignment, BLOSUM62), matches lysines column by column, and compares predictions: conserved lysines predicted in both / only one / neither, lysines
     present in only one protein, and the per-modification score difference
Outputs (folder multilygan_style_output/): training_report.tsv; <name>_ptm_predictions.tsv (each lysine, probability and call for each modification); with two proteins also
  comparison_aligned_lysines.tsv, comparison_summary.tsv, and the figures sites_per_modification.png, ptm_map_<A>_vs_<B>.png, aligned_lysine_scatter.png, ptm_difference_map.png (each also as .pdf)
Usage:   python 99_lysine_ptm_predict_compare.py eGP.FASTA mGP.FASTA
Options (environment variables):  PTMS=Ace,Sumo,Ubiq (subset: Ace Glyca Malon Meth Succ Sumo Ubiq)  N_TRAIN=15000 (positives per modification, same number of negatives)  N_TREES=200  FPR=0.05  SEED=1
Needs: python 3, numpy, scipy, scikit-learn, joblib, matplotlib   (pip install numpy scipy scikit-learn joblib matplotlib)
The probabilities come from a model trained on a 50:50 sample, so they are scores, not real-world probabilities - use the calls and the ranking.
"""
import os, sys, time, urllib.request, itertools
import numpy as np

AA = "ACDEFGHIKLMNPQRSTVWYX"
AIDX = {c: i for i, c in enumerate(AA)}
PTM_NAMES = {"Ace": "acetylation", "Glyca": "glycation", "Malon": "malonylation", "Meth": "methylation", "Succ": "succinylation", "Sumo": "SUMOylation", "Ubiq": "ubiquitination"}
BASE = "https://raw.githubusercontent.com/Lab-Xu/MultiLyGAN/main/Data/{}8.mat"
PTMS = [p.strip() for p in os.environ.get("PTMS", ",".join(PTM_NAMES)).split(",") if p.strip()]
N_TRAIN = int(os.environ.get("N_TRAIN", "15000")); N_TREES = int(os.environ.get("N_TREES", "200")); FPR = float(os.environ.get("FPR", "0.05")); SEED = int(os.environ.get("SEED", "1"))
DATA = "multilygan_data"; OUT = "multilygan_style_output"; W = 8   # 17-mer window = 8 + K + 8
os.makedirs(DATA, exist_ok=True); os.makedirs(OUT, exist_ok=True)
t0 = time.time()
def say(*a): print(f"[{(time.time() - t0) / 60:5.1f} min]", *a, flush=True)

# ---------------------------------------------------------------- sequences ----
def read_fasta(path):
    recs, name, seq = [], None, []
    for line in open(path):
        line = line.strip()
        if not line: continue
        if line.startswith(">"):
            if name is not None: recs.append((name, "".join(seq)))
            name, seq = line[1:].strip(), []
        else: seq.append("".join(c for c in line.upper() if c.isalpha()))
    if name is not None: recs.append((name, "".join(seq)))
    return recs

def windows(seq):
    out = []
    for i, c in enumerate(seq):
        if c == "K":
            a, b = i - W, i + W + 1
            out.append((i + 1, "X" * max(0, -a) + seq[max(0, a):b] + "X" * max(0, b - len(seq))))
    return out

def encode(wins):
    n = len(wins)
    arr = np.array([[AIDX.get(c, 20) for c in w] for w in wins], dtype=np.int16)            # n x 17
    onehot = np.eye(21, dtype=np.float32)[arr].reshape(n, -1)
    feats = [onehot]
    for k in (0, 1):                                                                         # CKSAAP: pairs of residues k apart
        a, b = arr[:, : 17 - k - 1], arr[:, k + 1:]
        idx = (np.arange(n)[:, None] * 441 + a * 21 + b).ravel()
        cnt = np.bincount(idx, minlength=n * 441).reshape(n, 441).astype(np.float32) / a.shape[1]
        feats.append(cnt)
    return np.hstack(feats)

# ---------------------------------------------------------------- training ----
def load_mat(p):
    path = os.path.join(DATA, f"{p}8.mat")
    if not os.path.exists(path) or os.path.getsize(path) < 10000:
        say(f"downloading {p}8.mat"); urllib.request.urlretrieve(BASE.format(p), path)
    import scipy.io as sio
    m = sio.loadmat(path)
    c = lambda k: np.array([str(x).strip() for x in np.ravel(m[k])])
    return c("Seq_pos"), c("ID_pos"), c("Seq_neg"), c("ID_neg")

def train(p):
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.metrics import roc_auc_score
    import joblib
    mp = os.path.join(DATA, f"rf_{p}_{N_TRAIN}_{N_TREES}_{SEED}.joblib")
    if os.path.exists(mp):
        d = joblib.load(mp); say(f"{p}: model loaded from cache"); return d
    sp, ip, sn, in_ = load_mat(p)
    ok = lambda s: np.array([len(x) == 17 and x[8] == "K" for x in s])
    kp, kn = ok(sp), ok(sn); sp, ip, sn, in_ = sp[kp], ip[kp], sn[kn], in_[kn]
    rng = np.random.default_rng(SEED)
    ids = np.unique(np.concatenate([ip, in_])); rng.shuffle(ids); test_ids = set(ids[: max(1, len(ids) // 5)])           # split by protein
    def pick(s, i, test, n):
        m = np.array([(x in test_ids) == test for x in i]); j = np.flatnonzero(m)
        return s[rng.choice(j, min(n, len(j)), replace=False)]
    trp, trn = pick(sp, ip, False, N_TRAIN), pick(sn, in_, False, N_TRAIN)
    tep, ten = pick(sp, ip, True, 4000), pick(sn, in_, True, 4000)
    n = min(len(trp), len(trn)); trp, trn = trp[:n], trn[:n]
    X = encode(list(trp) + list(trn)); y = np.r_[np.ones(n), np.zeros(n)]
    say(f"{p} ({PTM_NAMES[p]}): training on {n} + {n} windows ({X.shape[1]} features)")
    rf = RandomForestClassifier(n_estimators=N_TREES, min_samples_leaf=5, n_jobs=-1, random_state=SEED).fit(X, y)
    sc_p, sc_n = rf.predict_proba(encode(list(tep)))[:, 1], rf.predict_proba(encode(list(ten)))[:, 1]
    auc = roc_auc_score(np.r_[np.ones(len(sc_p)), np.zeros(len(sc_n))], np.r_[sc_p, sc_n])
    thr = float(np.quantile(sc_n, 1 - FPR)); sens = float(np.mean(sc_p >= thr))
    say(f"{p}: held-out proteins AUC {auc:.3f}; threshold {thr:.3f} (FPR {FPR:.0%}) gives sensitivity {sens:.1%}")
    d = {"rf": rf, "thr": thr, "auc": auc, "sens": sens, "n_train": n, "n_pos_total": int(kp.sum()), "n_neg_total": int(kn.sum())}
    joblib.dump(d, mp, compress=3); return d

# ---------------------------------------------------------------- alignment ----
B62 = """A R N D C Q E G H I L K M F P S T W Y V
A 4 -1 -2 -2 0 -1 -1 0 -2 -1 -1 -1 -1 -2 -1 1 0 -3 -2 0
R -1 5 0 -2 -3 1 0 -2 0 -3 -2 2 -1 -3 -2 -1 -1 -3 -2 -3
N -2 0 6 1 -3 0 0 0 1 -3 -3 0 -2 -3 -2 1 0 -4 -2 -3
D -2 -2 1 6 -3 0 2 -1 -1 -3 -4 -1 -3 -3 -1 0 -1 -4 -3 -3
C 0 -3 -3 -3 9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
Q -1 1 0 0 -3 5 2 -2 0 -3 -2 1 0 -3 -1 0 -1 -2 -1 -2
E -1 0 0 2 -4 2 5 -2 0 -3 -3 1 -2 -3 -1 0 -1 -3 -2 -2
G 0 -2 0 -1 -3 -2 -2 6 -2 -4 -4 -2 -3 -3 -2 0 -2 -2 -3 -3
H -2 0 1 -1 -3 0 0 -2 8 -3 -3 -1 -2 -1 -2 -1 -2 -2 2 -3
I -1 -3 -3 -3 -1 -3 -3 -4 -3 4 2 -3 1 0 -3 -2 -1 -3 -1 3
L -1 -2 -3 -4 -1 -2 -3 -4 -3 2 4 -2 2 0 -3 -2 -1 -2 -1 1
K -1 2 0 -1 -3 1 1 -2 -1 -3 -2 5 -1 -3 -1 0 -1 -3 -2 -2
M -1 -1 -2 -3 -1 0 -2 -3 -2 1 2 -1 5 0 -2 -1 -1 -1 -1 1
F -2 -3 -3 -3 -2 -3 -3 -3 -1 0 0 -3 0 6 -4 -2 -2 1 3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4 7 -1 -1 -4 -3 -2
S 1 -1 1 0 -1 0 0 0 -1 -2 -2 0 -1 -2 -1 4 1 -3 -2 -2
T 0 -1 0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1 1 5 -2 -2 0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1 1 -4 -3 -2 11 2 -3
Y -2 -2 -2 -3 -2 -1 -2 -3 2 -1 -1 -2 -1 3 -3 -2 -2 2 7 -1
V 0 -3 -3 -3 -1 -2 -2 -3 -3 3 1 -2 1 -1 -2 -2 0 -3 -1 4"""
_rows = [r.split() for r in B62.split("\n")]; _cols = _rows[0]
SUB = {(r[0], _cols[j]): int(v) for r in _rows[1:] for j, v in enumerate(r[1:])}

def align(a, b, gap=-6):                       # global (Needleman-Wunsch) alignment, linear gap; returns aligned strings
    n, m = len(a), len(b)
    S = np.zeros((n + 1, m + 1), dtype=np.int32); T = np.zeros((n + 1, m + 1), dtype=np.int8)
    S[:, 0] = np.arange(n + 1) * gap; S[0, :] = np.arange(m + 1) * gap; T[1:, 0] = 1; T[0, 1:] = 2
    sc = np.array([[SUB.get((x, y), -1) for y in b] for x in a], dtype=np.int32)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            d, u, l = S[i - 1, j - 1] + sc[i - 1, j - 1], S[i - 1, j] + gap, S[i, j - 1] + gap
            if d >= u and d >= l: S[i, j], T[i, j] = d, 0
            elif u >= l: S[i, j], T[i, j] = u, 1
            else: S[i, j], T[i, j] = l, 2
    i, j, ra, rb = n, m, [], []
    while i > 0 or j > 0:
        t = T[i, j]
        if i > 0 and j > 0 and t == 0: ra.append(a[i - 1]); rb.append(b[j - 1]); i -= 1; j -= 1
        elif i > 0 and (t == 1 or j == 0): ra.append(a[i - 1]); rb.append("-"); i -= 1
        else: ra.append("-"); rb.append(b[j - 1]); j -= 1
    return "".join(ra[::-1]), "".join(rb[::-1])

# ---------------------------------------------------------------- main ----
def main():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    files = [f for f in sys.argv[1:] if not f.startswith("-")]
    if not files: sys.exit("Usage: python 99_lysine_ptm_predict_compare.py eGP.FASTA [mGP.FASTA]")
    prots = []
    for f in files:
        recs = read_fasta(f)
        if not recs: sys.exit(f"No sequence in {f}")
        prots.append((os.path.splitext(os.path.basename(f))[0], recs[0][0].split()[0], recs[0][1]))
        say(f"{f}: {len(recs[0][1])} aa, {recs[0][1].count('K')} lysines" + (f" (only the first of {len(recs)} records is used)" if len(recs) > 1 else ""))
    models = {p: train(p) for p in PTMS}
    with open(f"{OUT}/training_report.tsv", "w") as t:
        t.write("modification\tname\tpositives_available\tnegatives_available\ttrained_on_each_class\theldout_protein_AUC\tcall_threshold\tsensitivity_at_threshold\tFPR\n")
        for p, d in models.items(): t.write(f"{p}\t{PTM_NAMES[p]}\t{d['n_pos_total']}\t{d['n_neg_total']}\t{d['n_train']}\t{d['auc']:.3f}\t{d['thr']:.3f}\t{d['sens']:.3f}\t{FPR}\n")
    P = {}                                                                                   # protein -> {pos: {ptm: (prob, call)}}
    for name, rec, seq in prots:
        wl = windows(seq); X = encode([w for _, w in wl]); P[name] = {}
        cols = {}
        for p, d in models.items(): cols[p] = d["rf"].predict_proba(X)[:, 1]
        with open(f"{OUT}/{name}_ptm_predictions.tsv", "w") as t:
            t.write("record\tposition\twindow\t" + "\t".join(f"P_{p}" for p in PTMS) + "\t" + "\t".join(f"call_{p}" for p in PTMS) + "\tn_calls\n")
            for k, (pos, w) in enumerate(wl):
                pr = {p: float(cols[p][k]) for p in PTMS}; ca = {p: pr[p] >= models[p]["thr"] for p in PTMS}; P[name][pos] = {p: (pr[p], ca[p]) for p in PTMS}
                t.write(f"{rec}\t{pos}\t{w}\t" + "\t".join(f"{pr[p]:.3f}" for p in PTMS) + "\t" + "\t".join("1" if ca[p] else "0" for p in PTMS) + f"\t{sum(ca.values())}\n")
        say(f"{name}: calls per modification (lysines predicted at FPR {FPR:.0%}): " + ", ".join(f"{p} {sum(P[name][q][p][1] for q in P[name])}/{len(wl)}" for p in PTMS))
    # ---- per-protein counts figure (and comparison, if two proteins) ----
    names = [p[0] for p in prots]; cmap = ["#1B7837", "#2B6CB0", "#E08214", "#762A83"]
    fig, ax = plt.subplots(figsize=(1.1 * len(PTMS) + 2.5, 4)); wd = 0.8 / len(names)
    for i, nm in enumerate(names):
        v = [sum(P[nm][q][p][1] for q in P[nm]) for p in PTMS]
        ax.bar(np.arange(len(PTMS)) + i * wd, v, wd, label=f"{nm} ({len(P[nm])} K)", color=cmap[i % 4])
        for x, y in zip(np.arange(len(PTMS)) + i * wd, v): ax.text(x, y + 0.1, str(y), ha="center", fontsize=7)
    ax.set_xticks(np.arange(len(PTMS)) + wd * (len(names) - 1) / 2); ax.set_xticklabels([PTM_NAMES[p] for p in PTMS], rotation=30, ha="right", fontsize=8)
    ax.set_ylabel(f"lysines predicted (FPR {FPR:.0%})"); ax.legend(frameon=False, fontsize=8); ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("Predicted modified lysines per modification", loc="left", fontweight="bold", fontsize=10)
    for e in ("png", "pdf"): fig.savefig(f"{OUT}/sites_per_modification.{e}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    if len(prots) < 2:
        say("one protein: predictions written; give a second FASTA to get the comparison"); return
    (A, ra, sa), (B, rb, sb) = prots[0], prots[1]
    al, bl = align(sa, sb); ncol = len(al)
    ident = sum(x == y and x != "-" for x, y in zip(al, bl)) / max(1, sum(x != "-" and y != "-" for x, y in zip(al, bl)))
    say(f"alignment: {ncol} columns, {100 * ident:.0f}% identity over aligned residues (a global alignment of divergent proteins is approximate)")
    pa = pb = 0; rows = []
    for c, (x, y) in enumerate(zip(al, bl)):
        if x != "-": pa += 1
        if y != "-": pb += 1
        ka, kb = x == "K", y == "K"
        if ka or kb:
            rows.append(dict(col=c + 1, pos_A=pa if x != "-" else 0, pos_B=pb if y != "-" else 0, res_A=x, res_B=y, status="conserved K" if ka and kb else (f"K only in {A}" if ka else f"K only in {B}")))
    with open(f"{OUT}/comparison_aligned_lysines.tsv", "w") as t:
        t.write(f"alignment_column\tpos_{A}\tpos_{B}\tresidue_{A}\tresidue_{B}\tstatus\t" + "\t".join(f"P_{A}_{p}\tP_{B}_{p}\tdiff_{p}\tcategory_{p}" for p in PTMS) + "\n")
        for r in rows:
            line = [str(r["col"]), str(r["pos_A"]), str(r["pos_B"]), r["res_A"], r["res_B"], r["status"]]
            for p in PTMS:
                a = P[A].get(r["pos_A"], {}).get(p); b = P[B].get(r["pos_B"], {}).get(p)
                ca, cb = bool(a and a[1]), bool(b and b[1])
                r.setdefault("v", {})[p] = (a[0] if a else np.nan, b[0] if b else np.nan, ca, cb)
                cat = ("both" if ca and cb else f"{A} only" if ca else f"{B} only" if cb else "neither") if (a and b) else (f"{A} lysine, predicted" if ca else f"{A} lysine, not predicted") if a else (f"{B} lysine, predicted" if cb else f"{B} lysine, not predicted")
                r.setdefault("cat", {})[p] = cat
                line += [f"{a[0]:.3f}" if a else "", f"{b[0]:.3f}" if b else "", f"{(a[0] - b[0]):.3f}" if a and b else "", cat]
            t.write("\t".join(line) + "\n")
    with open(f"{OUT}/comparison_summary.tsv", "w") as t:
        t.write(f"modification\t{A}_lysines_predicted\t{B}_lysines_predicted\tconserved_K_pairs\tpredicted_in_both\t{A}_only\t{B}_only\tneither\tscore_correlation_conserved_K\tconserved_pairs_with_gain_or_loss\n")
        for p in PTMS:
            cons = [r for r in rows if r["status"] == "conserved K"]; v = [r["v"][p] for r in cons]
            both = sum(1 for a, b, ca, cb in v if ca and cb); ao = sum(1 for a, b, ca, cb in v if ca and not cb); bo = sum(1 for a, b, ca, cb in v if cb and not ca); ne = len(v) - both - ao - bo
            xs, ys = np.array([a for a, b, _, _ in v]), np.array([b for a, b, _, _ in v]); cor = float(np.corrcoef(xs, ys)[0, 1]) if len(v) > 2 else float("nan")
            t.write(f"{p} ({PTM_NAMES[p]})\t{sum(P[A][q][p][1] for q in P[A])}\t{sum(P[B][q][p][1] for q in P[B])}\t{len(cons)}\t{both}\t{ao}\t{bo}\t{ne}\t{cor:.2f}\t{ao + bo}\n")
    say("comparison tables written; drawing figures")
    # ---- figure: per-modification maps of both proteins ----
    cm = plt.get_cmap("viridis"); fig, axes = plt.subplots(2, 1, figsize=(13, 0.55 * len(PTMS) * 2 + 2), sharex=False)
    for ax, (nm, seq) in zip(axes, [(A, sa), (B, sb)]):
        for r, p in enumerate(PTMS):
            for pos, d in P[nm].items():
                pr, ca = d[p]; ax.scatter(pos, r, s=26 if ca else 10, c=[cm(min(1, pr))], edgecolors="black" if ca else "none", linewidths=0.8, zorder=3 if ca else 2)
        ax.set_yticks(range(len(PTMS))); ax.set_yticklabels([PTM_NAMES[p] for p in PTMS], fontsize=8); ax.set_xlim(0, len(seq) + 1); ax.set_ylim(-0.7, len(PTMS) - 0.3); ax.invert_yaxis()
        ax.set_title(f"{nm}: every lysine, colour = score, black ring = called (FPR {FPR:.0%})", loc="left", fontsize=9, fontweight="bold"); ax.spines[["top", "right"]].set_visible(False)
    axes[1].set_xlabel("residue position (each protein in its own numbering)")
    sm = plt.cm.ScalarMappable(cmap=cm); fig.colorbar(sm, ax=axes, fraction=0.015, pad=0.01, label="RF score")
    for e in ("png", "pdf"): fig.savefig(f"{OUT}/ptm_map_{A}_vs_{B}.{e}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    # ---- figure: conserved lysine scatter ----
    cons = [r for r in rows if r["status"] == "conserved K"]
    nc = min(4, len(PTMS)); nr = int(np.ceil(len(PTMS) / nc)); fig, axs = plt.subplots(nr, nc, figsize=(3.2 * nc, 3.2 * nr), squeeze=False)
    for ax, p in zip(axs.ravel(), PTMS):
        xs = np.array([r["v"][p][0] for r in cons]); ys = np.array([r["v"][p][1] for r in cons]); thr = models[p]["thr"]
        col = ["#C53030" if (r["v"][p][2] != r["v"][p][3]) else ("#1B7837" if r["v"][p][2] else "#9AA0A6") for r in cons]
        ax.scatter(xs, ys, c=col, s=14); ax.axvline(thr, color="grey", lw=0.6, ls="--"); ax.axhline(thr, color="grey", lw=0.6, ls="--"); ax.plot([0, 1], [0, 1], color="grey", lw=0.5)
        for r, x, y in zip(cons, xs, ys):
            if r["v"][p][2] != r["v"][p][3] and abs(x - y) > 0.15: ax.text(x, y, f"K{r['pos_A']}/K{r['pos_B']}", fontsize=5.5)
        ax.set_title(PTM_NAMES[p], fontsize=9); ax.set_xlabel(f"{A} score", fontsize=8); ax.set_ylabel(f"{B} score", fontsize=8); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.spines[["top", "right"]].set_visible(False)
    for ax in axs.ravel()[len(PTMS):]: ax.axis("off")
    fig.suptitle(f"Conserved lysines ({len(cons)} aligned K pairs): score in {A} vs {B}  (green = called in both, red = called in one, grey = neither)", fontsize=9, x=0.01, ha="left")
    for e in ("png", "pdf"): fig.savefig(f"{OUT}/aligned_lysine_scatter.{e}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    # ---- figure: difference map along the alignment ----
    fig, ax = plt.subplots(figsize=(13, 0.5 * len(PTMS) + 2)); M = np.full((len(PTMS), ncol), np.nan)
    for r in rows:
        for k, p in enumerate(PTMS):
            a, b, _, _ = r["v"][p]
            M[k, r["col"] - 1] = (a if not np.isnan(a) else 0) - (b if not np.isnan(b) else 0)
    im = ax.imshow(np.ma.masked_invalid(M), aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1, interpolation="nearest")
    ax.set_yticks(range(len(PTMS))); ax.set_yticklabels([PTM_NAMES[p] for p in PTMS], fontsize=8); ax.set_xlabel("alignment column")
    ax.set_title(f"Score difference at lysines ({A} minus {B}); a lysine absent from one protein counts as 0 there", loc="left", fontsize=9, fontweight="bold")
    fig.colorbar(im, ax=ax, fraction=0.015, pad=0.01, label=f"{A} score - {B} score")
    for e in ("png", "pdf"): fig.savefig(f"{OUT}/ptm_difference_map.{e}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    say(f"done. Everything is in ./{OUT}/")

if __name__ == "__main__":
    main()
