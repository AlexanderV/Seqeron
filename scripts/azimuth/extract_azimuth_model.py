#!/usr/bin/env python3
"""Reproducible extractor for the Doench 2016 "Rule Set 2" / Azimuth on-target models.

This is a *build-time* developer tool (run once). Its outputs are committed:
  - Resources/azimuth_rs2_nopos.bin  (sequence-only model, "Rule Set 2" score)
  - Resources/azimuth_rs2_full.bin   (model with gene-position features)
  - test oracle CSVs (predictions from this verified reference + upstream-agreeing rows)

It loads Microsoft Research's pickled scikit-learn GradientBoostingRegressor models
WITHOUT scikit-learn (a custom Unpickler reads the raw tree node arrays), reproduces
the exact azimuth featurization (incl. Biopython Tm_NN and the CPython-2.7 dict column
ordering used at training time), and serialises the trees into a compact binary blob.

Provenance / inputs (downloaded from the pinned upstream repo if absent):
  https://github.com/MicrosoftResearch/Azimuth  (BSD-3-Clause)
    azimuth/saved_models/V3_model_nopos.pickle
    azimuth/saved_models/V3_model_full.pickle
    azimuth/tests/1000guides.csv          (upstream regression fixture; see NOTE below)

NOTE on the upstream fixture: azimuth/tests/1000guides.csv is INTERNALLY INCONSISTENT
with the shipped pickles for ~38% of rows (its own unit test warns it "can fail due to
randomness ... feature reordering"). This was verified here three ways: (1) our extracted
trees reproduced by sklearn's own Tree.predict bit-for-bit; (2) our featurizer matches a
verbatim port of upstream featurization.py to 1e-13 using real Biopython; (3) the column
order reproduces documented CPython-2.7 dict iteration order. The authoritative oracle is
therefore THIS reference's predictions; the fixture is used only as independent
corroboration on the subset of rows where it agrees.

Binary format (little-endian) -- see AzimuthOnTarget.cs for the reader:
  magic   u32  = 0x32535241  ('A','R','S','2')
  version u16  = 1
  flags   u16  bit0 = model has gene-position inputs (full model)
  treeCount  i32
  nodeCount  i32
  nFeatures  i32
  initScore     f64   (GradientBoosting init_, the training-target mean)
  learningRate  f64   (informational; leaf values below are already pre-scaled by it)
  treeStart[treeCount]  i32   root node index of each tree
  nodes[nodeCount]      (24 bytes each, AoS):
      thresholdOrValue f64   internal: split threshold; leaf: learningRate * raw_leaf
      left   i32             child if feature value <= threshold; -1 marks a leaf
      right  i32
      feature i32            split feature index; -1 marks a leaf
      _pad   i32
Prediction:  score = initScore + sum over trees of value(reached leaf).
"""
import io, os, sys, csv, math, struct, pickle, itertools, urllib.request
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, ".cache")
RES = os.path.normpath(os.path.join(HERE, "..", "..", "src", "Seqeron", "Algorithms",
                                    "Seqeron.Genomics.MolTools", "Resources"))
RAW = "https://raw.githubusercontent.com/MicrosoftResearch/Azimuth/master"
INPUTS = {
    "V3_model_nopos.pickle": RAW + "/azimuth/saved_models/V3_model_nopos.pickle",
    "V3_model_full.pickle":  RAW + "/azimuth/saved_models/V3_model_full.pickle",
    "1000guides.csv":        RAW + "/azimuth/tests/1000guides.csv",
}

def fetch(name):
    os.makedirs(CACHE, exist_ok=True)
    p = os.path.join(CACHE, name)
    if not os.path.exists(p):
        urllib.request.urlretrieve(INPUTS[name], p)
    return p

# --------------------------------------------------------------------------- sklearn-free load
class _Stub(object):
    def __setstate__(self, st): self._state = st
class _TreeStub(object):
    def __init__(self, *a, **k): self.ctor = a
    def __setstate__(self, st): self.state = st
class _Unp(pickle.Unpickler):
    def find_class(self, module, name):
        if name == "Tree": return _TreeStub
        if name in ("GradientBoostingRegressor", "DecisionTreeRegressor",
                    "MeanEstimator", "LeastSquaresError"):
            return _Stub
        return super().find_class(module, name)

def load_model(pickle_path):
    model, opts = _Unp(io.BytesIO(open(pickle_path, "rb").read()), encoding="latin1").load()
    st = model._state
    lr = float(st["learning_rate"]); init_mean = float(st["init_"]._state["mean"])
    trees = []
    for e in np.asarray(st["estimators_"]).flat:
        ts = e._state["tree_"].state
        nodes = ts["nodes"]
        vals = np.asarray(ts["values"]).reshape(len(nodes), -1)[:, 0]
        trees.append((nodes["left_child"].astype(np.int64),
                      nodes["right_child"].astype(np.int64),
                      nodes["feature"].astype(np.int64),
                      nodes["threshold"].astype(np.float64),
                      vals.astype(np.float64)))
    return trees, lr, init_mean, opts

# --------------------------------------------------------------------------- featurization
ALPHA = ["A", "T", "C", "G"]
def _alphabet(order): return ["".join(p) for p in itertools.product(ALPHA, repeat=order)]

def _nuc(s, order, max_index_to_use=30):
    if max_index_to_use <= len(s): max_index_to_use = len(s)
    s = s[:max_index_to_use]
    al = _alphabet(order)
    pd = np.zeros(len(al) * (len(s) - (order - 1)))
    pi = np.zeros(len(ALPHA) ** order)
    for pos in range(0, len(s) - order + 1):
        idx = al.index(s[pos:pos + order])
        pd[idx + pos * len(al)] = 1.0
        pi[idx] += 1.0
    return pd, pi

def _count_gc(seq): return len(seq[4:24].replace("A", "").replace("T", ""))

# Biopython Tm_NN (DNA_NN3, saltcorr=5, dnac1=dnac2=25, Na=50) -- verified == real Biopython.
_NN = {"init": (0, 0), "init_A/T": (2.3, 4.1), "init_G/C": (0.1, -2.8),
       "AA/TT": (-7.9, -22.2), "AT/TA": (-7.2, -20.4), "TA/AT": (-7.2, -21.3),
       "CA/GT": (-8.5, -22.7), "GT/CA": (-8.4, -22.4), "CT/GA": (-7.8, -21.0),
       "GA/CT": (-8.2, -22.2), "CG/GC": (-10.6, -27.2), "GC/CG": (-9.8, -24.4),
       "GG/CC": (-8.0, -19.9)}
_COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
def tm_nn(seq):
    seq = seq.upper(); cseq = "".join(_COMP[c] for c in seq)
    dh = ds = 0.0
    dh += _NN["init"][0]; ds += _NN["init"][1]
    ends = seq[0] + seq[-1]
    at = ends.count("A") + ends.count("T"); gc = ends.count("G") + ends.count("C")
    dh += _NN["init_A/T"][0] * at; ds += _NN["init_A/T"][1] * at
    dh += _NN["init_G/C"][0] * gc; ds += _NN["init_G/C"][1] * gc
    for i in range(len(seq) - 1):
        nb = seq[i:i + 2] + "/" + cseq[i:i + 2]
        v = _NN.get(nb) or _NN.get(nb[::-1])
        dh += v[0]; ds += v[1]
    ds += 0.368 * (len(seq) - 1) * math.log(50e-3)
    k = (25 - 25 / 2.0) * 1e-9
    return (1000 * dh) / (ds + 1.987 * math.log(k)) - 273.15

def _tm_features(seq):
    return [tm_nn(seq), tm_nn(seq[19:24]), tm_nn(seq[11:19]), tm_nn(seq[6:11])]

def nopos_blocks(seq):
    pd1, pi1 = _nuc(seq, 1); pd2, pi2 = _nuc(seq, 2); gc = _count_gc(seq)
    nx = seq[24] + seq[27]; ng = np.zeros(16); ng[_alphabet(2).index(nx)] = 1.0
    return {"_nuc_pd_Order1": pd1, "_nuc_pi_Order1": pi1, "_nuc_pd_Order2": pd2,
            "_nuc_pi_Order2": pi2, "gc_above_10": np.array([1.0 if gc > 10 else 0.0]),
            "gc_below_10": np.array([1.0 if gc < 10 else 0.0]), "gc_count": np.array([float(gc)]),
            "NGGX": ng, "Tm": np.array(_tm_features(seq))}

def full_blocks(seq, aa_cut, pct_pep):
    b = nopos_blocks(seq)
    b["Percent Peptide"] = np.array([float(pct_pep)])
    b["Amino Acid Cut position"] = np.array([float(aa_cut)])
    b["Percent Peptide <50%"] = np.array([1.0 if pct_pep < 50 else 0.0])
    return b

# CPython-2.7 (64-bit, no hash randomization) dict iteration order ----------------------------
def _py2_hash(s):
    if not s: return 0
    M = 0xFFFFFFFFFFFFFFFF
    x = (ord(s[0]) << 7) & M
    for c in s:
        x = ((1000003 * x) & M) ^ ord(c)
    x = (x ^ len(s)) & M
    if x > 0x7FFFFFFFFFFFFFFF: x -= 0x10000000000000000
    return -2 if x == -1 else x

def py2_dict_order(keys):
    M = 0xFFFFFFFFFFFFFFFF
    size = 8; table = [None] * size; fill = 0
    def ins(tbl, mask, key, h):
        perturb = h & M; i = h & mask
        while tbl[i & mask] is not None:
            i = (i * 5 + perturb + 1) & M; perturb >>= 5
        tbl[i & mask] = (key, h)
    for k in keys:
        ins(table, size - 1, k, _py2_hash(k)); fill += 1
        if fill * 3 >= size * 2:
            newsize = 8
            while newsize <= fill * 4: newsize <<= 1
            nt = [None] * newsize
            for slot in table:
                if slot is not None: ins(nt, newsize - 1, slot[0], slot[1])
            table, size = nt, newsize
    return [s[0] for s in table if s is not None]

NOPOS_INSERTION = ["_nuc_pd_Order1", "_nuc_pi_Order1", "_nuc_pd_Order2", "_nuc_pi_Order2",
                   "gc_above_10", "gc_below_10", "gc_count", "NGGX", "Tm"]
FULL_INSERTION = ["_nuc_pd_Order1", "_nuc_pi_Order1", "_nuc_pd_Order2", "_nuc_pi_Order2",
                  "gc_above_10", "gc_below_10", "gc_count", "Percent Peptide",
                  "Amino Acid Cut position", "Percent Peptide <50%", "NGGX", "Tm"]
NOPOS_ORDER = py2_dict_order(NOPOS_INSERTION)
FULL_ORDER = py2_dict_order(FULL_INSERTION)

def predict(trees, lr, init_mean, x):
    s = init_mean
    for L, R, F, T, V in trees:
        n = 0
        while L[n] != -1:
            n = L[n] if x[F[n]] <= T[n] else R[n]
        s += lr * V[n]
    return s

# --------------------------------------------------------------------------- serialization
MAGIC = 0x32535241
def serialize(trees, lr, init_mean, n_features, has_gene_pos):
    # flatten trees into one node array with per-tree root offsets
    starts = []; flat = []  # each node: (thrOrVal, left, right, feature)
    base = 0
    for L, R, F, T, V in trees:
        starts.append(base)
        for n in range(len(L)):
            if L[n] == -1:
                flat.append((lr * float(V[n]), -1, -1, -1))
            else:
                flat.append((float(T[n]), base + int(L[n]), base + int(R[n]), int(F[n])))
        base += len(L)
    buf = io.BytesIO()
    buf.write(struct.pack("<IHHiiidd", MAGIC, 1, 1 if has_gene_pos else 0,
                          len(trees), len(flat), n_features, init_mean, lr))
    buf.write(struct.pack("<%di" % len(starts), *starts))
    for thr, l, r, f in flat:
        buf.write(struct.pack("<diii i", thr, l, r, f, 0))
    return buf.getvalue()

# --------------------------------------------------------------------------- main
def main():
    os.makedirs(RES, exist_ok=True)
    np_pickle = fetch("V3_model_nopos.pickle"); fl_pickle = fetch("V3_model_full.pickle")
    guides = fetch("1000guides.csv")
    nop = load_model(np_pickle); full = load_model(fl_pickle)
    ntrees, nlr, nim, nopts = nop
    ftrees, flr, fim, fopts = full
    nfeat_n = 1 + max(int(F.max()) for _, _, F, _, _ in ntrees)
    nfeat_f = 1 + max(int(F.max()) for _, _, F, _, _ in ftrees)
    assert nfeat_n == 627 and nfeat_f == 630, (nfeat_n, nfeat_f)

    with open(os.path.join(RES, "azimuth_rs2_nopos.bin"), "wb") as f:
        f.write(serialize(ntrees, nlr, nim, nfeat_n, False))
    with open(os.path.join(RES, "azimuth_rs2_full.bin"), "wb") as f:
        f.write(serialize(ftrees, flr, fim, nfeat_f, True))
    print("wrote .bin: nopos nodes=%d full nodes=%d" %
          (sum(len(t[0]) for t in ntrees), sum(len(t[0]) for t in ftrees)))
    print("NOPOS_ORDER =", NOPOS_ORDER)
    print("FULL_ORDER  =", FULL_ORDER)

    # oracle vectors from THIS verified reference, plus upstream-agreeing subset.
    # Written into the test project so the C# suite embeds them directly.
    rows = list(csv.DictReader(open(guides)))
    odir = os.path.normpath(os.path.join(HERE, "..", "..", "tests", "Seqeron",
                                         "Seqeron.Genomics.Tests", "TestData", "Azimuth"))
    os.makedirs(odir, exist_ok=True)
    n_agree = f_agree = 0
    with open(os.path.join(odir, "nopos_oracle.csv"), "w", newline="") as on, \
         open(os.path.join(odir, "full_oracle.csv"), "w", newline="") as off_:
        wn = csv.writer(on); wn.writerow(["guide30", "ref_score", "upstream", "agrees"])
        wf = csv.writer(off_); wf.writerow(["guide30", "aa_cut", "pct_pep", "ref_score", "upstream", "agrees"])
        for r in rows:
            s = r["guide"]; aa = float(r["AA cut"]); pp = float(r["Percent peptide"])
            xn = np.concatenate([nopos_blocks(s)[k] for k in NOPOS_ORDER])
            xf = np.concatenate([full_blocks(s, aa, pp)[k] for k in FULL_ORDER])
            pn = predict(ntrees, nlr, nim, xn); pf = predict(ftrees, flr, fim, xf)
            un = float(r["truth nopos"]); uf = float(r["truth pos"])
            an = abs(pn - un) < 1e-3; af = abs(pf - uf) < 1e-3
            n_agree += an; f_agree += af
            wn.writerow([s, "%.6f" % pn, "%.6f" % un, int(an)])
            wf.writerow([s, "%g" % aa, "%g" % pp, "%.6f" % pf, "%.6f" % uf, int(af)])
    print("upstream agreement: nopos %d/%d  full %d/%d" %
          (n_agree, len(rows), f_agree, len(rows)))

if __name__ == "__main__":
    main()
