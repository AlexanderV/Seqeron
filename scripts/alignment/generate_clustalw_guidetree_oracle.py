#!/usr/bin/env python3
"""ALIGN-MULTI-001 (audit A1-12) oracle: the neighbour-joining guide tree of ClustalW 2.1 and the order of its progressive steps.

Seqeron API: SequenceAligner.BuildClustalWNeighborJoiningGuideTree (internal; canonical PhylogeneticAnalyzer.BuildTreeFromMatrix
NeighborJoining + ClustalW's printing / Tree::reRoot / Tree::createSets steps) and MultipleAlign*(..., MsaOptions { GuideTree =
MsaGuideTree.NeighborJoining }).

Reference (run here): the installed `clustalw` 2.1 (Debian clustalw 2.1+lgpl-7) `-TYPE=DNA -ALIGN -MAXDIV=0 -NEWTREE=x.dnd` (MAXDIV=0: no delayed divergent sequences, so every step is printed):
  * x.dnd = the unrooted neighbour-joining guide tree it writes (branch lengths with 5 decimals), built from the pairwise
    percent identities of its "slow" pairwise stage (F51);
  * the "Group n: Sequences: s" lines of the progressive stage = the sizes of the merges of the ROOTED tree in post-order
    (Tree::reRoot + Tree::groupSeqs of src/tree/Tree.cpp) — the order in which profiles are joined.
Both are recorded verbatim; nothing here models the binary.

Output TSV (tests/.../TestData/Alignment/clustalw_guidetree_oracle.tsv):
  seqs(;)  dnd(one line)  stepSizes(,)
Usage: python3 scripts/alignment/generate_clustalw_guidetree_oracle.py OUTFILE      (needs `clustalw` on PATH)
"""
import os, random, re, subprocess, sys, tempfile


def mutate(rng, s, p):
    out = []
    for c in s:
        r = rng.random()
        if r < p:
            out.append(rng.choice("ACGT"))
        elif r < p * 1.25:
            continue
        elif r < p * 1.5:
            out.append(c + rng.choice("ACGT"))
        else:
            out.append(c)
    return "".join(out) or "A"


def make_set(rng, kind):
    n = rng.randint(3, 11)
    if kind == "family":
        base = "".join(rng.choice("ACGT") for _ in range(rng.randint(30, 90)))
        return [mutate(rng, base, rng.choice([0.05, 0.15, 0.3])) for _ in range(n)]
    if kind == "twoclade":
        base = "".join(rng.choice("ACGT") for _ in range(rng.randint(40, 80)))
        a, b = mutate(rng, base, 0.3), mutate(rng, base, 0.3)
        return [mutate(rng, a if i % 2 == 0 else b, 0.05) for i in range(n)]
    if kind == "duplicates":
        base = "".join(rng.choice("ACGT") for _ in range(rng.randint(20, 60)))
        pool = [mutate(rng, base, 0.2) for _ in range(max(2, n // 2))]
        return [rng.choice(pool) for _ in range(n)]
    if kind == "twoletter":
        return ["".join(rng.choice("AC") for _ in range(rng.randint(12, 30))) for _ in range(n)]
    if kind == "random":
        return ["".join(rng.choice("ACGT") for _ in range(rng.randint(15, 70))) for _ in range(n)]
    if kind == "tiny":
        return ["".join(rng.choice("ACGT") for _ in range(rng.randint(3, 9))) for _ in range(rng.randint(2, 7))]
    raise ValueError(kind)


def run_clustalw(seqs):
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, "a.fa")
        with open(fa, "w") as f:
            for i, s in enumerate(seqs):
                f.write(f">s{i}\n{s}\n")
        out = subprocess.run(
            ["clustalw", f"-INFILE={fa}", "-TYPE=DNA", "-ALIGN", "-MAXDIV=0", f"-OUTFILE={d}/a.aln", f"-NEWTREE={d}/a.dnd"],
            capture_output=True, text=True, check=True).stdout
        dnd = open(os.path.join(d, "a.dnd")).read()
    steps = [int(m.group(1)) for m in re.finditer(r"^Group\s+\d+: Sequences:\s+(\d+)", out, re.M)]
    return re.sub(r"\s+", "", dnd), steps


def main(path):
    rng = random.Random(20261118)
    kinds = ["family"] * 5 + ["twoclade"] * 3 + ["duplicates"] * 4 + ["twoletter"] * 2 + ["random"] * 3 + ["tiny"] * 3
    rows = []
    for t in range(420):
        kind = kinds[t % len(kinds)]
        seqs = make_set(rng, kind)
        if len(seqs) < 2:
            continue
        dnd, steps = run_clustalw(seqs)
        assert len(steps) == len(seqs) - 1, (seqs, steps)
        rows.append((";".join(seqs), dnd, ",".join(map(str, steps))))
    with open(path, "w") as f:
        f.write("seqs\tdnd\tstepSizes\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    print(len(rows), "rows ->", path)


if __name__ == "__main__":
    main(sys.argv[1])
