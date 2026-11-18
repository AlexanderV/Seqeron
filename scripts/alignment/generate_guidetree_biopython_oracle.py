#!/usr/bin/env python3
"""ALIGN-MULTI-001 (audit A1-12) oracle: UPGMA (scipy average linkage) and neighbour-joining (Biopython 1.88) trees on random distance matrices.

Seqeron API: PhylogeneticAnalyzer.BuildTreeFromMatrix(taxa, matrix, TreeMethod.UPGMA | NeighborJoining) — the canonical builders the
guide tree of MultipleAlignProgressive / Iterative / Consistency now calls (SequenceAligner.BuildProgressiveGuideTree).

Reference (run here):
  * NJ: Biopython 1.88 Bio.Phylo.TreeConstruction.DistanceTreeConstructor().nj(...) (unrooted, trifurcating centre), full-precision lengths.
  * UPGMA: scipy.cluster.hierarchy.linkage(method='average') — the size-weighted (proportional) average of Sokal & Michener 1958 =
    UPGMA — converted to a rooted ultrametric tree (branch length = (height(parent) - height(child)) / 2, leaves height 0).
    NOT Biopython's DistanceTreeConstructor.upgma(): its update line is `dm[min_j, k] = (dm[min_i, k] + dm[min_j, k]) / 2`, the plain
    mean without cluster sizes, i.e. WPGMA; it agrees with UPGMA only until the first cluster of size >= 2 is averaged again.
Matrices are random real symmetric (6 decimals), so no ties occur and the topology is unique.

Output TSV (tests/.../TestData/Alignment/biopython_guidetree_oracle.tsv):
  n  matrix(row-major, ';' between rows, ',' between entries)  upgmaNewick  njNewick
Usage: python3 scripts/alignment/generate_guidetree_biopython_oracle.py OUTFILE
"""
import io, random, sys

import numpy as np
from Bio import Phylo
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor


def newick(tree):
    s = io.StringIO()
    Phylo.write(tree, s, "newick", branch_length_only=False, plain=False)
    return s.getvalue().strip()


def scipy_upgma_newick(m):
    n = len(m)
    root = to_tree(linkage(squareform(np.array(m), checks=False), method="average"))

    def fmt(node, parent_height):
        h = 0.0 if node.is_leaf() else node.dist / 2
        bl = parent_height - h
        label = f"t{node.id}" if node.is_leaf() else ""
        if node.is_leaf():
            return f"{label}:{float(bl)!r}"
        return f"({fmt(node.get_left(), h)},{fmt(node.get_right(), h)}):{float(bl)!r}"

    h = root.dist / 2
    return f"({fmt(root.get_left(), h)},{fmt(root.get_right(), h)});"


def main(path):
    rng = random.Random(20261118)
    rows = []
    for t in range(400):
        n = rng.randint(3, 14)
        m = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                m[i][j] = m[j][i] = round(rng.uniform(0.05, 2.0), 6)
        names = [f"t{i}" for i in range(n)]
        dm = DistanceMatrix(names, [m[i][: i + 1] for i in range(n)])
        c = DistanceTreeConstructor()
        up = scipy_upgma_newick(m)
        nj = newick(c.nj(dm))
        rows.append((str(n), ";".join(",".join(repr(x) for x in r) for r in m), up, nj))
    with open(path, "w") as f:
        f.write("n\tmatrix\tupgma\tnj\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    print(len(rows), "rows ->", path)


if __name__ == "__main__":
    main(sys.argv[1])
