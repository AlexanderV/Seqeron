"""Oracle for audit round 1, WP5 (A2-1 coverage cut-off, A2-2 tip clipping): literal transcription of Velvet 1.2.10
(dzerbino/velvet, fetched 2026-11-25 with a neutral User-Agent):

  src/concatenatedPreGraph.c  isEligibleTip / clipTips_pg      (multiplicity tip clipping, cutoff 2 * wordlength)
  src/correctedGraph.c        getTipLength / clipTipsHard      (hard tip clipping, cutoff 2 * wordlength)
  src/graphStats.c            removeLowCoverageNodesAndDenounceDubiousReads (two passes; the quotient
                              getTotalCoverage(node) / getNodeLength(node) is an INTEGER division of two
                              Coordinate = int64_t values compared with the double cut-off),
                              estimated_cov (length-weighted median node coverage; cut-off auto = /2)
  src/run2.c                  order: clipTips_pg -> correctGraph (tour bus [not modelled], clipTipsHard) ->
                              estimated_cov -> cov cut-off -> clipTipsHard
  src/preGraph.c              createPreArc_pg: arc multiplicity = number of (k+1)-mer occurrences of the read set,
                              a node and its twin share one arc (so both strands add to the same multiplicity)
  src/concatenatedGraph.c / graphReConstruction.c: node coverage = number of k-mer occurrences in the reads,
                              a node and its twin share it.

The node model is Velvet's bidirected one: node pair i has the two oriented ids +i and -i (the twin); arcs(id)
are the oriented successors, arcs(-id) those of the twin.  Nodes are the unitigs obtained with the node-centric
chain merge of wp4_oracle (not the vertex-path construction used by the C# code); the pair order (Velvet's node
index) is the ordinal order of min(sequence, reverse complement).  Single-strand mode (no twins) uses the reverse
direction view as -id.  Deviation from Velvet that is shared with the C# unit: arcs are the (k-1)-overlaps between
unitigs (not only read-consecutive k-mers), their multiplicity is the (k+1)-mer count of the reads (0 if none).

Usage: python3 -I wp5_oracle.py OUT.json [cases] [seed]   -> case file for the C# harness (see Evidence).
"""
import json
import random
import sys
from collections import defaultdict

sys.path.insert(0, __file__.rsplit('/', 1)[0])
import wp4_oracle as w4  # noqa: E402

revc = w4.revc


def occurrences(reads, k, n):
    c = defaultdict(int)
    for r in reads:
        r = w4.velvet_fold(r)
        for i in range(len(r) - n + 1):
            c[r[i:i + n]] += 1
    return c


def contigs_velvet_hairpin(ks, k, rc):
    """wp4_oracle.contigs_velvet (node-centric chain merge) with Velvet's hairpin rule: concatenatePreGraph_pg does
    not merge across a hairpin arc (isLoop_pg: an arc between a node and its own twin), i.e. x -> revc(x)."""
    succ = {x: sorted(y for y in ks if y[:-1] == x[1:]) for x in ks}
    pred = defaultdict(list)
    for x, ys in succ.items():
        for y in ys:
            pred[y].append(x)
    mergeable = lambda a, b: (len(succ[a]) == 1 and succ[a][0] == b and len(pred[b]) == 1
                              and not (rc and b == revc(a)))
    seen, out = set(), []
    for x in sorted(ks):
        if x in seen:
            continue
        s, guard = x, {x}
        while len(pred[s]) == 1 and mergeable(pred[s][0], s) and pred[s][0] not in guard:
            s = pred[s][0]
            guard.add(s)
        chain, c = [s], s
        seen.add(s)
        while len(succ[c]) == 1 and mergeable(c, succ[c][0]) and succ[c][0] not in seen:
            c = succ[c][0]
            chain.append(c)
            seen.add(c)
        out.append(chain[0] + ''.join(y[-1] for y in chain[1:]))
    return out


class Graph:
    def __init__(self, kset, k, rc, F, F1):
        self.k, self.rc, self.F, self.F1 = k, rc, F, F1
        seqs = contigs_velvet_hairpin(set(kset), k, rc)
        if rc:
            reps = sorted({min(s, revc(s)) for s in seqs})
        else:
            reps = sorted(set(seqs))
        self.seq = {}          # signed id -> sequence of that orientation
        self.alive = set()
        for i, s in enumerate(reps, 1):
            self.seq[i] = s
            if rc:
                self.seq[-i] = revc(s)
            self.alive.add(i)
        self.cov = {}
        for i in self.alive:
            s = self.seq[i]
            tot = 0
            for j in range(len(s) - k + 1):
                km = s[j:j + k]
                tot += F.get(km, 0)
                if rc and revc(km) != km:
                    tot += F.get(revc(km), 0)
            self.cov[i] = tot

    # --- structure
    def length(self, i):
        return len(self.seq[abs(i)]) - self.k + 1

    def _live(self, x):
        return abs(x) in self.alive

    def arcs(self, x):
        """Oriented successors of oriented id x (Velvet getArc list)."""
        k = self.k
        if self.rc or x > 0:
            tail = self.seq[x][-(k - 1):]
            cands = [y for i in self.alive for y in ((i, -i) if self.rc else (i,))]
            return sorted(y for y in cands if self.seq[y][:k - 1] == tail)
        # single strand, reverse view of node |x|: predecessors, negated
        head = self.seq[-x][:k - 1]
        return sorted(-p for p in self.alive if self.seq[p][-(k - 1):] == head)

    def count(self, x):
        return len(self.arcs(x))

    def mult(self, x, y):
        """Multiplicity of arc x -> y."""
        k = self.k
        if self.rc:
            w = self.seq[x][-k:] + self.seq[y][k - 1]
            return self.F1.get(w, 0) + self.F1.get(revc(w), 0)
        a, b = (x, y) if x > 0 else (-y, -x)       # real arc a -> b
        w = self.seq[a][-k:] + self.seq[b][k - 1]
        return self.F1.get(w, 0)

    def destroy(self, x):
        self.alive.discard(abs(x))

    # --- clipTips_pg
    def is_eligible_tip(self, index, cutoff):
        current = -index if self.rc or index > 0 else -index
        # NB: Velvet starts the walk at the twin of `index`.
        total, active, steps = 0, None, 0
        if self.arcs(index):
            return False
        while current is not None and self.count(-current) < 2 and self.count(current) < 2:
            steps += 1
            if steps > 10 * (len(self.alive) + 2):
                return False
            total += self.length(current)
            a = self.arcs(current)
            active = (current, a[0]) if a else None
            current = a[0] if a else None
        if total >= cutoff:
            return False
        if current is None:
            return True
        if self.count(-current) < 2:
            return False
        if self.mult(*active) == 1:
            return True
        mx = 0
        for y in self.arcs(-current):
            # arc of -current: -current -> y, i.e. the arc y' -> current with y' = -y
            mx = max(mx, self.mult(-current, y))   # twin arc of -y -> current: same multiplicity
        return mx > self.mult(*active)

    def clip_tips_pg(self):
        cutoff = 2 * self.k
        modified = True
        while modified:
            modified = False
            for index in sorted(self.alive):
                if index not in self.alive:
                    continue
                if self.is_eligible_tip(index, cutoff) or self.is_eligible_tip(-index, cutoff):
                    self.destroy(index)
                    modified = True

    # --- clipTipsHard
    def tip_length(self, node):
        current = -node
        if self.count(current) > 1:
            return self.length(node)
        length, steps = 0, 0
        while current is not None and self.count(-current) < 2 and self.count(current) < 2:
            steps += 1
            if steps > 10 * (len(self.alive) + 2):
                break
            length += self.length(current)
            a = self.arcs(current)
            current = a[0] if a else None
        return length

    def clip_tips_hard(self):
        cutoff = 2 * self.k
        modified = True
        while modified:
            modified = False
            for index in sorted(self.alive):
                if index not in self.alive:
                    continue
                if not self.arcs(index) and self.tip_length(index) < cutoff:
                    self.destroy(index)
                    modified = True
                elif not self.arcs(-index) and self.tip_length(-index) < cutoff:
                    self.destroy(index)
                    modified = True

    # --- coverage
    def remove_low_coverage(self, min_cov):
        for i in sorted(self.alive):
            if self.cov[i] // self.length(i) < min_cov:     # Coordinate / Coordinate
                self.destroy(i)

    def estimated_cov(self):
        nodes = sorted(self.alive, key=lambda i: self.cov[i] / self.length(i))
        if not nodes:
            return 0.0
        half = sum(self.length(i) for i in nodes) // 2
        s = 0
        for i in nodes:
            s += self.length(i)
            if s >= half:
                return self.cov[i] / self.length(i)
        return -1.0

    def kmers(self):
        out = set()
        for i in self.alive:
            s = self.seq[i]
            for j in range(len(s) - self.k + 1):
                out.add(s[j:j + self.k])
                if self.rc:
                    out.add(revc(s[j:j + self.k]))
        return out


def clean(reads, k, rc, clip, cov_cutoff, auto):
    F = occurrences(reads, k, k)
    F1 = occurrences(reads, k, k + 1)
    ks = w4.kmer_set(reads, k, rc)

    def build(s):
        return Graph(s, k, rc, F, F1)

    g = build(ks)
    if clip:
        g.clip_tips_pg()
        g.clip_tips_hard()
        ks = g.kmers()
        g = build(ks)
    cutoff = cov_cutoff
    if cutoff <= 0 and auto:
        cutoff = g.estimated_cov() / 2
    if cutoff > 0:
        g.remove_low_coverage(cutoff)
        g = build(g.kmers())
        g.remove_low_coverage(cutoff)
        if clip:
            g = build(g.kmers())
            g.clip_tips_hard()
    return g.kmers()


def expected(case):
    reads, k, rc, minc = case['reads'], case['k'], case['rc'], case['minContig']
    ks = clean(reads, k, rc, case['clip'], case['cov'], case['auto'])
    # contigs of the cleaned k-mer set exactly as wp4_oracle.dbg_expected does it
    a = w4.contigs_ba3m(ks)
    cyc = {c for c in a if c[:k - 1] == c[len(c) - (k - 1):] and w4.is_isolated_cycle(c, k, ks)}
    contigs = [c for c in w4.canon_pairs([c for c in a if c not in cyc], rc) if len(c) >= minc]
    cycles = sorted({(w4.cycle_class(c, k, rc), len(c)) for c in cyc if len(c) >= minc})
    ckm = set()
    for c in contigs + [c for c in cyc if len(c) >= minc]:
        for i in range(len(c) - k + 1):
            ckm.add(c[i:i + k])
            if rc:
                ckm.add(revc(c[i:i + k]))
    assembled = sum(1 for r in reads
                    if len(r) >= k and any(w4.velvet_fold(r)[i:i + k] in ckm for i in range(len(r) - k + 1)))
    return {"contigs": contigs, "cycles": [[c, n] for c, n in cycles], "assembled": assembled, "total": len(reads)}


def planted_case(rng):
    """Genome with planted errors: ground truth = the error-free genome (reads of one coverage level plus
    erroneous reads that create tips / low-coverage branches)."""
    k = rng.choice([5, 7, 9])  # odd: Velvet requires an odd word length
    n = rng.randint(30, 90)
    genome = ''.join(rng.choice('ACGT') for _ in range(n))
    reads, rl = [], rng.randint(k + 4, k + 14)
    cov = rng.randint(3, 6)
    for _ in range(cov * n // rl + 2):
        s = rng.randint(0, max(0, n - rl))
        r = genome[s:s + rl]
        if rng.random() < 0.25 and len(r) > 2:
            p = rng.randrange(len(r))
            r = r[:p] + rng.choice([c for c in 'ACGT' if c != r[p]]) + r[p + 1:]
        reads.append(r if rng.random() < 0.7 else revc(r))
    return genome, k, reads


def random_case(rng):
    mode = rng.random()
    if mode < 0.6:
        genome, k, reads = planted_case(rng)
    else:
        k = rng.choice([3, 5, 7])
        reads = [''.join(rng.choice('ACGT') for _ in range(rng.randint(k, k + 12))) for _ in range(rng.randint(1, 8))]
    rc = rng.random() < 0.7
    return {"reads": reads, "k": k, "rc": rc, "minContig": rng.choice([1, 1, k, 2 * k]),
            "clip": rng.random() < 0.75, "cov": rng.choice([0, 0, 1.5, 2, 2.5, 3, 4.5]),
            "auto": rng.random() < 0.3}


if __name__ == '__main__':
    out = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 20260525)
    cases = []
    for _ in range(n):
        c = random_case(rng)
        c['exp'] = expected(c)
        cases.append(c)
    json.dump(cases, open(out, 'w'))
    print(len(cases), 'cases;', sum(1 for c in cases if c['exp']['contigs'] == []), 'with no contig')
