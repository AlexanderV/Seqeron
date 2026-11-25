"""Oracles for audit round 1, WP4 (A2-4, A2-5, A3-2, A3-1, A8-2, A5-4a).

References transcribed (fetched 2026-11-25 with a neutral User-Agent):
- Velvet 1.2.10 (dzerbino/velvet) src/readSet.c `velvetifySequence` (case-insensitive A/C/G/T kept, every
  other symbol -> 'N') and src/tightString.c `writeNucleotide` ('N' -> Adenine); Zerbino & Birney 2008
  (each k-mer recorded with its reverse complement; chains of 1-out/1-in arcs merged).
  DBG oracle: two independent contig constructions over the folded distinct k-mer set (A: Rosalind BA3M
  MaximalNonBranchingPaths over DeBruijn(k-mers); B: Velvet node-centric chain merge), strand pairs
  {U, revcomp(U)} reported once as the ordinally smaller string (the Seqeron contract).
- Musket 1.1 kmer.cpp `Kmer::encode` (non-ACGT -> 0 = A): spectrum + k-mer lookups see A; a base no
  correction touches keeps its input symbol (musket_oracle.py is the stage 1-3 transcription).
- QUAST quast_libs/qutils.py `correct_seq`, `split_by_ns`; basic_stats.py `GC_content` formula;
  N50.py `NG50_and_LG50`, `au_metric(numlist, reference_length)`.
- BWA bwaseqio.c `bwa_trim_read` with bwtaln.h `BWA_MIN_RDLEN 35`.
- scaffold_oracle.py (Huson et al. Alg. 5.1, ABySS overlap collapse) with a parametric unknown-gap length.

Usage: python3 -I wp4_oracle.py OUT.json [cases_per_op] [seed]
Writes [{op, ..inputs.., expected}]; run the C# side with the same case file and compare (see
docs/Evidence/ASSEMBLY-DBG-001-Evidence.md "R1-WP4" for the harness and the numbers).
"""
import json
import random
import re
import sys
from collections import defaultdict
from fractions import Fraction

sys.path.insert(0, __file__.rsplit('/', 1)[0])
import gaps_oracle as go          # noqa: E402  (QUAST split_by_ns verbatim)
import musket_oracle as mo        # noqa: E402
import scaffold_oracle as so      # noqa: E402

COMP = str.maketrans('ACGT', 'TGCA')


def revc(x):
    return x.translate(COMP)[::-1]


# ---------------------------------------------------------------- DBG (A2-4, A2-5)
def velvet_fold(read):
    """velvetifySequence (upper-case, non-ACGT -> N) + writeNucleotide (N -> Adenine)."""
    return ''.join(c if c in 'ACGT' else 'A' for c in read.upper())


def kmer_set(reads, k, rc):
    s = set()
    for r in reads:
        r = velvet_fold(r)
        for i in range(len(r) - k + 1):
            s.add(r[i:i + k])
    if rc:
        s |= {revc(x) for x in s}
    return s


def contigs_ba3m(ks):
    g = defaultdict(list)
    for p in sorted(ks):
        g[p[:-1]].append(p[1:])
    indeg = defaultdict(int)
    outdeg = defaultdict(int)
    nodes = set(g)
    for v, ws in g.items():
        outdeg[v] += len(ws)
        for w in ws:
            indeg[w] += 1
            nodes.add(w)
    one = lambda v: indeg[v] == 1 and outdeg[v] == 1
    paths, used = [], set()
    for v in sorted(nodes):
        if not one(v) and outdeg[v] > 0:
            for w in g[v]:
                path = [v, w]
                while one(w):
                    used.add(w)
                    w = g[w][0]
                    path.append(w)
                paths.append(path)
    for v in sorted(nodes):
        if one(v) and v not in used:
            cyc, w = [v], v
            while True:
                used.add(w)
                w = g[w][0]
                cyc.append(w)
                if w == v:
                    break
            paths.append(cyc)
    return [p[0] + ''.join(n[-1] for n in p[1:]) for p in paths]


def contigs_velvet(ks):
    """Node-centric: nodes = k-mers, arc x->y iff x[1:] == y[:-1]; merge A->B when A has one outgoing
    arc and B one ingoing arc (Velvet 'Simplification'). Rotation-invariant cycle handling."""
    succ = {x: sorted(y for y in ks if y[:-1] == x[1:]) for x in ks}
    pred = defaultdict(list)
    for x, ys in succ.items():
        for y in ys:
            pred[y].append(x)
    mergeable = lambda a, b: len(succ[a]) == 1 and succ[a][0] == b and len(pred[b]) == 1
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


def cycle_class(c, k, rc):
    """An isolated cycle is spelled from its smallest node, and its reverse-complement cycle from another
    start, so cycles are compared as rotation classes (of the word and, with rc, of its reverse complement)."""
    w = c[:len(c) - (k - 1)]
    rots = lambda x: min(x[i:] + x[:i] for i in range(len(x)))
    return min(rots(w), rots(revc(w))) if rc else rots(w)


def canon_pairs(contigs, rc):
    if not rc:
        return sorted(contigs)
    return sorted({min(c, revc(c)) for c in contigs})


def is_isolated_cycle(c, k, ks):
    """True when every k-mer of the spelled path has in- and out-degree 1 in the k-mer graph (a closed loop
    of 1-in-1-out nodes: no k-mer outside the cycle overlaps it)."""
    kms = [c[i:i + k] for i in range(len(c) - k + 1)]
    for x in kms:
        if sum(1 for y in ks if y[:-1] == x[1:]) != 1 or sum(1 for y in ks if y[1:] == x[:-1]) != 1:
            return False
    return True


def dbg_expected(reads, k, min_contig, rc):
    ks = kmer_set(reads, k, rc)
    a = contigs_ba3m(ks)
    b = contigs_velvet(ks)
    # A and B agree on every contig except the spelling start of an isolated cycle (reported in "agree").
    cyc = {c for c in a if c[:k - 1] == c[len(c) - (k - 1):] and is_isolated_cycle(c, k, ks)}
    ca = canon_pairs([c for c in a if c not in cyc], rc)
    cb = canon_pairs([c for c in b if not (c[:k - 1] == c[len(c) - (k - 1):] and is_isolated_cycle(c, k, ks))], rc)
    cyc_b = {c for c in b if c[:k - 1] == c[len(c) - (k - 1):] and is_isolated_cycle(c, k, ks)}
    agree = (ca == cb) and ({cycle_class(c, k, rc) for c in cyc} == {cycle_class(c, k, rc) for c in cyc_b})
    contigs = [c for c in canon_pairs([c for c in a if c not in cyc], rc) if len(c) >= min_contig]
    cycles = sorted({(cycle_class(c, k, rc), len(c)) for c in cyc if len(c) >= min_contig})
    ckm = set()
    for c in contigs + [c for c in cyc if len(c) >= min_contig]:
        for i in range(len(c) - k + 1):
            ckm.add(c[i:i + k])
            if rc:
                ckm.add(revc(c[i:i + k]))
    assembled = sum(1 for r in reads
                    if len(r) >= k and any(velvet_fold(r)[i:i + k] in ckm for i in range(len(r) - k + 1)))
    return {"contigs": contigs, "cycles": [[c, n] for c, n in cycles], "assembled": assembled,
            "total": len(reads), "agree": agree}


# ---------------------------------------------------------------- QUAST (A3-2, A3-1)
def correct_seq_quast(seq):
    """QUAST qutils.correct_seq (verbatim) after mapping the symbols QUAST would reject to N."""
    corr = seq.upper()
    dic = {'M': 'N', 'K': 'N', 'R': 'N', 'Y': 'N', 'W': 'N', 'S': 'N', 'V': 'N', 'B': 'N', 'H': 'N', 'D': 'N'}
    pat = "(%s)" % "|".join(map(re.escape, dic.keys()))
    corr = re.sub(pat, lambda m: dic[m.group()], corr)
    corr = re.sub(r'[^ACGTN]', 'N', corr)     # QUAST: file skipped; Seqeron: lenient (seqtk gap reading)
    return corr


def ng50_and_lg50(numlist, reference_length, percentage=50.0):
    s = reference_length
    limit = reference_length * (100.0 - percentage) / 100.0
    lg = 0
    for l in numlist:
        s -= l
        lg += 1
        if s <= limit:
            return l, lg
    return None, None


def au_metric(numlist, reference_length=None):
    try:
        assert len(numlist) > 0
        denum = float(reference_length) if reference_length else float(sum(numlist))
        assert denum > 0.0
        return float(sum([n ** 2 for n in numlist])) / denum
    except AssertionError:
        return None


def stats_expected(seqs, min_gap, min_contig):
    corr = [correct_seq_quast(s) for s in seqs]
    lengths = sorted((len(s) for s in seqs), reverse=True)
    total = sum(lengths)
    n50, l50 = ng50_and_lg50(lengths, total, 50.0)
    n90, l90 = ng50_and_lg50(lengths, total, 90.0)
    gc = acgt = 0
    gaps_n = gap_len = 0
    for c in corr:
        acgt += len(c) - c.count('N')
        gc += c.count('G') + c.count('C')
        runs = [(m.start(), m.end() - m.start()) for m in re.finditer(r'N+', c)]
        gaps_n += len(runs)
        gap_len += sum(r[1] for r in runs)
    gap_list, contigs, extracted = [], [], []
    for s, c in zip(seqs, corr):
        gap_list += [[m.start(), m.end() - 1, m.end() - m.start()] for m in re.finditer(r'N+', c)
                     if m.end() - m.start() >= min_gap]
        pieces = go.quast_broken_pieces_exact(c, min_gap)
        contigs.append([[a, b - 1] for a, b in pieces])
        extracted += [s[a:b] for a, b in pieces if b - a >= max(min_contig, 1) or (min_contig == 0 and b > a)]
    return {"total": total, "n50": n50 or 0, "l50": l50 or 0, "n90": n90 or 0, "l90": l90 or 0,
            "gc": gc / acgt if acgt else 0.0, "gaps": gaps_n, "gapLen": gap_len,
            "noGaps": total - gap_len, "gapList": gap_list, "contigs": contigs, "extracted": extracted}


# ---------------------------------------------------------------- BWA (A8-2)
def bwa_trim_read(trim_qual, qual, min_rdlen=35):
    """bwaseqio.c bwa_trim_read: returns the new read length (p->len)."""
    n = len(qual)
    s, mx, max_l = 0, 0, n
    if trim_qual < 1:
        return n
    l = n - 1
    while l >= min_rdlen:
        s += trim_qual - (ord(qual[l]) - 33)
        if s < 0:
            break
        if s > mx:
            mx, max_l = s, l
        l -= 1
    return max_l


# ---------------------------------------------------------------- random cases
def rand_seq(rng, n, alphabet, weights=None):
    return ''.join(rng.choices(alphabet, weights=weights, k=n))


def gen(n, seed):
    rng = random.Random(seed)
    cases = []
    dna_noise = "ACGTacgtNnRYKMSWBDHVU-*x"
    wts = [20] * 4 + [5] * 4 + [3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    for _ in range(n):                                       # DBG
        k = rng.choice([3, 4, 5, 7, 9])
        genome = rand_seq(rng, rng.randint(20, 80), "ACGT")
        reads = []
        for _ in range(rng.randint(1, 8)):
            a = rng.randint(0, len(genome) - 1)
            r = genome[a:a + rng.randint(k, 30)]
            if rng.random() < 0.3:
                r = revc(r)
            r = ''.join(rng.choices(dna_noise, weights=wts, k=1)[0] if rng.random() < 0.06 else c for c in r)
            reads.append(r)
        rc = rng.random() < 0.7
        mc = rng.choice([1, 1, 5, 12])
        e = dbg_expected(reads, k, mc, rc)
        cases.append({"op": "dbg", "reads": reads, "k": k, "minContig": mc, "rc": rc, "exp": e})
    for _ in range(max(1, n // 4)):                          # CORRECT (Musket, N -> A)
        genome = rand_seq(rng, rng.randint(60, 140), "ACGT")
        k = rng.choice([5, 7, 9])
        reads = []
        for _ in range(rng.randint(12, 30)):
            a = rng.randint(0, len(genome) - k)
            r = list(genome[a:a + rng.randint(k, 40)])
            for i in range(len(r)):
                x = rng.random()
                if x < 0.03:
                    r[i] = rng.choice("ACGT")
                elif x < 0.05:
                    r[i] = rng.choice("NnRYKMSWBDHVU-*")
                elif x < 0.07:
                    r[i] = r[i].lower()
            reads.append(''.join(r))
        w = rng.choice([2, 3])
        cases.append({"op": "correct", "reads": reads, "k": k, "w": w,
                      "exp": {"reads": musket_fold_correct(reads, k, w)}})
    alpha = "ACGTacgtNnRYKMSWBDHVU-*"
    awts = [20] * 8 + [6, 2] + [1] * 13
    for _ in range(n):                                       # STATS / GAPS
        seqs = []
        for _ in range(rng.randint(1, 4)):
            length = rng.randint(0, 60)
            s, i = [], 0
            while i < length:
                if rng.random() < 0.3:
                    run = rng.choice([1, 1, 2, 3, 5, 9, 10, 11, 15, 25])
                    s.extend(rng.choice("NnRYKMSWBDHVU-*") for _ in range(run))
                    i += run
                else:
                    s.append(rng.choice("ACGTacgt"))
                    i += 1
            seqs.append(''.join(s))
        mg = rng.choice([1, 2, 3, 5, 10, 11])
        mc = rng.choice([0, 1, 2, 3, 5, 8])
        cases.append({"op": "stats", "seqs": seqs, "minGap": mg, "minContig": mc, "asGap": True,
                      "exp": stats_expected(seqs, mg, mc)})
    for _ in range(n):                                       # auN / auNG
        lens = [rng.choice([0, 1, 5, 50, 1000, 10 ** 6, rng.randint(1, 3 * 10 ** 9 // 1000)])
                for _ in range(rng.randint(0, 12))]
        g = rng.choice([None, rng.randint(1, 5 * 10 ** 7), sum(lens) or 1, 3 * 10 ** 9])
        exact = None if not lens or (g is None and sum(lens) == 0) else \
            float(Fraction(sum(x * x for x in lens), g if g else sum(lens)))
        cases.append({"op": "aun", "lengths": lens, "genome": g,
                      "exp": {"v": 0.0 if exact is None else exact,
                              "quast": au_metric(lens, g) if lens else None}})
    for _ in range(n):                                       # BWA trim
        reads = []
        for _ in range(rng.randint(1, 6)):
            L = rng.choice([0, 1, 10, 34, 35, 36, 40, 60, 100])
            q = ''.join(chr(33 + rng.choice([2, 5, 10, 15, 20, 30, 38, rng.randint(0, 41)])) for _ in range(L))
            reads.append(['A' * L, q])
        tq = rng.choice([0, 1, 5, 10, 15, 20, 30])
        cases.append({"op": "trim", "reads": reads, "q": tq, "floor": 35,
                      "exp": {"kept": [bwa_trim_read(tq, r[1]) for r in reads]}})
    for _ in range(max(1, n // 2)):                          # Scaffold unknownGapLength
        nc = rng.randint(2, 5)
        contigs = [rand_seq(rng, rng.randint(8, 25), "ACGT") for _ in range(nc)]
        links = []
        for _ in range(rng.randint(1, 6)):
            a, b = rng.sample(range(nc), 2)
            links.append([a, b, rng.choice([-30, -4, 0, 0, 5, 12]), rng.random() < .3, rng.random() < .3, rng.randint(1, 4)])
        ugl = rng.choice([1, 7, 100])
        so.UNKNOWN = ugl
        cases.append({"op": "scaffold", "contigs": contigs, "links": links, "ugl": ugl,
                      "exp": {"scaffolds": so.scaffold(contigs, [tuple(l) for l in links])}})
    return cases


def musket_fold_correct(reads, k, w):
    """Musket verbatim encoding: k-mers see A for every non-ACGT; a base no correction assigns keeps its symbol."""
    raw = [r.upper() for r in reads]
    folded = [velvet_fold(r) for r in raw]
    sp = mo.Spectrum(folded, k, True)
    out = []
    for r, f in zip(raw, folded):
        touched = set()
        res = para_correct_tracked(sp, w, f, touched)
        out.append(''.join(res[i] if (i in touched) else r[i] for i in range(len(r))))
    return out


def para_correct_tracked(sp, w, read, touched):
    """musket_oracle.para_correct with assignment tracking (identical control flow)."""
    seq = list(read)
    for _ in range(2):
        ef, corr = mo.two_sides(sp, w, seq, True)
        if ef:
            return ''.join(seq)
        if not corr:
            break
        for i, b in corr:
            seq[i] = mo.DEC[b]
            touched.add(i)
    for nerr in range(1, 5):
        pp = 1
        for _ in range(2):
            ef, corr = mo.one_side(sp, w, seq, nerr, 4 - nerr + 1, pp & 1, True)
            pp += 1
            if ef:
                return ''.join(seq)
            if not corr:
                break
            for i, b in corr:
                seq[i] = mo.DEC[b]
                touched.add(i)
        mv, votes = mo.vote(sp, w, seq)
        if mv == 0:
            return ''.join(seq)
        for i, b in mo.fix_by_votes(votes, mv, 3, sp.k, len(seq)):
            seq[i] = mo.DEC[b]
            touched.add(i)
    return ''.join(seq)


if __name__ == '__main__':
    out = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 1500
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 20261125
    cs = gen(n, seed)
    json.dump(cs, open(out, 'w'))
    disagree = sum(1 for c in cs if c['op'] == 'dbg' and not c['exp']['agree'])
    print(f"{len(cs)} cases; DBG method A vs B disagree on {disagree} (rotation of isolated cycles / rc pairs)")
