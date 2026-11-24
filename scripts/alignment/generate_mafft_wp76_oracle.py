#!/usr/bin/env python3
"""Builds the oracle of the options that WP76 of batch B09 ported: every case is one run of the installed script `mafft` (Debian 7.505-1), repeated twice (a case whose runs differ is dropped); the oracle holds return code, standard output,
standard error and every file the run wrote.

Groups (>= 200 sets each, scale with MC_SCALE):
  dir         --adjustdirection / --adjustdirectionaccurately (makedirectionlist -m, setdirection): nucleotide families whose members are partly reverse complemented, upper / lower case, N; protein sets ("Not necessary!"); with the options of the script
  dir-add     the same with --add / --addfragments / --addfull / --addlong / --addtotop (makedirectionlist -I n -F)
  dir-merge   the same with --merge (setdirection -H n: the directions inside a subalignment must agree)
  anchors     --anchors file (the external anchors of disttbfast: ExtAnch, Falign_givenanchors) with --skipanchorsremoterthan
Usage: python3 scripts/alignment/generate_mafft_wp76_oracle.py [group ...]   env MC_SCALE (default 1), MC_OUT (JSONL path), MC_TSV (TSV, gzip when it ends in .gz), MC_THREADS
"""
import os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mafft_command_oracle_lib as L
import generate_mafft_command_oracle as C
import generate_mafft_wp75_oracle as W

SCALE = float(os.environ.get('MC_SCALE', '1'))
OUT = os.environ.get('MC_OUT', '/tmp/mafft_wp76_oracle.jsonl')
THREADS = int(os.environ.get('MC_THREADS', '8'))
COMP = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')


def revcomp(s):
    return s.translate(COMP)[::-1]


def dna_family(rng, k, length):
    base = ''.join(rng.choice('ACGT') for _ in range(length))
    r = rng.choice([0.03, 0.08, 0.15])
    ind = rng.choice([0.0, 0.02, 0.05])
    return base, [L.mut(rng, base, 'ACGT', r, ind) for _ in range(k)]


def messy(rng, s):
    s = list(s)
    if rng.random() < 0.25:
        s = [c.lower() if rng.random() < 0.5 else c for c in s]
    if rng.random() < 0.15:
        for _ in range(rng.randint(1, 3)):
            s[rng.randrange(len(s))] = 'N' if rng.random() < 0.7 else 'n'
    return ''.join(s)


DIR_EXTRAS = [[], [], [], ['--reorder'], ['--anysymbol'], ['--retree', '1'], ['--maxiterate', '2'], ['--localpair'], ['--clustalout'], ['--thread', '2'], ['--6merpair'], ['--nuc'], ['--globalpair'], ['--op', '1.0'], ['--inputorder'],
              ['--genafpair'], ['--quiet'], ['--maxiterate', '1', '--localpair'], ['--phylipout'], ['--nofft'], ['--treeout', '--distout'], ['--preservecase']]


def dir_case(rng, k, group, mode):
    protein = rng.random() < 0.08
    cut0 = None
    length = rng.choice([40, 80, 160, 320])
    n = rng.choice([2, 3, 4, 6, 9, 14, 25]) if group != 'dir-merge' else rng.choice([4, 6, 9, 14, 25])
    if protein:
        base = ''.join(rng.choice('ACDEFGHIKLMNPQRSTVWY') for _ in range(length))
        seqs = [L.mut(rng, base, 'ACDEFGHIKLMNPQRSTVWY', 0.15, 0.02) for _ in range(n)]
    else:
        _, seqs = dna_family(rng, n, length)
        if group == 'dir-merge':
            # two subalignments (consecutive members), each in one direction: the table is consistent unless a member is turned (a quarter of the sets)
            # (every subalignment is a real alignment of the installed mafft: the script stops on unaligned subalignments before any direction is looked at)
            cut0 = rng.randint(1, n - 1)
            flip = [rng.random() < 0.5, rng.random() < 0.5]
            seqs = [messy(rng, s) for s in seqs]
            rows = []
            for part in (seqs[:cut0], seqs[cut0:]):
                rc_, aln, _, _ = L.run_once(['--quiet'], L.fasta(part))
                rows += [''.join(x.split('\n')[1:]) for x in aln.split('>')[1:]] if rc_ == 0 else part
            seqs = [revcomp(s) if flip[0 if i < cut0 else 1] else s for i, s in enumerate(rows)]
            if rng.random() < 0.25:
                j = rng.randrange(n)
                seqs[j] = revcomp(seqs[j])
        else:
            seqs = [messy(rng, revcomp(s) if rng.random() < 0.4 else s) for s in seqs]
    names = ['seq%d' % i if rng.random() < 0.85 else 'a very long sequence name that goes on and on %d' % i for i in range(n)]
    ex = list(rng.choice(DIR_EXTRAS))
    args = [mode] + ex
    files = {}
    text = L.fasta(seqs, names)
    if group == 'dir-add':
        nnew = rng.randint(1, min(4, n - 1)) if n > 2 else 1
        old, new = seqs[:-nnew], seqs[-nnew:]
        kind = rng.choice(['--add', '--addfragments', '--addfull', '--addlong', '--addtotop'])
        rc, aln, err, f = L.run_once(['--quiet', '--adjustdirection'], L.fasta(old, names[:len(old)]))
        text = aln if rc == 0 and aln.count('>') == len(old) else L.fasta(old, names[:len(old)])
        if kind == '--add' or rng.random() < 0.5:
            new = [s[rng.randrange(max(1, len(s) // 3)):] for s in new] if kind != '--add' else new
        files['new.fa'] = L.fasta(new, ['new%d' % i for i in range(len(new))])
        args = [mode, kind, 'new.fa'] + ex
    elif group == 'dir-merge':
        # subalignments: consecutive groups; the table lists the members of every group (1-origin), one per line
        cut = cut0 if cut0 is not None else rng.randint(1, n - 1)
        groups = [list(range(1, cut + 1)), list(range(cut + 1, n + 1))]
        files['tab.txt'] = '\n'.join(' '.join(map(str, g)) for g in groups if g) + '\n'
        args = [mode, '--merge', 'tab.txt'] + ex
    return ('%s-%04d' % (group, k), group, args, text, files)


ANCHOR_EXTRAS = [[], [], [], ['--skipanchorsremoterthan', '5'], ['--skipanchorsremoterthan', '50'], ['--retree', '1'], ['--retree', '2'], ['--nofft'], ['--6merpair'], ['--reorder'], ['--clustalout'], ['--thread', '2'], ['--anysymbol'],
                 ['--op', '1.0'], ['--bl', '45'], ['--treeout'], ['--maxiterate', '2'], ['--localpair'], ['--oneiteration']]


def anchor_case(rng, k):
    """--anchors: stretches of the sequences of a family that are identical (found with difflib), a random part of them, with overlaps and errors added."""
    import difflib
    protein = rng.random() < 0.5
    alpha = 'ACDEFGHIKLMNPQRSTVWY' if protein else 'ACGT'
    length = rng.choice([40, 80, 160, 300])
    n = rng.choice([2, 3, 4, 6, 9])
    base = ''.join(rng.choice(alpha) for _ in range(length))
    seqs = [L.mut(rng, base, alpha, rng.choice([0.05, 0.12, 0.25]), rng.choice([0.01, 0.04])) for _ in range(n)]
    lines = []
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.3:
                continue
            m = difflib.SequenceMatcher(None, seqs[i], seqs[j], autojunk=False)
            for blk in m.get_matching_blocks():
                if blk.size >= rng.choice([3, 5, 8]) and rng.random() < 0.7:
                    take = rng.randint(max(2, blk.size // 2), blk.size)
                    off = rng.randint(0, blk.size - take)
                    a, b = blk.a + off, blk.b + off
                    lines.append('%d %d %d %d %d %d %d' % (i + 1, j + 1, a + 1, a + take, b + 1, b + take, rng.choice([5, 10, 10, 20, 40, take * 5])))
    rng.shuffle(lines)
    if rng.random() < 0.25 and lines:
        # an overlapping / crossing / far anchor
        f = rng.choice(lines).split()
        f[2], f[3] = str(max(1, int(f[2]) - rng.randint(0, 4))), str(min(len(seqs[int(f[0]) - 1]), int(f[3]) + rng.randint(0, 4)))
        f[4], f[5] = str(max(1, int(f[4]) - rng.randint(0, 4))), str(min(len(seqs[int(f[1]) - 1]), int(f[5]) + rng.randint(0, 4)))
        lines.append(' '.join(f))
    if rng.random() < 0.06:
        lines.append(rng.choice(['9 10 1 3 1 3 5', '2 1 1 3 1 3 5', '1 2 1 3 1 1 5', '1 1 1 3 1 3 5']))
    if rng.random() < 0.15:
        lines.insert(rng.randrange(len(lines) + 1), '# a comment')
    files = {'anchors.txt': '\n'.join(lines) + ('\n' if lines else '')}
    args = ['--anchors', 'anchors.txt'] + list(rng.choice(ANCHOR_EXTRAS))
    text = L.fasta(seqs, ['s%d' % i for i in range(n)])
    return ('anchors-%04d' % k, 'anchors', args, text, files)


def build(rng, want):
    out = []

    def add(group, fn, n):
        if not want or group in want:
            for k in range(n):
                out.append(fn(k))

    n = int(210 * SCALE)
    add('dir', lambda k: dir_case(rng, k, 'dir', rng.choice(['--adjustdirection', '--adjustdirectionaccurately'])), n)
    add('dir-add', lambda k: dir_case(rng, k, 'dir-add', rng.choice(['--adjustdirection', '--adjustdirectionaccurately'])), n)
    add('dir-merge', lambda k: dir_case(rng, k, 'dir-merge', rng.choice(['--adjustdirection', '--adjustdirectionaccurately'])), n)
    add('anchors', lambda k: anchor_case(rng, k), n)
    return out


def run_case(case):
    res, bad = C.run_case(case)
    if res is not None and 'Segmentation fault' in res['err']:
        return None, case[0]
    return res, bad


def main():
    rng = random.Random(7600116)
    want = set(sys.argv[1:])
    allcases = build(rng, want)
    L.run_all(allcases, run_case, OUT, os.environ.get('MC_TSV'), THREADS)


if __name__ == '__main__':
    main()
