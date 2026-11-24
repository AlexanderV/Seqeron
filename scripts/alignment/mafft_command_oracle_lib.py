#!/usr/bin/env python3
"""Shared helpers of the MAFFT command-line oracle generators (B09 WP70): deterministic sequence sets and one run of the installed `mafft` script
(Debian 7.505-1) on a directory that holds the input file and the extra files the options name.  Every run is repeated; a case whose two runs differ is reported."""
import gzip, hashlib, json, os, random, subprocess, sys, tempfile
from concurrent.futures import ThreadPoolExecutor

P = 'ACDEFGHIKLMNPQRSTVWY'
N = 'ACGT'
INPUT_NAME = 'in.fa'


def mut(rng, s, alpha, r, ind):
    out = []
    for c in s:
        x = rng.random()
        if x < ind / 2:
            continue
        if x < ind:
            out.append(c)
            out.append(rng.choice(alpha))
            continue
        out.append(rng.choice(alpha) if rng.random() < r else c)
    return ''.join(out) or s[:1] or alpha[0]


def family(rng, alpha, n, k):
    base = ''.join(rng.choice(alpha) for _ in range(n))
    r = rng.choice([0.05, 0.15, 0.3])
    ind = rng.choice([0.02, 0.05, 0.1])
    return [mut(rng, base, alpha, r, ind) for _ in range(k)]


def pick_alpha(rng):
    return P if rng.random() < 0.5 else N


def small_set(rng):
    return family(rng, pick_alpha(rng), rng.choice([10, 30, 60, 120]), rng.randint(3, 25))


def pair_set(rng):
    return family(rng, pick_alpha(rng), rng.choice([20, 60, 150]), 2)


def mid_set(rng):
    return family(rng, pick_alpha(rng), rng.choice([15, 40]), rng.randint(30, 60))


def tiny_set(rng):
    a = pick_alpha(rng)
    return [''.join(rng.choice(a) for _ in range(rng.randint(1, 14))) for _ in range(rng.randint(3, 9))]


def odd_set(rng):
    a = pick_alpha(rng)
    seqs = family(rng, a, rng.choice([20, 40, 80]), rng.randint(3, 12))
    extra = 'XBZJOU*#.@x' if a == P else 'NRYKMSWBDHVXnryx'
    out = []
    for s in seqs:
        s = list(s)
        for _ in range(rng.randrange(0, 4)):
            s[rng.randrange(len(s))] = rng.choice(extra)
        if rng.random() < 0.4:
            s = [c.lower() if rng.random() < 0.5 else c for c in s]
        for _ in range(rng.randrange(0, 3)):
            s.insert(rng.randrange(len(s) + 1), '-')
        out.append(''.join(s))
    return out


def gapped_set(rng):
    return [s[:3] + '-' * rng.randint(0, 3) + s[3:] for s in small_set(rng)]


def fasta(seqs, names=None, width=None):
    out = []
    for i, s in enumerate(seqs):
        nm = names[i] if names else 's%d' % i
        if width:
            body = '\n'.join(s[j:j + width] for j in range(0, len(s), width)) if s else ''
        else:
            body = s
        out.append('>%s\n%s\n' % (nm, body))
    return ''.join(out)


def esc(s):
    return s.replace('\\', '\\\\').replace('\n', '\\n').replace('\r', '\\r').replace('\b', '\\b').replace('\t', '\\t')


def run_once(args, text, extra=None, binary='mafft'):
    """One run of the script in a fresh directory: returns (rc, stdout, stderr, files) where files are all files that the run created or changed
    (name -> text), except the inputs."""
    extra = extra or {}
    with tempfile.TemporaryDirectory() as td:
        open(os.path.join(td, INPUT_NAME), 'wb').write(text.encode('utf-8'))
        before = {INPUT_NAME: text}
        for k, v in extra.items():
            open(os.path.join(td, k), 'wb').write(v.encode('utf-8'))
            before[k] = v
        env = dict(os.environ, TMPDIR=td)
        p = subprocess.run([binary] + args + [INPUT_NAME], capture_output=True, cwd=td, env=env)
        files = {}
        for f in sorted(os.listdir(td)):
            data = open(os.path.join(td, f), 'rb').read().decode('utf-8', errors='replace')
            if f in before and before[f] == data:
                continue
            files[f] = data
        return p.returncode, p.stdout.decode('utf-8', errors='replace'), p.stderr.decode('utf-8', errors='replace'), files


def run_checked(args, text, extra=None):
    """The run repeated: (result, deterministic)."""
    a = run_once(args, text, extra)
    b = run_once(args, text, extra)
    return a, a == b


def run_all(allcases, run_case, out_jsonl, tsv=None, threads=8):
    """Runs every case (run_case: case -> (record or None, id)), writes the JSONL and, when `tsv` is given, the TSV the tests read (gzip when the name ends in .gz): id, group, args, input, files, rc, out, sha256(err), err (kept for
    short texts and every 25th), outfiles.  Returns the list of dropped (non-deterministic) ids."""
    dropped = []
    ft = (gzip.open(tsv, 'wt', encoding='utf-8', newline='\n') if tsv.endswith('.gz') else open(tsv, 'w', encoding='utf-8', newline='\n')) if tsv else None
    full = 0
    done = 0
    with ThreadPoolExecutor(threads) as ex, open(out_jsonl, 'w') as fo:
        for res, bad in ex.map(run_case, allcases):
            done += 1
            if res is None:
                dropped.append(bad)
                continue
            fo.write(json.dumps(res) + '\n')
            if ft:
                sha = '-' if res['nostderr'] else hashlib.sha256(res['err'].encode('utf-8')).hexdigest()
                keep = (not res['nostderr']) and len(res['err']) < 200000 and (full % 25 == 0 or len(res['err']) < 400)
                full += 1
                cols = [res['id'], res['group'], json.dumps(res['args']), res['input'], json.dumps(res['files']), str(res['rc']), res['out'], sha, res['err'] if keep else '', json.dumps(res['outfiles'])]
                ft.write('\t'.join(esc(c) for c in cols) + '\n')
            if done % 200 == 0:
                print(done, '/', len(allcases), file=sys.stderr)
    if ft:
        ft.close()
    print('cases', len(allcases), 'written', len(allcases) - len(dropped), 'non-deterministic (dropped)', len(dropped), dropped[:10], file=sys.stderr)
    return dropped
