#!/usr/bin/env python3
"""ALIGN-MULTI-001 (audit A2-27, B09 F120) oracle: the EMBOSS 6.6.0 sequence-set alignment formats of more than two rows (fasta, msf, clustal, mega, meganon, nexus, nexusnon, phylip, phylipnon,
selex, treecon), byte for byte.

Seqeron API: EmbossAlignmentWriter.Format(format, rows, options).   Reference (run here): the installed `seqret -sequence set.fa -osformat <format> -auto -stdout` on gapped FASTA sets of 1 - 9 rows
(DNA / IUPAC / protein; rows of different length; shared or distinct descriptions; names of 1 - 50 characters), with a faked run date (scripts/alignment/faketime.c).
Output TSV (tests/.../TestData/Alignment/emboss_seqout_oracle.tsv):  id  format  type  fasta  date  file      (type = N or P as seqret declared it: the MSF header says which)
Usage: python3 scripts/alignment/generate_emboss_seqout_oracle.py OUTFILE [setsPerFormat] [seed]      (needs EMBOSS seqret and gcc on PATH)
"""
import datetime, os, random, re, subprocess, sys, tempfile

FORMATS = ['fasta', 'msf', 'clustal', 'mega', 'meganon', 'nexus', 'nexusnon', 'phylip', 'phylipnon', 'selex', 'treecon']


def esc(t):
    return t.replace('\\', '\\\\').replace('\n', '\\n').replace('\t', '\\t')


def make_set(rnd):
    typ = rnd.choice('np')
    alpha = 'ACGT' if typ == 'n' else 'ACDEFGHIKLMNPQRSTVWY'
    if rnd.random() < 0.25:
        alpha += 'RYNKM' if typ == 'n' else 'BZX'
    n = rnd.choice([1, 2, 2, 3, 3, 4, 5, 7, 9])
    length = rnd.choice([1, 5, 12, 33, 50, 51, 60, 61, 100, 130])
    base = ''.join(rnd.choice(alpha) for _ in range(length))
    desc_mode = rnd.choice(['none', 'same', 'distinct'])
    shared = ' '.join(rnd.choice(['alpha', 'beta', 'gene', 'kinase', 'x1']) for _ in range(rnd.randint(1, 3)))
    names, seen = [], set()
    rows = []
    for i in range(n):
        while True:
            nm = rnd.choice(['s', 'seq', 'My_seq', 'sp', 'x.y']) + ''.join(rnd.choice('abcdefghijklmnopqrstuvwxyz0123456789_.-') for _ in range(rnd.choice([0, 1, 4, 9, 14, 25, 45])))
            if nm not in seen:
                seen.add(nm); break
        row = list(base)
        for j in range(len(row)):
            if rnd.random() < 0.08: row[j] = rnd.choice(alpha)
            elif rnd.random() < 0.08: row[j] = '-'
        r = ''.join(row)
        if rnd.random() < 0.3: r = r[:max(1, len(r) - rnd.randint(1, 6))]
        if rnd.random() < 0.15: r = r.lower()
        if set(r) == {'-'}: r = 'A' + r[1:]
        if r[0] == '-' and rnd.random() < 0.5: r = rnd.choice(alpha) + r[1:]
        d = '' if desc_mode == 'none' else (shared if desc_mode == 'same' else ' '.join(rnd.choice(['alpha', 'beta', 'gene', 'kinase', 'x1', 'zeta']) for _ in range(rnd.randint(1, 3))))
        rows.append((nm, d, r))
    return rows


def main(outfile, per, seed):
    rnd = random.Random(seed)
    out = []
    n = 0
    with tempfile.TemporaryDirectory() as d:
        shim = os.path.join(d, 'faketime.so')
        subprocess.run(['gcc', '-shared', '-fPIC', '-o', shim, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'faketime.c')], check=True)
        for case in range(per):
            rows = make_set(rnd)
            fasta = ''.join(f'>{nm}' + (f' {de}' if de else '') + f'\n{r}\n' for nm, de, r in rows)
            open(os.path.join(d, 'm.fa'), 'w').write(fasta)
            back = subprocess.run(['seqret', '-sequence', 'm.fa', '-osformat', 'fasta', '-auto', '-stdout'], cwd=d, capture_output=True).stdout.decode('latin-1')
            if ''.join(back.split('\n')) != ''.join(re.sub(r'(?m)^(>[^ \n]*).*$', r'\1', fasta).split('\n')) and re.sub(r'\s+', '', ''.join(l for l in back.split('\n') if not l.startswith('>'))) != ''.join(r for _, _, r in rows):
                continue    # seqret converts symbols on reading (e.g. X -> N in a nucleotide set): not the writer's business
            epoch = rnd.randint(978307200, 2051222400)
            env = dict(os.environ, TZ='UTC', LD_PRELOAD=shim, FAKE_EPOCH=str(epoch))
            for fmt in FORMATS:
                res = subprocess.run(['seqret', '-sequence', 'm.fa', '-osformat', fmt, '-auto', '-stdout'], cwd=d, env=env, capture_output=True)
                text = res.stdout.decode('latin-1')
                if res.returncode != 0 or not text or re.search(r'[^\x00-\x7f]', text):
                    continue
                m = re.search(r'Type: ([NP])', text)
                if fmt == 'msf':
                    typ = m.group(1)
                else:
                    t2 = subprocess.run(['seqret', '-sequence', 'm.fa', '-osformat', 'msf', '-auto', '-stdout'], cwd=d, env=env, capture_output=True).stdout.decode('latin-1')
                    typ = re.search(r'Type: ([NP])', t2).group(1)
                date = datetime.datetime.fromtimestamp(epoch, datetime.timezone.utc).replace(tzinfo=None)
                out.append('\t'.join([str(n), fmt, typ, esc(fasta), date.strftime('%Y-%m-%dT%H:%M:%S'), esc(text)]))
                n += 1
    with open(outfile, 'w', encoding='utf-8') as f:
        f.write('id\tformat\ttype\tfasta\tdate\tfile\n' + '\n'.join(out) + '\n')
    print(len(out), 'rows')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 110, int(sys.argv[3]) if len(sys.argv) > 3 else 20261120)
