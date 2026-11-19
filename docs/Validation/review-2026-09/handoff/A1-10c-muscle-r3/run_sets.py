#!/usr/bin/env python3
"""ROUND 3: runs MUSCLE 3.8.1551 (/opt/muscle3/muscle3) on every set of inputs/sets.tsv and writes out/muscle3_r3_oracle.tsv (never aborts: a failing set is recorded).

Output columns: id kind type options rc names rows tree1 tree2 dist1 dist2 weights usedtree stdout
  rc      = exit status of muscle (FAILED<n> / TIMEOUT / NOOUT when no output file was written)
  names   = the FASTA labels of the written rows, ',' joined, in the order muscle wrote them
  rows    = the written rows (';' joined, in written order); for -maketree the tree file text; for -spscore / -ppscore empty
  tree1 / tree2 / dist1 / dist2 / weights = the text of -tree1 / -tree2 / -distmx1 / -distmx2 / -computeweights files, usedtree = the tree handed to -usetree
  stdout  = stdout + stderr (trimmed to 400 chars, '-quiet' runs keep only the program's own messages)
Every tab / newline inside a field is written as \\t / \\n.
"""
import os, subprocess, tempfile

BIN = os.environ.get('MUSCLE3', '/opt/muscle3/muscle3')
HERE = os.path.dirname(os.path.abspath(__file__))
MATRICES = {'@blosum62': 'blosum62.txt', '@blosum62nostar': 'blosum62nostar.txt', '@pam250': 'pam250.txt'}
FILE_KEYS = {'tree1': '-tree1', 'tree2': '-tree2', 'distmx1': '-distmx1', 'distmx2': '-distmx2', 'computeweights': '-computeweights'}


def dec(b):
    """Console bytes -> text without ever failing: the program prints non-UTF-8 bytes (e.g. 0xff in the 'Matrix is not symmetrical' warning); they are written as \\xNN."""
    return b.decode('utf-8', errors='backslashreplace')


def esc(s): return s.replace('\\', '\\\\').replace('\t', '\\t').replace('\n', '\\n').replace('\r', '')


def write_fasta(path, rows, prefix):
    with open(path, 'w') as f:
        for i, r in enumerate(rows):
            f.write(f'>{prefix}{i}\n{r}\n')


def read_fasta(path):
    names, rows, name = [], {}, None
    for line in open(path, errors='replace'):
        line = line.rstrip('\n')
        if line.startswith('>'):
            name = line[1:].strip(); names.append(name); rows[name] = ''
        elif name:
            rows[name] += line.strip()
    return names, [rows[n] for n in names]


def parse_options(opts, d, typ):
    """-> (argv list, {file key: path}, wants_own_tree)"""
    args, files, own = [], {}, False
    for kv in opts.split(',') if opts else []:
        k, _, v = kv.partition('=')
        if k in FILE_KEYS:
            files[k] = os.path.join(d, k + '.txt'); args += [FILE_KEYS[k], files[k]]
        elif v == '@own':
            own = True; args += ['-' + k, os.path.join(d, 'own.nwk')]
        elif v in MATRICES:
            args += ['-' + k, os.path.join(HERE, 'inputs', 'matrices', MATRICES[v])]
        elif v == '':
            args += ['-' + k]
        else:
            args += ['-' + k, v]
    return args, files, own


def run_one(sid, kind, typ, data, opts):
    with tempfile.TemporaryDirectory() as d:
        out = os.path.join(d, 'o.fa')
        args, files, own = parse_options(opts, d, typ)
        base = [BIN, '-seqtype', typ]
        usedtree = ''
        if kind == 'align':
            seqs = data.split(';')
            write_fasta(os.path.join(d, 'x.fa'), seqs, 's')
            if own:
                pre = subprocess.run(base + ['-in', os.path.join(d, 'x.fa'), '-clusteronly', '-tree1', os.path.join(d, 'own.nwk'), '-quiet'], capture_output=True, cwd=d, timeout=600)
                if not os.path.exists(os.path.join(d, 'own.nwk')):
                    return [sid, kind, typ, opts, 'FAILED(pre)', '', '', '', '', '', '', '', '', esc(dec(pre.stdout + pre.stderr)[:400])]
                usedtree = open(os.path.join(d, 'own.nwk'), errors='replace').read()
            cmd = base + ['-in', os.path.join(d, 'x.fa'), '-out', out, '-quiet'] + args
        elif kind in ('profile', 'sw', 'ppscore'):
            a, b = data.split('|')
            write_fasta(os.path.join(d, 'a.fa'), a.split(';'), 'a'); write_fasta(os.path.join(d, 'b.fa'), b.split(';'), 'b')
            flag = {'profile': '-profile', 'sw': '-sw', 'ppscore': '-ppscore'}[kind]
            cmd = base + [flag, '-in1', os.path.join(d, 'a.fa'), '-in2', os.path.join(d, 'b.fa'), '-out', out, '-quiet'] + args
        elif kind == 'profdb':
            a, b = data.split('|')
            write_fasta(os.path.join(d, 'a.fa'), a.split(';'), 'a'); write_fasta(os.path.join(d, 'b.fa'), b.split(';'), 'b')
            cmd = base + ['-profdb', '-in1', os.path.join(d, 'a.fa'), '-in2', os.path.join(d, 'b.fa'), '-out', out, '-quiet'] + args
        elif kind == 'refine':
            write_fasta(os.path.join(d, 'x.fa'), data.split(';'), 's')
            cmd = base + ['-refine', '-in', os.path.join(d, 'x.fa'), '-out', out, '-quiet'] + args
        elif kind == 'refinew':
            write_fasta(os.path.join(d, 'x.fa'), data.split(';'), 's')
            cmd = base + ['-refinew', '-in', os.path.join(d, 'x.fa'), '-out', out, '-quiet'] + args
        elif kind == 'maketree':
            write_fasta(os.path.join(d, 'x.fa'), data.split(';'), 's')
            cmd = base + ['-maketree', '-in', os.path.join(d, 'x.fa'), '-out', out, '-quiet'] + args
        elif kind == 'spscore':
            write_fasta(os.path.join(d, 'x.fa'), data.split(';'), 's')
            cmd = base + ['-spscore', os.path.join(d, 'x.fa'), '-quiet'] + args
        else:
            return [sid, kind, typ, opts, 'UNKNOWNKIND', '', '', '', '', '', '', '', '', '']
        try:
            p = subprocess.run(cmd, capture_output=True, cwd=d, timeout=600)
        except subprocess.TimeoutExpired:
            return [sid, kind, typ, opts, 'TIMEOUT', '', '', '', '', '', '', '', usedtree, '']
        msg = esc(dec(p.stdout + p.stderr).strip()[:400])
        names, rows, text = [], [], ''
        if os.path.exists(out) and os.path.getsize(out) > 0:
            if kind == 'maketree':
                text = open(out, errors='replace').read()
            else:
                try:
                    names, rows = read_fasta(out)
                except Exception as e:
                    msg += ' | read error ' + str(e)
        rc = str(p.returncode)
        if p.returncode != 0:
            rc = f'FAILED{p.returncode}'
        elif kind not in ('spscore', 'ppscore') and not rows and not text and 'computeweights' not in files and 'clusteronly' not in opts:
            rc = 'NOOUT'
        def f(k):
            return esc(open(files[k], errors='replace').read()) if k in files and os.path.exists(files[k]) else ''
        return [sid, kind, typ, opts, rc, ','.join(names), ';'.join(rows) if rows else esc(text), f('tree1'), f('tree2'), f('distmx1'), f('distmx2'), f('computeweights'), esc(usedtree), msg]


def main():
    os.makedirs('out', exist_ok=True)
    n = failed = 0
    with open(os.path.join(HERE, 'inputs', 'sets.tsv')) as f, open('out/muscle3_r3_oracle.tsv', 'w') as o:
        f.readline()
        o.write('id\tkind\ttype\toptions\trc\tnames\trows\ttree1\ttree2\tdist1\tdist2\tweights\tusedtree\tstdout\n')
        for line in f:
            sid, kind, typ, data, opts = line.rstrip('\n').split('\t')
            try:
                row = run_one(sid, kind, typ, data, opts)
            except Exception as e:
                row = [sid, kind, typ, opts, 'RUNNER_ERROR', '', '', '', '', '', '', '', '', esc(str(e)[:300])]
            n += 1
            failed += not row[4].isdigit()
            o.write('\t'.join(row) + '\n'); o.flush()
    print(f'{n} sets run; non-zero or no output: {failed}')


if __name__ == '__main__':
    main()
