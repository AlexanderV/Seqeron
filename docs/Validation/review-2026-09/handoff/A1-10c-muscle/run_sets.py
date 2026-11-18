#!/usr/bin/env python3
"""Runs MUSCLE 3.8.1551 (/opt/muscle3/muscle3) on every set of inputs/sets.tsv and writes out/muscle3_oracle.tsv.

Output: id  type  options  order(ids of the rows in the order MUSCLE wrote them, comma separated)  rows(;)  (rows in INPUT order, aligned, upper case as MUSCLE writes them)
A set on which MUSCLE fails (non-zero exit / no output) is written with rows = FAILED and the message; the count is printed.
"""
import os, subprocess, tempfile

BIN = os.environ.get('MUSCLE3', '/opt/muscle3/muscle3')
FLAG = {'maxiters': '-maxiters', 'maxtrees': '-maxtrees', 'gapopen': '-gapopen', 'gapextend': '-gapextend', 'center': '-center', 'distance1': '-distance1',
        'cluster1': '-cluster1', 'cluster2': '-cluster2', 'sueff': '-sueff'}


def run(typ, seqs, opts):
    with tempfile.TemporaryDirectory() as d:
        fa, out = os.path.join(d, 'x.fa'), os.path.join(d, 'o.fa')
        with open(fa, 'w') as f:
            for i, s in enumerate(seqs):
                f.write(f'>s{i}\n{s}\n')
        args = [BIN, '-in', fa, '-out', out, '-quiet', '-seqtype', typ]
        for kv in opts.split(',') if opts else []:
            k, _, v = kv.partition('=')
            if k == 'model':
                args.append('-' + v)
            else:
                args += [FLAG[k], v]
        p = subprocess.run(args, capture_output=True, text=True, cwd=d)
        if p.returncode != 0 or not os.path.exists(out):
            return None, None, (p.stdout + p.stderr).strip().replace('\n', ' | ')[:300]
        rows, order, name = {}, [], None
        for line in open(out):
            line = line.rstrip('\n')
            if line.startswith('>'):
                name = line[1:].strip()
                rows[name] = ''
                order.append(int(name[1:]))
            elif name:
                rows[name] += line.strip()
        return [rows[f's{i}'] for i in range(len(seqs))], order, ''


def main():
    os.makedirs('out', exist_ok=True)
    failed = 0
    with open('inputs/sets.tsv') as f, open('out/muscle3_oracle.tsv', 'w') as o:
        f.readline()
        o.write('id\ttype\toptions\torder\trows\n')
        for line in f:
            sid, typ, seqs, opts = line.rstrip('\n').split('\t')
            rows, order, msg = run(typ, seqs.split(';'), opts)
            if rows is None:
                failed += 1
                o.write(f'{sid}\t{typ}\t{opts}\t\tFAILED {msg}\n')
            else:
                o.write(f'{sid}\t{typ}\t{opts}\t{",".join(map(str, order))}\t{";".join(rows)}\n')
    print('sets run; failed:', failed)


if __name__ == '__main__':
    main()
