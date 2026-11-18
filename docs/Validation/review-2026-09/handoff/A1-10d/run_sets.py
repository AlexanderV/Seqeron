#!/usr/bin/env python3
"""Runs the patched clustalw 2.1 on every set of inputs/sets.tsv in four modes and writes out/end_gaps_oracle.tsv.

modes (environment variables read by the patched `-ENDGAPS` handler, see Dockerfile): control (none), useendgaps (CW_USEENDGAPS=1),
endgappen (CW_ENDGAPPEN=1), both. Every run passes -ENDGAPS so that the patched handler executes.
Output: id  type  mode  rows(;)     + out/control_check.txt (control rows equal the committed binary rows? count)
"""
import os, subprocess, sys, tempfile

FLAGS = {'gapopen': 'GAPOPEN', 'gapext': 'GAPEXT', 'transweight': 'TRANSWEIGHT', 'maxdiv': 'MAXDIV', 'noweights': 'NOWEIGHTS', 'negative': 'NEGATIVE',
         'quick': 'QUICKTREE', 'ktuple': 'KTUPLE', 'window': 'WINDOW', 'topdiags': 'TOPDIAGS', 'pairgap': 'PAIRGAP', 'pwgapopen': 'PWGAPOPEN',
         'pwgapext': 'PWGAPEXT', 'matrix': 'MATRIX', 'pwmatrix': 'PWMATRIX', 'nopgap': 'NOPGAP', 'nohgap': 'NOHGAP', 'novgap': 'NOVGAP',
         'gapdist': 'GAPDIST', 'hgapresidues': 'HGAPRESIDUES', 'endgaps': 'ENDGAPS'}
MODES = {'control': {}, 'useendgaps': {'CW_USEENDGAPS': '1'}, 'endgappen': {'CW_ENDGAPPEN': '1'}, 'both': {'CW_USEENDGAPS': '1', 'CW_ENDGAPPEN': '1'}}
BIN = os.environ.get('CLUSTALW', '/opt/clustalw/clustalw')


def run(typ, seqs, opts, env):
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, 'x.fa')
        with open(fa, 'w') as f:
            for i, s in enumerate(seqs):
                f.write(f'>s{i}\n{s}\n')
        out = os.path.join(d, 'o.fa')
        args = [BIN, f'-INFILE={fa}', f'-TYPE={typ}', '-ALIGN', '-OUTPUT=FASTA', '-OUTORDER=INPUT', f'-OUTFILE={out}', '-ENDGAPS']
        for kv in opts.split(',') if opts else []:
            k, _, v = kv.partition('=')
            if k == 'endgaps':
                continue
            args.append(f'-{FLAGS[k]}' + (f'={v}' if v else ''))
        subprocess.run(args, capture_output=True, text=True, cwd=d, env=dict(os.environ, **env))
        rows, name = {}, None
        for line in open(out):
            line = line.rstrip('\n')
            if line.startswith('>'):
                name = line[1:].strip(); rows[name] = ''
            elif name:
                rows[name] += line.strip()
    return [rows[f's{i}'] for i in range(len(seqs))]


def main():
    os.makedirs('out', exist_ok=True)
    same = total = 0
    with open('inputs/sets.tsv') as f, open('out/end_gaps_oracle.tsv', 'w') as o:
        f.readline()
        o.write('id\ttype\tmode\trows\n')
        for line in f:
            sid, typ, seqs, opts, control = line.rstrip('\n').split('\t')
            seqs = seqs.split(';')
            for mode, env in MODES.items():
                rows = run(typ, seqs, opts, env)
                o.write(f'{sid}\t{typ}\t{mode}\t{";".join(rows)}\n')
                if mode == 'control':
                    total += 1
                    same += rows == control.split(';')
    with open('out/control_check.txt', 'w') as f:
        f.write(f'control rows equal the committed binary rows: {same} of {total}\n')
    print('control equal', same, 'of', total)


if __name__ == '__main__':
    main()
