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


def run(typ, seqs, opts, env, tag='x'):
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
        p = subprocess.run(args, capture_output=True, text=True, cwd=d, env=dict(os.environ, **env))
        rows, name = {}, None
        raw = open(out).read() if os.path.exists(out) else None
        for line in (raw or '').splitlines():
            if line.startswith('>'):
                name = line[1:].strip(); rows[name] = ''
            elif name:
                rows[name] += line.strip()
        got = [rows.get(f's{i}') for i in range(len(seqs))]
        if any(g is None for g in got):
            # keep everything needed to diagnose (the self-built binary may differ from the installed one)
            dbg = os.path.join('out', 'debug', f'{tag}')
            os.makedirs(dbg, exist_ok=True)
            with open(os.path.join(dbg, 'info.txt'), 'w') as f:
                f.write('args: ' + ' '.join(args) + '\nenv: ' + repr(env) + f'\nrc: {p.returncode}\n--- stdout ---\n{p.stdout}\n--- stderr ---\n{p.stderr}\n--- output file ({"missing" if raw is None else "present"}) ---\n{raw}\n--- input ---\n' + open(fa).read())
            return None
    return got


def preflight():
    # Written first, so that even a total failure leaves a diagnosis: version, a trivial DNA run, a trivial protein run.
    os.makedirs('out', exist_ok=True)
    with open('out/preflight.txt', 'w') as f:
        for label, typ, seqs in (('DNA', 'DNA', ['ACGTACGTACGTAAAC', 'ACGTACGAACGTAAAC', 'ACGTTCGTACGTAAAC']),
                                 ('PROTEIN', 'PROTEIN', ['MKVLAAGIVGLLLAQ', 'MKVLAGIVGLLAQ', 'MKILAAGVVGLLLAQ'])):
            with tempfile.TemporaryDirectory() as d:
                fa = os.path.join(d, 'x.fa')
                open(fa, 'w').write(''.join(f'>s{i}\n{s}\n' for i, s in enumerate(seqs)))
                args = [BIN, f'-INFILE={fa}', f'-TYPE={typ}', '-ALIGN', '-OUTPUT=FASTA', '-OUTORDER=INPUT', f'-OUTFILE={d}/o.fa', '-ENDGAPS']
                p = subprocess.run(args, capture_output=True, text=True, cwd=d)
                f.write(f'=== {label}: {" ".join(args)}\nrc={p.returncode}\n--- stdout ---\n{p.stdout}\n--- stderr ---\n{p.stderr}\n--- files in dir: {sorted(os.listdir(d))}\n')
                for fn in sorted(os.listdir(d)):
                    if fn.endswith(('.fa', '.dnd')) and fn != 'x.fa':
                        f.write(f'--- {fn} ---\n' + open(os.path.join(d, fn)).read() + '\n')


def main():
    os.makedirs('out', exist_ok=True)
    same = total = failed = 0
    preflight()
    with open('inputs/sets.tsv') as f, open('out/end_gaps_oracle.tsv', 'w') as o:
        f.readline()
        o.write('id\ttype\tmode\trows\n')
        for line in f:
            sid, typ, seqs, opts, control = line.rstrip('\n').split('\t')
            seqs = seqs.split(';')
            for mode, env in MODES.items():
                rows = run(typ, seqs, opts, env, f'{sid}_{mode}')
                if rows is None:
                    failed += 1
                    o.write(f'{sid}\t{typ}\t{mode}\tFAILED\n')
                    if failed <= 3:
                        print('FAILED', sid, mode, '(details in out/debug/)', flush=True)
                    if mode == 'control':
                        total += 1
                    continue
                o.write(f'{sid}\t{typ}\t{mode}\t{";".join(rows)}\n')
                if mode == 'control':
                    total += 1
                    same += rows == control.split(';')
    with open('out/control_check.txt', 'w') as f:
        f.write(f'control rows equal the committed binary rows: {same} of {total}; failed runs: {failed}\n')
    print('control equal', same, 'of', total)


if __name__ == '__main__':
    main()
