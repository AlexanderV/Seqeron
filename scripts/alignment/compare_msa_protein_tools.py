#!/usr/bin/env python3
"""INFORMATIONAL comparison (audit A1-14, B09 F61): sum-of-pairs score of protein alignments made by the locally installed MSA programs
(muscle, mafft, kalign, clustalo, clustalw, t_coffee) and by the Seqeron methods, rescored under ONE objective (BLOSUM62, natural affine gaps with the
EMBOSS needle numbers gapopen 10 / gapextend 0.5 = GapOpen -9.5, GapExtend -0.5; Altschul 1989).

The Seqeron rows are produced by the independent Python implementations of generate_msa_matrix_oracle.py (they equal the C# rows on all 238 oracle sets,
so these are the Seqeron numbers). The other programs optimise their own objectives (profile scores, consistency, HMMs), so a lower SP score is not an
error - this table only shows where the repository's methods stand. Nothing is asserted.
Usage: python3 scripts/alignment/compare_msa_protein_tools.py [SETS=20]
"""
import importlib.util, os, random, subprocess, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('mo', os.path.join(HERE, 'generate_msa_matrix_oracle.py'))
mo = importlib.util.module_from_spec(spec); spec.loader.exec_module(mo)

SC = mo.Sc('BLOSUM62', -9.5, -0.5)


def fasta(path, seqs):
    open(path, 'w').write(''.join(f'>s{i}\n{s}\n' for i, s in enumerate(seqs)))


def read_fasta(path):
    rows, cur = {}, None
    for l in open(path):
        l = l.strip()
        if l.startswith('>'): cur = l[1:].split()[0]; rows[cur] = ''
        elif cur is not None: rows[cur] += l
    return rows


def run(tool, seqs, d):
    fa, out = os.path.join(d, 'in.fa'), os.path.join(d, 'out.fa')
    fasta(fa, seqs)
    cmds = {
        'muscle': ['muscle', '-align', fa, '-output', out],
        'mafft': ['sh', '-c', f'mafft --quiet --auto {fa} > {out}'],
        'kalign': ['kalign', '-i', fa, '-o', out, '-f', 'fasta'],
        'clustalo': ['clustalo', '-i', fa, '-o', out, '--force', '--outfmt=fa'],
        'clustalw': ['clustalw', '-TYPE=PROTEIN', '-ALIGN', '-OUTPUT=FASTA', f'-INFILE={fa}', f'-OUTFILE={out}', '-OUTORDER=INPUT'],
        't_coffee': ['t_coffee', '-seq', fa, '-output', 'fasta_aln', '-outfile', out, '-quiet', '-n_core', '1'],
    }
    try:
        subprocess.run(cmds[tool], capture_output=True, cwd=d, timeout=300, env=dict(os.environ, HOME=d))
        rows = read_fasta(out)
        res = [rows[f's{i}'].upper().replace('.', '-') for i in range(len(seqs))]
        if [r.replace('-', '') for r in res] != seqs: return None
        return res
    except Exception:
        return None


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    rnd = random.Random(20261124)
    totals, wins = {}, {}
    sets = []
    for t in range(n):
        base = mo.rseq(rnd, rnd.randint(45, 70), mo.AA)
        k = rnd.randint(4, 6)
        sets.append([mo.mutate(rnd, base, rnd.choice([0.1, 0.25, 0.4]), mo.AA) for _ in range(k)])
    for seqs in sets:
        scores = {}
        with tempfile.TemporaryDirectory() as d:
            for tool in ('muscle', 'mafft', 'kalign', 'clustalo', 'clustalw', 't_coffee'):
                rows = run(tool, seqs, d)
                if rows: scores[tool] = mo.sp_affine(rows, SC)
        st = mo.star(seqs, SC, 'affine'); pr, it = mo.progressive_iterative(seqs, SC, 'affine'); co = mo.consistency(seqs, SC, 'affine')
        for name, rows in (('seqeron-star', st), ('seqeron-progressive', pr), ('seqeron-iterative', it), ('seqeron-consistency', co)):
            scores[name] = mo.sp_affine(rows, SC)
        scores['seqeron-exact-bound(sum pair optima)'] = sum(mo.bio_optimum(seqs[i], seqs[j], SC, 'affine') for i in range(len(seqs)) for j in range(i + 1, len(seqs)))
        best_tool = max(v for k_, v in scores.items() if not k_.startswith('seqeron'))
        for k_, v in scores.items():
            totals.setdefault(k_, []).append(mo.Fraction(v, SC.scale))
            wins[k_] = wins.get(k_, 0) + (1 if v >= best_tool else 0)
    print(f'{n} random protein families (k = 4..6, length 45-70, divergence 10-40 %); SP score under BLOSUM62 / gapopen 10 / gapextend 0.5 (natural affine)')
    for k_, v in sorted(totals.items(), key=lambda kv: -sum(kv[1])):
        print(f'{k_:42s} mean SP {float(sum(v)) / len(v):9.1f}   sets >= best external program: {wins[k_]:3d}/{len(v)}   (sets scored: {len(v)})')


if __name__ == '__main__':
    main()
