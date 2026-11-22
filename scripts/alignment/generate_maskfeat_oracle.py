#!/usr/bin/env python3
"""A5-8c3g (B09 F249) oracle: EMBOSS 6.6.0 `maskfeat -auto -sequence in.dat [-sformat1 F] -type TYPES [-tolower | -maskchar C] [-supper1 | -slower1] -outseq o.fa -osformat2 fasta`: the sequence with the regions of the
features of the wanted types masked (a character, or lower case). Inputs: entries with feature tables in the formats whose readers carry the features (embl, genbank, swiss, nbrf, gff3, refseqp, gff2) made as in
generate_seqret_oracle.py, FASTA files with no features, damaged files; types: single types, wildcards (repeat*, *UTR), lists separated by blanks / commas / semicolons / bars, the default.
JSONL.GZ rows: id, fmt, text, spec, args, rc, fa, err
Usage: python3 scripts/alignment/generate_maskfeat_oracle.py OUT.jsonl.gz [n=1000] [seed]      (needs EMBOSS 6.6.0)
"""
import gzip, json, os, random, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_seqret_oracle as S  # noqa: E402
import generate_diffseq_features_oracle as D2  # noqa: E402

FMTS = ['embl', 'genbank', 'swiss', 'nbrf', 'gff3', 'refseqp', 'gff2', 'fasta', 'embl', 'genbank', 'gff3', 'swiss']
TYPES = ['repeat*', 'gene', 'CDS', 'mRNA', 'exon', 'intron', 'misc_feature', 'source', '*UTR', 'repeat_region', 'tRNA rRNA', 'gene,CDS', 'exon;intron', 'gene|CDS', 'misc*', 'g?ne', '*', 'binding_site', 'variation',
         'polypeptide_region', 'polypeptide_domain', 'domain', 'region', 'Region', 'CDS gene mRNA', 'x', 'repeat_unit', 'databank_entry', '*_region', 'sequence_feature', '']


def main(out, n, seed):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        work = os.path.join(d, 'w'); os.mkdir(work)
        gen = os.path.join(d, 'g'); os.mkdir(gen)
        while len(rows) < n:
            fmt = FMTS[len(rows) % len(FMTS)]
            text, how = S.make_input(rnd, gen, fmt)
            if rnd.random() < 0.15: text, how = S.damaged(rnd, text), how + '+damaged'
            spec = 'in.dat'
            args = ['-sformat1', fmt] if rnd.random() < 0.6 else []
            if rnd.random() < 0.85:
                t = rnd.choice(TYPES)
                if t == '' and rnd.random() < 0.5: t = rnd.choice(TYPES[:-1])
                args += ['-type', t]
            if rnd.random() < 0.25: args.append('-tolower')
            elif rnd.random() < 0.5: args += ['-maskchar', rnd.choice(['N', 'X', 'n', 'x', '-', '*', '#', ' ', '.', '?', 'Q', ''])]
            if rnd.random() < 0.1: args.append(rnd.choice(['-supper1', '-slower1']))
            for nme in os.listdir(work): os.remove(os.path.join(work, nme))
            open(os.path.join(work, 'in.dat'), 'w', newline='\n', encoding='latin-1').write(text)
            hang = False
            try:
                import resource
                def limit():
                    resource.setrlimit(resource.RLIMIT_FSIZE, (4 << 20, 4 << 20))
                p = subprocess.run(['maskfeat', '-auto', '-sequence', spec] + args + ['-outseq', 'o.fa', '-osformat2', 'fasta'], cwd=work, capture_output=True, timeout=5, preexec_fn=limit)
                rc, err = p.returncode, p.stderr.decode('latin-1')
                if rc == -25: hang = True
            except subprocess.TimeoutExpired:
                rc, err, hang = -9, '', True
                subprocess.run(['pkill', '-9', 'maskfeat'])
            q = os.path.join(work, 'o.fa')
            fa = open(q, encoding='latin-1', newline='').read() if os.path.exists(q) and not hang else ''
            if err.startswith('Mask off features in a sequence\n'): err = err[len('Mask off features in a sequence\n'):]
            rows.append({'id': len(rows), 'fmt': fmt, 'how': how, 'text': text, 'spec': spec, 'args': args, 'rc': rc, 'fa': fa, 'err': err, 'hang': hang})
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    from collections import Counter
    print(len(rows), 'rows;', sorted(Counter((r['fmt'], 'hang' if r['hang'] else 'fa' if r['fa'] else 'none') for r in rows).items()))


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 1000, int(sys.argv[3]) if len(sys.argv) > 3 else 20261625)
