#!/usr/bin/env python3
"""A5-8c3h (B09 F263, WP58) oracle: EMBOSS 6.6.0 `twofeat` (emboss/twofeat.c: pairs of features of a sequence in a required relation, reported)
    twofeat -auto -sequence in.dat -outfile o.rep [-sformat1 F] [-sbegin1 N] [-send1 N] [-asource S] [-atype T] [-asense 0|+|-] [-aminscore X] [-amaxscore X] [-atag T] [-avalue V] [-bsource ..] [-btype ..] [-bsense ..] ...
             [-overlap A|O|NO|NW|AW|BW] [-minrange N] [-maxrange N] [-rangetype N|L|R|F] [-sense A|S|O] [-order A|AB|BA] [-twoout] [-typeout T] [-rformat2 FORMAT] [-raccshow2] [-rdesshow2] ... [-rmaxall2 N] [-rmaxseq2 N]
Inputs: EMBL entries with CDS / mRNA / exon / gene / misc features (generate_coderet_oracle.py: joins, complements, orders, partial ends, tags, several entries per file), the entries with feature tables of the reader oracle (generate_seqret_oracle.py),
damaged files; every report format of the library. JSONL.GZ rows: id, fmt, how, text, args, cmd, epoch, rc, out, err, hang, memerr
Usage: python3 scripts/alignment/generate_twofeat_oracle.py OUT.jsonl.gz [n=1000] [seed]      (needs EMBOSS 6.6.0, valgrind, gcc)
"""
import gzip, json, os, random, re, shutil, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_seqret_oracle as S  # noqa: E402
import generate_coderet_oracle as C  # noqa: E402
import generate_featreport_oracle as F  # noqa: E402

FMTS = ['own', 'own', 'own', 'own', 'own', 'embl', 'genbank', 'gff3', 'gff2', 'swiss', 'nbrf', 'refseqp', 'own', 'own']
TYPES = ['CDS', 'mRNA', 'exon', 'gene', 'misc_feature', 'repeat_region', 'intron', 'tRNA', '*', 'CDS exon', 'mRNA|gene', 'misc*', 'source', 'variation', 'binding_site', 'domain', 'region', 'polypeptide_region', '*UTR', 'x']
RF = ['table', 'table', 'table', 'simple', 'gff', 'gff3', 'excel', 'embl', 'genbank', 'dasgff', 'debug', 'srs', 'listfile', 'regions', 'feattable', 'seqtable', 'nametable', 'tagseq', 'motif', 'summary', 'draw', 'restrict', 'dbmotif', 'diffseq', 'pir', 'swiss']


def main(out, n, seed):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        work = os.path.join(d, 'w'); os.mkdir(work)
        gen = os.path.join(d, 'g'); os.mkdir(gen)
        shim = os.path.join(d, 'faketime.so')
        subprocess.run(['gcc', '-shared', '-fPIC', '-o', shim, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'faketime.c')], check=True)
        while len(rows) < n:
            fmt = FMTS[len(rows) % len(FMTS)]
            if fmt == 'own':
                text, how = ''.join(C.embl_entry(rnd, rnd.choice(['X%d' % rnd.randint(1, 99), 'AB%06d' % rnd.randint(1, 999999), 'Seq_%d' % rnd.randint(1, 9)])) for _ in range(rnd.randint(1, 3))), 'own'
                fmt = 'embl'
            else:
                text, how = S.make_input(rnd, gen, fmt)
            if rnd.random() < 0.06: text, how = S.damaged(rnd, text), how + '+damaged'
            args = ['-auto', '-sequence', 'in.dat', '-outfile', 'o.rep']
            if rnd.random() < 0.1: args += ['-sformat1', fmt]
            if rnd.random() < 0.06: args += ['-sbegin1', str(rnd.randint(1, 20))]
            if rnd.random() < 0.06: args += ['-send1', str(rnd.randint(25, 80))]
            if rnd.random() < 0.03: args += ['-sreverse1']
            if rnd.random() < 0.65: args += ['-atype', rnd.choice(TYPES)]
            if rnd.random() < 0.65: args += ['-btype', rnd.choice(TYPES)]
            if rnd.random() < 0.07: args += ['-asource', rnd.choice(['*', 'EMBL', 'embl', 'x', 'GenBank'])]
            if rnd.random() < 0.07: args += ['-bsource', rnd.choice(['*', 'EMBL', 'embl', 'x', 'GenBank'])]
            if rnd.random() < 0.15: args += ['-asense', rnd.choice(['0', '+', '-'])]
            if rnd.random() < 0.15: args += ['-bsense', rnd.choice(['0', '+', '-'])]
            if rnd.random() < 0.05: args += ['-aminscore', rnd.choice(['0', '1', '1.5', '-1', '100'])]
            if rnd.random() < 0.05: args += ['-amaxscore', rnd.choice(['0', '1', '2.7', '100', '-1'])]
            if rnd.random() < 0.05: args += ['-bminscore', rnd.choice(['0', '1', '1.5'])]
            if rnd.random() < 0.05: args += ['-bmaxscore', rnd.choice(['0', '1', '2.7', '100'])]
            if rnd.random() < 0.08: args += ['-atag', rnd.choice(['*', 'gene', 'note', 'translation', 'x', 'gene|note'])]
            if rnd.random() < 0.08: args += ['-btag', rnd.choice(['*', 'gene', 'note', 'translation', 'x'])]
            if rnd.random() < 0.05: args += ['-avalue', rnd.choice(['*', 'g1', 'a*', 'x'])]
            if rnd.random() < 0.05: args += ['-bvalue', rnd.choice(['*', 'g1', 'a*', 'x'])]
            if rnd.random() < 0.5: args += ['-overlap', rnd.choice(['A', 'O', 'NO', 'NW', 'AW', 'BW'])]
            if rnd.random() < 0.3: args += ['-minrange', str(rnd.choice([0, 1, 5, 20, -3, 50]))]
            if rnd.random() < 0.3: args += ['-maxrange', str(rnd.choice([0, 10, 30, 100, 500, -2]))]
            if rnd.random() < 0.3: args += ['-rangetype', rnd.choice(['N', 'L', 'R', 'F'])]
            if rnd.random() < 0.25: args += ['-sense', rnd.choice(['A', 'S', 'O'])]
            if rnd.random() < 0.1: args += ['-order', rnd.choice(['A', 'AB', 'BA'])]
            two = rnd.random() < 0.2
            if two: args.append('-twoout')
            elif rnd.random() < 0.1: args += ['-typeout', rnd.choice(['gene', 'CDS', 'misc_feature', 'repeat_region', 'x', 'exon'])]
            if rnd.random() < 0.9: args += ['-rformat2', rnd.choice(RF)]
            if rnd.random() < 0.1: args.append('-raccshow2')
            if rnd.random() < 0.1: args.append('-rdesshow2')
            if rnd.random() < 0.06: args += ['-rscoreshow2', 'N']
            if rnd.random() < 0.06: args += ['-rstrandshow2', 'N']
            if rnd.random() < 0.08: args.append('-rusashow2')
            if rnd.random() < 0.08: args += ['-rmaxall2', str(rnd.randint(1, 4))]
            if rnd.random() < 0.08: args += ['-rmaxseq2', str(rnd.randint(1, 4))]
            args = F.report_file_args(rnd, args, 2)
            F.clean(work)
            open(os.path.join(work, 'in.dat'), 'w', newline='\n', encoding='latin-1').write(text)
            os.mkdir(os.path.join(work, 'sub'))
            epoch = rnd.randint(978307200, 2051222400)
            rc, err, hang = F.run('twofeat', work, args, shim, epoch)
            outfile, outtext = F.report_output(work, hang)
            if err.startswith('Find neighbouring pairs of features in sequence(s)\n'): err = err[len('Find neighbouring pairs of features in sequence(s)\n'):]
            memerr = (not hang) and (rc != 0 or len(rows) % 5 == 0) and F.valgrind('twofeat', work, args)
            rows.append({'id': len(rows), 'fmt': fmt, 'how': how, 'text': text, 'args': args, 'cmd': F.entries_of(args), 'epoch': epoch, 'rc': rc, 'out': outtext, 'outfile': outfile, 'err': err, 'hang': hang, 'memerr': memerr})
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    from collections import Counter
    print(len(rows), 'rows;', sorted(Counter((r['fmt'], 'hang' if r['hang'] else 'died' if r['rc'] else 'out' if r['out'] else 'none') for r in rows).items()))


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 1000, int(sys.argv[3]) if len(sys.argv) > 3 else 20265809)
