"""Oracle sections of the Bio.Align readers msf, emboss, hhr and tabular (B09 A2-27d): the input texts come from the installed EMBOSS 6.6.0 (needle / water / stretcher / matcher / seqret),
HH-suite 3.3.0 (hhalign), BLAST+ and fasta36 where a program writes the format, plus hand-made variants and refusals; the expected value is what `Bio.Align.parse` gives."""
import os
import random
import re
import shutil
import subprocess
import tempfile

AA = 'ACDEFGHIKLMNPQRSTVWY'
DNA = 'ACGT'


def run(cmd, cwd, timeout=60):
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout)
    return r.stdout


def fixdate(text):
    return re.sub(r'(# Rundate: ).*', r'\1Sat 21 Nov 2026 19:00:25', text)


def rseq(rng, alphabet, n):
    return ''.join(rng.choice(alphabet) for _ in range(n))


def related(rng, alphabet, base, p):
    out = []
    for c in base:
        r = rng.random()
        if r < p / 2:
            continue
        if r < p:
            out.append(rng.choice(alphabet))
            continue
        out.append(c)
        if rng.random() < p / 4:
            out.append(rng.choice(alphabet))
    return ''.join(out)


# ---------------------------------------------------------------------------- msf

def msf_checksum_line(names, seqs, header_kind, rng, gap='.'):
    pass


def gcg_msf(rng, names, seqs, kind='N', dot='.', lead='  stdout', check_word='CompCheck:', block=50, group=10, coords=True, pileup=False, weights=None, lenmode='residues'):
    """A GCG MSF text built like EMBOSS / MUSCLE output."""
    n = len(seqs)
    cols = len(seqs[0])
    out = []
    out.append('PileUp' if pileup else '!!%s_MULTIPLE_ALIGNMENT 1.0' % ('NA' if kind == 'N' else 'AA'))
    out.append('')
    out.append('%s MSF: %d Type: %s 21/11/26 %s %d ..' % (lead, cols, kind, check_word, rng.randrange(1000, 9999)))
    out.append('')
    for i in range(n):
        if lenmode == 'residues':
            ln = len(seqs[i].replace('-', ''))
        else:
            ln = len(seqs[i])
        w = weights[i] if weights else 1.0
        out.append('  Name: %-10s Len: %d  Check: %d Weight: %.2f' % (names[i], ln, rng.randrange(0, 9999), w))
    out.append('')
    out.append('//')
    out.append('')
    for s0 in range(0, cols, block):
        if coords:
            out.append('%12d %*d' % (s0 + 1, block + (block // group) - 1 + 1, min(s0 + block, cols)) if False else '%11d%s%d' % (s0 + 1, ' ' * 3, min(s0 + block, cols)))
        for i in range(n):
            chunk = seqs[i][s0:s0 + block]
            groups = [chunk[g:g + group] for g in range(0, len(chunk), group)]
            out.append('%-10s %s' % (names[i], ' '.join(g.replace('-', dot) for g in groups)))
        out.append('')
    return '\n'.join(out) + '\n'


def random_msa_rows(rng, alphabet, n, cols):
    base = rseq(rng, alphabet, cols)
    rows = []
    for _ in range(n):
        r = ''.join(c if rng.random() > 0.2 else '-' for c in base)
        rows.append(r)
    return rows


def run_msf(ctx):
    rng = ctx.rng
    add_read = ctx.add_read
    tmp = tempfile.mkdtemp()
    try:
        for i in range(70):
            alphabet = rng.choice([DNA, AA])
            n = rng.randrange(2, 6)
            cols = rng.randrange(8, 130)
            rows = random_msa_rows(rng, alphabet, n, cols)
            if rng.random() < 0.3:
                rows = [r.lower() if rng.random() < 0.3 else r for r in rows]
            names = ['seq%d' % (k + 1) if rng.random() < 0.5 else rng.choice(['sp|P%05d' % rng.randrange(99999), 'Q%d_HUMAN' % rng.randrange(9999), 'gi%d' % rng.randrange(99999)]) + '_%d' % k for k in range(n)]
            # real EMBOSS output
            fa = ''.join('>%s\n%s\n' % (names[k], rows[k]) for k in range(n))
            open(os.path.join(tmp, 'a.fa'), 'w').write(fa)
            out = run(['seqret', '-sequence', 'a.fa', '-osformat', 'msf', '-auto', '-stdout'], tmp)
            add_read('msf', out)
            if i % 3 == 0:
                add_read('msf', out.replace('CompCheck:', 'Check:'))
        for i in range(90):
            alphabet = rng.choice([DNA, AA])
            n = rng.randrange(1, 6)
            cols = rng.randrange(5, 140)
            rows = random_msa_rows(rng, alphabet, n, cols)
            if n == 1:
                rows = [rows[0].replace('-', '')]
                cols = len(rows[0])
            names = ['s%d' % (k + 1) for k in range(n)]
            kind = 'N' if alphabet == DNA else 'P'
            text = gcg_msf(rng, names, rows, kind=kind, dot=rng.choice(['.', '-', '~']), lead=rng.choice(['  stdout', '', 'x.msf', '  ']), check_word=rng.choice(['Check:', 'CompCheck:']),
                           block=rng.choice([50, 20, 60]), group=rng.choice([10, 5]), coords=rng.random() < 0.7, pileup=rng.random() < 0.4,
                           weights=[round(rng.uniform(0.1, 2), 2) for _ in range(n)] if rng.random() < 0.5 else None, lenmode=rng.choice(['residues', 'residues', 'columns']))
            add_read('msf', text)
        for i in range(12):
            alphabet = rng.choice([DNA, AA])
            rows = random_msa_rows(rng, alphabet, 3, 40)
            names = ['a', 'b', 'c']
            base = gcg_msf(rng, names, rows, kind='N' if alphabet == DNA else 'P', dot='-', lead='x', check_word='Check:')
            # damaged variants
            add_read('msf', base.replace('//', ''))
            add_read('msf', base.split('//')[0] + '//\n')
            add_read('msf', base.replace('Name: b ', 'Name: a '))
            add_read('msf', base.replace('Type: N', 'Type: X').replace('Type: P', 'Type: X'))
            add_read('msf', base.replace('!!NA_MULTIPLE_ALIGNMENT 1.0', 'Header').replace('!!AA_MULTIPLE_ALIGNMENT 1.0', 'Header'))
            add_read('msf', base + 'trailing text\n')
            add_read('msf', re.sub(r'MSF: \d+', 'MSF: 7', base))
            add_read('msf', re.sub(r'MSF: \d+', 'MSF: abc', base))
            add_read('msf', base.replace('\n//\n\n', '\n//\n'))
            add_read('msf', base.rstrip().rsplit('\n', 3)[0] + '\n')
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    add_read('msf', '')
    add_read('msf', '\n')


# ---------------------------------------------------------------------------- emboss

def run_emboss(ctx):
    rng = ctx.rng
    add_read = ctx.add_read
    tmp = tempfile.mkdtemp()
    try:
        count = 0
        for i in range(150):
            prot = rng.random() < 0.5
            alphabet = AA if prot else DNA
            base = rseq(rng, alphabet, rng.randrange(25, 120))
            a = related(rng, alphabet, base, rng.choice([0.05, 0.2, 0.4]))
            b = related(rng, alphabet, base[rng.randrange(0, 10):], rng.choice([0.05, 0.2, 0.4]))
            nseq = rng.choice([1, 1, 1, 2])
            fa = ''.join('>a%d first%s\n%s\n' % (k, ' x' if k else '', a if k == 0 else related(rng, alphabet, base, 0.3)) for k in range(nseq))
            fb = ''.join('>b%d\n%s\n' % (k, b if k == 0 else related(rng, alphabet, base, 0.3)) for k in range(nseq))
            open(os.path.join(tmp, 'a.fa'), 'w').write(fa)
            open(os.path.join(tmp, 'b.fa'), 'w').write(fb)
            prog = rng.choice(['needle', 'water', 'stretcher', 'needle', 'water', 'matcher'])
            fmt = rng.choice(['pair', 'srspair', 'simple', 'pair', 'srspair'])
            cmd = [prog, '-asequence', 'a.fa', '-bsequence', 'b.fa', '-aformat3', fmt, '-auto', '-stdout']
            if prog in ('needle', 'water'):
                cmd += ['-gapopen', str(rng.choice([8, 10, 12])), '-gapextend', str(rng.choice([0.5, 1, 2]))]
            if prog == 'stretcher':
                cmd += ['-gapopen', str(rng.choice([8, 12])), '-gapextend', str(rng.choice([2, 4]))]
            if prog == 'matcher':
                cmd += ['-alternatives', str(rng.choice([1, 2, 3]))]
            if not prot and rng.random() < 0.3:
                cmd += ['-sreverse2']
            if rng.random() < 0.25 and prog in ('needle', 'water'):
                cmd += ['-nobrief']
            if rng.random() < 0.2:
                cmd += ['-awidth3', str(rng.choice([30, 50, 80, 100]))]
            try:
                out = run(cmd, tmp)
            except Exception:  # noqa: BLE001
                continue
            if not out.strip():
                continue
            add_read('emboss', fixdate(out))
            count += 1
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    # variants of the last texts
    ex = ctx.last_rows('emboss', 12)
    for text in ex:
        add_read('emboss', text.replace('# Aligned_sequences: 2', '# Aligned_sequences: 3'))
        add_read('emboss', text.replace('# Length:', '# Zength:'))
        add_read('emboss', text[:len(text) // 2])
        add_read('emboss', text.replace('# Score:', '# Score: x'))
    add_read('emboss', '')
    add_read('emboss', 'not emboss\n')
    add_read('emboss', '########################################\n')
    add_read('emboss', '########################################\n# Program: needle\n########################################\n')


# ---------------------------------------------------------------------------- hhr

def write_a3m(path, entries):
    with open(path, 'w') as fh:
        for name, seq in entries:
            fh.write('>%s\n%s\n' % (name, seq))


def run_hhr(ctx):
    rng = ctx.rng
    add_read = ctx.add_read
    tmp = tempfile.mkdtemp()
    try:
        for i in range(110):
            n_q = rng.randrange(1, 5)
            base = rseq(rng, AA, rng.randrange(40, 160))
            q = [('query%d %s' % (i, rng.choice(['desc one', 'protein x', 'a'])), base)]
            for k in range(n_q - 1):
                q.append(('q%d' % k, related(rng, AA, base, 0.25)))
            write_a3m(os.path.join(tmp, 'q.a3m'), q)
            t_entries = []
            for k in range(rng.randrange(1, 4)):
                lo = rng.randrange(0, max(1, len(base) // 3))
                hi = rng.randrange(len(base) * 2 // 3, len(base) + 1)
                t_entries.append(('target%d %s' % (k, rng.choice(['some description', 'hypothetical protein', 'x y z'])), related(rng, AA, base[lo:hi], rng.choice([0.1, 0.3, 0.5]))))
            write_a3m(os.path.join(tmp, 't.a3m'), t_entries[:1])
            cmd = ['hhalign', '-i', 'q.a3m', '-t', 't.a3m', '-o', 'out.hhr', '-v', '0']
            if rng.random() < 0.3:
                cmd += ['-glob']
            if rng.random() < 0.3:
                cmd += ['-ssm', '0']
            try:
                run(cmd, tmp)
                text = open(os.path.join(tmp, 'out.hhr')).read()
            except Exception:  # noqa: BLE001
                continue
            text = re.sub(r'(Date\s+).*', r'\1Sat Nov 21 18:57:54 2026', text)
            text = re.sub(r'(Command\s+).*', r'\1hhalign -i q.a3m -t t.a3m -o out.hhr -v 0', text)
            add_read('hhr', text)
            if i % 6 == 0:
                add_read('hhr', text.replace('Done!', 'Done!\nextra'))
                add_read('hhr', text.replace('Done!\n', ''))
            if i % 9 == 0:
                add_read('hhr', text.replace('Match_columns', 'Match_colums'))
                add_read('hhr', text.replace('\nNo 1\n', '\nNo 2\n'))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    add_read('hhr', '')
    add_read('hhr', 'Query x\n\n')
    add_read('hhr', 'Query x\nMatch_columns 5\n\n No Hit\n')


# ---------------------------------------------------------------------------- tabular

BLAST_FIELD_SETS = [
    'std', 'std btop', 'qacc sacc qstart qend sstart send btop', 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore btop',
    'std qlen slen', 'std qseq sseq', 'qacc.ver sacc.ver qstart qend sstart send btop qlen', 'std score', 'std nident positive gaps ppos',
    'qseqid sseqid qstart qend sstart send btop staxids stitle', 'std qframe sframe', 'std frames', 'qacc sacc evalue bitscore',
]


def run_tabular(ctx):
    rng = ctx.rng
    add_read = ctx.add_read
    tmp = tempfile.mkdtemp()
    try:
        for i in range(140):
            prog = rng.choice(['blastn', 'blastn', 'tblastn', 'blastp', 'blastx', 'tblastx'])
            q_prot = prog in ('blastp', 'tblastn')
            db_prot = prog in ('blastp', 'blastx')
            base_dna = rseq(rng, DNA, rng.randrange(150, 400))
            prot_of = lambda s: ''.join(rng.choice(AA) for _ in range(len(s) // 3))  # noqa: E731
            if prog == 'blastn':
                db = [('subj%d' % k, related(rng, DNA, base_dna, 0.1) if k == 0 else rseq(rng, DNA, 300)) for k in range(3)]
                q = [('qry%d' % i, base_dna[rng.randrange(0, 40):])]
                dbtype = 'nucl'
            elif prog == 'blastp':
                bp = rseq(rng, AA, rng.randrange(60, 130))
                db = [('subj%d' % k, related(rng, AA, bp, 0.2) if k == 0 else rseq(rng, AA, 100)) for k in range(3)]
                q = [('qry%d desc here' % i, bp)]
                dbtype = 'prot'
            elif prog == 'tblastn':
                bp = rseq(rng, AA, rng.randrange(60, 120))
                codons = {a: ''.join(rng.choice(DNA) for _ in range(3)) for a in AA}
                # use the standard code through a translation of a random DNA instead: simply take protein from translated DNA
                from Bio.Seq import Seq
                dna = base_dna
                prot = str(Seq(dna[:len(dna) // 3 * 3]).translate())
                db = [('subj%d' % k, related(rng, DNA, dna, 0.05) if k == 0 else rseq(rng, DNA, 300)) for k in range(2)]
                q = [('qry%d' % i, prot[rng.randrange(0, 10):].replace('*', 'A'))]
                dbtype = 'nucl'
            elif prog == 'blastx':
                from Bio.Seq import Seq
                prot = str(Seq(base_dna[:len(base_dna) // 3 * 3]).translate()).replace('*', 'A')
                db = [('subj%d' % k, prot if k == 0 else rseq(rng, AA, 100)) for k in range(2)]
                q = [('qry%d' % i, base_dna)]
                dbtype = 'prot'
            else:
                db = [('subj%d' % k, related(rng, DNA, base_dna, 0.05) if k == 0 else rseq(rng, DNA, 300)) for k in range(2)]
                q = [('qry%d' % i, base_dna)]
                dbtype = 'nucl'
            write_a3m(os.path.join(tmp, 'db.fa'), db)
            write_a3m(os.path.join(tmp, 'q.fa'), q)
            run(['makeblastdb', '-in', 'db.fa', '-dbtype', dbtype, '-out', 'db', '-parse_seqids'], tmp)
            fields = rng.choice(BLAST_FIELD_SETS)
            cmd = [prog, '-query', 'q.fa', '-db', 'db', '-outfmt', '7 ' + fields, '-evalue', '10']
            if prog == 'blastn' and rng.random() < 0.5:
                cmd += ['-strand', rng.choice(['plus', 'minus', 'both'])]
            if prog in ('blastn',):
                cmd += ['-task', 'blastn']
            try:
                out = run(cmd, tmp)
            except Exception:  # noqa: BLE001
                continue
            out = re.sub(r'# Database: .*', '# Database: db', out)
            if out.strip():
                add_read('tabular', out)
        for i in range(60):
            # FASTA m8CB / m8CC
            dna = rng.random() < 0.5
            alphabet = DNA if dna else AA
            base = rseq(rng, alphabet, rng.randrange(60, 200))
            write_a3m(os.path.join(tmp, 'f1.fa'), [('fq%d desc' % i, base)])
            write_a3m(os.path.join(tmp, 'f2.fa'), [('fs%d' % k, related(rng, alphabet, base, 0.15) if k == 0 else rseq(rng, alphabet, 120)) for k in range(3)])
            prog = rng.choice(['ssearch36', 'fasta36', 'glsearch36'])
            cmd = [prog, '-m', rng.choice(['8CB', '8CC']), '-q', 'f1.fa', 'f2.fa']
            if dna:
                cmd.insert(1, '-n')
            try:
                out = run(cmd, tmp)
            except Exception:  # noqa: BLE001
                continue
            if out.strip():
                add_read('tabular', out)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    add_read('tabular', '')
    add_read('tabular', 'no header\n')
    add_read('tabular', '# BLASTN 2.12.0+\n# Query: q\n# Database: d\n# Fields: query id, bogus\n# 1 hits found\nq\tx\n# BLAST processed 1 queries\n')
