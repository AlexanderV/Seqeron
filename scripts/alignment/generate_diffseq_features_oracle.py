#!/usr/bin/env python3
"""A5-8c3b (B09 F218) oracle: EMBOSS 6.6.0 `diffseq` on GenBank / EMBL entries WITH FEATURES (the `first_feature` / `second_feature` notes of the 7-tag report, the report in other formats, both GFF3 files), alone and with the
sequence qualifiers (-sbegin / -send / -sreverse: the feature table is reversed and trimmed with the sequence), mixed with FASTA inputs, -sformat.
Features: random types (EMBL keys, aliases, SO names, unknown keys, other case), locations (a..b, complement, join, order, complement(join), single positions, a^b, <a, >b, (a.b)..c, a..(b.c), remote entries), qualifiers (quoted text over
several lines, unquoted values, no value, db_xref, protein_id, translation, codon_start, unknown names, other case).
JSONL.GZ rows: id, nuc, file1, file2 (texts), format1, format2 (the names used in the file names a.gb / a.embl / a.fa), wordsize, glob, rq (the qualifier text), cmd, epoch, report, feat1, feat2, error, rc.
Usage: python3 scripts/alignment/generate_diffseq_features_oracle.py OUT.jsonl.gz [n=300] [seed]      (needs EMBOSS 6.6.0 and gcc)
"""
import gzip, json, os, random, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_diffseq_oracle as D  # noqa: E402

TYPES = ['gene', 'CDS', 'mRNA', 'exon', 'intron', 'misc_feature', 'tRNA', 'rRNA', 'repeat_region', 'promoter', 'misc_RNA', 'variation', 'STS', "5'UTR", "3'UTR", 'source', 'Gene', 'cds', 'binding_site', 'sequence_conflict', 'zzz_type', 'enhancer', 'terminator', 'misc_binding', 'protein_bind']
WORDS = ['alpha', 'beta', 'kinase', 'putative', 'hypothetical protein', 'ribosomal', 'subunit', 'transport', 'regulator', 'of', 'the', 'binding', 'domain', 'zinc', 'finger']


def wrap(text, width):
    """wrap on blanks at width characters (a word longer than width is cut)"""
    out, cur = [], ''
    for w in text.split(' '):
        while len(w) > width:
            if cur: out.append(cur); cur = ''
            out.append(w[:width]); w = w[width:]
        if not cur: cur = w
        elif len(cur) + 1 + len(w) <= width: cur += ' ' + w
        else: out.append(cur); cur = w
    if cur: out.append(cur)
    return out


def locations(rnd, n):
    def rng():
        a = rnd.randint(1, n); b = rnd.randint(a, min(n, a + rnd.randint(0, 60)))
        return a, b
    kind = rnd.random()
    if kind < 0.03:
        return '%d..%d' % (rnd.randint(1, n), n + rnd.randint(1, 20))
    if kind < 0.3:
        a, b = rng(); return '%d..%d' % (a, b)
    if kind < 0.45:
        a, b = rng(); return 'complement(%d..%d)' % (a, b)
    if kind < 0.65:
        parts = sorted(rng() for _ in range(rnd.randint(2, 4)))
        loc = ','.join('%d..%d' % p for p in parts)
        return 'join(%s)' % loc
    if kind < 0.75:
        parts = sorted(rng() for _ in range(rnd.randint(2, 3)))
        return 'complement(join(%s))' % ','.join('%d..%d' % p for p in parts)
    if kind < 0.8:
        parts = [rng() for _ in range(rnd.randint(2, 3))]
        return 'order(%s)' % ','.join('%d..%d' % p for p in parts)
    if kind < 0.86:
        return str(rnd.randint(1, n))
    if kind < 0.89:
        a = rnd.randint(1, n - 1); return '%d^%d' % (a, a + 1)
    if kind < 0.93:
        a, b = rng(); return rnd.choice(['<%d..%d', '%d..>%d', '<%d..>%d']) % (a, b)
    if kind < 0.96:
        a, b = rng(); return '(%d.%d)..%d' % (a, min(n, a + 3), b)
    parts = sorted(rng() for _ in range(2))
    return 'join(%d..%d,J00194.1:%d..%d)' % (parts[0][0], parts[0][1], parts[1][0], parts[1][1])


def qualifiers(rnd):
    quals = []
    for _ in range(rnd.choice([0, 1, 2, 3, 4, 6])):
        r = rnd.random()
        if r < 0.15: quals.append(('gene', '"%s"' % rnd.choice(['abc', 'recA', 'dnaK', 'x1'])))
        elif r < 0.25: quals.append(('product', '"%s"' % ' '.join(rnd.choice(WORDS) for _ in range(rnd.randint(1, 5)))))
        elif r < 0.4: quals.append(('note', '"%s"' % ' '.join(rnd.choice(WORDS) for _ in range(rnd.randint(1, 30)))))
        elif r < 0.47: quals.append(('db_xref', '"%s"' % rnd.choice(['taxon:9606', 'GI:12345', 'GeneID:77', 'badvalue'])))
        elif r < 0.52: quals.append(('protein_id', '"%s"' % rnd.choice(['AAA12345.1', 'CAB1234567.2', 'bad_id', 'ENSP00000123456'])))
        elif r < 0.57: quals.append(('codon_start', rnd.choice(['1', '2', '3'])))
        elif r < 0.6: quals.append(('transl_table', '11'))
        elif r < 0.66: quals.append(('translation', '"%s"' % ''.join(rnd.choice('ACDEFGHIKLMNPQRSTVWY') for _ in range(rnd.randint(5, 150)))))
        elif r < 0.7: quals.append(('EC_number', '"%s"' % rnd.choice(['1.1.1.1', '3.4.21.4', '2.7.-.-'])))
        elif r < 0.74: quals.append(('evidence', rnd.choice(['experimental', 'not_experimental'])))
        elif r < 0.77: quals.append(('pseudo', None))
        elif r < 0.8: quals.append(('standard_name', '"%s"' % rnd.choice(['std1', 'alpha factor'])))
        elif r < 0.83: quals.append(('function', '"%s"' % ' '.join(rnd.choice(WORDS) for _ in range(rnd.randint(1, 4)))))
        elif r < 0.86: quals.append(('number', str(rnd.randint(1, 9))))
        elif r < 0.89: quals.append(('locus_tag', '"%s"' % rnd.choice(['b0001', 'TAG_00012', 'x'])))
        elif r < 0.92: quals.append(('Gene', '"%s"' % rnd.choice(['upper', 'Case'])))
        elif r < 0.96: quals.append((rnd.choice(['foo', 'unknown_tag', 'myqual']), '"%s"' % rnd.choice(['bar', 'baz qux'])))
        elif r < 0.98: quals.append(('organism', '"Homo sapiens"'))
        elif r < 0.99: quals.append(rnd.choice([('note', '"say ""hello"" there"'), ('note', '""'), ('note', '"a/b c/d"'), ('gene', '"x\\y"'), ('note', 'unquoted'), ('Note', '"Mixed Case"'), ('note', '"line one"')]))
        else: quals.append(('allele', '"%s"' % rnd.choice(['a1', 'A2'])))
    return quals


def feature_lines(rnd, fmt, n):
    """the feature table lines of an entry"""
    lines = []
    for _ in range(rnd.choice([0, 1, 2, 3, 5, 8])):
        ty = rnd.choice(TYPES); loc = locations(rnd, n)
        quals = qualifiers(rnd)
        if fmt == 'genbank':
            lines.append('     %-16s%s' % (ty, loc) if len(ty) < 16 else '     %s %s' % (ty, loc))
            for name, val in quals:
                text = '/%s' % name + ('' if val is None else '=%s' % val)
                for w in wrap(text, 58): lines.append(' ' * 21 + w)
        else:
            lines.append('FT   %-16s%s' % (ty, loc) if len(ty) < 16 else 'FT   %s %s' % (ty, loc))
            for name, val in quals:
                text = '/%s' % name + ('' if val is None else '=%s' % val)
                for w in wrap(text, 58): lines.append('FT' + ' ' * 19 + w)
    return lines


def entry(rnd, fmt, name, desc, seq, features):
    n = len(seq)
    if fmt == 'genbank':
        head = 'LOCUS       %s %d bp    DNA     linear   UNK 01-JAN-2000\nDEFINITION  %s\nACCESSION   %s\nVERSION     %s.1\nKEYWORDS    .\n' % (name, n, desc, name, name)
        if rnd.random() < 0.9:
            head += 'SOURCE      synthetic\n  ORGANISM  synthetic\n'
            if rnd.random() < 0.6: head += '            Eukaryota; Metazoa;\n            Chordata.\n'
            # the line after the SOURCE block is skipped by the program: a FEATURES line there loses the features
            if rnd.random() < 0.8: head += 'REFERENCE   1  (bases 1 to %d)\n  AUTHORS   Smith,J.\n  TITLE     A title\n  JOURNAL   Nature 1 (1), 1-2 (2000)\n' % n
            if rnd.random() < 0.3: head += 'COMMENT     a comment\n            continued.\n'
        o = head + 'FEATURES             Location/Qualifiers\n'
        o += '     source          1..%d\n' % n
        o += ''.join(l + '\n' for l in features)
        o += 'ORIGIN\n'
        for i in range(0, n, 60):
            o += '%9d %s\n' % (i + 1, ' '.join(seq[i + j:i + j + 10].lower() for j in range(0, 60, 10) if seq[i + j:i + j + 10]))
        return o + '//\n' + ('LOCUS       OTHER 10 bp    DNA     linear   UNK 01-JAN-2000\nORIGIN\n        1 acgtacgtac\n//\n' if rnd.random() < 0.05 else '')
    o = 'ID   %s; SV 1; linear; genomic DNA; STD; UNC; %d BP.\nXX\nAC   %s;\nXX\nDE   %s\nXX\nFH   Key             Location/Qualifiers\nFH\n' % (name, n, name, desc)
    o += 'FT   source          1..%d\n' % n
    o += ''.join(l + '\n' for l in features)
    o += 'XX\nSQ   Sequence %d BP; 0 A; 0 C; 0 G; 0 T; 0 other;\n' % n
    for i in range(0, n, 60):
        chunk = seq[i:i + 60]
        o += '     %-66s%9d\n' % (' '.join(chunk[j:j + 10].lower() for j in range(0, 60, 10) if chunk[j:j + 10]), min(i + 60, n))
    return o + '//\n'


def main(out, n, seed):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        shim = os.path.join(d, 'faketime.so')
        subprocess.run(['gcc', '-shared', '-fPIC', '-o', shim, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'faketime.c')], check=True)
        while len(rows) < n:
            nuc, a, b = True, '', ''
            while not a or not b or len(a) < 40:
                nuc, a, b = D.case(rnd)
                nuc = True
                if not all(c in 'ACGTacgtN' for c in a + b):
                    a = ''.join(rnd.choice('ACGT') for _ in range(rnd.randint(60, 300))); b = D.mutate(rnd, a, 'ACGT', rnd.randint(1, 6))
            n1 = rnd.choice(['SEQA', 'AB123456', 'gb_a', 'X56789']); n2 = rnd.choice(['SEQB', 'sample2', 'chr7', 'Y99999'])
            fmt1 = rnd.choice(['genbank', 'genbank', 'embl']); fmt2 = rnd.choice(['genbank', 'embl', 'embl', 'fasta'])
            desc1 = rnd.choice(['first sequence', 'a gene region', 'x']); desc2 = rnd.choice(['second sequence', 'other', 'y'])
            f1 = entry(rnd, fmt1, n1, desc1, a, feature_lines(rnd, fmt1, len(a)))
            f2 = D.fasta(n2, desc2, b) if fmt2 == 'fasta' else entry(rnd, fmt2, n2, desc2, b, feature_lines(rnd, fmt2, len(b)))
            ws = rnd.choice([3, 4, 5, 6, 8, 10]); glob = rnd.random() < 0.3
            rq = []
            if rnd.random() < 0.35:
                for k, length in ((1, len(a)), (2, len(b))):
                    if rnd.random() < 0.7: rq += ['-sbegin%d' % k, str(rnd.choice([rnd.randint(1, max(1, length // 3)), -rnd.randint(1, max(1, length // 2))]))]
                    if rnd.random() < 0.7: rq += ['-send%d' % k, str(rnd.choice([rnd.randint(max(1, length // 2), length), -rnd.randint(1, max(1, length // 3))]))]
                    if rnd.random() < 0.35: rq.append('-sreverse%d' % k)
            if rnd.random() < 0.2: rq += ['-sformat1', {'genbank': 'genbank', 'embl': 'embl'}[fmt1]]
            rq += ['-rformat', rnd.choice(['diffseq'] * 6 + ['table', 'embl', 'gff', 'excel', 'tagseq', 'simple', 'feattable'])]
            if rnd.random() < 0.15: rq += ['-raccshow', '-rdesshow']
            if rnd.random() < 0.1: rq.append('-rusashow')
            ext = {'genbank': 'gb', 'embl': 'embl', 'fasta': 'fa'}
            fa, fb = 'a.' + ext[fmt1], 'b.' + ext[fmt2]
            args = [fa, fb, '-wordsize', str(ws)] + (['-globaldifferences'] if glob else []) + rq + ['-outfile', 'o.rep', '-aoutfeat', 'a.gff', '-boutfeat', 'b.gff', '-auto']
            for nme in ('o.rep', 'a.gff', 'b.gff'):
                if os.path.exists(os.path.join(d, nme)): os.remove(os.path.join(d, nme))
            open(os.path.join(d, fa), 'w', newline='\n').write(f1); open(os.path.join(d, fb), 'w', newline='\n').write(f2)
            epoch = rnd.randint(978307200, 2051222400)
            env = dict(os.environ, LD_PRELOAD=shim, FAKE_EPOCH=str(epoch))
            try:
                p = subprocess.run(['diffseq'] + args, cwd=d, capture_output=True, errors='replace', env=env, timeout=20)
            except subprocess.TimeoutExpired:
                continue
            def read(nme):
                q = os.path.join(d, nme)
                return open(q, encoding='latin-1', newline='').read() if os.path.exists(q) else ''
            if p.returncode != 0 and not read('o.rep'): continue
            entries = ['[-asequence] ' + fa, '[-bsequence] ' + fb, '-wordsize %d' % ws] + (['-globaldifferences'] if glob else [])
            i = 0
            while i < len(rq):
                if i + 1 < len(rq) and (not rq[i + 1].startswith('-') or rq[i + 1].lstrip('-').isdigit()): entries.append(rq[i] + ' ' + rq[i + 1]); i += 2
                else: entries.append(rq[i]); i += 1
            entries += ['-outfile o.rep', '-aoutfeat a.gff', '-boutfeat b.gff', '-auto']
            rows.append({'id': len(rows), 'nuc': nuc, 'file1': f1, 'file2': f2, 'name1': fa, 'name2': fb, 'wordsize': ws, 'glob': glob, 'rq': ' '.join(rq), 'cmd': entries, 'epoch': epoch,
                         'report': read('o.rep'), 'feat1': read('a.gff'), 'feat2': read('b.gff'), 'error': p.stderr.strip() or None, 'rc': p.returncode})
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    print(len(rows), 'rows', sum(1 for r in rows if 'Feature:' in r['report'] or 'first_feature' in r['report']), 'with features in the report')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 300, int(sys.argv[3]) if len(sys.argv) > 3 else 20261224)
