#!/usr/bin/env python3
"""B09 WP60 oracle: sequence DATABASES whose access method reads an INDEX built by the EMBOSS 6.6.0 indexers, queried through the installed `seqret` / `entret`:
    method emboss   : dbxfasta, dbxflat (b+tree files *.ent *.xid *.xac *.xsv *.xde *.xkw *.xtx and their *.p* parameter files; ajindex.c, ajtextdb.c textAccessEmboss)
    method emblcd   : dbifasta, dbiflat (EMBL CD-ROM index files division.lkp entrynam.idx acnum.trg acnum.hit ...; ajtextdb.c textAccessEmblcd)
    method gcg      : dbigcg (the same index files over GCG .ref / .seq pairs; ajseqdb.c seqAccessGcg)
    method embossgcg: dbxgcg (ajseqdb.c seqAccessEmbossGcg)
    seqret -auto [-feature] -sequence SPEC [-sformat1 F] ... -outseq o.fa -osformat2 fasta [-offormat2 gff3 -ofname2 o.gff]
The files of a row are the database files and the index files (as produced by the indexer: latin-1 text), the configuration (`.embossrc`) and the query; the reference is run with EMBOSSRC set to the work directory.
JSONL.GZ rows: id, how (method/indexer/format), definitions, files ({path: text}), spec, args, feature, epoch, rc, fa, gff, err, hang, memerr, program ('seqret' or 'entret'), entret (the file entret wrote), log (the indexer log)
Usage: python3 scripts/alignment/generate_dbindex_oracle.py OUT.jsonl.gz COUNT SEED METHOD      (METHOD: emboss | emblcd | gcg | embossgcg ; needs EMBOSS 6.6.0, gcc)
"""
import gzip, json, os, random, re, shutil, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_seqret_oracle as G  # noqa: E402
import generate_usa_oracle as U  # noqa: E402

WORDS = ['kinase', 'protein', 'alpha', 'beta', 'receptor', 'binding', 'human', 'mouse', 'putative', 'zinc', 'finger', 'transcription', 'factor', 'ribosomal', 'subunit', 'hypothetical']
ORGS = ['Homo sapiens', 'Mus musculus', 'Escherichia coli', 'Saccharomyces cerevisiae', 'Arabidopsis thaliana', 'Drosophila melanogaster']
KEYS = ['KINASE', 'TRANSFERASE', 'ATP-BINDING', 'RIBOSOMAL PROTEIN', 'ZINC-FINGER', 'MEMBRANE', 'HYDROLASE']


def seq_of(rnd, nuc, n=None):
    n = n or rnd.choice([12, 30, 61, 120, 200])
    return ''.join(rnd.choice('ACGT' if nuc else 'ACDEFGHIKLMNPQRSTVWY') for _ in range(n))


def names(rnd, count):
    pool = ['ABC', 'XYZ', 'seq%d', 'AB%06d', 'Q%05d', 'P%05d', 'ID_%d', 'hsp%d', 'GENE%d', 'a%d', 'T%d']
    out = []
    while len(out) < count:
        p = rnd.choice(pool)
        s = p % rnd.randint(1, 99999) if '%' in p else p + str(rnd.randint(1, 99))
        if s.lower() not in [x.lower() for x in out]:
            out.append(s)
    return out


def fasta_entries(rnd, fmt, nuc, count):
    """(text, entries) in the header form of the dbxfasta id format"""
    out, ents = '', []
    for nm in names(rnd, count):
        acc = rnd.choice(['XX%06d' % rnd.randint(1, 999999), 'P%05d' % rnd.randint(1, 99999), 'A%d' % rnd.randint(1, 9999)])
        sv = acc + '.%d' % rnd.randint(1, 9)
        gi = str(rnd.randint(1000, 99999999))
        desc = ' '.join(rnd.choice(WORDS) for _ in range(rnd.randint(0, 5)))
        if fmt == 'simple': head = nm + (' ' + desc if desc else '')
        elif fmt == 'idacc': head = nm + rnd.choice([' %s', ' (%s)']) % acc + (' ' + desc if desc else '')
        elif fmt == 'idsv': head = nm + rnd.choice([' %s', ' (%s)']) % sv + (' ' + desc if desc else '')
        elif fmt == 'gcgid': head = 'db:' + nm + (' ' + desc if desc else '')
        elif fmt == 'gcgidacc': head = 'db:%s %s' % (nm, acc) + (' ' + desc if desc else '')
        elif fmt == 'dbid': head = 'db %s' % nm + (' ' + desc if desc else '')
        else: head = 'gi|%s|gb|%s|%s' % (gi, sv, nm) + (' ' + desc if desc else '')
        s = seq_of(rnd, nuc)
        out += '>%s\n%s\n' % (head, '\n'.join(s[j:j + rnd.choice([60, 70, 80])] for j in range(0, len(s), 60)))
        ents.append({'id': nm, 'acc': acc, 'sv': sv, 'gi': gi, 'des': desc})
    return out, ents


def flat_entry(rnd, fmt, nm, nuc):
    acc = rnd.choice(['XX%06d' % rnd.randint(1, 999999), 'P%05d' % rnd.randint(1, 99999), 'A%d' % rnd.randint(1, 9999)])
    acc2 = 'Q%05d' % rnd.randint(1, 99999)
    ver = rnd.randint(1, 5)
    gi = rnd.randint(1000, 9999999)
    desc = ' '.join(rnd.choice(WORDS) for _ in range(rnd.randint(1, 5))).capitalize() + '.'
    org = rnd.choice(ORGS)
    kws = rnd.sample(KEYS, rnd.randint(0, 3))
    s = seq_of(rnd, nuc)
    n = len(s)
    lines = []
    if fmt == 'embl':
        lines.append('ID   %s; SV %d; linear; genomic DNA; STD; HUM; %d BP.' % (nm, ver, n)); lines.append('XX')
        lines.append('AC   %s;%s' % (acc, ' ' + acc2 + ';' if rnd.random() < 0.4 else '')); lines.append('XX')
        lines.append('SV   %s.%d' % (acc, ver)); lines.append('XX')
        lines.append('DE   ' + desc); lines.append('XX')
        if kws: lines.append('KW   ' + '; '.join(kws) + '.'); lines.append('XX')
        lines.append('OS   ' + org); lines.append('OC   Eukaryota; Metazoa.'); lines.append('XX')
        lines.append('SQ   Sequence %d BP; %d A; %d C; %d G; %d T; 0 other;' % (n, s.count('A'), s.count('C'), s.count('G'), s.count('T')))
        low = s.lower()
        for i in range(0, n, 60):
            ch = low[i:i + 60]
            lines.append('     %-66s%9d' % (' '.join(ch[j:j + 10] for j in range(0, len(ch), 10)), min(n, i + 60)))
        lines.append('//')
    elif fmt == 'swiss':
        lines.append('ID   %s     STANDARD;      PRT;   %d AA.' % (nm, n)); lines.append('AC   %s;%s' % (acc, ' ' + acc2 + ';' if rnd.random() < 0.4 else ''))
        lines.append('DT   01-JAN-1990 (Rel. 1, Created)'); lines.append('DE   ' + desc)
        lines.append('OS   ' + org + '.'); lines.append('OC   Eukaryota; Metazoa.')
        if kws: lines.append('KW   ' + '; '.join(kws) + '.')
        lines.append('SQ   SEQUENCE   %d AA;  %d MW;  0 CRC64;' % (n, n * 110))
        for i in range(0, n, 60):
            ch = s[i:i + 60]
            lines.append('     ' + ' '.join(ch[j:j + 10] for j in range(0, len(ch), 10)))
        lines.append('//')
    else:   # gb, refseq
        lines.append('LOCUS       %-16s %d bp    DNA     linear   PRI 01-JAN-1990' % (nm, n))
        lines.append('DEFINITION  ' + desc)
        lines.append('ACCESSION   %s%s' % (acc, ' ' + acc2 if rnd.random() < 0.3 else ''))
        lines.append('VERSION     %s.%d  GI:%d' % (acc, ver, gi))
        if kws: lines.append('KEYWORDS    ' + '; '.join(kws) + '.')
        lines.append('SOURCE      ' + org); lines.append('  ORGANISM  ' + org); lines.append('            Eukaryota; Metazoa.')
        lines.append('FEATURES             Location/Qualifiers')
        lines.append('     source          1..%d' % n)
        lines.append('ORIGIN')
        low = s.lower()
        for i in range(0, n, 60):
            ch = low[i:i + 60]
            lines.append('%9d %s' % (i + 1, ' '.join(ch[j:j + 10] for j in range(0, len(ch), 10))))
        lines.append('//')
    return '\n'.join(lines) + '\n', {'id': nm, 'acc': acc, 'sv': '%s.%d' % (acc, ver), 'gi': str(gi), 'des': desc, 'org': org, 'key': kws}


def run(prog, args, cwd, env, timeout=60):
    try:
        p = subprocess.run([prog] + args, cwd=cwd, capture_output=True, timeout=timeout, env=env)
        return p.returncode, p.stdout.decode('latin-1'), p.stderr.decode('latin-1')
    except subprocess.TimeoutExpired:
        return -9, '', ''


def collect(work, skip=()):
    files = {}
    for root, _, ns in os.walk(work):
        for n in sorted(ns):
            rel = os.path.relpath(os.path.join(root, n), work)
            if rel in skip or rel == '.embossrc':
                continue
            files[rel] = open(os.path.join(root, n), 'rb').read().decode('latin-1')
    return files


def clean(work):
    for nme in os.listdir(work):
        p = os.path.join(work, nme)
        shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)


QUERIES = ['{db}:{id}', '{db}:{id}', '{db}:{id}', '{db}:{idw}', '{db}:*', '{db}-id:{id}', '{db}-acc:{acc}', '{db}:acc:{acc}', '{db}-sv:{sv}', '{db}-des:{word}', '{db}:des:{word}', '{db}-des:{word}*', '{db}:{{{id}|{id2}}}',
           '{db}:{id},{id2}', '{db}:{id}|{id2}', '{db}-id:{id}&acc:{acc}', '{db}-id:{id}|acc:{acc}', '{db}-id:{id}!acc:{acc}', '{db}-acc:{acc}^id:{id}', '{db}:{id}[2:9]', '{db}:{id}[3:20:r]', '{DB}:{id}', '{db}:{ID}', '{db}:nosuchid',
           '{db}-org:{org}', '{db}-org:{orgw}', '{db}-key:{key}', '{db}-key:{keyw}', '{db}-gi:{gi}', '{db}:{tailw}', '{db}-des:*{word}*', '{db}-id:{idw}&des:{word}', '{db}:{acc}', '{db}:{sv}']


def build(rnd, method, work, shim, files, defs):
    """one database: returns (info dict) or None"""
    kind = 'N' if rnd.random() < 0.75 else 'P'
    nuc = kind == 'N'
    name = rnd.choice(['tdb', 'mydb', 'db1', 'seqs', 'a_b', 'swdb', 'embl2', 'tdb', 'mydb'] + (['Prot_db'] if rnd.random() < 0.1 else []))
    dbdir, idxdir = rnd.choice(['db', 'data/sub']), rnd.choice(['idx', 'index/x', 'db'])
    nfiles = rnd.choice([1, 1, 1, 2, 3])
    pagesize = rnd.choice([2048, 2048, 2048, 4096, 1024])
    compressed = rnd.random() < 0.6
    fields_all = ['id', 'acc', 'sv', 'des']
    ents, fmtname, indexer, ifmt, ext = [], '', '', '', ''
    texts = {}
    if method in ('emboss', 'emblcd') and rnd.random() < 0.45:
        indexer = 'fasta'
        ifmt = rnd.choice(['simple', 'idacc', 'idsv', 'gcgid', 'gcgidacc', 'dbid', 'ncbi']) if method == 'emboss' else rnd.choice(['simple', 'idacc', 'idsv', 'gcgid', 'gcgidacc', 'dbid', 'ncbi'])
        fmtname, ext = 'fasta', 'fasta'
        for i in range(nfiles):
            t, e = fasta_entries(rnd, ifmt, nuc, rnd.choice([1, 2, 3, 5, 8, 15]))
            texts['%s/f%d.fasta' % (dbdir, i)] = t
            ents += e
    else:
        indexer = 'flat'
        ifmt = rnd.choice(['embl', 'swiss', 'gb', 'refseq']) if method == 'emboss' else rnd.choice(['embl', 'swiss', 'genbank'])
        fmtname = {'embl': 'embl', 'swiss': 'swiss', 'gb': 'genbank', 'genbank': 'genbank', 'refseq': 'refseq'}[ifmt]
        nuc = ifmt != 'swiss'
        kind = 'N' if nuc else 'P'
        ext = 'dat'
        for i in range(nfiles):
            t = ''
            for nm in names(rnd, rnd.choice([1, 2, 3, 5, 8, 15])):
                et, e = flat_entry(rnd, ifmt, nm, nuc)
                t += et; ents.append(e)
            texts['%s/f%d.dat' % (dbdir, i)] = t
    if len({e['id'].lower() for e in ents}) != len(ents):
        return None
    if indexer == 'fasta':
        fields = ['id'] + rnd.sample(['acc', 'sv', 'des'], rnd.randint(0, 3))
    else:
        fields = ['id'] + rnd.sample(['acc', 'sv', 'des', 'key', 'org'], rnd.randint(0, 5))
    resname = name
    res = 'RES %s [ type: Index\n  idlen: %d\n  acclen: %d\n  svlen: %d\n  keylen: %d\n  deslen: %d\n  orglen: %d\n]\n' % (resname, rnd.choice([15, 15, 20]), rnd.choice([15, 20]), rnd.choice([15, 20]), rnd.choice([15, 20]), rnd.choice([15, 25, 30]), rnd.choice([15, 25, 30]))
    for nme in os.listdir(work):
        p = os.path.join(work, nme)
        shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)
    os.makedirs(os.path.join(work, idxdir), exist_ok=True)
    for pth, t in texts.items():
        os.makedirs(os.path.dirname(os.path.join(work, pth)), exist_ok=True)
        open(os.path.join(work, pth), 'w', newline='\n', encoding='latin-1').write(t)
    cfg = 'SET PAGESIZE %d\nSET CACHESIZE %d\n' % (pagesize, rnd.choice([50, 200])) + res
    env = dict(os.environ, EMBOSSRC=work, EMBOSS_RCHOME='N', LD_PRELOAD=shim, FAKE_EPOCH='1700000000')
    open(os.path.join(work, '.embossrc'), 'w').write(cfg)
    attrs = ['type: %s' % rnd.choice([kind, 'Nucleotide' if nuc else 'Protein']), 'format: %s' % fmtname, 'directory: "%s"' % dbdir, 'indexdirectory: "%s"' % idxdir, 'filename: "*.%s"' % ext]
    if rnd.random() < 0.3:
        attrs.append('release: "1.0"')
    if rnd.random() < 0.15:
        attrs.append('caseidmatch: %s' % rnd.choice(['Y', 'N']))
    if rnd.random() < 0.1:
        attrs.append('hasaccession: %s' % rnd.choice(['Y', 'N']))
    dbmethod = method
    if method == 'emboss':
        prog = 'dbxfasta' if indexer == 'fasta' else 'dbxflat'
        args = ['-auto', '-dbname', name, '-dbresource', resname, '-idformat', ifmt if indexer == 'fasta' else ifmt.upper() if ifmt != 'genbank' else 'GB', '-directory', dbdir, '-filenames', '*.' + ext, '-fields', ','.join(fields),
                '-compressed' if compressed else '-nocompressed', '-outfile', 'index.log', '-indexoutdir', idxdir]
        if rnd.random() < 0.2: args += ['-release', '1.0']
        attrs.append('method: emboss')
        attrs.append('fields: "%s"' % ' '.join(fields)) if rnd.random() < 0.7 else None
    else:
        prog = 'dbifasta' if indexer == 'fasta' else 'dbiflat'
        args = ['-auto', '-dbname', name, '-idformat', ifmt if indexer == 'fasta' else ifmt.upper(), '-directory', dbdir, '-filenames', '*.' + ext, '-fields', ','.join(f for f in fields if f != 'id'), '-outfile', 'index.log', '-indexoutdir', idxdir]
        attrs.append('method: emblcd')
        attrs.append('fields: "%s"' % ' '.join(f for f in fields if f != 'id')) if rnd.random() < 0.7 else None
    attrs = [a for a in attrs if a]
    rnd.shuffle(attrs)
    defs = cfg + 'DBNAME %s [\n  %s\n]\n' % (name, '\n  '.join(attrs))
    open(os.path.join(work, '.embossrc'), 'w').write(defs)
    viacompress = method == 'emboss' and compressed and rnd.random() < 0.25       # an uncompressed index compressed by dbxcompress, field by field
    if viacompress:
        args = [('-nocompressed' if a == '-compressed' else a) for a in args]
    rc, out, err = run(prog, args, work, env)
    if rc != 0:
        return None
    if viacompress:
        for f in fields:
            rc2, _, _ = run('dbxcompress', ['-auto', '-dbname', name, '-field', f, '-outfile', 'cmp.log', '-indexdir', idxdir], work, env)
            if rc2 != 0:
                return None
        if os.path.exists(os.path.join(work, 'cmp.log')): os.remove(os.path.join(work, 'cmp.log'))
    for nme in ('index.log',):
        if os.path.exists(os.path.join(work, nme)): os.remove(os.path.join(work, nme))
    return {'name': name, 'ents': ents, 'defs': defs, 'files': collect(work), 'fields': fields, 'fmt': fmtname, 'how': '%s/%s/%s%s%s' % (method, prog, ifmt, '/c' if compressed else '', '+dbxcompress' if viacompress else '')}


def main(out, count, seed, method):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        shim = os.path.join(d, 'faketime.so')
        subprocess.run(['gcc', '-shared', '-fPIC', '-o', shim, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'faketime.c')], check=True)
        work = os.path.join(d, 'w'); os.mkdir(work)
        made = 0
        while made < count:
            info = build(rnd, method, work, shim, {}, '')
            if not info:
                continue
            ents = info['ents']
            e1, e2 = rnd.choice(ents), rnd.choice(ents)
            allw = [w for e in ents for w in e['des'].lower().replace('.', '').split()] or ['kinase']
            keys = [k for e in ents for k in e.get('key', [])] or ['KINASE']
            qs = [x for x in QUERIES if rnd.random() < 0.1 or all((f in info['fields']) for f in re.findall(r'[-:&|!^](id|acc|sv|des|org|key|gi)[:]', x.replace('{db}', 'db')) if f != 'gi') and not ('{acc}' in x and 'acc' not in info['fields'] and False)]
            q = rnd.choice(qs).format(db=info['name'], DB=info['name'].upper(), id=e1['id'], ID=e1['id'].upper(), id2=e2['id'], idw=e1['id'][:max(1, len(e1['id']) // 2)] + '*', acc=e1['acc'], sv=e1['sv'], gi=e1['gi'],
                                          word=rnd.choice(allw), org=e1.get('org', 'Homo sapiens'), orgw=e1.get('org', 'Homo')[:4] + '*', key=rnd.choice(keys), keyw=rnd.choice(keys)[:3] + '*', tailw='*' + e1['id'][-2:])
            spec = q
            if rnd.random() < 0.1:
                spec = '%s::%s' % (rnd.choice([info['fmt'], 'fasta', 'embl']), q)
            args = []
            if rnd.random() < 0.2: args += ['-sformat1', rnd.choice([info['fmt'], 'fasta', 'embl'])]
            if rnd.random() < 0.15: args += ['-sbegin1', rnd.choice(['2', '5', '-6'])]
            if rnd.random() < 0.15: args += ['-send1', rnd.choice(['8', '-3', '40'])]
            if rnd.random() < 0.1: args += ['-sreverse1']
            if rnd.random() < 0.08: args += ['-firstonly']
            if rnd.random() < 0.05: args += ['-sdbname1', 'xdb']
            feature = rnd.random() < 0.3
            program = 'seqret' if rnd.random() < 0.8 else 'entret'
            epoch = rnd.randint(978307200, 2051222400)
            env = dict(os.environ, LD_PRELOAD=shim, FAKE_EPOCH=str(epoch), EMBOSSRC=work, EMBOSS_RCHOME='N')
            if program == 'seqret':
                cmd = ['-auto'] + (['-feature'] if feature else []) + ['-sequence', spec] + args + ['-outseq', 'o.fa', '-osformat2', 'fasta'] + (['-offormat2', 'gff3', '-ofname2', 'o.gff'] if feature else [])
            else:
                feature = False
                cmd = ['-auto', '-sequence', spec] + [a for a in args if a == '-firstonly'] + ['-outfile', 'o.fa']
            rc, _, err = run(program, cmd, work, env, timeout=10)
            hang = rc in (-9, -25)
            if 'ajSysCreate' in err:
                continue
            rd = lambda n: open(os.path.join(work, n), encoding='latin-1', newline='').read() if os.path.exists(os.path.join(work, n)) and not hang else ''
            if err.startswith('Read and write (return) sequences\n'):
                err = err[len('Read and write (return) sequences\n'):]
            if err.startswith('Retrieve sequence entries from flatfile databases and files\n'):
                err = err[len('Retrieve sequence entries from flatfile databases and files\n'):]
            memerr = False
            rows.append({'id': len(rows), 'how': info['how'], 'definitions': info['defs'], 'files': info['files'], 'spec': spec, 'args': args if program == 'seqret' else [a for a in args if a == '-firstonly'], 'feature': feature, 'epoch': epoch, 'rc': rc,
                         'fa': rd('o.fa'), 'gff': rd('o.gff'), 'err': err, 'hang': hang, 'memerr': memerr, 'program': program})
            made += 1
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows:
            fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    from collections import Counter
    print(len(rows), 'rows;', sorted(Counter((r['how'].split('/')[1], 'hang' if r['hang'] else 'died' if r['rc'] else 'fa' if r['fa'] else 'none') for r in rows).items()))


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 20266010, sys.argv[4] if len(sys.argv) > 4 else 'emboss')
