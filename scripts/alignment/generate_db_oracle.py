#!/usr/bin/env python3
"""A5-8c3h (B09 F256, WP58) oracle: sequence DATABASES defined in an EMBOSS configuration file (`.embossrc`: `DBNAME name [ type: ... format: ... method: direct directory: ... filename: ... fields: ... ]`, `SET var value`, `$var`)
and queried by `db:id`, `db-acc:value`, `db:id1,id2`, `db:{a|b}`, `db:field:value`, `db:*`, `db:wild*` through the installed EMBOSS 6.6.0 `seqret` (ajquery.c `ajQuerystrParseRead`, ajnam.c `ajNamDbData` / `ajNamDbQuery`,
ajtextdb.c `textAccessDirect`):
    EMBOSSRC=work seqret -auto [-feature] -sequence SPEC [-sformat1 F] [-sbegin1 N] [-send1 N] [-sreverse1] -outseq o.fa -osformat2 fasta [-offormat2 gff3 -ofname2 o.gff]
The database files are made from the inputs of the reader oracle (writer outputs of seqret, entries with feature tables). The method `direct` reads the files of the directory that match `filename` (not those of `exclude`) in the order of the
directory listing, which depends on the file system: a database of several files is used only for queries whose result does not depend on that order (the row is run with the files created in two orders and dropped when the outputs differ).
JSONL.GZ rows: id, how, definitions (the text of the configuration file), files ({path: text}), spec, args, feature, epoch, rc, fa, gff, err, hang, memerr.
Usage: python3 scripts/alignment/generate_db_oracle.py OUT.jsonl.gz COUNT [seed]      (needs EMBOSS 6.6.0)
"""
import fnmatch, gzip, json, os, random, re, shutil, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_seqret_oracle as G  # noqa: E402
import generate_usa_oracle as U  # noqa: E402

DBFMT = [('fasta', 'N', 'fa'), ('fasta', 'P', 'fa'), ('embl', 'N', 'embl'), ('genbank', 'N', 'gb'), ('swiss', 'P', 'sw'), ('nbrf', 'N', 'pir'), ('gcg', 'N', 'gcg'), ('codata', 'N', 'cod'), ('ig', 'N', 'ig'), ('fastq', 'N', 'fq'), ('gff3', 'N', 'gff')]


def db_text(rnd, fmt, kind):
    """text of a database file of format fmt (several entries)"""
    if fmt == 'fasta':
        return U.multi_fasta(rnd, kind == 'N')
    for _ in range(4):
        t, how = G.make_input(rnd, GEN, fmt)
        if t.strip():
            return t
    return U.multi_fasta(rnd, kind == 'N')


GEN = None


def definition(name, fmt, kind, directory, filename, rnd):
    attrs = ['type: %s' % rnd.choice([kind, 'Nucleotide' if kind == 'N' else 'Protein', 'Sequence']), 'format: %s' % fmt, 'method: direct', 'directory: "%s"' % directory]
    attrs.append('filename: "%s"' % filename)
    if rnd.random() < 0.3:
        attrs.append('fields: "%s"' % rnd.choice(['sv des org key', 'des', 'sv gi des', 'acc sv des org']))
    if rnd.random() < 0.15:
        attrs.append('caseidmatch: %s' % rnd.choice(['Y', 'N']))
    if rnd.random() < 0.15:
        attrs.append('hasaccession: %s' % rnd.choice(['Y', 'N']))
    if rnd.random() < 0.2:
        attrs.append('release: "1.0"')
        attrs.append('comment: "test database"')
    if rnd.random() < 0.1:
        attrs.append('exclude: "skip*"')
    rnd.shuffle(attrs)
    return 'DBNAME %s [\n  %s\n]\n' % (name, '\n  '.join(attrs))


QUERIES = ['{db}:{id}', '{db}:{id}', '{db}:{idw}', '{db}:*', '{db}-id:{id}', '{db}-acc:{acc}', '{db}:acc:{acc}', '{db}-des:{word}', '{db}:des:{word}', '{db}:{{{id}|{id2}}}', '{db}:{id},{id2}', '{db}:{id}|{id2}', '{db}-sv:{sv}',
           '{db}:{id}[2:9]', '{db}:{id}[3:20:r]', '{DB}:{id}', '{db}:{ID}', '{db}:nosuchid', '{db}', '{db}:', '{db}-nofield:x', '{db}-org:Homo*', '{db}-key:kinase*', '{db}:*{tail}']


def main(out, count, seed):
    global GEN
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        shim = os.path.join(d, 'faketime.so')
        subprocess.run(['gcc', '-shared', '-fPIC', '-o', shim, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'faketime.c')], check=True)
        work = os.path.join(d, 'w'); os.mkdir(work)
        GEN = os.path.join(d, 'g'); os.mkdir(GEN)
        made = 0
        while made < count:
            for nme in os.listdir(work):
                p = os.path.join(work, nme)
                shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)
            files, defs, dbs = {}, 'SET tdir db\n' if rnd.random() < 0.4 else '', []
            ndb = rnd.choice([1, 1, 2, 3])
            for k in range(ndb):
                fmt, kind, ext = rnd.choice(DBFMT)
                name = rnd.choice(['tdb', 'mydb', 'Prot_db', 'db1', 'seqs', 'tdb2', 'a_b'])
                name = '%s%d' % (name, k) if k else name
                if name in [x['name'] for x in dbs]:
                    continue
                dirname = rnd.choice(['db%d' % k, 'data%d/sub' % k])
                nfiles = rnd.choice([1, 1, 1, 2, 3])
                fnames = ['%s%d.%s' % (rnd.choice(['f', 'part', 'x']), i, ext) for i in range(nfiles)]
                if nfiles > 1 and rnd.random() < 0.2:
                    fnames.append('skip1.%s' % ext)
                filename = rnd.choice(['*.%s' % ext, '*', fnames[0], 'f*.%s' % ext, 'part?.%s' % ext]) if nfiles > 1 else rnd.choice(['*.%s' % ext, fnames[0], '*'])
                texts = {}
                for fn in fnames:
                    texts[os.path.join(dirname, fn)] = db_text(rnd, fmt, kind)
                    if rnd.random() < 0.1 and nfiles == 1:
                        texts[os.path.join(dirname, fn)] = G.damaged(rnd, texts[os.path.join(dirname, fn)])
                directory = dirname
                dtext = definition(name, fmt, kind, directory, filename, rnd)
                defs += dtext
                excl = 'skip*' in dtext
                matching = {k: v for k, v in texts.items() if fnmatch.fnmatchcase(os.path.basename(k), filename) and not (excl and fnmatch.fnmatchcase(os.path.basename(k), 'skip*'))}
                files.update(texts)
                dbs.append({'name': name, 'fmt': fmt, 'kind': kind, 'texts': matching or texts, 'multi': len(matching) > 1})
            if not dbs:
                continue
            db = rnd.choice(dbs)
            # ids, accessions, words of the entries of the database file(s)
            ids, accs, svs, words = [], [], [], []
            for t in db['texts'].values():
                for m in re.finditer(r'^[>@]\s*(\S+)|^ID\s+(\S+?)[; ]|^LOCUS\s+(\S+)|^ENTRY\s+(\S+)|^>[A-Z0-9]{2};(\S+)', t, re.M):
                    ids.append(next(g for g in m.groups() if g).rstrip(';'))
                accs += re.findall(r'^AC\s+([A-Z0-9_.]+);|^ACCESSION\s+(\S+)', t, re.M) and [a for tup in re.findall(r'^AC\s+([A-Z0-9_.]+);|^ACCESSION\s+(\S+)', t, re.M) for a in tup if a]
                words += re.findall(r'[A-Za-z]{4,}', t[:400])
            ids = ids or ['a']
            accs = accs or ['XX000001']
            id1, id2 = rnd.choice(ids), rnd.choice(ids)
            queries = QUERIES
            if db['multi']:
                # the method direct reads the files in the order of the directory listing (the file system's): only a query for an identifier that occurs once is independent of it
                queries = ['{db}:{id}', '{db}-id:{id}', '{DB}:{id}']
                counts = {}
                for t in db['texts'].values():
                    for m in re.finditer(r'^[>@]\s*(\S+)|^ID\s+(\S+?)[; ]|^LOCUS\s+(\S+)|^ENTRY\s+(\S+)|^>[A-Z0-9]{2};(\S+)', t, re.M):
                        k = next(g for g in m.groups() if g).rstrip(';').lower()
                        counts[k] = counts.get(k, 0) + 1
                ids = [i for i in ids if counts.get(i.lower(), 0) == 1] or ids
            q = rnd.choice(queries).format(db=db['name'], DB=db['name'].upper(), id=id1, ID=id1.upper(), id2=id2, idw=id1[:max(1, len(id1) // 2)] + '*', acc=rnd.choice(accs), sv=rnd.choice(accs) + '.1', word=rnd.choice(words or ['desc']), tail=id1[-2:])
            spec = q
            if rnd.random() < 0.12:
                spec = '%s::%s' % (rnd.choice([db['fmt'], 'fasta', 'embl']), q)
            lstfile = None
            if rnd.random() < 0.12:
                files['l1.lst'] = '%s\n%s\n' % (q, '%s:%s' % (db['name'], rnd.choice(ids)))
                spec = '@l1.lst'
            args = []
            if rnd.random() < 0.25: args += ['-sformat1', rnd.choice([db['fmt'], 'fasta', 'embl'])]
            if rnd.random() < 0.2: args += ['-sbegin1', rnd.choice(['2', '5', '-6'])]
            if rnd.random() < 0.2: args += ['-send1', rnd.choice(['8', '-3', '40'])]
            if rnd.random() < 0.15: args += ['-sreverse1']
            if rnd.random() < 0.08: args += ['-sid1', rnd.choice(['a', 'zz'])]
            if rnd.random() < 0.1: args += ['-firstonly']
            if rnd.random() < 0.05: args += ['-sdbname1', 'xdb']
            feature = rnd.random() < 0.6
            epoch = rnd.randint(978307200, 2051222400)
            env = dict(os.environ, LD_PRELOAD=shim, FAKE_EPOCH=str(epoch), EMBOSSRC=work, EMBOSS_RCHOME='N')
            cmd = ['-auto'] + (['-feature'] if feature else []) + ['-sequence', spec] + args + ['-outseq', 'o.fa', '-osformat2', 'fasta'] + (['-offormat2', 'gff3', '-ofname2', 'o.gff'] if feature else [])

            def run(order):
                for nme in os.listdir(work):
                    p = os.path.join(work, nme)
                    shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)
                open(os.path.join(work, '.embossrc'), 'w', newline='\n').write(defs)
                items = list(files.items())
                if order:
                    items.reverse()
                for pth, t in items:
                    full = os.path.join(work, pth)
                    os.makedirs(os.path.dirname(full) or work, exist_ok=True)
                    open(full, 'w', newline='\n', encoding='latin-1').write(t)
                rc, err, hang = G.seqret(work, cmd, env=env)
                rd = lambda n: open(os.path.join(work, n), encoding='latin-1', newline='').read() if os.path.exists(os.path.join(work, n)) and not hang else ''
                return rc, err, hang, rd('o.fa'), rd('o.gff')

            r1 = run(False)
            rc, err, hang, fa, gff = r1
            if db['multi'] and db['fmt'] in ('gff3', 'gff2') and feature:
                continue                    # a GFF file without features read with its features: the library dies, and which file is read first decides it
            if db['multi'] and not hang:
                # the result must not depend on the order of the files: exactly one of the files of the database may hold the entry (the query on that file alone is a file query)
                hits = 0
                bad = False
                for pth in db['texts']:
                    rcv, errv, hangv = G.seqret(work, ['-auto', '-sequence', pth, '-sformat1', db['fmt'], '-outseq', 'chk.fa', '-osformat2', 'fasta'], env=env)
                    if rcv != 0 or hangv or 'Error' in errv or 'Warning' in errv:
                        bad = True              # a file that does not read cleanly in the format of the database makes the end of the read depend on the order of the files
                if bad:
                    continue
                qid = re.sub(r'^[^:\-]+[:\-](?:id:)?', '', q)
                for pth in db['texts']:
                    rc1, err1, hang1 = G.seqret(work, ['-auto', '-sequence', '%s:%s' % (pth, qid.split('[')[0]), '-sformat1', db['fmt'], '-outseq', 'chk.fa', '-osformat2', 'fasta'], env=env)
                    cp = os.path.join(work, 'chk.fa')
                    if rc1 == 0 and os.path.exists(cp) and os.path.getsize(cp) > 0:
                        hits += 1
                    if os.path.exists(cp): os.remove(cp)
                if hits != 1:
                    continue
            if err.startswith('Read and write (return) sequences\n'):
                err = err[len('Read and write (return) sequences\n'):]
            if 'ajSysCreate' in err:
                continue
            rows.append({'id': len(rows), 'how': db['fmt'], 'definitions': defs, 'files': files, 'spec': spec, 'args': args, 'feature': feature, 'epoch': epoch, 'rc': rc, 'fa': fa, 'gff': gff, 'err': err, 'hang': hang})
            made += 1
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows:
            fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    from collections import Counter
    print(len(rows), 'rows;', sorted(Counter('hang' if r['hang'] else 'died' if r['rc'] else 'fa' if r['fa'] else 'none' for r in rows).items()))


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 20265802)
