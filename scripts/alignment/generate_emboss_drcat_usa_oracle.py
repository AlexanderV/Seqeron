#!/usr/bin/env python3
"""B09 WP67 (F309): oracle of the USAs `db:id` of the databases of the data resource catalogue (DRCAT.dat) for seqret / entret of EMBOSS 6.6.0 (queryDbfind -> ajResourceGetDbdata: the database is found in the catalogue when one of its query lines has a
format the sequence reader knows; the entries are then fetched from the URL of the line). The reference runs without network access: a fetch fails ("Unable to read sequence"). With network access (the user's run, handoff/WP62-drcat-lookup) the fetch of some
databases does not return within 60 s: such rows depend on the network and are not part of this oracle (see the F309 entry of B09.md).
JSONL.GZ rows as emboss_prog_oracle_lib.py (prog, args, files, rc, out, err, outfiles, hang, epoch).
Usage: python3 scripts/alignment/generate_emboss_drcat_usa_oracle.py OUT.jsonl.gz [n=600] [seed]      (needs EMBOSS 6.6.0)"""
import gzip, json, os, random, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import emboss_prog_oracle_lib as L  # noqa: E402


def names():
    out = []
    for line in open('/usr/share/EMBOSS/data/DRCAT.dat', encoding='latin-1'):
        if line.startswith('ID      ') or line.startswith('IDalt   '):
            out.append(line.split(None, 1)[1].strip())
    return out


def main():
    out = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 600
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 20262601
    rnd = random.Random(seed)
    pool = names()
    forms = ['%s:X1', '%s:J00231', '%s-id:X1', '%s-acc:X1', '%s:X*', '%s:{X1}', '%s-id:{X1}', '%s:X1,X2', '%s::fasta', 'fasta::%s:X1', '%s']
    rows = []
    with tempfile.TemporaryDirectory() as d:
        work = os.path.join(d, 'w')
        os.mkdir(work)
        shim = L.clock_shim(d)
        for k in range(n):
            nm = rnd.choice(pool)
            if rnd.random() < 0.15: nm = nm.lower()
            elif rnd.random() < 0.1: nm = nm.upper()
            if rnd.random() < 0.05: nm = rnd.choice(['nosuchdb', 'embl_x', 'x', 'edam', 'go', 'taxon', 'drcat', 'chebi', 'swo', 'so', 'pw', 'ro', 'eco', 'srs', 'ensembl', 'biomart', 'das', 'mrs', 'entrez', 'dbfetch'])
            usa = rnd.choice(forms) % nm
            prog = rnd.choice(['seqret', 'seqret', 'seqret', 'entret'])
            args = ['-sequence', usa]
            if prog == 'seqret': args += ['-outseq', 'stdout']
            else: args += ['-outfile', 'stdout']
            files = {}
            if rnd.random() < 0.04: files[nm] = '>a\nACGT\n'          # a file of the name of the database
            epoch = random.Random(seed * 1000003 + k).randint(978307200, 2051222400)
            env = dict(os.environ, LD_PRELOAD=shim, FAKE_EPOCH=str(epoch), TZ='UTC')
            rc, o, e, outfiles, hang = L.run(prog, args, files, work, env=env)
            rc2, o2, e2, outfiles2, hang2 = L.run(prog, args, files, work, env=env)
            if hang or (rc, o, e, outfiles) != (rc2, o2, e2, outfiles2):
                continue
            rows.append({'id': len(rows), 'prog': prog, 'args': args, 'files': files, 'rc': rc, 'out': o, 'err': e, 'outfiles': outfiles, 'hang': hang, 'epoch': epoch})
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows:
            fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    from collections import Counter
    print(len(rows), 'rows;', Counter((r['rc'], r['err'].split('\n')[0][:60]) for r in rows).most_common(12))


if __name__ == '__main__':
    main()
