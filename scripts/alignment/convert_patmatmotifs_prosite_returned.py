#!/usr/bin/env python3
"""B09 WP67 (F306): converts the returned run of the Docker package docs/Validation/review-2026-09/handoff/WP62-patmatmotifs-prosite (patmatmotifs 6.6.0 on the real PROSITE release, 85 rows: standard output, error stream, exit status,
written files, valgrind verdict) into the test data
    tests/.../TestData/Alignment/emboss_prosite_returned_oracle.jsonl.gz     one JSON row per run: id, args, files (the sequence files of the run), rc, out, err, outfiles, epoch (the run date of the report header)
    tests/.../TestData/Alignment/emboss_prosite_release_data.json.gz         the files of the PROSITE data directory that the runs read: prosite.lines of the entries that the reports name (and the six pruned names) and the documentation of the entries of
                                                                             the -full runs (PROSITE release 2026_04, SIB Swiss Institute of Bioinformatics, CC BY-NC-ND 4.0; the documentation texts are reproduced unchanged)
The data files are made by prosextract from the release with the checksums of out/environment.txt (this script checks them); the files are read from a directory made by `prosextract` of that release (REAL_PROSITE_DIR=.../PROSITE).
Usage: REAL_PROSITE_DIR=DIR/PROSITE python3 -I convert_patmatmotifs_prosite_returned.py UNZIPPED_ARCHIVE_DIR OUT_DIR"""
import calendar, gzip, hashlib, json, os, re, sys, time

KNOWN_INPUTS = {'one.fa', 'two.fa', 'five.fa', 'many.fa', 'damaged.fa', 'empty.fa', 'rows.txt'}
PRUNED = ['myristyl', 'asn_glycosylation', 'camp_phospho_site', 'pkc_phospho_site', 'ck2_phospho_site', 'tyr_phospho_site']
DATADIR = '/usr/share/EMBOSS/data/PROSITE/'


def read(path):
    return open(path, 'rb').read().decode('latin-1')


def main(src, out):
    prosite = os.environ['REAL_PROSITE_DIR']
    env = read(os.path.join(src, 'out', 'environment.txt'))
    assert 'ad366c80e518d2e919725d5d164652573fd95d52b98a49ea3c50e4317ee4b7a6' in env, 'prosite.dat of another release'
    rows_dir = os.path.join(src, 'out', 'rows')
    rows, motifs, full_motifs = [], set(), set()
    for name in sorted(os.listdir(rows_dir)):
        d = os.path.join(rows_dir, name)
        args = read(os.path.join(d, 'args.txt')).split()
        files = {f: read(os.path.join(d, f)) for f in sorted(os.listdir(d)) if f in KNOWN_INPUTS and f != 'rows.txt'}
        # only the input of the run: the sequence file named by -sequence
        seq = args[args.index('-sequence') + 1] if '-sequence' in args else None
        files = {k: v for k, v in files.items() if k == seq}
        outfiles = {f: read(os.path.join(d, f)) for f in sorted(os.listdir(d)) if f not in KNOWN_INPUTS and f not in ('args.txt', 'status.txt', 'stderr.txt', 'stdout.txt', 'valgrind.txt', 'valgrind_status.txt')}
        rc = int(read(os.path.join(d, 'status.txt')))
        assert int(read(os.path.join(d, 'valgrind_status.txt'))) != 99, 'valgrind reports an error in ' + name        # --error-exitcode=99; any other status is the exit status of the program
        out_text = read(os.path.join(d, 'stdout.txt'))
        err_text = read(os.path.join(d, 'stderr.txt'))
        report = outfiles.get('out.rpt', '') or out_text
        m = re.search(r'# Rundate: (\w+ \d+ \w+ \d+ [\d:]+)', report) or re.search(r'# Rundate: (\w+ \d+ \w+ \d+ [\d:]+)', ''.join(outfiles.values()))
        epoch = calendar.timegm(time.strptime(m.group(1), '%a %d %b %Y %H:%M:%S')) if m else None
        for text in [report] + list(outfiles.values()):
            motifs.update(re.findall(r'^Motif = (\S+)', text, re.M))
            full_motifs.update(re.findall(r'^# Motif: (\S+)', text, re.M))
        row = {'id': len(rows), 'args': args, 'files': files, 'rc': rc, 'out': out_text, 'err': err_text, 'outfiles': outfiles}
        if epoch is not None:
            row['epoch'] = epoch
        rows.append(row)
    lines = read(os.path.join(prosite, 'prosite.lines'))
    blocks = [b for b in lines.split('\n\n') if b.strip()]
    keep, docs = [], {}
    for b in blocks:
        head = b.split('\n')[0].split(' ')
        nm, ac = head[0], head[1]
        if nm in motifs or nm in full_motifs or nm.lower() in PRUNED:
            keep.append(b)
            if nm in full_motifs:
                docs[DATADIR + ac] = read(os.path.join(prosite, ac))
    data = {DATADIR + 'prosite.lines': '\n\n'.join(keep) + '\n\n'}
    data.update(docs)
    with gzip.open(os.path.join(out, 'emboss_prosite_returned_oracle.jsonl.gz'), 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows:
            fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    with gzip.open(os.path.join(out, 'emboss_prosite_release_data.json.gz'), 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        json.dump(data, fh, ensure_ascii=True)
    print(len(rows), 'rows;', len(keep), 'entries;', len(docs), 'documentation files')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
