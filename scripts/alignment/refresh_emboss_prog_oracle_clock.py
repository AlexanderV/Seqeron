#!/usr/bin/env python3
"""B09 WP62 (F281): re-runs every row of an emboss_prog_<program>_oracle.jsonl.gz (the program, the arguments and the input files of the row) under the clock shim faketime.c with a fixed date per row ('epoch', UTC) and writes the
file again with the new results and the epoch: the oracles of WP61 were recorded on the day of their generation, so the rows that hold the date of the run (genbank, embl, ... headers) could only pass on that day.
Usage: python3 scripts/alignment/refresh_emboss_prog_oracle_clock.py PROGRAM FILE.jsonl.gz      (needs EMBOSS 6.6.0)"""
import gzip, json, os, random, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import emboss_prog_oracle_lib as L  # noqa: E402


def main(prog, path):
    rows = [json.loads(line) for line in gzip.open(path, 'rt', encoding='utf-8')]
    changed = 0
    with tempfile.TemporaryDirectory() as d:
        work = os.path.join(d, 'w')
        os.mkdir(work)
        shim = L.clock_shim(d)
        for r in rows:
            epoch = random.Random(0x5eed * 1000003 + r['id']).randint(978307200, 2051222400)
            env = dict(os.environ, LD_PRELOAD=shim, FAKE_EPOCH=str(epoch), TZ='UTC')
            rc, o, e, outfiles, hang = L.run(prog, r['args'], r['files'], work, env=env)
            if (rc, hang) != (r['rc'], r['hang']) and not (r['hang'] or r['rc'] < 0):
                print('row', r['id'], 'status changed', r['rc'], '->', rc, file=sys.stderr)
            if (o, e, outfiles) != (r['out'], r['err'], r['outfiles']):
                changed += 1
            r.update(rc=rc if not (r['hang'] or r['rc'] < 0) else r['rc'], out=o if not (r['hang'] or r['rc'] < 0) else r['out'], err=e if not (r['hang'] or r['rc'] < 0) else r['err'], outfiles=outfiles if not (r['hang'] or r['rc'] < 0) else r['outfiles'], epoch=epoch)
    with gzip.open(path, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows:
            fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    print(prog, len(rows), 'rows,', changed, 'with another result')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
