#!/usr/bin/env python3
"""A5-8c3h (B09 F257, WP58) oracle: EMBOSS 6.6.0 `entret` (emboss/entret.c: the full text of every entry the sequence readers of ajseqread.c read, `ajSeqGetEntryS` = the text the readers store with `ajTextinStoreReadline`)
    entret -auto -sequence SPEC [-sformat1 F] [-sbegin1 ...] [-firstonly] -outfile o.txt
SPEC is `in.dat`, `F::in.dat`, a query on the file (`in.dat:id`, `in.dat:wild*`), a list file, a range (the range does not change the text of the entry). The inputs are the texts of the reader oracle in every input format
(generate_seqret_oracle.py: writer outputs of seqret, entries with feature tables, hand-made texts, damaged files).
JSONL.GZ rows: id, fmt, how, text, files, spec, args, epoch, rc, out, err, hang, memerr.
Usage: python3 scripts/alignment/generate_entret_oracle.py OUT.jsonl.gz PER_FORMAT [seed] [format ...]      (needs EMBOSS 6.6.0)
"""
import gzip, json, os, random, shutil, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_seqret_oracle as G  # noqa: E402

FORMATS = [f for f in G.FORMATS if f not in ('ace', 'consed', 'hennig86', 'unknown')]


def main(out, per, seed, formats):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        work = os.path.join(d, 'w'); os.mkdir(work)
        gen = os.path.join(d, 'g'); os.mkdir(gen)
        for fmt in formats:
            made = 0
            while made < per:
                text, how = G.make_input(rnd, gen, fmt)
                if rnd.random() < 0.15: text, how = G.damaged(rnd, text), how + '+damaged'
                if rnd.random() < 0.05: text = text.rstrip('\n')
                spec, args = 'in.dat', []
                r = rnd.random()
                if r < 0.55: args = ['-sformat1', fmt]
                elif r < 0.75: spec = '%s::in.dat' % fmt
                if rnd.random() < 0.12: spec += rnd.choice(['[2:9]', ':a', ':*', ':seq*', ':%d' % rnd.randint(0, 30)])
                if rnd.random() < 0.1: args += ['-firstonly']
                if rnd.random() < 0.05: args += ['-sbegin1', '2']
                files = {'in.dat': text}
                for nme in os.listdir(work):
                    os.remove(os.path.join(work, nme))
                open(os.path.join(work, 'in.dat'), 'w', newline='\n', encoding='latin-1').write(text)
                cmd = ['-auto', '-sequence', spec] + args + ['-outfile', 'o.txt']
                hang = False
                proc = subprocess.Popen(['entret'] + cmd, cwd=work, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                try:
                    _, errb = proc.communicate(timeout=5)
                    rc, err = proc.returncode, errb.decode('latin-1')
                except subprocess.TimeoutExpired:
                    proc.kill(); proc.communicate(); rc, err, hang = -9, '', True
                if err.startswith('Retrieve sequence entries from flatfile databases and files\n'): err = err[len('Retrieve sequence entries from flatfile databases and files\n'):]
                op = os.path.join(work, 'o.txt')
                outtext = open(op, encoding='latin-1', newline='').read() if os.path.exists(op) and not hang else ''
                memerr = False
                if fmt in G.VG_FORMATS and not hang and shutil.which('valgrind'):
                    vgargs = [a.replace('o.txt', 'ovg.txt') for a in cmd]
                    try:
                        vg = subprocess.run(['valgrind', '-q', '--error-exitcode=99', 'entret'] + vgargs, cwd=work, capture_output=True, timeout=120)
                        import re
                        memerr = vg.returncode == 99 or re.search(rb'^==\d+== (Invalid|Conditional|Use of|Syscall|Mismatched|Source and destination|Process terminating)', vg.stderr, re.M) is not None
                    except subprocess.TimeoutExpired:
                        memerr = True
                rows.append({'id': len(rows), 'fmt': fmt, 'memerr': memerr, 'how': how, 'files': files, 'spec': spec, 'args': args, 'rc': rc, 'out': outtext, 'err': err, 'hang': hang})
                made += 1
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows:
            fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    from collections import Counter
    print(len(rows), 'rows;', sorted(Counter((r['fmt'], 'hang' if r['hang'] else 'died' if r['rc'] else 'out' if r['out'] else 'none') for r in rows).items())[:12])


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 20265803, sys.argv[4:] or FORMATS)
