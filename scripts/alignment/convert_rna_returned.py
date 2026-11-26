#!/usr/bin/env python3
"""B09 WP83 (item A2-19c-rna): converts the returned archive of the Docker package docs/Validation/review-2026-09/handoff/A2-19c-rna (returned/A2-19c-rna-out.zip = run 3 of 2026-11-26, with CONTRAfold 2.02 and DAFS 0.0.3 built; run 2 stays as .run2.zip: the sets of inputs/sets.tsv run twice through `mafft <options> in.fa`
of the EXTENDED build of MAFFT 7.505 (mccaskillwrap, mxscarnamod; installed at /opt/mafft-ext) with ViennaRNA 2.5.1, linux/amd64) into the test oracle tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/mafft_rna_oracle.tsv.gz
(the format of every mafft oracle, see SequenceAligner_MafftProgram_Tests.Load: id <TAB> group <TAB> args (JSON) <TAB> input <TAB> files (JSON) <TAB> rc <TAB> stdout <TAB> sha256(stderr) <TAB> stderr (always: the texts are short enough) <TAB> outfiles (JSON); escaped: \\ \\n \\t \\r).
Groups: mccaskill (rc 0), mccaskill-stop (rc 1: protein input, "nuc only"), rnaalifold (rc 0), rnaalifold-abort (the glibc heap check aborts tbfast on a protein input; only the process number differs between the repetitions: normalised),
scarnapair-stop (rc 1: `tbfast _ ... -s` reads hat4, which only --mccaskill / --contrafold / --dafs write), dafspair-stop (rc 1: `pairlocalalign -G` reads hat4: the same file is missing),
dafs-error (rc 1: `mccaskillwrap -G` runs `dafs --mafft-out ...` and the DAFS of the package, the public release 0.0.3, has no such option: "ERROR IN mccaskill_mea"), contrafold (rc 0: contrafoldwrap + CONTRAfold 2.02 inference + the RNA progressive step).
Only deterministic runs are taken (the two repetitions equal after the normalisation of the process numbers); the rest is listed in the report. The archive is untrusted data: it is parsed, never executed.
Usage: python3 -I scripts/alignment/convert_rna_returned.py <returned.zip>"""
import gzip, hashlib, json, os, re, sys, tempfile, zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', '..', 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment')


def esc(s):
    return s.replace('\\', '\\\\').replace('\n', '\\n').replace('\t', '\\t').replace('\r', '\\r')


def norm_err(t):
    return re.sub(r'(mafft: line \d+:)\s+\d+ (Aborted|Segmentation fault)', r'\1 P \2', t)


def group_of(args, run):
    opts = ' '.join(args)
    err = run['stderr']
    if 'Aborted' in err:
        return 'rnaalifold-abort'
    if '--dafspair' in args and run['rc'] != 0 and 'Cannot open hat4' in err:
        return 'dafspair-stop'
    if '--dafs' in args and run['rc'] != 0 and 'ERROR IN mccaskill_mea' in err:
        return 'dafs-error'
    if '--contrafold' in args and run['rc'] == 0:
        return 'contrafold'
    if run['rc'] != 0 and 'Cannot open hat4' in err:
        return 'scarnapair-stop'
    if run['rc'] != 0 and 'nuc only' in err:
        return 'mccaskill-stop'
    if '--mccaskill' in args:
        return 'mccaskill'
    if '--rnaalifold' in args:
        return 'rnaalifold'
    raise SystemExit('unclassified set: ' + opts)


def main(archive):
    rows, left, counts = [], [], {}
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(archive) as z:
            for n in z.namelist():
                assert not n.startswith('/') and '..' not in n, n
            z.extractall(td)
        out = os.path.join(td, 'out')
        for k, opts, rc, det, secs in json.load(open(os.path.join(out, 'index.json'), encoding='utf-8')):
            rec = json.load(open(os.path.join(out, 'sets', k + '.json'), encoding='utf-8'))
            a, b = rec['run1'], rec['run2'] or rec['run1']
            same = all((norm_err(a[x]) if x == 'stderr' else a[x]) == (norm_err(b[x]) if x == 'stderr' else b[x]) for x in ('rc', 'stdout', 'stderr', 'files'))
            if a['rc'] == 'timeout' or not same:
                left.append((k, opts, 'timeout' if a['rc'] == 'timeout' else 'differs'))
                continue
            args = [x for x in rec['args'] if x != 'in.fa']
            err = norm_err(a['stderr'])
            g = group_of(args, a)
            counts[g] = counts.get(g, 0) + 1
            outfiles = {n: t for n, t in a['files'].items() if n != 'in.fa'}
            cols = [k, g, json.dumps(args), rec['input'], '{}', str(a['rc']), a['stdout'], hashlib.sha256(err.encode('utf-8')).hexdigest(), err, json.dumps(outfiles)]
            rows.append('\t'.join(esc(c) for c in cols))
    with gzip.open(os.path.join(DATA, 'mafft_rna_oracle.tsv.gz'), 'wt', encoding='utf-8', newline='\n') as f:
        f.write('\n'.join(rows) + '\n')
    print('rows', len(rows), counts, 'left out', left, file=sys.stderr)


if __name__ == '__main__':
    main(sys.argv[1])
