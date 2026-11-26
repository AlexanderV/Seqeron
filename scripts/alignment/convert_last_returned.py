#!/usr/bin/env python3
"""B09 WP77 (item A2-19c-last): converts the returned archive of the Docker package docs/Validation/review-2026-09/handoff/A2-19c-last (returned/A2-19c-last-out.zip: the sets of inputs/sets.tsv run through `mafft --lastpair` / `--lastmultipair`
[--addfragments | --addfull | --addlong | --addtotop new.fa] of MAFFT 7.505-1 with the real LAST programs behind logging wrappers) into two test oracles in tests/Seqeron/Seqeron.Genomics.Last.Tests/TestData/Alignment (the test project of the GPL assembly that holds LAST):
  mafft_last_oracle.tsv.gz   the format of every mafft oracle (SequenceAligner_MafftProgram_Tests.Load): id <TAB> group <TAB> args (JSON) <TAB> input <TAB> files (JSON: new.fa) <TAB> rc <TAB> stdout <TAB> sha256(stderr) <TAB> stderr (only when short,
                             every 25th) <TAB> outfiles (JSON). Groups: last-plain (--lastpair), last-plain-thread (the thread number of the progress lines is a race: normalised to T), last-crash (the segmentation fault of pairlocalalign -R), last-add (with --addfragments / --addfull / --addlong / --addtotop), last-protein (the stop "Not yet supported"), last-stop (the stops of the script).
                             Only deterministic runs (the two repetitions of the package equal) with a return code are taken; a set that timed out or differs between its two runs is listed in the report and left out.
  last_calls_oracle.tsv.gz   one row per call of lastdb / lastal: set id <TAB> index <TAB> program <TAB> arguments (JSON list) <TAB> input files (JSON: name -> text; _dbd / _db<k>d / _q<k> / _scoringmatrixforlast) <TAB> stdout (the raw output of lastal: the
                             `-f 0` tabular format that readlastresx / readlastresx_group read) <TAB> stderr <TAB> rc
(escaped: \\ \\n \\t \\r in every column.) The archive is untrusted data: it is parsed, never executed.
Usage: python3 -I scripts/alignment/convert_last_returned.py <returned.zip> [--report]"""
import gzip, hashlib, io, json, os, re, sys, tempfile, zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', '..', 'tests', 'Seqeron', 'Seqeron.Genomics.Last.Tests', 'TestData', 'Alignment')


def esc(s):
    return s.replace('\\', '\\\\').replace('\n', '\\n').replace('\t', '\\t').replace('\r', '\\r')


def norm_err(t):
    """The random names of the temporary directories and the process number of a crash are not part of the program's behaviour (run 1 of the package called every set non-deterministic because of them)."""
    t = re.sub(r'mafft\.[A-Za-z0-9]{6,12}', 'mafft.X', t)
    t = re.sub(r'/tmp/tmp[a-z0-9_]{6,10}', '/tmp/tmpX', t)
    return re.sub(r'(mafft: line \d+:)\s+\d+ (Segmentation fault|Aborted)', r'\1 P \2', t)


def norm_thread(t):
    """The number of the thread that prints a progress line is a race of the program (run 2: sets 15, 30, 63, 73, 81, 129, 191 differ only there and in the process numbers of the crashes)."""
    return re.sub(r'\((thread|threads) +\d+', r'(\1 T', re.sub(r'by thread +\d+', 'by thread T', t))


def group_of(rec, run):
    if 'Segmentation fault' in run['stderr'] or 'Aborted' in run['stderr']:
        return 'last-crash'          # pairlocalalign -R overflows its heap (valgrind in out/valgrind): the outcome is the glibc's, kept apart

    opts = rec['options'].split()
    if '--lastmultipair' in opts and rec['mode'] == 'plain':
        return 'last-stop'
    if run['rc'] != 0 and 'Not yet supported' in run['stderr']:
        return 'last-protein'
    if rec['mode'] == 'add':
        return 'last-add'
    return 'last-plain'


def main(archive, report=False):
    oracle, calls, left_out = [], [], []
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(archive) as z:
            for n in z.namelist():
                assert not n.startswith('/') and '..' not in n, n
            z.extractall(td)
        out = os.path.join(td, 'out')
        full = 0
        for line in open(os.path.join(out, 'index.tsv'), encoding='utf-8'):
            k = line.split('\t')[0]
            rec = json.load(gzip.open(os.path.join(out, 'sets', k + '.json.gz'), 'rt', encoding='utf-8'))
            run = rec['run']
            run2 = rec.get('run2') or run
            thread_only = False
            if rec['deterministic'] is False and run['rc'] != 'timeout' and run2['rc'] != 'timeout':
                same = (run['rc'], run['stdout'], run['files']) == (run2['rc'], run2['stdout'], run2['files'])
                if same and norm_thread(norm_err(run['stderr'])) == norm_thread(norm_err(run2['stderr'])):
                    thread_only = True
                    rec['deterministic'] = True
            if '--thread' in rec['options'].split() and norm_thread(norm_err(run['stderr'])) != norm_err(run['stderr']):
                thread_only = True   # which thread takes which pair is a race: it is not part of the behaviour even when both repetitions agree
            if rec['deterministic'] is False or run['rc'] == 'timeout':
                left_out.append((k, 'nondeterministic' if rec['deterministic'] is False else 'timeout'))
                continue
            args = rec['options'].split()
            files = {'new.fa': rec['new.fa']} if rec['new.fa'] else {}
            err = norm_err(run['stderr'])
            if thread_only:
                err = norm_thread(err)
            sha = hashlib.sha256(err.encode('utf-8')).hexdigest()
            keep = len(err) < 200000 and (full % 25 == 0 or len(err) < 400)
            full += 1
            outfiles = {n: t for n, t in run['files'].items() if n not in ('in.fa', 'new.fa')}
            cols = [k, group_of(rec, run), json.dumps(args), rec['in.fa'], json.dumps(files), str(run['rc']), run['stdout'], sha, err if keep else '', json.dumps(outfiles)]
            if thread_only:
                cols[1] += '-thread' if cols[1] != 'last-crash' else ''
            oracle.append('\t'.join(esc(c) for c in cols))
            for i, c in enumerate(run['calls']):
                it = c['items']
                argl = it.get('args', '').split('\n')[:-1] if it.get('args') else []
                # only the files that the call names (the package logs every file of the directory: the same databases again for every call)
                inputs = {n[3:]: t for n, t in it.items() if n.startswith('in.') and n[3:] in argl}
                cols = [k, str(i), c['prog'], json.dumps(it.get('args', '').split('\n')[:-1] if it.get('args') else []), json.dumps(inputs), it.get('stdout', ''), it.get('stderr', ''), it.get('rc', '').strip()]
                calls.append('\t'.join(esc(x) for x in cols))
    for name, rows in (('mafft_last_oracle.tsv.gz', oracle), ('last_calls_oracle.tsv.gz', calls)):
        with gzip.open(os.path.join(DATA, name), 'wt', encoding='utf-8', newline='\n') as f:
            f.write('\n'.join(rows) + '\n')
    print('mafft runs', len(oracle), 'lastdb / lastal calls', len(calls), 'left out', len(left_out), left_out[:10], file=sys.stderr)


if __name__ == '__main__':
    main(sys.argv[1], '--report' in sys.argv)
