#!/usr/bin/env python3
"""The corruption rows of work package S2.4: every damaged container of a corpus (written by the explicit test CorruptionCorpusS24Tests.WriteTheCorpus) is opened and used by the runner Corpus.csproj built against each build under test (the
   S1 head, the S2.4 head), in Checked and in Trusted mode, one process per file and mode, with a watchdog. Prints the counts per build and mode (refused, answered, invalid, other, crash, hang) and the number of files on which the two builds differ.
   usage: run_corpus_s24.py <runner dir of the base> <runner dir of the new build> <corpus dir> <out.tsv> [workers]"""
import collections, concurrent.futures, os, subprocess, sys

base, new, corpus, out = sys.argv[1:5]
workers = int(sys.argv[5]) if len(sys.argv) > 5 else 4
env = dict(os.environ, PATH=os.environ['PATH'] + ':/root/.dotnet', DOTNET_ROOT='/root/.dotnet', DOTNET_CLI_TELEMETRY_OPTOUT='1', DOTNET_NOLOGO='1')
files = sorted(f for f in os.listdir(corpus) if f.endswith('.sesa'))


def run(runner, path, mode):
    try:
        p = subprocess.run(['dotnet', os.path.join(runner, 'Corpus.dll'), path, mode], capture_output=True, text=True, timeout=25, env=env)
    except subprocess.TimeoutExpired:
        return 'hang'
    if p.returncode != 0:
        return 'crash rc=%d' % p.returncode
    line = p.stdout.strip().splitlines()[-1] if p.stdout.strip() else 'no-output'
    return line


def one(args):
    f, mode = args
    path = os.path.join(corpus, f)
    return f, mode, run(base, path, mode), run(new, path, mode)


jobs = [(f, m) for f in files for m in ('Checked', 'Trusted')]
rows = []
with concurrent.futures.ThreadPoolExecutor(workers) as ex:
    for r in ex.map(one, jobs):
        rows.append(r)
with open(out, 'w') as o:
    o.write('file\tmode\tbase\tnew\n')
    for f, m, a, b in rows:
        o.write('%s\t%s\t%s\t%s\n' % (f, m, a, b))


def klass(x):
    return x.split(' ')[0] + (' ' + x.split(' ')[1] if x.startswith('crash') else '')


for mode in ('Checked', 'Trusted'):
    for name, idx in (('base', 2), ('new', 3)):
        c = collections.Counter(klass(r[idx]) for r in rows if r[1] == mode)
        print('%-8s %-5s %s' % (mode, name, dict(sorted(c.items()))))
    diff = [r for r in rows if r[1] == mode and r[2] != r[3]]
    print('%-8s files on which the two builds differ: %d of %d' % (mode, len(diff), len([r for r in rows if r[1] == mode])))
    for r in diff[:10]:
        print('   ', r)
