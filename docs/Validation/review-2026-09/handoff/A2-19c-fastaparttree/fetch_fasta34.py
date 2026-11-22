#!/usr/bin/env python3
"""Finds the source of FASTA 3.4 (fasta34), builds it and installs `fasta34` into /usr/local/bin. Neutral User-Agent.
Order: a tarball named fasta34*.tar.gz / fasta3*.tar.gz in the directory given as argument (./inputs), then the upstream archive locations below.
Prints the sha256 of the tarball used and of the built binary (returned with the results)."""
import glob, hashlib, os, shutil, subprocess, sys, tarfile, urllib.request

URLS = [
    'https://fasta.bioch.virginia.edu/wrpearson/fasta/fasta34/fasta-34.26.5.tar.gz',
    'https://faculty.virginia.edu/wrpearson/fasta/fasta34/fasta-34.26.5.tar.gz',
    'ftp://ftp.virginia.edu/pub/fasta/fasta3.tar.gz',
    'https://github.com/wrpearson/fasta34/archive/refs/heads/master.tar.gz',
]


def get(url):
    return urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'Seqeron-review/1.0'}), timeout=120).read()


local = sorted(glob.glob(os.path.join(sys.argv[1], 'fasta3*.tar.gz')) + glob.glob(os.path.join(sys.argv[1], 'fasta-3*.tar.gz')))
data = None
if local:
    data = open(local[0], 'rb').read()
    print('using', local[0])
else:
    for u in URLS:
        try:
            data = get(u)
            print('fetched', u)
            break
        except Exception as e:  # noqa: BLE001
            print('failed', u, e)
if data is None:
    sys.exit('no fasta34 source: place fasta34*.tar.gz into inputs/')
open('fasta34.tar.gz', 'wb').write(data)
print('sha256', hashlib.sha256(data).hexdigest())
with tarfile.open('fasta34.tar.gz') as t:
    t.extractall('.')
src = [d for d in os.listdir('.') if os.path.isdir(d)][0]
os.chdir(os.path.join(src, 'src'))
mk = 'Makefile.linux64_sse2' if os.path.exists('../make/Makefile.linux64_sse2') else None
subprocess.check_call('make -f ../make/Makefile.linux64 all' if os.path.exists('../make/Makefile.linux64') else 'make all', shell=True)
for name in ('fasta34', 'fasta3', 'fasta34_t', 'fasta3_t'):
    for base in ('.', '../bin'):
        p = os.path.join(base, name)
        if os.path.exists(p):
            shutil.copy(p, '/usr/local/bin/fasta34')
            print('installed', p, hashlib.sha256(open(p, 'rb').read()).hexdigest())
            sys.exit(0)
sys.exit('fasta34 binary not found after the build')
