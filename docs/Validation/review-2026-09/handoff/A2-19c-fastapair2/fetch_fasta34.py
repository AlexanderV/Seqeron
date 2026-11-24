#!/usr/bin/env python3
"""Builds `fasta34` (FASTA 3.4, release 34.26.5 of 2007-04-29, W. R. Pearson) and installs it into /usr/local/bin. Neutral User-Agent.
Source: inputs/fasta-34.26.5.shar.Z (the file https://fasta.bioch.virginia.edu/wrpearson/fasta/fasta33-35/fasta-34.26.5.shar.Z, sha256 d1611406c868a1f4c02773853a7be74b080426061f2445348ffade12ad7bfd02,
kept in the package), else the same URL. The shell archive is NOT executed: it is parsed here (each `sed 's/^X//' << 'SHAR_EOF' > 'name'` section is written to a file).
The build is `make -f Makefile.linux fasta34` with a gcc line that accepts the 2007 C. The build log and the sha256 of the binary are returned with the results."""
import hashlib, os, re, shutil, subprocess, sys, urllib.request

URL = 'https://fasta.bioch.virginia.edu/wrpearson/fasta/fasta33-35/fasta-34.26.5.shar.Z'
SHA = 'd1611406c868a1f4c02773853a7be74b080426061f2445348ffade12ad7bfd02'

path = os.path.join(sys.argv[1], 'fasta-34.26.5.shar.Z')
if os.path.exists(path):
    data = open(path, 'rb').read()
    print('using', path)
else:
    data = urllib.request.urlopen(urllib.request.Request(URL, headers={'User-Agent': 'Seqeron-review/1.0'}), timeout=120).read()
    print('fetched', URL)
print('sha256 of the archive', hashlib.sha256(data).hexdigest(), '(expected %s)' % SHA)
text = data
if data[:2] in (b'\x1f\x9d', b'\x1f\x8b'):  # compress(1) .Z or gzip: gzip -d reads both
    text = subprocess.run(['gzip', '-dc'], input=data, capture_output=True, check=True).stdout
lines = text.split(b'\n')
pat = re.compile(rb"^sed 's/\^X//' << 'SHAR_EOF' > '([^']+)' &&$")
os.makedirs('src', exist_ok=True)
i = n = 0
while i < len(lines):
    m = pat.match(lines[i])
    if m:
        name = m.group(1).decode()
        assert '..' not in name and not name.startswith('/'), name
        body = []
        i += 1
        while lines[i] != b'SHAR_EOF':
            body.append(lines[i][1:] if lines[i].startswith(b'X') else lines[i])
            i += 1
        open(os.path.join('src', name), 'wb').write(b'\n'.join(body) + (b'\n' if body else b''))
        n += 1
    i += 1
print('extracted', n, 'files')
os.chdir('src')
subprocess.check_call(['make', '-f', 'Makefile.linux', 'CC=gcc -g -O -std=gnu89 -fcommon -fno-strict-aliasing -w', 'fasta34'])
shutil.copy('fasta34', '/usr/local/bin/fasta34')
print('installed fasta34', hashlib.sha256(open('fasta34', 'rb').read()).hexdigest())
