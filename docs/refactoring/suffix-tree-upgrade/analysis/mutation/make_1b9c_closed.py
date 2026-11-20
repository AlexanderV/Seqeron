#!/usr/bin/env python3
"""Work package 1b.9c: the lists for the second run of every mutant that a new test closes (1b9c_closed.py: CLOSED), against the final code and the final tests ("the kill must reproduce"), one list per source file:
   1b9c_c_<d|e>_<file>_listed.json. The mutants are taken from the lists of 1b.7d and 1b.7e (the same single-site edit); a mutant whose pattern the production seams of this package changed is written again in CLOSED_PATTERNS.
   They are run with MUTANT_FILTER_ONLY set to the new test classes (the kill must come from a test of this package):  run_1b9c_closed.sh.
   usage: make_1b9c_closed.py"""
import glob, json, os, re
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '../../../../..'))
exec(open(os.path.join(HERE, '1b9c_closed.py')).read())
# the mutants whose edit moved with a seam of this package (the same mutation, written against the code as it is now)
CLOSED_PATTERNS = {
    ('e', 'P63'): ('FreeBlock(_cur);\n', ''),
    ('e', 'P64'): ('FreeBlock(_end);\n', ''),
    ('d', 'M18'): ('if (received != end - start) throw new InvalidDataException', 'if (received > end - start) throw new InvalidDataException'),  # the check became CheckWindow
    ('e', 'B124'): ('if (actual != planned) throw', 'if (actual < planned) throw'),  # the check became CheckContainerLength
}
pool = {}
for f in sorted(glob.glob(os.path.join(HERE, '1b7d_*_listed.json'))):
    if re.match(r'1b7d_r\d', os.path.basename(f)):
        continue
    for m in json.load(open(f)):
        pool[('d', m['id'])] = m
for f in sorted(glob.glob(os.path.join(HERE, '1b7e_r[12]_*_listed.json'))):
    for m in json.load(open(f)):
        pool[('e', m['id'])] = m
by_file = {}
for key in sorted(CLOSED):
    m = dict(pool[key])
    if key in CLOSED_PATTERNS:
        m['old'], m['new'] = CLOSED_PATTERNS[key]
    src = open(os.path.join(ROOT, m['file'])).read()
    if src.count(m['old']) != 1:
        print('STALE', key, src.count(m['old']), m['desc'])
    by_file.setdefault((key[0], m['file']), []).append({'id': m['id'], 'desc': m['desc'], 'file': m['file'], 'old': m['old'], 'new': m['new'], 'filter': '', 'closed_by': CLOSED[key]})
for (k, file), ms in by_file.items():
    name = '1b9c_c_%s_%s_listed.json' % (k, os.path.splitext(os.path.basename(file))[0].lower())
    json.dump(ms, open(os.path.join(HERE, name), 'w'), indent=1)
    print(name, [m['id'] for m in ms])
