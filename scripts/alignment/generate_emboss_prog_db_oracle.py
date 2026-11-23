#!/usr/bin/env python3
"""B09 WP67 (F309) oracle of the EMBOSS 6.6.0 data resource programs on the standard configuration (emboss.standard, EDAM.obo, DRCAT.dat of the installation): drfinddata, dbtell, showdb
    <program> -auto ARGS ...      every row: arguments, standard output, error stream, exit status, written files (see emboss_prog_oracle_lib.py)
Usage: python3 scripts/alignment/generate_emboss_prog_db_oracle.py PROGRAM OUT.jsonl.gz [n] [seed]      (needs EMBOSS 6.6.0)"""
import os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import emboss_prog_oracle_lib as L  # noqa: E402

DATA = '/usr/share/EMBOSS/'


def standard_dbs():
    return re.findall(r'^DBNAME\s+(\S+)', open(DATA + 'emboss.standard', encoding='latin-1').read(), re.M)


def drcat_ids():
    return [m.group(1).strip() for m in re.finditer(r'^ID\s+(.*)$', open(DATA + 'data/DRCAT.dat', encoding='latin-1').read(), re.M)]


def edam_words():
    ids, words = [], []
    for ln in open(DATA + 'data/EDAM.obo', encoding='latin-1'):
        if ln.startswith('id: '):
            ids.append(ln[4:].strip())
        elif ln.startswith('name: '):
            words += re.findall(r'[A-Za-z0-9]+', ln[6:])
    return ids, sorted(set(words))


IDS, WORDS = edam_words()
DBS = standard_dbs()
DRCAT = drcat_ids()


def case(rnd, s):
    r = rnd.random()
    return s.upper() if r < 0.1 else s.capitalize() if r < 0.2 else s


def gen_dbtell(rnd, i, d):
    r = rnd.random()
    if r < 0.45: name = case(rnd, rnd.choice(DBS))
    elif r < 0.8: name = case(rnd, rnd.choice(DRCAT))
    elif r < 0.85: name = rnd.choice(['', ' ', 'nosuch', 'drcat:x', 'edam*', 'x y', 'taxon ', 'edam_', 'EDAM_DATA', '*'])
    else: name = rnd.choice(DRCAT).lower()[:rnd.randint(1, 5)]
    args = ['-database', name]
    if rnd.random() < 0.3: args.append('-full')
    r = rnd.random()
    if r < 0.4: args += ['-outfile', 'stdout']
    elif r < 0.55: args += ['-outfile', 'out.txt']
    if rnd.random() < 0.03: args += ['-server', '']
    if rnd.random() < 0.03: args = []
    return args, {}


def gen_showdb(rnd, i, d):
    args = []
    flags = ['html', 'protein', 'nucleic', 'sequence', 'feature', 'text', 'taxonomy', 'resource', 'assembly', 'obo', 'xml', 'standard', 'user', 'include', 'full', 'access', 'methods', 'taxscope', 'fields', 'numfields',
             'aliases', 'numaliases', 'examples', 'numexamples', 'defined', 'release', 'only', 'heading', 'type', 'id', 'query', 'all', 'comment']
    p = rnd.choice([0.05, 0.15, 0.3])
    for f in flags:
        if rnd.random() < p:
            args.append(('-' if rnd.random() < 0.55 else '-no') + f)
    units = [[x] for x in args]
    if rnd.random() < 0.15:
        units.append(['-database', rnd.choice(DBS + [x.upper() for x in DBS[:3]] + ['nosuch', 'embl', ''])])
    r = rnd.random()
    if r < 0.4: units.append(['-outfile', 'stdout'])
    elif r < 0.5: units.append(['-outfile', 'out.txt'])
    rnd.shuffle(units)
    return [x for u in units for x in u], {}


def gen_drfinddata(rnd, i, d):
    r = rnd.random()
    if r < 0.3: q = rnd.choice(WORDS)
    elif r < 0.5: q = rnd.choice(IDS).split(':')[-1]
    elif r < 0.6: q = rnd.choice(IDS)
    elif r < 0.7: q = rnd.choice(WORDS)[:rnd.randint(3, 5)] + '*'
    elif r < 0.8: q = ','.join(rnd.choice(WORDS) for _ in range(rnd.randint(2, 3)))
    elif r < 0.9: q = rnd.choice(WORDS).upper()
    else: q = rnd.choice(['', 'nosuchword', 'x y', 'protein|sequence', '?', 'pro?ein', 'edam*', '0001', 'format', 'data'])
    args = ['-query', q]
    if rnd.random() < 0.25: args.append('-nosubclasses')
    if rnd.random() < 0.1: args.append('-subclasses')
    if rnd.random() < 0.1: args.append('-sensitive')
    r = rnd.random()
    if r < 0.5: args += ['-outfile', 'stdout']
    elif r < 0.65: args += ['-outfile', 'out.res']
    if rnd.random() < 0.3: args += ['-oformat', rnd.choice(['drcat', 'basic', 'list', 'DRCAT', 'nosuch', 'text'])]
    if rnd.random() < 0.05: args += ['-odirectory', rnd.choice(['', 'sub', 'nodir'])]
    return args, {}


GENS = {'dbtell': gen_dbtell, 'showdb': gen_showdb, 'drfinddata': gen_drfinddata}

if __name__ == '__main__':
    prog = sys.argv.pop(1)
    L.main(prog, GENS[prog], 300)
