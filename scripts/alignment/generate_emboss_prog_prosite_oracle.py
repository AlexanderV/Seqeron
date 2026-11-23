#!/usr/bin/env python3
"""B09 WP67 (F306) oracle of the PROSITE programs of EMBOSS 6.6.0: patmatmotifs (the PROSITE database prosite.lines / documentation files of prosextract searched with a protein sequence), patmatdb (one PROSITE pattern against
protein sequences) and prosextract (prosite.dat / prosite.doc -> PROSITE/prosite.lines + one documentation file per accession):
    <program> -auto ARGS ...      every row: arguments, input files, standard output, error stream, exit status, written files (see emboss_prog_oracle_lib.py)
The data of a row are files of the row: the patterns come from the PATTERN entries of a PROSITE release (prosite.dat of ftp.expasy.org: the lines ID / AC / DE / PA of the entries are used, the documentation texts are made up here), the
program runs without the PROSITE data of the installation (a private mount over its directory, see emboss_prog_oracle_lib.launcher) and with the data directory WORK/edata/ (EMBOSS_DATA) when the row says so.
Usage: PROSITE_DAT=prosite.dat python3 scripts/alignment/generate_emboss_prog_prosite_oracle.py PROGRAM OUT.jsonl.gz [n=300] [seed]      (needs EMBOSS 6.6.0, the file prosite.dat, unshare)"""
import os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import emboss_prog_oracle_lib as L  # noqa: E402
import generate_seqret_oracle as S  # noqa: E402
import generate_emboss_prog_seq_oracle as G  # noqa: E402
import generate_emboss_prog_wp63_oracle as W  # noqa: E402

AA = 'ACDEFGHIKLMNPQRSTVWY'
PRUNED = ['MYRISTYL', 'ASN_GLYCOSYLATION', 'CAMP_PHOSPHO_SITE', 'PKC_PHOSPHO_SITE', 'CK2_PHOSPHO_SITE', 'TYR_PHOSPHO_SITE']
_ENTRIES = None


def load_entries():
    """the PATTERN entries of prosite.dat: dicts id, ac, de (the lines), pa (the pattern lines, without the trailing dot)"""
    global _ENTRIES
    if _ENTRIES is not None:
        return _ENTRIES
    text = open(os.environ.get('PROSITE_DAT', 'prosite.dat'), encoding='latin-1').read()
    out = []
    for block in text.split('\n//\n'):
        lines = block.split('\n')
        idl = [l for l in lines if l.startswith('ID')]
        if not idl or not idl[0].rstrip().endswith('PATTERN.'):
            continue
        e = {'id': idl[0].split()[1].rstrip(';'), 'ac': '', 'de': [], 'pa': []}
        for l in lines:
            if l.startswith('AC'): e['ac'] = l.split()[1].rstrip(';')
            elif l.startswith('DE'): e['de'].append(l[5:].rstrip())
            elif l.startswith('PA'): e['pa'].append(l[5:].rstrip())
        if e['pa']:
            out.append(e)
    _ENTRIES = out
    return out


def pattern_text(e):
    return ''.join(e['pa']).rstrip('.')


def prosite_to_regexp(s):
    """embPatPrositeToRegExp (nucleus/embpat.c)"""
    c = ' '.join(s.upper().split()) if False else s.upper().strip()
    t = ''
    isnt = isct = False
    if '<' in c: isnt = True
    if '>' in c: isct = True
    if isnt: t += '^'
    i = 0
    while i < len(c):
        ch = c[i]
        if ch == '?': t += '[^#]'; i += 1; continue
        if ch == 'X': t += '[^BJOUXZ]'; i += 1; continue
        if ch in '([{':
            close = {'(': ')', '[': ']', '{': '}'}[ch]
            j = c.find(close, i + 1)
            if j < 0: raise ValueError('unmatched')
            body = ''.join(x for x in c[i + 1:j] if x not in '<>')
            t += {'(': '{' + body + '}', '[': '[' + body + ']', '{': '[^' + body + ']'}[ch]
            i = j + 1
            continue
        if ch in AA: t += ch; i += 1; continue
        if ch in ' -<>': i += 1; continue
        raise ValueError('char')
    if isct: t += '$'
    return t


def lines_block(e):
    """the entry as prosextract writes it into PROSITE/prosite.lines"""
    de = e['de'][0].split('.')[0].split()[0] if e['de'] else ''
    pa = pattern_text(e)
    return '%s %s\n %s\n %s\n^%s\n\n' % (e['id'], e['ac'], de, pa, prosite_to_regexp(pa))


def instance(rnd, pa):
    """a string that matches the PROSITE pattern pa"""
    toks = pa.split('-')
    out = []
    for t in toks:
        rep = 1
        m = re.match(r'^(.*?)\((\d+)(?:,(\d+))?\)$', t)
        if m:
            t = m.group(1)
            rep = rnd.randint(int(m.group(2)), int(m.group(3) or m.group(2)))
        t = t.strip('<>')
        if t.lower() == 'x' or t == '':
            pool = AA
        elif t.startswith('['):
            pool = [c for c in t[1:-1] if c in AA] or list(AA)
        elif t.startswith('{'):
            pool = [c for c in AA if c not in t[1:-1]]
        else:
            pool = [c for c in t if c in AA] or list(AA)
        out += [rnd.choice(pool) for _ in range(rep)]
    return ''.join(out)


def doc_text(rnd, e):
    """a made-up documentation text of an entry (the real texts of prosite.doc are copyrighted and are not used)"""
    lines = ['Documentation of %s (%s), a made-up text for the test.' % (e['ac'], e['id']), '']
    for _ in range(rnd.choice([0, 1, 3, 6])):
        lines.append(' '.join(rnd.choice(S.WORDS) for _ in range(rnd.randint(3, 14))) + rnd.choice(['', '', ' ', '\t']))
        if rnd.random() < 0.3: lines.append('')
    lines.append('-Consensus pattern: ' + pattern_text(e))
    return '\n'.join(lines) + '\n'


def choose(rnd, k):
    ents = load_entries()
    pruned = [e for e in ents if e['id'] in PRUNED or e['id'].startswith(('TYR_PHOSPHO_SITE', 'MYRISTYL'))]
    short = [e for e in ents if len(pattern_text(e)) < 40]
    pick = []
    if rnd.random() < 0.5: pick += rnd.sample(pruned, min(len(pruned), rnd.randint(1, 4)))
    pick += rnd.sample(short, min(len(short), k))
    if rnd.random() < 0.3: pick += rnd.sample(ents, rnd.randint(1, 3))
    rnd.shuffle(pick)
    seen, out = set(), []
    for e in pick:
        if e['id'] not in seen:
            seen.add(e['id']); out.append(e)
    return out


def weird_lines(rnd, text):
    """a hand-edited prosite.lines: CR LF, trailing blanks, missing accession, extra blank lines, a line of another shape"""
    out = []
    for ln in text.split('\n'):
        r = rnd.random()
        if r < 0.04: ln += ' '
        elif r < 0.07: ln += '\r'
        elif r < 0.12 and ln and ln[0] not in ' ^': ln = ln.split(' ')[0]
        out.append(ln)
        if rnd.random() < 0.02: out.append('')
    return '\n'.join(out)


def seq_with_sites(rnd, entries, n):
    ents = entries or []
    parts = [''.join(rnd.choice(AA) for _ in range(rnd.randint(2, 60)))]
    for _ in range(rnd.choice([0, 1, 2, 3, 5]) if ents else 0):
        parts.append(instance(rnd, pattern_text(rnd.choice(ents))))
        parts.append(''.join(rnd.choice(AA) for _ in range(rnd.randint(1, 40))))
    s = ''.join(parts)
    if rnd.random() < 0.15: s = s.lower()
    if rnd.random() < 0.05: s = s[:3]
    return s


def gen_patmatmotifs(rnd, i, d):
    ents = choose(rnd, rnd.choice([2, 4, 8, 15, 25]))
    text = ''.join(lines_block(e) for e in ents)
    if rnd.random() < 0.1: text = weird_lines(rnd, text)
    noacc = rnd.random() < 0.04                                       # name lines without an accession: -full opens the documentation file named by the directory PROSITE/
    if noacc: text = '\n'.join(ln.split(' ')[0] if ln and ln[0] not in ' ^' else ln for ln in text.split('\n'))
    files = {}
    in_edata = rnd.random() < 0.3
    prefix = 'edata/PROSITE/' if in_edata else 'PROSITE/'
    r = rnd.random()
    if r < 0.93: files[prefix + 'prosite.lines'] = text
    elif r < 0.97: files['PROSITE/prosite.lines'] = ''                       # an empty data file
    # r >= 0.97: no data file at all
    if rnd.random() < 0.85:
        for e in ents:
            if rnd.random() < 0.85: files[prefix + e['ac']] = doc_text(rnd, e)
    n = rnd.choice([1, 1, 2, 3])
    items = []
    for k in range(n):
        items.append((S.rand_name(rnd) + ('' if k == 0 else str(k)), S.rand_desc(rnd), seq_with_sites(rnd, ents, 100)))
    sfmt = rnd.choice(['fasta'] * 6 + ['gcg', 'embl', 'swiss', 'pir', 'raw', 'genbank', 'ncbi'])
    stext, fmt = G.seq_file(rnd, d, items, sfmt)
    files['in.dat'] = stext
    args = ['-sequence', 'in.dat'] + G.seq_quals(rnd, fmt, '', 0.1)
    if noacc or rnd.random() < 0.4: args.append('-full')
    if rnd.random() < 0.4: args.append('-noprune')
    if rnd.random() < 0.03: args.append(rnd.choice(['-full=x', '-prune=maybe', '-full', '-nofull']))
    args += W.report_quals(rnd, 0.35)
    args += W.report_file(rnd, 'patmatmotifs')
    if rnd.random() < 0.04: args += ['-sequence', 'in.dat']
    if 'sub' in ' '.join(args): files['sub/keep'] = 'x'
    return args, files, {'edata': True}


def fix_prosextract_input(rnd, ents, weird):
    """prosite.dat and prosite.doc of a few entries"""
    dat, doc = [], []
    for k, e in enumerate(ents):
        if weird and rnd.random() < 0.15:
            dat.append('ID   %s; MATRIX.\nAC   %s;\nDE   %s\n//' % (e['id'], e['ac'], ' '.join(e['de'])))
            continue
        lines = ['ID   %s; PATTERN.' % e['id'], 'AC   %s;' % e['ac'], 'DT   01-APR-1990 CREATED; 01-APR-1990 DATA UPDATE; 01-APR-1990 INFO UPDATE.']
        lines += ['DE   ' + l for l in e['de']] if e['de'] else []
        pa = e['pa']
        lines += ['PA   ' + l for l in pa]
        lines += ['CC   /SITE=1,phosphorylation;', 'CC   /VERSION=1;', 'DO   PDOC%05d;' % (k + 1)]
        dat.append('\n'.join(lines) + '\n//')
        body = doc_text(rnd, e).rstrip('\n')
        if weird and rnd.random() < 0.15 and k + 1 < len(ents):
            doc.append('{PDOC%05d}\n{%s; %s}\n{%s; %s}\n{BEGIN}\n%s\n{END}' % (k + 1, e['ac'], e['id'], ents[k + 1]['ac'], ents[k + 1]['id'], body))
        else:
            doc.append('{PDOC%05d}\n{%s; %s}\n{BEGIN}\n%s\n{END}' % (k + 1, e['ac'], e['id'], body))
    sep = rnd.choice(['\n', '\n', '\n\n'])
    return 'CC   made-up header\n//\n' + '\n'.join(dat) + '\n', '{PDOC00000}\n{BEGIN}\nmade-up header of the documentation\n{END}\n' + sep.join(doc) + '\n'


def gen_prosextract(rnd, i, d):
    ents = choose(rnd, rnd.choice([1, 2, 4, 8, 15]))
    weird = rnd.random() < 0.35
    dat, doc = fix_prosextract_input(rnd, ents, weird)
    if weird and rnd.random() < 0.3: dat = weird_lines(rnd, dat)
    if weird and rnd.random() < 0.3: doc = weird_lines(rnd, doc)
    sub = rnd.choice(['', '', 'sub/', 'sub'])
    files = {sub.rstrip('/') + ('/' if sub else '') + 'prosite.dat': dat, sub.rstrip('/') + ('/' if sub else '') + 'prosite.doc': doc}
    if rnd.random() < 0.05: del files[next(iter(files))]
    elif rnd.random() < 0.05: del files[list(files)[-1]]
    args = ['-prositedir', rnd.choice(['./', '.', ''] if not sub else [sub, 'sub', 'sub/'])]
    if rnd.random() < 0.05: args = ['-prositedir', rnd.choice(['nodir/', 'nodir'])]
    if rnd.random() < 0.03: args = []
    return args, files, {'edata': True}


def pcre_tail(rnd):
    """pieces that mean something to the regular expression library"""
    return rnd.choice(['(2)', '(1,3)', '(0)', '(0,2)', '(2,)', '(,2)', '(a)', '(2,1)', '()', '(70000)', '[A-C]', '[C-A]', '[^A]', '[\\d]', '[\\w]', '[\\n]', '[\\x41]', '[a-c]', '[A-]', '[-A]', '[A\\]', '[[:alpha:]', '[A[]', '[^]', '[]', '{A}', '{AB}', '{}', '?', '[AB](2)(3)', 'x(2)(3)', '<', '>', '[<A]', '(<2)', 'xx', 'x-', '-x', '--', ' ', '  A ', 'BJO', 'U', 'Z', 'ab', '*', '+', '.', '|', '\\', '[\\p{L}]', '[\\s]', '[\\.]', '[#]', '[^#]'])


def gen_patmatdb(rnd, i, d):
    ents = load_entries()
    e = rnd.choice(ents)
    pa = pattern_text(e)
    r = rnd.random()
    if r < 0.5: motif = pa
    elif r < 0.6: motif = pa.lower()
    elif r < 0.75: motif = rnd.choice(['A-x(2)-C', '[ST]-x-[RK]', 'N-{P}-[ST]-{P}', '<M', 'K>', '<M-x(1,3)-[ST]>', 'C-x(2,4)-C', 'G-x(2)-G', 'A', 'AC', 'x-A', 'L-x(6)-L', 'x(3)', '[AC]', '{A}'])
    elif r < 0.85: motif = pa + '-' + pcre_tail(rnd)
    elif r < 0.93: motif = rnd.choice([pcre_tail(rnd), pcre_tail(rnd) + rnd.choice('ACDEFG'), rnd.choice('ACDEFG') + pcre_tail(rnd), 'A' + pcre_tail(rnd) + pcre_tail(rnd)])
    else: motif = rnd.choice(['A-x(2', '[AC', '{A', 'A-1', 'A-B-J', 'A@', 'A,C', 'xyz', '', 'A', '(2)', '<>', '-', '[<>]', '(<>)'])
    n = rnd.choice([1, 1, 2, 3, 4])
    items = [(S.rand_name(rnd) + ('' if k == 0 else str(k)), S.rand_desc(rnd), seq_with_sites(rnd, [e], 100)) for k in range(n)]
    text, fmt = G.seq_file(rnd, d, items, rnd.choice(['fasta'] * 6 + ['gcg', 'embl', 'swiss', 'pir', 'raw', 'genbank']))
    files = {'in.dat': text}
    args = ['-sequence', 'in.dat', '-motif', motif] + G.seq_quals(rnd, fmt, '', 0.1)
    args += W.report_quals(rnd, 0.35)
    args += W.report_file(rnd, 'patmatdb')
    if 'sub' in ' '.join(args): files['sub/keep'] = 'x'
    return args, files, ({'vg': True} if rnd.random() < 0.15 else {})


GENS = {'patmatmotifs': gen_patmatmotifs, 'patmatdb': gen_patmatdb, 'prosextract': gen_prosextract}

if __name__ == '__main__':
    prog = sys.argv.pop(1)
    L.main(prog, GENS[prog], 300)
