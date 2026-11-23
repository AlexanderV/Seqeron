#!/usr/bin/env python3
"""B09 WP67 (F306) oracle of the regular expression engine of the library of EMBOSS 6.6.0 (ajRegComp / ajRegExecC = pcre_compile / pcre_exec of the bundled PCRE 7, libepcre.so.7, no options): the expressions that the PROSITE conversion
(embPatPrositeToRegExp) can write (letters, [^BJOUXZ], [^#], classes and braces with any text of a pattern, the counted repeats of x(n,m), < > anchors) and expressions of the PCRE syntax in general, each searched in several subjects.
TSV (gzip): pattern TAB subject TAB result; a text is written with the escapes \\\\ \\t \\n \\r and \\xHH for the other control and non-ASCII characters; the result is `M offset length`, `-` (no match) or `E position message` (the library refuses the expression:
`ajRegComp` returns NULL and the programs that use it die by a segmentation fault).
Usage: python3 scripts/alignment/generate_emboss_pcre_oracle.py OUT.tsv.gz [n=5000] [seed]      (needs /usr/lib/emboss/lib/libepcre.so.7 of the EMBOSS 6.6.0 package)"""
import ctypes, gzip, random, sys

L = ctypes.CDLL('/usr/lib/emboss/lib/libepcre.so.7')
L.epcre_compile.restype = ctypes.c_void_p
L.epcre_compile.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.POINTER(ctypes.c_char_p), ctypes.POINTER(ctypes.c_int), ctypes.c_void_p]
L.pcre_exec.restype = ctypes.c_int
L.pcre_exec.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.c_int]

AA = 'ACDEFGHIKLMNPQRSTVWY'
BODY = list('ABCDEFGHIKLMNPQRSTVWY') + list('abcxz') + ['-', '-', '^', '$', '.', '*', '+', '?', '(', ')', '|', '[', ':', '{', '}', ',', '0', '1', '2', '3', '9', ' ', '#', '=', '!', '<', '>', '_', '/', '@', '~'] + \
    ['\\d', '\\w', '\\s', '\\D', '\\W', '\\S', '\\b', '\\B', '\\x41', '\\x', '\\0', '\\12', '\\1', '\\E', '\\Q', '\\n', '\\t', '\\.', '\\-', '\\]', '\\[', '\\\\', '\\A', '\\z', '\\Z', '\\h', '\\v', '\\R', '\\K', '\\C', '\\e', '\\a', '\\f', '\\cA', '\\c', '\\p', '\\X', '\\N', '\\u', '\\g', '\\8', '\\ ', '[:alpha:]', '[:^digit:]', '[:foo:]', '[.a.]', '[=a=]', '(?:', '(?=', '(?!', '(?>', '(?#', '(?<n>', '(?P<n>', '(?<=', '(?<!', '(?i', '(?i)', '(?z', '(?R', '(?1', '(?|', '(?C', '(?(', "(?'n'", '(?P', '(?<', '\\k<n>', '\\k{n}', '\\g1', '\\g{1}', '\\g{-1}', '\\g<1>', '\\g{0}', 'A-C', 'C-A', 'A-', '0-9']


def body(rnd, n):
    return ''.join(rnd.choice(BODY) for _ in range(n))


def token(rnd):
    r = rnd.random()
    if r < 0.22: return rnd.choice(AA)
    if r < 0.30: return rnd.choice(['[^BJOUXZ]', '[^#]', '[ST]', '[RK]', '[^P]', '[AG]', '[DE]', '[ILVM]'])
    if r < 0.38: return rnd.choice(['{2}', '{1,3}', '{2,}', '{0,2}', '{0}', '{3,3}', '{2,1}', '{70000}', '{,2}', '{a}', '{}', '{ 2}', '{2 }', '{2,a}'])
    if r < 0.50: return '[' + ('^' if rnd.random() < 0.3 else '') + body(rnd, rnd.randint(0, 4)) + ']'
    if r < 0.58: return '{' + body(rnd, rnd.randint(0, 4)) + '}'
    if r < 0.63: return '[^' + body(rnd, rnd.randint(0, 3)) + ']'
    if r < 0.68: return rnd.choice(['^', '$', '^', '$'])
    if r < 0.74: return rnd.choice(['.', '*', '+', '?', '|', '(', ')', '(?:', '(?=', '(?!', '(?>', '*?', '+?', '??', '*+', '++', '?+', '{2}?', '{2}+'])
    if r < 0.82: return rnd.choice(['\\d', '\\w', '\\s', '\\b', '\\B', '\\D', '\\x41', '\\0', '\\1', '\\2', '\\.', '\\A', '\\z', '\\Z', '\\h', '\\Q', '\\E', '\\Qab\\E', '\\n', '\\cA', '\\R', '\\K', '\\G', '\\C', '\\e', '\\y', '\\12', '\\18', '\\8', '\\08', '\\x{41}', '\\x{4}', '\\x{', '\\g1', '\\k', '\\g{2}', '\\k<n>', '(a)\\g1', '(a)\\1', '(?<n>a)\\k<n>', '(a)(b)\\g{-2}'])
    return body(rnd, 1)


def pattern(rnd):
    n = rnd.choice([1, 2, 3, 3, 4, 5, 6, 8])
    t = ''.join(token(rnd) for _ in range(n))
    if rnd.random() < 0.2: t = '^' + t
    if rnd.random() < 0.1: t += '$'
    return t


def subjects(rnd, pat):
    out = ['', ''.join(rnd.choice(AA) for _ in range(rnd.randint(1, 10)))]
    # strings made from the characters of the pattern
    lit = [c for c in pat if c.isalnum() or c in '{}-,#'] or ['A']
    for _ in range(7):
        out.append(''.join(rnd.choice(lit) for _ in range(rnd.randint(1, 12))))
    out.append(''.join(rnd.choice(lit + list(AA[:6])) for _ in range(rnd.randint(1, 12))))
    out.append(rnd.choice(['A{2}', '{a}', 'A{,2}', 'AAA', 'ab\n', 'a-c', '\n', 'A\nB', ' a ', '0', 'a|b']))
    return out


def esc(s):
    o = []
    for ch in s:
        c = ord(ch)
        if ch == '\\': o.append('\\\\')
        elif ch == '\t': o.append('\\t')
        elif ch == '\n': o.append('\\n')
        elif ch == '\r': o.append('\\r')
        elif c < 32 or c > 126: o.append('\\x%02X' % c)
        else: o.append(ch)
    return ''.join(o)


def run(p, s):
    err = ctypes.c_char_p(); pos = ctypes.c_int()
    h = L.epcre_compile(p.encode('latin-1'), 0, ctypes.byref(err), ctypes.byref(pos), None)
    if not h: return 'E %d %s' % (pos.value, err.value.decode())
    ov = (ctypes.c_int * 90)()
    n = L.pcre_exec(h, None, s.encode('latin-1'), len(s), 0, 0, ov, 90)
    if n == -1: return '-'
    if n < 0: return None
    return 'M %d %d' % (ov[0], ov[1] - ov[0])


def main():
    out = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 5000
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 20262301
    rnd = random.Random(seed)
    seen = set()
    rows = []
    while len(seen) < n:
        p = pattern(rnd)
        if p in seen or '\x00' in p: continue
        seen.add(p)
        for s in subjects(rnd, p):
            r = run(p, s)
            if r is None: continue                         # recursion limit and the like: no row
            rows.append('%s\t%s\t%s' % (esc(p), esc(s), r))
            if r[0] == 'E': break                          # one row of an expression the library refuses
    with gzip.open(out, 'wt', encoding='ascii', newline='\n', compresslevel=9) as fh:
        fh.write('\n'.join(rows) + '\n')
    from collections import Counter
    print(len(rows), 'rows;', Counter(r.rsplit('\t', 1)[1][0] for r in rows))


if __name__ == '__main__':
    main()
