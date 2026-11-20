#!/usr/bin/env python3
"""A4-3 (B09 F153): which characters the EMBOSS 6.6.0 programs drop from an input sequence.

For every ported program that reads plain sequences (needle water stretcher matcher merger megamerger seqmatchall fuzznuc fuzzpro fuzztran etandem equicktandem dotmatcher palindrome) and for the two
alignment readers (consambig showalign) the script inserts each character 1..255 (except newline, carriage return and '>') into the middle of a sequence of the program, runs the installed binary and
classifies the character by comparing the program's output with (a) the output without the character ('drop'), (b) with N / X in its place ('toN'), (c) with '-' (alignment readers: 'gap').
Result recorded in docs/Validation/review-2026-09/B09.md (F153): for the plain programs 199 characters are 'drop' for every program, sequence and type (nucleotide, protein; first and second sequence), i.e. everything
except the 52 ASCII letters and '*'; the alignment readers (type gapany) turn `. ~ -` into gaps and read `? # +` as N / X, every other non-letter is dropped.
Usage: python3 scripts/alignment/measure_emboss_reader_chars.py [program ...]      (needs the EMBOSS 6.6.0 programs on PATH; takes a few minutes)
"""
import subprocess, re, os, sys, tempfile, collections
D = tempfile.mkdtemp()
def w(name, text): open(os.path.join(D,name),'w').write(text)
def sh(args):
    p = subprocess.run(args, cwd=D, capture_output=True, text=True, errors='replace')
    return p
def norm(s):
    out=[]
    for l in s.splitlines():
        if re.search(r'Rundate|Commandline|Program:|Report_file|# Align_format|Date|\d{4}$', l) and not l.startswith('Seq'): 
            if re.search(r'Rundate|Commandline|Date', l): continue
        out.append(l)
    return '\n'.join(out)
def two(prog, extra=()):
    def f(a, b):
        w('a.fa','>a\n'+a+'\n'); w('b.fa','>b\n'+b+'\n')
        p = sh([prog,'-asequence','a.fa','-bsequence','b.fa','-auto','-stdout','-aformat3','pair']+list(extra))
        return norm(p.stdout) + '|ERR:' + (p.stderr.strip()[-80:] if p.returncode else '')
    return f
def matcher(a,b):
    w('a.fa','>a\n'+a+'\n'); w('b.fa','>b\n'+b+'\n')
    p = sh(['matcher','a.fa','b.fa','-auto','-stdout','-aformat3','pair']); return norm(p.stdout)+'|ERR:'+(p.stderr.strip()[-80:] if p.returncode else '')
def merger(a,b):
    w('a.fa','>a\n'+a+'\n'); w('b.fa','>b\n'+b+'\n')
    p = sh(['merger','-asequence','a.fa','-bsequence','b.fa','-outfile','o.txt','-outseq','o.fa','-auto'])
    o = open(os.path.join(D,'o.fa')).read() if os.path.exists(os.path.join(D,'o.fa')) else ''
    return norm(o) + '|ERR:'+(p.stderr.strip()[-80:] if p.returncode else '')
def megam(a,b):
    w('a.fa','>a\n'+a+'\n'); w('b.fa','>b\n'+b+'\n')
    p = sh(['megamerger','-asequence','a.fa','-bsequence','b.fa','-wordsize','4','-outfile','o.txt','-outseq','o.fa','-auto'])
    o = open(os.path.join(D,'o.fa')).read() if os.path.exists(os.path.join(D,'o.fa')) else ''
    r = open(os.path.join(D,'o.txt')).read() if os.path.exists(os.path.join(D,'o.txt')) else ''
    return norm(o)+norm(r) + '|ERR:'+(p.stderr.strip()[-80:] if p.returncode else '')
def seqmatchall(a,b):
    w('a.fa','>a\n'+a+'\n>b\n'+b+'\n')
    p = sh(['seqmatchall','a.fa','-wordsize','3','-auto','-stdout'])
    return norm(p.stdout)+'|ERR:'+(p.stderr.strip()[-80:] if p.returncode else '')
def fuzz(prog, pat, extra=()):
    def f(a,b=None):
        w('a.fa','>a\n'+a+'\n')
        p = sh([prog,'-sequence','a.fa','-pattern',pat,'-auto','-stdout','-rformat','gff']+list(extra))
        return norm(p.stdout)+'|ERR:'+(p.stderr.strip()[-80:] if p.returncode else '')
    return f
def tandem(prog, extra=()):
    def f(a,b=None):
        w('a.fa','>a\n'+a+'\n')
        p = sh([prog,'-sequence','a.fa','-outfile','r.txt','-origfile','o.txt','-auto']+list(extra))
        r = open(os.path.join(D,'r.txt')).read() if os.path.exists(os.path.join(D,'r.txt')) else ''
        o = open(os.path.join(D,'o.txt')).read() if os.path.exists(os.path.join(D,'o.txt')) else ''
        return norm(r)+norm(o)+'|ERR:'+(p.stderr.strip()[-80:] if p.returncode else '')
    return f
def dotm(a,b):
    w('a.fa','>a\n'+a+'\n'); w('b.fa','>b\n'+b+'\n')
    for f in os.listdir(D):
        if f.startswith('dm'): os.remove(os.path.join(D,f))
    p = sh(['dotmatcher','a.fa','b.fa','-windowsize','3','-threshold','10','-graph','data','-goutfile','dm','-auto'])
    outs=''
    for f in sorted(os.listdir(D)):
        if f.startswith('dm'): outs+=open(os.path.join(D,f)).read()
    return norm(outs)+'|ERR:'+(p.stderr.strip()[-80:] if p.returncode else '')
def pal(a,b=None):
    w('a.fa','>a\n'+a+'\n')
    p = sh(['palindrome','-sequence','a.fa','-minpallen','3','-maxpallen','20','-gaplimit','0','-nummismatches','0','-overlap','-outfile','r.txt','-auto'])
    r = open(os.path.join(D,'r.txt')).read() if os.path.exists(os.path.join(D,'r.txt')) else ''
    return norm(r)+'|ERR:'+(p.stderr.strip()[-80:] if p.returncode else '')
# name: (func, base1 nuc, base2 nuc, base1 prot, base2 prot, insertion idx offset)
NUC1='ACGTTAGCATGCAACGTTAGGCTA'; NUC2='ACGTTAGCATGCTACGTTAGGCTA'
PRO1='ACDEFGHIKLMNPQRSTVWYACDEFGHIK'; PRO2='ACDEFGHIKLMNPSRSTVWYACDEFGHIK'
progs = {
 'needle': (two('needle'),'both'),
 'water': (two('water'),'both'),
 'stretcher': (two('stretcher'),'both'),
 'matcher': (matcher,'both'),
 'merger': (merger,'both'),
 'megamerger': (megam,'nuc'),
 'seqmatchall': (seqmatchall,'both'),
 'fuzznuc': (fuzz('fuzznuc','ACG'),'nuc'),
 'fuzzpro': (fuzz('fuzzpro','DEF'),'prot'),
 'fuzztran': (fuzz('fuzztran','D',['-frame','1']),'nuc'),
 'etandem': (tandem('etandem',['-minrepeat','2','-maxrepeat','4','-threshold','8']),'nuc'),
 'equicktandem': (tandem('equicktandem',['-maxrepeat','4','-threshold','8']),'nuc'),
 'dotmatcher': (dotm,'both'),
 'palindrome': (pal,'nuc'),
}
def plain_sweep(sel):
    for name in sel:
        f, kind = progs[name]
        for which in (['nuc','prot'] if kind=='both' else [kind]):
            b1,b2 = (NUC1,NUC2) if which=='nuc' else (PRO1,PRO2)
            if name in ('etandem','equicktandem'): b1='ACACACACACGTGTGTGTGTAAGCTAGCTAGCTAGC'; b2=b1
            if name=='palindrome': b1='AAGCTTGGATCCAAGCTTAGGCCTAATTCCGGAATTGC'; b2=b1
            for seqno in ([1,2] if name in ('needle','water','stretcher','matcher','merger','megamerger','seqmatchall','dotmatcher') else [1]):
                ref = f(b1,b2)
                groups = collections.defaultdict(list)
                for o in range(1,256):
                    c=chr(o)
                    if c in '\n\r>': continue
                    mid = len(b1)//2
                    d1 = b1[:mid]+c+b1[mid:]; d2 = b2[:mid]+c+b2[mid:]
                    r = f(d1,b2) if seqno==1 else f(b1,d2)
                    if r==ref: cls='drop'
                    else:
                        alt = ('N' if which=='nuc' else 'X')
                        a1 = b1[:mid]+alt+b1[mid:]; a2=b2[:mid]+alt+b2[mid:]
                        r2 = f(a1,b2) if seqno==1 else f(b1,a2)
                        cls = 'toN' if r==r2 else ('ERR' if 'ERR:Error' in r or r.endswith('ERR:') is False and 'ERR:' in r and r.split('|ERR:')[1] else 'other')
                    groups[cls].append(o)
                print(name, which, 'seq', seqno, {k:''.join(chr(x) if 33<=x<127 else '\\x%02x'%x for x in v) if len(v)<80 else '%d chars'%len(v) for k,v in groups.items()})

def cons(rows):
    w('a.fa',''.join('>r%d\n%s\n'%(i,r) for i,r in enumerate(rows)))
    p = sh(['consambig','-sequence','a.fa','-auto','-outseq','o.fa'])
    o = open(os.path.join(D,'o.fa')).read() if os.path.exists(os.path.join(D,'o.fa')) else ''
    return norm(o)+'|ERR:'+(p.stderr.strip()[-80:] if p.returncode else '')
def show(rows):
    w('a.fa',''.join('>r%d\n%s\n'%(i,r) for i,r in enumerate(rows)))
    p = sh(['showalign','-sequence','a.fa','-auto','-outfile','o.txt','-show','A'])
    o = open(os.path.join(D,'o.txt'),errors='replace').read() if os.path.exists(os.path.join(D,'o.txt')) else ''
    return norm(o)+'|ERR:'+(p.stderr.strip()[-80:] if p.returncode else '')


def alignment_sweep():
    for name, f in (('consambig', cons), ('showalign', show)):
        for which, rows in (('nuc', ['ACGTTAGCAT', 'ACGTTAGGAT', 'ACCTTAGCAT']), ('prot', ['ACDEFGHIKL', 'ACDEFGHIKM', 'ACDEYGHIKL'])):
            ref = f(rows); g = collections.defaultdict(str)
            for o in range(1, 128):
                c = chr(o)
                if c in '\n\r>' or c.isalpha(): continue
                d = [r[:5] + c + r[5:] for r in rows]
                r = f(d)
                if r == ref: cls = 'drop'
                elif r == f([r_[:5] + '-' + r_[5:] for r_ in rows]): cls = 'gap'
                elif r == f([r_[:5] + 'N' + r_[5:] for r_ in rows]): cls = 'N'
                elif r == f([r_[:5] + 'X' + r_[5:] for r_ in rows]): cls = 'X'
                elif r == f([r_[:5] + '*' + r_[5:] for r_ in rows]): cls = 'star'
                else: cls = 'other'
                g[cls] += c if c.isprintable() else '\\x%02x' % o
            print(name, which, dict(g))


if __name__ == '__main__':
    sel = sys.argv[1:] or list(progs)
    plain_sweep([x for x in sel if x in progs])
    if not sys.argv[1:] or any(x in ('consambig', 'showalign') for x in sys.argv[1:]):
        alignment_sweep()
