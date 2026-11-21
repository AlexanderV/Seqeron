#!/usr/bin/env python3
"""A6-2 (B09 F177) oracle: how the EMBOSS 6.6.0 ACD library resolves the TEXT of sequence qualifiers that follows the sequence arguments of a program (`needle a.fa b.fa <text> ...`).

For each of the 16 ported programs random texts are generated from a grammar of everything a user can write - numbered / unnumbered names, unique and ambiguous abbreviations, `-` / `--` / `/` / no dash (with `=`), any letter case,
`=value` and separate values, `-no` negations, boolean values, odd integers - and the program is run with EMBOSS_ACDLOG=1; the ACD processing log (`<program>.acdlog`: `set qualifier -sbegin[2] (param 0) = 3`, `set next qualifier ...`)
is the program's own statement of what every token of the text set, so the oracle is the program's resolution, not an assumption about it.
TSV: id  program  text  expected, where expected is
  ERROR     the program died while resolving the text (unknown / ambiguous qualifier, value required, `no` prefix with a value) or refused a value of a sequence qualifier (Bad value for '-sbegin' ...)
  REFUSED   the program accepted a token that is not a qualifier of a sequence input (a parameter, -stdout, -gapopen ...): the port takes sequence qualifiers only and refuses it
  `name[pnum]=value;...`   the FINAL value of every sequence qualifier the text set (the last assignment of each item), sorted by name and number
  `?name[pnum]=value;...`  the same, but the program died AFTER resolving the text for a reason that hides whether the values are valid (it could not read the sequence for -sformat gcg, the type check ...): only the resolution is compared
Usage: python3 scripts/alignment/generate_emboss_acd_args_oracle.py OUTFILE [per_program=130] [seed]   (needs the EMBOSS 6.6.0 programs on PATH)
"""
import os, random, re, subprocess, sys, tempfile, collections

TWO = {'needle', 'water', 'stretcher', 'matcher', 'supermatcher', 'needleall', 'merger', 'megamerger', 'dotmatcher'}
NUC = 'ACGT'
SEQ_NAMES = ['sbegin', 'send', 'sreverse', 'sask', 'snucleotide', 'sprotein', 'slower', 'supper', 'scircular', 'squick', 'sformat', 'iquery', 'ioffset', 'sdbname', 'sid', 'ufo', 'fformat', 'fopenfile']
INT_NAMES = {'sbegin', 'send', 'ioffset'}
STR_NAMES = {'sformat', 'iquery', 'sdbname', 'sid', 'ufo', 'fformat', 'fopenfile'}
CORE = ['sbegin', 'send', 'sreverse', 'slower', 'supper', 'scircular', 'snucleotide', 'sprotein']
DECOYS = ['gapopen', 'stdout', 'auto', 'sequence', 'asequence', 'outfile', 'pattern', 'brief', 'wordsize', 'complement', 'datafile', 'aformat', 'rformat']

TAIL = {
    'needle': ['-auto', '-stdout', '-aformat3', 'pair'], 'water': ['-auto', '-stdout', '-aformat3', 'pair'], 'stretcher': ['-auto', '-stdout', '-aformat3', 'pair'],
    'matcher': ['-auto', '-stdout', '-aformat3', 'pair'], 'supermatcher': ['-wordlen', '3', '-auto', '-stdout', '-aformat3', 'pair'], 'needleall': ['-auto', '-stdout', '-aformat3', 'pair'],
    'merger': ['-outfile', 'o.txt', '-outseq', 'o.fa', '-osformat', 'fasta', '-auto'], 'megamerger': ['-wordsize', '8', '-outseq', 'o.fa', '-osformat', 'fasta', '-outfile', 'o.txt', '-auto'],
    'seqmatchall': ['-wordsize', '4', '-outfile', 's.out', '-aformat', 'match', '-auto'], 'fuzznuc': ['-pattern', 'ACG', '-rformat', 'gff', '-stdout', '-auto'],
    'fuzzpro': ['-pattern', 'A-C', '-rformat', 'gff', '-stdout', '-auto'], 'fuzztran': ['-pattern', 'L-x-S', '-frame', 'F', '-table', '0', '-rformat', 'gff', '-stdout', '-auto'],
    'etandem': ['-outfile', 'o.rep', '-origfile', 'o.orig', '-auto'], 'equicktandem': ['-outfile', 'o.rep', '-origfile', 'o.orig', '-auto'],
    'dotmatcher': ['-windowsize', '4', '-threshold', '20', '-graph', 'data', '-goutfile', 'dm', '-auto'], 'distmat': ['-outfile', 'o.txt', '-auto'],
}
PROGRAMS = list(TAIL)


def pos(rnd, n):
    c = rnd.random()
    if c < 0.2: return 0
    if c < 0.55: return rnd.randint(1, n)
    if c < 0.75: return -rnd.randint(1, n)
    return rnd.choice([n + 3, -n - 3, 1, -1])


def int_value(rnd):
    c = rnd.random()
    if c < 0.75: return str(pos(rnd, 30))
    return rnd.choice(['+3', '007', '-0', '+0', 'abc', '1.5', '', '1e3', '0x10', '2147483648', '-2147483649', '9223372036854775807', '9223372036854775808', '99999999999999999999', '-9223372036854775808', '5-', '--5', '٣', '010', ' 3'.strip()])


BOOL_TOKENS = ['y', 'Y', 'yes', 'YES', 'n', 'N', 'no', 'No', 't', 'T', 'true', 'TRUE', 'f', 'F', 'false', '0', '1']
BOOL_EQ = BOOL_TOKENS + ['maybe', 'tr', 'ye', '2', '0.5', '1e2', '00', '01', '-1', '+0', '0x1', 'nan', 'inf', '1e999', 'yEs', '']


def mangle_case(rnd, s):
    c = rnd.random()
    if c < 0.8: return s
    if c < 0.9: return s.upper()
    return ''.join(ch.upper() if rnd.random() < 0.5 else ch for ch in s)


def item(rnd, two, names=None, decoys=True, strays=True):
    """one qualifier of the text as a list of tokens; names = the sequence qualifiers to draw from, decoys = also draw program qualifiers / parameters (which the port refuses)"""
    names = names or SEQ_NAMES
    name = rnd.choice(CORE) if rnd.random() < 0.8 else (rnd.choice(names) if (not decoys or rnd.random() < 0.8) else rnd.choice(DECOYS))
    is_seq = name in SEQ_NAMES
    is_bool = is_seq and name not in INT_NAMES and name not in STR_NAMES
    spelled = name
    c = rnd.random()
    if c < 0.30 and len(name) > 1: spelled = name[:rnd.randint(1, len(name) - 1)]
    elif c < 0.32: spelled = name[:-1] + rnd.choice('xyz') + name[-1:]
    neg = ''
    if is_bool and rnd.random() < 0.22: neg = 'no'
    number = ''
    c = rnd.random()
    if c < 0.40: number = ''
    elif c < 0.68: number = '1'
    elif c < 0.92: number = '2' if two else '1'
    elif c < 0.96: number = '3'
    elif c < 0.98: number = '0'
    else: number = '01'
    dash = rnd.choices(['-', '--', '/', ''], [80, 8, 6, 6])[0]
    form = rnd.random()
    if is_bool:
        r = rnd.random()
        eq = r >= 0.78
        if not strays and dash == '' and not eq: dash = '-'              # a word without a dash and without '=' is a parameter of the program, not a qualifier
        qual = dash + mangle_case(rnd, neg + spelled + number)
        if r < 0.55 or (not strays and neg and not eq): return [qual]    # -noqual takes no value: a following token would be a parameter
        if r < 0.78: return [qual, rnd.choice(BOOL_TOKENS)]
        return [qual + '=' + rnd.choice(BOOL_EQ)]
    if not strays and dash == '' and form >= 0.22: dash = '-'
    qual = dash + mangle_case(rnd, neg + spelled + number)
    value = int_value(rnd) if name in INT_NAMES else rnd.choice(['fasta', 'x', 'fastq', 'gcg', 'abc', '3', '-3', '-fasta', 'sbegin', '-sbegin', 'a=b'])
    if not is_seq: value = rnd.choice(['5', '10', 'x', '4'])
    if form < 0.22 or value == '': return [qual + '=' + value]           # an empty word on the command line would be an empty argument
    if form < 0.30: return [qual]                         # a value-taking qualifier with no value (the next token is whatever follows)
    return [qual, value]


def make_text(rnd, two):
    out = []
    for _ in range(rnd.choice([1, 1, 2, 2, 3, 3, 4, 5])):
        out += item(rnd, two)
    if rnd.random() < 0.03: out.append(rnd.choice(['stray', '-', '--', 'x=1']))
    return out


def seq_files(d, prog, rnd):
    def dna(n): return ''.join(rnd.choice(NUC) for _ in range(n))
    prot = lambda n: ''.join(rnd.choice('ACDEFGHIKLMNPQRSTVWY') for _ in range(n))
    n = rnd.randint(30, 60)
    base = prot(n) if prog == 'fuzzpro' else dna(n)
    files = {'a.fa': '>a\n' + base + '\n', 'b.fa': '>b\n' + (base[3:] + base[:3]) + '\n'}
    if prog in ('seqmatchall', 'distmat', 'etandem', 'equicktandem', 'fuzznuc', 'fuzzpro', 'fuzztran'):
        files['s.fa'] = '>s1\n' + base + '\n>s2\n' + (base[2:] + base[:2]) + '\n'
        files['a.fa'] = files['s.fa'] if prog in ('seqmatchall', 'distmat') else files['a.fa']
    for f, t in files.items():
        with open(os.path.join(d, f), 'w') as h: h.write(t)
    return ['a.fa', 'b.fa'] if prog in TWO else (['s.fa'] if prog in ('etandem', 'equicktandem') else ['a.fa'])


def classify(prog, text, d, positionals, p, seqset):
    """expected from the acdlog of the first run and the status / messages of a second run without logging (fuzztran crashes after the resolution when it logs)."""
    log = os.path.join(d, prog + '.acdlog')
    if not os.path.exists(log): return None
    lines = open(log, errors='replace').read().split('\n')
    # blocks: 'argv[i] <tok> + argv[i+1] ...' followed by set lines until the next block
    first = len(positionals) + 1
    last = first + len(text) - 1
    block = None
    sets = collections.OrderedDict()
    seen = set()
    params_in_text = False
    for ln in lines:
        m = re.match(r'argv\[(\d+)\] <', ln)
        if m: block = int(m.group(1)); continue
        if block is None or block < first or block > last: continue
        seen.add(block)
        m = re.match(r'(set qualifier|set next qualifier) -(\w+)\[(\d+)\] \(param \d+\) = (.*)$', ln)
        if m: sets[(m.group(2), int(m.group(3)))] = m.group(4)
        if re.match(r'Parameter \d+: ', ln) or re.match(r'Parameter \(\d+\) ', ln): params_in_text = True
        if re.match(r'Parameter \d+: ', ln): params_in_text = True
    # a token consumed as a value has no block: the blocks seen must cover every argv index that is not the value of the previous qualifier
    if p.returncode != 0:
        err = p.stderr
        if re.search(r"Unknown qualifier|Value required|'no' prefix|Died: .*ambiguous|Too many parameters|Bad value for '-(sbegin|send|sreverse|sask|snucleotide|sprotein|slower|supper|scircular|squick|sformat|iquery|ioffset|sdbname|sid|ufo|fformat|fopenfile)'", err):
            return 'ERROR'
    if params_in_text: return 'REFUSED'
    if any(k[0] not in SEQ_NAMES for k in sets): return 'REFUSED'
    body = ';'.join('%s[%d]=%s' % (k[0], k[1], v) for k, v in sorted(sets.items()))
    return ('?' if p.returncode != 0 else '') + body


def main(out, per, seed):
    rnd = random.Random(seed)
    rows = []
    stats = collections.Counter()
    with tempfile.TemporaryDirectory() as d:
        for prog in PROGRAMS:
            made = 0
            tries = 0
            while made < per and tries < per * 4:
                tries += 1
                two = prog in TWO
                text = make_text(rnd, two)
                if any(('\t' in t or '\n' in t) for t in text): continue
                positionals = seq_files(d, prog, rnd)
                for f in os.listdir(d):
                    if f.endswith('.acdlog'): os.remove(os.path.join(d, f))
                env = dict(os.environ, EMBOSS_ACDLOG='1')
                try:
                    p = subprocess.run([prog] + positionals + text + TAIL[prog], cwd=d, capture_output=True, text=True, env=env, timeout=60)
                    plain = subprocess.run([prog] + positionals + text + TAIL[prog], cwd=d, capture_output=True, text=True, timeout=60)
                except subprocess.TimeoutExpired:
                    continue
                exp = classify(prog, text, d, positionals, plain, None)
                if exp is None: continue
                if exp == 'ERROR' and any(t.startswith(('-', '/')) and False for t in text): pass
                if exp == 'ERROR' and stats[(prog, 'ERROR')] > per * 0.55: continue
                if exp == 'REFUSED' and stats[(prog, 'REFUSED')] > per * 0.25: continue
                rows.append('\t'.join([str(len(rows)), prog, ' '.join(text), exp]))
                stats[(prog, exp if exp in ('ERROR', 'REFUSED') else 'OK')] += 1
                made += 1
    with open(out, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('id\tprogram\ttext\texpected\n' + '\n'.join(rows) + '\n')
    print(len(rows), 'rows', dict(collections.Counter(k[1] for k in stats.elements())))


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 130, int(sys.argv[3]) if len(sys.argv) > 3 else 20261220)
