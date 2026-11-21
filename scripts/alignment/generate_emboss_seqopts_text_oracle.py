#!/usr/bin/env python3
"""A6-2 (B09 F177) oracle, effect level: the EMBOSS 6.6.0 programs run with the sequence qualifiers written as command-line TEXT after the sequence arguments (`needle a.fa b.fa -sbegin 3 -sre2 ...`: unnumbered, numbered,
abbreviated, `=` forms, any case; see generate_emboss_acd_args_oracle.py for the grammar) and what they REPORT (the same observables as emboss_seqopts_oracle.tsv, produced by the same program runners of
generate_emboss_seqopts_oracle.py). The port resolves the text with EmbossSequenceQualifiers.ParseInputs(program, text) and must produce the same report, or refuse exactly the texts the program refuses (ERROR).
TSV: id  program  kind  params  seq1  seq2  seq3  text  expected      (params / seqs / expected as in emboss_seqopts_oracle.tsv)
Usage: python3 scripts/alignment/generate_emboss_seqopts_text_oracle.py OUTFILE [per_program=40] [seed]   (needs the EMBOSS 6.6.0 programs on PATH)
"""
import os, random, sys, tempfile, collections

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_emboss_seqopts_oracle as old          # noqa: E402
import generate_emboss_acd_args_oracle as acd         # noqa: E402

MARK = '\x00'
EFFECT_NAMES = acd.CORE + ['sask', 'squick']
state = {'n': 0, 'text': [], 'two': False}


def text_item(rnd, two):
    """acd.item restricted to the sequence qualifiers that act on a residue string (no program options, no file-format values)."""
    return acd.item(rnd, two, EFFECT_NAMES, decoys=False, strays=False)


def make_text(rnd, two):
    out = []
    for _ in range(rnd.choice([1, 1, 2, 2, 3, 3, 4])):
        out += text_item(rnd, two)
    return out


def fake_rand_q(rnd, n, kind, plain=False):
    state['n'] += 1
    return {'_t': state['text'] if state['n'] == 1 else []}


def fake_qargs(q, sfx):
    return [MARK + t for t in q['_t']]


class TextRunner(old.Runner):
    def run(self, args):
        text = [a[1:] for a in args if a.startswith(MARK)]
        rest = [a for a in args if not a.startswith(MARK)]
        prog, tail = rest[0], rest[1:]
        positionals, others = [], []
        i = 0
        while i < len(tail):
            t = tail[i]
            if t in ('-asequence', '-bsequence', '-sequence') and i + 1 < len(tail):
                positionals.append((0 if t != '-bsequence' else 1, tail[i + 1])); i += 2; continue
            if not t.startswith('-') and not others:
                positionals.append((len(positionals), t)); i += 1; continue
            others.append(t); i += 1
        positionals.sort()
        return super().run([prog] + [p for _, p in positionals] + text + others)


def main(out, per, seed):
    rnd = random.Random(seed)
    old.rand_q = fake_rand_q
    old.qargs = fake_qargs
    rows = []
    errors = collections.Counter(); counts = collections.Counter()
    with tempfile.TemporaryDirectory() as d:
        rn = TextRunner(d)
        for name, fn, kinds in old.PROGRAMS:
            made = attempts = 0
            two = name in acd.TWO
            while made < per and attempts < per * 8:
                attempts += 1
                kind = kinds[made % len(kinds)] if len(kinds) > 1 else kinds
                text = make_text(rnd, two)
                state['n'] = 0; state['text'] = text
                res = fn(rn, rnd, kind)
                if res is None: continue
                params, seqs, _qa, _qb, expected = res
                if any('\t' in s or '\n' in s for s in seqs): continue
                if expected == 'ERROR' and errors[name] > per * 0.4: continue
                rows.append('\t'.join([str(len(rows)), name, kind, params] + seqs + [' '.join(text), expected]))
                made += 1
                if expected == 'ERROR': errors[name] += 1
            counts[name] = made
    with open(out, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('id\tprogram\tkind\tparams\tseq1\tseq2\tseq3\ttext\texpected\n' + '\n'.join(rows) + '\n')
    print(len(rows), 'rows', dict(counts), 'errors', dict(errors))


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 40, int(sys.argv[3]) if len(sys.argv) > 3 else 20261221)
