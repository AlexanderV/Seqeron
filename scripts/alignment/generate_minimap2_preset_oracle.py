#!/usr/bin/env python3
"""A5-9 (B09 F173) oracle: the `-x` presets of minimap2 2.26 (`options.c` `mm_set_opt`) up to the end of the chaining stage, as PAF TEXT (columns 1 - 12 and the tags tp cm s1 s2 dv rl of every line of `minimap2 -x PRESET ref.fa query.fa`).

JSON lines: {"preset": "...", "extra": [...], "names": [...], "targets": [...], "qname": "...", "query": "...", "paf": "<the program's stdout>"}
Presets (>= 250 cases each): map-ont, map-pb, map-hifi, ava-ont, ava-pb, asm5, asm10, asm20, sr, splice, splice:hq (the last two are recorded to measure how far the chaining stage reaches; the test
compares the presets the port reproduces and states which it does not). Cases: a reference of 1 - 3 random sequences (1.5 - 14 kb; for the ava presets several reads and the query among them, so that the
self hit is skipped), the query a mutated (1 - 15 %: substitutions, deletions, insertions), possibly reverse-complemented, possibly partial copy of a reference stretch, with repeats (a second copy elsewhere, sometimes
inverted), homopolymer runs (the HPC presets), 80 - 250 bp reads for sr, spliced queries (exons of a gene joined) for splice.
Usage: python3 scripts/alignment/generate_minimap2_preset_oracle.py [per_preset=220] > minimap2_preset_oracle.jsonl      (needs minimap2 2.26 on PATH)
"""
import json, random, subprocess, sys, tempfile

rnd = random.Random(20261210)
COMP = str.maketrans('ACGT', 'TGCA')


def rseq(n, hp=False):
    s = []
    while len(s) < n:
        c = rnd.choice('ACGT')
        s.extend([c] * (rnd.randint(3, 9) if hp and rnd.random() < 0.15 else 1))
    return ''.join(s[:n])


def revcomp(s): return s.translate(COMP)[::-1]


def mutate(s, rate):
    out = []
    for c in s:
        r = rnd.random()
        if r < rate * 0.5: out.append(rnd.choice('ACGT'))
        elif r < rate * 0.75: pass
        elif r < rate: out.append(c); out.append(rnd.choice('ACGT'))
        else: out.append(c)
    return ''.join(out)


def make(preset):
    hp = preset in ('map-pb', 'ava-pb', 'map-hifi')
    if preset in ('ava-ont', 'ava-pb'):
        n = rnd.randint(3, 7)
        base = rseq(rnd.randint(3000, 9000), hp)
        reads = []
        for i in range(n):
            a = rnd.randint(0, len(base) // 3); b = rnd.randint(2 * len(base) // 3, len(base))
            r = mutate(base[a:b], rnd.choice([0.02, 0.08, 0.14]))
            reads.append(revcomp(r) if rnd.random() < 0.4 else r)
        names = ['read%d' % i for i in range(n)]
        qi = rnd.randrange(n)
        return names, reads, names[qi], reads[qi]
    if preset == 'sr':
        ref = rseq(rnd.randint(4000, 20000))
        refs = [ref]
        if rnd.random() < 0.4: refs.append(ref[:rnd.randint(2000, 4000)] + rseq(300) + ref[2000:rnd.randint(4000, 6000)])
        a = rnd.randint(0, len(ref) - 300)
        q = mutate(ref[a:a + rnd.randint(80, 250)], rnd.choice([0.0, 0.01, 0.03, 0.06]))
        if rnd.random() < 0.5: q = revcomp(q)
        return ['ref%d' % i for i in range(len(refs))], refs, 'rd%d' % rnd.randint(0, 99999), q
    if preset.startswith('splice'):
        gene = rseq(rnd.randint(3000, 9000))
        exons = []; pos = rnd.randint(0, 200)
        for _ in range(rnd.randint(2, 5)):
            ln = rnd.randint(80, 300)
            exons.append(gene[pos:pos + ln]); pos += ln + rnd.randint(100, 1200)
        q = mutate(''.join(exons), rnd.choice([0.0, 0.02, 0.05]))
        if rnd.random() < 0.4: q = revcomp(q)
        return ['gene'], [gene], 'tx%d' % rnd.randint(0, 9999), q
    ln = rnd.randint(1500, 14000)
    ref = rseq(ln, hp)
    refs = [ref]
    if rnd.random() < 0.3: refs.append(rseq(rnd.randint(1000, 8000), hp))
    a = rnd.randint(0, ln // 2); b = rnd.randint(a + 700, min(ln, a + 7000))
    seg = ref[a:b]
    rate = {'map-ont': rnd.choice([0.03, 0.08, 0.15]), 'map-pb': rnd.choice([0.05, 0.12]), 'map-hifi': rnd.choice([0.002, 0.01, 0.03]),
            'asm5': rnd.choice([0.002, 0.005, 0.02]), 'asm10': rnd.choice([0.01, 0.05]), 'asm20': rnd.choice([0.03, 0.1, 0.17])}.get(preset, 0.05)
    q = mutate(seg, rate)
    if rnd.random() < 0.35:                      # a repeat elsewhere (a secondary hit), sometimes inverted
        c = mutate(seg[:len(seg) // 2], rate * 1.5)
        refs[0] = ref + rseq(rnd.randint(200, 3000)) + (revcomp(c) if rnd.random() < 0.5 else c)
    if rnd.random() < 0.4: q = revcomp(q)
    if rnd.random() < 0.15: q = q[:rnd.randint(len(q) // 3, len(q))]
    return ['ref%d' % i for i in range(len(refs))], refs, 'q%d' % rnd.randint(0, 99999), q


def main():
    per = int(sys.argv[1]) if len(sys.argv) > 1 else 220
    with tempfile.TemporaryDirectory() as d:
        for preset in ['map-ont', 'map-pb', 'map-hifi', 'ava-ont', 'ava-pb', 'asm5', 'asm10', 'asm20', 'sr', 'splice', 'splice:hq']:
            made = 0; tries = 0
            while made < per and tries < per * 5:
                tries += 1
                names, targets, qname, query = make(preset)
                open(d + '/r.fa', 'w').write(''.join('>%s\n%s\n' % (n, t) for n, t in zip(names, targets)))
                open(d + '/q.fa', 'w').write('>%s\n%s\n' % (qname, query))
                extra = []                      # -F (max fragment length, effective with sr) and -G (max gap on the reference = -r b,B) after the preset, as the program applies them
                c = rnd.random()
                if preset == 'sr' and c < 0.3: extra = ['-F', str(rnd.choice([150, 300, 500, 2000]))]
                elif c < 0.2: extra = ['-G', str(rnd.choice([500, 2000, 10000, 50000]))]
                p = subprocess.run(['minimap2', '-x', preset] + extra + [d + '/r.fa', d + '/q.fa'], capture_output=True, text=True)
                if p.returncode != 0: continue
                made += 1
                print(json.dumps({'preset': preset, 'extra': extra, 'names': names, 'targets': targets, 'qname': qname, 'query': query, 'paf': p.stdout}))


if __name__ == '__main__':
    main()
