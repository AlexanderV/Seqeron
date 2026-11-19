#!/usr/bin/env python3
"""ALIGN-MULTI-001 (audit A2-14c, B09 F106) oracle: ClustalW 2.1 secondary-structure gap-penalty masks of the profile commands.

Seqeron API: ClustalWProfileFileReader.ParseClustal / ParseGde -> SequenceAligner.RunClustalWProfile / RunClustalWSequences (ClustalWRunOptions.ProfileStructure1 / 2, Structure) -> ClustalWAlignmentWriter.Write.

Reference (run here): the installed `clustalw` 2.1 (Debian 2.1+lgpl-7) `-PROFILE1=a.aln -PROFILE2=b.aln -PROFILE | -SEQUENCES` with CLUSTAL-format profile files (written by `clustalw -ALIGN` itself, then given `!SS_name` structure lines
or `!GM_name` gap-penalty-mask lines before the sequence lines of every block; some files carry both, some lose the line in a later block: the program then reports "bad secondary structure format") or GDE-format files
(`"SS_name` entries), -NOSECSTR1 / -NOSECSTR2, -SECSTROUT=STRUCTURE|MASK|BOTH|NONE, -HELIXGAP / -STRANDGAP / -LOOPGAP / -TERMINALGAP / -HELIXENDIN / -HELIXENDOUT / -STRANDENDIN / -STRANDENDOUT, -OUTPUT=CLUSTAL|GDE|FASTA|PIR
and the alignment option variants of the earlier ClustalW oracles; the whole output file (rows + structure lines) must be equal.
Output TSV (tests/.../TestData/Alignment/clustalw_structure_oracle.tsv):  id  command(profile|sequences)  type  inputformat(clustal|gde)  file1  file2  args  format  output
(files and output with \\n / \\t / \\\\ escaped; args exclude the file names / -TYPE / -OUTPUT)
Usage: python3 scripts/alignment/generate_clustalw_structure_oracle.py OUTFILE [sets] [seed]      (needs `clustalw` on PATH)
"""
import importlib.util, os, random, re, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + '.py'))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


dna = load('generate_clustalw_msa_full_oracle')
prot = load('generate_clustalw_protein_msa_oracle')
FLAGS = dict(prot.FLAGS)
FLAGS['transweight'] = 'TRANSWEIGHT'


def esc(t):
    return t.replace('\\', '\\\\').replace('\n', '\\n').replace('\t', '\\t')


def cw(args, cwd):
    return subprocess.run(['clustalw'] + args, capture_output=True, text=True, cwd=cwd, timeout=30)


def structure_string(rnd, L, kind):
    if kind == 'gm':
        return ''.join(rnd.choice('123456789') if rnd.random() < 0.9 else rnd.choice('0.-') for _ in range(L))
    out = []
    while len(out) < L:
        r = rnd.random()
        k = rnd.randint(1, 12)
        if r < 0.3:
            seg = [rnd.choice('Aa')] * k
            if rnd.random() < 0.4: seg[-1] = '$'
        elif r < 0.55:
            seg = [rnd.choice('Bb')] * k
            if rnd.random() < 0.4: seg[-1] = '%'
        else:
            seg = [rnd.choice('.-~C?')] * k
        out.extend(seg)
    return ''.join(out[:L])


def add_clustal_structure(rnd, text, lines_spec):
    """lines_spec: list of (prefix '!SS_' / '!GM_', name, full string); inserts one line per spec before the first sequence line of every block"""
    lines = text.split('\n')
    out = [lines[0]]
    i = 1
    pos = 0
    while i < len(lines):
        line = lines[i]
        if line.strip() == '' or line[0] == ' ':
            out.append(line); i += 1; continue
        # first line of a block: its residue token length gives the block width
        tokens = line.split()
        width = len(tokens[1]) if len(tokens) > 1 else 0
        for prefix, name, full, drop_from in lines_spec:
            if drop_from is not None and pos >= drop_from: continue
            chunk = full[pos:pos + width]
            out.append(f'{prefix}{name}  {chunk}')
        while i < len(lines) and lines[i].strip() != '' and lines[i][0] != ' ':
            out.append(lines[i]); i += 1
        pos += width
    return '\n'.join(out)


def aligned_profile(d, tag, names, seqs, typ, rnd, structure):
    """the CLUSTAL-format file of an alignment (from the program), with structure lines; returns (text, columns)"""
    fa = os.path.join(d, tag + '0.fa')
    with open(fa, 'w') as f:
        for n, s in zip(names, seqs): f.write(f'>{n}\n{s}\n')
    out = os.path.join(d, tag + '0.aln')
    if len(seqs) == 1:
        cw([f'-INFILE={tag}0.fa', f'-TYPE={typ.upper()}', '-CONVERT', '-OUTPUT=CLUSTAL', f'-OUTFILE={tag}0.aln'], d)
    else:
        cw([f'-INFILE={tag}0.fa', f'-TYPE={typ.upper()}', '-ALIGN', '-OUTPUT=CLUSTAL', '-OUTORDER=INPUT', f'-OUTFILE={tag}0.aln'], d)
    text = open(out).read()
    L = 0
    for line in text.split('\n')[1:]:
        t = line.split()
        if t and line[0] != ' ' and len(t) > 1 and re.fullmatch(r'[A-Za-z\-]+', t[1]):
            pass
    # alignment length = sum of the residue-token lengths of the first sequence over the blocks
    first = names[0]
    for line in text.split('\n')[1:]:
        if line.startswith(first + ' ') or line.startswith(first[:30] + ' '):
            t = line.split()
            if len(t) > 1: L += len(t[1])
    return text, L


def gde_text(d, tag, typ, names, seqs):
    fa = os.path.join(d, tag + '0.fa')
    with open(fa, 'w') as f:
        for n, s in zip(names, seqs): f.write(f'>{n}\n{s}\n')
    out = os.path.join(d, tag + '0.gde')
    if len(seqs) == 1:
        cw([f'-INFILE={tag}0.fa', f'-TYPE={typ.upper()}', '-CONVERT', '-OUTPUT=GDE', f'-OUTFILE={tag}0.gde'], d)
    else:
        cw([f'-INFILE={tag}0.fa', f'-TYPE={typ.upper()}', '-ALIGN', '-OUTPUT=GDE', '-OUTORDER=INPUT', f'-OUTFILE={tag}0.gde'], d)
    return open(out).read()


def case(rnd, i):
    typ = 'dna' if i % 2 == 0 else 'protein'
    command = 'profile' if (i // 2) % 2 == 0 else 'sequences'
    m = dna if typ == 'dna' else prot
    dna.rnd.seed(rnd.randint(0, 2 ** 31)); prot.rnd.seed(rnd.randint(0, 2 ** 31))
    maxlen = rnd.choice([30, 50, 80, 120, 160])
    seqs = dna.make_set(maxlen, 'ACGT', 2, rnd.choice([4, 6, 9])) if typ == 'dna' else prot.make_set(maxlen, 2, rnd.choice([4, 6, 9]))
    if len(seqs) < 2: seqs = seqs + seqs[:1]
    names = [f's{j:03d}' for j in range(len(seqs))]
    cut = rnd.randint(1, len(seqs) - 1)
    informat = 'clustal' if rnd.random() < 0.8 else 'gde'
    with tempfile.TemporaryDirectory() as d:
        files = []
        for tag, nn, ss in (('a', names[:cut], seqs[:cut]), ('b', names[cut:], seqs[cut:])):
            has_structure = rnd.random() < (0.9 if tag == 'a' else 0.6) or (command == 'sequences' and len(ss) == 1 and rnd.random() < 0.8)
            if informat == 'clustal':
                text, L = aligned_profile(d, tag, nn, ss, typ, rnd, None)
                if has_structure and L > 0:
                    sname = rnd.choice([f'ss{tag}', 'helix_strand_mask', 'x' * rnd.randint(1, 40)])
                    mode = rnd.random()
                    spec = []
                    if mode < 0.62: spec = [('!SS_', sname, structure_string(rnd, L, 'ss'), None)]
                    elif mode < 0.8: spec = [('!GM_', sname, structure_string(rnd, L, 'gm'), None)]
                    elif mode < 0.9: spec = [('!SS_', sname, structure_string(rnd, L, 'ss'), None), ('!GM_', sname + '2', structure_string(rnd, L, 'gm'), None)]
                    elif mode < 0.95: spec = [('!GM_', sname, structure_string(rnd, L, 'gm'), None), ('!SS_', sname + '2', structure_string(rnd, L, 'ss'), None)]
                    else: spec = [('!SS_', sname, structure_string(rnd, L, 'ss'), rnd.randint(1, L))]   # structure line missing from a later block
                    text = add_clustal_structure(rnd, text, spec)
            else:
                text = gde_text(d, tag, typ, nn, ss)
                L = max(len(x) for x in ss) + 0
                if has_structure:
                    # the aligned length: from the alignment itself (first sequence of the GDE file)
                    body = [l for l in text.split('\n') if l and l[0] not in '#%"']
                    L = len(''.join(body)) // max(1, len(nn))
                    full = structure_string(rnd, L, 'ss' if rnd.random() < 0.75 else 'gm')
                    tag2 = 'SS_' if set(full) - set('123456789') else 'GM_'
                    if tag2 == 'GM_' and rnd.random() < 0.3: tag2 = 'SS_'
                    sname = rnd.choice([f'ss{tag}', 'helix_mask'])
                    width = 60
                    text = text.rstrip('\n') + '\n"' + tag2 + sname + '\n' + '\n'.join(full[k:k + width] for k in range(0, L, width)) + '\n'
            files.append(text)
        with open(os.path.join(d, 'a.in'), 'w') as f: f.write(files[0])
        with open(os.path.join(d, 'b.in'), 'w') as f: f.write(files[1])
        o = {}
        m.rnd.seed(rnd.randint(0, 2 ** 31))
        o = m.options()
        for k in ('ktuple', 'window', 'topdiags', 'pairgap', 'novgap', 'endgaps', 'hgapresidues'): o.pop(k, None)
        args = [f'-{FLAGS[k]}' + ('' if v is None else f'={v}') for k, v in o.items()]
        if rnd.random() < 0.15: args.append('-NOSECSTR1')
        if rnd.random() < 0.15: args.append('-NOSECSTR2')
        if rnd.random() < 0.5: args.append('-SECSTROUT=' + rnd.choice(['STRUCTURE', 'MASK', 'BOTH', 'NONE']))
        for flag, lo, hi in (('HELIXGAP', 1, 9), ('STRANDGAP', 1, 9), ('LOOPGAP', 1, 9), ('TERMINALGAP', 1, 9), ('HELIXENDIN', 0, 3), ('HELIXENDOUT', 0, 3), ('STRANDENDIN', 0, 3), ('STRANDENDOUT', 0, 3)):
            if rnd.random() < 0.22: args.append(f'-{flag}={rnd.randint(lo, hi)}')
        fmt = rnd.choice(['CLUSTAL', 'CLUSTAL', 'CLUSTAL', 'GDE', 'GDE', 'FASTA', 'PIR'])
        ext = {'CLUSTAL': 'aln', 'GDE': 'gde', 'FASTA': 'fasta', 'PIR': 'pir'}[fmt]
        flag = '-PROFILE' if command == 'profile' else '-SEQUENCES'
        r = cw(['-PROFILE1=a.in', '-PROFILE2=b.in', flag, f'-TYPE={typ.upper()}', f'-OUTPUT={fmt}', f'-OUTFILE=o.{ext}'] + args, d)
        out = os.path.join(d, 'o.' + ext)
        text = open(out).read() if os.path.exists(out) and os.path.getsize(out) > 0 else ''
        return [i, command, typ, informat, esc(files[0]), esc(files[1]), ' '.join(args), fmt, esc(text)], bool(text)


def main(outfile, n, seed):
    rnd = random.Random(seed)
    rows = []
    i = 0
    while len(rows) < n:
        row, ok = case(rnd, len(rows) + 0)
        if ok: rows.append(row)
        i += 1
    with open(outfile, 'w') as f:
        f.write('id\tcommand\ttype\tinputformat\tfile1\tfile2\targs\tformat\toutput\n')
        for r in rows: f.write('\t'.join(str(x) for x in r) + '\n')
    from collections import Counter
    print(len(rows), 'sets', Counter((r[1], r[2], r[3]) for r in rows))


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 1200, int(sys.argv[3]) if len(sys.argv) > 3 else 20261123)
