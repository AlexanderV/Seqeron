#!/usr/bin/env python3
"""Oracle generator for the T-Coffee methods that run an EXTERNAL program (audit A2-9 / A2-10 / A2-11, B09 F114+).

Runs the installed T-Coffee program (Debian t-coffee 13.45.0.4846264+really13.41.0.28bdc39+dfsg-1, `t_coffee`) on random protein sets together with
the external programs it drives, and writes `tcoffee_external_oracle.tsv`:

  feature <TAB> flags <TAB> seqs (a|b|c) <TAB> rows (aligned FASTA rows in input order, or FAIL:<rc>) <TAB> extra (compact JSON, see below)

  TMalign_pair / mustang_pair   (`-method TMalign_pair -template_file tpl`; every sequence may carry a `_P_` template = a CA-only PDB file that this script
                                 generates deterministically: a random-walk backbone, related members share the backbone with indels, substitutions, a rigid motion
                                 and coordinate noise)
        extra = {"pdb": [PDB sequence of every template ("" = no template), ...],
                 "out": {"TMalign_pair": {"i,j": <what the external program wrote for the pair of templates i<j, exactly the text t_coffee reads>, ...}, "mustang_pair": {...}}}
                (only the structure methods named in the flags)
        for TMalign the LAST EIGHT lines of its standard output (t_coffee's glue script keeps `tail -4`), for mustang the file results.afasta.
  blastp_msa                    (`-method blastp_msa`: the sequences are formatted as a local BLAST database (legacy `formatdb` = makeblastdb) and searched
                                 against themselves with `blastall -p blastp -m7`)
        extra = {"xmlz": base64 of the zlib-compressed BLAST XML (-m7 / -outfmt 5) the program read,
                 "lib": for the flags `-method blastp_msa` alone: the library `-lib_only -out_lib` writes ("a:p:b:q:weight:count;...", 0-based sequences, 1-based residues)}

The external programs run on the SAME files the program hands them (the templates are copied byte for byte), so that a port of t_coffee's side of the interface can be fed
with their output and must reproduce every row. Every run uses `-n_core 1 -no_warning` and a private HOME whose `.t_coffee/plugins/linux` links the installed helpers (the program
extends PATH / MAFFT_BINARIES with that directory, see SequenceAligner.TCoffeeLibrary.cs and docs): TMalign, mustang, legacy_blast.pl (= /usr/bin/legacy_blast of ncbi-blast+-legacy).

usage: generate_tcoffee_external_oracle.py OUT_DIR [FEATURE ...]   (env SCALE, SEED)
"""
import base64, json, math, os, random, subprocess, sys, tempfile, zlib
from concurrent.futures import ThreadPoolExecutor

PROTEIN = 'ACDEFGHIKLMNPQRSTVWY'
AA3 = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'}


def setup_plugins(d):
    plug = d + '/.t_coffee/plugins/linux'
    os.makedirs(plug, exist_ok=True)
    for src in ('/usr/lib/mafft/bin', '/usr/lib/mafft/lib/mafft'):
        if os.path.isdir(src):
            for fn in os.listdir(src):
                if not os.path.exists(plug + '/' + fn): os.symlink(src + '/' + fn, plug + '/' + fn)
    if os.path.exists('/usr/bin/muscle3') and not os.path.exists(plug + '/muscle'): os.symlink('/usr/bin/muscle3', plug + '/muscle')
    if os.path.exists('/usr/bin/legacy_blast') and not os.path.exists(plug + '/legacy_blast.pl'): os.symlink('/usr/bin/legacy_blast', plug + '/legacy_blast.pl')
    for tool in ('TMalign', 'mustang'):
        if os.path.exists('/usr/bin/' + tool) and not os.path.exists(plug + '/' + tool): os.symlink('/usr/bin/' + tool, plug + '/' + tool)
    if os.environ.get('SAP_REAL') and os.path.exists('/opt/sap/bin/sap') and not os.path.exists(plug + '/sap'): os.symlink('/opt/sap/bin/sap', plug + '/sap')   # the hand-off image A2-12 installs the REAL SAP


# ----- synthetic structures ---------------------------------------------------------------------------------------------------------------------------------------

def walk(n, rng):
    """Persistent random walk with 3.8 A between consecutive CA atoms."""
    p = [(0.0, 0.0, 0.0)]
    d = (1.0, 0.0, 0.0)
    for _ in range(n - 1):
        v = [d[k] + rng.gauss(0, 0.7) for k in range(3)]
        l = math.sqrt(sum(x * x for x in v)) or 1.0
        d = tuple(x / l for x in v)
        p.append(tuple(p[-1][k] + 3.8 * d[k] for k in range(3)))
    return p


def helix(n, rng):
    return [(2.3 * math.cos(math.radians(100 * i)), 2.3 * math.sin(math.radians(100 * i)), 1.5 * i) for i in range(n)]


def rigid(pts, rng, noise):
    a, b, c = [rng.uniform(0, 2 * math.pi) for _ in range(3)]
    ca, sa, cb, sb, cc, sc = math.cos(a), math.sin(a), math.cos(b), math.sin(b), math.cos(c), math.sin(c)
    R = [[cb * cc, -cb * sc, sb], [ca * sc + sa * sb * cc, ca * cc - sa * sb * sc, -sa * cb], [sa * sc - ca * sb * cc, sa * cc + ca * sb * sc, ca * cb]]
    t = [rng.uniform(-20, 20) for _ in range(3)]
    return [tuple(sum(R[i][j] * q[j] for j in range(3)) + t[i] + rng.gauss(0, noise) for i in range(3)) for q in pts]


def pdb_text(seq, pts, chain='A'):
    # t_coffee's is_pdb_file needs "\nHEADER" or "\nSEQRES" and "\nATOM": the first line must not be the only keyword line
    out = ['REMARK   1 SYNTHETIC CA TRACE', 'HEADER    SYNTHETIC                               01-JAN-00   0SYN']
    for i in range(0, len(seq), 13):
        out.append('SEQRES %3d %s %4d  %s' % (i // 13 + 1, chain, len(seq), ' '.join(AA3[c] for c in seq[i:i + 13])))
    for i, (c, q) in enumerate(zip(seq, pts)):
        out.append('ATOM  %5d  CA  %3s %s%4d    %8.3f%8.3f%8.3f  1.00 20.00           C' % (i + 1, AA3[c], chain, i + 1, q[0], q[1], q[2]))
    out += ['TER', 'END']
    return '\n'.join(out) + '\n'


def make_structure_set(rng, krange, lrange, p_none):
    """Returns (sequences, pdb sequences or '' per sequence, coordinates or None per sequence)."""
    k = rng.randint(*krange)
    n = rng.randint(*lrange)
    base = ''.join(rng.choice(PROTEIN) for _ in range(n))
    bpts = walk(n, rng) if rng.random() < 0.7 else helix(n, rng)
    seqs, pdbseqs, coords = [], [], []
    for _ in range(k):
        if rng.random() < 0.8:
            s = list(base); pts = list(bpts)
            for _ in range(rng.randint(0, max(1, n // 5))):
                op = rng.choice('sid'); p = rng.randrange(len(s))
                if op == 's': s[p] = rng.choice(PROTEIN)
                elif op == 'i':
                    s.insert(p, rng.choice(PROTEIN))
                    q = pts[p]; pts.insert(p, (q[0] + rng.uniform(-1.5, 1.5), q[1] + rng.uniform(-1.5, 1.5), q[2] + rng.uniform(-1.5, 1.5)))
                elif len(s) > 8: del s[p]; del pts[p]
        else:
            m = rng.randint(max(8, lrange[0] - 2), lrange[1])
            s = [rng.choice(PROTEIN) for _ in range(m)]; pts = walk(m, rng)
        seq = ''.join(s)
        seqs.append(seq)
        if rng.random() < p_none:
            pdbseqs.append(''); coords.append(None); continue
        ps = list(seq)
        r = rng.random()
        if r < 0.15:                                  # a structure whose sequence differs from the sequence (substitutions): the template is accepted at > 35 % identity / > 50 % coverage
            for _ in range(rng.randint(1, max(1, len(ps) // 3))): ps[rng.randrange(len(ps))] = rng.choice(PROTEIN)
        pp = list(pts)
        if rng.random() < 0.1 and len(ps) > 12:      # a short structure: low coverage
            cut = rng.randint(len(ps) // 3, len(ps) // 2)
            ps = ps[:cut]; pp = pp[:cut]
        pdbseqs.append(''.join(ps))
        coords.append(rigid(pp[:len(ps)], rng, rng.choice([0.0, 0.3, 0.8])))
    return seqs, pdbseqs, coords


# ----- running the program ------------------------------------------------------------------------------------------------------------------------------------------

def read_fasta_rows(path, names):
    rows, cur = {}, None
    for l in open(path):
        l = l.strip()
        if l.startswith('>'): cur = l[1:].split()[0]; rows[cur] = ''
        elif cur: rows[cur] += l
    return '|'.join(rows[n] for n in names)


def run_structure(feature, flags, seqs, pdbseqs, coords):
    names = [f'q{i:03d}' for i in range(len(seqs))]
    with tempfile.TemporaryDirectory() as d:
        setup_plugins(d)
        self_templates = '3dcoffee' in flags       # -mode 3dcoffee = `-template_file SELF_P_`: the templates stand in the headers of the sequence file itself
        files = []
        with open(d + '/tpl.fa', 'w') as tpl:
            for i, (n, ps, pts) in enumerate(zip(names, pdbseqs, coords)):
                if pts is None: files.append(None); continue
                path = f'{d}/p{i}.pdb'
                with open(path, 'w') as g: g.write(pdb_text(ps, pts))
                files.append(path)
                tpl.write(f'>{n} _P_ {path}\n')
        with open(d + '/s.fa', 'w') as f:
            for n, s, path in zip(names, seqs, files): f.write(f'>{n} _P_ {path}\n{s}\n' if self_templates and path else f'>{n}\n{s}\n')
        out = {}
        # sap_pair and `-mode 3dcoffee` (= -in Msap_pair): SAP is not installed, so the program switches to TMalign ("SAP is not installed, TMalign will be used instead")
        want = {'TMalign_pair': any(k in flags for k in ('TMalign_pair', 'sap_pair', '3dcoffee')), 'mustang_pair': 'mustang_pair' in flags}
        for tool in ('TMalign_pair', 'mustang_pair'):
            if not want[tool]: continue
            out[tool] = {}
            for i in range(len(seqs)):
                for j in range(i + 1, len(seqs)):
                    if files[i] is None or files[j] is None: continue
                    if tool == 'TMalign_pair':
                        r = subprocess.run(['TMalign', files[i], files[j]], capture_output=True, text=True, encoding='latin-1', cwd=d)
                        if r.returncode == 0: out[tool][f'{i},{j}'] = '\n'.join(r.stdout.split('\n')[-9:])        # the last eight lines + the final newline
                    else:
                        wd = f'{d}/m{i}_{j}'; os.makedirs(wd)
                        r = subprocess.run(['mustang', '-i', files[i], files[j], '-F', 'fasta'], capture_output=True, text=True, encoding='latin-1', cwd=wd)
                        if r.returncode == 0 and os.path.exists(wd + '/results.afasta'): out[tool][f'{i},{j}'] = open(wd + '/results.afasta').read()   # the glue script keeps the result only when mustang exits with status 0 (it segfaults on some pairs after writing the file)
        env = dict(os.environ, HOME=d)
        cmd = ['t_coffee', '-seq', d + '/s.fa'] + ([] if self_templates else ['-template_file', d + '/tpl.fa']) + ['-n_core', '1', '-no_warning', '-outfile', d + '/o.fa', '-output', 'fasta_aln'] + flags.split()
        r = subprocess.run(cmd, capture_output=True, text=True, encoding='latin-1', cwd=d, env=env)
        try: rows = read_fasta_rows(d + '/o.fa', names)
        except (FileNotFoundError, KeyError): rows = f'FAIL:{r.returncode}'
        extra = {'pdb': pdbseqs, 'out': out}
        cache = d + '/.t_coffee/cache'
        if os.environ.get('SAP_REAL') and os.path.isdir(cache):   # SAP's raw output: t_coffee keeps <template1>_<template2>.sap_results in its cache
            extra['sap'] = {fn: open(cache + '/' + fn, encoding='latin-1').read() for fn in sorted(os.listdir(cache)) if fn.endswith('.sap_results')}
    return rows, extra


def run_plain(flags, seqs):
    names = [f'q{i:03d}' for i in range(len(seqs))]
    with tempfile.TemporaryDirectory() as d:
        setup_plugins(d)
        with open(d + '/s.fa', 'w') as f:
            for n, s in zip(names, seqs): f.write(f'>{n}\n{s}\n')
        env = dict(os.environ, HOME=d)
        cmd = ['t_coffee', '-seq', d + '/s.fa', '-n_core', '1', '-no_warning', '-outfile', d + '/o.fa', '-output', 'fasta_aln'] + flags.split()
        r = subprocess.run(cmd, capture_output=True, text=True, encoding='latin-1', cwd=d, env=env)
        try: return read_fasta_rows(d + '/o.fa', names)
        except (FileNotFoundError, KeyError): return f'FAIL:{r.returncode}'


def run_blast(flags, seqs):
    names = [f'q{i:03d}' for i in range(len(seqs))]
    with tempfile.TemporaryDirectory() as d:
        setup_plugins(d)
        with open(d + '/s.fa', 'w') as f:
            for n, s in zip(names, seqs): f.write(f'>{n}\n{s}\n')
        env = dict(os.environ, HOME=d)
        cmd = ['t_coffee', '-seq', d + '/s.fa', '-n_core', '1', '-no_warning', '-outfile', d + '/o.fa', '-output', 'fasta_aln'] + flags.split()
        r = subprocess.run(cmd, capture_output=True, text=True, encoding='latin-1', cwd=d, env=env)
        try: rows = read_fasta_rows(d + '/o.fa', names)
        except (FileNotFoundError, KeyError): rows = f'FAIL:{r.returncode}'
        lib = None
        if flags == '-method blastp_msa':     # the primary library the program writes (`-lib_only -out_lib`): "a:p:b:q:weight:count;..." (0-based sequences, 1-based residues) as in the library oracle of the program
            subprocess.run(['t_coffee', '-seq', d + '/s.fa', '-n_core', '1', '-no_warning', '-method', 'blastp_msa', '-lib_only', '-out_lib', d + '/o.lib'], capture_output=True, text=True, encoding='latin-1', cwd=d, env=env)
            lib, cur = [], None
            try:
                for l in open(d + '/o.lib'):
                    l = l.strip()
                    if l.startswith('#'): t = l[1:].split(); cur = (int(t[0]), int(t[1]))
                    elif l and not l.startswith('!') and cur is not None:
                        t = l.split()
                        if len(t) >= 4 and t[0].isdigit(): lib.append(f'{cur[0] - 1}:{t[0]}:{cur[1] - 1}:{t[1]}:{t[2]}:{t[3]}')
            except FileNotFoundError:
                lib = None
            lib = ';'.join(lib) if lib is not None else None
        # what the program's glue runs: makeblastdb -dbtype prot on the sequences, blastall -p blastp -m7 of the sequences against them
        subprocess.run(['makeblastdb', '-in', 's.fa', '-dbtype', 'prot', '-logfile', '/dev/null'], capture_output=True, cwd=d, env=env)
        x = subprocess.run(['blastall', '-i', 's.fa', '-d', 's.fa', '-m7', '-p', 'blastp'], capture_output=True, text=True, encoding='latin-1', cwd=d, env=env)
        extra = {'xmlz': base64.b64encode(zlib.compress(x.stdout.encode('latin-1'), 9)).decode('ascii')}
        if lib is not None: extra['lib'] = lib
        return rows, extra


PLAN = {
    'TMalign_pair': [('-method TMalign_pair', 70), ('-method TMalign_pair,slow_pair', 30), ('-method proba_pair,TMalign_pair', 20), ('-method TMalign_pair -relax_lib 0', 20), ('-method TMalign_pair -tree_mode upgma', 20)],
    'mustang_pair': [('-method mustang_pair', 70), ('-method mustang_pair,slow_pair', 30), ('-method proba_pair,mustang_pair', 20), ('-method mustang_pair,TMalign_pair', 20), ('-method mustang_pair -relax_lib 0', 20)],
    # the installed program has no SAP (sap_pair and -mode 3dcoffee run TMalign instead), no DaliLite and no fugueali (it stops: dali_pair "ERROR: UNSPECIFIED", fugue_pair "The Program fugueali Needed by T-COFFEE Could not be found")
    'sap_pair': [('-method sap_pair', 60), ('-method sap_pair,slow_pair', 30), ('-method proba_pair,sap_pair', 20), ('-method sap_pair -relax_lib 0', 20), ('-method sap_pair -tree_mode upgma', 20)],
    # hand-off image A2-12 only (SAP_REAL=1, /opt/sap/bin/sap): the same sets with the REAL SAP; extra["sap"] = the *.sap_results files of the cache
    'sap_real': [('-method sap_pair', 100), ('-method sap_pair,slow_pair', 30)],
    # -mode 3dcoffee also looks the profile templates of the sequence NAMES up on the network (-netfile qNNN: no network here, ~15 s per sequence): only a few sets
    'sap_3dcoffee': [('-mode 3dcoffee', 8)],
    'dali_pair': [('-method dali_pair', 40), ('-method TMalign_pair,dali_pair', 20)],
    'fugue_pair': [('-method fugue_pair', 40), ('-method mustang_pair,fugue_pair', 20)],
    'blastp_msa': [('-method blastp_msa', 80), ('-method blastp_msa,slow_pair', 30), ('-method proba_pair,blastp_msa', 20), ('-method blastp_msa -relax_lib 0', 20), ('-method blastp_msa -tree_mode upgma', 20)],
}


def job(a):
    feature, flags, idx, seed = a
    rng = random.Random(seed * 1000003 + idx * 7919 + sum(map(ord, feature + flags)) * 104729)
    if feature == 'blastp_msa':
        base = ''.join(rng.choice(PROTEIN) for _ in range(rng.randint(30, 70)))
        seqs = []
        for _ in range(rng.randint(2, 8)):
            if rng.random() < 0.7:
                s = list(base)
                for _ in range(rng.randint(0, max(1, len(base) // 4))):
                    op = rng.choice('sid'); p = rng.randrange(len(s))
                    if op == 's': s[p] = rng.choice(PROTEIN)
                    elif op == 'i': s.insert(p, rng.choice(PROTEIN))
                    elif len(s) > 12: del s[p]
                seqs.append(''.join(s))
            else:
                seqs.append(''.join(rng.choice(PROTEIN) for _ in range(rng.randint(25, 70))))
        rows, extra = run_blast(flags, seqs)
    else:
        seqs, pdbseqs, coords = make_structure_set(rng, (2, 6), (14, 40), 0.0 if '3dcoffee' in flags else 0.15)   # a sequence without a template makes the program look for a PDB entry of that NAME on the network (-netfile)
        rows, extra = run_structure(feature, flags, seqs, pdbseqs, coords)
    return '\t'.join([feature, flags, '|'.join(seqs), rows, json.dumps(extra, separators=(',', ':'))])


def main():
    out = sys.argv[1]
    feats = sys.argv[2:] or list(PLAN)
    scale = float(os.environ.get('SCALE', '1')); seed = int(os.environ.get('SEED', '7'))
    os.makedirs(out, exist_ok=True)
    jobs = [(f, fl, i, seed) for f in feats for fl, n in PLAN[f] for i in range(max(1, int(n * scale)))]
    with ThreadPoolExecutor(max_workers=os.cpu_count() or 4) as ex:
        res = list(ex.map(job, jobs))
    path = os.path.join(out, 'tcoffee_external_oracle.tsv')
    new = not os.path.exists(path)
    with open(path, 'a') as f:
        if new: f.write('feature\tflags\tseqs\trows\textra\n')
        for r in res: f.write(r + '\n')


if __name__ == '__main__':
    main()
