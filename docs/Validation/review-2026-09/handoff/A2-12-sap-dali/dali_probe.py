#!/usr/bin/env python3
"""Records what DaliLite v5 (/opt/dali/DaliLite.v5/bin) does with the synthetic CA-only structures of generate_tcoffee_external_oracle.py (24 pairs): the import step, the pairwise alignment and every
raw output file. Never raises: every outcome is written to OUT_DIR/dali_probe.tsv (pair <TAB> step <TAB> return code <TAB> first 2000 characters of stdout+stderr)."""
import os, random, subprocess, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_tcoffee_external_oracle as G

out = sys.argv[1]
BIN = '/opt/dali/DaliLite.v5/bin'
rows = []
rng = random.Random(7)
done = 0
while done < 24:
    seqs, pdbseqs, coords = G.make_structure_set(rng, (2, 4), (20, 40), 0.0)
    with tempfile.TemporaryDirectory() as d:
        files = []
        for i, (ps, pts) in enumerate(zip(pdbseqs, coords)):
            p = f'{d}/p{i}.pdb'
            open(p, 'w').write(G.pdb_text(ps, pts))
            files.append(p)
        for i in range(len(files)):
            for j in range(i + 1, len(files)):
                if done >= 24: break
                env = dict(os.environ, PATH=BIN + ':' + os.environ.get('PATH', ''))
                steps = [('import1', [BIN + '/import.pl', '--pdbfile', files[i], '--pdbid', 'aaaa', '--dat', d, '--clean']),
                         ('import2', [BIN + '/import.pl', '--pdbfile', files[j], '--pdbid', 'bbbb', '--dat', d, '--clean']),
                         ('dali', [BIN + '/dali.pl', '--cd1', 'aaaaA', '--cd2', 'bbbbA', '--dat1', d, '--dat2', d, '--outfmt', 'summary,alignments,equivalences', '--clean'])]
                for name, cmd in steps:
                    try:
                        r = subprocess.run(cmd, capture_output=True, text=True, encoding='latin-1', cwd=d, env=env, timeout=120)
                        rows.append(f'{done}\t{name}\t{r.returncode}\t' + (r.stdout + r.stderr)[:2000].replace('\n', '\\n').replace('\t', ' '))
                    except Exception as e:                                  # noqa: BLE001 - the probe records everything
                        rows.append(f'{done}\t{name}\texception\t{e}')
                done += 1
open(os.path.join(out, 'dali_probe.tsv'), 'w').write('pair\tstep\trc\toutput\n' + '\n'.join(rows) + '\n')
