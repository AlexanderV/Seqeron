#!/usr/bin/env python3
"""One heavy comparison of the installed `mafft` 7.505-1 with the port on a case in which `nodepair` DIVIDES a node into jobs (B09 WP73, A2-19c-s4c-lhlimit): the node must have n0 * n1 * lenav^2 > MINBLOCKLEN2 (1e9; mltaln9.c
recalcpairs_para4), which needs two groups of 36 sequences of 1350 nucleotides; `--lhlimit 20` keeps the pairs of only the longest sequences so that the real run takes about a minute.  The case is not in the test
suite (the port needs 80 seconds on it); the result is recorded in docs/Validation/review-2026-09/B09.md (F336).

Usage: python3 scripts/alignment/verify_mafft_lhlimit_divided.py   (needs `mafft` 7.505 on the PATH and the built tests project: `dotnet build Seqeron.sln -c Debug`)
Prints the three lines of the program's standard error that name the division ("numjob=..., nnodesdivided=...", "Divided n nodes to m jobs") and whether return code, standard output, standard error and files of the real program and of the port are equal.
"""
import json, os, random, subprocess, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mafft_command_oracle_lib as L

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
HARNESS = '''using Seqeron.Genomics.Alignment;
var d = args[0];
var files = new Dictionary<string, string>();
foreach (var f in Directory.GetFiles(d)) files[Path.GetFileName(f)] = File.ReadAllText(f, System.Text.Encoding.Latin1);
var r = MafftProgram.Run(args.Skip(1).ToArray(), files);
Console.Out.Write(System.Text.Json.JsonSerializer.Serialize(new { rc = r.ExitCode, outp = r.Stdout, err = r.Stderr, files = r.Files }));
'''


def main():
    rng = random.Random(4)
    n, ln = 36, 1350
    b1 = ''.join(rng.choice(L.N) for _ in range(ln))
    b2 = L.mut(rng, b1, L.N, 0.5, 0.0)
    seqs = [L.mut(rng, b1, L.N, 0.03, 0.0) for _ in range(n)] + [L.mut(rng, b2, L.N, 0.03, 0.0) for _ in range(n)]
    text = L.fasta(seqs)
    args = ['--large', '--globalpair', '--lhlimit', '20']
    real = L.run_once(args, text)
    for line in real[2].split('\n'):
        if 'Divided' in line or 'numjob' in line or 'sizeav' in line:
            print(line)
    with tempfile.TemporaryDirectory() as td:
        proj = os.path.join(td, 'h')
        os.makedirs(proj)
        open(os.path.join(proj, 'Program.cs'), 'w').write(HARNESS)
        open(os.path.join(proj, 'h.csproj'), 'w').write('<Project Sdk="Microsoft.NET.Sdk"><PropertyGroup><OutputType>Exe</OutputType><TargetFramework>net10.0</TargetFramework><ImplicitUsings>enable</ImplicitUsings><Nullable>enable</Nullable></PropertyGroup>'
                                                         '<ItemGroup><ProjectReference Include="%s" /></ItemGroup></Project>' % os.path.join(ROOT, 'src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/Seqeron.Genomics.Alignment.csproj'))
        subprocess.run(['dotnet', 'build', '-c', 'Release', '-o', os.path.join(proj, 'out'), proj], check=True, capture_output=True)
        data = os.path.join(td, 'data')
        os.makedirs(data)
        open(os.path.join(data, L.INPUT_NAME), 'w').write(text)
        p = subprocess.run(['dotnet', os.path.join(proj, 'out', 'h.dll'), data] + args + [L.INPUT_NAME], capture_output=True, text=True, check=True)
        mine = json.loads(p.stdout)
    print('equal: rc', real[0] == mine['rc'], 'stdout', real[1] == mine['outp'], 'stderr', real[2] == mine['err'], 'files', real[3] == mine['files'])


if __name__ == '__main__':
    main()
