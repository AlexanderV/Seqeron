# Hand-written mutants of the 1b.4a prototype code (Links/*.cs), one at a time: apply, build the Esa tests, run the Links tests (100 s timeout = hang counts as killed), restore.
# Usage: python3 -u mutation_p1_runner.py [mutant ids]   (paths are those of the session machine). Result of the run: mutation_p1_result.txt
import subprocess, sys, shutil, os
ROOT='/home/user/Seqeron'
L=ROOT+'/apps/SuffixTree.Esa.Bench/Links/'
M=[
 ('H1','LcpRangeHelper.cs','            if (_lcp[j] < bound) return j;\n        }\n\n        int idx = pos >> _logFanout;\n        for (int k = 1; k <= _levels.Length; k++)\n        {\n            int[] level = _levels[k - 1];\n            for (int b = idx - 1','            if (_lcp[j] <= bound) return j;\n        }\n\n        int idx = pos >> _logFanout;\n        for (int k = 1; k <= _levels.Length; k++)\n        {\n            int[] level = _levels[k - 1];\n            for (int b = idx - 1'),
 ('H3','LcpRangeHelper.cs','for (int j = hi & ~(f - 1); j <= hi; j++)','for (int j = (hi & ~(f - 1)) + 1; j <= hi; j++)'),
 ('H4','LcpRangeHelper.cs','return NextLess(lo, min + 1, ref cnt);','return NextLess(lo, min, ref cnt);'),
 ('H5','LcpRangeHelper.cs','        int position = b * f;\n        while (true)','        int position = (b * f) + 1;\n        while (true)'),
 ('H6','LcpRangeHelper.cs','            l = (l >> _logFanout) + 1;','            l = (l >> _logFanout);'),
 ('H7','LcpRangeHelper.cs','if (level[b] < bound) return RightmostBelow(k, b, bound, ref cnt);','if (level[b] <= bound) return RightmostBelow(k, b, bound, ref cnt);'),
 ('K1','LinkKernels.cs','            if (lcp[l] < t) break;\n            l--;','            if (lcp[l] <= t) break;\n            l--;'),
 ('K2','LinkKernels.cs','        int r = b + 1;\n        while (r <= n)','        int r = b;\n        while (r <= n)'),
 ('K3','LinkKernels.cs','while (boundary >= 0 && boundary <= a)','while (boundary >= 0 && boundary < a)'),
 ('K4','LinkKernels.cs','h.ArgMin(a + 1, b, ref cnt)','h.ArgMin(a, b, ref cnt)'),
 ('K5','LinkKernels.cs','h.NextLess(b + 1, t, ref cnt)','h.NextLess(b, t, ref cnt)'),
 ('K6','LinkKernels.cs','return new LinkTarget(l, r - 1, Math.Max(left, right), left > right);','return new LinkTarget(l, r - 1, Math.Max(left, right), left >= right);'),
 ('K7','LinkKernels.cs','c.Isa![c.Sa[v.Rb] + 1]','c.Isa![c.Sa[v.Rb]]'),
 ('K8','LinkKernels.cs','target = new LinkTarget(cl, cr, ell, cLast);','target = new LinkTarget(cl, cr, ell, !cLast);'),
 ('D1','MemDriver.cs','_text[r - 1] != qsym[start - 1]','_text[r - 1] == qsym[start - 1]'),
 ('D2','MemDriver.cs','for (int k = hi + 1; k <= _n; k++)','for (int k = hi + 1; k < _n; k++)'),
 ('D3','MemDriver.cs','Depth = l.V.Depth - 1, Last','Depth = l.V.Depth, Last'),
 ('D4','MemDriver.cs','            if (l.Len == l.C.Depth)\n            {\n                l.V = l.C;\n                l.HasEdge = false;\n            }\n            return true;','            if (l.Len > l.C.Depth)\n            {\n                l.V = l.C;\n                l.HasEdge = false;\n            }\n            return true;'),
 ('D5','MemDriver.cs','            if (v < minLength) break;\n            if (v < run) run = v;\n            k--;','            if (v <= minLength) break;\n            if (v < run) run = v;\n            k--;'),
 ('D6','MemDriver.cs','if (r == 0 || start == 0 ||','if (start == 0 ||'),
 ('N1','LinkNodes.cs','nodeDepthAll[emissionOfPush[eParentPush[e]]] + 1','nodeDepthAll[emissionOfPush[eParentPush[e]]] + 2'),
 ('T1','LinkTables.cs','psi[r] = isa[sa[r] + 1];','psi[r] = isa[sa[r]];'),
]
only=sys.argv[1:] 
env=dict(os.environ, PATH=os.environ['PATH']+':/root/.dotnet', DOTNET_ROOT='/root/.dotnet', DOTNET_CLI_TELEMETRY_OPTOUT='1', SEQERON_BENCH_DIR='/tmp/stbench')
for mid,f,old,new in M:
    if only and mid not in only: continue
    path=L+f
    src=open(path).read()
    if src.count(old)!=1:
        print(mid,'APPLY-FAIL count',src.count(old)); continue
    open(path,'w').write(src.replace(old,new))
    try:
        b=subprocess.run(['dotnet','build','tests/SuffixTree/SuffixTree.Esa.Tests','-c','Release','--nologo','-v','q'],cwd=ROOT,env=env,capture_output=True,text=True)
        if b.returncode!=0:
            print(mid,'BUILD-FAIL',(b.stdout+b.stderr)[-300:].replace('\n',' ')); continue
        try:
            t=subprocess.run(['dotnet','test','tests/SuffixTree/SuffixTree.Esa.Tests','-c','Release','--nologo','--no-build','--filter','FullyQualifiedName~SuffixTree.Esa.Tests.Links'],cwd=ROOT,env=env,capture_output=True,text=True,timeout=100)
        except subprocess.TimeoutExpired:
            subprocess.run(['pkill','-f','SuffixTree.Esa.Tests.dll'])
            print(mid,'KILLED (hang: the tests did not finish in 100 s)',flush=True)
            continue
        out=t.stdout
        res='KILLED' if t.returncode!=0 else 'SURVIVED'
        failed=[l.strip() for l in out.splitlines() if l.strip().startswith('Failed ')]
        print(mid,res,len(failed),'failing:',', '.join(x.split()[1] for x in failed[:3]),flush=True)
    finally:
        open(path,'w').write(src)
sys.stdout.flush()
