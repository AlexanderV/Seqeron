import subprocess, os, sys, glob, statistics, json, re
SP=os.environ['SP']
env=dict(os.environ, PATH=os.environ['PATH']+':/root/.dotnet', DOTNET_ROOT='/root/.dotnet', DOTNET_CLI_TELEMETRY_OPTOUT='1')
DLL='/home/user/Seqeron/apps/SuffixTree.Esa.Bench/bin/Release/net10.0/SuffixTree.Esa.Bench.dll'
def libsais(f, reps):
    out=[]; sa=lcp=None
    for i in range(reps):
        o=subprocess.run([SP+'/libsais_work/libsais_drv2',f,'hash'],capture_output=True,text=True).stdout
        m=dict(x.split('=') for x in o.split())
        out.append(float(m['total'])); sa=m['sa_sha']; lcp=m['lcp_sha']
    return out,sa,lcp
def bench(f, lcpk, reps, sak='pointer'):
    res=[]; hs=None
    for i in range(reps):
        o=subprocess.run(['dotnet',DLL,'build','--input',f,'--sa-kernel',sak,'--lcp-kernel',lcpk],capture_output=True,text=True,env=env).stdout
        d={}
        for l in o.splitlines():
            p=l.split('\t')
            if len(p)>=4 and p[0]=='R': d[p[2]]=p[3]
        res.append(d)
    return res
files=sys.argv[2:]
reps=int(sys.argv[1])
rows=[]
print('file\tn\tlibsais_total_med\tlibsais_sa\tdirect_total_med\tdirect_sa\tdirect_lcp\tp3_total_med\tp3_lcp\tref_lcp(1run)\tkind\trss_Bpc\thash_ok')
for f in files:
    name=os.path.basename(f).replace('.txt','')
    lt,sa,lcp=libsais(f,reps)
    d=bench(f,'direct',reps); p=bench(f,'pointer',reps)
    r=bench(f,'reference',1)
    med=lambda L,k: statistics.median(float(x[k]) for x in L)
    ok=all(x['sa_sha256']==sa and x['lcp_sha256']==lcp for x in d+p+r)
    n=d[0]['input_chars']
    print(f"{name}\t{n}\t{statistics.median(lt):.3f}\t-\t{med(d,'build_total_s'):.3f}\t{med(d,'build_sa_s'):.3f}\t{med(d,'build_lcp_s'):.3f}\t{med(p,'build_total_s'):.3f}\t{med(p,'build_lcp_s'):.3f}\t{float(r[0]['build_lcp_s']):.3f}\t{d[0]['lcp_kind']}\t{float(d[0]['peak_rss_bytes_per_char']):.1f}\t{int(ok)}",flush=True)
