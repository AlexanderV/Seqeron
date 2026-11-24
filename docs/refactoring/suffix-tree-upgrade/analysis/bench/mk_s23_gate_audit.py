# builds results_s23_gate_rows_audit.txt from the raw rows of sessions 7 - 7e (scratch directory /tmp/stbench/s23; the raw rows are not kept: the file is the record)
import statistics as st, collections, os, subprocess, sys
B='/tmp/stbench/s23'
out=[]
out.append("# S2.3 gate rows after the audit (nav, query): the S1 head (base) against the head of the first audit commit 646f78c0 (new2), the head of the build before it a563e755 (new) and the head after the out-of-line bridges 0bfb2114 (new3). Fresh processes of run_1b9_gate.sh, interleaved by round.")
out.append("# The libsais drift control ran before and after every session (results_s23_drift.tsv): 1.0 - 1.3 of the reference 0.328 s, i.e. the machine was slower than at the S1 head in every session; the variants are interleaved, so the ratio is paired, but the rounds of one session differ by +-12 % (q25-q75), and the same build read 1.00 and 1.15 in two sessions: the figures below are good to about +-7 %.")
out.append("# sessions: 7 (20 rounds: base, new2: query dna_10m, nav mem, nav anchors), 7b (30 rounds: anchors; base, new, new2), 7c (30: anchors; base, new, new3), 7d (60: anchors, the order rotated every round; base, new, new2, new3), 7e (40: anchors and mem, rotated; base, new, new3)")
out.append("")
out.append("## session 7: base against new2 (20 rounds), every metric of query dna_10m and nav dna_10m")
txt=subprocess.run(['python3',f'{B}/agg_gate.py',f'{B}/gate_base/raw7.tsv',f'{B}/gate_new2/raw7.tsv'],capture_output=True,text=True).stdout
out.append(txt.rstrip())
out.append("")
out.append("## the nav steady and first-call figures per session: median of the variant / median of base (rounds in brackets)")
out.append("session\tmetric\tvariant\tvariant/base\trounds")
def load(v,s):
    p=f'{B}/gate_{v}/raw{s}.tsv'
    d=collections.defaultdict(list)
    if not os.path.exists(p): return None
    for l in open(p):
        x=l.rstrip('\n').split('\t')
        if len(x)>6 and x[4]=='nav_esa' and x[5] in ('anchors_steady_median_s','mem_steady_median_s','anchors_first_call_s','mem_first_call_s'):
            d[x[5]].append(float(x[6]))
    return d
sess={'7':('base','new2'),'7b':('base','new','new2'),'7c':('base','new','new3'),'7d':('base','new','new2','new3'),'7e':('base','new','new3')}
tot=collections.defaultdict(list)
for s in ('7','7b','7c','7d','7e'):
    data={v:load(v,s) for v in sess[s]}
    for m in ('anchors_steady_median_s','mem_steady_median_s','anchors_first_call_s','mem_first_call_s'):
        if m not in data['base']: continue
        b=st.median(data['base'][m])
        for v in sess[s][1:]:
            r=st.median(data[v][m])/b
            out.append(f"{s}\t{m}\t{v}\t{r:.3f}\t{len(data[v][m])}")
            tot[(v,m)].append((r,len(data[v][m])))
out.append("")
out.append("## all sessions: the median of the session ratios and the ratio weighted by rounds")
out.append("variant\tmetric\tsessions\tmedian of the ratios\tweighted by rounds\tmin\tmax")
for (v,m),l in sorted(tot.items()):
    r=[x for x,_ in l]; n=[k for _,k in l]
    out.append(f"{v}\t{m}\t{len(l)}\t{st.median(r):.3f}\t{sum(x*k for x,k in l)/sum(n):.3f}\t{min(r):.3f}\t{max(r):.3f}")
open(sys.argv[1],'w').write('\n'.join(out)+'\n')
print(len(out),'lines')
