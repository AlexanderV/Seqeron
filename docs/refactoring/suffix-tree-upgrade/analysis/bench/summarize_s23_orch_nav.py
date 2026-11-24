# summary of run_s23_orch_nav_ab.sh rows (arms base, base2 = null, new): quantiles of the per-process times and the paired ratios
import re,sys,statistics as st
def q(xs,f): xs=sorted(xs); return xs[int(f*(len(xs)-1))]
a={};s={}
for line in open(sys.argv[1]):
    p=line.rstrip('\n').split('\t'); t='\t'.join(p[2:])
    m=re.search(r'warm_median_s\t[0-9.]+\ts\tmedian of \d+ warm calls \([^)]*\) \[([0-9.]+) \.\. ([0-9.]+)\]',t)
    m2=re.search(r'steady_median_s\t([0-9.]+)\ts',t)
    if m: a.setdefault(p[1],{})[int(p[0])]=float(m.group(1))
    if m2: s.setdefault(p[1],{})[int(p[0])]=float(m2.group(1))
for name,d in (('min over warm calls',a),('steady median',s)):
    print(name,{k:len(v) for k,v in d.items()})
    for k in ('base','base2','new'): v=list(d[k].values()); print(f"  {k:6s} min={min(v):.4f} q25={q(v,.25):.4f} med={st.median(v):.4f} q75={q(v,.75):.4f}")
    for f in (0,.1,.25,.5):
        print(f"  q{f}: new/base={q(list(d['new'].values()),f)/q(list(d['base'].values()),f):.3f}  null base2/base={q(list(d['base2'].values()),f)/q(list(d['base'].values()),f):.3f}")
    r=[d['new'][i]/d['base'][i] for i in d['base'] if i in d['new']]; r0=[d['base2'][i]/d['base'][i] for i in d['base'] if i in d['base2']]
    print('  paired median new/base %.3f (q25 %.3f q75 %.3f)  null %.3f (q25 %.3f q75 %.3f)'%(st.median(r),q(r,.25),q(r,.75),st.median(r0),q(r0,.25),q(r0,.75)))
