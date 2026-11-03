import csv, glob, collections, re
d = collections.defaultdict(dict)
for p in glob.glob("func/K*.tsv")+["audit/K0_func.tsv"]:
    if p.endswith(".extra.tsv"): continue
    rd=csv.reader(open(p,encoding="utf-8"),delimiter="\t"); next(rd)
    for r in rd:
        d[r[0]][r[1]] = (int(r[2]), r[3], r[4], r[5])
names = {"K0":"Seqeron","K1":"sdsl","K2":"SeqAn3","K3":"GenomeT","K4":"MUMmer4","K5":"libsais","K6":"Go","K7":"rust-bio","K8":"Python","K9":".NET","K10":"BWA/Bt2","K11":"r-index","K12":"ropebwt3"}
ids = list(names)
F = [f"F{i:02d}" for i in range(1,82)]
for k in ids:
    miss=[f for f in F if f not in d[k]]
    if miss: print("MISSING",k,miss)
def cov(k, mode="norm"):
    s=0
    for f in F:
        v=d[k][f][0]
        if mode=="strict": v = 2 if v==2 else 0
        if mode=="lenient": v = 2 if v>=1 else 0
        s+=v
    return s/162
print("tool        norm  strict lenient  U-cells")
for k in sorted(ids,key=lambda k:-cov(k)):
    print(f"{k:<4}{names[k]:<9} {cov(k)*100:5.1f}  {cov(k,'strict')*100:5.1f}  {cov(k,'lenient')*100:5.1f}   {sum(1 for f in F if d[k][f][1]=='U')}")
groups = {"A Construction":F[0:11],"B Exact":F[11:24],"C Patterns":F[24:30],"D Approx":F[30:35],"E Matching":F[35:44],"F Repeats/stats":F[44:53],"G Navigation":F[53:65],"H Persist/ops":F[65:75],"I Alphabets":F[75:81]}
print("\ngroup coverage % (Seqeron | best other (who))")
for g,fs in groups.items():
    def gc(k): return sum(d[k][f][0] for f in fs)/(2*len(fs))*100
    best=max((k for k in ids if k!="K0"),key=gc)
    print(f"{g:<17} {gc('K0'):5.1f} | {gc(best):5.1f} ({names[best]})")
print("\nFunctions where Seqeron < best competitor (gap = best-Seq; n2 = #competitors scoring 2; n12 = #scoring>=1)")
rows=[]
for f in F:
    s=d["K0"][f][0]
    others={k:d[k][f][0] for k in ids if k!="K0"}
    best=max(others.values())
    if best>s:
        rows.append((best-s, sum(1 for v in others.values() if v==2), sum(1 for v in others.values() if v>=1), f, s, best))
rows.sort(key=lambda r:(-r[0],-r[1],-r[2]))
func_names = {}
for line in open("FUNCTIONS.md",encoding="utf-8"):
    m=re.match(r"(F\d\d) (.*)",line.strip())
    if m: func_names[m.group(1)]=m.group(2)
for gap,n2,n12,f,s,b in rows:
    print(f"{f} seq={s} best={b} n2={n2:2d} n>=1={n12:2d}  {func_names[f]}")
print("\nSeqeron ahead (score > every competitor):", [f for f in F if d['K0'][f][0]>max(d[k][f][0] for k in ids if k!='K0')])
