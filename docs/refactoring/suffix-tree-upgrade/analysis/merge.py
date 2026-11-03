import csv, glob, sys, collections
W = dict(C1=3,C2=2,C3=3,C4=3,C5=2,C6=2,C7=3,C8=2,C9=2,C10=3,C11=3,C12=3,C13=3,C14=2,C15=2,C16=2,C17=1,C18=2,C19=2,C20=2,C21=2,C22=1,C23=3,C24=1)
CAPS = [f"C{i}" for i in range(1,25)]
def load(paths):
    d = collections.defaultdict(dict)
    for p in paths:
        with open(p, encoding="utf-8") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                d[r["competitor_id"]][r["cap_id"]] = (int(r["score"]), r["flag"], r["evidence"], r["url"])
    return d
paths = [p for p in glob.glob("results/K*.tsv") if not p.endswith(("meta.tsv","survey.tsv"))] + []
d = load(paths)
for k,v in load(["audit/bench_overrides.tsv"]).items():
    d[k].update(v)
names = {"K0":"Seqeron","K1":"sdsl","K2":"SeqAn3","K3":"GenomeTools","K4":"MUMmer4","K5":"libsais","K6":"Go","K7":"rust-bio","K8":"Python","K9":".NET","K10":"BWA/Bt2","K11":"r-index","K12":"ropebwt3"}
ids = ["K0"]+[f"K{i}" for i in range(1,13)]
for k in ids:
    miss = [c for c in CAPS if c not in d[k]]
    if miss: print("MISSING", k, miss)
def wmean(k, ws=W):
    return sum(d[k][c][0]*ws[c] for c in CAPS)/sum(ws.values())
def umean(k): return sum(d[k][c][0] for c in CAPS)/24
print("cap  w  " + " ".join(f"{names[k][:6]:>6}" for k in ids))
for c in CAPS:
    print(f"{c:<4} {W[c]}  " + " ".join(f"{d[k][c][0]}{d[k][c][1]:<1}".rjust(6) for k in ids))
print()
for k in sorted(ids, key=lambda k:-wmean(k)):
    u = sum(1 for c in CAPS if d[k][c][1]=="U")
    print(f"{k:<4}{names[k]:<12} weighted {wmean(k):.2f}  unweighted {umean(k):.2f}  U-cells {u}")
