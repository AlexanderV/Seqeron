"""Projection model for the suffix-tree upgrade plan.
Inputs: results/K*.tsv + audit/K0_audit.tsv (property matrix, 0-4, weights W), func/K*.tsv + audit/K0_func.tsv (function matrix, 0-2).
Stages set TARGET cells for Seqeron (max with current). A target is reached only if the stage's acceptance gate passes.
"""
import csv, glob, collections, itertools, json, sys

W = dict(C1=3,C2=2,C3=3,C4=3,C5=2,C6=2,C7=3,C8=2,C9=2,C10=3,C11=3,C12=3,C13=3,C14=2,C15=2,C16=2,C17=1,C18=2,C19=2,C20=2,C21=2,C22=1,C23=3,C24=1)
CAPS = list(W)
F = [f"F{i:02d}" for i in range(1, 82)]
NAMES = {"K0":"Seqeron","K1":"sdsl","K2":"SeqAn3","K3":"GenomeTools","K4":"MUMmer4","K5":"libsais","K6":"Go","K7":"rust-bio","K8":"Python","K9":".NET","K10":"BWA/Bt2","K11":"r-index","K12":"ropebwt3"}
IDS = list(NAMES)

def load_c():
    d = collections.defaultdict(dict)
    for p in glob.glob("results/K*.tsv") + ["audit/K0_audit.tsv"]:
        if p.endswith(("meta.tsv", "survey.tsv")): continue
        rd = csv.reader(open(p, encoding="utf-8"), delimiter="\t"); next(rd)
        for r in rd: d[r[0]][r[1]] = int(r[2])
    rd = csv.reader(open("audit/bench_overrides.tsv", encoding="utf-8"), delimiter="\t"); next(rd)
    for r in rd: d[r[0]][r[1]] = int(r[2])   # measured cells win
    return d

def load_f():
    d = collections.defaultdict(dict)
    for p in glob.glob("func/K*.tsv") + ["audit/K0_func.tsv"]:
        if p.endswith(".extra.tsv"): continue
        rd = csv.reader(open(p, encoding="utf-8"), delimiter="\t"); next(rd)
        for r in rd: d[r[0]][r[1]] = int(r[2])
    return d

def cmean(d, k): return sum(d[k][c] * W[c] for c in CAPS) / sum(W.values())
def fcov(d, k): return sum(d[k][f] for f in F) / 162

# ---- stage targets (id -> {cell: target}); C = property cells, F = function cells -------------------------------
# Values are what a strict auditor would grant IF the stage gate passes. Dependencies: S2,S3,S4,S5,S7 need S1.
STAGES = {
 "S0": dict(name="Safety, honest baseline, perf harness", cost=1, deps=[],
    C={"C23": 3}, F={"F69": 2}),
 "S1": dict(name="Enhanced-suffix-array backend (SA-IS + LCP + child table) behind ISuffixTree", cost=8, deps=[],
    C={"C10": 2, "C11": 2, "C12": 3, "C8": 3},
    F={"F07":2,"F08":2,"F09":2,"F10":2,"F17":2,"F54":2,"F55":2,"F56":2,"F57":1,"F58":2,"F59":2,"F60":1,"F61":2,"F62":2,"F63":2,"F65":2,"F43":2,"F44":1,"F53":2,"F52":2,"F19":2,"F15":2}),
 "S2": dict(name="Genome scale: int64 positions, 2-bit DNA, v7 format, file/stream build", cost=3, deps=["S1"],
    C={"C13": 3, "C19": 3},
    F={"F79":2,"F81":2,"F78":2,"F03":2}),
 "S3": dict(name="Both strands, ms[], SMEM/seeds, collections, batch queries, public chaining", cost=3, deps=["S1"],
    C={"C5": 3, "C6": 3, "C7": 3},
    F={"F02":2,"F22":2,"F23":2,"F29":2,"F35":2,"F36":2,"F39":2,"F40":2,"F50":2,"F75":2}),
 "S4": dict(name="Approximate / degenerate search on the index (mismatch, edit, IUPAC, N, bidirectional)", cost=3, deps=["S1"],
    C={"C2": 3, "C3": 3},
    F={"F31":2,"F32":2,"F33":2,"F34":2,"F18":2,"F20":2}),
 "S5": dict(name="Parallel construction (SA-IS / LCP / child table)", cost=3, deps=["S1"],
    C={"C15": 2, "C10": 3}, F={}),
 "S6": dict(name="Online append, merge, segment delete", cost=1, deps=[],
    C={"C17": 2}, F={"F04":2,"F05":1,"F06":1}),
 "S7": dict(name="Compressed FM-index (BWT + rank/select + sampled locate)", cost=8, deps=["S1"],
    C={"C18": 3, "C11": 4}, F={"F11":2,"F18":2,"F20":2}),
 "S8": dict(name="Packaging and proof: NuGet, dotnet tool, MCP sessions, published benchmarks", cost=1, deps=[],
    C={"C22": 3, "C24": 3, "C23": 4}, F={"F70":2,"F71":2}),
}

def apply(base_c, base_f, chosen):
    c = dict(base_c); f = dict(base_f)
    for s in chosen:
        for k, v in STAGES[s]["C"].items(): c[k] = max(c[k], v)
        for k, v in STAGES[s]["F"].items(): f[k] = max(f[k], v)
    return c, f

def valid(chosen):
    return all(all(dep in chosen for dep in STAGES[s]["deps"]) for s in chosen)

if __name__ == "__main__":
    dc, df = load_c(), load_f()
    bc, bf = dc["K0"], df["K0"]
    best_c = max((k for k in IDS if k != "K0"), key=lambda k: cmean(dc, k))
    best_f = max((k for k in IDS if k != "K0"), key=lambda k: fcov(df, k))
    tc = cmean(dc, best_c); tf = fcov(df, best_f)
    print(f"best competitor property mean: {NAMES[best_c]} {tc:.2f}; function coverage: {NAMES[best_f]} {tf*100:.1f}%")
    print(f"Seqeron now: property {cmean(dc,'K0'):.2f}, function {fcov(df,'K0')*100:.1f}%\n")
    order = ["S0","S1","S2","S3","S4","S5","S6","S7","S8"]
    chosen = []
    print("cumulative (in numeric order)")
    for s in order:
        chosen.append(s)
        c, f = apply(bc, bf, chosen)
        pm = sum(c[k] * W[k] for k in CAPS) / sum(W.values()); fc = sum(f[x] for x in F) / 162
        print(f"+{s} {STAGES[s]['name'][:58]:<58} cost {STAGES[s]['cost']:>2}  prop {pm:.2f}  func {fc*100:5.1f}%")
    print("\nstand-alone gain of each stage (given S1 where it depends on it)")
    for s in order:
        sel = [s] + [d for d in STAGES[s]["deps"]]
        c, f = apply(bc, bf, sel)
        pm = sum(c[k] * W[k] for k in CAPS) / sum(W.values()); fc = sum(f[x] for x in F) / 162
        c0, f0 = apply(bc, bf, STAGES[s]["deps"])
        pm0 = sum(c0[k] * W[k] for k in CAPS) / sum(W.values()); fc0 = sum(f0[x] for x in F) / 162
        print(f"{s} dprop {pm-pm0:+.2f}  dfunc {100*(fc-fc0):+5.1f} pp  cost {STAGES[s]['cost']}")
    print("\ncheapest valid subsets beating the best competitor by margin (prop +0.4, func +15 pp)")
    res = []
    for r in range(1, len(order) + 1):
        for sub in itertools.combinations(order, r):
            if not valid(sub): continue
            c, f = apply(bc, bf, sub)
            pm = sum(c[k] * W[k] for k in CAPS) / sum(W.values()); fc = sum(f[x] for x in F) / 162
            if pm >= tc + 0.4 and fc >= tf + 0.15:
                res.append((sum(STAGES[s]["cost"] for s in sub), sub, pm, fc))
    res.sort()
    for cost, sub, pm, fc in res[:6]:
        print(cost, "+".join(sub), f"prop {pm:.2f} func {fc*100:.1f}%")

    # ---- robustness to the weights: 2000 random weightings (weights 1..3 per capability, uniformly) ----
    import random
    random.seed(7)
    def wins(c, rnd_w):
        me = sum(c[k]*rnd_w[k] for k in CAPS)/sum(rnd_w.values())
        others = max(sum(dc[o][k]*rnd_w[k] for k in CAPS)/sum(rnd_w.values()) for o in IDS if o != "K0")
        return me > others, me - others
    for label, sel in (("now", []), ("S0+S1+S3+S4+S5+S8", ["S0","S1","S3","S4","S5","S8"]), ("all but S7", ["S0","S1","S2","S3","S4","S5","S6","S8"]), ("all", order)):
        c, f = apply(bc, bf, sel)
        ok = 0; gaps = []
        for _ in range(2000):
            rw = {k: random.choice([1,2,3]) for k in CAPS}
            w, g = wins(c, rw); ok += w; gaps.append(g)
        gaps.sort()
        print(f"weights robustness [{label}]: Seqeron first in {ok/20:.1f}% of 2000 weightings; worst gap {gaps[0]:+.2f}, median {gaps[1000]:+.2f}")
