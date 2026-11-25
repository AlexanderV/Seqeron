"""Compare the C# harness output with wp4_oracle.py expectations: python3 -I compare_wp4.py cases.json got.json"""
import json, sys
sys.path.insert(0, __file__.rsplit("/", 1)[0])
cases = json.load(open(sys.argv[1])); got = json.load(open(sys.argv[2]))
bad = {}
for c, g in zip(cases, got):
    e = c['exp']; op = c['op']
    ok = True
    if op == 'dbg':
        k, rc = c['k'], c['rc']
        import wp4_oracle as wo
        cyc_classes = {(cl, n) for cl, n in e['cycles']}
        rest, got_cyc = [], set()
        for x in g['contigs']:
            if x[:k - 1] == x[len(x) - (k - 1):] and (wo.cycle_class(x, k, rc), len(x)) in cyc_classes:
                got_cyc.add((wo.cycle_class(x, k, rc), len(x)))
            else:
                rest.append(x)
        ok = sorted(rest) == sorted(e['contigs']) and got_cyc == cyc_classes and g['assembled'] == e['assembled'] and g['total'] == e['total']
    elif op == 'correct':
        ok = g['reads'] == e['reads']
    elif op == 'stats':
        ok = all((abs(g[k] - e[k]) < 1e-12 if k == 'gc' else g[k] == e[k]) for k in ('total', 'n50', 'l50', 'n90', 'l90', 'gc', 'gaps', 'gapLen', 'noGaps'))
        ok = ok and g['gapList'] == e['gapList'] and g['contigs'] == e['contigs'] and g['extracted'] == e['extracted']
    elif op == 'aun':
        ok = abs(g['v'] - e['v']) <= 1e-15 * max(1, abs(e['v'])) and (e['quast'] is None or abs(g['v'] - e['quast']) <= 1e-9 * max(1, abs(e['quast'])))
    elif op == 'trim':
        ok = g['kept'] == e['kept']   # minLength 0 keeps every read, in order
    elif op == 'scaffold':
        ok = g['scaffolds'] == e['scaffolds']
    bad.setdefault(op, [0, 0])
    bad[op][0] += 1
    if not ok:
        bad[op][1] += 1
        if bad[op][1] <= 2:
            print('MISMATCH', op, json.dumps({k: v for k, v in c.items() if k != 'exp'})[:300], '\n exp', json.dumps(e)[:300], '\n got', json.dumps(g)[:300])
print({k: f"{v[0]-v[1]}/{v[0]} match" for k, v in bad.items()})
