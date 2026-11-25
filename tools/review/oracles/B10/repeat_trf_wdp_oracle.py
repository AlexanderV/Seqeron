import random, json, sys
def wdp_score(s, c, M=2, X=-7, G=-7):
    n, p = len(s), len(c); best = 0
    prev = [0]*p
    for i in range(n):
        cur = [0]*p
        for j in range(p):
            sc = M if s[i] == c[j] else X
            v = max(0, prev[(j-1) % p] + sc, prev[j] + G)
            cur[j] = v
        for _ in range(2):  # close the cyclic deletion (gap in sequence) loop
            for j in range(p):
                cur[j] = max(cur[j], cur[(j-1) % p] + G)
        best = max(best, max(cur)); prev = cur
    return best
def mutate(u, copies, sub, indel, rng):
    out = []
    for _ in range(copies):
        for ch in u:
            r = rng.random()
            if r < sub: out.append(rng.choice([x for x in 'ACGT' if x != ch]))
            elif r < sub + indel/2: pass
            elif r < sub + indel: out.append(ch); out.append(rng.choice('ACGT'))
            else: out.append(ch)
    return ''.join(out)
def gen(rng):
    p = rng.randint(2, 25); u = ''.join(rng.choice('ACGT') for _ in range(p))
    copies = rng.randint(4, 25); sub = rng.choice([0, .03, .08, .15]); indel = rng.choice([0, .02, .05])
    rep = mutate(u, copies, sub, indel, rng)
    left = ''.join(rng.choice('ACGT') for _ in range(rng.randint(10, 80)))
    right = ''.join(rng.choice('ACGT') for _ in range(rng.randint(10, 80)))
    return dict(seqs=[left + rep + right], unit=u, planted=(len(left), len(left)+len(rep)-1))
if __name__ == '__main__':
    rng = random.Random(int(sys.argv[1])); cs = [gen(rng) for _ in range(int(sys.argv[2]))]
    json.dump(cs, open(sys.argv[3], 'w'))
