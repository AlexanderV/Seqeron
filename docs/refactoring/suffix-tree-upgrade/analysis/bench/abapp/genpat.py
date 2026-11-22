#!/usr/bin/env python3
"""Pattern files for the in-process A/B client (abapp): lines 'P <present, count>', 'L <present, locate>', 'A <absent, count>'.
usage: genpat.py <text file> <alphabet or @text> <len P> <len L> <len A> <n P> <n L> <n A> <seed> <out>"""
import sys
M = (1 << 64) - 1
class Mix:
    def __init__(self, s): self.s = s & M
    def nxt(self):
        self.s = (self.s + 0x9E3779B97F4A7C15) & M
        z = self.s
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & M
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & M
        return z ^ (z >> 31)
    def below(self, n): return self.nxt() % n
text_path, alpha, lp, ll, la, np_, nl, na, seed, out = sys.argv[1:11]
lp, ll, la, np_, nl, na, seed = map(int, (lp, ll, la, np_, nl, na, seed))
text = open(text_path, encoding='latin-1').read().rstrip('\n')
if alpha == '@text': alpha = ''.join(sorted(set(text)))
r = Mix(seed)
with open(out, 'w') as f:
    for _ in range(np_): f.write('P ' + text[(p := r.below(len(text) - lp)):p + lp] + '\n')
    for _ in range(nl): f.write('L ' + text[(p := r.below(len(text) - ll)):p + ll] + '\n')
    for _ in range(na): f.write('A ' + ''.join(alpha[r.below(len(alpha))] for _ in range(la)) + '\n')
