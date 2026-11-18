#!/usr/bin/env python3
"""Round 3 of the verifier audit: the probe scheme of long pairs changed after round 2 (tail probes), the width-rule function was renamed, and the stack-bound option gained an alternative: the mutants of those sites against the final code.
   usage: make_e8_round3.py <repo root>   (writes e8_round3_outputverifier_listed.json and e8_round3_e8_listed.json)"""
import json
V = 'apps/SuffixTree.EmProto/Verify/OutputVerifier.cs'
O = []
def m(id, desc, old, new):
    O.append({'id': id, 'desc': desc, 'old': old, 'new': new, 'file': V})
m('C18', 'no probes of the matching part of a long pair', 'for (int k = 0; k < o.Probes && span > 0; k++)', 'for (int k = 0; k < 0 && span > 0; k++)')
m('C19', 'the tail probes start two symbols before the end of the claimed match', 'long off = k < o.Probes / 2 ? stored - 1 - k :', 'long off = k < o.Probes / 2 ? stored - 2 - k :')
m('C19b', 'no tail probes (random probes only)', 'long off = k < o.Probes / 2 ? stored - 1 - k :', 'long off = k < 0 ? stored - 1 - k :')
m('C19c', 'a third of the probes are tail probes', 'long off = k < o.Probes / 2 ? stored - 1 - k :', 'long off = k < o.Probes / 3 ? stored - 1 - k :')
m('C19d', 'the random probes are all at the first offset above the cap', '% (ulong)span);', '% 1UL);')
m('C23', 'a pair left over in a block is not reported', 'if (lcp.HasEscape && lcp.PairCount > 0 && pair != lcp.LowerBound(end)) totals.Lcp.Add(', 'if (lcp.HasEscape && lcp.PairCount > 0 && pair < 0) totals.Lcp.Add(')
m('C24', 'the width rule takes 17 times the escapes', 'if (16 * escapes8 <= n) return 1;', 'if (17 * escapes8 <= n) return 1;')
m('C25', 'the width rule for U16 takes 9 times the escapes', 'if (8 * escapes16 <= n) return 2;', 'if (9 * escapes16 <= n) return 2;')
m('C26', 'the width rule takes 15 times the escapes', 'if (16 * escapes8 <= n) return 1;', 'if (15 * escapes8 <= n) return 1;')
m('C27', 'the first-probe loop condition ignores the span', 'k < o.Probes && span > 0; k++)', 'k < o.Probes; k++)')
json.dump(O, open('e8_round3_outputverifier_listed.json', 'w'), indent=1)
E = [{'id': 'H05', 'desc': 'the stack cap is given in entries instead of Mi entries', 'old': ': o.Long("stack-cap-m", 64) << 20,', 'new': ': o.Long("stack-cap-m", 64),', 'file': 'apps/SuffixTree.EmProto/Experiments/E8.cs'},
     {'id': 'H08', 'desc': 'the default number of probes is 0', 'old': 'Probes = o.Int("probes", 16),', 'new': 'Probes = o.Int("probes", 0),', 'file': 'apps/SuffixTree.EmProto/Experiments/E8.cs'}]
json.dump(E, open('e8_round3_e8_listed.json', 'w'), indent=1)
print(len(O) + len(E), 'mutants for round 3')
