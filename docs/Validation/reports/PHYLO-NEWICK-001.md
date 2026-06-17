# Validation Report: PHYLO-NEWICK-001 — Newick format I/O (N-ary / multifurcating refactor)

- **Validated:** 2026-06-17   **Area:** Phylogenetics
- **Canonical method(s):** `PhylogeneticAnalyzer.ParseNewick(string)`, `PhylogeneticAnalyzer.ToNewick(PhyloNode, bool)`, `PhyloNode.Children`/`Left`/`Right` (`src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs`)
- **Refactor under test:** commit c4f0190 — `PhyloNode` now stores an ordered `List<PhyloNode> Children`; `Left`/`Right` are convenience accessors over the first two; Newick parses & writes genuine multifurcations `(A,B,C);`; the previous throw-on->2-children guard removed.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN (one code defect found and completely fixed this session)

## Stage A — Description (the grammar)

### Sources opened this session

1. **Wikipedia, "Newick format"** — <https://en.wikipedia.org/wiki/Newick_format>
   BNF productions (verbatim):
   - `Tree → Subtree ";"`
   - `Subtree → Leaf | Internal`
   - `Leaf → Name`
   - `Internal → "(" BranchSet ")" Name`
   - `BranchSet → Branch | Branch "," BranchSet`
   - `Branch → Subtree Length`
   - `Name → empty | string`
   - `Length → empty | ":" number`

   The recursive `BranchSet → Branch | Branch "," BranchSet` permits an arbitrary-length
   comma-separated list of branches under one internal node. A node with **≥3** children is a
   **multifurcation / polytomy**. `(A,B,C);` = a single internal node (the root) with three leaf
   children A, B, C.

2. **Felsenstein / Olsen, PHYLIP "The Newick tree format"** — <https://phylipweb.github.io/phylip/newicktree.html>
   Verbatim: *"Interior nodes are represented by a pair of matched parentheses. Between them are
   representations of the nodes that are immediately descended from that node, separated by commas."*;
   *"Trees can be multifurcating at any level."*; *"Branch lengths can be incorporated … by putting a
   real number … after a node and preceded by a colon."*; *"The tree ends with a semicolon."*
   Worked multifurcating example: `(B:6.0,(A:5.0,C:3.0,E:4.0):5.0,D:11.0);` — an interior node with
   **three** immediate descendants (B, an interior node, D).

### Grammar checks

- **Multifurcation is a first-class grammar feature** — yes, via the recursive `BranchSet` (Wikipedia)
  and "multifurcating at any level" (Felsenstein/Olsen). Not an extension; it is the canonical format.
- **Closing parenthesis is mandatory** — `Internal → "(" BranchSet ")" Name`. Every opening `(` must
  be matched by a `)`. Unbalanced parentheses are not a valid `Tree`.
- **Branch length attachment** — `Length → ":" number`, attached per branch; decimal point uses `.`
  (locale-independent).
- **Terminator** — `Tree → Subtree ";"`.

### Hand-derived expected results (from the grammar, not the code)

| Input | Parse (derived) | Serialize (no lengths) |
|-------|-----------------|------------------------|
| `(A,B,C);` | root, 3 leaf children [A,B,C] in order | `(A,B,C);` |
| `(A,B,(C,D,E));` | root 3 children: A, B, internal{C,D,E}; 5 leaves | `(A,B,(C,D,E));` |
| `(A:0.1,B:0.2,C:0.3);` | 3 children, BL 0.1/0.2/0.3 | `(A:0.1000,B:0.2000,C:0.3000);` (with lengths, F4) |
| `(A,B,C,D,E,F,G);` | 7 children | `(A,B,C,D,E,F,G);` |
| `(B:6.0,(A:5.0,C:3.0,E:4.0):5.0,D:11.0);` | root 3 children: B, internal{A,C,E}, D; 5 leaves | `(B,(A,C,E),D);` |
| `(A,B);` (binary) | root, 2 children A,B | `(A,B);` |

**Stage A verdict: PASS** — grammar productions, multifurcation definition, branch-length syntax,
and the mandatory closing parenthesis are confirmed against two authoritative sources.

## Stage B — Implementation

### Code path reviewed

- `ParseNewick` (`PhylogeneticAnalyzer.cs:632`) → `ParseNewickRecursive` (`:662`): on `(`, loops
  `BranchSet` parsing each child + optional `:length`, separated by `,`; then handled `)`; then
  optional internal-node label. `ToNewick`/`ToNewickRecursive` (`:566`/`:576`): emits every child
  comma-separated regardless of count (genuine N-ary write), F4 branch lengths, valid unquoted labels.

### Independent verification (executed probe + tests)

A standalone probe parsed each grammar case and printed child count / leaf list / re-serialization.
All multifurcation cases matched the hand-derived table exactly:

```
(A,B,C);                                -> children=3 leaves=[A,B,C]     write=(A,B,C);
(A,B,(C,D,E));                          -> children=3 leaves=[A,B,C,D,E] write=(A,B,(C,D,E));
(A:0.1,B:0.2,C:0.3);                    -> children=3 (BL 0.1/0.2/0.3)   write(F4)=(A:0.1000,B:0.2000,C:0.3000);
(A,B,C,D,E,F,G);                        -> children=7 leaves=[A..G]      write=(A,B,C,D,E,F,G);
(B:6.0,(A:5.0,C:3.0,E:4.0):5.0,D:11.0); -> children=3 leaves=[B,A,C,E,D] write=(B,(A,C,E),D);
```

- **No 3rd-child dropping.** The `ParseNewick_MultifurcatingRoot_ParsesThreeChildren` test asserts
  `Children.Count == 3` **and** `Children[0..2].Name == A/B/C` **and** the ordered leaf set — a
  regression that dropped/truncated the 3rd child (the old behaviour threw) would fail it. Confirmed.
- **Binary trees unchanged.** `(A,B);`, `(A:0.1,B:0.2);`, `((A,B),(C,D));`, `A;`, and the
  full-format `((A:0.1,B:0.2)AB:0.3,(C:0.4,D:0.5)CD:0.6)Root;` round-trip exactly as before
  (regression tests green).

### Defect found and fixed

**DEFECT (code, silent acceptance of malformed input):** `ParseNewickRecursive` never required the
matching `)`. The closing-paren handling was `if (pos < len && newick[pos] == ')') pos++;` with no
`else`. Consequently malformed, unbalanced inputs were **silently accepted**:

| Malformed input | Old (buggy) behaviour | Correct (per grammar) |
|-----------------|-----------------------|-----------------------|
| `(A,B`   (unclosed, no `;`) | accepted as `(A,B)` | reject (no matching `)`) |
| `(A,B;`  (unclosed, with `;`) | accepted as `(A,B)` | reject |
| `((A,B);` (missing inner `)`) | accepted as degenerate `((A,B))` — **wrong topology** | reject |

The pre-existing trailing-garbage guard (`ParseNewick:655`) only fires on input remaining *after* the
root subtree, so an unclosed `(` (which leaves nothing trailing) slipped through. `(A,B));` (extra `)`)
and `(A,B);extra` (trailing junk) already threw correctly.

**Fix** (`PhylogeneticAnalyzer.cs`, in `ParseNewickRecursive` after the descendant loop): require the
`)`; throw `FormatException("…unbalanced parentheses — an opening '(' has no matching ')'…")` when it
is absent. This makes the parser enforce `Internal → "(" BranchSet ")" Name` and matches the existing
strictness for trailing garbage and surplus `)`.

### Tests added (4, strict)

In `PhylogeneticAnalyzer_NewickIO_Tests.cs` (Malformed-Input region):
- `ParseNewick_ExtraCloseParen_Throws` — `(A,B));`
- `ParseNewick_UnclosedParen_Throws` — `(A,B` and `(A,B;` (asserts "unbalanced" message)
- `ParseNewick_MissingInnerCloseParen_Throws` — `((A,B);`

### Test-quality gate

- Exact round-trip strings (`(A,B,C);`, `(A:0.1000,B:0.2000,C:0.3000);`), exact child counts, and
  ordered leaf sets — all derived from the grammar, not echoed from code output.
- The no-child-dropped test would fail if any of the 3 children were dropped.
- **Mutation check:** reverting the fix (restoring the lenient `)` handling) makes
  `ParseNewick_UnclosedParen_Throws` and `ParseNewick_MissingInnerCloseParen_Throws` **fail 2/2** —
  the new tests genuinely pin the fix; no green-washing.

**Stage B verdict: PASS** (after fix).

## Verdict & follow-ups

- **Stage A:** PASS. **Stage B:** PASS. **End-state:** ✅ CLEAN.
- **Defect:** unbalanced-parentheses silent acceptance — found, root-caused, fixed in code, locked by
  4 strict tests, mutation-verified. Logged in `FINDINGS_REGISTER.md`; ledger Phase-3 row #8.
- **Full suite:** 6776 passed / 0 failed; build 0 errors (4 pre-existing NUnit2007 warnings only in
  unrelated `ApproximateMatcher_EditDistance_Tests.cs`).
- No further follow-ups.
