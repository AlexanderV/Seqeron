# Validation Report: PHYLO-NEWICK-001 — Newick format I/O (N-ary / multifurcating)

- **Validated:** 2026-06-24   **Area:** Phylogenetics
- **Canonical method(s):** `PhylogeneticAnalyzer.ParseNewick(string)`, `PhylogeneticAnalyzer.ToNewick(PhyloNode, bool)`, `PhyloNode.Children` (with `Left`/`Right` convenience accessors over the first two) — `src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs`
- **Under test (re-confirmation):** commit c4f0190 (N-ary multifurcating tree model) + commit 650f18d1 (mandatory closing `)`, reject silent unbalanced-paren acceptance).
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN (no defect found this session; prior fix re-confirmed independently)

## Stage A — Description (the grammar)

### Sources opened this session (fetched live, not memory)

1. **Wikipedia, "Newick format"** — <https://en.wikipedia.org/wiki/Newick_format>
   BNF productions (fetched verbatim):
   - `Tree → Subtree ";"`
   - `Subtree → Leaf | Internal`
   - `Leaf → Name`
   - `Internal → "(" BranchSet ")" Name`
   - `BranchSet → Branch | Branch "," BranchSet`
   - `Branch → Subtree Length`
   - `Name → empty | string`
   - `Length → empty | ":" number`

   The recursive `BranchSet → Branch | Branch "," BranchSet` permits an arbitrary-length,
   comma-separated list of branches under one internal node. A node with **≥3** children is a
   **multifurcation / polytomy**. `(A,B,C);` is a single internal node (the root) with three leaf
   children A, B, C. The production `Internal → "(" BranchSet ")" Name` makes the closing `)`
   **mandatory** — an opened descendant list that is never closed is not a valid `Tree`.

2. **Felsenstein / Olsen, PHYLIP "The Newick tree format"** — <https://phylipweb.github.io/phylip/newicktree.html>
   Fetched verbatim: *"Trees can be multifurcating at any level."*; *"Interior nodes are
   represented by a pair of matched parentheses. Between them are representations of the nodes
   that are immediately descended from that node, separated by commas."*; *"Branch lengths can be
   incorporated into a tree by putting a real number, with or without decimal point, after a node
   and preceded by a colon."*; *"The tree ends with a semicolon."*

### Grammar checks

- **Multifurcation is a first-class grammar feature** — confirmed via the recursive `BranchSet`
  (Wikipedia) and "multifurcating at any level" (Felsenstein/Olsen). Not an extension; it is the
  canonical format.
- **Closing parenthesis mandatory** — `Internal → "(" BranchSet ")" Name`; "matched parentheses"
  (Olsen). Unbalanced parentheses are not a valid tree.
- **Branch length** — `Length → ":" number`, attached per branch; decimal point `.`
  (locale-independent).
- **Terminator** — `Tree → Subtree ";"`.

### Hand-derived expected results (from the grammar, not the code)

| Input | Parse (derived) | Serialize (no lengths) |
|-------|-----------------|------------------------|
| `(A,B,C);` | root, 3 leaf children [A,B,C] in order | `(A,B,C);` |
| `(A:0.1,B:0.2,C:0.3);` | 3 children, BL 0.1/0.2/0.3 | `(A,B,C);` / with-lengths `(A:0.1000,B:0.2000,C:0.3000);` |
| `(A,B,(C,D,E));` | root 3 children: A, B, internal{C,D,E}; 5 leaves | `(A,B,(C,D,E));` |
| `(A,B);` (binary) | root, 2 children A,B | `(A,B);` |
| `(B:6.0,(A:5.0,C:3.0,E:4.0):5.0,D:11.0);` | root 3 children: B, internal{A,C,E}, D; 5 leaves | `(B,(A,C,E),D);` |
| `A;` | single leaf node | `A;` |
| `(A,B`, `(A,B;`, `((A,B);` | **malformed** — unbalanced parens, reject | — |
| `(A,B));`, `(A,B);extra` | **malformed** — surplus `)` / trailing garbage, reject | — |
| `""`, `"   "` | empty / no tree, reject | — |

**Stage A verdict: PASS** — grammar productions, the multifurcation definition, branch-length
syntax, the mandatory closing parenthesis, and the semicolon terminator are confirmed against two
authoritative sources.

## Stage B — Implementation

### Code path reviewed

- `ParseNewick` (`PhylogeneticAnalyzer.cs:657`): null/empty/whitespace → `ArgumentException`
  (`:659`); trims, strips trailing `;`; recurses; handles optional root branch length (`:671`);
  rejects unconsumed trailing input as `FormatException` (`:680`).
- `ParseNewickRecursive` (`:687`): on `(`, loops the `BranchSet` parsing each child + optional
  `:length`, separated by `,`, appending every child to `node.Children` (`:718`) — genuine N-ary.
  After the loop it **requires** the matching `)` (`:734-743`), throwing
  `FormatException("…unbalanced parentheses — an opening '(' has no matching ')'…")` when absent.
  A `MaxParseDepth=1000` guard rejects pathological nesting with a catchable `FormatException`
  (`:693`) instead of a `StackOverflowException`.
- `ToNewick`/`ToNewickRecursive` (`:577`/`:587`): null node → `""` (`:579`); emits every child
  comma-separated regardless of count (genuine N-ary write, `:599-609`), `F4` branch lengths,
  trailing `;` (`:583`), internal labels only when a valid unquoted Newick label (`:617`,
  `IsValidUnquotedNewickLabel :628`).

### Independent verification (executed probe + tests)

A standalone probe (against the actual compiled assembly) parsed each grammar case and printed
child count / leaf list / re-serialisation. All results matched the hand-derived table exactly:

```
(A,B,C);                                -> children=3 leaves=[A,B,C]     write=(A,B,C);     writeL=(A:0.0000,B:0.0000,C:0.0000);
(A:0.1,B:0.2,C:0.3);                    -> children=3 leaves=[A,B,C]     write=(A,B,C);     writeL=(A:0.1000,B:0.2000,C:0.3000);
(A,B,(C,D,E));                          -> children=3 leaves=[A,B,C,D,E] write=(A,B,(C,D,E));
(A,B);                                  -> children=2 leaves=[A,B]       write=(A,B);
(B:6.0,(A:5.0,C:3.0,E:4.0):5.0,D:11.0); -> children=3 leaves=[B,A,C,E,D] write=(B,(A,C,E),D);
A;                                      -> children=0 leaves=[A]         write=A;

(A,B    -> REJECT FormatException "unbalanced parentheses"
(A,B;   -> REJECT FormatException "unbalanced parentheses"
((A,B); -> REJECT FormatException "unbalanced parentheses"
(A,B)); -> REJECT FormatException "unexpected trailing input"
(A,B);extra -> REJECT FormatException "unexpected trailing input"
""      -> REJECT ArgumentException "Newick string is empty."
"   "   -> REJECT ArgumentException "Newick string is empty."
```

- **Multifurcation parsed faithfully, no child dropping.** `(A,B,C);` yields exactly 3 children in
  order; `(A,B,(C,D,E));` yields 5 leaves with the nested polytomy intact; the PHYLIP worked
  example `(B:6.0,(A:5.0,C:3.0,E:4.0):5.0,D:11.0);` yields 3 root children / 5 leaves.
- **Byte-for-byte round-trips:** `(A,B,C);` → `(A,B,C);`; `(A:0.1,B:0.2,C:0.3);` → with-lengths
  `(A:0.1000,B:0.2000,C:0.3000);` (F4); `(A,B,(C,D,E));` → `(A,B,(C,D,E));`.
- **Malformed input rejected:** unbalanced parens (`(A,B`, `(A,B;`, `((A,B);`) → `FormatException`
  (unbalanced); surplus `)` and trailing garbage → `FormatException` (trailing); empty/whitespace
  → `ArgumentException`. The closing-`)` enforcement (`:734-743`) is present and active.

### Test quality audit

`PhylogeneticAnalyzer_NewickIO_Tests.cs` — 36 tests, all green. The multifurcation region
(`:638`) asserts exact `Children.Count`, ordered child names, exact per-child branch lengths, and
byte-for-byte round-trip strings (`(A,B,C);`, `(A:0.1000,B:0.2000,C:0.3000);`, `(A,B,(C,D,E));`) —
derived from the grammar, not echoed from code. The malformed-input region (`:725`) pins each
rejection: `ParseNewick_UnclosedParen_Throws` (`(A,B`, `(A,B;`, asserts "unbalanced" message),
`ParseNewick_MissingInnerCloseParen_Throws` (`((A,B);`), `ParseNewick_ExtraCloseParen_Throws`
(`(A,B));`), `ParseNewick_TrailingGarbage_Throws` (`(A,B);extra`). Empty/null/whitespace assert
`ArgumentException`. The assertions are exact-value, not "no-throw" tautologies.

### Findings / defects

None. The prior session's fix (mandatory closing `)`) is present in `PhylogeneticAnalyzer.cs` and
independently re-confirmed: the unbalanced-paren cases that the old lenient parser silently
accepted now throw `FormatException`. No code was changed this session.

**Stage B verdict: PASS.**

## Verdict & follow-ups

- **Stage A:** PASS. **Stage B:** PASS. **End-state:** ✅ CLEAN.
- **Code changed:** none (re-confirmation only).
- **Full suite:** 18208 passed / 0 failed (1 benchmark skipped); build 0 warnings / 0 errors.
- No follow-ups.
