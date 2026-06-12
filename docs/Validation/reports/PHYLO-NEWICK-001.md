# Validation Report: PHYLO-NEWICK-001 — Newick format tree I/O (parse + serialize)

- **Validated:** 2026-06-12   **Area:** Phylogenetic
- **Canonical method(s):** `PhylogeneticAnalyzer.ToNewick(PhyloNode, bool)`, `PhylogeneticAnalyzer.ParseNewick(string)`
  - Source: `src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs:471` (ToNewick), `:542` (ParseNewick)
  - Tests: `tests/Seqeron/Seqeron.Genomics.Tests/PhylogeneticAnalyzer_NewickIO_Tests.cs`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS (binary-only by design; multifurcation and trailing garbage are now rejected with a clear `FormatException` — silent truncation eliminated; full N-ary support remains a documented scope boundary)

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia, "Newick format"** (https://en.wikipedia.org/wiki/Newick_format) — fetched. Confirms the formal grammar verbatim:
  ```
  Tree     → Subtree ";"
  Subtree  → Leaf | Internal
  Leaf     → Name
  Internal → "(" BranchSet ")" Name
  BranchSet→ Branch | Branch "," BranchSet
  Branch   → Subtree Length
  Name     → empty | string
  Length   → empty | ":" number
  ```
  Confirms the eight canonical example trees, including `(A,B,(C,D));` (leaf names only), `(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);` (popular format), `(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;` (full format). Confirms: whitespace prohibited only within numbers and unquoted strings (allowed elsewhere); quoted labels use single quotes with doubled-quote escaping; underscore→blank in unquoted strings; comments in `[]`.
- **Olsen (1990), "Interpretation of Newick's 8:45 Tree Format"** (https://phylipweb.github.io/phylip/newick_doc.html) — fetched. Confirms:
  ```
  tree            ::= descendant_list [root_label] [':' branch_length] ';'
  descendant_list ::= '(' subtree {',' subtree} ')'
  subtree         ::= descendant_list [internal_node_label] [':' branch_length]
                    | leaf_label [':' branch_length]
  ```
  Confirms the root may carry both a label and a branch length, and that unquoted labels may NOT contain blanks, parentheses, square brackets, single quotes, colons, semicolons, or commas.
- **PHYLIP / Felsenstein** (http://evolution.genetics.washington.edu/phylip/newicktree.html) — old URL now serves only a redirect notice; canonical content is mirrored at phylipweb.github.io and was cross-checked via the Olsen doc and the Evidence file's PHYLIP example table (`A;`, `(B,(A,C,E),D);`, multifurcation, single taxon).

### Grammar rules verified
- Tree terminates with `;` (N1). ✓
- Nested parentheses express internal structure; siblings separated by `,`. ✓
- Node labels optional; optional `:branch_length` after a label/subtree (N5). ✓
- Internal node label allowed after `)` (N8). ✓
- Root may have a label and a trailing `:branch_length` (Olsen). ✓
- Branch-length numbers are locale-independent, `.` decimal separator (N9). ✓

### Edge-case semantics (sourced)
- Empty/whitespace-only/null input is not a tree → reject (Wikipedia: a tree requires at least `Subtree ";"`). ✓
- Single leaf `A;` is a valid tree (PHYLIP example). ✓
- Missing `;` is technically malformed (grammar requires it); lenient stripping is an implementation choice, documented. ✓

### Independent cross-check (hand-constructed)
- `(A,B,(C,D));` — full grammar: a root with **three** children A, B, and the subtree (C,D); 4 leaves.
- `(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);` — same topology with branch lengths 0.1/0.2 on A/B, 0.3/0.4 on C/D, 0.5 on the (C,D) clade.
- `((A,B),(C,D));` — root with two children, each a cherry; 4 leaves; serialization reproduces the structure.

### Findings / divergences
Description is correct and faithfully sourced. The TestSpec/Evidence explicitly and correctly document the implemented subset (binary trees, no quoted labels, no `[]` comments, no underscore→blank, F4 precision) with spec references. **Stage A: PASS.**

## Stage B — Implementation

### Code path reviewed
- `ToNewick` (`PhylogeneticAnalyzer.cs:471`) → `ToNewickRecursive` (`:481`); internal-label gate `IsValidUnquotedNewickLabel` (`:527`).
- `ParseNewick` (`:542`) → `ParseNewickRecursive` (`:565`); `ParseLabel` (`:620`); `ParseNumber` (`:634`).

### Formula/grammar realised correctly? (evidence — hand traces vs code)
- **`(A,B);`** → Trim, strip `;`. pos=0 `(`; Left=leaf `A`; skip `,`; Right=leaf `B`; skip `)`. Root internal, Left=A, Right=B. Matches M05. ✓
- **`(A:0.1,B:0.2);`** → `ParseNumber` reads `0.1`→Left.BL, `0.2`→Right.BL via `double.TryParse(..., NumberStyles.Float, InvariantCulture)`. Matches M06. ✓
- **`((A,B),(C,D));`** → root internal; Left=`(A,B)` internal; Right=`(C,D)` internal; 4 leaves. Matches M07. ✓
- **`(A:0.1,B:0.2):0.0;`** → after the subtree, `ParseNewick:556` sees `:` and assigns root BL 0.0 (Olsen root branch length). Matches M14. ✓
- **`(A:0.1,B:0.2)Root:0.5;`** → internal-name branch (`:603`) reads `Root`, then root `:` reads 0.5. Matches M15. ✓
- **Serialization** uses `F4` + `InvariantCulture` (`:496`/`:506`), `.` decimal separator (N9). Internal labels emitted only if `IsValidUnquotedNewickLabel` (Olsen prohibited-char rule). UPGMA/NJ auto-names like `(Human,Chimp)` contain metacharacters → correctly suppressed (M18). ✓
- **Round-trip** `((A:0.1,B:0.2)AB:0.3,(C:0.4,D:0.5)CD:0.6)Root;` → serialize→parse preserves topology, leaf+internal names, and all six branch lengths within ±0.0001 (F4 gives ±0.00005). Matches M16/S03. ✓

### Edge cases in code
- Empty/null/whitespace → `string.IsNullOrWhiteSpace` guard throws `ArgumentException` (`:544`). M11/M12/S05a. ✓
- `ToNewick(null)` → returns `""` (`:473`). M13. ✓
- Missing `;` → `EndsWith(";")` check is skipped, parses leniently (`:548`). S06. ✓
- Single taxon `A;` → leaf branch of `ParseNewickRecursive`. C01. ✓
- Trailing whitespace → handled by `Trim()` (`:547`). S05b. ✓

### Cross-verification table recomputed vs code
All MUST (M01–M18), SHOULD (S01–S06), COULD (C01) tests recomputed by trace and executed: 29 test methods, all pass.

### Variant/delegate consistency
`ToNewick` and `ParseNewick` are the only public surface; `GetLeaves`/`CalculateTreeLength` are consistent helpers. No `*Fast`/delegate variants.

### Numerical robustness
`ParseNumber` accepts digits, `.`, `+`, `-`, `e`, `E` (scientific notation) and falls back to 0 on parse failure; `InvariantCulture` prevents locale issues. F4 output is bounded; no overflow/div-by-zero.

### Test quality audit
Tests assert exact sourced values (leaf names, branch lengths 0.1/0.2/0.3/0.4/0.5/0.6, internal names, root BL), not just "no-throw" or tautologies. Topology checked via sibling-set equivalence. Deterministic. Real and adequate.

### Findings / defects
No defect against the unit's stated (binary-tree) contract. One **documented scope divergence from the full grammar**, confirmed by trace:

- **Multifurcation silently truncated.** `ParseNewickRecursive` reads exactly one Left child, one comma, one Right child, then expects `)`. For a 3-child node such as the Wikipedia/PHYLIP example `(A,B,(C,D));`, after parsing Left=A and Right=B the cursor sits on `,` (not `)`); the `)`-skip and root-`:` checks both no-op and the remaining `,(C,D))` is **dropped silently** (no throw, no error). The retained binary subset is correct; extra children of any node are lost.
  - This is exactly the **"Binary trees only"** limitation documented in `TestSpecs/PHYLO-NEWICK-001.md` and the Evidence file, with the correct spec citation (`BranchSet → Branch | Branch "," BranchSet` permits N children). UPGMA/NJ — the only producers consumed by this library — emit strictly bifurcating trees, so the round-trip contract holds for all data the system generates. Recorded as PASS-WITH-NOTES rather than a defect.

Not-implemented (also documented, out of scope, not exercised by current algorithms): quoted labels, `[]` comments, underscore→blank. F4 branch-length precision (±0.00005) is a deliberate, documented formatting choice.

## Fix applied (2026-06-12)

The "multifurcation silently truncated" divergence (Stage B note) has been eliminated. The parser now **rejects** input it cannot faithfully represent in the binary `PhyloNode` model, instead of silently dropping children:

- **Multifurcation throws.** `ParseNewickRecursive` (`PhylogeneticAnalyzer.cs`), after reading the Left and Right children and their branch lengths, now checks for a comma at the same level. A comma there means a third (or further) child — a multifurcation — and the parser throws `FormatException("Multifurcating (non-binary) Newick trees are not supported by the binary tree model; ... at position {pos}.")`. This fires for both top-level (`(A,B,C);`) and nested (`(A,B,(C,D,E));`) multifurcations.
- **Trailing garbage throws.** `ParseNewick`, after consuming the root subtree and optional root branch length (the terminal `;` and surrounding whitespace are stripped beforehand), now throws `FormatException("Malformed Newick string: unexpected trailing input ... at position {pos}.")` if any input remains unconsumed (e.g. `(A,B);extra`).
- **Exception type:** `FormatException` (malformed-format convention); the pre-existing empty/null guard keeps its `ArgumentException`.
- **No change to successful parses.** Valid binary trees parse exactly as before: `(A,B);`, `(A:0.1,B:0.2);`, `((A,B),(C,D));`, single leaf `A;`, internal labels, root branch length/label, and serialize→parse round-trips (topology + branch lengths within ±0.0001) all unchanged.
- **Safe for internal pipelines.** Confirmed by grep: `ParseNewick` is only called on external input (MCP tools) and in tests. UPGMA/NJ build strictly bifurcating `PhyloNode`s directly (Left/Right) and `ToNewickRecursive` emits only binary output, so no internal round-trip can produce a multifurcating string — the new throw cannot break them.
- **Scope boundary unchanged.** Full N-ary Newick support (a child-list `PhyloNode` model and refactored consumers: UPGMA/NJ/RF/MRCA/serialize) remains explicitly out of scope.

Tests added to `PhylogeneticAnalyzer_NewickIO_Tests.cs`: multifurcating-root throws, multifurcating-nested-subtree throws, trailing-garbage throws, plus regression tests asserting valid binary trees and round-trip still parse. Newick filter: 34 passed. Full suite: 4470 passed, 0 failed.

## Verdict & follow-ups
- **Stage A: PASS** — grammar and examples verified against Wikipedia + Olsen (1990); description faithfully sourced.
- **Stage B: PASS** — implementation correctly realises the binary-tree Newick subset (parse, serialize, round-trip, root branch length, internal-label suppression); multifurcating input and trailing garbage are now **rejected with a clear `FormatException`** rather than silently truncated. Silent data loss eliminated.
- **State: CLEAN** — fix applied conservatively (no tree-model refactor); full suite green (4470 passed).
- **Follow-up (optional, out of scope):** if multifurcating Newick input ever becomes required, extend `ParseNewickRecursive` to a child-list model (and `PhyloNode` to N children) so non-binary trees parse without loss.
