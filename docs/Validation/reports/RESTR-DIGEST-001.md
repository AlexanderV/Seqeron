# Validation Report: RESTR-DIGEST-001 — Restriction Digest Simulation (circular + linear fragment math)

- **Validated:** 2026-06-24   **Area:** MolTools
- **Re-validation scope:** the circular-molecule digest overload added in commit `293e139`
  (`Digest(DnaSequence, MoleculeTopology, params string[])` + `MoleculeTopology { Linear, Circular }`)
  and the unchanged linear path. Fragment-count + exact fragment-length math verified by hand for
  BOTH topologies and the 0-site / 1-site edge cases.
- **Canonical method(s):** `RestrictionAnalyzer.Digest(DnaSequence, MoleculeTopology, params string[])`
  (`RestrictionAnalyzer.cs:350`), private `DigestCircular` (`:365`); control:
  `Digest(DnaSequence, params string[])` (`:259`).
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ✅ PASS
- **End-state:** ✅ CLEAN (no code defect, no new code change this session; tests confirmed green)

## Stage A — Description

### Sources opened & confirmed THIS session
- **WebSearch — standard fragment-count rule (Quora summary + Vanderbilt restriction-mapping guide,
  Addgene Plasmids 101):** *"For circular DNA you get as many fragments as sites present in the
  molecule"* and *"If the DNA is linear, then the number of cuts is one less than the number of
  fragments"* ⇒ **circular: k sites → k fragments; linear: k sites → k+1 fragments.** The geometric
  argument (N cuts in a ring → N arcs; N cuts in a segment → N+1 pieces) is the same rule.
- **Prior session sources (re-confirmed, unchanged):** Univ. of Illinois MolBio "If you cut a circle
  once, you get one linear fragment... cut it a second time to get 2 linear fragments" (1 cut
  linearizes, does not separate); EMBOSS `restrict -plasmid` (origin-spanning sites/fragments must be
  considered for circular; default is linear); Wikipedia EcoRI cleavage **`G↓AATTC`** ⇒ a site
  starting at index `i` cuts at `i+1` (DB `CutPositionForward = 1`).

### Fragment-count + length rule (exact)
- **Linear:** k distinct forward-strand cut sites → **k+1** fragments. Cut boundary list is
  `[0, …cuts, len]`; fragment i = `[cuts[i], cuts[i+1])`. Lengths sum to `len`.
- **Circular:** k distinct forward-strand cut sites → **k** fragments. Non-origin fragment i =
  `[cuts[i], cuts[i+1])`; the single origin-spanning fragment = `[lastCut, len) + [0, firstCut)`,
  length **`(len − lastCut) + firstCut`**. Lengths sum to `len`.
- **Edge:** 0 sites circular → 1 uncut full-length circular fragment; 1 site circular → 1
  full-length linearized fragment (the k=1 case of the origin-span loop). 0 sites linear → 1
  whole-sequence fragment.

### Independent hand cross-check (EcoRI cut = site_start + 1; lengths sum to molecule length)

| Example | Seq (len) | Sites→cuts | Topology | Fragments (start, seq, len) | Count | Sum |
|---|---|---|---|---|---|---|
| 2-site | `GAATTCAAAGAATTCAAA` (18) | {0,9}→{1,10} | circular | (1,"AATTCAAAG",9); origin(10,"AATTCAAA"+"G"="AATTCAAAG",9) | 2 | 18 |
| 2-site | same | same | linear | (0,"G",1);(1,…,9);(10,"AATTCAAA",8) | 3 | 18 |
| 1-site | `AAAGAATTCAAA` (12) | {3}→{4} | circular | origin(4,"AATTCAAA"+"AAAG"="AATTCAAAAAAG",12) | 1 | 12 |
| 1-site | same | same | linear | (0,"AAAG",4);(4,"AATTCAAA",8) | 2 | 12 |
| 3-site | `GAATTCAAAGAATTCAAAGAATTCAAA` (27) | {0,9,18}→{1,10,19} | circular | (1,…,9);(10,…,9); origin(19,"AATTCAAA"+"G",9) | 3 | 27 |
| 2-site non-repeat | `CGAATTCTTTGAATTCAA` (18) | {1,10}→{2,11} | circular | (2,"AATTCTTTG",9); origin(11,"AATTCAA"+"CG"="AATTCAACG",9) | 2 | 18 |

All hand derivations: counts follow linear k+1 / circular k; origin-span length
`(len−lastCut)+firstCut`; every row's fragment lengths sum exactly to the molecule length.

### Source-quality note
A lay aggregator answer claiming "3 sites → 4 fragments" for a circle is incorrect (conflates the
linear formula). Authoritative sources (Illinois MolBio, Vanderbilt, Addgene, EMBOSS) and the
geometric argument all give circular **k→k**. The implementer's rule is correct.

### Findings / divergences
None.

## Stage B — Implementation

### Code path reviewed
- **Linear** (`RestrictionAnalyzer.cs:259-318`): collects distinct **forward-strand** cut positions
  (avoids palindromic double-count), builds boundary list `{0} ∪ cuts ∪ {len}`, emits
  `[cuts[i], cuts[i+1])` → exactly k+1 fragments; 0 cuts → one whole fragment (`:281`).
- **Topology overload** (`:350`): identical arg-guards; `Linear` delegates byte-for-byte to the
  original method; `Circular` → `DigestCircular`.
- **`DigestCircular`** (`:365`): distinct forward-strand cuts normalised `cut % length` (origin/end
  identity per EMBOSS `-plasmid`); 0 cuts → one uncut full-length fragment (`:387`); else k fragments
  — non-origin `[start, nextCut)` (`:414-417`), origin-spanning `tail(seq[start..n]) + head(seq[0..nextCut])`
  with `length = (n − start) + nextCut` (`:419-425`). On a circle both flanks carry an enzyme.

### Formula realised correctly?
Yes. `length = (n − start) + nextCut` is exactly `(len − lastCut) + firstCut`; join order
`tail + head` (last-cut→end, then start→first-cut) is correct. k=1 is the loop's single iteration with
`nextCut = cuts[0]` ⇒ one full-length linearized fragment.

### Cross-verification vs code
All six hand rows above match the actual code output (confirmed via the test suite, which asserts
exact start/length/sequence/sum for the 0/1/2/3-site circular cases and the linear control counts).

### Variant/delegate consistency
`Linear` topology delegates to the original `Digest`; default overload parity holds. Linear path
untouched — all pre-existing linear/summary/map/compatibility tests remain green.

### Test quality audit
- 5 circular tests are discriminating: same seq+enzyme yields different counts under linear vs
  circular (3≠2, 2≠1), so a wrong topology branch fails them; exact lengths/sequences/start asserted;
  0/1/2/3 sites covered; sum invariant checked; arg-guards covered.
- The non-tandem-repeat hardening test (`…OriginSpanningSequenceIsUnique_NotTandemRepeat`,
  `CGAATTCTTTGAATTCAA`) asserts the origin-span "AATTCAACG" appears nowhere as a contiguous substring,
  so only a correct `tail+head` join can produce it; guards against the prior tandem-repeat
  self-confirmation. Verified present and passing.

### Findings / defects
None.

## Verdict & follow-ups
- **Stage A:** ✅ PASS — linear k+1 and circular k confirmed against Quora/Vanderbilt/Addgene/Illinois
  MolBio/EMBOSS; origin-span `(len−lastCut)+firstCut`, 0/1/k semantics sourced and hand-verified.
- **Stage B:** ✅ PASS — code realises both rules exactly; filtered run **149 passed, 0 failed**
  (`FullyQualifiedName~RestrictionAnalyzer`), build 0 warnings/0 errors.
- **End-state:** ✅ CLEAN. No code change this session; no open follow-ups.
