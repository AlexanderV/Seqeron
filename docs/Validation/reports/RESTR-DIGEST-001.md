# Validation Report: RESTR-DIGEST-001 — Restriction Digest Simulation (circular-molecule enhancement)

- **Validated:** 2026-06-17   **Area:** MolTools
- **Re-validation scope:** the NEW circular-molecule digest added in commit `293e139`
  (`Digest(DnaSequence, MoleculeTopology, params string[])` overload + `MoleculeTopology { Linear, Circular }`),
  plus confirmation the unchanged linear path is unaffected.
- **Canonical method(s):** `RestrictionAnalyzer.Digest(DnaSequence, MoleculeTopology, params string[])`,
  private `DigestCircular`; control: `Digest(DnaSequence, params string[])`.
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ✅ PASS
- **End-state:** ✅ CLEAN (no code defect; one test-quality hardening added)

## Stage A — Description

### Sources opened & retrieved THIS session
- **Univ. of Illinois MolBio, "Restriction Digestion/Gel Electrophoresis Assignment 1"** —
  https://www.life.illinois.edu/molbio/geldigest/assign1.html — verbatim: *"If you cut a circle
  once, you get one linear fragment. You must cut it a second time to get 2 linear fragments."*
  ⇒ circular molecule: **k cut sites → k fragments** (1 cut linearizes, does not separate).
- **EMBOSS `restrict` documentation (`-plasmid` option)** —
  https://emboss.sourceforge.net/apps/cvs/emboss/apps/restrict.html — verbatim: *"If this is set
  then this allows searches for restriction enzyme recognition site and cut positions that span
  the end of the sequence to be considered."* ⇒ circular digests must consider origin-spanning
  (wrap-around) sites/fragments; default mode is linear.
- **Addgene, "Plasmids 101: How to Verify Your Plasmid Using a Restriction Digest"** (search-confirmed)
  — a single cutter linearizes the plasmid into one fragment; number of fragments = number of cut sites.
- **Wikipedia, "EcoRI"** — https://en.wikipedia.org/wiki/EcoRI — verbatim cleavage **`G↓AATTC`**
  (cut after the first G, 1 nt into the recognition site). Validates the DB entry
  `("EcoRI","GAATTC",1,5,…)`, i.e. `CutPositionForward = 1` ⇒ a site starting at index `i` cuts at `i+1`.
- **WebSearch** (fragment-count rule): the standard rule is linear **k+1** vs circular **k**.

### Fragment-count rule (the exact rule recorded)
- **Linear:** k distinct forward-strand cut sites → **k+1** fragments (unchanged).
- **Circular:** k distinct forward-strand cut sites → **k** fragments. The origin-spanning fragment
  joins the last-cut→end piece with the start→first-cut piece, of length
  **`(SequenceLength − lastCut) + firstCut`**. 0 sites → one full-length uncut circular fragment;
  1 site → one full-length linearized fragment.

### Source-quality note
A search aggregator (toppr/brainly) snippet claimed "3 sites → 4 fragments" for a circle; this is a
known-incorrect lay answer that conflates the linear formula. The authoritative Illinois MolBio
source and the geometric argument both give circular **k→k**. The implementer's rule is correct.

### Independent cross-check (hand-derived, NOT lifted from repo tests)
EcoRI cut = (site start) + 1.

| Example | Seq (len) | Sites→cuts | Circular fragments (start, seq, len) | Sum | Linear count (control) |
|---|---|---|---|---|---|
| 2-site | `GAATTCAAAGAATTCAAA` (18) | {0,9}→{1,10} | (1,"AATTCAAAG",9); origin (10,"AATTCAAA"+"G"="AATTCAAAG",9) | 18 | 3 |
| 1-site | `AAAGAATTCAAA` (12) | {3}→{4} | origin (4,"AATTCAAA"+"AAAG"="AATTCAAAAAAG",12) | 12 | 2 |
| 3-site | `GAATTCAAAGAATTCAAAGAATTCAAA` (27) | {0,9,18}→{1,10,19} | (1,…,9);(10,…,9); origin (19,"AATTCAAA"+"G",9) | 27 | 4 |
| 2-site (validator, non-repeat) | `CGAATTCTTTGAATTCAA` (18) | {1,10}→{2,11} | (2,"AATTCTTTG",9); origin (11,"AATTCAA"+"CG"="AATTCAACG",9) | 18 | 3 |

All four hand-derivations match the code output. Origin-span length `(len−lastCut)+firstCut` confirmed.

### Findings / divergences
None. Rule and edge semantics (0/1/k) are sourced and correct.

## Stage B — Implementation

### Code path reviewed
`RestrictionAnalyzer.cs:350-437`. The topology overload (`:350`) validates args identically to the
linear overload and delegates `Linear` to `Digest(DnaSequence, params string[])` (byte-for-byte
unchanged default behaviour). `DigestCircular` (`:365`) collects distinct **forward-strand** cut
positions (same convention as linear, avoids palindromic double-count), normalises each `cut %
length` (origin/end identity, per EMBOSS `-plasmid`), returns: 0 sites → one uncut full-length
fragment (`:387`); else **k** fragments, each cut→next-cut cyclically, with the last (origin-spanning)
fragment built as `tail + head` of length `(n−start)+nextCut` (`:419-426`).

### Formula realised correctly?
Yes. `length = (n − start) + nextCut` is exactly `(len − lastCut) + firstCut`; the join order
`tail + head` (last-cut→end then start→first-cut) is correct. 1 site is the k=1 case of the loop:
`nextCut = cuts[0] = start` ⇒ a single origin-spanning fragment of full length.

### Cross-verification vs code
The four-row table above was recomputed against actual code output (via the test run) — all counts,
start positions, lengths, and sequences match. `Digest(seq, "EcoRI")` (default) equals
`Digest(seq, MoleculeTopology.Linear, "EcoRI")` (asserted in the suite).

### Variant/delegate consistency
`Linear` topology delegates to the original method; default overload parity asserted. Linear path
untouched (all pre-existing linear/summary/map/compatibility tests still green).

### Test quality audit
- Implementer's 5 circular tests: discriminating linear-vs-circular counts on the SAME seq+enzyme
  (3≠2, 2≠1) — a test passing against the linear impl on circular input would fail these; exact
  lengths/sequences/start asserted; 0/1/2/3 sites covered; sum invariant checked; arg-guards. Strong.
- **Gap identified:** the implementer's circular examples are tandem repeats where every fragment
  (including the origin-spanning one) reads "AATTCAAAG", so the wrap-around *sequence* assertion is
  partly self-confirming — a subtly wrong slice could coincidentally still read "AATTCAAAG".
- **Hardening added (test-only, 0 code change):**
  `Digest_Circular_TwoSites_OriginSpanningSequenceIsUnique_NotTandemRepeat` on the non-repetitive
  18 bp circle `CGAATTCTTTGAATTCAA` (cuts {2,11}). Origin-span = "AATTCAACG", which appears nowhere
  as a contiguous substring of the molecule (explicitly asserted), so it can only arise from a
  correct `tail+head` join; also asserts the two fragments are not equal, linear=3 vs circular=2,
  and sum=18.
- **Mutation check:** swapping the join to `head + tail` fails the new test AND the 3 implementer
  discriminating tests; reverted. Test-quality gate: PASS.

### Findings / defects
No code defect. One Stage-B test-quality gap (tandem-repeat self-confirmation) closed with a strict
sourced test.

## Verdict & follow-ups
- **Stage A:** ✅ PASS — circular k→k, origin-span `(len−lastCut)+firstCut`, 0/1/k semantics all
  confirmed against Illinois MolBio, EMBOSS `-plasmid`, Addgene, Wikipedia EcoRI.
- **Stage B:** ✅ PASS — code realises the rule exactly; tests hardened; full unfiltered suite
  **6769 passed, 0 failed**, build 0 errors.
- **End-state:** ✅ CLEAN. No open follow-ups.
