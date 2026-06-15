# Validation Report: ANNOT-REPEAT-001 — Repetitive Element Detection and Classification

- **Validated:** 2026-06-15   **Area:** Annotation
- **Canonical method(s):** `GenomeAnnotator.FindRepetitiveElements(string, int minRepeatLength, int minCopies)`, `GenomeAnnotator.ClassifyRepeat(string, IReadOnlyDictionary<string,string>)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS (defect found and fully fixed this session)

## Stage A — Description

### Sources opened this session (retrieved, not just cited)
1. **Wikipedia — Tandem repeat** (https://en.wikipedia.org/wiki/Tandem_repeat). Confirmed verbatim definition "A pattern of one or more nucleotides is repeated and the repetitions are directly adjacent to each other"; worked example `ATTCG ATTCG ATTCG` (the spec's M1 dataset); microsatellite/STR < 10 nt; minisatellite "between 10 and 60 nucleotides are repeated".
2. **NIH / genome.gov Genetics Glossary — Tandem Repeat** (https://www.genome.gov/genetics-glossary/Tandem-Repeat). Confirmed the **"two or more"** / head-to-tail minimum: "A sequence of two or more DNA bases that is repeated numerous times in a head-to-tail manner on a chromosome." This is the authoritative basis for M2 and for `minCopies >= 2`.
3. **Wikipedia — Inverted repeat** (https://en.wikipedia.org/wiki/Inverted_repeat). Confirmed "a single stranded sequence of nucleotides followed downstream by its reverse complement"; spacer "any length including zero"; gap-0 ⇒ palindrome; worked example `5'-TTACGnnnnnnCGTAA-3'` with revcomp(`TTACG`)=`CGTAA` (the spec's M4 dataset).
4. **Hampson et al. (2021) IUPACpal, BMC Bioinformatics 22:51** (https://pmc.ncbi.nlm.nih.gov/articles/PMC7866733/). Confirmed formal IR definition WW̄ᴿ; gapped IR WGW̄ᴿ with |G| ≥ 0; perfect (ungapped, 0 mismatch) vs imperfect (δ_H ≤ k); detection params = min/max arm length, max gap, max mismatches.
5. **RepeatMasker documentation** (referenced; class vocabulary): SINE, LINE, LTR, DNA, Satellite, Simple_repeat, Low_complexity, RNA, Unclassified/Unknown; classification by best homology match to a library entry.

### Formula / definition check
- **Tandem repeat:** head-to-tail adjacent copies of a primitive unit, ≥ 2 copies. Matches sources 1 + 2 exactly. Span is `copies × unitLen`; 0-based, end-exclusive (INV-1/INV-2). ✔
- **Inverted repeat:** WGW̄ᴿ, right arm = reverse complement of left arm, gap |G| ≥ 0, gap-0 = palindrome (INV-3). Matches sources 3 + 4 exactly. ✔
- **Hand-checked reverse complements:** revcomp(`GAA`)=`TTC` (so `GAATTC`, EcoRI, is a gap-0 IR/palindrome — M3); revcomp(`TTACG`)=`CGTAA` (so `TTACGAAAAAACGTAA` is a gap-6 IR — M4); revcomp(`GGATCC`)=`GGATCC` (BamHI palindrome). All confirmed. ✔
- **Classification vocabulary** (INV-4): all labels are from the RepeatMasker source class set. ✔

### Edge-case semantics
- Single copy is not a tandem repeat (source 2 "two or more"). ✔
- Non-primitive units collapse to primitive period (source 1 corner case). ✔
- Gap-0 IR = palindrome (sources 3, 4). ✔
- No library match → Unclassified/Unknown, with a documented simple-repeat (STR 1–6 bp) fallback by motif size (source 1 + RepeatMasker). ✔

### Independent cross-check (numbers)
- M1 `ATTCGATTCGATTCG`: unit `ATTCG` (5 bp), 3 copies, span [0,15) — matches Wikipedia worked example.
- M4 `TTACGAAAAAACGTAA`: left `TTACG`[0,5), gap `AAAAAA`[5,11), right `CGTAA`[11,16); revcomp(`TTACG`)=`CGTAA` — matches Wikipedia IR example.

**Stage A findings:** none. The description is biologically and mathematically correct. PASS.
(The one documented simplification — exact substring containment instead of Smith-Waterman/Repbase homology — is a Framework limitation explicitly flagged in Assumption 1, not an error. It was tightened this session, see Stage B.)

## Stage B — Implementation

### Code path reviewed
- `GenomeAnnotator.FindRepetitiveElements` / `FindRepetitiveElementsCore` / `FindTandemRepeats` / `IsPrimitive` / `FindInvertedRepeats` — `src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs:652-792`.
- `GenomeAnnotator.ClassifyRepeat` / `IsSimpleRepeat` — same file, `:808-869`.

### Formula realised correctly?
- **Tandem repeats:** scans unit lengths 1..min(len/minCopies, 60); skips non-primitive units (`IsPrimitive`); counts adjacent equal blocks; enforces `minCopies` and `minRepeatLength`; reports left-maximal arrays de-duplicated by span. Verified by direct execution: `ATTCGATTCGATTCG`→[0,15) `ATTCG`×3; `AAAAAA`→[0,6) primitive `A`. ✔ (It also reports legitimate overlapping arrays such as `TTCGA`×2 at [1,11) — these are genuine tandem repeats per the definition, not false positives; tests filter by exact span.)
- **Inverted repeats:** for each left arm computes `GetReverseComplementString` and searches the right arm within gap 0..50. Verified: `GAATTC`→[0,6) palindrome; `TTACGAAAAAACGTAA`→[0,16) `TTACG[AAAAAA]CGTAA`. INV-3 holds for all three hand-checked enzymes. ✔
- **Input contract:** null→`ArgumentNullException`, `minCopies<2`→`ArgumentOutOfRangeException`, `minRepeatLength<1`→throws, empty→empty. Verified. ✔

### DEFECT FOUND AND FIXED — `ClassifyRepeat` bidirectional containment
**Root cause.** Matching used `query.Contains(element) || element.Contains(query)`. The reverse direction (`element.Contains(query)`) let a trivially short query inherit a class purely because a longer consensus contains those letters. Direct execution before the fix:
```
ClassifyRepeat("A",  {Alu:20bp, Satellite:6bp})  -> "SINE/Alu"   (WRONG: "A" is not an Alu)
ClassifyRepeat("GGCCGGG", {Alu:20bp})            -> "SINE/Alu"   (WRONG: fragment of consensus)
```
A 1-bp query matches almost any non-trivial library. This contradicts RepeatMasker, which **screens the query for occurrences of known library elements** (element found *within* the query) and reports the best match — it does not classify a query that is merely a substring of a consensus.

**Why uncaught.** No MUST test depended on the reverse direction; no existing test exercised a short or sub-consensus query. The behaviour would also pass a deliberately-wrong implementation (code-echo blind spot).

**Fix (this session).** Matching made one-directional (`element ⊆ query`); XML doc updated to the screen-query-for-elements model. After the fix:
```
ClassifyRepeat("A", ...)        -> "Unknown"     ClassifyRepeat("GGCCGGG", ...) -> "Unknown"
ClassifyRepeat("CCCCGGCCGGGCGCGGTGGCTCACGGGG", ...) -> "SINE/Alu"   (M5 still correct)
```
Evidence Assumption 1 and TestSpec Assumption 1 corrected to drop "either direction / is contained by".

### Cross-verification table recomputed vs code (post-fix)
| Case | Input | Expected (source) | Actual |
|------|-------|-------------------|--------|
| M1 tandem | `ATTCGATTCGATTCG` | [0,15) `ATTCG`×3 | [0,15) ✔ |
| M3 IR gap0 | `GAATTC` | [0,6) `GAATTC`, revcomp arms | [0,6) ✔ |
| M4 IR gap6 | `TTACGAAAAAACGTAA` | [0,16) `TTACG[AAAAAA]CGTAA` | [0,16) ✔ |
| S1 primitive | `AAAAAA` | [0,6) primitive `A` | [0,6) ✔ |
| M5 classify | query ⊇ 20bp Alu | `SINE/Alu` | `SINE/Alu` ✔ |
| M6 classify | `CACACACA` | `Simple_repeat` | `Simple_repeat` ✔ |
| M7 classify | random 21bp | `Unknown` | `Unknown` ✔ |
| M7b (new) | `"A"` | `Unknown` | `Unknown` ✔ |
| M7c (new) | `GGCCGGG` (Alu fragment) | `Unknown` | `Unknown` ✔ |
| M5b (new) | query ⊇ Alu + Satellite | longest=`SINE/Alu` | `SINE/Alu` ✔ |

### Variant/delegate consistency
Two independent public methods; no `*Fast`/delegate variants. Case-insensitivity verified (`gaattc`→`GAATTC` IR).

### Test-quality audit (HARD gate)
- **Before:** 13 tests, all happy/edge but **no test exercised the `ClassifyRepeat` reverse-containment direction or a short query** → the defect was invisible; M7-style tests would pass a wrong impl.
- **After:** added M7b, M7c, M5b, M6b (ClassifyRepeat) and an INV-3 reverse-complement arm test (FindRepetitiveElements). All expectations trace to external sources (RepeatMasker screen-for-elements model; hand-computed reverse complements). No weakened assertions, no widened tolerances, no skips. Exact-value asserts throughout.
- **Honest green:** full **unfiltered** suite `Failed: 0, Passed: 6566` (1 pre-existing benchmark skipped, unrelated). `dotnet build` 0 errors; the 4 build warnings are pre-existing in an unrelated file (`ApproximateMatcher_EditDistance_Tests.cs`) — the changed files (`GenomeAnnotator.cs`, this test file) build warning-free.

**Test-quality gate result: PASS.**

## Verdict & follow-ups
- **Stage A:** PASS. **Stage B:** PASS (one real defect — `ClassifyRepeat` bidirectional containment misclassifying short queries — found and **completely fixed** this session: code + tests + Evidence + TestSpec).
- **End-state: CLEAN.** Algorithm fully functional; full suite green.
- Logged as finding **A17** in `FINDINGS_REGISTER.md`; ledger row **58**.
