# Validation Report: ASSEMBLY-CONSENSUS-001 — Consensus Computation

- **Validated:** 2026-06-15   **Area:** Assembly
- **Canonical method(s):** `SequenceAssembler.ComputeConsensus(IReadOnlyList<string> alignedReads, double threshold = 0.5, char ambiguous = 'N')`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm

1. **Biopython `Bio.Align.AlignInfo.SummaryInfo.dumb_consensus` (tag `biopython-179`)** — fetched
   the raw source this session
   (https://raw.githubusercontent.com/biopython/biopython/biopython-179/Bio/Align/AlignInfo.py).
   Confirmed verbatim:
   - Signature: `def dumb_consensus(self, threshold=0.7, ambiguous="X", require_multiple=False):`
   - `con_len = self.alignment.get_alignment_length()` (consensus spans full alignment length).
   - Per-column tally skips gaps: counted only when `record.seq[n] != "-" and record.seq[n] != "."`;
     `num_atoms` increments per non-gap residue.
   - Bounds check `if n < len(record.seq):` guards records shorter than `con_len`.
   - Per-column reset of `atom_dict={}`, `num_atoms=0`, `max_atoms=[]`, `max_size=0`.
   - Max set: `if atom_dict[atom] > max_size: max_atoms=[atom]; max_size=...` /
     `elif atom_dict[atom] == max_size: max_atoms.append(atom)`.
   - Decision rule (verbatim):
     ```python
     if require_multiple and num_atoms == 1:
         consensus += ambiguous
     elif (len(max_atoms) == 1) and ((float(max_size) / float(num_atoms)) >= threshold):
         consensus += max_atoms[0]
     else:
         consensus += ambiguous
     ```
   This matches every claim in the TestSpec/Evidence/algorithm doc.

2. **EMBOSS `cons`** — plurality cut-off "below which there is no consensus" corroborates the
   threshold→ambiguous semantics. (Reference accepted as cited; not re-fetched this session, used
   only for the conceptual plurality cut-off, which Biopython already pins exactly.)

3. **Wikipedia "Consensus sequence"** — definition ("most frequent residues … at each position in
   a sequence alignment") and IUPAC `N` = any base. Used for definition + default symbol rationale.

### Formula check

Commit predicate `unique-max ∧ num_atoms>0 ∧ (max_size/num_atoms ≥ threshold)`, else ambiguous —
matches Biopython's `(len(max_atoms)==1) and (max_size/num_atoms >= threshold)` with the all-gap
short-circuit. Inclusive `≥` matches the source. Gap set `{'-','.'}` matches. Length = longest
read matches `con_len`.

### Edge-case semantics check

Empty list → `""` (no columns); null → `ArgumentNullException`; all-gap column (`num_atoms==0`)
→ ambiguous, no div-by-zero; tie → ambiguous; sub-threshold → ambiguous; ragged reads span longest
read. All have a defined, sourced expected behaviour. INV-01..INV-05 are genuine properties.

### Independent cross-check (numbers)

Ran the actual reference implementation **Biopython 1.85** `dumb_consensus` on the test datasets
(gap-padded to equal length so `-` is skipped). Outputs (ambiguous `N`, threshold per case):

| Case | Reads | thr | Biopython output |
|------|-------|-----|------------------|
| M1 | ACGT,ACGT,ACGT | 0.5 | `ACGT` |
| M2 | ACGT,ACGT,ACGT,TCGT | 0.7 | `ACGT` |
| M3 | AC,AC,TC | 0.7 | `NC` |
| M4 | A,G | 0.5 | `N` |
| M5 | A-GT,A.GT,ACGT | 0.5 | `ACGT` |
| M6 | ACGT,ACG | 0.5 | `ACGT` |
| M7 | A-T,A-T | 0.5 | `ANT` |
| M9a | A,A,A,A,T | 0.7 | `A` |
| M9b | A,A,T | 0.7 | `N` |
| C1 | A,G (ambiguous=X) | 0.5 | `X` |

Biopython's own deprecation-notice example independently confirms the rule: `ACGT/ATGT/ATGT` →
`ANGT` at default 0.7 (col1 T=2/3≈0.667 < 0.7 → N). All ten outputs equal the TestSpec expected
values; every expected value is therefore externally sourced, not a code echo.

### Findings / divergences

Two **documented default-value divergences** from Biopython, both presentation-only and fully
reachable via parameters (so they do not change the decision rule):

- `threshold` default `0.5` (true simple-majority / EMBOSS plurality) vs Biopython's documented
  `0.7`. Pass `threshold: 0.7` to reproduce Biopython exactly (verified by M2/M3/M9).
- `ambiguous` default `'N'` (DNA IUPAC any-base) vs Biopython's `'X'` (protein). Pass `ambiguous:'X'`
  for Biopython parity (verified by C1).

These are registered assumptions in the spec/Evidence and are the reason for **PASS-WITH-NOTES**
rather than PASS. The `require_multiple` flag is intentionally not implemented (documented;
default `False` is the standard rule under test), which is out of scope and not a defect.

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs:861-926`.

### Formula realised correctly? (evidence)

- Null guard `ArgumentNullException.ThrowIfNull` (L866); empty → `""` (L868). ✓
- `length` = max read length (L871-875) = `con_len`. ✓
- Per-column reset of `counts`/`numAtoms` (L881-882); bounds check `pos >= read.Length → continue`
  (L886) = Biopython `n < len(record.seq)`. ✓
- Gap skip `c == '-' || c == '.'` after `ToUpperInvariant` (L888-889). ✓
- Tie detection via `maxCount` (count of residues sharing the max), L896-912 — equivalent to
  `len(max_atoms)`. ✓
- Commit predicate `maxCount==1 && numAtoms>0 && (double)maxSize/numAtoms >= threshold` (L918-920),
  else `ambiguous` (L922). The `numAtoms>0` short-circuits before the division, matching
  Biopython's `len(max_atoms)==1` guard for all-gap columns. ✓

### Cross-verification table recomputed vs code

The full unfiltered suite passes (below); all twelve unit tests assert exactly the Biopython-derived
values in the table above. Hand-recomputation of each case matches both the code and Biopython 1.85.

### Variant/delegate consistency

Only one canonical method. The MCP wrapper `AlignmentTools.ComputeConsensus`
(`src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs:348`) delegates to the canonical
method with defaults plus an equal-length precondition; no divergent consensus logic. Out of scope
for this unit.

### Test quality audit (HARD gate)

- **Sourced, not code-echoes:** every expected value re-derived from the Biopython rule and confirmed
  against Biopython 1.85 output (table above). A deliberately-wrong implementation (e.g. `MaxBy`
  tie-break, or `>` instead of `≥`) would fail M3/M4/M9.
- **No green-washing:** all assertions are exact `Is.EqualTo` full-string checks — no
  Greater/AtLeast/Contains/ranges, no widened tolerances, no skipped/ignored tests.
- **Coverage:** all five invariants and all decision-rule branches exercised — unanimous (M1),
  threshold-above (M2/M9a), threshold-below (M3/M9b), tie (M4), gap-skip for both `-` and `.` (M5),
  ragged length (M6), all-gap column (M7), empty list (M8), null (S1), lowercase normalization (S2),
  custom ambiguous symbol (C1). Both defaults exercised. Old weak `#region ComputeConsensus` block
  was already removed from `SequenceAssemblerTests.cs` (L222 marker).
- **Honest green:** full unfiltered suite `Failed: 0, Passed: 6532`; `dotnet build` 0 errors. The 4
  build warnings are in an unrelated file (`ApproximateMatcher_EditDistance_Tests.cs`), not changed
  here.

**Test-quality gate: PASS.**

### Findings / defects

None. No code or test change was required.

## Verdict & follow-ups

- **Stage A:** PASS-WITH-NOTES (two documented, parameter-reachable default divergences from Biopython).
- **Stage B:** PASS.
- **End-state:** ✅ CLEAN — no defect found; implementation and tests faithfully realise the
  externally-validated Biopython `dumb_consensus` decision rule.
- **Cross-check:** Biopython 1.85 reference implementation run this session; 10/10 datasets match.
- No follow-ups.
