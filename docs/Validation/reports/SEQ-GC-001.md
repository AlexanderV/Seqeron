# Validation Report: SEQ-GC-001 — GC Content Calculation

- **Validated:** 2026-06-24   **Area:** Composition
- **Canonical method(s):**
  `SequenceExtensions.CalculateGcContent(ReadOnlySpan<char>)` (percentage 0–100),
  `CalculateGcFraction(ReadOnlySpan<char>)` (fraction 0–1);
  **new opt-in overload** `CalculateGcFraction(ReadOnlySpan<char>, GcAmbiguityMode {Remove,Ignore,Weighted})`
  plus its `string` overload (commit `6e900e92`);
  delegates `CalculateGcContentFast`/`CalculateGcFractionFast`, `DnaSequence.GcContent()`.
- **Stage A verdict:** 🟡 PASS-WITH-NOTES
- **Stage B verdict:** ✅ PASS

## Stage A — Description

**Sources opened (retrieved live, not from memory).**
- Wikipedia *GC-content*: percentage formula `(G+C)/(A+T+G+C) × 100%`.
- **Biopython `Bio/SeqUtils/__init__.py` (master, retrieved 2026-06-24** from
  `https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py`).
  Exact `gc_fraction` source:
  ```python
  gc = sum(seq.count(x) for x in "CGScgs")
  if ambiguous == "remove":
      length = gc + sum(seq.count(x) for x in "ATWUatwu")
  else:
      length = len(seq)
  if ambiguous == "weighted":
      gc += sum((seq.count(x)+seq.count(x.lower()))*_gc_values[x] for x in "BDHKMNRVXY")
  if length == 0: return 0
  return gc / length
  ```
  with `_gc_values`: G=C=S=1.0, A=T=U=W=0.0, M=R=Y=K=X=N=0.5, V=B=2/3, H=D=1/3.

**Formula check.** The default fraction = `(G+C)/(A+T+G+C(+U))` matches Wikipedia and the
Biopython `remove` denominator for standard nucleotides. The new `GcAmbiguityMode` overload
reproduces the three Biopython modes verbatim (numerator counts C,G,S; `Remove` length excludes
all other ambiguity codes; `Ignore`/`Weighted` length = `len(seq)`; `Weighted` adds the per-code
mean GC value for B,D,H,K,M,N,R,V,X,Y). Empty/zero-length → 0.

**Edge-case semantics.** All sourced: empty → 0; no valid nucleotides → 0; ambiguity codes per
the selected mode; U treated as a valid non-GC nucleotide; case-insensitive.

**Independent cross-check (exact numbers, re-implemented Biopython algorithm in Python).**

| Input | Mode | Biopython `gc_fraction` |
|-------|------|--------------------------|
| `GCAT` | remove | 0.5 |
| `GCSW` | remove | 0.75 |
| `GCATN` | remove | 0.5 |
| `GCATN` | ignore | 0.4 |
| `GCATN` | weighted | 0.5 |
| `VH` | weighted | 0.5 |
| `GCVBHD` | weighted | 0.6666666666666666 |
| `ATWWS` | remove | 0.2 |
| `NNNN` | ignore | 0.0 |
| `acgtSWn` | weighted | 0.5 |
| `""` | remove | 0 |

**Finding (note, not a defect).** The TestSpec (`tests/TestSpecs/SEQ-GC-001.md`, dated 2026-02-14)
still describes only the legacy default (A/T/G/C/U) behaviour and the historical S/W-exclusion
divergence; it predates the `6e900e92` opt-in overload and does not yet document the new
`GcAmbiguityMode` surface. The *code* is correct and the divergence is now resolved by opt-in
parity, so this is a documentation-staleness note → PASS-WITH-NOTES, not a failure.

## Stage B — Implementation

**Code path reviewed.**
`src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs`:
- `CalculateGcContent` :28–58 (percentage), `CalculateGcFraction` :81–111 (default fraction) —
  unchanged from prior validation; loop counts G/C into numerator+denominator, A/T/U into
  denominator only, all else skipped; `validCount==0 → 0`.
- **New** `GcAmbiguityMode` enum :136–155; weighted constants :160–164; overload
  `CalculateGcFraction(ReadOnlySpan<char>, GcAmbiguityMode)` :181–214; `string` overload :220–223;
  `WeightedGcValue` :225–233.

**Formula realised correctly.** The overload: numerator `gc` gets +1.0 for G/C/S; `strongWeakCount`
(Remove denominator) increments for A/T/G/C/S/W/U; `totalLength = sequence.Length` is the
Ignore/Weighted denominator; in `Weighted` mode the `default` branch adds `WeightedGcValue` for
the 10 IUPAC codes (V=B=2/3, H=D=1/3, M/R/Y/K/X/N=0.5, everything else 0). `length>0 ? gc/length : 0`.
This is byte-for-byte the Biopython algorithm: S/W are handled by the explicit switch cases (not
the weighted dict), exactly as Biopython keeps S/W out of its `"BDHKMNRVXY"` weighted loop.
`char.ToUpperInvariant` gives the case-insensitivity Biopython gets from the `count(x)+count(x.lower())`.

**Cross-verification recomputed against the code** — the 94 GC tests were executed (below); every
value in the table above reproduces. Traced the non-tested extras by hand through the switch:
`GCVBHD` weighted = (1+1 + 2/3+2/3+1/3+1/3)/6 = 4/6 = 0.6667 ✓; `ATWWS` remove = 1/(1+1+2+1)=0.2 ✓;
`acgtSWn` weighted = (1+1 + 1[S] + 0.5[n])/7 = 2.5? — recheck: gc base = c,g,S = 3; +n·0.5 = 3.5;
len 7 → 0.5 ✓ (matches Python reference).

**Variant/delegate consistency.** `*Fast` and `DnaSequence.GcContent()` forward to the canonical
Span methods; the new `string` overload guards null/empty → 0 then forwards to the Span overload.
No independent arithmetic to diverge.

**Numerical robustness.** Integer/`double` counts; denominator guarded (`length>0`, `validCount==0`,
`IsEmpty`/null). No overflow for realistic lengths.

**Test quality audit.**
- `ConventionCompatibility_OptIn_Tests.cs` :25–112 — 7 evidence-based cases for the new overload
  citing the exact Biopython lines: Remove unambiguous, Remove S/W semantics, Remove other-code
  exclusion, Ignore full-length denominator, Weighted N=0.5, Weighted V=2/3 & H=1/3, empty/null → 0,
  and a default-overload-unchanged guard. Assertions are exact sourced values via `Is.EqualTo(...).Within(Tol)`.
- `SequenceExtensions_CalculateGcContent_Tests.cs` (467 lines) — exact-value cases for the legacy
  methods incl. Biopython `GDVV`=100, `CCTGNN`=75, `GCAU`=50, `NNNNN`=0, plus biological references.
- Filtered run: **94 passed, 0 failed** (`~CalculateGcContent | ~ConventionCompatibility_OptIn`).

**Findings / defects.** None. The opt-in overload faithfully reproduces Biopython `gc_fraction` for
all three ambiguity modes; default behaviour is unchanged (verified by the dedicated guard test).

## Verdict & follow-ups

- **Stage A: 🟡 PASS-WITH-NOTES** — biology/maths correct and now matched to Biopython by opt-in
  parity; only note is that `tests/TestSpecs/SEQ-GC-001.md` is stale (predates the `GcAmbiguityMode`
  overload and still frames S/W as an unavoidable divergence). Recommend (non-blocking) adding a
  §1.6 to the TestSpec describing the three modes.
- **Stage B: ✅ PASS** — implementation realises the validated formula and the Biopython modes
  exactly; all cross-check numbers reproduce; delegates consistent; tests assert exact sourced values.
- **No code defects logged. No code changed this session.**
