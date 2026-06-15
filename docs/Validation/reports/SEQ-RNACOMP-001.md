# Validation Report: SEQ-RNACOMP-001 — RNA-specific Complement (per-base, IUPAC-complete)

- **Validated:** 2026-06-16   **Area:** Composition
- **Canonical method(s):** `SequenceExtensions.GetRnaComplementBase(char)` (`src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs:175`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

---

## Stage A — Description

### Sources opened & what they confirm (retrieved this session)

1. **Biopython `Bio/Data/IUPACData.py`** — https://raw.githubusercontent.com/biopython/biopython/master/Bio/Data/IUPACData.py (WebFetch, 2026-06-16). Verbatim:
   `ambiguous_rna_complement = {"A":"U","C":"G","G":"C","U":"A","M":"K","R":"Y","W":"W","S":"S","Y":"R","K":"M","V":"B","H":"D","D":"H","B":"V","X":"X","N":"N"}`
   `ambiguous_dna_complement = {"A":"T", ... ,"X":"X","N":"N"}` (identical except T↔U in the base alphabet).
2. **Biopython `Bio/Seq.py`** — https://raw.githubusercontent.com/biopython/biopython/master/Bio/Seq.py (WebFetch, 2026-06-16). Confirms `ambiguous_rna_complement["T"] = ambiguous_rna_complement["U"]` (so **T → A**), the docstring statement **"Any T in the sequence is treated as a U"**, and the worked example `Seq("CGAUT").complement_rna()` → `Seq('GCUAA')`.
3. **Biopython 1.80 `Bio.Seq` docs** — https://biopython.org/docs/1.80/api/Bio.Seq.html (WebFetch, 2026-06-16). Verbatim doctest: **`complement_rna("ACGTUacgtuXYZxyz")` → `'UGCAAugcaaXRZxrz'`** and `reverse_complement_rna("ACGTUacgtuXYZxyz")` → `'zrxZRXaacguAACGU'`. Lowercase preserved; X, Z pass through; Y→R.
4. **Bioinformatics.org SMS IUPAC table** — https://www.bioinformatics.org/sms/iupac.html (WebFetch, 2026-06-16). Code→bases→complement: R(A,G)→Y, Y(C,T)→R, S(G,C)→S, W(A,T)→W, K(G,T)→M, M(A,C)→K, B(C,G,T)→V, D(A,G,T)→H, H(A,C,T)→D, V(A,C,G)→B, N→N. U→A. (U substitutes for T in RNA.)
5. **NC-IUB 1984 (Cornish-Bowden, NAR 13(9):3021)** — originating standard for ambiguity codes; encoded by sources 1 and 4.

### Formula check

The "formula" is a fixed IUPAC lookup. The complement of an ambiguity set S is the set of complements of S's members. Hand-verified each ambiguity code against its member set with Watson–Crick A↔U, C↔G:
- R={A,G}→{U,C}=Y ✓; Y={C,U}→{G,A}=R ✓; S={G,C}→{C,G}=S ✓; W={A,U}→{U,A}=W ✓;
- K={G,U}→{C,A}=M ✓; M={A,C}→{U,G}=K ✓; B={C,G,U}→{G,C,A}=V ✓; D={A,G,U}→{U,C,A}=H ✓;
- H={A,C,U}→{U,G,A}=D ✓; V={A,C,G}→{U,G,C}=B ✓; N→N ✓.
All eleven match sources 1 and 4 exactly.

### Edge-case semantics

- **T → A** (T treated as U): sourced verbatim from source 2/3.
- **U → A, A → U** (RNA alphabet emits U, never T): source 1, 3.
- **Non-IUPAC pass-through** (gaps `-`/`.`, digits, Z, space, x, z): source 3 shows X, Z and lowercase passing through unchanged.
- **Case**: Biopython *preserves* input case; the repo *uppercases recognized bases* (documented Assumption #1, mirrors DNA sibling SEQ-COMP-001). This is the only divergence and affects letter casing only, not the complement identity.

### Independent cross-check (numbers)

`complement_rna("ACGTUacgtuXYZxyz")` (Biopython doc 1.80, verbatim) = `UGCAAugcaaXRZxrz`. Per-base:
A→U, C→G, G→C, T→A, U→A | a→u, c→g, g→c, t→a, u→a | X→X, Y→R, Z→Z | x→x, y→r, z→z. Hand-derived per-base mapping agrees with both the Biopython table (source 1) and the IUPAC table (source 4). Also `Seq("CGAUT").complement_rna()`=`GCUAA` independently confirms T→A.

### Findings / divergences

Stage A is sound. The description (algorithm doc, Evidence, TestSpec) accurately reflects all retrieved sources, including the single documented case-normalization divergence. **Stage A: PASS.**

---

## Stage B — Implementation

### Code path reviewed

`GetRnaComplementBase(char)` — `SequenceExtensions.cs:175-194`, a `switch` expression (case-insensitive arms, AggressiveInlining).

### Formula realised correctly?

The switch maps exactly: A/a→U, U/u→A, G/g→C, C/c→G, T/t→A, R/r→Y, Y/y→R, S/s→S, W/w→W, K/k→M, M/m→K, B/b→V, D/d→H, H/h→D, V/v→B, N/n→N, `_`→passthrough. This is **byte-for-byte identical** to Biopython `ambiguous_rna_complement` (source 1) plus the T→A rule (source 2/3), modulo the repo's uppercase-on-recognized convention. X is not an explicit arm but falls through `_ => nucleotide`, yielding X→X — coincides with Biopython's explicit X→X.

### Cross-verification table recomputed vs code

| Input | Source-expected (Biopython/IUPAC) | Repo (uppercase conv.) | Code output |
|-------|-----------------------------------|------------------------|-------------|
| A,U,C,G | U,A,G,C | U,A,G,C | U,A,G,C ✓ |
| T,t | A,A | A,A | A,A ✓ |
| R,Y,S,W,K,M,B,D,H,V,N | Y,R,S,W,M,K,V,H,D,B,N | same | same ✓ |
| a,u,c,g (lc) | u,a,g,c (Biopython) | U,A,G,C (repo) | U,A,G,C ✓ |
| `"ACGTUacgtuXYZxyz"` | `UGCAAugcaaXRZxrz` (Biopython) | `UGCAAUGCAAXRZxRz` | `UGCAAUGCAAXRZxRz` ✓ |
| `-`,`.`,`5`,`Z`,`z`,space | passthrough | passthrough | passthrough ✓ |

All values trace to externally retrieved sources, not to the code.

### Variant/delegate consistency

Single canonical method; no `*Fast`/instance variants. Contrasted with `GetComplementBase` (DNA, line 138): A→T vs A→U — consistent with the RNA-vs-DNA distinction (source 3: `complement` "TGC" vs `complement_rna` "UGC").

### Test quality audit (HARD gate)

File: `SequenceExtensions_GetRnaComplementBase_Tests.cs` (9 tests).
- **Sourced, not code-echoes:** every expected value (M1–M6, S1) is an exact literal traceable to Biopython/IUPAC sources retrieved this session. A deliberately-wrong impl (e.g. A→T, or T preserved, or any ambiguity swap) would fail M1/M2/M3/M5/M6.
- **No green-washing:** all known-value assertions use `Is.EqualTo` exact equality; M5 asserts the full exact string. No widened tolerances, no skips, no Greater/Contains where exact is known.
- **Coverage:** standard pairing, T→A (upper+lower), all 11 ambiguity codes, lowercase→uppercase, full-alphabet worked string, RNA≠DNA, non-IUPAC pass-through (incl. gaps/digit/Z/space/lowercase z), involution property (S2), no-T-emitted invariant (C1). All Stage-A branches and edge cases exercised. The sole public method has full coverage; no error/throw cases exist by contract (no `char` input throws).
- **C1** uses `Is.Not.EqualTo('T')` — appropriate as a supplementary invariant (INV-2); the exact values it relates to are already locked by M1/M3, so this is not a weakened substitute for an exact assertion.
- **Honest green:** full unfiltered suite **Failed: 0, Passed: 6573** (1 unrelated benchmark `MFE_Benchmark_AllScenarios` explicitly skipped). `dotnet build` 0 errors; the 4 build warnings are pre-existing NUnit2007 warnings in an unrelated file (`ApproximateMatcher_EditDistance_Tests.cs`), not in this unit's changed files (no files were changed this session).

**Test-quality gate: PASS.**

### Findings / defects

None. Code, tests, and description are mutually consistent and all match the external authoritative sources. **Stage B: PASS.**

---

## Verdict & follow-ups

- **Stage A: PASS.** Description matches Biopython `ambiguous_rna_complement`, the bioinformatics.org IUPAC table, and Biopython doctest examples; T→A and case-normalization divergence correctly documented.
- **Stage B: PASS.** Implementation is byte-for-byte the validated table; tests assert exact sourced values and cover all branches/edge cases.
- **Test-quality gate: PASS.**
- **End-state: ✅ CLEAN.** No defect found; full suite green (6573/0).
- No follow-ups; no findings-register entry required.
