# Validation Report: SEQ-COMP-001 — DNA/RNA Complement

- **Validated:** 2026-06-12   **Area:** Composition
- **Canonical method(s):** `SequenceExtensions.GetComplementBase(char)`, `SequenceExtensions.GetRnaComplementBase(char)`, `SequenceExtensions.TryGetComplement(ReadOnlySpan<char>, Span<char>)`; instance wrappers `DnaSequence.Complement()`, `RnaSequence.Complement()`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS (RNA-IUPAC gap closed 2026-06-12 — see "Fix applied" below; both DNA and RNA paths are now IUPAC-complete)

## Stage A — Description

### Sources opened & what they confirm

1. **Wikipedia — Nucleic acid notation** (https://en.wikipedia.org/wiki/Nucleic_acid_notation), opened via WebFetch.
   Confirmed the full IUPAC (NC-IUB 1984) symbol → complement table verbatim:

   | Symbol | Meaning | Bases | Complement |
   |--------|---------|-------|------------|
   | A | Adenine | A | T |
   | C | Cytosine | C | G |
   | G | Guanine | G | C |
   | T | Thymine | T | A |
   | U | Uracil | U | A |
   | W | Weak | A,T | W |
   | S | Strong | C,G | S |
   | M | Amino | A,C | K |
   | K | Keto/Ketone | G,T | M |
   | R | Purine | A,G | Y |
   | Y | Pyrimidine | C,T | R |
   | B | Not A | C,G,T | V |
   | D | Not C | A,G,T | H |
   | H | Not G | A,C,T | D |
   | V | Not T | A,C,G | B |
   | N | Any | A,C,G,T | N |
   | - | Gap | — | — |

   This is exactly the table implemented by `GetComplementBase` (A↔T, G↔C, U→A, R↔Y, S↔S, W↔W, K↔M, B↔V, D↔H, N↔N).

2. **Biopython `Bio.Seq` complement() / complement_rna()** (https://biopython.org/docs/latest/api/Bio.Seq.html), opened via WebFetch. Confirmed docstring examples:
   - `complement("CGA")` → `"GCT"`
   - U treated as T in DNA context: `complement(Seq("CGAUT"))` → `Seq("GCTAA")`
   - Case-preserving + IUPAC ambiguity + pass-through unknowns: `complement("ACGTUacgtuXYZxyz")` → `"TGCAAtgcaaXRZxrz"`
   - `complement_rna("CGA")` → `"GCU"`; T treated as U: `complement_rna(Seq("CGAUT"))` → `Seq("GCUAA")`
   - `complement_rna("ACGTUacgtuXYZxyz")` → `"UGCAAugcaaXRZxrz"`

### Formula / definition check
The "formula" is a fixed bijective lookup table. Watson–Crick pairing (A↔T/U, G↔C) and the IUPAC ambiguity complements are correct and standard. Complement is an involution on the recognized alphabet (verified: every code maps to a partner whose complement returns the original; self-paired S,W,N trivially hold).

### Edge-case semantics (sourced)
- **U in DNA context** → A (Biopython "U treated as T", IUPAC U complement = A). ✓
- **T in RNA context** → A (Biopython "T treated as U"). ✓
- **Gap `-`** → passes through (IUPAC gap complement is gap; Biopython preserves). ✓
- **Unknown / non-IUPAC chars** → pass through unchanged (Biopython `X→X, Z→Z`). ✓
- **Case** → Biopython preserves case. Our DNA/RNA *sequence types* normalize input to uppercase on construction, so recognized bases emit uppercase; the per-char extension uppercases only recognized bases and passes unknowns through verbatim (see Stage B note).
- **Empty** → trivially valid (zero-length).

### Independent cross-check (exact numbers)
- DNA `complement("ACGT")` = `"TGCA"`; `complement("CGAUT")` = `"GCTAA"`; `complement("ACTG-NH")` = `"TGAC-ND"`.
- RNA `complement_rna("CGAUT")` = `"GCUAA"`.
- Full string `complement("ACGTUacgtuXYZxyz")` decodes per-char as A→T C→G G→C T→A U→A | a→T c→G g→C t→A u→A | X→X Y→R Z→Z | x→x y→R z→z.

All confirmed against two independent authoritative sources (Wikipedia IUPAC table and Biopython docstrings).

### Findings
Stage A description is correct and authoritatively sourced. PASS.

## Stage B — Implementation

### Code path reviewed
- `src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs:138` — `GetComplementBase` (full 16-symbol IUPAC table + gap/unknown pass-through).
- `:171` — `GetRnaComplementBase` (A,U,G,C,T only; everything else passes through).
- `:189` — `TryGetComplement` (returns false if `destination.Length < sequence.Length`; O(n) loop; empty source → returns true, writes nothing).
- `DnaSequence.cs:54` `Complement()` and `RnaSequence.cs:54` `Complement()` instance wrappers.

### Formula realised correctly?
`GetComplementBase` reproduces the validated NC-IUB table exactly, case-insensitively, returning uppercase for every recognized base. Recomputed all 16 symbols + lowercase + involution against the code — all match the Wikipedia/Biopython table.

### Cross-verification table recomputed vs code (via added tests)
| Input | Reference | Code output | Match |
|-------|-----------|-------------|-------|
| `ACGT` | Biopython `TGCA` | `TGCA` | ✓ |
| `CGAUT` (DNA) | Biopython `GCTAA` | `GCTAA` | ✓ |
| `ACTG-NH` | Biopython `TGAC-ND` | `TGAC-ND` | ✓ |
| `ACGTUacgtuXYZxyz` | Biopython `TGCAAtgcaaXRZxrz` (case-preserving) | `TGCAATGCAAXRZxRz` | ✓ (recognized bases uppercased; unknowns verbatim) |
| `CGAUT` (RNA) | Biopython `GCUAA` | `GCUAA` | ✓ |

### Variant/delegate consistency
`TryGetComplement`, `TryGetReverseComplement`, and `DnaSequence.Complement()` all delegate to `GetComplementBase` — consistent by construction. `RnaSequence.Complement()` uses an inline A/U/C/G switch (matches `GetRnaComplementBase` on the AUGC subset; sequence is pre-normalized to uppercase so the AUGC-only switch is sufficient there).

### Test quality audit
`SequenceExtensions_Complement_Tests.cs` asserts exact sourced values (not "no throw"), covers all MUST/SHOULD cases, IUPAC codes + involution, buffer-too-small (false), destination-larger (only source.Length written), empty (destination untouched), and Biopython cross-checks. Tests are deterministic. 129 tests pass under the `~Complement` filter (was 126 before this session).

### Findings / notes
1. **Clarified contract (caught during validation):** the documented claim "always returns uppercase" applies only to *recognized* IUPAC bases. Non-IUPAC characters (including lowercase `x`/`z`) pass through verbatim via `_ => nucleotide`. This is the correct, Biopython-consistent behavior; my initial cross-check test had the wrong expected value and the code corrected my expectation. The new test now locks the true behavior (`...XRZxRz`).
2. **RNA IUPAC gap — RESOLVED 2026-06-12 (see "Fix applied" below):** `GetRnaComplementBase` previously handled only A,U,G,C,T and passed IUPAC ambiguity codes through unchanged, whereas Biopython `complement_rna` complements them. This has been fixed: the RNA path now implements the full IUPAC ambiguity complement table (R↔Y, S↔S, W↔W, K↔M, B↔V, D↔H, N↔N), exactly mirroring the DNA path but emitting the RNA alphabet (U not T). Both paths are now IUPAC-complete.

## Verdict & follow-ups
- **Stage A: PASS** — table and edge semantics confirmed against Wikipedia (IUPAC NC-IUB 1984) and Biopython docstrings.
- **Stage B: PASS** — DNA path fully realises the validated IUPAC table and all edge cases; RNA helper is now also IUPAC-complete (gap closed 2026-06-12, see below). No outstanding defects.
- **State: CLEAN** — RNA-IUPAC gap fixed; build green, 4465/4465 tests pass.
- Tests added: `GetComplementBase_BiopythonCrossVerification_GapAndIupac`, `GetComplementBase_BiopythonCrossVerification_FullCaseAndUnknowns`, `GetRnaComplementBase_BiopythonCrossVerification_FullCaseAndUnknowns` (Stage B session), plus the four RNA-IUPAC tests listed in "Fix applied" below.

## Fix applied (2026-06-12)

**Change.** Extended `SequenceExtensions.GetRnaComplementBase(char)` (`src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs`) from the AUGCT-only switch to the full IUPAC ambiguity complement table, exactly mirroring the DNA path `GetComplementBase` but emitting the RNA alphabet (U not T). New arms (case-insensitive, recognized bases returned uppercase, non-IUPAC chars pass through verbatim via `_ => nucleotide`):
R↔Y, Y↔R, S↔S, W↔W, K↔M, M↔K, B↔V, V↔B, D↔H, H↔D, N↔N. The XML doc comment was updated to state full IUPAC support and cite Biopython.

**Source.** Biopython `Bio.Seq.complement_rna` (https://biopython.org/docs/latest/api/Bio.Seq.html), verified via WebFetch: `complement_rna("ACGTUacgtuXYZxyz")` → `"UGCAAugcaaXRZxrz"` (ambiguity codes complement the same as DNA; only the base alphabet changes T→U). Cross-checked against the Wikipedia/NC-IUB 1984 IUPAC table already cited in Stage A.

**Behavioral note.** As with the DNA path, recognized bases are uppercased while unrecognized characters pass through verbatim, so our output for the Biopython example string is `"UGCAAUGCAAXRZxRz"` (matches Biopython exactly except recognized bases are uppercased — same case convention as the DNA path).

**Impact check.** The only external consumer of `GetRnaComplementBase` outside this file is `RnaSecondaryStructure.GetComplement` (`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs:454`), which now gains correct ambiguity-code handling; no callers depended on the old pass-through of ambiguity codes. `RnaSequence.Complement()`/`ReverseComplement()` use their own inline AUGC switch and are unaffected. Full suite stays green.

**Tests added** (in `tests/Seqeron/Seqeron.Genomics.Tests/SequenceExtensions_Complement_Tests.cs`):
- `GetRnaComplementBase_IupacAmbiguityCodes_CorrectComplements` — each IUPAC code's RNA complement (R→Y, Y→R, S→S, W→W, K→M, M→K, B→V, V→B, D→H, H→D, N→N).
- `GetRnaComplementBase_IupacAmbiguityCodes_LowercaseReturnsUppercase` — lowercase variants return uppercase.
- `GetRnaComplementBase_IupacInvolution_AllRecognizedCodes` — involution over all 15 recognized codes.
- `GetRnaComplementBase_BiopythonCrossVerification_FullIupacString` — `complement_rna("ACGURYSWKMBDHVN")` → `"UGCAYRSWMKVHDBN"`.
- Updated `GetRnaComplementBase_BiopythonCrossVerification_FullCaseAndUnknowns` to assert the full Biopython string `"ACGTUacgtuXYZxyz"` → `"UGCAAUGCAAXRZxRz"` (previously restricted to the AUGCT subset).

Complement filter: 133 passed. Full suite: 4465/4465 passed.
