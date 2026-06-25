# Validation Report: TRANS-CODON-001 — Codon → amino-acid translation tables

- **Validated:** 2026-06-24   **Area:** Translation
- **Canonical method(s):** `GeneticCode.Translate`, `IsStartCodon`, `IsStopCodon`, `GetCodonsForAminoAcid`, `GetByTableNumber` (+ static codes Standard/VertebrateMitochondrial/YeastMitochondrial/BacterialPlastid)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **NCBI Genetic Codes** (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) — fetched live. Retrieved the verbatim `AAs`/`Starts`/`Base1`/`Base2`/`Base3` strings for tables 1, 2, 3, 11. These match the strings recorded in `docs/Evidence/TRANS-CODON-001-Evidence.md` character-for-character.

### Method: decode-then-diff (not hand transcription)
Rather than trusting a hand-typed 64-codon table, I decoded the NCBI strings programmatically (Base1/Base2/Base3 → codon, with T→U; `AAs[i]` → amino acid; `Starts[i]=='M'` → start; `AAs[i]=='*'` → stop) and diffed the result against the codon tables produced by the code's own construction logic (Standard table + per-table overrides). Each table = exactly 64 codons.

NCBI-decoded start/stop sets:
| Table | Starts (decoded) | Stops (decoded) |
|-------|------------------|-----------------|
| 1 Standard | AUG, CUG, UUG | UAA, UAG, UGA |
| 2 Vert. Mito | AUA, AUC, AUG, AUU, GUG | AGA, AGG, UAA, UAG |
| 3 Yeast Mito | AUA, AUG, GUG | UAA, UAG |
| 11 Bact./Plastid | AUA, AUC, AUG, AUU, CUG, GUG, UUG | UAA, UAG, UGA |

### Key cross-check entries (by hand, against NCBI)
- Table 1: AUG=M (start), UAA/UAG/UGA=* (stops), 6-fold Leu {UUA,UUG,CUU,CUC,CUA,CUG}, 6-fold Ser {UCU,UCC,UCA,UCG,AGU,AGC}, 6-fold Arg {CGU,CGC,CGA,CGG,AGA,AGG}.
- Table 2 differences from std: AGA/AGG=* (not R), AUA=M (not I), UGA=W (not *) — confirmed.
- Table 3 differences from std: CUN=T (CUU/CUC/CUA/CUG, not L), AUA=M, UGA=W — confirmed.
- Table 11: identical codon→AA table to Standard; 7 alternative starts per Starts line — confirmed.

### Edge-case semantics
Length≠3 / empty / null → ArgumentException; non-IUPAC nucleotides (X, Z, digits) → ArgumentException; IUPAC ambiguity codons (N, R, Y, …) → 'X'. The IUPAC→'X' behavior matches Biopython/EMBOSS ambiguous-codon handling and is the established convention for this library.

### Findings / divergences
Minor description nuance (not a defect): the Evidence/TestSpec "Known Failure Modes" table lists `NNN` under "Unknown codon → Exception", but the implemented and tested behavior is `NNN → 'X'` (valid IUPAC, untranslatable). Code and tests are internally consistent and follow the documented IUPAC-X rule in the source comments; only the older Evidence prose under-specifies it. No correctness impact.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Core/GeneticCode.cs`:
- `Translate` (L63–79): length guard, `ToUpperInvariant().Replace('T','U')` normalization, dict lookup, IUPAC→'X' fallback, else throw.
- `CreateStandardCodonTable` (L230–319): 64 entries.
- Per-table overrides: Vert. Mito L180–183 (AUA=M, UGA=W, AGA=*, AGG=*); Yeast Mito L199–204 (CUN=T, AUA=M, UGA=W); Bacterial L214–225 (table unchanged, only starts differ).
- Start/stop sets L168–169, 187–188, 208–209, 221–222.

### Formula realised correctly?
Diff of code-constructed tables vs NCBI-decoded tables: **0 differences** for all four tables; all 64 codons present in each. Hardcoded start/stop sets in the source equal the NCBI-decoded sets exactly (verified above).

### Cross-verification vs code (tests run)
`GeneticCodeTests`: **50 passed, 0 failed** (dotnet test, net10.0). Includes the full-64-codon NCBI comparison (`Translate_CompleteStandardCodonTable_MatchesNcbi`), all three stops, AUG=M, exact degeneracy sets (M, W, L×6, S×6, R×6, I×3, stop×3), all table-2/3/11 differences, exact start/stop set counts and membership, DNA/RNA normalization, case-insensitivity, and invalid-input exceptions.

### Variant/delegate consistency
`GetByTableNumber(1/2/3/11)` returns the same singletons used by the static properties (tested). `IsStartCodon`/`IsStopCodon` read the same normalized sets `Translate` is consistent with. `GetCodonsForAminoAcid` is the inverse of the codon table (case-insensitive AA input).

### Test quality audit
Assertions check exact sourced values and exact set membership/counts, not tautologies; deterministic; cover every Stage-A edge case. Strong.

## Verdict & follow-ups
PASS / PASS. State CLEAN — no defect. Tables 1, 2, 3, 11 independently re-derived from live NCBI strings and matched the code with zero diffs. Only follow-up (optional, non-blocking): update the Evidence "Known Failure Modes" row so `NNN` reads "→ 'X' (IUPAC ambiguous)" to match the implemented/tested behavior. No code change made.
