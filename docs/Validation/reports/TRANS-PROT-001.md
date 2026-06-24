# Validation Report: TRANS-PROT-001 — DNA/mRNA → Protein Translation

- **Validated:** 2026-06-24   **Area:** Translation
- **Canonical method(s):** `Translator.Translate(DnaSequence|RnaSequence|string, geneticCode, frame, toFirstStop)`, `Translator.TranslateSixFrames(DnaSequence, geneticCode)`, `Translator.FindOrfs(DnaSequence, geneticCode, minLength, searchBothStrands)`
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Core/Translator.cs` (codon table in `GeneticCode.cs`)
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/TranslatorTests.cs`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES

## Stage A — Description

### Sources opened & what they confirm
- **Rosalind PROT** (https://rosalind.info/problems/prot/) — translate an mRNA codon-by-codon 5'→3'. Stop codons (UAA, UAG, UGA) terminate translation and are **not** emitted; sequence after the stop is ignored.
  - Sample input: `AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA`
  - Sample output: `MAMAPRTEINSTRING`
- **Biopython `Bio.Seq.translate`** (biopython.org/docs Bio.Seq) — confirms the exact convention pair the library implements:
  - `stop_symbol` defaults to `'*'`.
  - `to_stop` defaults to `False` → translation continues past in-frame stops, each rendered as `'*'`.
  - `to_stop=True` → terminates at the first in-frame stop, and the stop symbol is **not** appended.
  - Ambiguous codons (`NNN`, `GCN`, `TAN`) → `'X'` (unknown amino acid), not an error.
- **Wikipedia: Reading frame** — single strand has 3 frames (offset 0/1/2); dsDNA has 6 (3 forward + 3 on the reverse complement read 5'→3').
- **NCBI Genetic Codes** (Evidence doc, re-checked counts) — Standard Table 1, Vert. Mito Table 2, Yeast Mito Table 3, Bacterial Table 11 mappings + start/stop sets.

### Semantics confirmed
- **Direction/reading:** triplet codons 5'→3', exact-table lookup. ✓
- **Stop handling (dual convention):** `toFirstStop=false` (default, = Biopython default) emits `'*'` and continues; `toFirstStop=true` (= Rosalind / Biopython `to_stop=True`) terminates and excludes the stop. ✓
- **Frames:** `Translate` uses 0-based {0,1,2}; `TranslateSixFrames` uses ±1/±2/±3 (forward offsets 0/1/2; negative = reverse-complement offsets 0/1/2, the Biopython/EMBOSS "alternative" convention documented in the method's `<remarks>`). ✓
- **T→U normalisation:** DNA `T` mapped to `U` before lookup (in both `TranslateSequence` and `GeneticCode.Translate`). ✓
- **Trailing incomplete codon:** loop `i + 3 <= len` drops any 1–2 nt remainder. ✓ (Biopython six_frame_translations / EMBOSS transeq triplet truncation.)
- **Ambiguous codon → 'X':** `GeneticCode.Translate` returns `'X'` when a codon is all-IUPAC but not in the table; throws for truly invalid characters. ✓ (matches Biopython.)
- **Start handling:** plain `Translate` translates the whole frame (no AUG requirement); the AUG-gated path is `FindOrfs`. Both are standard and cleanly separated. ✓

### Independent cross-checks (hand-computed)
- Rosalind PROT: `AUG·GCC·AUG·GCG·CCC·AGA·ACU·GAG·AUC·AAU·AGU·ACC·CGU·AUU·AAC·GGG·UGA` → `MAMAPRTEINSTRING` (`*` at UGA terminates). ✓
- Internal stop, default convention: `AUG·UAA·GCU·UAG·GGG` → `M*A*G`; `toFirstStop=true` → `M`. ✓
- Trailing partial codon: `AUG·GCU·CC` → `MA` (2-nt remainder dropped). ✓
- Ambiguous codon: `AUG·NNN·GCN·UAA` → `MXX*` (string overload). ✓
- Standard codon table contains all 64 codons (3 stops), matching NCBI Table 1.

### Findings / divergences
No biological divergence. The dual stop convention is internally consistent and the `toFirstStop=true` / default paths exactly match Rosalind and Biopython respectively.

## Stage B — Implementation

### Code path reviewed
- `Translator.TranslateSequence` (`Translator.cs:138-162`) — frame validation (`<0||>2` → `ArgumentOutOfRangeException`), `T→U`, codon loop `for(i=frame; i+3<=len; i+=3)`, `geneticCode.Translate`, `if(toFirstStop && aa=='*') break;` else append.
- `Translator.TranslateSixFrames` (`Translator.cs:115-136`) — frames +1/+2/+3 = forward offsets 0/1/2; -1/-2/-3 = reverse-complement offsets 0/1/2.
- `Translator.FindOrfsInSequence` (`Translator.cs:164-229`) — per-frame AUG-gated scan; emits on stop or end-of-sequence; minLength filter; reverse reported with negative frame; inclusive end at last base of stop codon.
- `GeneticCode.Translate` (`GeneticCode.cs:63-89`) — exact table lookup; all-IUPAC codon → `'X'`; otherwise throws.

### Formula realised correctly?
Yes — codon-by-codon 5'→3' with the validated stop semantics, frame offsets, T→U, and triplet truncation.

### Cross-verification recomputed vs code (tests run, all pass)
| Input | Convention | Expected | Code |
|-------|-----------|----------|------|
| `AUGGCC…GGGUGA` (Rosalind) | `toFirstStop:true` | `MAMAPRTEINSTRING` | ✓ |
| same | default | `MAMAPRTEINSTRING*` | ✓ |
| `AUGUAAGCUUAGGGG` | default | `M*A*G` | ✓ |
| `AUGUAAGCUUAGGGG` | `toFirstStop:true` | `M` | ✓ |
| `AUGGCUCC` | default | `MA` (2-nt remainder dropped) | ✓ |
| `AUGNNNGCNUAA` (string overload) | default | `MXX*` | ✓ |
| Insulin B CDS (UniProt P01308) | default | `FVNQHLCGSHLVEALYLVCGERGFFYTPKT` | ✓ |

### Variant/delegate consistency
- `TranslateSixFrames[+1]` == `Translate(frame:0)`; `[-(f+1)]` == `Translate(revComp, frame:f)` — covered by existing tests. ✓
- DNA/RNA/string overloads delegate to one `TranslateSequence`. ✓

### Test quality audit
Added 3 lock tests for previously implicit Stage-A edge cases:
`Translate_AmbiguousIupacCodon_ProducesX` (→ `MXX*`), `Translate_InternalStop_DefaultConventionEmitsStarAndContinues` (→ `M*A*G` / `M`), `Translate_TrailingPartialCodon_IsDropped` (→ `MA`). All assert exact sourced values. Translator filter: 35 passed.

### Findings / notes (non-defects)
- **N1 — Ambiguity path reachable only via the `string` overload.** The typed `DnaSequence`/`RnaSequence` constructors reject IUPAC ambiguity codes (only A/C/G/U[/T] are valid bases), so the `'X'` (unknown amino acid) translation path in `GeneticCode.Translate` is **unreachable** through the typed overloads — they throw at construction. The `Translate(string)` overload does not pre-validate the alphabet, so it is the only entry point where an ambiguous codon yields `'X'`. This is internally consistent (strict typed sequences vs. permissive raw-string translation matching Biopython) and is now documented by the lock test, hence PASS-WITH-NOTES rather than a defect.
- **N2 — Trailing partial codon silently dropped.** Matches EMBOSS transeq / Biopython six_frame_translations. (Recent Biopython additionally warns on partial codons for `cds=True`; not applicable here.) Defensible, documented convention.

## Verdict & follow-ups
- **Stage A: PASS** — semantics confirmed against Rosalind PROT, Biopython `Seq.translate`, Wikipedia Reading frame, and NCBI genetic-code tables.
- **Stage B: PASS-WITH-NOTES** — code faithfully realises the validated description; two non-defect notes recorded (ambiguity reachable only via string overload; trailing partial codon dropped). No code change needed.
- **State: CLEAN** — no defect; added 3 lock tests for under-covered Stage-A edge cases.
- **Tests:** Translator filter 35 passed; full `Seqeron.Genomics.Tests` suite 18211 passed / 0 failed.
