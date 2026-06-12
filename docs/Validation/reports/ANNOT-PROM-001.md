# Validation Report: ANNOT-PROM-001 — Prokaryotic Promoter Motif Detection (-10 / -35 boxes)

- **Validated:** 2026-06-12   **Area:** Annotation
- **Canonical method(s):** `GenomeAnnotator.FindPromoterMotifs(string dnaSequence)`
  → `IEnumerable<(int position, string type, string sequence, double score)>`
  (`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs:461`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN

---

## Stage A — Description

### Sources opened & what they confirm

| Source | What it confirms |
|--------|------------------|
| Wikipedia "Promoter (genetics)" | -35 box consensus **TTGACA**; -10 box consensus **TATAAT**; optimal spacer **17 bp** between -35 and -10; per-position occurrence frequencies — **-35:** T 69%, T 79%, G 61%, A 56%, C 54%, A 54%; **-10:** T 77%, A 76%, T 60%, A 61%, A 56%, T 82%; coordinate convention = positions upstream of the transcription start site are **negative numbers counting back from -1**. |
| Wikipedia "Pribnow box" | -10 element consensus **TATAAT**, centered ~10 bp upstream of the initiation site; discovered by **David Pribnow (March 1975, PNAS)**; alternative E. coli frequencies (Harley et al.) T 82%, A 89%, T 52%, A 59%, A 49%, T 89%. |
| Pribnow (1975) PNAS 72(3):784-788 | Original identification of the -10 ("Pribnow box") at an early T7 RNA-polymerase binding site. |
| Harley & Reynolds (1987) NAR 15(5):2343-2361 | E. coli promoter compilation underlying the per-position nucleotide occurrence probabilities. |

### Definition / convention check
- **Consensus correct & not swapped:** TTGACA → "-35 box", TATAAT → "-10 box". Confirmed against both Wikipedia articles. No TATATA/TATAAT confusion; no -10/-35 swap.
- **Coordinate convention:** authoritative convention is distance upstream of TSS as negative. The implementation does **not** report -10/-35 distances; its `position` field is a **0-based index into the input string** (documented in spec M01/M07 and tests). This is a substring-search/scan tool, not a TSS-anchored predictor — consistent and clearly scoped.
- **Spacer window (17 ± 1 bp, 15–19):** confirmed in literature; the spec explicitly lists "spacing validation — not enforced" as a deliberate deviation (independent motif search). Documented, not a hidden defect.
- **Mismatch tolerance:** real promoters tolerate 2–3 mismatches; the implementation uses **exact** matching to consensus substrings (full 6-mer + prefix-5/suffix-5/prefix-4 variants). Documented deviation (would require a PWM/HMM).

### Scoring formula
`score = sum(matched-position probabilities) / sum(all 6 consensus probabilities)`.
- -35 denominator = 69+79+61+56+54+54 = **373** (3.73).
- -10 denominator = 77+76+60+61+56+82 = **412** (4.12).

### Worked example & independent cross-check (hand-computed)

Promoter example: `TTGACA …(17 bp)… TATAAT`. Both consensus 6-mers are detected as full hits (score 1.000); spacing is reported via positions (not enforced). Variant scores recomputed by hand:

| Variant | Positions | Calc | Score |
|---------|-----------|------|-------|
| TTGACA | 1–6 | 373/373 | 1.000 ✓ |
| TTGAC  | 1–5 | (69+79+61+56+54)/373 = 319/373 | 0.855 ✓ |
| TGACA  | 2–6 | (79+61+56+54+54)/373 = 304/373 | 0.815 ✓ |
| TTGA   | 1–4 | (69+79+61+56)/373 = 265/373 | 0.710 ✓ |
| TATAAT | 1–6 | 412/412 | 1.000 ✓ |
| TATAA  | 1–5 | (77+76+60+61+56)/412 = 330/412 | 0.801 ✓ |
| ATAAT  | 2–6 | (76+60+61+56+82)/412 = 335/412 | 0.813 ✓ |
| TATA   | 1–4 | (77+76+60+61)/412 = 274/412 | 0.665 ✓ |

All eight scores match the spec table and the code constants exactly.

### Stage A findings
Consensus sequences, per-position probabilities, spacer window, and coordinate convention all match authoritative sources. The two divergences (no spacing enforcement; exact-substring instead of mismatch-tolerant matching) are explicitly declared in the TestSpec "Deviations from Literature" section and are legitimate scoping choices for a motif-scan primitive. **Stage A: PASS.**

---

## Stage B — Implementation

### Code path reviewed
`GenomeAnnotator.FindPromoterMotifs` — `GenomeAnnotator.cs:461-515`.

- Motif tables (`minus35Motifs`, `minus10Motifs`) hold the full 6-mer plus prefix-5, suffix-5, prefix-4 with the literature-derived scores listed above — values match Stage-A hand computation exactly.
- Input is upper-cased via `ToUpperInvariant()` → case-insensitive (M07).
- Each motif scanned with `for (i = 0; i <= len - motif.Length; i++)` and exact substring compare; every occurrence yielded with its 0-based `position`. Loop bound is safe for empty / short input (no scan when `motif.Length > len`), so empty string and all-C sequences yield nothing (M05/M06).
- `-35 box` ↔ TTGACA, `-10 box` ↔ TATAAT — labels correct, not swapped.

### Cross-verification recomputed vs code
The Stage-A worked-example table was reproduced directly from the code's constant arrays; all 8 variant scores agree. Edge cases traced: no-motif → empty; empty → empty; one box only → only that box's variants; multiple same-type → all positions reported (M08: TTGACA at 0 and 11).

### Variant/delegate consistency
Single canonical static method; no `*Fast` variant or instance delegate to cross-check.

### Test quality audit
`GenomeAnnotator_PromoterMotif_Tests.cs` (20 methods → 28 cases with TestCases). Assertions check **exact** sourced values: positions, type strings, and scores to 0.001, plus monotonic ordering (S02) and exact variant decomposition (S01). Edge cases M05/M06/M07/M08 covered. Not tautological. Real-world C01 verifies both boxes with ~18 bp spacing.

### Stage B findings
Implementation faithfully realises the validated description. **Stage B: PASS.**

---

## Verdict & follow-ups

- **Stage A: PASS** — consensus (-35 TTGACA / -10 TATAAT), probabilities, 17 bp spacer, and upstream-negative coordinate convention confirmed against Wikipedia (Promoter genetics / Pribnow box), Pribnow (1975), Harley & Reynolds (1987).
- **Stage B: PASS** — code constants and labels match sources; all 8 variant scores recomputed and verified; edge cases handled.
- **End-state: ✅ CLEAN** — no defect found. Documented deviations (spacing not enforced; exact-match not PWM) are declared scope limits of a motif-scan primitive, not errors.
- **Tests:** `--filter ~Promoter` → 28 passed / 0 failed. Full suite `Seqeron.Genomics.Tests` → **4461 passed, 0 failed**. No code changed.
