# Validation Report: ANNOT-PROM-001 — Prokaryotic Promoter Motif Detection (-10 / -35 boxes)

- **Validated:** 2026-06-24   **Area:** Annotation
- **Canonical method(s):** `GenomeAnnotator.FindPromoterMotifs(string dnaSequence)`
  → `IEnumerable<(int position, string type, string sequence, double score)>`
  (`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs:584-638`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN

---

## Stage A — Description

### Sources opened & what they confirm

| Source | What it confirms (fetched live, not from memory) |
|--------|---------------|
| Wikipedia "Promoter (genetics)" | −35 box consensus **TTGACA**; −10 box consensus **TATAAT**; optimal spacer **17 bp** ("affects promoter strength by up to 600-fold"). Per-position frequencies — **−35:** T 69%, T 79%, G 61%, A 56%, C 54%, A 54%; **−10:** T 77%, A 76%, T 60%, A 61%, A 56%, T 82%. |
| Wikipedia "Pribnow box" | −10 element consensus **TATAAT**, "centered roughly ten base pairs upstream from the site of initiation of transcription"; alternative E. coli frequencies (Harley et al.) T 82%, A 89%, T 52%, A 59%, A 49%, T 89%. |
| Pribnow (1975) PNAS 72(3):784-788 | Original identification of the −10 ("Pribnow box"). (citation cross-referenced via Wikipedia primary refs) |
| Harley & Reynolds (1987) NAR 15(5):2343-2361 | E. coli promoter compilation underlying the per-position nucleotide occurrence probabilities. |

### Definition / convention check
- **Consensus correct & not swapped:** TTGACA → "-35 box", TATAAT → "-10 box". Confirmed against both Wikipedia articles. No −10/−35 swap; TATAAT (not TATATA) correct.
- **Per-position probabilities match exactly:** the values hard-coded in the source (69/79/61/56/54/54 for −35; 77/76/60/61/56/82 for −10) match the Wikipedia "Promoter (genetics)" frequency tables verbatim.
- **Coordinate convention:** authoritative convention is distance upstream of TSS as a negative number. The implementation's `position` field is a **0-based index into the input string** (per spec M01/M07 + tests). This is a substring-scan primitive, not a TSS-anchored predictor — consistently scoped and documented.
- **Spacer window (17 bp optimal):** confirmed in both sources; spec explicitly lists "spacing validation — not enforced" as a deliberate deviation. Documented, not a hidden defect.
- **Mismatch tolerance:** real promoters tolerate 2–3 mismatches; implementation uses **exact** matching to consensus substrings (full 6-mer + prefix-5/suffix-5/prefix-4). Declared deviation (a mismatch-tolerant detector would need a full PWM/HMM).

### Scoring formula
`score = sum(matched-position probabilities) / sum(all 6 consensus probabilities)`.
- −35 denominator = 69+79+61+56+54+54 = **373** (3.73).
- −10 denominator = 77+76+60+61+56+82 = **412** (4.12).

### Independent cross-check (hand-recomputed this session, Python)

| Variant | Positions | Calc | Code/Spec | Recomputed |
|---------|-----------|------|-----------|------------|
| TTGACA | 1–6 | 373/373 | 1.000 | 1.000 ✓ |
| TTGAC  | 1–5 | (69+79+61+56+54)/373 | 0.855 | 0.855 ✓ |
| TGACA  | 2–6 | (79+61+56+54+54)/373 | 0.815 | 0.815 ✓ |
| TTGA   | 1–4 | (69+79+61+56)/373 | 0.710 | 0.710 ✓ |
| TATAAT | 1–6 | 412/412 | 1.000 | 1.000 ✓ |
| TATAA  | 1–5 | (77+76+60+61+56)/412 | 0.801 | 0.801 ✓ |
| ATAAT  | 2–6 | (76+60+61+56+82)/412 | 0.813 | 0.813 ✓ |
| TATA   | 1–4 | (77+76+60+61)/412 | 0.665 | 0.665 ✓ |

All eight scores match the spec table and the code constants exactly.

### Stage A findings
Consensus sequences, per-position probabilities, spacer window, and the upstream-of-TSS convention all match authoritative external sources fetched live. The two divergences (no spacing enforcement; exact-substring rather than mismatch-tolerant matching) are explicitly declared in the TestSpec "Deviations from Literature" section and are legitimate scoping choices for a motif-scan primitive. **Stage A: PASS.**

---

## Stage B — Implementation

### Code path reviewed
`GenomeAnnotator.FindPromoterMotifs` — `GenomeAnnotator.cs:584-638`.

- Motif tables (`minus35Motifs`, `minus10Motifs`) hold the full 6-mer plus prefix-5, suffix-5, prefix-4 with the literature-derived scores above — values match Stage-A hand computation exactly.
- Input is upper-cased via `ToUpperInvariant()` → case-insensitive (M07).
- Each motif scanned with `for (i = 0; i <= seq.Length - motif.Length; i++)` and exact `Substring` compare; every occurrence yielded with its 0-based `position`. Loop bound is safe for empty / short input (`motif.Length > len` ⇒ no iterations), so empty and all-C sequences yield nothing (M05/M06).
- `-35 box` ↔ TTGACA, `-10 box` ↔ TATAAT — labels correct, not swapped.

### Cross-verification recomputed vs code
The 8-row variant table was reproduced directly from the code's constant arrays; all scores agree. Edge cases traced: no-motif → empty; empty → empty; one box only → only that box's variants; multiple same-type → all positions reported (M08: TTGACA at 0 and 11).

### Variant/delegate consistency
Single canonical static method; no `*Fast` variant or instance delegate to cross-check.

### Test quality audit
`GenomeAnnotator_PromoterMotif_Tests.cs` (20 methods → 28 cases via `[TestCase]`). Assertions check **exact** sourced values: positions, type strings, and scores to ±0.001, plus monotonic ordering (S02) and exact variant decomposition (S01: TATAAT → positions 0,0,1,0). Edge cases M05/M06/M07/M08 covered. Not tautological. Real-world C01 verifies both boxes with ~18 bp spacing at exact positions 4 and 28.

### Stage B findings
Implementation faithfully realises the validated description. **Stage B: PASS.**

---

## Verdict & follow-ups

- **Stage A: PASS** — consensus (−35 TTGACA / −10 TATAAT), per-position probabilities, 17 bp optimal spacer, and upstream-negative coordinate convention all confirmed against Wikipedia (Promoter genetics / Pribnow box), with Pribnow (1975) and Harley & Reynolds (1987) as the cited primary refs.
- **Stage B: PASS** — code constants and labels match sources; all 8 variant scores recomputed (Python) and verified; edge cases handled.
- **End-state: ✅ CLEAN** — no defect found. Documented deviations (spacing not enforced; exact-match not PWM) are declared scope limits of a motif-scan primitive, not errors.
- **Tests:** `--filter ~GenomeAnnotator_PromoterMotif_Tests` → **20 passed / 0 failed** (28 cases). No code changed.
