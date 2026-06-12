# Validation Report: ANNOT-GENE-001 — Prokaryotic Gene Prediction (Shine-Dalgarno RBS + ORF)

- **Validated:** 2026-06-12   **Area:** Annotation
- **Canonical method(s):** `GenomeAnnotator.PredictGenes(dna, minOrfLength, prefix)`, `GenomeAnnotator.FindRibosomeBindingSites(dna, upstreamWindow, minDistance, maxDistance)` (supported by `FindOrfs`)
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs`
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/GenomeAnnotator_Gene_Tests.cs`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES
- **End state:** CLEAN (no defect found; no code change required)

---

## Stage A — Description

### Sources opened & what they confirm

| Source | Confirms |
|--------|----------|
| Wikipedia "Shine-Dalgarno sequence" | Six-base consensus = **AGGAGG** (E. coli AGGAGG**U**); anti-SD at 16S rRNA 3' terminus = **5'-YACCUCCUUA-3'** (complementary to AGGAGGU); SD "generally located around **8 bases upstream** of the start codon AUG". |
| Chen et al. (1994) NAR 22(23):4953-4957 (PMC523762 / academic.oup.com) | "Measurements of protein expression demonstrated an **optimal aligned spacing of 5 nt** for both series." Aligned spacing "corresponds naturally to the spacing between the 3'-end of the 16S rRNA and the P-site" — i.e. measured from the **3' end of the SD core** to the initiation codon. |
| Wikipedia "Ribosome-binding site" / "Gene prediction" | RBS = purine-rich region ~5' of the initiation codon directing the ribosome; prokaryotic gene finding is ORF-based with signal (RBS/SD) detection; prokaryotic start codons ATG/GTG/TTG; stop codons TAA/TAG/TGA. |
| Shine & Dalgarno (1975) | Original identification of the AGGAGGU polypurine tract complementary to 16S rRNA 3' end. |

### Defined-semantics check

- **SD consensus AGGAGG** — confirmed canonical (the AGGAGGU complementary to anti-SD CCUCCU / 5'-ACCUCCUUA). Correct.
- **SD variants** {AGGAGG, GGAGG, AGGAG, GAGG, AGGA} — all are contiguous substrings of the AGGAGG consensus core (degenerate/truncated SDs are biologically real; weaker SDs use a subset of the core). Sourced as "substrings of consensus" — acceptable.
- **Spacer window 4–15 bp, optimal 5 nt** — optimal 5 nt is Chen (1994); the broader 4–15 bp functional window is the TestSpec's sourced design choice (Wikipedia's "~8 bases upstream" sits squarely inside it). Reasonable and sourced.
- **Aligned-spacing reference point** — measured from the **3' (last) base of the SD core to the start codon**, NOT from the 5' end. This matches Chen (1994). (Key check requested by protocol: no wrong-end measurement, no anti-SD/SD confusion.)
- **Gene = RBS + downstream ORF** — `PredictGenes` is ORF-based (ATG/GTG/TTG … stop, both strands, ≥ minOrfLength aa); `FindRibosomeBindingSites` locates SD signals upstream of forward-strand ORFs within the spacer window. Correct model for prokaryotes (no introns).

### Worked examples (hand-computed)

`CreateSequenceWithSd` layout: `C×10 (padding) + SD + C×spacer + ATG…TAA`. SD found at genomic index 10; ORF start at `10 + len(SD) + spacer`. Code computes `distanceToStart = orfStart − genomicPos − len(SD)` = spacer.

- **Positive (with SD):** AGGAGG, spacer 8 → ORF start at index 24, `distanceToStart = 24 − 10 − 6 = 8` ∈ [4,15] → **predicted RBS (score 6/6 = 1.0)**. ✓
- **Optimal (Chen 1994):** AGGAGG, spacer 5 → `distanceToStart = 5` → detected. ✓
- **Multiple upstream (C3):** `C×10 + AGGAGG + C×4 + GAGG + C×5 + ATG…`; ATG at 29. AGGAGG spacing = 29−10−6 = 13; GAGG spacing = 29−20−4 = 5 → both reported. ✓ (hand-verified with a script)
- **Negative (no SD):** ORF with non-purine-rich upstream (e.g. poly-C/A) → no motif match → **no RBS** (lower confidence / not reported). An SD with no downstream ORF (`FindRibosomeBindingSites_NoOrfs_ReturnsEmpty`) → empty. ✓

### Findings / divergences
None biological. The 4–15 bp window is a sourced design parameter rather than a single literature constant; optimal 5 nt and the AGGAGG consensus are exactly sourced.

---

## Stage B — Implementation

### Code path reviewed
- `FindRibosomeBindingSites` — `GenomeAnnotator.cs:246-284`
- `PredictGenes` — `GenomeAnnotator.cs:289-321`
- `FindOrfs` / `FindOrfsInFrame` — `GenomeAnnotator.cs:90-212`
- Start/stop codon sets — `:52-63`

### Realised correctly?
- SD consensus + variants at `:253` = `{AGGAGG, GGAGG, AGGAG, GAGG, AGGA}` — matches description; comment "(binds to 3' end of 16S rRNA)" correctly states anti-SD complementarity (no SD/anti-SD confusion).
- **Aligned spacing** at `:272`: `distanceToStart = orf.Start − genomicPos − motif.Length`. For a motif of length L at `genomicPos`, its 3' base is at `genomicPos+L−1`; the spacer count between that base and `orf.Start` is exactly `orf.Start − genomicPos − L`. This is the **correct 3'-end-to-start-codon measurement** per Chen (1994) — measured from the right end. ✓
- Window filter `:274`: `distanceToStart ∈ [minDistance, maxDistance]` (default 4..15). ✓
- Score `:276`: `motif.Length / 6.0` (normalised to consensus length 6). AGGAGG → 1.0. ✓
- `PredictGenes`: Type "CDS" (`:317`), strand from `IsReverseComplement` (`:303`), sequential IDs `{prefix}_{n:D4}` (`:313`), `protein_length = ProteinSequence.TrimEnd('*').Length` (`:308`) = (End−Start)/3 − 1 (excludes stop). ✓

### Cross-verification vs code (tests executed)
| Case | Expected | Result |
|------|----------|--------|
| AGGAGG spacer 8 detected | yes, exact "AGGAGG", score 1.0 | PASS |
| Optimal spacer 5 (Chen) | detected | PASS |
| Min 4 / Max 15 / Beyond 16 | detect / detect / not detect | PASS |
| Too close (2) | full AGGAGG filtered | PASS |
| Shorter variants GGAGG/AGGAG/GAGG/AGGA | detected | PASS |
| Multiple upstream (spacing 13 & 5) | both | PASS |
| protein_length = (End−Start)/3 − 1 | exact | PASS |
| GTG/TTG starts, both strands, overlaps, long ORF | as specified | PASS |

### Edge cases
Empty/null → empty (`FindOrfs` guards null/empty); ATG-only / no-stop → empty; no-start → empty; SD too close (<4) / too far (>15) → filtered; multiple candidate starts handled by `FindOrfsInFrame` accumulating all starts before a stop.

### Findings / notes (PASS-WITH-NOTES, not defects)
1. **RBS only on forward-strand ORFs** — `FindRibosomeBindingSites` filters `o => !o.IsReverseComplement` (`:257`), so it never reports SD signals for reverse-strand genes. This is a scope limitation of the RBS helper (PredictGenes itself does cover both strands); it does not affect any tested/contracted behaviour and is not a correctness error.
2. **4–15 bp window is a design parameter**, not a single literature constant (optimal 5 nt is the sourced point value). Documented in the TestSpec.

---

## Verdict & follow-ups
- **Stage A: PASS** — AGGAGG consensus, anti-SD 5'-YACCUCCUUA-3', ~8 nt typical distance, and optimal 5 nt aligned spacing all confirmed against authoritative sources.
- **Stage B: PASS-WITH-NOTES** — implementation faithfully realises the validated description; aligned spacing is measured from the SD 3' end to the start codon (correct), worked examples reproduce, all tests pass. Notes are scope observations, not defects.
- **End state: CLEAN** — no defect; no code or test changes made.
- **Tests:** `--filter FullyQualifiedName~Gene` → 490 passed, 0 failed. Full suite → **4461 passed, 0 failed** (matches baseline).
