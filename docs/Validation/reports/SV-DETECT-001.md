# Validation Report: SV-DETECT-001 — Structural Variant Detection from Paired-End Mapping signatures

- **Validated:** 2026-06-15   **Area:** StructuralVar
- **Canonical method(s):** `StructuralVariantAnalyzer.ClassifySV`, `StructuralVariantAnalyzer.DetectSVs`, `StructuralVariantAnalyzer.FindDiscordantPairs` (+ helper `IsConcordantOrientation`, `ClusterDiscordantPairs`)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS-WITH-NOTES (after in-session fix of a real defect)

## Stage A — Description

### Sources opened this session (and what they confirm)

1. **Medvedev, Stanciu & Brudno (2009), Nat Methods 6(11s):S13–S20** — author PDF `cs.toronto.edu/~brudno/nmeth_review09.pdf`, text extracted locally with `pypdf` and read verbatim:
   - Deletion: *"A mate pair that spans an isolated deletion event maps to the corresponding regions of the reference, but the mapped distance is greater than the insert size."* → DEL = span **greater** than insert size. ✓
   - Insertion: *"Conversely, if the event is an insertion, then the distance is smaller."* → INS = span **smaller**. ✓
   - Inversion: *"A mate pair that spans either (but not both) of its breakpoints will map to the reference with the orientation of the read, lying within the inversion, flipped."* → INV = same-strand (flipped). ✓
   - Linking/translocation: *"Other linking signatures can connect regions that are arbitrarily distant or even on different chromosomes."* → cross-chromosome = translocation. ✓
   - Corner case: *"the basic insertion signature does not appear when the size of the insertion is greater than the insert size of the sequenced fragment, and it does not indicate the inserted sequence itself."* ✓
2. **BreakDancer README** (`raw.githubusercontent.com/genome/breakdancer/master/README`), fetched this session:
   - `-c` default **3**; *"upper: mean + std * threshold … lower: mean − std * threshold"*; *"the upper and the lower separation threshold would be: mean + 3 std and mean − 3 std respectively."* → bounds = μ ± c·σ, default c = 3. ✓
   - `-r` = *"minimum number of read pairs required to establish a connection [2]"* → default 2. ✓
   - SV codes DEL/INS/INV/ITX/CTX confirmed. ✓
3. **cureffi.org (FR proper-pair / SAM FLAG 0x02)**, fetched this session: FR is the proper orientation; *"If they instead align RF, FF or RR, that's a problem."* — i.e. **RF is abnormal, alongside FF/RR.**
4. **DELLY (Rausch et al. 2012, Bioinformatics 28(18):i333)** and **SVXplorer (Kumar et al. 2020, PLoS Comput Biol 16(4):e1007737)**, fetched this session: an **FR** cluster ⇒ deletion candidate; an **RF** cluster ⇒ **tandem-duplication** candidate; FF/RR ⇒ inversion. DELLY: tandem duplications are *"paired-ends where the first and second read changed their relative order but kept the alignment strand"* (= RF / everted).

### Formula check
- μ ± c·σ span cutoff with default c = 3 matches BreakDancer verbatim. Inclusive bounds (discordant iff strictly outside) is a reasonable, sourced convention.
- DEL/INS/INV/CTX signature mapping matches Medvedev 2009 verbatim.

### Edge-case semantics
- Insertion > fragment size invisible to span (sourced). Below-support clusters dropped (BreakDancer -r). Empty → empty, null → throw. All sourced/defined.

### Findings / divergences (Stage A)
- **DEFECT (now corrected in the description):** the original description (SV_Detection.md §2.1/§4.1/§5.3), TestSpec (item 6, INV-2 area, S4) and Evidence all claimed **RF is concordant** ("FR/RF both proper, point inward"). This is **wrong** for the short-insert FR library the unit models: RF (outward-facing / everted) is the basic **tandem-duplication** signature and is discordant — confirmed independently by DELLY, LUMPY/Manta, SVXplorer, and even the unit's own cited cureffi/BWA source ("RF … that's a problem"). RF is "proper" only for opposite-orientation mate-pair libraries, which is not this contract. The description, TestSpec and Evidence were corrected this session and DELLY + SVXplorer added as sources 5–6.
- Stage A is otherwise correct. Verdict **PASS-WITH-NOTES** (one sourced correction).

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/StructuralVariantAnalyzer.cs`:
`FindDiscordantPairs` (L154), `IsConcordantOrientation` (L204), `ClassifySV` (L225), `DetectSVs` (L265), `ClusterDiscordantPairs` (L286).

### Formula realised correctly?
- Span bounds `μ ± c·σ`, strict `<`/`>` → exactly-at-bound is concordant. ✓ Matches BreakDancer.
- ClassifySV order inter-chr→INV→(now RF→DUP)→DEL→INS→Complex. ✓
- Min-support gate and cluster sweep correct. ✓

### Defect found and fixed (the RF orientation bug)
- **Before:** `IsConcordantOrientation` returned true for both FR **and RF**, so an RF same-chromosome pair within span bounds was never flagged discordant, and `ClassifySV` had no Duplication branch (an everted pair would have fallen to ComplexRearrangement). `SVType.Duplication` existed in the enum but was unreachable from the PEM path.
- **Fix:** `IsConcordantOrientation` → concordant iff **FR only** (`strand1=='+' && strand2=='-'`; the data model stores the upstream mate first). Added an RF branch to `ClassifySV` (`strand1=='-' && strand2=='+'` → `SVType.Duplication`). No other caller relied on the old behaviour (`StructuralVariantAnalyzerTests.cs` uses FR for "normal" and FF for "abnormal"; unaffected).

### Cross-verification table (recomputed vs the fixed code, μ=400, σ=50, c=3 ⇒ bounds [250,550])

| Input (chr1,str1,chr2,str2,span) | Source expectation | ClassifySV / discordant | Match |
|---|---|---|---|
| chr1,+,chr1,−,5000 | DEL (span>550) — Medvedev | Deletion | ✓ |
| chr1,+,chr1,−,100 | INS (span<250) — Medvedev | Insertion | ✓ |
| chr1,+,chr1,+,400 | INV (same strand) — Medvedev/cureffi | Inversion | ✓ |
| chr1,−,chr1,+,400 | **DUP (RF/everted)** — DELLY/SVXplorer | **Duplication; discordant** | ✓ (was wrong before fix) |
| chr1,+,chr2,−,400 | CTX (inter-chr) — Medvedev | Translocation | ✓ |
| chr1,+,chr2,+,400 | CTX precedence — A1 | Translocation | ✓ |
| chr1,+,chr1,−,400 | concordant FR | not discordant | ✓ |
| chr1,+,chr1,−,250 | bound inclusive | not discordant | ✓ |
| chr1,+,chr1,−,249 | outside → INS | discordant; Insertion | ✓ |
| chr1,+,chr1,−,550 | bound inclusive | not discordant | ✓ |
| 3× chr1,+,chr1,−,5000 (minSupport 2) | 1 DEL, support 3 — BreakDancer -r | 1 Deletion, SupportingReads=3 | ✓ |
| 1× DEL pair (minSupport 2) | no SV | empty | ✓ |

### Variant/delegate consistency
`DetectSVs`→`FindDiscordantPairs`→`ClusterDiscordantPairs`→`CreateSVFromCluster`→`ClassifySV` all share the same cutoff/orientation logic; the fix to `IsConcordantOrientation` and `ClassifySV` propagates consistently to `DetectSVs`. ✓

### Test-quality audit (HARD gate)
- **Green-washed test fixed:** S4 (`FindDiscordantPairs_RfOrientationWithinBounds_NotDiscordant`) asserted RF is *not* discordant — a defect that locked in the wrong behaviour. Rewritten to `…_IsDiscordantDuplication`: asserts the RF pair **is** flagged discordant and classifies as `Duplication`, citing DELLY/SVXplorer.
- **Coverage gap closed:** added **M9** `ClassifySV_RfEvertedOrientationSameChr_ReturnsDuplication` — the previously-unreachable Duplication branch now has an exact-value test.
- All other tests assert exact `SVType` values / exact counts / exact `SupportingReads`, cover all branches (TRANS, INV, DUP, DEL, INS, both boundaries, min-support below/meets, empty, null). No weakened assertions, no widened tolerances, no skips.
- **Honest green:** full unfiltered suite `dotnet test` = **Failed: 0, Passed: 6485** (pre-change baseline 6484 + 1 new test M9; S4 rewritten in place); the SV-DETECT-001 fixture itself = 15 passed. `dotnet build` 0 errors, no new warnings in changed files.

### Findings / defects (Stage B)
- One real algorithm defect (RF mis-classed as concordant; no Duplication branch) — **completely fixed** this session in code, tests, description, Evidence and TestSpec. Logged as FINDINGS_REGISTER A14.

## Verdict & follow-ups
- **Stage A:** PASS-WITH-NOTES (RF-concordance claim corrected against DELLY/SVXplorer/cureffi).
- **Stage B:** PASS-WITH-NOTES (RF/Duplication defect fixed; tests strengthened).
- **End-state:** ✅ CLEAN — defect fully fixed; algorithm fully functional; suite green (6485, Failed: 0).
- Out-of-scope by design (documented, not defects): split-read insertions > fragment, copy-paste/interspersed duplications and overlapping FR+RF complex clusters fold to `ComplexRearrangement`; linear cluster sweep rather than BreakDancer connection scoring.
