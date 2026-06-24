# Validation Report: SPLICE-ACCEPTOR-001 — Acceptor Site Detection (3' splice site)

- **Validated:** 2026-06-24   **Area:** Splicing
- **Canonical method(s):** `SpliceSitePredictor.FindAcceptorSites(sequence, minScore, includeNonCanonical)`; internal `ScoreAcceptorSite`, `ScoreU12AcceptorSite`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — RNA splicing** (https://en.wikipedia.org/wiki/RNA_splicing, fetched 2026-06-24):
  3' splice site (acceptor) consensus given as **`Y-rich-N-C-A-G-[cut]-G`** — the intron
  terminates with an "almost invariant AG". GT-AG rule confirmed (intron starts GU, ends AG).
  Branch point consensus **`Y-U-R-A-C`**, 20–50 nt upstream of the acceptor. PPT (C/U-rich) sits
  between the branch point and the AG.
- **Wikipedia — Polypyrimidine tract** (https://en.wikipedia.org/wiki/Polypyrimidine_tract,
  fetched 2026-06-24): pyrimidine-rich (especially uracil), **15–20 bp**, located **5–40 bp**
  before the 3' end of the intron.
- TestSpec/Evidence cite **Shapiro & Senapathy (1987)** for position frequencies
  (−2 A=100%, −1 G=100%, −3 C≈65–70%, +0 exonic G≈50%, upstream PPT C+U≈70–80%);
  **Burge et al. (1999)** for PPT-quality scoring; **Hall & Padgett (1994)/Jackson (1991)** for
  the U12 YCCAC variant. These are consistent with the fetched references.

### Feature/coordinate check
- Invariant **AG** at the intron 3' end: recognized via `upper[i]=='A' && upper[i+1]=='G'`.
- **Acceptor = AG, not donor GT/GU**: donor path uses GU/GC/AU; acceptor path uses AG; U12 acceptor
  uses AC. No dinucleotide confusion.
- **Polypyrimidine tract**: modeled — pyrimidine (C/U) count in window `[pos−15, pos−3)` contributes
  to the score (`pptScore/12*2`). Matches the documented C/U-rich PPT upstream of AG.
- **AcceptorPwm consensus `(Y)nNCAG|G`**: −2 A=1.00, −1 G=1.00, −3 C=0.70, +0 G=0.50, upstream
  pyrimidine-enriched (U=0.50/C=0.30) — matches Shapiro & Senapathy frequencies.
- **Branch point (YURAC)**: modeled separately in `BranchPointPwm`/`FindBranchPoints`, consumed by
  intron prediction — a distinct element, correctly not folded into the acceptor score itself.
- **Coordinate convention**: 0-based scan; reported `Position = i+1` = index of the terminal **G of
  AG** = last intronic nucleotide (the `...CAG|G` cut point). Self-consistent, documented (INV-5).

### Independent cross-check (numbers)
Hand-recomputed M1 `UUUUUUUUUUUUUUUUCAGGG` (AG at i=17; PWM pos = i+2+offset = 19+offset):
- PPT window `[2,14)` = 12×U → 12/12 → contribution **2.0**
- PWM log2 terms over the 8 offsets = **10.21362**; raw total = **12.21362**, count=8
- normalized = `(12.21362/9 + 2)/4` = **0.839267** → matches spec/test `0.8393` exactly.

Also traced by hand: M10 (two AG at pos 18 & 22), S1 (U12 YCCAC pos 17 → 3.5/3.5 = 1.0),
S3 (AG at i=15 → pos 16). All consistent with code and tests.

### Findings / divergences (Stage A notes)
- **NOTE:** The acceptor *score* omits a branch-point term; BPS is scored separately. Biologically
  sound (BPS is a distinct element) and not claimed by the spec. Not a defect.
- **NOTE:** Normalization `(score/(count+1)+2)/4` is a documented heuristic mapping to [0,1], not a
  literature formula. The spec asserts behavioural properties (range, monotonicity, ordering), which
  the tests (M5/M6/S5) verify. Acceptable design decision.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs`
- `FindAcceptorSites` lines 209–250 (guard len<20; scan from i=15; AG/AC recognition; `Position=i+1`;
  motif; confidence; U12 gated by `includeNonCanonical`)
- `ScoreAcceptorSite` lines 308–343 (PPT count + PWM at `position+2+offset`; normalization to [0,1])
- `ScoreU12AcceptorSite` lines 386–430 (YCCAC consensus + PPT; /3.5 normalization)

### Formula realised correctly? (evidence)
Yes. The PWM offset is applied as `position + 2 + offset`, so offset −2 lands on the A and −1 on the
G of AG, +0 on the first exonic base — correct splice-site-relative alignment. M1 reproduced by hand
to 6 dp (0.839267) matching the code/test value.

### Cross-verification table recomputed vs code
| Case | Expected | Recomputed/Observed | Match |
|------|----------|---------------------|-------|
| M1 canonical AG, pos | 18 | 18 (i+1, i=17) | ✅ |
| M1 score | 0.8393 | 0.839267 (hand) | ✅ |
| M2 no AG | empty | empty | ✅ |
| M3/M4 guards | empty | empty (null/empty, len<20) | ✅ |
| M5 strong>weak PPT | strong>0.7, weak<0.7 | passes | ✅ |
| M10 two AG | pos 18 & 22, first higher | passes | ✅ |
| S1 U12 YCCAC | pos 17, score 1.0 | 3.5/3.5=1.0 (hand) | ✅ |
| S3 position | 16 (AG at i=15) | passes | ✅ |
| C1 AG before scan start | empty | empty (scan starts i=15) | ✅ |

### Variant/delegate consistency
DNA(T)↔RNA(U) equivalence and case-insensitivity verified (`Replace('T','U')` + `ToUpperInvariant`).
U12 AC path gated by `includeNonCanonical`. Donor path uses GU/GC/AU — no acceptor/donor swap.

### Test quality audit
Canonical file `SpliceSitePredictor_AcceptorSite_Tests.cs` — 17 tests asserting exact positions,
exact scores (0.8393, 1.0), exact counts, type, [0,1] ranges, ordering, threshold subset, and an
independent PPT-fraction helper. Not tautological; deterministic. Covers all Stage-A edge cases.
Ran filtered: **17 passed / 0 failed**.

### Findings / defects
None. No off-by-one, no dinucleotide confusion, no missing spec-claimed feature.

## Verdict & follow-ups
- **Stage A:** PASS-WITH-NOTES — biology and consensus confirmed against current Wikipedia + cited
  primary literature; two documented design notes (BPS scored separately; heuristic normalization).
- **Stage B:** PASS — code faithfully realises the validated description; M1 and S1 reproduced by hand.
- **State:** CLEAN — no defect found; no code changes required.
</content>
</invoke>
