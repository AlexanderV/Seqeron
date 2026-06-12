# Validation Report: SPLICE-ACCEPTOR-001 — Acceptor Site Detection (3' splice site)

- **Validated:** 2026-06-12   **Area:** Splicing
- **Canonical method(s):** `SpliceSitePredictor.FindAcceptorSites(sequence, minScore, includeNonCanonical)`; internal `ScoreAcceptorSite`, `ScoreU12AcceptorSite`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — RNA splicing** (https://en.wikipedia.org/wiki/RNA_splicing): the 3' splice site
  "terminates the intron with an almost invariant AG sequence"; full consensus given as
  **Y-rich-N-C-A-G | G** (cut before the final exonic G). GU-AG rule confirmed: intron begins GU,
  ends AG. Branch point consensus **Y-U-R-A-C**, 20–50 nt upstream of the acceptor.
- **Wikipedia — Polypyrimidine tract** (https://en.wikipedia.org/wiki/Polypyrimidine_tract):
  pyrimidine-rich (especially U), **15–20 bp** long, located **5–40 bp** before the 3' end of the intron.
- TestSpec / Evidence cite **Shapiro & Senapathy (1987)** for the position frequencies
  (pos −2 A=100%, pos −1 G=100%, pos −3 C≈65–70%, pos 0 exonic G≈50%, upstream PPT C+U≈70–80%);
  **Burge et al. (1999)** for PPT-quality scoring; **Hall & Padgett (1994)/Jackson (1991)** for the
  U12 YCCAC variant. These match the cited literature.

### Feature/coordinate check
- Invariant **AG** at intron 3' end: present (recognition `upper[i]=='A' && upper[i+1]=='G'`).
- **Polypyrimidine tract**: modeled — pyrimidine count in window [pos−15, pos−3) contributes to score.
- **Branch point** (YNYURAC with branch A): modeled separately in `BranchPointPwm`/`FindBranchPoints`
  and consumed by intron prediction (not by the acceptor score itself, which is acceptable — the
  acceptor signal is AG + upstream PPT; the BPS is a distinct element scored separately).
- **AcceptorPwm consensus (Y)nNCAG|G**: positions −2 A=1.0, −1 G=1.0, −3 C=0.70, 0 G=0.50 — match
  Shapiro & Senapathy frequencies.
- **Coordinate convention**: 0-based scan; reported `Position = i+1` = index of the terminal G of AG
  (last intronic nucleotide). Self-consistent and documented (INV-5).
- **Donor vs acceptor not confused**: donor uses GT/GU, acceptor uses AG, U12 acceptor uses AC — all correct.

### Independent cross-check (numbers)
Hand-computed M1 worked example `UUUUUUUUUUUUUUUUCAGGG` (AG at index 17–18):
- PPT window [2,14) = 12 U → 12/12 → contribution 2.0
- PWM log2 terms over offsets {−15..0} relative to splice site (pos = i+2+offset): raw sum = 12.2136 over count 8
- normalized = (12.2136/9 + 2)/4 = **0.83927** → matches spec/test `0.8393` exactly.

### Findings / divergences (Stage A notes)
- **NOTE:** The acceptor *score* itself does not include a branch-point term; BPS is scored in a
  separate function. This is biologically sound (BPS is a distinct element) and the spec does not
  claim the acceptor score must embed it. Not a defect.
- **NOTE:** Normalization `(score/(count+1)+2)/4` is a documented heuristic mapping to [0,1], not a
  literature formula. Behavioural properties (range, monotonicity, ordering) are what the spec asserts,
  and those are verified by tests M5/M6/S5. Acceptable design decision.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs`
- `FindAcceptorSites` lines 209–250 (scan from i=15, AG/AC recognition, position i+1, motif, confidence)
- `ScoreAcceptorSite` lines 308–343 (PPT count + PWM at `position+2+offset`, normalization)
- `ScoreU12AcceptorSite` lines 386–430 (YCCAC consensus + PPT, /3.5 normalization)

### Formula realised correctly? (evidence)
Yes. The 2026-03-16 off-by-2 PWM-alignment fix is in place (`position + 2 + offset`), so PWM offset
−2 lands on the A and −1 on the G of the AG. Traced M1 by hand → 0.83927, matching code/test.

### Cross-verification table recomputed vs code
| Case | Expected | Recomputed/Observed | Match |
|------|----------|---------------------|-------|
| M1 canonical AG, pos | 18 | 18 (i+1, i=17) | ✅ |
| M1 score | 0.8393 | 0.83927 (hand) | ✅ |
| M2 no AG | empty | empty | ✅ |
| M3/M4 guards | empty | empty (len<20, null/empty) | ✅ |
| M5 strong>weak PPT | strong>0.7, weak<0.7 | passes | ✅ |
| M10 two AG | pos 18 & 22, first higher | passes | ✅ |
| S1 U12 YCCAC | pos 17, score 1.0 | passes | ✅ |
| S3 position | 16 (AG at idx15) | passes | ✅ |

### Variant/delegate consistency
DNA(T)↔RNA(U) equivalence and case-insensitivity verified (T→U + ToUpperInvariant). U12 path gated
by `includeNonCanonical`. Donor path uses GT/GU correctly — no acceptor/donor dinucleotide swap.

### Test quality audit
Canonical file `SpliceSitePredictor_AcceptorSite_Tests.cs` — 17 tests asserting exact positions,
exact scores (0.8393, 1.0), exact counts, type, range, ordering, and an independent PPT helper.
Not tautological. Filter `~AcceptorSite` = 23 passed (17 canonical + 6 donor-class name matches).

### Findings / defects
None. No off-by-one, no dinucleotide confusion, no missing spec-claimed feature.

## Verdict & follow-ups
- **Stage A:** PASS-WITH-NOTES (biology correct; acceptor score omits BPS term by design and uses a
  heuristic normalization — both acceptable and documented, not defects).
- **Stage B:** PASS — code faithfully realises the validated description; M1 reproduced to 4 dp.
- **State:** CLEAN — no defect found; no code changes required. Full suite 4486 passed / 0 failed.
