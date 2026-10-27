# Evidence Artifact: CRISPR-OFF-001

**Test Unit ID:** CRISPR-OFF-001
**Algorithm:** Off-target analysis (`CrisprDesigner.FindOffTargets`, `CalculateSpecificityScore`,
`CalculateMitHitScore`, `CalculateMitSpecificityScore` (both overloads), `CalculateSaCas9HitScore`,
`CalculateCfdScore`, `CalculateCfdSpecificityScore` (both overloads), `CfdAggregation`)
**Date Collected:** 2026-10-14 (campaign `docs/Validation/review-2026-09`, batch B08)

---

## Every URL opened (2026-10-14, re-verification pass)

All of these were fetched with `curl` through this container in this session; byte counts are the
downloaded files' own sizes. Nothing below was taken from memory, and no value in this document was
read off the C# implementation.

| URL | Status | Size | Used for |
|---|---|---|---|
| `https://raw.githubusercontent.com/maximilianh/crisporWebsite/master/crispor.py` | 200 | 375 382 B | every line citation below |
| `https://raw.githubusercontent.com/maximilianh/crisporWebsite/master/crisporEffScores.py` | 200 | — | checked; contains no off-target scoring, not cited |
| `https://raw.githubusercontent.com/maximilianh/crisporWebsite/master/bin/src/pairwise-library-screen/predictSingle.py` | 200 | 2 567 B | `SaCas9Scorer`, `mm_dict`, the position remap |
| `…/bin/src/pairwise-library-screen/models/rstan_mult_21_coeffs.txt` | 200 | 25 771 B | the 252 transcribed coefficients + the 21-nt oracles |
| `…/models/rstan_mult_20_coeffs.txt` | 200 | 25 053 B | the reference script's 20-nt assert |
| `…/models/rstan_mult_22_coeffs.txt` | 200 | 27 119 B | confirming tables ship for lengths 20-24 |
| `…/models/rstan_mult_23_coeffs.txt` | 200 | 28 586 B | same |
| `…/models/rstan_mult_24_coeffs.txt` | 200 | 29 738 B | same |
| `https://raw.githubusercontent.com/maximilianh/crisporWebsite/master/CFD_Scoring/mismatch_score.pkl` | 200 | 6 965 B | `mm_scores["rG:dA,1"] = 1.0` for the end-to-end CFD rows |
| `https://raw.githubusercontent.com/maximilianh/crisporWebsite/master/CFD_Scoring/pam_scores.pkl` | 200 | 344 B | `pam_scores["GG"]=1.0`, `["AG"]=0.25925925899999996`, `["TT"]=0.0` |
| `https://api.github.com/repos/maximilianh/crisporWebsite/contents/bin/src` | **403** | — | negative result: the directory cannot be listed from here |
| `…/bin/src/pairwise-library-screen/README` | **404** | — | negative result (see the SaCas9 section) |

Hsu et al. 2013 and Doench et al. 2016 are cited for the *origin* of the models only; their numeric
content is taken from the reference implementation's own code and shipped tables, listed above.

---

## Online Sources

### CRISPOR — reference implementation source (`crispor.py`, master)

**URL:** https://raw.githubusercontent.com/maximilianh/crisporWebsite/master/crispor.py
**Accessed:** 2026-10-14 (`curl`, 375 382 bytes; read verbatim, line numbers below are of that file)
**Authority rank:** 3 (reference implementation — Concordet & Haeussler, NAR 46:W242 (2018))

**Key extracted points (file:line):**

| Lines | What |
|---|---|
| **L889-891** | `pamIsCpf1(pam)` — the Cpf1/Cas12a PAM list (`TNN TTN TTTN TYCV TATV TTTV TTTR ATTN TTTA TCTA TCCA CCCA YTTV TTYN`). It does **not** contain CasX's `TTCN`, which has its own `pamIsCasX` at **L897-899**. |
| **L901-903** | `pamIsSaCas9(pam)` — `pam.split("-")[0] in ["NNGRRT", "NNNRRT"]`. |
| **L1872-1880** | `annotateOfftargets` docstring: "All 0MM matches in this range are ignored from scoring ('ontargets')". |
| **L1920-1929** | The on-target skip, and it is **narrower than the docstring suggests**: `if editDist==0 and repCount==0 and ((chrom==inChrom and start >= inStart and end <= inEnd) or (inChrom=='' and foundOneOntarget==False)): … continue`. Only 0MM hits inside the caller-supplied locus are dropped (or, in "auto-ontarget" mode with no locus, just the first 0MM hit met — `foundOneOntarget` then blocks the branch). A second perfect match elsewhere in the genome is scored as an off-target worth `calcHitScore = 100`; the aggregates now reproduce both modes (*Divergences* 10), and `repCount` (L1906-1910) is `maxOccurrences` (*Divergences* 11). |
| **L1931-1936** | `guideNoPam = guideSeq[:len(guideSeq)-len(pam)]` (PAM stripped before MIT scoring); `if len(otSeqNoPam) == 19: otSeqNoPam = "A" + otSeqNoPam` (the 19-nt pad is applied here too). |
| **L1939-1942** | `if isCpf1: mitScore = 0.0 ; cfdScore = 0.0` — "Cpf1 has no off-target scores yet". |
| **L1943-1945** | `elif isSaCas9: mitScore = calcSaHitScore(guideNoPam, otSeqNoPam) ; cfdScore = -1` — SaCas9 uses its own model and is excluded from CFD. |
| **L1947-1958** | `else: mitScore = calcHitScore(...)`; then `if pam == "NGG" and otSeq[-2:] != "GG": mitScore *= 0.2` (the alternative-PAM downweight). |
| **L1960-1962** | `mitOtScores.append(mitScore)`; `if cfdScore != -1: cfdScores.append(cfdScore)` — the CFD sentinel exclusion. |
| **L1986-1988** | `if pamIsCpf1(pam): guideScore = -1 ; guideCfdScore = -1`. |
| **L1994-1999** | `guideScore = calcMitGuideScore(sum(mitOtScores))`; `guideCfdScore = calcCfdGuideScore(sum(cfdScores))` under `doCfdFix`, else `calcMitGuideScore(sum(cfdScores))`. |
| **L2025-2040** | `calcSaHitScore` — delegates to `predictSingle.SaCas9Scorer` and returns `100.0 * saScorer.calcScore(...)`, "to be compatible with the MIT score, has to be in the range 0-100 for the MIT aggregate guide specificity score". |
| **L2045** | `hitScoreM` — the published Hsu weight vector, `[0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583]`, indexed 5'→3' (index 19 = PAM-proximal). |
| **L2047-2067** | `calcHitScore`'s own length handling: `len > 20` → last 20 nt of both ("S. aureus requires 21bp long guides. We fudge by using only last 20bp"); `len == 19` → left-pad both with `"A"` ("first pos has no weight anyways"); `len < 19` → `matrixStart = 20 - len`, `maxDist = len - 1` with the comment "for shorter guides, I'm not sure if this score makes sense anymore, we force things". |
| **L2098-2104** | `calcMitGuideScore(hitSum) = int(round(100/(100+hitSum) * 100))`. Seqeron returns the **unrounded** value. |
| **L2106-2109** | `calcCfdGuideScore(hitSum) = 100. * 100.0 / hitSum`, "suggested by Nicholas Parkinson". |
| **L2149-2185** | `calcCfdScore` — `re.search('[^ATCG]')` screen on both sequences → sentinel `-1`; `pam = off[-2:]`, `sg = off[:20]`; `if doCfdFix: cfd_score = cfd_score * 100.0` (**L2181-2182**). |
| **L340, L7045-7047** | `doCfdFix = False` by default; set to `True` only by the URL argument `fixCfd`. |
| **L4724-4729** | `showPamWarning` — the user-facing Cpf1 banner. See *Divergences* below: its prose contradicts the code. |

### CRISPOR — SaCas9 model (`predictSingle.py` + trained coefficients)

**URLs (both HTTP 200, downloaded 2026-10-14):**
- https://raw.githubusercontent.com/maximilianh/crisporWebsite/master/bin/src/pairwise-library-screen/predictSingle.py (2 567 bytes)
- https://raw.githubusercontent.com/maximilianh/crisporWebsite/master/bin/src/pairwise-library-screen/models/rstan_mult_21_coeffs.txt (25 771 bytes; also fetched `rstan_mult_20_coeffs.txt` 25 053 B and `rstan_mult_22_coeffs.txt` 27 119 B)

**Authority rank:** 3 (the trained model distributed with the reference tool).
**Primary source:** Tycko J, Wainberg M, Marinov GK, et al., "Pairwise library screen systematically
interrogates *Staphylococcus aureus* Cas9 specificity in human cells", *Nat Commun* **9**, 2962 (2018),
doi:10.1038/s41467-018-05391-2 — the paper CRISPOR's `calcSaHitScore` cites.

**Key extracted points:**

1. `SaCas9Scorer.__init__` reads `models/rstan_mult_<guideLen>_coeffs.txt` (tables shipped for guide
   lengths 20-24) and keys the `Median` column by
   `guideLenDict[guide_len - Position + 1][Mismatch]` — i.e. the shipped `Position` counts from the
   **3'/PAM-proximal** end and is remapped to a 1-based index from the **5'** end.
2. `calcScore` multiplies the coefficients of all mismatched positions, skipping `mm_pos == 1`
   ("Mismatches at the 5' position are not truly determined by the model"), and refuses `-`
   ("This model does not support gaps/bulges").
3. `mm_dict` maps the raw `(targetBase, guideBase)` pair to the coefficient id: `('A','C') → "dT:rC"`,
   `('A','T') → "dT:rU"`. So `d` is the **complement of the off-target DNA base** and `r` is the
   **guide base with T written as U** — the same convention as CFD, written in the opposite order.

**What was tried before concluding this model is obtainable** (the previous review had recorded it as
"a trained model not obtainable as published parameters"): the GitHub contents API
(`api.github.com/repos/maximilianh/crisporWebsite/contents/bin/src`) returns **403** from this
container, and `bin/src/pairwise-library-screen/README` **404**s, so the directory cannot be listed.
The files themselves are nonetheless served directly by `raw.githubusercontent.com` once the exact
paths from `predictSingle.py` (`models/rstan_mult_%s_coeffs.txt`) are used. They are plain TSV, so
`pandas` is not required to read them. **Conclusion: the model is fully obtainable and was
implemented** (finding F14) — campaign rule 2 applies, and the earlier "not obtainable" note is
withdrawn.

### Doench 2016 CFD matrices

The matrices themselves were transcribed and cross-checked in an earlier unit (recorded in
`docs/Evidence/CRISPR-GUIDE-001-Evidence.md` / the CFD block of `CrisprDesigner.cs`:
`mismatch_score.pkl` 240 entries and `pam_scores.pkl` 16 entries, decoded and cross-checked
element-by-element across `maximilianh/crisporWebsite` and `bm2-lab/iGWOS`, 240/240 + 16/16
identical) and are **unchanged by this unit**, which only adds the guide-level aggregation on top.

**URLs re-opened this session** (to source the per-hit CFD values used in the end-to-end rows below,
rather than reading them off the C# tables; both HTTP 200, `pickle.load`-decoded):

- https://raw.githubusercontent.com/maximilianh/crisporWebsite/master/CFD_Scoring/mismatch_score.pkl (6 965 bytes, 240 entries) — `mm_scores["rG:dA,1"] = 1.0`
- https://raw.githubusercontent.com/maximilianh/crisporWebsite/master/CFD_Scoring/pam_scores.pkl (344 bytes, 16 entries) — `pam_scores["GG"] = 1.0`, `pam_scores["AG"] = 0.25925925899999996`, `pam_scores["TT"] = 0.0`

Primary source: Doench, Fusi, Sullender, Hegde, et al., *Nat Biotechnol* **34**, 184-191 (2016),
PMID 26780180.

### Primary literature for the heuristic's motivation

- Hsu, Scott, Weinstein, et al., *Nat Biotechnol* **31**, 827-832 (2013), PMID 23873081 — PAM-proximal
  8-12 bp defines specificity; the position-weight model aggregated by crispr.mit.edu.
- Fu, Foden, Khayter, et al., *Nat Biotechnol* **31**, 822-826 (2013), PMID 23792628 — off-target
  frequency and truncated gRNAs.

These motivate the *direction* of `CalculateSpecificityScore`'s seed weighting only. **No source was
found for its numbers** (5 / 2, the 12-nt seed width, the linear sum, the subtraction from 100); it is
declared in the XML doc and in the `crispr_specificity_score` tool description as this library's own
heuristic, with the published aggregates offered alongside.

---

## Reference Cross-Checks (numbers)

All values below come from independent Python runs of verbatim transcriptions of the reference
functions, driven against the downloaded coefficient tables and pickles — never read off the C# code.
The harness was a throwaway script run under `python3 -I` from this session's scratchpad; it is not
checked in, because it is a literal copy-paste of the reference bodies listed above and reproducing it
needs nothing but the URLs in *Online Sources*:

- `calcHitScore`, `calcMitGuideScore` (with the `int(round(...))` dropped, see *Divergences* 6) and
  `calcCfdGuideScore` copied verbatim from `crispor.py` L2045-2109;
- `SaCas9Scorer.__init__` / `calcScore` copied from `predictSingle.py`, with `pandas.read_table`
  replaced by a plain tab-split read of the same `Position` / `Mismatch` / `Median` columns and the
  same `guideLenDict[guide_len - Position + 1]` remap;
- the CFD mismatch / PAM tables read straight out of the reference's own
  `CFD_Scoring/mismatch_score.pkl` (240 entries) and `pam_scores.pkl` (16 entries) with `pickle.load`.

The 12 × 21 SaCas9 coefficient table in `CrisprDesigner.cs` was additionally diffed **mechanically**
against `rstan_mult_21_coeffs.txt` (parse the `Median` column, apply the same position remap, extract
the C# arrays by regex, compare cell by cell): **252 values compared, 0 mismatching cells, maximum
absolute difference 0.0.**

### `calcHitScore` with the >20-nt fudge (21-nt pairs)

| guide | off-target | reference | note |
|---|---|---|---|
| `ACGTACGTACGTACGTACGTA` | `ACGTACGTACGTACGTACGTA` | `100.0` | identical |
| `ACGTACGTACGTACGTACGTA` | `TCGTACGTACGTACGTACGTA` | `100.0` | mismatch at index 0 falls outside the last-20 window |
| `ACGTACGTACGTACGTACGTA` | `ACGTACGTACGTACGTACGTC` | `41.7` | window index 19, W = 0.583 |
| `ACGTACGTACGTACGTACGTA` | `ACGTCCGTACGTACGTACGTC` | `6.389516129032259` | two mismatches; all three score terms active |

### `calcSaHitScore` (Tycko 2018, 21-nt table) — the reference script's own asserted values

| guide | off-target | `calcScore` | `calcSaHitScore` (×100) |
|---|---|---|---|
| `GCAACCACAAACCCACGAGGG` | `GCAACCACAAACCCACGAGGG` | `1.0` | `100.0` |
| `GCAACCACAAACCCACGAGGG` | `ACAAACACATACCCACAAGGA` | `0.041293256150042444` | `4.1293256150042446` |
| `GGGTGAGTGAGTGTGTGCGTG` | `GAGCGAGTGGGTGTGTGCGTG` | `0.13441983202750601` | `13.4419832027506` |

(The 20-nt table's own assert, `calcScore("TGTGGGTGAGTGTGTGCGTG","TGTGAGTGAGTGTGTGCGTG") =
0.314432049304524`, also reproduced exactly; that table is not transcribed into C# — see
*Divergences*.)

Coefficient-orientation guard (same mismatch type `dT:rC`, guide `C` vs off-target `A`):
index 1 → `54.6734106002513`, index 20 → `10.569373315278801`. A table transcribed without the
`guide_len - Position + 1` remap swaps these.

### `calcMitGuideScore` (unrounded)

| Σ hit scores | reference |
|---|---|
| `0` | `100.0` |
| `{2.5, 10.0}` → 12.5 | `88.88888888888889` |
| `60.5` | `62.30529595015576` |
| `{60.5, 41.7}` → 102.2 | `49.45598417408507` |

### CFD aggregation (over RAW per-hit CFD values in [0, 1])

| raw CFD scores | `MitStyle` | `ParkinsonNormalized` |
|---|---|---|
| `{0.5, 0.25}` | `99.25558312655087` | `133.33333333333334` |
| `{1.0}` | `99.00990099009901` | `100.0` |
| `{}` / `{0.0}` | `100.0` | **undefined** — division by zero (the reference raises `ZeroDivisionError`) |

`ParkinsonNormalized` is `100 / Σ raw`, **not** `100 × 100 / Σ raw`: the reference's `doCfdFix` flag
applies *two* changes together — `calcCfdScore` rescales every per-hit CFD by ×100 (L2181-2182) **and**
the aggregate becomes `100 × 100 / Σ` (L2106-2109) — so over raw inputs one factor of 100 cancels:
`100 × 100 / (100 × Σ raw) = 100 / Σ raw`.

### End-to-end genome-level cross-checks (locked in tests)

| Input | Result | Model |
|---|---|---|
| guide `GACGCATAAAGATGAGACGC`, genome `TACGCATAAAGATGAGACGCAGG` | MIT `50.0`, heuristic `98.0`, CFD(MitStyle) `99.00990099009901` | one PAM-distal off-target |
| guide `ACGTACGTACGTACGTACGTA`, genome `ACGTACGTACGTACGTACGTCAAGAAT`, SaCas9 | MIT `95.91285612347565` (per-hit `4.2613097364785695`) | Tycko 2018 |
| guide `ACGTACGTACGTACGTACGTA`, genome `ACGTCCGTACGTACGTACGTCAAGAAT`, SaCas9 | MIT `97.54818977171612` (per-hit `2.5134348817970387`) | Tycko 2018, two mismatches |
| guide `GACGCATAAAGATGAGACGC`, genome `CTTGCGTCTCATCTTTATGCGTA`, SpCas9-NAG | CFD(MitStyle) `99.74141115651946`, CFD(Parkinson) `385.71428610000004` | reverse-strand hit, PAM `AAG` → CFD `mm_scores["rG:dA,1"] (1.0) × pam_scores["AG"] (0.25925925899999996)`, both read from the reference pickles this session |

The last row doubles as the reverse-strand-PAM convention guard: if the PAM were reported
reverse-complemented (`CTT` → key `TT`, score `0.0`) the per-hit CFD would collapse to 0 and the
aggregate would be exactly `100`.

---

## Divergences recorded (not defects)

1. **SaCas9 coefficient tables for guide lengths 20, 22, 23, 24 are not transcribed.** The reference
   ships all five. Only the 21-nt table is in `CrisprDesigner.cs`, because 21 nt is the guide length of
   this library's only SaCas9 system (`CrisprSystemType.SaCas9`), and `FindOffTargets` rejects a guide
   that does not match its system's length — so no other length is reachable. Other lengths are
   **refused** by `CalculateSaCas9HitScore` rather than approximated. Adding a table is a pure
   data change if a 20/22/23/24-nt SaCas9 variant is ever added.
2. **The alternative-PAM ×0.2 MIT downweight (L1947-1958) is implemented, over a PAM-annotated hit
   set (B08 WP8 / F50).** CRISPOR runs one wider BWA scan whose NGG query also returns non-GG
   (NAG/NGA) sites and downweights those: `if pam == "NGG" and otSeq[-2:] != "GG": mitScore *= 0.2`
   ("alternative PAMs represent only ~10% of all cleaveage events... so that these off-targets are not
   ranked among the top but still appear in the list somewhat"). A *single-system* Seqeron scan cannot
   trigger it — an `SpCas9` scan yields only `NGG` hits (`otSeq[-2:] == "GG"`, condition false) and an
   `SpCas9_NAG` scan has `pam != "NGG"` (condition also false), which is why the previous round
   recorded it as unreachable — so the rule is offered on the two entry points that CAN see several
   PAMs: `CalculateMitSpecificityScore(string, IEnumerable<PamAnnotatedOffTarget>, CrisprSystemType)`
   (any PAM-annotated hit set, from the new public `FindOffTargetsWithPam`) and
   `CalculateMitSpecificityScoreWithAlternativePams`, which scans the main PAM plus every alternative
   PAM of the reference's `offtargetPams` table (**L243-253**) that this library offers as a system.
   Three properties of the reference rule are reproduced exactly and locked by test: it is gated on the
   *guide's* PAM being literally `NGG` (so an `SpCas9-NAG` guide is never downweighted); it tests the
   hit's last two bases of `protospacer + PAM`, not membership of the alternative table (any non-`*GG`
   PAM is downweighted); and **only** the MIT score is downweighted — the reference computes the hit's
   CFD score afterwards with no 0.2 factor, so `CalculateCfdSpecificityScore` is untouched.
   B08 round 3 G1 (F80) gives the CFD aggregate the same wider search: `annotateOfftargets` computes
   `cfdScore = calcCfdScore(guideSeq, otSeq)` for **every** hit of the NAG/NGA search (L1905-1990, no 0.2
   factor), after the pre-filter of L4104-4122 that tests the **MIT** `calcHitScore` against
   `ALTPAMMINSCORE`, so a guide whose only off-target sits on `NAG` is not perfectly specific.
   `CalculateCfdSpecificityScoreWithAlternativePams` shares the hit collection (alternative-PAM list,
   scannable systems, main/alternative dedup, MIT-score cutoff) with the MIT method and sums the raw CFD
   values. Oracle (crispor.py L2043-2185 exec'd verbatim in Python, pickled CFD tables, `pam_scores["AG"]`
   = 0.25925925899999996): NAG-only genome -> CFD 0.18518518492592592 -> **99.81515711670932** (was 100);
   mixed NGG/NAG/NGA genome -> 99.24515917350747 (cutoff 1.0), 99.20727233740708 (cutoff 0).
   B08 round 2 E1 (F66/F67) closed two of the former sub-divergences. (a) `SpCas9_NGA` (`NGA`, 20 nt)
   and `SaCas9_NNGRRN` (`NNGRRN`, 21 nt) were appended to `CrisprSystemType`, so the alternative-PAM
   aggregate now scans `NGG → NAG, NGA` and `NNGRRT → NNGRRN` (before, the `NGA`/`NNGRRN` scans were
   silently dropped and an SaCas9 guide's alt-PAM aggregate equalled the plain one); a site matching both
   `NNGRRT` and `NNGRRN` is counted once (CRISPOR runs one search), and SaCas9 `NNGRRN` hits are NOT
   downweighted (the 0.2 rule is gated on `pam=="NGG"`). Oracle: guide `GACGCATAAAGATGAGACGC`, one
   1-mismatch hit `GACGCTTAAAGATGAGACGC` + `CGA` -> `calcHitScore` 60.5 (crispor.py code re-run) ->
   `100/(100+0.2x60.5)x100 = 89.20606601248885` (was 100.0). (b) The pipeline-level cutoff
   (`ALTPAMMINSCORE = 1.0`, **L318**, applied at **L4108-4122**, `MINSCORE = 0.0` for main PAM hits;
   `--minAltPamScore`, L8439-8441) is the optional `minAlternativePamHitScore` of
   `CalculateMitSpecificityScoreWithAlternativePams`; the DEFAULT is `1.0`, the
   reference's own (B08 round 3 G6; it was `0` = no filtering before — pass `0` for that). The reference tests the
   alternative list first, so for SaCas9 the cutoff also reaches `NNGRRT` main-PAM hits; reproduced.
   The original single-system overload `CalculateMitSpecificityScore(string, DnaSequence, int,
   CrisprSystemType)` keeps its behaviour unchanged (locked by test).
3. **PAM-before-target systems are refused outright**, including CasX. For Cpf1/Cas12a this follows the
   reference's code exactly (`mitScore = cfdScore = 0.0`, `guideScore = guideCfdScore = -1`).
   CasX (`TTCN`) is **not** in the reference's `pamIsCpf1` list (L889-891) — it has its own
   `pamIsCasX` (L897-899) — so CRISPOR routes a CasX guide into the `else` branch and scores it with
   `calcHitScore`. **That path is demonstrably mis-oriented**, traced through the code this session:
   `pamIsCasX` sets `pamIsFirst = True` (**L483-487**); `concatGuideAndPam` therefore builds the scored
   string as `pamPlusSeq + pamSeq + guideSeq` — PAM at the **5'** end (**L5543-5548**); but
   `annotateOfftargets` strips the PAM with `guideSeq[:len(guideSeq)-len(pam)]` (**L1932-1933**), which
   removes bases from the **3'** end. What reaches `calcHitScore` thus still carries the 4-nt `TTCN`
   PAM at its 5' end and has lost the 3'-most 4 protospacer bases, while `hitScoreM` is indexed
   PAM-proximal-last — so every one of the 20 position weights lands on the wrong base. Reproducing
   that would be reproducing a reference bug, so this library refuses instead. Returning the reference's own sentinel `-1` was rejected because it would break
   the documented `[0, 100]` range; returning `100` would be unsafe (a guide whose off-targets were
   never scored would look perfectly specific — exactly defect F9).
4. **CRISPOR's own Cpf1 prose contradicts its code.** `showPamWarning` (L4724-4729) says "there is no
   off-target ranking algorithm available in the literature, to our knowledge. **We use Hsu and CFD
   scores below for off-target ranking**, but they were developed for spCas9." The scoring path
   (L1939-1942, L1986-1988) zeroes both scores for Cpf1 PAMs. The **code** is taken as the reference
   behaviour; the banner is recorded here so the discrepancy is not re-litigated.
5. **The `calcHitScore` branch for guides shorter than 19 nt is not implemented.** The reference's own
   comment is "for shorter guides, I'm not sure if this score makes sense anymore, we force things",
   and the branch is unreachable here: every `CrisprSystem.GuideLength` is ≥ 20 and `FindOffTargets`
   rejects a mismatched guide. The 19-nt `A`-pad **is** implemented, because `annotateOfftargets`
   applies it as well (L1934-1936).
6. **The aggregate is returned unrounded**; `calcMitGuideScore` rounds to an integer
   (`int(round(score*100))`). Callers may round as they wish.
7. **Ambiguity codes cannot reach the off-target scan.** `DnaSequence` validates A/C/G/T at
   construction (`DnaSequence.ValidateSequence` → `SequenceExtensions.IndexOfInvalidDna`), so the
   reference's `[^ATCG] → cfdScore = -1` exclusion inside `CalculateCfdSpecificityScore` is a
   **defensive mirror**, not a live parity fix. It is kept so the reference's exclusion rule is
   encoded rather than implied, and a test pins the `DnaSequence` rejection that makes it unreachable.
7a. **Raw-string guide alphabet (B08 round 3 G2).** A guide is a raw string, so `DnaSequence` validation does not cover it. The off-target scan, the aggregates built on it, `CalculateMitHitScore` and `GetGrafMotif` normalise it as F70 does for `EvaluateGuideRna`: ASCII fold, `U` read as `T`, every other character an `ArgumentException`. Without it a `N`x20 guide mismatched every site (specificity 100) and an RNA-spelled guide reported its own on-target as a 5-mismatch off-target. The reference `calcHitScore` itself compares characters without validating; rejecting is this library's choice (same class as F9/F70). The CFD and SaCas9 per-hit scorers keep their A/C/G/T check without the `U` reading (the CFD reference screens `[^ATCG]`).
8. **Bulges / indels are out of scope** for every model here; the reference's SaCas9 model refuses them
   explicitly and CFD is defined only for a substitution-only 20-nt:20-nt comparison.
9. **`CalculateCfdSpecificityScore` refuses SaCas9; the reference returns a vacuous 100 for it.**
   The per-hit side is identical — the reference sets `cfdScore = -1` for every SaCas9 hit
   (L1943-1945). But its guide level then does *not* emit a sentinel: those hits are all dropped by
   `if cfdScore != -1` (L1961-1962), so `guideCfdScore = calcMitGuideScore(sum([])) = 100`. Only Cpf1
   gets `guideCfdScore = -1` (L1986-1988). Reporting a perfect 100 for a guide whose off-targets were
   never CFD-scored is exactly the silent-perfect-score failure mode this unit fixed elsewhere
   (defect F9), so this library throws instead. CFD is in any case defined only for a 20-nt
   protospacer and the reference provides no length fallback for it.
10. **The on-target is ONE hit (or the hits inside an interval), not every exact match (B08 round 3 G7 / F81, corrected
   in round 4 H5).** The reference's skip (L1917-1928) is scoped to the caller-supplied on-target locus - or, with
   no locus ("auto-ontarget": `inChrom=='' and foundOneOntarget==False`), to the single first 0MM hit - so a
   *second* perfect match elsewhere in the genome is scored as a full `calcHitScore = 100` off-target.
   `FindOffTargets` (contract: the 1..max-mismatch enumeration) still drops every exact match. The four
   CRISPOR-faithful genome aggregates take an optional `onTargetStart`/`onTargetEnd` (0-based, end-exclusive,
   forward strand; footprint = protospacer + PAM; strand not compared) and skip exactly
   `editDist == 0 and start >= inStart and end <= inEnd`. **Without an interval they now reproduce
   auto-ontarget** (round 3 had declined it - "first follows the aligner's output order, so no deterministic
   equivalent exists" - which the code contradicted: the scan is fully deterministic, and exact matches of one
   PAM class score alike, so the choice cannot change the value; only the alternative-PAM aggregates could
   differ, and there the first exact match of the main-PAM scan is the on-target while exact alternative-PAM
   matches are kept). Oracle: CRISPOR's own `annotateOfftargets`, **executed** (the genuine function body
   loaded by AST from `crispor.py` sha256 `16c3a54b...`, with `makeAlnStr`/`isInPar` stubbed and the CFD matrices
   read from their JSON conversion): N perfect sites, no interval -> guideScore 100, 50, 20, 2, 0 and
   guideCfdScore 100, 99, 96, 67, 17 for N = 1, 2, 5, 50, 500 (integer-rounded; unrounded 100/N and
   100/(100+N-1)*100); duplicated site + one 1-mismatch hit -> 38 / 98 (38.38771593090211 / 98.3146067418492);
   exact NGG + NAG + NGA -> 71 / 100 (71.42857142857143 / 99.67237321836325). With an interval `[0,23)` -> MIT
   `50.0`, CFD `99.00990099009901`; `[0,50)` -> 100; a footprint only partly inside (`[1,23)`) is not skipped.
11. **Repeat guides - `repCount` (round 4 H4).** The reference suppresses a guide wholesale when any hit has
   `totalAlnCount > MAXOCC` (60000, L269) - or `totalAlnCount > 1 and isRep` (L1909-1910): `guideScore = 0`,
   `guideCfdScore = 0`, an empty off-target list, and no on-target skip (L1924 requires `repCount==0`). The old
   code scored a guide with 500 perfect sites 100 where the reference gives 0 (and the same 100 for 2, 5, 50).
   `totalAlnCount` is BWA's `X0 + X1` for the guide alone, "not looking at the PAM" (L4081): **verified against
   BWA 0.7.17** with the reference's own command line (`bwa aln -o 0 -m 2000000 -n k -k k -N -l 20`, `bwa samse
   -n 60000`) on 40 random genomes with planted 0-5-mismatch copies on both strands, every fifth guide a
   palindrome: X0+X1 equals the number of guide-length windows within k mismatches on either strand in 40/40.
   `isRep` is BWA omitting the alignments from the XA tag, i.e. the same condition (`samse -n` is `MAXOCC`), so the
   rule is one comparison. The four aggregates take `maxOccurrences` (default 60000) and return 0 above it; the
   enumerations have no score to zero and still list their hits. The alignment count is PAM-agnostic, which is why it
   is a separate window scan and not the PAM-adjacent hit count.

---

## Unsourced element declared (not fixed, cannot be)

`CalculateSpecificityScore` — `max(0, 100 − Σ(5 per seed mismatch + 2 per distal mismatch))`. No
source was found for the weights, the 12-nt seed width, the linear sum or the subtraction from 100.
Per the campaign's rule for unsourced convenience metrics it is **kept and declared honestly** (XML
doc, algorithm doc, `crispr_specificity_score` MCP description, this file) and the two published
alternatives in the same class are offered in each of those places and exposed over MCP as
`crispr_mit_specificity_score` / `crispr_cfd_specificity_score`. A test
(`CalculateSpecificityScore_HeuristicAndMitAggregate_AreDistinct`) pins that the heuristic and the
published aggregate are visibly different numbers on the same input, so one cannot be mistaken for
the other.

---

## D19: DNA/RNA-bulge enumeration oracle (2026-10-25)

**Primary sources.** Bae, Park & Kim, *Bioinformatics* 30:1473–1475 (2014),
https://doi.org/10.1093/bioinformatics/btu048, describes Cas-OFFinder's search.
The reference implementation's 3.0 README and release supply the DNA/RNA-bulge
input and output semantics:
https://github.com/snugel/cas-offinder/blob/develop/README.md and
https://github.com/snugel/cas-offinder/releases/tag/3.0.0b3.

**Executed reference.** Official Windows x86-64 `3.0.0b3` release archive;
`cas-offinder.exe` prints `Cas-OFFinder 3.0.0`. It ran on the available
Intel Arc OpenCL GPU (`G` mode); no CPU OpenCL device was present. The exact
archive SHA-256 is `29C9AE74197C1EC11A609D030A647BF252785A03E2D56BEC3C33AA09DE5A09BC`. The exact
PowerShell commands were:

```powershell
$dir = (Resolve-Path '.git/d19-oracle').Path
& "$dir/cas-offinder/build/cas-offinder.exe" "$dir/input.txt" G "$dir/output.txt"
& "$dir/cas-offinder/build/cas-offinder.exe" "$dir/input-mismatch.txt" G "$dir/output-mismatch.txt"
& "$dir/cas-offinder/build/cas-offinder.exe" "$dir/input-two.txt" G "$dir/output-two.txt"
```

The first input file contained these three lines (the first line is the absolute
path `D:\Prototype\Bio\Seqeron\.git\d19-oracle\toy.fa` on the oracle machine):

```text
D:\Prototype\Bio\Seqeron\.git\d19-oracle\toy.fa
NNNNNNNNNNNNNNNNNNNNNGG 2 2
GGCCGACCTGTCGCTGACGCTGG 0
```

Full `toy.fa`:

```fasta
>toy
TTTTTTTTTTGGCCGACCATGTCGCTGACGCTGGTTTTTTTTTTTTTTTGGCCGACCGTCGCTGACGCTGGTTTTTTTTTTTTTTTCCAGCGTCAGCGACATGGTCGGCCTTTTTTTTTTTTTTTCCAGCGTCAGCGACGGTCGGCCTTTTTTTTTT
```

Full raw `output.txt` (tabs between fields):

```text
##Generated by Cas-OFFinder 3.0.0
#Id	Bulge Type	crRNA	DNA	Chromosome	Location	Direction	Mismatches	Bulge Size
0	DNA	GGCCGACC-TGTCGCTGACGCTGG	GGCCGACCATGTCGCTGACGCTGG	toy	10	+	0	1
0	DNA	GGCCGACC-TGTCGCTGACGCTGG	GGCCGACCATGTCGCTGACGCTGG	toy	87	-	0	1
0	RNA	GGCCGACCTGTCGCTGACGCTGG	GGCCGACC-GTCGCTGACGCTGG	toy	49	+	0	1
0	RNA	GGCCGACCTGTCGCTGACGCTGG	GGCCGACC-GTCGCTGACGCTGG	toy	128	-	0	1
```

The second input file used `toy-mismatch.fa`, DNA bulge 1/RNA bulge 0,
and one mismatch:

```text
D:\Prototype\Bio\Seqeron\.git\d19-oracle\toy-mismatch.fa
NNNNNNNNNNNNNNNNNNNNNGG 1 0
GGCCGACCTGTCGCTGACGCTGG 1
```

```fasta
>toy_mm
TTTTTTTTTTGGCCGACCATGTCGCAGACGCTGGTTTTTTTTTT
```

Full raw `output-mismatch.txt`:

```text
##Generated by Cas-OFFinder 3.0.0
#Id	Bulge Type	crRNA	DNA	Chromosome	Location	Direction	Mismatches	Bulge Size
0	DNA	GGCCGACC-TGTCGCTGACGCTGG	GGCCGACCATGTCGCaGACGCTGG	toy_mm	10	+	1	1
```

The third input tests the allowed two-base ceiling:

```text
D:\Prototype\Bio\Seqeron\.git\d19-oracle\toy-two.fa
NNNNNNNNNNNNNNNNNNNNNGG 2 2
GGCCGACCTGTCGCTGACGCTGG 0
```

```fasta
>toy_two
TTTTTTTTTTGGCCGACCATTGTCGCTGACGCTGGTTTTTTTTTTTTTTTGGCCGACCTCGCTGACGCTGGTTTTTTTTTT
```

Full raw `output-two.txt`:

```text
##Generated by Cas-OFFinder 3.0.0
#Id	Bulge Type	crRNA	DNA	Chromosome	Location	Direction	Mismatches	Bulge Size
0	DNA	GGCCGACC--TGTCGCTGACGCTGG	GGCCGACCATTGTCGCTGACGCTGG	toy_two	10	+	0	2
0	RNA	GGCCGACCTGTCGCTGACGCTGG	GGCCGACC--TCGCTGACGCTGG	toy_two	50	+	0	2
0	RNA	GGCCGACCTGTCGCTGACGCTGG	GGCCGACCT--CGCTGACGCTGG	toy_two	50	+	0	2
```

**Comparison contract.** The reference reports `Location` 10/87/49/128 in the
four rows above. The physical forward-strand starts of their complete
protospacer-plus-PAM footprints in the pasted FASTA are 10/86/49/125.
The reverse-strand reference locations are therefore shifted by +1 (DNA)
and +3 (RNA) relative to those starts in this release. The
[reference output code](https://github.com/snugel/cas-offinder/blob/3.0.0b3/cas-offinder.cpp#L395-L446)
prints `loci - m_chrpos[idx] + offset` after padding the search pattern for
the maximum DNA-bulge size; `offset = 0` in the reverse-PAM branches. We
retain the physically checked FASTA coordinates instead of reproducing the
reference's padded scan offset. Seqeron's
`BulgedOffTargetSite.Position` keeps the existing `OffTargetSite.Position`
convention: the forward-strand PAM coordinate. Its `TargetStart` is the
forward-strand protospacer start. Thus the same four oracle sites map to
`(Position, TargetStart, IsForwardStrand)` = `(31,10,true)`, `(86,89,false)`,
`(68,49,true)`, `(125,128,false)`. The aligned strings in the API omit the
three PAM bases; the gap offset is 8 in the first five oracle rows. The reference's
lowercase `a` marks the one substitution in the last row; the API retains
uppercase sequence plus `Mismatches = 1`.

**D19 decision record.** I used the official Cas-OFFinder 3.0.0b3
executable rather than the older separate `cas-offinder-bulge` wrapper,
because this release has native DNA/RNA-bulge support; this Windows host
exposed an OpenCL GPU but no CPU device, so the exact invocation used `G`.
This is a separate alignment enumerator and record,
with no changes to the existing substitution-only scan or scoring chain.
It reports each admissible gap offset and includes zero-mismatch bulged sites,
but drops zero-mismatch ungapped sites because the API has no on-target locus.
The 0..2 bound is this public API's explicit practical limit, not a claim
that Cas-OFFinder itself rejects larger bulges. The reverse-coordinate
difference just described is a deliberate departure from this reference
release's reported `Location`, while the sequence alignment and mismatch
fields match exactly. The two-base RNA target has two equivalent gap offsets;
both are reported, as the reference does. No MCP wrapper was added.

## D19 follow-up: admissible gap offsets, PAM on the hit, on-target interval (B08 round 4, H2/H3/H26)

**What the D19 oracle could not show.** The eight real-binary rows above all sit on sites that exist *only* as bulged
alignments, so they cannot expose an alignment the reference does **not** generate. The first scan reported two such
classes: a DNA bulge at the gap offset that rebuilds the ungapped query, and an RNA bulge whose removed run ends at the
last guide base. A perfect site came back as `Dna/1@0`, `Dna/2@0`, `Rna/1@0`, `Rna/2@0` and never once as itself.

**Source (Cas-OFFinder 3.0.0b3 `cas-offinder.cpp`, tag `3.0.0b3`, read in full for `parseInput` and the output block).**
With `s0` the guide padded with the PAM's `N` columns, `first`/`last` = `find_first_not_of('N')`/`find_last_not_of('N')`
and `dnab`/`rnab` the bulge ceilings:

| Branch | Loop | Gap offsets (string coordinates) |
|---|---|---|
| DNA, PAM after | `j = first .. last+1`, `compare = N^(dnab-i) s0[:j] N^i s0[j:]` | `j = 0` equals `tmp = N^dnab s0` -> type `X`, size 0 |
| DNA, PAM before | same with `compare = s0[:j] N^i s0[j:] N^(dnab-i)` | `j = last+1` equals `tmp = s0 N^dnab` -> type `X`, size 0 |
| RNA | `j = first .. last-i`, removes `s0[j:j+i]` | the last guide base is never removable |

The older `cas-offinder-bulge` v1.2b wrapper is stricter still (`range(1, len(target))` for DNA, `range(1, len(target)-size)`
for RNA). Seqeron follows 3.0.0b3, as its documentation always said.

**Oracle.** `scripts/cas-offinder/cof3_bulge_oracle.py` is a line-by-line port of those loops and of the output
formatting (both PAM geometries), with the OpenCL matcher replaced by its documented semantics. **Validated against the
real binary:** on the toy genomes of the D19 section it reproduces all 8 recorded rows exactly (4 exact, 1 with a
substitution, 3 two-base rows including the two equivalent RNA offsets). It is the oracle for the sweep below; it is not
a new run of the binary.

**Differential (this session, 2026-10-27).** 8000 random cases (800 x 10 seeds; four systems SpCas9 NGG/20,
SaCas9 NNGRRT/21, Cas12a TTTV/23, CasX TTCN/20; 0-4 substitutions; DNA and RNA ceilings 0-2; planted perfect, substituted,
DNA-bulged and RNA-bulged sites on both strands), exact exclusion of the 0-mismatch ungapped row applied to the oracle:

| | oracle rows | Seqeron only | oracle only |
|---|---:|---:|---:|
| before the fix (seed 1, 600 cases) | 9632 | **1317** (SpCas9/SaCas9 DNA@0: 406; Cas12a/CasX DNA@L: 399; RNA at windowLength: 512) | 0 |
| after the fix (10 seeds x 800 cases) | 129 221 | **0** | **0** |

(The earlier audit's randomised sweep reported 70 extra rows against 3.0.0b3 over 8 SpCas9 genomes; the larger and
four-system sweep finds the same two classes and, with Cas12a/CasX, the mirrored DNA collapse at the other end.)

**Sequence ends - a recorded difference, not a defect.** The tool pads its pattern with `dnab` wildcard columns, so a PAM
whose padded window crosses a sequence end is never matched: with sites planted at the start (seed 11, 600 cases) the oracle
lists 12 116 rows and Seqeron reports those plus 578 more, none missing. Those alignments lie entirely inside the sequence,
so they are reported.

**H26.** The reference prints the PAM as the last columns of the crRNA and DNA fields; `BulgedOffTargetSite` dropped it, so
no hit could reach a PAM-dependent score without re-deriving the PAM from coordinates. It now carries `Pam` (read on the
protospacer strand, like `PamAnnotatedOffTarget.Pam`). The optional on-target interval has CRISPOR's semantics
(`annotateOfftargets`, crispor.py L1920-1929, the same footprint rule as the aggregates): an exact ungapped match is
dropped only when its protospacer+PAM footprint lies inside it; with none given, every exact match is dropped as before (this method's contract; the specificity aggregates follow auto-ontarget since H5).
