# Evidence Artifact: CRISPR-GUIDE-001

**Test Unit ID:** CRISPR-GUIDE-001
**Algorithm:** Guide RNA design & on-target evaluation (`CrisprDesigner.DesignGuideRnas`, `EvaluateGuideRna`,
`GetCutSite`, `CountSelfComplementaryStems`, `GetGrafMotif`, `CalculateOnTargetDoench2014`,
`CalculateOnTargetRuleSet2` / `AzimuthRuleSet2`, `GuideRnaParameters` / `GuideRnaCandidate`)
**Date Collected:** 2026-10-01/02 (campaign `docs/Validation/review-2026-09`, batch B08)

---

## Online Sources

### CRISPOR — reference implementation source (`crispor.py`, `crisporEffScores.py`, master)

**URL:** https://raw.githubusercontent.com/maximilianh/crisporWebsite/master/crispor.py
**Accessed:** 2026-10-01 (`curl`; read verbatim)
**Authority rank:** 3 (reference implementation — Concordet & Haeussler, NAR 2018)

**Key extracted points:**

1. **Cas9 cleavage position.** The expected cleavage position is "-3bp 5' of the PAM site"; the cut feature is
   built as `startFt = start - 3` on the plus strand and as `ftSeq + "---"` on the minus strand, i.e. the three
   protospacer bases adjacent to the PAM. → forward `Position - 3`, reverse `Position + pamLen + 2`
   (forward-strand coordinates, as produced by `FindPamSites`).
2. **Cas12a/Cpf1 cleavage position.** Staggered cut "after the 18th base on the non-targeted strand which has
   the TTTV PAM motif" (Zetsche et al. 2015, Cell 163:759, Fig. 3). → forward `Position + pamLen + 18`,
   reverse `Position - 19`.
3. **Rule Set 2 context window.** `crisporEffScores.calcAllScores` builds the Doench-2016 input with
   `trimSeqs(seqs, -24, 6)` around the PAM start = 4 nt 5' flank + 20-nt protospacer + 3-nt NGG PAM + 3 nt
   3' flank (30 nt, protospacer strand); windows containing N are not scored.
4. **Ranking.** `mergeGuideInfo` sorts the guide table `reverse=True` by the selected score column
   (descending by score, best first).
5. **Graf 2019 motifs** (`crisporEffScores.getGrafType`): TT-motif when the guide ends in `TTC`/`TTT`, or its
   last four bases are only T/C with ≥2 T, or contain `TT` with ≥3 T or ≥1 C; GCC-motif when it ends in
   `[AGT]GCC` or `GCCT`. Applied to NGG systems only and reported as a warning ("Inefficient"), with no change
   to any score. Primary source: Graf et al. 2019, Cell Reports 26:1098–1103.

### CHOPCHOP — reference implementation source (`chopchop.py`)

**URL:** https://raw.githubusercontent.com/Carl-labhub/chopchop-clone/master/chopchop.py
(mirror of the Bitbucket `valenlab/chopchop` repository; 144 791 bytes)
**Accessed:** 2026-10-02 (`curl`)
**Sources tried first and unreachable:** `https://bitbucket.org/valenlab/chopchop/raw/HEAD/chopchop.py` (404),
`https://raw.githubusercontent.com/valenlab/chopchop/master/chopchop.py` (404),
`https://raw.githubusercontent.com/MicrobialDarkMatter/chopchop/master/chopchop.py` (404).
**Authority rank:** 3 (reference implementation — Labun et al., Nucleic Acids Res 44:W272 (2016); 47:W171 (2019))

**Key extracted points (line numbers in the fetched file):**

1. `STEM_LEN = 4` (L213).
2. `selfComp(fwd, backbone)` (L1882-1893) and `Guide.calcSelfComplementarity` (L378-394, L588-604):
   `rvs = revComp(fwd)`, `L = len(fwd) - STEM_LEN - 1`, and for every window `fwd[i:i+STEM_LEN]` with
   `gccontent(...) >= 0.5`, `folding += 1` when the window occurs in `rvs[0:(L-i)]` **or** inside any backbone
   region. The PAM is excluded: `fwd = self.guideSeq[len(PAM):]  # Do not include PAM motif in folding calculations`.
3. `gccontent` (L963-968) is the plain G/C fraction of the window (case-insensitive).
4. **Backbone orientation** (L2857-2858): `tmp = args.backbone.strip().split(",")`;
   `args.backbone = [str(Seq(el).reverse_complement()) for el in tmp]` — the caller supplies the scaffold on the
   guide's strand (CLI help, L2794: "comma-separated list, same strand as guide") and the tool
   reverse-complements it internally.
5. **Filter default** (L2828): `-filterSelfCompMax` default `-1`, "no filter"; the filter is applied as
   `if add and filterSelfCompMax != -1: … if folding > filterSelfCompMax` (L1814-1820, L1854-1860).
6. **GC window** `GC_LOW = 40` / `GC_HIGH = 70` — the same acceptance window as `GuideRnaParameters.Default`.

---

## Reference Cross-Checks (numbers)

### Self-complementarity stem count (CHOPCHOP `selfComp`)

The reference algorithm was re-implemented verbatim in Python from the source above and driven against the C#
`CrisprDesigner.CountSelfComplementaryStems`:

| Guide | Backbone argument | Reference | Seqeron |
|---|---|---|---|
| `ACGGACTAGCCTACGTACGT` | — (guide-internal only) | 0 | 0 |
| `ACGGACTAGCCTACGTACGT` | `AGGCTAGTCCGT` (standard scaffold, guide strand) | 8 (`ACGG CGGA GGAC GACT CTAG TAGC AGCC GCCT`) | 8 |
| `ACGGACTAGCCTACGTACGT` | `ACGGACTAGCCT` (already reverse-complemented — *not* what the tool expects) | 1 (`CTAG`, which matches in both orientations) | 1 |
| `ACGTACGTACGTACGTACGT` | — / standard scaffold | 10 / 10 | 10 / 10 |

### Nearest-neighbour Tm used by the Azimuth features

`AzimuthRuleSet2.MeltingTemp` (inlined DNA_NN3 table, salt-correction method 5, dnac1 = dnac2 = 25 nM,
Na⁺ = 50 mM — the parameters azimuth passes to Biopython `MeltingTemp.Tm_NN`) was compared with the canonical
`ThermoConstants.CalculateNearestNeighborTm` over **821 840** comparisons: 200 000 random 30-mers × the four
azimuth windows `[0,30) [19,24) [11,19) [6,11)`, plus every 2..7-mer exhaustively — **0 bitwise differences**.
Both reproduce Biopython 1.88 exactly, e.g.

| Sequence | Biopython `Tm_NN(..., nn_table=DNA_NN3)` | Seqeron |
|---|---|---|
| `ACGTACGTACGTACGTACGTACGTAGGACG` | 62.55983015288416 | 62.55983015288416 |
| `TAGG` | -67.90079099878676 | -67.90079099878676 |
| `GCGCGCGC` | 36.45901389113385 | 36.45901389113385 |

The duplicated table was therefore deleted in favour of the canonical call (closes the B07 cross-batch request).

---

## Divergences recorded (not defects)

- The library's legacy `GuideRnaCandidate.SelfComplementarityScore` (normalised complementary-pair fraction with
  an `L²` denominator and a 0.3 threshold) has no published source. It is kept for backward compatibility and
  documented as the older metric; the sourced CHOPCHOP stem count (`SelfComplementaryStems`) is now reported
  alongside it, and only the latter can filter.
- CHOPCHOP adds `folding × SCORE['FOLDING']` to its own composite score. Seqeron reports the stem count and
  filters on it only when the caller sets `MaxSelfComplementaryStems` (CHOPCHOP's own CLI default is no filter),
  because Seqeron's 0..100 composition score is not CHOPCHOP's score and mixing the two scales would be unsourced.
