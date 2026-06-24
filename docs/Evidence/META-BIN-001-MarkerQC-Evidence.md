# Evidence Artifact: META-BIN-001 (marker-gene completeness/contamination addendum)

**Test Unit ID:** META-BIN-001
**Algorithm:** Genome binning — CheckM-style single-copy marker-gene completeness & contamination
**Date Collected:** 2026-06-25

> This addendum covers the marker-gene quality metrics added on top of the existing TNF/coverage
> binning + proxy metrics (documented in `META-BIN-001-Evidence.md`). The proxy completeness
> (length ratio) and contamination (GC stddev) and the default `BinContigs` path are unchanged.

---

## Online Sources

### CheckM paper — Parks, Imelfort, Skennerton, Hugenholtz & Tyson (2015)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC4484387/ (open-access PMC of *Genome Res* 25:1043)
**Accessed:** 2026-06-25
**Authority rank:** 1 (peer-reviewed paper)

**Retrieval:** WebSearch "CheckM completeness contamination formula collocated marker sets Parks 2015";
WebFetch of the PMC article, prompt asking for the verbatim completeness/contamination equations
from the Methods section ("Estimation of completeness, contamination, and strain heterogeneity").

**Key Extracted Points:**

1. **Completeness (Eq. 1):** the average, over all collocated marker sets `M`, of the recovered
   fraction of each set: `Σ_{s∈M} |s ∩ G_M| / |s| / |M|`, "where *s* is a set of collocated marker
   genes; *M* is the set of all collocated marker sets *s*; and *G_M* is the set of marker genes
   identified in a genome." (verbatim variable definitions quoted from the article).
2. **Contamination (Eq. 2):** `Σ_{s∈M} Σ_{g∈s} C_g / |s| / |M|`, "where *C_g* is *N* − 1 for a gene
   *g* identified *N* ≥ 1 times, and 0 for a missing gene." (verbatim).
3. **Why marker SETS:** "Marker genes that are consistently collocated within a lineage do not
   provide independent evidence of a genome's quality, so collocated marker genes were grouped into
   marker sets in order to further refine estimates of genome quality." Completeness/contamination
   are averaged per-set and then over the number of sets (reduces the influence of correlated
   markers).
4. **Marker basis:** completeness ≈ fraction of unique single-copy genes (SCGs) present;
   contamination ≈ how many SCGs are present in multiple copies (only one copy expected per genome).

### CheckM reference implementation — `MarkerSet.genomeCheck` (Ecogenomics/CheckM)

**URL:** https://raw.githubusercontent.com/Ecogenomics/CheckM/master/checkm/markerSets.py
**Accessed:** 2026-06-25
**Authority rank:** 3 (reference implementation of the authoritative tool)

**Retrieval:** WebFetch of the raw `markerSets.py`, asking for the loop that computes completeness
and contamination from marker sets + a marker-hit-count map.

**Key Extracted Points (verbatim arithmetic):**

```python
comp = 0.0; cont = 0.0
for ms in self.markerSet:
    present = 0; multiCopy = 0
    for marker in ms:
        count = len(hits.get(marker, []))
        if count == 1:   present += 1
        elif count > 1:  present += 1; multiCopy += (count - 1)
    comp += float(present) / len(ms)
    cont += float(multiCopy) / len(ms)
percComp = 100 * comp / len(self.markerSet)
percCont = 100 * cont / len(self.markerSet)
```

This pins the exact arithmetic: a multi-copy marker counts **once** toward `present` (so
completeness uses `|s ∩ G_M|`, not the copy count) AND contributes `count − 1` to `multiCopy`;
both sums are divided by the set size `|s|`, averaged over `|M|`, and ×100. Matches Eqs. 1–2.

### Universal single-copy ribosomal markers — Xu et al. (2022)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC9411617/ ("A revisit to universal single-copy
genes in bacterial genomes")
**Accessed:** 2026-06-25
**Authority rank:** 1 (peer-reviewed)

**Retrieval:** WebSearch for universal single-copy ribosomal marker Pfams; WebFetch of the article.

**Key Extracted Points:**

1. "the vast majority of USCGs in every set are annotated with the functional category J
   (Translation, ribosomal structure, and biogenesis)" — ribosomal proteins are the canonical core
   of universal single-copy gene (USCG) marker sets.
2. Confirms ribosomal proteins (S/L families) are the standard universal single-copy markers used
   for genome completeness/lineage analysis.

### Pfam HMM profiles — EMBL-EBI InterPro API

**URL (example):** https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/PF00410/?annotation=hmm
**Accessed:** 2026-06-25
**Authority rank:** 5 (curated database)

**Retrieval:** `curl -L` of `entry/pfam/<ACC>/?annotation=hmm` for nine ribosomal Pfams; the
response is a gzip of the HMMER3/f ASCII profile, decompressed verbatim (no edits). Each parsed
header confirms `HMMER3/f [3.3 | Nov 2019]`, `ALPH amino`, and a `GA` gathering threshold line:

| Pfam | Name | LENG | GA1 (bits) |
|------|------|------|-----------|
| PF00318 | Ribosomal_S2  | 216 | 26.1 |
| PF00177 | Ribosomal_S7  | 148 | 24.0 |
| PF00410 | Ribosomal_S8  | 125 | 24.0 |
| PF00380 | Ribosomal_S9  | 121 | 22.5 |
| PF00338 | Ribosomal_S10 |  98 | 21.2 |
| PF00411 | Ribosomal_S11 | 110 | 25.1 |
| PF00203 | Ribosomal_S19 |  81 | 27.0 |
| PF00687 | Ribosomal_L1  | 198 | 32.2 |
| PF00297 | Ribosomal_L3  | 369 | 22.0 |

### Pfam licence — InterPro documentation

**URL:** https://interpro-documentation.readthedocs.io/en/latest/pfam.html
**Accessed:** 2026-06-25
**Authority rank:** 2 (official documentation)

**Key Extracted Point (verbatim):** "Pfam is freely available under the Creative Commons Zero
('CC0') licence." → public domain, freely redistributable; profiles are embedded as-is.

### True-positive protein — E. coli uS8 (RpsH), UniProt P0A7W7

**URL:** https://rest.uniprot.org/uniprotkb/P0A7W7.fasta
**Accessed:** 2026-06-25
**Authority rank:** 5 (curated database)

**Retrieval:** WebFetch of the UniProt FASTA. Sequence (130 aa), the canonical small ribosomal
subunit protein uS8 — the true positive for PF00410 (Ribosomal_S8):
`MSMQDPIADMLTRIRNGQAANKAAVTMPSSKLKVAIANVLKEEGFIEDFKVEGDTKPELELTLKYFQGKAVVESIQRVSRPGLRIYKRKDELPKVMAGLGIAVVSTSKGVMTDRAARQAGLGGEIICYVA`

---

## Documented Corner Cases and Failure Modes

### From the CheckM reference implementation

1. **Multi-copy marker:** counts once toward `present` (completeness) and `N−1` toward contamination.
2. **Missing marker:** `C_g = 0`; contributes 0 to both present and multiCopy.
3. **Empty marker set:** dividing by `len(ms)` would fail — empty sets must be excluded.
4. **No marker sets (|M|=0):** both metrics undefined; return 0 (avoid div-by-zero).

### From the Plan7 / HMMER scoring model

5. **Glocal full-length score:** this engine computes a glocal (whole-sequence) Viterbi log-odds
   score; HMMER's GA thresholds were tuned with HMMER's local model + composition (null2)
   corrections, so absolute bit scores differ from `hmmsearch`. The true-positive separation is
   nonetheless decisive (uS8 vs PF00410 = +176 bits; vs all 8 other ribosomal families < 0 bits).

---

## Test Datasets

### Dataset: Synthetic bin — hand-derived CheckM formula (pins Eqs. 1–2)

**Source:** Parks et al. (2015) Eqs. 1–2 / `MarkerSet.genomeCheck` (above), computed by hand.

Marker sets `M` (3 sets) and per-marker copy counts:

| Set | Markers | Counts | present_s | multiCopy_s | comp term | cont term |
|-----|---------|--------|-----------|-------------|-----------|-----------|
| s1 | {A, B} (|s|=2) | A=1, B=0 | 1 | 0 | 1/2 = 0.5 | 0/2 = 0 |
| s2 | {C, D, E} (|s|=3) | C=2, D=1, E=1 | 3 | 1 | 3/3 = 1.0 | 1/3 |
| s3 | {F} (|s|=1) | F=1 | 1 | 0 | 1/1 = 1.0 | 0/1 = 0 |

- `comp = 0.5 + 1.0 + 1.0 = 2.5` ⇒ **Completeness = 100 · 2.5 / 3 = 250/3 = 83.33333…%**
- `cont = 0 + 1/3 + 0 = 1/3` ⇒ **Contamination = 100 · (1/3) / 3 = 100/9 = 11.11111…%**
- MarkersPresent (distinct ≥1) = {A,C,D,E,F} = 5; MarkerSetCount = 3; MarkerCount = 6.

### Dataset: HMM detection integration

**Source:** UniProt P0A7W7 (E. coli uS8) + bundled PF00410 CC0 HMM (GA1 = 24 bits).

| Marker | Protein | Engine bit score | Threshold | Hit? |
|--------|---------|------------------|-----------|------|
| PF00410 (S8, true family) | uS8 | ≈ 176.38 | 24.0 | yes (count 1) |
| PF00177, PF00318, PF00380, PF00338, PF00411, PF00203, PF00687, PF00297 | uS8 | all < 0 | their GA1 | no |

⇒ `DetectMarkers([uS8], bundled)` gives PF00410 = 1, all others = 0. With each bundled marker as
its own singleton set (|M|=9), exactly one of nine sets is fully present ⇒ Completeness = 100/9 ≈
11.111%, Contamination = 0.

---

## Assumptions

1. **ASSUMPTION: Marker present iff Viterbi bit score ≥ GA1.** CheckM calls a marker present from
   `hmmsearch` hits above the model's gathering threshold. We reuse the Pfam GA1 per-sequence
   threshold as the bit-score gate against this repo's glocal Plan7 Viterbi score. The numeric gate
   is sourced (GA1 from the HMM file); the *choice of glocal Viterbi vs HMMER local+null2 scoring*
   is a documented engine difference (not a fabricated value), kept consistent with the existing
   `ProteinMotifFinder.FindDomainsByHmm`. The exact-formula tests do not depend on this — they feed
   counts directly.
2. **ASSUMPTION: bundled markers as singleton sets.** The bundled universal ribosomal families are
   distinct, non-collocated Pfams, so each is its own marker set |s|=1. CheckM's operon-based
   collocation grouping needs the full lineage DB; left as the residual.

---

## Recommendations for Test Coverage

1. **MUST Test:** the exact hand-derived synthetic-bin Completeness = 250/3% and Contamination =
   100/9% from `EstimateBinQualityFromMarkerCounts` — Evidence: CheckM Eqs. 1–2 / `genomeCheck`.
2. **MUST Test:** all-present single-copy ⇒ 100% complete, 0% contaminated; all-absent ⇒ 0/0;
   empty marker set ignored (no div-by-zero); |M|=0 ⇒ 0/0 — Evidence: corner cases 1–4.
3. **MUST Test:** HMM detection — bundled PF00410 detects E. coli uS8 (count 1) and no other
   bundled family does; end-to-end `EstimateBinQualityFromMarkers` over the 9 singleton sets gives
   Completeness = 100/9% — Evidence: HMM integration dataset.
4. **SHOULD Test:** contamination with a triplicated marker (`C_g = N−1 = 2`); null-argument guards.
5. **COULD Test:** caller-supplied `LoadMarkerHmms` round-trips a profile and keys it by accession.

---

## References

1. Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW (2015). CheckM: assessing the
   quality of microbial genomes recovered from isolates, single cells, and metagenomes.
   *Genome Research* 25:1043–1055. https://pmc.ncbi.nlm.nih.gov/articles/PMC4484387/
2. Ecogenomics/CheckM — `checkm/markerSets.py` (`MarkerSet.genomeCheck`).
   https://raw.githubusercontent.com/Ecogenomics/CheckM/master/checkm/markerSets.py
3. Xu L, et al. (2022). A revisit to universal single-copy genes in bacterial genomes.
   https://pmc.ncbi.nlm.nih.gov/articles/PMC9411617/
4. EMBL-EBI InterPro / Pfam — HMM profiles (CC0).
   https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/PF00410/?annotation=hmm ;
   licence: https://interpro-documentation.readthedocs.io/en/latest/pfam.html
5. UniProt P0A7W7 — Small ribosomal subunit protein uS8, *E. coli*.
   https://rest.uniprot.org/uniprotkb/P0A7W7.fasta

---

## Change History

- **2026-06-25**: Initial documentation of the marker-gene completeness/contamination addendum.
