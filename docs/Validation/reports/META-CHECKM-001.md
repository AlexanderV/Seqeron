# Validation Report: META-CHECKM-001 — CheckM Marker-Gene Completeness/Contamination

- **Validated:** 2026-06-25   **Area:** Metagenomics
- **Canonical method(s):** `EstimateBinQualityFromMarkerCounts`, `EstimateBinQualityFromMarkers`, `DetectMarkers`, `LoadBundledRibosomalMarkerHmms`, `LoadBundledBacterialMarkerHmms`, `LoadBundledArchaealMarkerHmms`, `Bundled{Ribosomal,Bacterial,Archaeal}MarkerSets`, `LoadMarkerHmms`
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ✅ PASS
- **State:** ✅ CLEAN

Source: `src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs` (region "CheckM-style marker-gene completeness / contamination", lines 1840–2231).
Tests: `tests/Seqeron/Seqeron.Genomics.Tests/MetagenomicsAnalyzer_MarkerGeneQuality_Tests.cs` (20 tests).

## Authoritative sources opened this session

1. **Parks DH et al. (2015) "CheckM: assessing the quality of microbial genomes recovered from
   isolates, single cells, and metagenomes." *Genome Res* 25:1043–1055** — paper text
   (genome.cshlp.org abstract/methods is paywalled; method prose obtained via the publisher
   abstract and corroborating summaries): completeness = unique single-copy genes present /
   unique single-copy genes expected; contamination = single-copy genes present in multiple
   copies (only one copy should be present per genome); collocated marker genes are grouped into
   **marker sets** to refine estimates ("36% of marker genes, on average, are grouped into a set
   with one or more other marker genes"); lineage-specific collocated marker sets give robust
   estimates across bacterial/archaeal lineages.
2. **Ecogenomics/CheckM `checkm/markerSets.py` — `MarkerSet.genomeCheck(hits, bIndividualMarkers)`**
   (raw GitHub, retrieved verbatim this session) — the exact arithmetic (see Stage A).
3. **Parks DH et al. (2018) "A standardized bacterial taxonomy…" *Nat Biotechnol* 36:996 (GTDB)**
   + GTDB release-214 Pfam HMM directory listing (data.gtdb.ecogenomic.org) — confirms the bundled
   bac120 Pfam-subset accessions (PF00380, PF00410, PF00466, PF01025, PF02576, PF03726) all appear
   in the GTDB Pfam marker HMM set; bac120 = 120 markers (Pfam + TIGRFAM), ar122 = 122 markers.
4. **Pfam = CC0** (InterPro documentation) — bundled HMMs are public-domain redistributable.

## Stage A — Description

### Formula (markerSets.py "Marker Set Mode", `bIndividualMarkers=False`), verbatim:

```python
comp = 0.0 ; cont = 0.0
for ms in self.markerSet:
    present = 0 ; multiCopy = 0
    for marker in ms:
        count = len(hits.get(marker, []))
        if count == 1:            present += 1
        elif count > 1:           present += 1 ; multiCopy += (count - 1)
    comp += float(present) / len(ms)
    cont += float(multiCopy) / len(ms)
percComp = 100 * comp / len(self.markerSet)
percCont = 100 * cont / len(self.markerSet)
```

i.e. for marker sets `M`, with `present_s = |s ∩ G_M|` (markers in `s` found ≥1×) and
`multiCopy_s = Σ_{g∈s}(N_g − 1)` over markers found `N_g ≥ 1` times:

- **Completeness = 100 · (1/|M|) · Σ_{s∈M} present_s/|s|**
- **Contamination = 100 · (1/|M|) · Σ_{s∈M} multiCopy_s/|s|**

This matches the paper's prose (per-set normalisation, then averaging across collocated sets;
contamination = excess copies, `N−1` per multi-copy marker). PASS.

### Edge-case semantics
- |M| = 0 → undefined (0/0); both metrics return 0 to avoid div-by-zero (defined, sourced as a guard).
- Empty set `s` (|s| = 0) → would divide by zero; excluded from |M| (matches CheckM, which never
  emits empty sets).
- Multi-copy marker still counts **once** toward `present` (completeness) and `N−1` toward
  contamination — exactly the markerSets.py `count > 1` branch.
- Contamination can exceed 100% under heavy duplication (no cap) — correct per formula.

### Independent cross-check (hand computation, NOT from repo)

Synthetic bin re-derived directly from the markerSets.py arithmetic:
`M = {{A,B},{C,D,E},{F}}`, copies `A=1,B=0,C=2,D=1,E=1,F=1`.

| set | present/|s| | multiCopy/|s| |
|-----|-------------|---------------|
| {A,B}   | 1/2 | 0/2 |
| {C,D,E} | 3/3 | 1/3 |
| {F}     | 1/1 | 0/1 |

`comp = 0.5+1+1 = 2.5 → 100·2.5/3 = 250/3 = `**`83.3333…%`**;
`cont = 0+1/3+0 = 1/3 → 100·(1/3)/3 = 100/9 = `**`11.1111…%`**.
Re-run of a standalone Python transcription of `genomeCheck` produced
`Completeness = 83.33333333333333`, `Contamination = 11.111111111111109` — confirms the
memory value independently of the repo.

**Stage A verdict: PASS.** Formula matches the reference implementation and the paper.

## Stage B — Implementation

### Code realises the formula
`EstimateBinQualityFromMarkerCounts` (lines 1986–2035) loops over sets, computes
`present` (count ≥ 1) and `multiCopy += count−1` (count > 1), accumulates `present/ids.Count`
and `multiCopy/ids.Count`, then `100·comp/usedSets`, `100·cont/usedSets`. This is the
markerSets.py "Marker Set Mode" exactly (`count==1→present`, `count>1→present & +N−1` is
logically identical to `count≥1→present`, `count>1→+N−1`). Empty/null sets skipped; |M|=0→0/0.

`DetectMarkers` (2047–2072) scores each protein with `Plan7ProfileHmm.ViterbiScore`, converts
nats→bits (`/ln2`), and counts a hit when `bits ≥ marker.BitScoreThreshold` (the Pfam GA1
gathering threshold) — the copy count per marker. `EstimateBinQualityFromMarkers` composes
`DetectMarkers` + `EstimateBinQualityFromMarkerCounts`.

### Independent HMM detection cross-check (pyhmmer 0.12.1, real HMMER engine)

`hmmsearch` of the **bundled** HMMs against the test proteomes, with `bit_cutoffs="gathering"`:

| HMM (bundled) | GA reported | E. coli uS8 (P0A7W7) | E. coli GrpE (P09372) |
|---------------|-------------|----------------------|------------------------|
| PF00410 Ribosomal_S8 | 24.0 | hit, 166.58 bits, included | no hit |
| PF01025 GrpE         | 25.8 | no hit               | hit, 163.42 bits, included |
| PF00380 Ribosomal_S9 | 22.5 | no hit               | no hit |

The C# `DetectMarkers` reproduces this exactly: `DetectMarkers_EcoliuS8_HitsOnlyPF00410` and
`DetectMarkers_EcoliGrpE_HitsOnlyPF01025InBac120` assert PF00410→1 / PF01025→1 with all other
bundled families →0. The GA thresholds pinned in the C# tests (PF00410 = 24.0, PF01025 = 25.8)
match pyhmmer's reported gathering cutoffs. (uS8 is the canonical single-copy ribosomal-protein
marker; detection confirmed against the reference HMMER pipeline.)

### Cross-verification recomputed vs code (test run)

| Scenario | Expected (source) | C# test | Result |
|----------|-------------------|---------|--------|
| Synthetic bin {A,B},{C,D,E},{F} | 250/3 %, 100/9 % | `M-FORMULA` | PASS |
| All single-copy present | 100 %, 0 % | `M-ALLPRESENT` | PASS |
| Nothing found | 0 %, 0 % | `M-ALLABSENT` | PASS |
| Triplicated singleton | 100 %, 200 % | `M-TRIPLICATE` | PASS |
| Empty set excluded from \|M\| | \|M\|=1, 100 %, 0 % | `M-EMPTYSET` | PASS |
| \|M\|=0 | 0 %, 0 % | `M-NOMARKERSETS` | PASS |
| uS8 only / 9 ribo sets | 100/9 %, 0 % | `H-COMPLETENESS` | PASS |
| 2× uS8 | 100/9 %, 100/9 % | `H-DUPLICATE` | PASS |
| GrpE only / 6 bac120 sets | 100/6 %, 0 % | `B-COMPLETENESS` | PASS |
| uS8 / 35 ar122 sets | 100/35 %, 0 % | `A-TRUEPOSITIVE` | PASS |

### Variant / loader consistency
- `LoadBundledRibosomalMarkerHmms` → 9 markers; PF00410 GA=24.0, LENG=125 (asserted).
- `LoadBundledBacterialMarkerHmms` → 6 bac120 Pfam markers {PF00380,PF00410,PF00466,PF01025,
  PF02576,PF03726}; PF01025 GA=25.8, LENG=165 (asserted). Matches GTDB Pfam listing.
- `LoadBundledArchaealMarkerHmms` → 35 ar122 Pfam markers (asserted equivalent to the listed set).
- `LoadMarkerHmms` (caller-supplied reader path) round-trips a profile, keyed by accession
  (PF00410.25), and detects uS8.
- Original 9-marker accessor preserved (regression guard `B-DEFAULTS-UNCHANGED`).

### Test-quality audit
All assertions trace to the CheckM formula / hand computation / pyhmmer (exact rationals like
250/3, 100/9, 100/6, 100/35; GA=24.0/25.8; LENG=125/165), not code echoes. Null-argument throws
covered for every public entry point. Stage-A edge cases (complete→100/0, missing lowers
completeness, duplication raises contamination, empty set, |M|=0, collocation sets, bac vs ar
loaders) are all covered. No green-washing observed.

Full suite (`dotnet test Seqeron.sln -c Debug`): Failed 0 (see ledger). Marker-quality fixture: 20/20 pass.

## Verdict & follow-ups

**Stage A PASS / Stage B PASS → State CLEAN.** The CheckM completeness/contamination formula is
implemented verbatim from `markerSets.py`, independently reconfirmed by hand (83.33%/11.11%), and
the real-marker detection path matches the reference HMMER engine (pyhmmer). No defect.

**Documented boundary (acceptable, not LIMITED):** the full per-lineage CheckM `checkm_data` DB
(lineage-specific collocated marker sets + reference genome tree for tree-based placement, plus
operon-based collocation grouping) and the **TIGRFAM-defined** bac120/ar122 markers (CC BY-SA 4.0)
are **caller-supplied** via `LoadMarkerHmms` — they are gated/non-redistributable, not a defect.
The CheckM **formula** and the **CC0 Pfam** domain-level markers are validated here; bundled sets
use one singleton collocated set per Pfam family (operon collocation requires the lineage DB).
