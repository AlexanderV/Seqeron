# Validation Report: ONCO-PLOIDY-001 — Tumor Ploidy Estimation + Whole-Genome-Doubling Detection

- **Validated:** 2026-06-24   **Area:** Oncology
- **Canonical method(s):**
  `OncologyAnalyzer.EstimatePloidy(IEnumerable<AlleleSpecificSegment>)`,
  `OncologyAnalyzer.DetectWholeGenomeDoubling(IEnumerable<AlleleSpecificSegment>, ReferenceGenome=GRCh38)`,
  `OncologyAnalyzer.DetectWholeGenomeDoublingFromSuppliedLength(...)` (legacy variant),
  `GetAutosomeLengths(ReferenceGenome)` / `GetAutosomalGenomeLength(ReferenceGenome)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** CLEAN

This is a re-validation after the limitations-campaign fix (commits c47e87f8 and 940995dd) which
re-based the WGD genome fraction on an embedded reference chromosome-size table (UCSC hg38/hg19,
Σchr1–22) instead of the supplied-segment total length, and added a `ReferenceGenome` selector and
table accessors. The prior report (2026-06-16) was PASS-WITH-NOTES because of exactly the
supplied-segment-length denominator that this fix has now resolved.

## Stage A — Description

### Sources opened & what they confirm (retrieved this session via web tools)

1. **facets-suite `copy-number-scores.R` `is_genome_doubled`** (raw GitHub, MSKCC, PMID 30013179 /
   Bielski 2018) — fetched verbatim this session:
   ```r
   is_genome_doubled = function(segs, chrom_info, treshold = 0.5) {
       autosomal_genome = sum(as.numeric(chrom_info$size[chrom_info$chr %in% 1:22]))
       frac_elevated_mcn = sum(as.numeric(segs$length[which(segs$mcn >= 2 & segs$chrom %in% 1:22)])) / autosomal_genome
       wgd = frac_elevated_mcn > treshold
       wgd
   }
   ```
   with `mcn = tcn - lcn` (from `parse_segs`). Confirms: denominator = reference autosomal genome
   length (Σ chrom 1–22 from a build-parameterised chromosome-size table, NOT supplied segments);
   numerator restricted to autosomes; key on **major** CN ≥ 2; strict `>` 0.5.

2. **Ensembl GRCh38 assembly REST** (`rest.ensembl.org/info/assembly/homo_sapiens`) — fetched this
   session; assembly **GRCh38.p14** (GCA_000001405.29). Spot-checked chromosome lengths:
   chr1 = 248,956,422; chr7 = 159,345,973; chr21 = 46,709,983; chr22 = 50,818,468; chrX = 156,040,895.

3. **Ensembl GRCh37 assembly REST** (`grch37.rest.ensembl.org/info/assembly/homo_sapiens`) — fetched
   this session; assembly **GRCh37.p13**. Spot-checked: chr1 = 249,250,621; chr7 = 159,138,663;
   chr21 = 48,129,895; chr22 = 51,304,566.

   (UCSC `hg38/hg19.chrom.sizes` were unreachable from this environment this session — repeated
   timeouts on `hgdownload.soe.ucsc.edu`. The embedded values were independently verified against the
   Ensembl GRCh38.p14 / GRCh37.p13 assembly REST endpoints, which carry the same authoritative
   chromosome lengths, plus a by-hand sum cross-check below.)

   Ploidy definition (Patchwork "average total copy number of all genomic segments weighted by segment
   length" → ψ = Σ(CN·L)/Σ(L), CN total = Major+Minor; ASCAT `ploidy = sum((nA+nB)*length)/sum(length)`;
   Van Loo PNAS n-scale, >2.7n aneuploidy) were re-validated against primary literature in the prior
   session and are unchanged by this fix; carried forward.

### Formula check
- Ploidy ψ = Σ((Major+Minor)·Length)/Σ(Length) — matches Patchwork prose and ASCAT code exactly.
- WGD = (Σ autosomal length where major CN ≥ 2) / (Σ chr1–22 reference length) **> 0.5** strictly —
  matches facets-suite `is_genome_doubled` verbatim (denominator now the reference autosomal genome,
  numerator autosome-restricted, major-CN key, strict operator).

### Edge-case semantics (all sourced)
- Empty ploidy set → Σ(L)=0 → reject (weighted mean undefined).
- Length ≤ 0 / negative CN → invalid input (reject).
- WGD against the fixed reference denominator: empty set → numerator 0 → fraction 0 → **false** (no
  throw); this differs from ploidy because the denominator is a constant, not Σ supplied length.
- WGD boundary exactly 0.5 → false (strict `>`).
- WGD keys on **major** CN ≥ 2 (a 1:1 / total-2 genome is NOT doubled; a 2:0 LOH or 2:2 IS).
- WGD numerator excludes sex chromosomes/contigs (`chrom %in% 1:22`).

### Independent cross-check (hand-computed this session)
| Item | Computation | Result |
|------|-------------|--------|
| GRCh38 Σ(chr1–22) | sum of 22 embedded values (Python) | **2,875,001,522** ✓ (= INV-5) |
| GRCh37 Σ(chr1–22) | sum of 22 embedded values (Python) | **2,881,033,286** ✓ (= INV-5) |
| GRCh38 half | 2,875,001,522 / 2 (even) | **1,437,500,761** (M9 exact-half value) |
| M1 ploidy | (2·100+4·100+3·50)/250 | **3.0** |
| M3 ploidy | (2·300+4·10)/310 = 640/310 | **2.0645…** (length-weighting, ≠ plain mean 3.0) |
| chr-size spot checks | Ensembl GRCh38.p14: chr1 248,956,422 / chr7 159,345,973 / chr21 46,709,983 / chr22 50,818,468; GRCh37.p13: chr1 249,250,621 / chr7 159,138,663 / chr21 48,129,895 / chr22 51,304,566 | **all match embedded table** |
| WGD just-over-half | (G/2)+1 = 1,437,500,762 / G > 0.5 | **true** |
| WGD exactly half | 1,437,500,761 / G = 0.5, not > 0.5 | **false** |
| WGD half−1 | 1,437,500,760 / G < 0.5 | **false** |

### Findings / divergences
- None. The Stage-A description (TestSpec + Evidence) is source-exact and the prior PASS-WITH-NOTES
  assumption (supplied-segment-length denominator) is correctly resolved: it is now the facets-suite
  reference autosomal genome length, with the legacy behaviour retained as a clearly-named variant.

**Verdict: PASS.**

## Stage B — Implementation

### Code path reviewed (`src/.../Seqeron.Genomics.Oncology/OncologyAnalyzer.cs`)
- `EstimatePloidy` — 6554-6578
- `DetectWholeGenomeDoubling` — 6608-6637 (reference-table overload, `genome = GRCh38` default)
- `DetectWholeGenomeDoublingFromSuppliedLength` — 6659-6687 (legacy variant)
- `GetAutosomeLengths` / `GetAutosomalGenomeLength` — 6481-6507
- `GRCh38AutosomeLengths` / `GRCh37AutosomeLengths` tables — 6411-6468
- `TryGetAutosomeNumber` (chr-prefix + 1..22 parse) — 6517-6533
- `ValidateSegment` (End>Start, CN≥0) — 2135-2151; `AlleleSpecificSegment.Length = End-Start`,
  `MajorCopyNumber` = larger allele — 1926-1936
- constants `WholeGenomeDoublingMajorCopyNumber = 2` (6377), `…FractionThreshold = 0.5` (6384)

### Formula realised correctly?
- `EstimatePloidy`: `weightedCopyNumberSum += (Major+Minor)*length`, `totalLength += length`, returns
  the quotient; rejects empty (`totalLength==0`). Matches ASCAT/Patchwork. (Unchanged by the fix.)
- `DetectWholeGenomeDoubling`: denominator = `GetAutosomalGenomeLength(genome)`; loops segments,
  validates, **skips non-autosomes** via `TryGetAutosomeNumber`, adds `Length` to `elevatedLength`
  when `MajorCopyNumber >= 2`; returns `(double)elevatedLength / autosomalGenomeLength > 0.5`. This is
  the facets-suite rule verbatim — reference denominator, autosome restriction, major-CN key, strict
  operator. Empty set ⇒ elevated 0 ⇒ false (no throw), matching the fixed-denominator semantics.
- `GetAutosomalGenomeLength` sums the 22-entry table with `long` accumulation (no overflow: ~2.9e9 ≪
  long max). `TryGetAutosomeNumber` strips a case-insensitive `chr`/`CHR` prefix then parses 1..22;
  rejects empty, sex chromosomes, contigs.
- Validation contract: null → ArgumentNullException; End≤Start or negative CN → ArgumentException —
  shared with ONCO-LOH-001/HRD-001.

### Cross-verification recomputed vs code (unit tests, 30/30 pass)
- M14/M15: embedded GRCh38/GRCh37 tables equal the spot-checked Ensembl values (and the full
  22-value arrays asserted in-test); M16: sums = 2,875,001,522 / 2,881,033,286 — match the by-hand
  Python sums above.
- M8/M9/M10 flip exactly at (G/2)±1 / (G/2) against the true GRCh38 denominator; M11 (all 1:1 → false)
  and S1 (2:0 LOH → true) lock major-vs-total; M12 (100 Mb amplified → false) proves the
  supplied-segment bias is gone; S3 (chrX/chrY excluded) and S4 (chr7 recognised) lock the autosome
  parser; M17 selects the hg19 denominator; M18 (empty → false); ploidy M1/M3/M4/S2/C1 recompute to the
  hand-derived exact values.

### Variant/delegate consistency
- `DetectWholeGenomeDoublingFromSuppliedLength` is the pre-fix Σ-supplied-length variant; L1 (60%→true),
  L2 (exactly 50%→false), L3 (empty→throw), L4 (null→throw) verify it. Property
  (`OncologyProperties.cs`) and combinatorial (`OncologyCombinatorialTests.cs`) WGD oracles were
  re-pointed at this legacy overload (their length-based oracles encode the supplied-length semantics);
  both files build with 0 warnings and the unit subset is green.

### Test quality audit (HARD gate)
- **Exact sourced values, not code-echoes:** ploidy assertions `Is.EqualTo(x).Within(1e-10)` with x
  hand-derived (3.0; 640/310; 2.0; 580/190); WGD assertions exact booleans at the strict 0.5 boundary
  against the true 2,875,001,522 bp denominator; table assertions are full 22-value arrays equal to the
  authoritative sizes.
- **Discriminating tests present:** M3 (640/310 ≠ plain-mean 3.0) catches an unweighted-mean bug;
  M9 (exactly half → false) catches a `>=` bug; M11 vs S1 catch a total-vs-major bug; M12 catches a
  supplied-length-denominator regression; S3 catches a missing autosome restriction.
- **No green-washing:** no `Greater`/`AtLeast`/`Contains`/widened tolerances/skipped tests; all
  edge/error branches asserted for both WGD overloads and ploidy (empty, non-positive length, negative
  CN, null).

### Findings / defects
- **No defect.** Code faithfully realises the Stage-A-validated formulas. The previously-noted
  supplied-length denominator is resolved (now the facets-suite reference autosomal genome).
- **Note (input contract, not a defect):** `MajorCopyNumber` is documented as the larger allele, so it
  equals facets-suite `mcn = tcn − lcn`. Acceptable, documented.
- **Note (synthetic test inputs):** boundary WGD tests use single segments far longer than any real
  chromosome (e.g. 1.43 Gb on "chr1"); the method does not bound segment length to chromosome size, so
  these are valid synthetic boundary probes — by design, not a defect.

**Verdict: PASS.**

## Verdict & follow-ups
- Stage A: **PASS**. Stage B: **PASS**.
- End-state: **CLEAN** — no defect; no code/test changes made this session.
- Unit tests: `OncologyAnalyzer_EstimatePloidy_Tests` 30/30 pass; build 0 warnings / 0 errors.
- No code changed, so the full unfiltered suite was not required to be re-run this session (the unit
  subset is green and the build is clean).
- No findings logged to FINDINGS_REGISTER (no defect).
