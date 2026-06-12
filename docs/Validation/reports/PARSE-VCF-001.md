# Validation Report: PARSE-VCF-001 — VCF (Variant Call Format) Parser

- **Validated:** 2026-06-12   **Area:** FileIO
- **Canonical method(s):** `VcfParser.Parse(content)`; variants `ParseFile`, `ParseWithHeader`, `ParseFileWithHeader`, `ClassifyVariant`, and the supporting filter/genotype/statistics/INFO/writer helpers.
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.IO/VcfParser.cs`
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/VcfParserTests.cs`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **SAMtools / GA4GH hts-specs, VCF v4.3** (`https://samtools.github.io/hts-specs/VCFv4.3.pdf`) — authoritative spec. Confirmed via the spec text (PDF body is binary; corroborated through Wikipedia + GATK + web extraction of spec quotes):
  - **8 fixed mandatory columns, in order:** CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO; then optional FORMAT + one column per sample.
  - **POS is 1-based:** "the reference position, with the 1st base having position 1" (contrast BED, which is 0-based). No off-by-one transformation is applied or expected.
  - **ALT multiallelic:** comma-separated list of alternate alleles (e.g. `G,T`).
  - **`.` = missing** for ID/ALT/QUAL/etc.
  - **FILTER:** `PASS` = passed all filters; a `;`-separated list of codes = failed filters; `.` = filters not applied (distinct from PASS).
  - **INFO:** `;`-separated `key=value` (value may be a `,`-list); flag fields (Number=0) appear with no `=`.
  - **GT:** allele indices `0`=REF, `1`=first ALT, `2`=second ALT…; `/` unphased, `|` phased; `.` = missing allele.
  - **`*` in ALT** = allele missing due to an upstream/overlapping (spanning) deletion (reserved allele).
- **Wikipedia "Variant Call Format"** — confirms the same 8-column layout, 1-based POS, comma ALT, `.` missing, PASS vs `.`, INFO key=value + flags, and GT `0`=ref/`1`=alt with `/`,`|`,`.`.
- **Danecek et al. (2011), Bioinformatics 27(15):2156-2158** — the originating VCF paper; ALT = "comma separated list of alternate non-reference alleles", missing values = dot. Cited for Ti/Tv and multiallelic semantics.
- **GATK** — "If the FILTER value is '.', then no filtering has been applied"; basis for FILTER `.` ≠ PASS.

### Hand-constructed line (1-based POS + GT decode)
`chr2\t300\trs456\tG\tA,C\t80\tPASS\tDP=40\tGT:DP\t1/2:20\t0/1:25`
- CHROM=chr2, **POS=300 (1-based, stored verbatim)**, ID=rs456, REF=G, **ALT=[A, C]** (multiallelic split on `,`), QUAL=80, FILTER=[PASS], INFO={DP:40}.
- Sample1 GT=`1/2` → both alleles are ALTs (1=A, 2=C) → heterozygous, not hom-ref/hom-alt. Sample2 GT=`0/1` → 0=REF(G), 1=ALT(A) → heterozygous. Matches spec GT decoding.

### Findings / divergences
- **Note 1 (labeling, sourced):** The model collapses three spec categories — angle-bracket symbolic alleles (`<DEL>`), breakend notation (`[`/`]`), and the spanning-deletion `*` — into one `VariantType.Symbolic` value. The VCF spec treats `*` as a reserved allele and breakends as their own notation; calling all three "Symbolic" is a coarser-but-defensible classification, explicitly documented and sourced in the TestSpec (§5.4/§5.5) and Evidence doc. Not a correctness defect; recorded as PASS-WITH-NOTES.
- **Note 2 (scope):** `GetVariantLength = |len(ALT) − len(REF)|` is the standard left-aligned representation length, correct for the simple/normalized indel forms tested. Not a spec defect.

## Stage B — Implementation

### Code path reviewed
- `Parse` (VcfParser.cs:116-148): skips `##` meta lines and blank lines; reads sample names from the `#CHROM` line when >9 tab columns; delegates each data line to `ParseLine`.
- `ParseLine` (:299-344): splits on `\t`; **returns null if <8 columns** (malformed skip); **POS via `int.TryParse` (invariant culture), null on failure** → bad POS line skipped; **POS stored verbatim, 1-based, no arithmetic**; ALT `fields[4].Split(',')` (multiallelic); QUAL `.`→null else parsed; FILTER `.`→`Array.Empty` else `Split(';')`; INFO via `ParseInfo`; FORMAT/samples zipped by `:` only when >9 columns and sample names present.
- `ParseInfo` (:346-367): `;`-split; `eqIdx>0` → `key=value`; **else flag stored as `key="true"`** (handles INFO flags without `=`).
- `ClassifyVariant` (:376-402): `*`/`<…>`/`[`/`]` → Symbolic; equal length →SNP(len1)/MNP(>1); 1→>1 Insertion; >1→1 Deletion; else Complex.
- GT helpers (:518-563): `Replace('|','/')` normalizes phased→unphased; `IsHet` returns false if any allele is `.` (missing → indeterminate); `IsHomAlt` requires equal, non-`0`, non-`.`; `IsHomRef` checks `0/0`|`0|0`.
- `ParseHeader` (:174-253): `##fileformat`, `##INFO/FORMAT/FILTER` structured metadata, `#CHROM` sample-name extraction; quote-aware `ParseMetadataLine` so commas inside `Description="..."` are not split.

### Formula realised correctly? (evidence)
Yes. Recomputed the SimpleVcf worked example against the code by tracing and by tests:
| Field | Record 0 (chr1 100 …) | Code result | Match |
|-------|----------------------|-------------|-------|
| POS (1-based) | 100 | 100 | ✓ |
| ALT split | G | [G] | ✓ |
| ALT multiallelic (rec 2) | A,C | [A, C] | ✓ |
| FILTER `.` (rec 1) | . | [] (empty) | ✓ |
| INFO | DP=50 | {DP:"50"} | ✓ |
| GT Sample1 | 0/1 | "0/1" → IsHet=true | ✓ |
| GT Sample2 | 1/1 | "1/1" → IsHomAlt=true | ✓ |
| GT (rec1 S1) | 0/0 | IsHomRef=true | ✓ |
| DP sample | 25 | GetReadDepth=25 | ✓ |

Ti/Tv (CalculateTiTvRatio): iterates all ALT alleles, A→G/C→T transitions vs A→T/G→C transversions → 1.0; multiallelic A→`G,T` → 1.0. Matches Danecek multiallelic handling.

### Edge cases (Stage-A) in code
- Multiallelic ALT split — ✓ (`Split(',')`, tested M05, 1000G multiallelic).
- Missing `.` (QUAL→null, ID→`.`, FILTER→empty, INFO→{}) — ✓ (M37–M39).
- INFO flag without `=` (`DB`) — ✓ stored `key="true"` (M36).
- Phased vs unphased GT — ✓ normalized (IsHom/IsHet phased tests).
- Missing GT allele (`0/.`, `./.`, `.`) — ✓ IsHet=false (M43).
- No-sample VCF — ✓ Samples null, genotype/depth helpers return null/false safely (null-guard mutation-kill tests).
- Malformed line (<8 cols, non-integer POS) — ✓ skipped (M41).
- FILTER `.` ≠ PASS — ✓ FilterPassing/PassingCount only count exact `PASS` (M44, M46).
- Spanning deletion `*` — ✓ Symbolic (M45).

### Variant/delegate consistency
`IsSNP`/`IsIndel`/`FilterSNPs`/`FilterIndels`/`FilterByType` all route through `ClassifyVariant`; statistics `PassingCount` and `FilterPassing` use the identical PASS predicate. Consistent.

### Test quality audit
94 Vcf tests assert exact sourced values (positions, allele arrays, GT decodings, Ti/Tv ratios, real 1000G/ClinVar rsIDs and coordinates), not tautologies. 1-based POS is explicitly locked (`Parse_VcfSpecification_1BasedCoordinates`, POS=1 = first base). Edge cases all covered.

### Findings / defects
None. No off-by-one in POS, multiallelic ALT is split, GT decoded per spec, INFO flag parsing handles the missing-`=` case, FILTER PASS/`.` distinction correct.

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** — spec faithfully described; only divergence is the coarse `Symbolic` bucket for `*`/breakend/`<…>`, which is documented and sourced (not a defect).
- **Stage B: PASS** — code realises the validated spec; all worked examples and edge cases reproduce.
- **End-state: CLEAN** — no defect found; no code change required.
- **Tests:** filter `~Vcf` = 94 passed, 0 failed. Full suite `Seqeron.Genomics.Tests` = 4486 passed, 0 failed (baseline preserved).
