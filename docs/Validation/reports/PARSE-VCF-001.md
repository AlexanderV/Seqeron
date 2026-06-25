# Validation Report: PARSE-VCF-001 — VCF (Variant Call Format) Parser

- **Validated:** 2026-06-24   **Area:** FileIO
- **Canonical method(s):** `VcfParser.Parse(content)`; variants `ParseFile`, `ParseWithHeader`, `ParseFileWithHeader`, `ClassifyVariant`, and the supporting filter/genotype/statistics/INFO/writer helpers.
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.IO/VcfParser.cs`
- **Test files:** `tests/Seqeron/Seqeron.Genomics.Tests/VcfParserTests.cs`, `VcfParser_MutationKillers_Tests.cs`, `ConventionCompatibility_OptIn_Tests.cs`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **SAMtools / GA4GH hts-specs, VCF v4.3** (`https://samtools.github.io/hts-specs/VCFv4.3.pdf`) — authoritative spec (PDF body is binary; confirmed through the Wikipedia + GATK extractions below, which quote the same conventions). 8 fixed mandatory columns in order CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO; POS 1-based ("1st base having position 1"); ALT comma-separated alternate alleles; `.` = missing; FILTER PASS vs `.`; INFO `;`-separated `key=value` with Number=0 flags; GT `0`=REF/`1`=first ALT…, `/` unphased / `|` phased / `.` missing; `*` = spanning-deletion allele.
- **Wikipedia "Variant Call Format"** (fetched 2026-06-24) — quoted directly: POS is "The 1-based position of the variation on the given sequence" (first base = position 1, not 0); the 8 mandatory columns in order CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO; ALT alternates comma-separated without spaces (`G,T`); `.` represents missing/unknown data throughout; PASS = passed all filters, other codes = failed filters; INFO = `;`-separated `<key>=<data>[,data]` key-value pairs.
- **Danecek et al. (2011), Bioinformatics 27(15):2156-2158** — the originating VCF paper; cited for multiallelic ALT and Ti/Tv semantics.
- **GATK** — "If the FILTER value is '.', then no filtering has been applied"; basis for FILTER `.` ≠ PASS.

### Coordinate convention (the focus of this re-validation)
- **VCF file POS is 1-based** (confirmed verbatim from sources above).
- **`VcfParser.VcfRecord.Pos` stores the file's 1-based POS verbatim** — there is **no** 0-based internal storage in the parser, so no off-by-one transform is applied or expected on parse or write.
- The **new `Variant.VcfPosition` accessor (commit 6e900e92)** lives on a *different* type: `VariantCaller.Variant` (`src/.../Seqeron.Genomics.Annotation/VariantCaller.cs:418-431`), whose `Position` is 0-based (internal/array convention). `VcfPosition => Position + 1` converts that 0-based internal coordinate to the VCF 1-based POS. This is the variant-caller's model, not the parser's — the two are unrelated except that both honour the same 1-based VCF convention at the file boundary.

### Hand-constructed line (1-based POS + multiallelic split + GT decode)
`chr2\t300\trs456\tG\tA,C\t80\tPASS\tDP=40\tGT:DP\t1/2:20\t0/1:25`
- CHROM=chr2, **POS=300 (1-based, stored verbatim)**, ID=rs456, REF=G, **ALT=[A, C]** (multiallelic split on `,`), QUAL=80, FILTER=[PASS], INFO={DP:40}.
- Sample1 GT=`1/2` → both alleles ALTs (1=A, 2=C) → het, not hom. Sample2 GT=`0/1` → 0=REF(G),1=ALT(A) → het. Matches spec GT decoding.

### Findings / divergences
- **Note 1 (labeling, sourced):** the model collapses angle-bracket symbolic alleles (`<DEL>`), breakend notation (`[`/`]`), and the spanning-deletion `*` into one `VariantType.Symbolic` value. Coarser-but-defensible, explicitly documented/sourced in the TestSpec (§5.4/§5.5) and Evidence doc. Not a correctness defect → PASS-WITH-NOTES.
- **Note 2 (scope):** `GetVariantLength = |len(ALT) − len(REF)|` is the standard left-aligned representation length, correct for the simple/normalized indel forms tested.

## Stage B — Implementation

### Code path reviewed
- `Parse` (VcfParser.cs:116-148): skips `##` meta and blank lines; reads sample names from `#CHROM` when >9 tab columns; delegates each data line to `ParseLine`.
- `ParseLine` (:299-344): splits on `\t`; **<8 columns → null (malformed skip)**; **POS via `int.TryParse` (invariant culture), null on failure**; **POS stored verbatim, 1-based, no arithmetic** (line 306, `out int pos`); ALT `fields[4].Split(',')` (multiallelic); QUAL `.`→null else parsed; FILTER `.`→`Array.Empty` else `Split(';')`; INFO via `ParseInfo`; FORMAT/samples zipped by `:` only when >9 columns and sample names present.
- `ParseInfo` (:346-367): `;`-split; `eqIdx>0` → `key=value`; else flag stored `key="true"`.
- `ClassifyVariant` (:376-402): `*`/`<…>`/`[`/`]` → Symbolic; equal length → SNP(len1)/MNP(>1); 1→>1 Insertion; >1→1 Deletion; else Complex.
- GT helpers (:518-563): `Replace('|','/')` normalizes phased→unphased; `IsHet` false if any allele `.`; `IsHomAlt` requires equal non-`0` non-`.`; `IsHomRef` checks `0/0`/`0|0`.
- `FormatRecord` (:766-794): writes `record.Pos.ToString()` **verbatim** (line 771) → parse→write round-trips the 1-based POS exactly, no off-by-one.

### `VcfPosition` accessor verification (commit 6e900e92)
- `VariantCaller.Variant` (VariantCaller.cs:418-431): `Position` 0-based; `VcfPosition => Position + 1` (line 430).
- `ToVcfLines` (:350-369) already emits `variant.Position + 1` for the POS column (line 367) — so the serialized POS column equals `VcfPosition` (both 1-based). Confirmed by `ConventionCompatibility_OptIn_Tests.VcfPosition_MatchesToVcfLinesPosColumn` (POS=42 for Position=41) and `VcfPosition_IsOneBasedOffsetOfInternalPosition` (Position=0 ⇒ VcfPosition=1).
- **Defaults unchanged:** `VcfPosition` is a *new computed property only*; it adds no field, mutates no constructor/positional layout, and does not touch `Position` or any existing serialization. The parser (`VcfRecord.Pos`) is untouched by the commit. Verified the parser still stores POS verbatim and the full Vcf suite is unchanged.

### Formula realised correctly? (evidence — recomputed vs code)
| Field | Record | Code result | Match |
|-------|--------|-------------|-------|
| POS (1-based) | chr1 100 | 100 | ✓ |
| ALT split | G | [G] | ✓ |
| ALT multiallelic | A,C | [A, C] | ✓ |
| FILTER `.` | . | [] | ✓ |
| INFO | DP=50 | {DP:"50"} | ✓ |
| GT Sample1 0/1 | — | IsHet=true | ✓ |
| GT Sample2 1/1 | — | IsHomAlt=true | ✓ |
| DP sample | 25 | GetReadDepth=25 | ✓ |
| Variant.VcfPosition | Position=41 | 42 (=ToVcfLines POS) | ✓ |

Ti/Tv (`CalculateTiTvRatio`): iterates all ALT alleles; A→G/C→T transitions vs A→T/G→C transversions; multiallelic counted per ALT — matches Danecek.

### Edge cases (Stage-A) in code
Multiallelic ALT split ✓ (M05); `.` missing for QUAL/ID/FILTER/INFO ✓ (M37–M39); INFO flag without `=` ✓ (M36); phased vs unphased GT ✓; missing GT allele → IsHet=false ✓ (M43); no-sample VCF ✓ (null-safe helpers); malformed line (<8 cols, non-integer POS) → skipped ✓ (M41); FILTER `.` ≠ PASS ✓ (M44, M46); spanning deletion `*` → Symbolic ✓ (M45).

### Variant/delegate consistency
`IsSNP`/`IsIndel`/`FilterSNPs`/`FilterIndels`/`FilterByType` all route through `ClassifyVariant`; `PassingCount` and `FilterPassing` share the identical exact-`PASS` predicate. Consistent. `Variant.VcfPosition` and `ToVcfLines` agree (both +1).

### Numerical robustness note (scope, not a defect)
POS is parsed into a 32-bit `int` (`int.TryParse`), max ≈ 2.15e9. Sufficient for all individual human chromosomes (chr1 ≈ 249 Mbp) and the tested datasets, but a single contig POS > 2^31−1 (e.g. some concatenated/large non-human references) would fail to parse and be silently skipped. The VCF spec does not mandate 64-bit; this is a documented scope boundary, not a divergence from any tested or sourced value.

### Test quality audit
Vcf tests assert exact sourced values (positions, allele arrays, GT decodings, Ti/Tv ratios, real 1000G/ClinVar rsIDs and coordinates), not tautologies; 1-based POS explicitly locked. The new opt-in tests lock `VcfPosition = Position+1` and its equality with the `ToVcfLines` POS column. Filtered run (`~Vcf` ∪ `~ConventionCompatibility_OptIn`) = **199 passed, 0 failed**.

### Findings / defects
None. No off-by-one in POS (parser stores verbatim; caller `VcfPosition`/`ToVcfLines` apply the single +1 consistently); multiallelic ALT split; GT decoded per spec; INFO flag handling correct; FILTER PASS/`.` distinction correct; defaults unchanged by commit 6e900e92.

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** — spec faithfully described; only divergence is the coarse `Symbolic` bucket for `*`/breakend/`<…>` (documented & sourced).
- **Stage B: PASS** — code realises the validated spec; the new `Variant.VcfPosition` accessor is `Position+1`, consistent with `ToVcfLines`, with defaults unchanged; parser keeps the file's 1-based POS verbatim.
- **End-state: CLEAN** — no defect found; no code change required. POS-32-bit and coarse-Symbolic are sourced scope notes, not defects.
- **Tests:** filtered Vcf+OptIn = 199 passed, 0 failed. No code touched → full-suite baseline preserved.
