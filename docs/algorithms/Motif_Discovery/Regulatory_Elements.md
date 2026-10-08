# Regulatory Element Scan

| Field | Value |
|-------|-------|
| Algorithm Group | Motif Discovery (Matching) |
| Test Unit ID | MOTIF-REGULATORY-001 |
| Related Projects | Seqeron.Genomics.Analysis, Seqeron.Genomics.Core |
| Implementation Status | Production |
| Last Reviewed | 2026-09-29 |

## 1. Overview

`FindRegulatoryElements` scans a DNA sequence for a fixed library of well-characterised regulatory consensus motifs (eukaryotic and prokaryotic promoter elements, translation-initiation signals, a polyadenylation signal, and several transcription-factor binding sites). Each library entry is matched as its published consensus string in IUPAC nucleotide code (degenerate positions kept, e.g. Kozak `R`, AP-1 `S`, κB `R/N/W/Y`, E-box `N`) and every occurrence is reported with its 0-based start position. The algorithm is specification-driven: detection is exact pattern matching against published consensus sequences, not a probabilistic score, so a hit means the input literally contains the cited consensus (or a member of its IUPAC family) [1][2][3][4][5][6][7][8][9].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Regulatory elements are short, conserved DNA sequences recognised by the transcription/translation machinery. The "consensus" of an element is the most-frequent base at each aligned position across many functional examples; individual instances may vary, but the consensus is the canonical signature used for recognition [2]. The library mixes eukaryotic core-promoter elements (TATA box, CCAAT box, GC box) [1][3], prokaryotic σ70 promoter hexamers (-10 Pribnow box, -35 box) [2], a bacterial ribosome-binding site (Shine-Dalgarno), the eukaryotic translation-initiation Kozak context [4], the poly(A) signal [5], and four transcription-factor binding sites (E-box, AP-1, NF-κB, CREB) [6][7][8][9].

### 2.2 Core Model

For a sequence `S` of length `n` and a consensus pattern `P` of length `m`, the scan reports every start index `i` with `0 <= i <= n - m` such that for all `j`, `S[i+j]` is in the IUPAC base set of `P[j]`. Plain bases match themselves; `N` matches any of A/C/G/T; the remaining IUPAC ambiguity codes match their defined subsets (NC-IUB 1985; canonical `IupacHelper.MatchesIupac`, identical to Biopython `Bio.SeqUtils.nt_search`). The library consensus strings are (5'→3'):

| Element | Pattern | Source |
|---------|---------|--------|
| TATA Box | `TATAAA` | Bucher (1990) [1] |
| CAAT Box | `CCAAT` | Bucher (1990) [1] |
| GC Box | `GGGCGG` | Dynan & Tjian (1983); Gidoni, Dynan & Tjian (1984) [3] |
| -10 Box (Pribnow) | `TATAAT` | Harley & Reynolds (1987) [2] |
| -35 Box | `TTGACA` | Harley & Reynolds (1987) [2] |
| Kozak | `GCCGCCRCCATGG` (= GCCGCC(A/G)CCATGG) | Kozak (1987) [4] |
| Shine-Dalgarno | `AGGAGG` | Shine & Dalgarno (1974) [10] |
| Poly(A) Signal | `AATAAA` | Proudfoot & Brownlee (1976) [5] |
| E-box | `CANNTG` | Massari & Murre (2000) [6] |
| AP-1 (TRE) | `TGASTCA` (= TGA(C/G)TCA) | Lee, Mitchell & Tjian (1987) [7]; Angel et al. (1987) [11] |
| NF-κB | `GGGRNWYYCC` | Gilmore (2006) [12]; includes Sen & Baltimore (1986) site GGGACTTTCC [8] |
| CREB | `TGACGTCA` | Montminy et al. (1986) [9] |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every reported element's matched `Sequence` has length equal to its `Pattern` and occurs at `Position` (0-based) in the input. | scan emits `S[i..i+m)` only at matching `i` |
| INV-02 | Each matched `Sequence` IUPAC-matches its `Pattern`. | per-position IUPAC membership test |
| INV-03 | The scan is exhaustive: all matching start indices `0 <= i <= n-m` are reported. | linear window over every offset |
| INV-04 | Consensus strings equal the published values; no fabricated constants. | each constant cites its primary source [1]–[9] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `DnaSequence` | required | DNA sequence to scan | non-null; A/C/G/T (case-insensitive, normalised to upper) |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Name` | `string` | Human-readable element name (e.g. "TATA Box", "AP-1"). |
| `Position` | `int` | 0-based start index of the occurrence in the sequence. |
| `Sequence` | `string` | The matched substring (length = pattern length). |
| `Pattern` | `string` | The consensus pattern that matched (IUPAC). |
| `Description` | `string` | Short biological role label. |

Returns an `IEnumerable<RegulatoryElement>`; elements are yielded library-entry by library-entry, in increasing position within each entry.

### 3.3 Preconditions and Validation

Null `sequence` → `ArgumentNullException`. Empty sequence → empty result. Matching is case-insensitive (`DnaSequence` normalises to uppercase). Coordinates are 0-based, inclusive start. Alphabet is DNA; the input must already be a valid `DnaSequence` (validated by its constructor).

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate `sequence` is non-null.
2. For each library entry `(Name, Pattern, Description)`:
3. Run the IUPAC degenerate scan of `Pattern` over the sequence.
4. For each match, yield a `RegulatoryElement` with the name, 0-based position, matched substring, pattern, and description.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The element library (§2.2 table) is the reference table; each consensus string is a named constant in `MotifFinder.KnownMotifs` carrying an inline source citation. IUPAC ambiguity codes follow the standard nucleotide code; the library uses `R` (Kozak −3, κB), `S` (AP-1 centre), `W`, `Y` (κB) and `N` (E-box, κB).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Full scan | O(n × r × m̄) | O(1) extra | n = sequence length, r = number of library entries (12), m̄ = mean pattern length; matches `O(n × r)` in the registry for bounded m |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [MotifFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs)

- `MotifFinder.FindRegulatoryElements(DnaSequence)`: scans the consensus library and yields `RegulatoryElement` records.
- `MotifFinder.KnownMotifs`: nested static class of source-cited consensus constants.
- `MotifFinder.FindDegenerateMotif(DnaSequence, string)`: the underlying IUPAC per-position scan reused for each entry.
- `MotifFinder.FindRegulatoryElements(DnaSequence, bool bothStrands)` → `StrandedRegulatoryElement(..., Strand)` (2026-09 follow-up, additive): with `bothStrands = true` the orientation-independent elements whose IUPAC pattern is not its own reverse complement — CAAT box, GC box, NF-κB — are also matched on the minus strand (pattern reverse-complemented via the canonical IUPAC `DnaSequence.GetReverseComplementString`, same `FindDegenerateMotif` path). Orientation-independent set (`OrientationIndependentRegulatoryElements`): CCAAT — "found in the forward or reverse orientation" (Mantovani 1998, NAR 26:1135); GC box — Sp1 binds GC boxes in both orientations, bidirectional SV40 transcription (Gidoni et al. 1985, Science 230:511); AP-1, NF-κB, E-box, CREB — enhancer elements, which act in either orientation (Banerji, Rusconi & Schaffner 1981, Cell 27:299). AP-1 `TGASTCA`, E-box `CANNTG`, CREB `TGACGTCA` are self-reverse-complementary (not rescanned); NF-κB `GGGRNWYYCC` is not (reverse complement `GGRRWNYCCC`). TATA, −10/−35, Kozak, Shine–Dalgarno and poly(A) are strand-specific (given strand only). Positions are forward window starts; `Sequence` is the site read 5'→3' on its own strand. Verified against Biopython `nt_search` on `seq` and `seq.reverse_complement()`.
- `MotifFinder.BucherPromoterMatrices` (`TataBox` POL012.1, `CapSignal` POL002.1 INR, `CcaatBox` POL004.1, `GcBox` POL003.1) and `MotifFinder.FindPromoterElementsByMatrix(DnaSequence, double falsePositiveRate, bool bothStrands = true)` ([MotifFinder.PromoterMatrices.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.PromoterMatrices.cs)): the Bucher (1990) count matrices as distributed in the JASPAR POLII collection (MEDLINE 2329577; identical in JASPAR 2014/2016/2018/2020, read from the pyjaspar 4.0.0 SQLite releases), turned into log2-odds PWMs with JASPAR pseudocounts (`PositionWeightMatrix.FromCounts` + `JasparPseudocounts`, = Biopython `motifs.read(f,"jaspar")` + `calculate_pseudocounts` + `.pssm`) and scanned through `ScanWithPwm` / `ScanWithPwmBothStrands` (both strands for CCAAT/GC box) at the per-matrix threshold `ScoreDistribution().ThresholdFpr(falsePositiveRate)`.

- `MotifFinder.FindSigma70Promoters(DnaSequence, maxMismatches35 = 2, maxMismatches10 = 2, minSpacer = 15, maxSpacer = 21, bothStrands = false)` → `Sigma70PromoterCandidate(Strand, Minus35Start, Minus35, Minus10Start, Minus10, Spacer, Mismatches35, Mismatches10, TotalMismatches, SpacerDeviation)` ([MotifFinder.Sigma70Promoters.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.Sigma70Promoters.cs), audit group D, part 2): bacterial σ70 −35/−10 pairing per Harley & Reynolds (1987) [2] — consensus TTGACA / TATAAT, spacer 15–21 bp with 17 ± 1 in 92 % of promoters (constants `Sigma70Minus35Consensus`, `Sigma70Minus10Consensus`, `Sigma70OptimalSpacer`). Every −35 hexamer within the mismatch limit (canonical `SequenceExtensions.HammingDistance`) is paired with every −10 hexamer within its limit at each allowed spacer; mismatch counts and |spacer − 17| are reported (no combined score is invented). Order: '+' then '−', −35 start, spacer; minus-strand boxes in forward coordinates, read 5'→3' on the minus strand. = independent Python brute force on 300 random sequences with planted boxes (358 candidates).
- `MotifFinder.PredictSigma70Promoters(DnaSequence, bothStrands = true, inVitro = false)` → `Sigma70PromoterPrediction` per TSS: port of the **Promoter Calculator v1.0** (La Fleur, Hossain & Salis 2022) [18] — the published, open σ70 model that scores a −35/−10 pair together with its spacer, UP element, extended −10, discriminator and initial transcribed region. Configurations UP (24) · 1 · −35 · spacer 15–20 · −10 · discriminator 6–10 · ITR 20; ΔG_total = ΔG_−10 + ΔG_−35 + ΔG_disc + ΔG_ITR + ΔG_ext−10 + ΔG_spacer + ΔG_UP + intercept (343 trained coefficients `free_energy_coeffs.npy` + intercept; quadratic spacer term 0.1463 s² − 4.9113 s + 41.119, minimal at s ≈ 17; dimer tables of `util.py`); minimum per TSS; rate = 42·exp(−β·ΔG_total) (β E. coli 1.636217004872062, in vitro 0.81632623). TSS follows the reference (minus strand: n − t). = the authors' Python (ported to Python 3, sklearn 1.9.1 `OneHotEncoder(sparse_output=False)`): 60 random sequences, 4634 per-TSS predictions, both strands and both β — every field bit-identical. lacUV5 example: best TSS 108, −35 TTTACA, 18-nt spacer, −10 TATAAT, ΔG −3.044998584365458, rate 6123.861140216731.
- `GenomeAnnotator.FindPromoterMotifs` (Annotation, batch B11) still reports unpaired −35/−10 hexamer prefixes; pairing is available here (cross-batch request recorded in B05.md).

### 5.2 Current Behavior

Each library entry is scanned independently via `FindDegenerateMotif`, so results are grouped by entry rather than globally sorted by position. The suffix tree (`SuffixTree.FindAllOccurrences`) was **not** used: the E-box pattern `CANNTG` is IUPAC-degenerate (not a single exact substring), and the library is a small, fixed set of short patterns scanned once over typically short promoter regions, so a single linear IUPAC scan is the correct and simplest algorithm. The exact-substring suffix tree would require enumerating all 16 E-box expansions and offers no asymptotic benefit at this scale.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- All 12 consensus strings copied from their sources (§2.2): TATAAA [1], CCAAT [1], GGGCGG [3], TATAAT/TTGACA [2], GCCGCCRCCATGG [4], AGGAGG [10], AATAAA [5], CANNTG [6], TGASTCA [7][11], GGGRNWYYCC [12], TGACGTCA [9].
- Exhaustive 0-based exact/IUPAC occurrence reporting (INV-01..INV-03).

**Intentionally simplified:**

- `FindRegulatoryElements(DnaSequence)` scans TATA / CCAAT / GC box as their core consensus strings (`TATAAA`, `CCAAT`, `GGGCGG`); weight-matrix detection of the Bucher (1990) elements is `FindPromoterElementsByMatrix` (above). Bucher's own cut-off values (e.g. −8.16 for the TATA box, WebSearch record of the paper) are defined on his smoothed natural-log weight scale; that transformation (smoothing constant, normalisation) was not obtainable (paper, EPD and JASPAR sites blocked), so thresholds are set by background false-positive rate instead.
- `FindRegulatoryElements(DnaSequence)` scans the given strand only (unchanged); reverse-orientation hits of the orientation-independent elements are available via `FindRegulatoryElements(sequence, bothStrands: true)`.

(2026-09 review: Kozak, AP-1 and NF-κB were previously scanned as single representative strings `GCCGCCACCATGG`, `TGACTCA`, `GGGACTTTCC`; they now use the published IUPAC consensus.)

**Not implemented:**

- Bucher's native weight scale and cut-offs (see above).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | AP-1 consensus corrected `TGAGTCA` → `TGACTCA` | Deviation (fix) | prior value reported wrong AP-1 sites and missed real ones | fixed | Lee, Mitchell & Tjian (1987) [7] |
| 2 | Added -10 (`TATAAT`) and -35 (`TTGACA`) prokaryotic hexamers | Deviation (addition) | prokaryotic promoters now detected | fixed | Harley & Reynolds (1987) [2] |
| 3 | Item 1 reconsidered (2026-09): `TGAGTCA` is the collagenase TRE and the reverse complement of `TGACTCA`; AP-1 is now the published consensus `TGASTCA` (both reported) | Fix | TREs written as TGAGTCA were missed | fixed | [7][11] |
| 4 | Kozak `GCCGCCACCATGG` → `GCCGCCRCCATGG`; NF-κB `GGGACTTTCC` → `GGGRNWYYCC` (published IUPAC consensus instead of one representative string) | Fix | G at −3 Kozak contexts and variant κB sites were missed | fixed | [4][12] |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null sequence | `ArgumentNullException` | contract |
| Empty sequence | empty result | no offset satisfies `0 <= i <= n-m` |
| Multiple occurrences of one element | all reported with their positions | INV-03 exhaustiveness |
| Degenerate E-box (`CACGTG`, `CAGCTG`, …) | matched as `CANNTG` | IUPAC `N` membership |
| Sequence with no consensus | empty result | no match |

### 6.2 Limitations

`FindRegulatoryElements` detects only the fixed library of consensus strings; it is not a general motif discovery method and does not score partial matches (use `FindPromoterElementsByMatrix` / PWM scanning) and scans strand-specific elements on the given strand only. The −35 and −10 hexamers are reported independently there; spacing-constrained σ70 pairing with mismatches is `FindSigma70Promoters`, and a quantitative σ70 promoter model (free energy, transcription rate) is `PredictSigma70Promoters` (Promoter Calculator v1.0, trained on E. coli; other organisms use the E. coli β).

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var seq = new DnaSequence("GGGTATAAAGGG");
var hits = MotifFinder.FindRegulatoryElements(seq).ToList();
// hits[0]: Name="TATA Box", Position=3, Sequence="TATAAA", Pattern="TATAAA"
// "AATGAGTCAGG" → AP-1 at 2, Sequence="TGAGTCA", Pattern="TGASTCA" (Biopython nt_search → [2])
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [MotifFinder_FindRegulatoryElements_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/MotifFinder_FindRegulatoryElements_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [MOTIF-REGULATORY-001-Evidence.md](../../../docs/Evidence/MOTIF-REGULATORY-001-Evidence.md)
- Related algorithms: [Shared_Motifs](Shared_Motifs.md), [Known_Motif_Search](../Motif_Analysis/Known_Motif_Search.md)

## 8. References

1. Bucher P. 1990. Weight matrix descriptions of four eukaryotic RNA polymerase II promoter elements derived from 502 unrelated promoter sequences. J Mol Biol 212(4):563-578. https://doi.org/10.1016/0022-2836(90)90223-9
2. Harley C.B., Reynolds R.P. 1987. Analysis of E. coli promoter sequences. Nucleic Acids Res 15(5):2343-2361. https://doi.org/10.1093/nar/15.5.2343
3. Dynan W.S., Tjian R. 1983. The promoter-specific transcription factor Sp1 binds to upstream sequences in the SV40 early promoter. Cell 35:79-87; Gidoni D., Dynan W.S., Tjian R. 1984. Multiple specific contacts between a mammalian transcription factor and its cognate promoters. Nature 312:409-413. (Lundin et al. 1994, previously cited here, concerns yeast MIG1, not Sp1.)
4. Kozak M. 1987. An analysis of 5'-noncoding sequences from 699 vertebrate messenger RNAs. Nucleic Acids Res 15(20):8125-8148. https://doi.org/10.1093/nar/15.20.8125
5. Proudfoot N.J., Brownlee G.G. 1976. 3' non-coding region sequences in eukaryotic messenger RNA. Nature 263:211-214. https://doi.org/10.1038/263211a0
6. Massari M.E., Murre C. 2000. Helix-loop-helix proteins: regulators of transcription in eucaryotic organisms. Mol Cell Biol 20(2):429-440. https://doi.org/10.1128/MCB.20.2.429-440.2000
7. Lee W., Mitchell P., Tjian R. 1987. Purified transcription factor AP-1 interacts with TPA-inducible enhancer elements. Cell 49(6):741-752. https://doi.org/10.1016/0092-8674(87)90612-X
8. Sen R., Baltimore D. 1986. Multiple nuclear factors interact with the immunoglobulin enhancer sequences. Cell 46(5):705-716. https://doi.org/10.1016/0092-8674(86)90346-6
9. Montminy M.R., Sevarino K.A., Wagner J.A., Mandel G., Goodman R.H. 1986. Identification of a cyclic-AMP-responsive element within the rat somatostatin gene. PNAS 83(18):6682-6686. https://doi.org/10.1073/pnas.83.18.6682
10. Shine J., Dalgarno L. 1974. The 3'-terminal sequence of Escherichia coli 16S ribosomal RNA: complementarity to nonsense triplets and ribosome binding sites. PNAS 71(4):1342-1346. https://doi.org/10.1073/pnas.71.4.1342
11. Angel P., Imagawa M., Chiu R., Stein B., Imbra R.J., Rahmsdorf H.J., Jonat C., Herrlich P., Karin M. 1987. Phorbol ester-inducible genes contain a common cis element recognized by a TPA-modulated trans-acting factor. Cell 49(6):729-739. https://doi.org/10.1016/0092-8674(87)90611-8
12. Gilmore T.D. 2006. Introduction to NF-κB: players, pathways, perspectives. Oncogene 25:6680-6684. https://doi.org/10.1038/sj.onc.1209954
13. Mantovani R. 1998. A survey of 178 NF-Y binding CCAAT boxes. Nucleic Acids Res 26(5):1135-1143. https://doi.org/10.1093/nar/26.5.1135 (CCAAT "found in the forward or reverse orientation"; WebSearch record)
14. Gidoni D., Kadonaga J.T., Barrera-Saldaña H., Takahashi K., Chambon P., Tjian R. 1985. Bidirectional SV40 transcription mediated by tandem Sp1 binding interactions. Science 230:511-517 (PMID 2996137; WebSearch record)
15. Banerji J., Rusconi S., Schaffner W. 1981. Expression of a β-globin gene is enhanced by remote SV40 DNA sequences. Cell 27:299-308 (enhancers act in either orientation; WebSearch record)
16. JASPAR POLII collection, matrices POL012.1 (TATA-Box), POL002.1 (INR), POL004.1 (CCAAT-box), POL003.1 (GC-box), MEDLINE 2329577 = [1]; read from pyjaspar 4.0.0 (PyPI) `JASPAR2014/2016/2018/2020.sqlite` (identical in all four releases).
17. Biopython 1.88 `Bio.motifs.jaspar` (`read`, `calculate_pseudocounts`), `Bio.motifs.matrix` (`search`, `distribution`), `Bio.motifs.thresholds` (`ScoreDistribution.threshold_fpr`).
18. La Fleur T.L., Hossain A., Salis H.M. 2022. Automated model-predictive design of synthetic promoters to control transcriptional profiles in bacteria. Nat Commun 13:5159. https://doi.org/10.1038/s41467-022-32829-5. Reference implementation opened: raw.githubusercontent.com/hsalis/SalisLabCode/master/Promoter_Calculator/ `Promoter_Calculator_v1_0.py`, `util.py`, `free_energy_coeffs.npy` (343 float64), `model_intercept.npy`, `README.md`.
