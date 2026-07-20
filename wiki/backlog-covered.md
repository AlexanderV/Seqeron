---
type: index
title: "Ingestion backlog — covered-via-concept history + resolution log"
tags: [meta, coverage]
created: 2026-07-18
updated: 2026-07-20
---

# Ingestion backlog — covered history

Split out of [[backlog]] on 2026-07-18 to keep that hub under the page-size cap.
This page holds the **chronological resolution-note history**. The static
**Covered via concept (done)** table was further split into
[[backlog-covered-table]] on 2026-07-20 (this page grows append-only, so the
static table was moved off it). The active Pending/Queued/Notes list stays in
[[backlog]].

## Resolution-note history

algorithm-doc pending backlog is now fully closed.
(Variants/Variant_Detection → [[germline-variant-calling-snp-indel]] resolved 2026-07-17, **closing the
Variants domain AND the entire per-domain algorithm-doc pending backlog** (RECONCILE/REUSE: the umbrella
`CallVariants` caller VARIANT-CALL-001 is the parent of the already-reconciled SNP (`SNP_Detection.md`)
and indel (`Indel_Detection.md`) facets and the core of this page; treated the spec as the canonical
PRIMARY spec (added first in `sources:`, ahead of the SNP/indel specs and the three Evidence docs) rather
than creating a redundant `variant-detection` page. Enriched a "The umbrella `CallVariants` caller" section:
the `CallVariants` / `CallVariantsFromAlignment` / `CalculateStatistics` entry points, the shared
Needleman–Wunsch `SequenceAligner.GlobalAlign` engine + `CallVariantsFromAlignment` pre-aligned string
entry point (`ArgumentException` on unequal lengths), the full five-member `VariantType` enum (only
SNP/Insertion/Deletion emitted; MNP/Complex reserved), the `VariantStatistics` record + per-kilobase
`VariantDensity` formula, and the `VcfPosition = Position+1` opt-in accessor. See `log.md` for full detail.
Contradictions: none.)
Variants/Variant_Annotation → [[variant-effect-annotation-vep]] resolved 2026-07-17 (RECONCILE/REUSE:
the variant-annotation spec VARIANT-ANNOT-001 — `VariantAnnotator.PredictFunctionalImpact` (per-variant
codon-translation consequence engine) / `Annotate` (batch, most-severe `VariantAnnotation` per variant) /
`GetImpactLevel` / `GetConsequenceRank` — is the VEP-style consequence concept **already synthesized** on
this Evidence-derived page (its VARIANT-ANNOT-001 Evidence was already in `sources:`, and the page already
carried the `OverlapConsequence` IMPACT/rank table, the codon-translation predicate set, and the worked
oracles); treated the spec as the canonical PRIMARY spec (added ahead of the two Evidence docs in
`sources:`) rather than creating a redundant `variant-annotation` page — kept distinct from the upstream
germline caller [[germline-variant-calling-snp-indel]] (VARIANT-CALL/SNP/INDEL) and the still-pending
umbrella `Variant_Detection.md` sibling. Enriched a "Method contract (VARIANT-ANNOT-001 algorithm spec)"
section: the four `VariantAnnotator` entry points + `PredictFunctionalImpact(variant, transcript,
referenceSequence, sequenceStart=1)` / `Annotate(variants, annotations, referenceSequence?, sequenceStart)`
signatures, the `FunctionalImpact` (`Consequence`/`Impact`/`CodonChange` `c.`/`AminoAcidChange` `p.`) output
+ most-severe `VariantAnnotation` INV-04, the input contract (1-based `Position`, forward-strand alleles,
`sequenceStart` window origin), the ArgumentException/ArgumentNullException + ambiguous-codon→X→
`coding_sequence_variant` + no-window→coarse-coding-term edges, O(E) / O(v×g) complexity + suffix-tree-N/A,
and the ASM-02 forward-strand-only (minus-strand + intron-split codons unhandled) + out-of-scope
SIFT/PolyPhen-invented-constants limitations. See `log.md` for full detail. Contradictions: none.)
Variants/SNP_Detection → [[germline-variant-calling-snp-indel]] resolved 2026-07-17 (RECONCILE/REUSE:
the SNP-detection spec VARIANT-SNP-001 — `VariantCaller.FindSnpsDirect` (canonical positional Hamming
enumerator over `string` inputs) / `FindSnps` (`DnaSequence` delegate = `CallVariants` filtered to
`VariantType.SNP`) — is the SNP facet already synthesized on this SNP/indel page (its VARIANT-SNP-001
Evidence was already in `sources:`, and the page carried a dedicated "SNP detection: FindSnps /
FindSnpsDirect and the Hamming-mismatch invariant" section); treated the spec as the canonical PRIMARY
spec, kept distinct from the sibling indel `Indel_Detection.md` (already reconciled) and the still-pending
umbrella `Variant_Detection.md` / `Variant_Annotation.md`; enriched the SNP section with a "Method
contract" paragraph (the two entry points incl. the `string`-input `FindSnpsDirect` vs `DnaSequence`
`FindSnps`, 0-based `Position`/`QueryPosition` + single-base REF/ALT INV-02/04, null/empty vs
ArgumentNullException contract, O(n)/O(1) vs O(n×m) complexity, suffix-tree-N/A). See `log.md` for full
detail. Contradictions: none.)
Variants/Indel_Detection → [[germline-variant-calling-snp-indel]] resolved 2026-07-17 (RECONCILE/REUSE:
the indel-detection spec VARIANT-INDEL-001 — `VariantCaller.FindInsertions` / `FindDeletions` /
`FindIndels` (union) filters over the aligned-column caller — is the indel facet already synthesized on
this SNP/indel page (its VARIANT-INDEL-001 Evidence was already in `sources:`); treated the spec as the
canonical PRIMARY spec, kept distinct from the still-pending SNP `SNP_Detection.md` and umbrella
`Variant_Detection.md` / `Variant_Annotation.md` siblings; enriched the "Indel detection" section with a
"Method contract" paragraph (the three entry points incl. the `FindIndels` union behind `find_indels`,
0-based `Position` INV-06 + `QueryPosition`, case-norm / null / empty contract, `ATGCAT→ATGTCAT` oracle,
O(n×m) NW cost, suffix-tree-N/A). See `log.md` for full detail. Contradictions: none.)
Translation/Six_Frame_Translation → [[genetic-code-translation]] resolved 2026-07-17, **closing the
Translation domain** (last pending doc) (RECONCILE/REUSE: the six-frame translation + START→STOP ORF spec
TRANS-SIXFRAME-001 — `Translator.TranslateSixFrames` (six frames keyed ±1…±3) / `FindOrfs` — is the
six-frame surface already synthesized on this page's `Translator` layer; treated the spec as the canonical
PRIMARY spec rather than creating a redundant `six-frame-translation` page, joining its already-reconciled
`Codon_Translation.md` / `Protein_Translation.md` siblings. Enriched a "Six-frame translation contract"
subsection: forward `+f` = offset `f−1` / reverse `−f` = reverse complement at offset `f−1` (Biopython
independent-offset convention), per-frame length ⌊(len−offset)/3⌋, the never-early-terminates `*`-rendering
behavior (does not honour `toFirstStop`) with one shared reverse complement, the full `OrfResult` field set
(0-based `StartPosition`, 0-based **inclusive** `EndPosition`, `Frame` ±1…±3, START-included/STOP-excluded
`Protein`, derived `NucleotideLength = End−Start+1` and `AminoAcidLength = Protein.Length`, INV-05/INV-06),
O(n) time / O(n) space (suffix tree N/A), and the `MAIVMGR*KGAR*` / `LSGTLSAAHYNGH` + `MKP` oracles).
Contradictions: none — spec, Evidence artifact, and impl agree on the ±1…±3 six-frame + START→STOP ORF
contract.)
Translation/Codon_Translation → [[genetic-code-translation]] resolved 2026-07-17 (RECONCILE/REUSE: the
single-codon lookup spec TRANS-CODON-001 — `GeneticCode.Translate` / `IsStartCodon` / `IsStopCodon` /
`GetCodonsForAminoAcid` / `GetByTableNumber` — is exactly the concept already synthesized on this
Evidence-derived page (four NCBI tables 1/2/3/11, T→U + case normalization, exactly-3-char contract);
treated the spec as the canonical PRIMARY spec rather than creating a redundant `codon-translation` page.
Key enrichment: the spec doc documents the IUPAC-ambiguity → `'X'` return as an **intentional
simplification** (§5.2/§5.3/§6.1), which **resolves** the source-vs-implementation discrepancy the page had
flagged against the older Evidence corner-case table (which expected an `ArgumentException` for `NNN`);
also added the full family-size degeneracy distribution (1/2/3/4/6-fold) and the O(1) per-codon complexity.
Kept distinct from the still-pending whole-CDS `Six_Frame_Translation.md` sibling,
both already described on this page's `Translator` layer).
Translation/Protein_Translation → [[genetic-code-translation]] resolved 2026-07-17 (RECONCILE/REUSE: the
whole-sequence framed-translation + ORF spec TRANS-PROT-001 — `Translator.Translate` (Dna/Rna/string
overloads) / `TranslateSixFrames` / `FindOrfs` — is exactly the `Translator` layer already synthesized on
this page above the single-codon `GeneticCode` lookup; treated the spec as the canonical PRIMARY spec rather
than creating a redundant `protein-translation` page. Enriched with a "Method contract (algorithm spec)"
section: the three entry-point signatures + defaults (`frame=0`, `toFirstStop=false`, `minLength=100` aa,
`searchBothStrands=true`), the null/empty contract (Dna/Rna overloads throw `ArgumentNullException`, string
overload returns empty, invalid frame → `ArgumentOutOfRangeException`, six-frame always returns all six keys),
invariants INV-01…INV-04 (frame-0 == six-frame `[1]`; length ≤ floor((len−frame)/3); exactly six keys;
`NucleotideLength == End−Start+1`), and the accepted deviations / not-implemented scope (ORFs run off the
sequence end still emitted, reverse-strand ORF coordinates in the RC scan frame needing external remap, no
nested-ORF reporting). Kept distinct from the still-pending six-frame `Six_Frame_Translation.md` sibling.)
Transcriptome/Expression_Quantification → [[expression-quantification]] resolved 2026-07-17, **closing the
Transcriptome domain** (last pending doc) (RECONCILE/REUSE: the TPM/FPKM/RPKM + quantile-normalization
quantification spec TRANS-EXPR-001 — `CalculateTPM` / `CalculateFPKM` / `QuantileNormalize` — is the concept
already synthesized on this Evidence-derived page; treated the spec as the canonical PRIMARY spec and enriched
a "Method contract (algorithm spec)" section with the three entry-point signatures + tuple/matrix I/O types,
the TPM-also-fills-FPKM-field (N = sample total raw count) and QN-uses-first-sample-length behaviors, the
explicit INV-04/INV-05 rank-order/permutation invariants, the tied-rank deviation fix, per-method complexity
(TPM O(n), FPKM O(1), QN O(s·n log n)), and the not-implemented TMM/median-of-ratios scope, kept distinct from
the downstream consumer [[differential-expression]], rather than creating a redundant page).
Transcriptome/Differential_Expression → [[differential-expression]] resolved 2026-07-17 (RECONCILE/REUSE:
the algorithm spec TRANS-DIFF-001 — `CalculateFoldChange` / `FindDifferentiallyExpressed` — is the
concept already synthesized on this Evidence-derived page; treated the spec as the canonical PRIMARY
spec and enriched a "Method contract (algorithm spec)" section with the two entry-point signatures +
defaults (alpha=0.05, log2FoldChangeThreshold=1.0), the control=Condition1 / treatment=Condition2
input contract, the already-normalized-input assumption, the six-field `DifferentialExpression` output
(incl. the `Regulation` Upregulated/Downregulated/Unchanged label), O(g·s + g·log g) complexity, and
the suffix-tree-N/A note, rather than creating a redundant page).
Transcriptome/Alternative_Splicing → [[alternative-splicing-psi]] resolved 2026-07-17 (RECONCILE/REUSE:
the algorithm spec TRANS-SPLICE-001 — `CalculatePSI` / `DetectAlternativeSplicing` — is the concept
already synthesized on this Evidence-derived page; treated the spec as the canonical PRIMARY spec and
enriched a "Method contract" section with the entry-point signatures + effective-length defaults, the
negative-read `ArgumentOutOfRangeException`, the decoupled `InclusionLevel = NaN` on detected events,
the multi-difference→SkippedExon fallback, the tolerant null/<2-isoform/identical-pair contract, and
O(g·k²·e) complexity, rather than creating a redundant `alternative-splicing` page).
StructuralVar/SV_Detection → [[discordant-pair-sv-detection]] resolved 2026-07-17, **closing the
StructuralVar domain** (last pending doc) (RECONCILE/REUSE: the paired-end-mapping (PEM) discordant-pair
SV caller, unit SV-DETECT-001, recognises SV-type signatures from read-pair span + orientation
(DEL=larger span, INS=smaller span, INV=same-strand FF/RR, DUP=everted RF, CTX=interchromosomal) and is
the concept already synthesized on this Evidence-derived page — reconciled the algorithm spec onto it as
its canonical PRIMARY spec rather than creating a redundant `sv-detection` page, kept distinct from the
sibling split-read [[breakpoint-detection-split-reads]] and read-depth [[read-depth-cnv-segmentation]]
units; enriched with a "Method contract (algorithm spec)" section carrying the `FindDiscordantPairs` /
`ClassifySV` / `DetectSVs` entry points + defaults (μ=400/σ=50/c=3.0σ/clusterDistance=500/minSupport=2),
the 8-tuple read-pair input + `StructuralVariant`(Type/Start/End/Length/SupportingReads/Quality) /
`SVType` output, the first-match classification order with the `ComplexRearrangement` fallback +
`maxInsertSize=10000` guard, the O(n log n) sort-dominated cost, and the ArgumentNullException /
inclusive-`μ±c·σ`-boundary contract. No new page; spec added to `sources:` ahead of the Evidence doc,
`source_commit`→a408154d, `updated`→2026-07-17. No `wiki/sources/` page — a spec, not an Evidence/Validation
report; hub [[algorithm-validation-evidence]] and graph edges unchanged, all already cover SV-DETECT-001);
StructuralVar/Copy_Number_Variation → [[read-depth-cnv-segmentation]] resolved 2026-07-17
(RECONCILE/REUSE: the Copy Number Variation spec, unit SV-CNV-001, describes the read-depth →
windowed-mean → log2 ratio → integer copy-number pipeline already synthesized on this Evidence-derived
concept — reconciled the algorithm spec onto that page as its canonical PRIMARY spec rather than
creating a redundant `copy-number-variation` page. Added a "Method contract (algorithm spec)" section
carrying the genuinely-distinct implementation surface: the `StructuralVariantAnalyzer.DetectCNV(depthData,
windowSize=100, referenceDepth=null, chromosome="chr1")` + `SegmentCopyNumber(logRatios, chromosome="chr1")`
entry points + `LogRatioToCopyNumber`/`OverallMedianNonZero` privates, the `CopyNumberSegment` output
contract (Start/End, LogRatio, CopyNumber=round(2·2^log2)≥0, BAlleleFrequency=NaN, ProbeCount), invariants
INV-01…INV-06, the exception/edge-case contract (null⇒ArgumentNullException, windowSize≤0⇒ArgumentOutOfRangeException,
trailing partial window dropped, zero-depth/NaN⇒no-call, ref≤0⇒no segments), O(n)/O(n·windowSize⁻¹) cost,
and the 100/50/150-depth CN 2/1/3 oracle. Verified the log2/CN arithmetic and window default against the
spec. Kept CNV distinct from the sibling SV_Detection (still pending) and the SV anchor
[[breakpoint-detection-split-reads]]. No new page; spec added to `sources:` ahead of the Evidence doc,
`source_commit`→941abc5b, `updated`→2026-07-17. No `wiki/sources/` page — a spec, not an Evidence/Validation
report; hub [[algorithm-validation-evidence]] and graph edges unchanged, all already cover SV-CNV-001);
Statistics/Sequence_Summary → [[base-composition]] resolved 2026-07-17, **closing the Statistics
domain** (last pending doc) (RECONCILE/REUSE: the top-level `SummarizeNucleotideSequence`
summary-record aggregator, unit SEQ-SUMMARY-001, bundles a sequence's length, composition, GC
content, Shannon entropy, linguistic complexity and melting temperature into one `SequenceSummary`
record by pure field-by-field delegation to already-validated per-metric methods — per its own
Evidence [[seq-summary-001-evidence]] it "adds no new concept," and every aggregated field already
has a home ([[base-composition]] itself for composition/GC, [[shannon-entropy]], [[linguistic-complexity]],
[[melting-temperature]]); reconciled the spec onto the composition page that already carried the
SEQ-SUMMARY-001 Evidence in `sources:` and described it as the SEQ-STATS-001 umbrella's aggregation
wrapper, rather than creating a redundant `sequence-summary` page. Added a "`SummarizeNucleotideSequence`
aggregator (SEQ-SUMMARY-001 primary spec)" section carrying the genuinely-distinct implementation
surface — the `SequenceStatistics.SummarizeNucleotideSequence(string)` entry point + `SequenceSummary`
record, the six field→sub-analyzer map (INV-01…INV-06), the summary's only branch
`useWallaceRule: S.Length < 14` (strict `<`, 14→GC/Marmur-Doty branch), the O(n)/O(n) cost dominated
by linguistic complexity (k≤6), suffix-tree-N/A, the degenerate zero summary / case-insensitive /
U+N-counted edge cases, and the `ATGCATGC` oracle. Verified the field list against the spec. No new
page; spec added to `sources:`, `source_commit`→f3f84d5b. No `wiki/sources/` page — a spec, not an
Evidence/Validation report; hub [[algorithm-validation-evidence]] and graph edges unchanged, all
already cover SEQ-SUMMARY-001);
Statistics/Molecular_Weight_Calculation → [[molecular-weight]] resolved 2026-07-17
(RECONCILE/REUSE: the Molecular Weight Calculation spec, unit SEQ-MW-001, describes the
average-isotopic protein/DNA/RNA mass scalar `MW = Σ monomer masses − (n−1)·water` already
synthesized on this Evidence-derived concept — reconciled the algorithm spec onto that page as its
canonical PRIMARY spec rather than creating a new `molecular-weight-calculation` page (the
`molecular-weight` slug already existed, paired with [[isoelectric-point]] as the ExPASy Compute
pI/Mw partner), and enriched it with a "Method contract (algorithm spec)" section carrying the
genuinely-distinct `SequenceStatistics.CalculateMolecularWeight(string)` (protein) /
`CalculateNucleotideMolecularWeight(string, isDna=true)` (DNA/RNA) entry points + Analysis-assembly
location, the named mass-table constants, the O(n)-time/O(1)-space single-pass cost, the named
invariants INV-01…INV-04, and the single-stranded-average-only scope note (monoisotopic /
double_stranded / circular are Biopython-only, not implemented here). Verified the average mass
tables + water 18.0153 Da against the spec. No new page; spec added to `sources:` ahead of the
Evidence doc, `source_commit`→d4ae46ad. No `wiki/sources/` page — a spec, not an Evidence/Validation
report; hub [[algorithm-validation-evidence]] and graph edges unchanged, all already cover SEQ-MW-001);
(Statistics/Isoelectric_Point → [[isoelectric-point]] resolved 2026-07-17
(RECONCILE/REUSE: the Isoelectric Point spec, unit SEQ-PI-001, describes the pI = net-charge-zero
pH found by bisecting the Henderson–Hasselbalch net-charge function over [0,14] already synthesized
on this Evidence-derived concept — reconciled the algorithm spec onto that page as its canonical
PRIMARY spec rather than creating a new `isoelectric-point` page (the slug already existed), and
enriched it with the genuinely-distinct `SequenceStatistics.CalculateIsoelectricPoint(string)`
public entry point + private `NetCharge` helper + Analysis-assembly location, the O(n)-time/O(1)-space
cost with ≈11 fixed bisection iterations, the named invariants INV-01…INV-04, the not-implemented
PTM/non-standard-residue pKa scope (defer to IPC/pIChemiSt), and the 2-D gel/IEF + ion-exchange
applications. No new page; spec added to `sources:`, `source_commit`→5088660a. No `wiki/sources/`
page — a spec, not an Evidence/Validation report; hub [[algorithm-validation-evidence]] and graph
edges unchanged, all already cover SEQ-PI-001);
(Statistics/Hydrophobicity_Analysis → [[hydrophobicity-gravy-and-profile]] resolved 2026-07-17
(RECONCILE/REUSE: the Hydrophobicity Analysis spec, unit SEQ-HYDRO-001, describes the
Kyte-Doolittle GRAVY scalar + sliding-window hydropathy profile already synthesized on this
Evidence-derived concept — reconciled the algorithm spec onto that page as its canonical PRIMARY
spec rather than creating a new `hydrophobicity-analysis` page, and enriched it with the
genuinely-distinct `SequenceStatistics.CalculateHydrophobicity(string)` /
`CalculateHydrophobicityProfile(string, windowSize=9)` entry points + Analysis-assembly location,
the **default windowSize=9** (surface; 19 for TM), the lazy `IEnumerable<double>` O(n·W)/O(1)
profile contract, and the named invariants INV-01…INV-05. No new page; spec added to `sources:`,
`source_commit`→1163a158. No `wiki/sources/` page — a spec, not an Evidence/Validation report;
hub [[algorithm-validation-evidence]] and graph edges unchanged, all already cover SEQ-HYDRO-001);
(Statistics/GC_Content_Profile → [[windowed-gc-profile-and-variance]] resolved 2026-07-17
(RECONCILE/REUSE: the GC Content Profile spec, unit SEQ-GC-PROFILE-001, describes the standalone
sliding-window GC%-content channel already synthesized as the composite page's standalone entry
point — reconciled the algorithm spec onto that concept as its canonical PRIMARY spec rather than
creating a new page, and enriched the standalone section with the genuinely-distinct
`SequenceStatistics.CalculateGcContentProfile(string, windowSize=100, stepSize=1, fraction=false)`
signature/location, its bare-`IEnumerable<double>` lazy output contract (no coordinates/skew/variance),
the distinct w=100/step=1 defaults vs the composite's 1000/100, the opt-in `fraction` [0,1] flag,
`O(W·windowSize)`/`O(1)`-streaming cost, and the overlap-aware oracle
`CalculateGcContentProfile("GGGAAATGCC",4,3) → [75,0,75]`. No new page; spec added to `sources:`,
`source_commit`→2beeab21. No `wiki/sources/` page — a spec, not an Evidence/Validation report; hub
and graph edges unchanged, all already cover SEQ-GC-PROFILE-001);
(Statistics/Entropy_Profile → [[entropy-profile]] resolved 2026-07-17
(NEW concept: the sliding-window **Shannon entropy profile** SEQ-ENTROPY-PROFILE-001 —
`SequenceStatistics.CalculateEntropyProfile(seq, windowSize=50, stepSize=1)` yields one
`H = −Σ pᵢ log₂ pᵢ` (bits) per window as a lazy `IEnumerable<double>`, the windowed consumer of the
scalar [[shannon-entropy]]'s **general-alphabet** kernel `SequenceStatistics.CalculateShannonEntropy`
(counts ALL letters, so `N`/degenerate count as themselves, protein windows exceed 2 bits, no T↔U
normalization). Created its own page rather than reconciling because it is a **genuinely-distinct
method** from the DNA-canonical [[windowed-sequence-complexity-profile]]
(`SequenceComplexity.CalculateWindowedComplexity`, `ComplexityPoint` = A/T/G/C Shannon **+** linguistic
complexity, W=64/step=10) — the Statistics-domain slug `entropy-profile` was already reserved as
pending on [[shannon-entropy]]. Moved SEQ-ENTROPY-PROFILE-001's home off the windowed page's
"second entry point" section (now a sibling pointer) onto the new page; added reciprocal typed edges
[[entropy-profile]]↔[[shannon-entropy]] / ↔[[windowed-sequence-complexity-profile]] / →[[k-mer-statistics]]
/ →[[test-unit-registry]]; repointed the SEQ-ENTROPY-PROFILE-001 verdict row in [[validation-verdicts]]
(was `?`); added the spec to `sources:` on the windowed + shannon-entropy pages (both now describe the
method). No `wiki/sources/` page — a spec, not an Evidence/Validation report; hub unchanged);
Statistics/Dinucleotide_Analysis → [[dinucleotide-relative-abundance]] resolved 2026-07-17
(RECONCILE/REUSE: the Dinucleotide Analysis spec, unit SEQ-DINUC-001, describes exactly the
already-synthesized dinucleotide-frequency / Karlin odds-ratio / codon-frequency engine — reconciled
onto the existing Evidence-derived concept as its canonical PRIMARY spec rather than creating the
placeholder `dinucleotide-analysis` slug; added a "Method contract (algorithm spec)" section carrying
the genuinely-distinct implementation detail — the three `SequenceStatistics` signatures
(`CalculateDinucleotideFrequencies` / `CalculateDinucleotideRatios` /
`CalculateCodonFrequencies(seq, readingFrame=0)`) + `IReadOnlyDictionary<string,double>` return, the
`{A,T,G,C,U}` dinuc / `{A,T,G,C}` codon alphabet filters, `O(n)`/`O(k≤25|64)` tabulation (suffix tree
N/A), the INV-01/INV-02 frequency-sum-to-1 invariants, and the single-strand-vs-strand-symmetrized-`ρ*`
deviation + not-implemented `0.78`/`1.23` classification; kept distinct from the CG special case
[[cpg-island-detection]] and the k-mer-diversity view [[k-mer-statistics]]; no `wiki/sources/` page —
a spec, not an Evidence report; hub unchanged);
Statistics/DNA_Thermodynamics → [[dna-duplex-nearest-neighbor-thermodynamics]] resolved 2026-07-17
(RECONCILE/REUSE: the nearest-neighbour DNA-duplex thermodynamics spec, unit SEQ-THERMO-001, describes
exactly the already-synthesized full-tuple ΔH°/ΔS°/ΔG°/Tm engine — reconciled onto the existing
Evidence-derived concept as its canonical PRIMARY spec rather than creating the placeholder
`dna-thermodynamics` slug; added a "Method contract (algorithm spec)" section carrying the
genuinely-distinct implementation detail — the `CalculateThermodynamics(seq, naConcentration=0.05,
primerConcentration=2.5e-7)` signature + defaults (C_T ÷ F = 4), the 4-tuple output rounding
(ΔH/ΔS/ΔG 2 dp, Tm 1 dp), `O(n)`/`O(1)` tabulation with both WC complements stored (suffix tree N/A),
the length-<2 ⇒ `(0,0,0,0)` / case-insensitive / non-ACGT-⇒0 guards, the fixed-F=4 simplification and
not-implemented mismatch/dangling-end/Mg²⁺ corrections, and the one-terminus→two-terminus initiation
fix; kept distinct from the `%GC` scalar Tm on [[melting-temperature]] and the 2004-parameter primer
engine [[primer-dimer-thermodynamics-tm]]; no `wiki/sources/` page — a spec, not an Evidence report;
hub unchanged);
Splicing/Gene_Structure_Prediction → [[gene-structure-prediction-intron-exon]] resolved 2026-07-17,
**closing the Splicing domain** (last pending doc) (REUSE: the composite intron/exon gene-structure
predictor spec reconciled onto the existing splicing-family composite page rather than creating the
placeholder `gene-structure-prediction` slug; enriched with a "Method contract (algorithm spec)"
section carrying the genuinely-distinct implementation content — the `PredictGeneStructure(seq,
minExonLength=30, minIntronLength=60, minScore=0.5)` + `PredictIntrons(seq, minIntronLength=60,
maxIntronLength=100000, minScore=0.5)` signatures (the latter's `maxIntronLength` **fixed to 100000**
inside `PredictGeneStructure`), the `FindDonorSites`/`FindAcceptorSites` **`minScore*0.8`** site
threshold + `includeNonCanonical=true` for both, intron length `acceptor.Position-donor.Position+1`,
the branch-point search window **`[acceptor.Position-50, acceptor.Position-18]`** at min branch score
**0.4**, combined score `(donor+acceptor+branch)/3` with branch else `(donor+acceptor)/2`, greedy
descending-score non-overlap selection, the `GeneStructure`/`Exon`/`Intron` output records, invariants
INV-01…INV-04, and the **SplicedSequence-vs-exon-record coordinate caveat** (spliced sequence is
defined by intron removal, `minExonLength` filters only the `Exons` list, so a short gap can be
dropped from `Exons` yet remain spliced-in); no `wiki/sources/` page — a spec, not an
Evidence/Validation report; hub unchanged). This closes Splicing (Donor + Acceptor + Gene_Structure
all resolved 2026-07-17);
Splicing/Donor_Site_Detection → [[splice-donor-site-prediction]] resolved 2026-07-17
(REUSE: the algorithm spec describes exactly the already-synthesized SPLICE-DONOR-001 5' donor
unit — reconciled onto the existing splicing-family donor page rather than creating the
placeholder `donor-site-detection` slug; enriched with a "Method contract (algorithm spec)"
section — `FindDonorSites(seq, minScore=0.5, includeNonCanonical=false)` + opt-in
`ScoreDonorMaxEnt`, the scan starting at **index 0** (vs the acceptor's 15), the `SpliceSite`
output contract (`Position = i` = the donor dinucleotide **start**, `Type` `Donor` for GU **and**
GC / `U12Donor` for AU, `Motif` via `GetMotifContext(seq,i,3,6)`, `Confidence`=clamp((score−0.5)/0.5)),
the `ScoreDonorSite` binary `DonorPwm` −3..+5 match-fraction internals + `ScoreU12DonorSite`
`AUAUCC`/6 scorer, and invariants INV-01…INV-04; **corrected** the guard claim from "shorter than
the 9-nt window → empty" to the actual **length < 6** guard (the default scorer normalizes over
in-bounds positions and does not need a full 9-nt window; only opt-in `ScoreDonorMaxEnt` requires
exactly 9 nt); no `wiki/sources/` page — a spec, not an Evidence/Validation report). Gene_Structure
is the last doc pending in Splicing);
Splicing/Acceptor_Site_Detection → [[splice-acceptor-site-prediction]] resolved 2026-07-17
(REUSE: the algorithm spec describes exactly the already-synthesized SPLICE-ACCEPTOR-001 unit —
reconciled onto the existing splicing-family anchor rather than creating the placeholder
`acceptor-site-detection` slug; enriched with the `SpliceSitePredictor` method signatures +
defaults, the `SpliceSite` output contract (`Position = i+1` = the `G`, `Type`
`Acceptor`/`U12Acceptor`, `Motif`/`Score`/`Confidence`), the canonical PPT window
`[position-15,position-3)`/12·2 + sparse `AcceptorPwm` offsets `-15,-10,-5,-4,-3,-2,-1,0`, the U12
`3.5` normalizer, and invariants INV-01…INV-07; **corrected** the branch-point report threshold
from `≥ 0.8` to `≥ minScore` **default 0.5** (spec §2.2 + `FindAcceptorBranchPoint` impl); no
`wiki/sources/` page — a spec, not an Evidence/Validation report). Gene_Structure remains
pending in Splicing);
Sequence_Composition/Shannon_Entropy → [[shannon-entropy]] resolved 2026-07-17, **closing the
Sequence_Composition domain** (last pending doc) (NEW concept: the base per-symbol Shannon entropy
`H = −Σ p·log₂ p` bits — a genuinely-distinct scalar member of the `SEQ-COMPLEX-*` complexity/entropy
family that had no home. Two entry points with deliberately different alphabets: canonical DNA
`SequenceComplexity.CalculateShannonEntropy` (only A/T/G/C, range [0,2]) and general-alphabet
`SequenceStatistics.CalculateShannonEntropy` (all letters, the field bundled by `SummarizeNucleotideSequence`);
kept DISTINCT from the higher-order k-mer k-entropy of [[k-mer-statistics]] (`CalculateKmerEntropy`,
SEQ-COMPLEX-KMER-001 — per-base entropy is its k=1 composition-only case) and from the per-window
Shannon channel of [[windowed-sequence-complexity-profile]] / the still-pending Statistics-domain
`entropy-profile` consumer; mapped the MCP pair `complexity_shannon`/`shannon_entropy` onto it and
repointed the SEQ-ENTROPY-001 verdict row in [[validation-verdicts]] (was `?`); no `wiki/sources/`
page — a spec, not an Evidence/Validation report);
Sequence_Composition/Sequence_Validation → [[sequence-validation]] resolved 2026-07-17 (NEW concept:
the standalone SEQ-VALID-001 strict DNA/RNA alphabet-membership unit — `SequenceExtensions.IsValidDna`/`IsValidRna(ReadOnlySpan<char>)`
+ `DnaSequence.TryCreate(string, out DnaSequence?)` — had no home; it is genuinely distinct from
[[fasta-parsing]]'s opt-in multi-alphabet `SequenceAlphabet` parse-time validator (StrictDna/IupacNucleotide/Rna/Protein)
and from the ambiguity-*interpreting* [[iupac-degenerate-matching]]; created as the strict go/no-go
character-membership primitive whose tally counterpart is [[base-composition]]; mapped the MCP predicate
pair `is_valid_dna`/`is_valid_rna` onto it and repointed the SEQ-VALID-001 verdict row in [[validation-verdicts]]
(was `?`); no `wiki/sources/` page — a spec, not an Evidence/Validation report);
Sequence_Composition/Replication_Origin_Prediction → [[replication-origin-cumulative-skew]] resolved
2026-07-17 (REUSE: the SEQ-REPLICATION-001 primary spec reconciles onto the existing cumulative-GC-skew
origin/terminus concept — enriched with the `PredictReplicationOrigin(DnaSequence|string)` entry points,
the `ReplicationOriginPrediction` record-struct output contract, INV-01…06, O(n)/O(1) single-pass
complexity, and the null/empty-vs-`ArgumentNullException` overload edge cases);
Sequence_Composition/RNA_Complement → [[rna-base-pairing]] resolved 2026-07-17 (REUSE: the
SEQ-RNACOMP-001 per-base IUPAC-complete RNA complement `GetRnaComplementBase` is already synthesized
in the RNA-base-pairing concept's "SEQ-family full-IUPAC RNA complement" section — enriched with the
spec's entry-point signature/location, O(1) `switch` complexity, X pass-through correction, and the
`"ACGTUacgtuXYZxyz"` worked oracle);
Sequence_Composition/Linguistic_Complexity → [[linguistic-complexity]] resolved 2026-07-17 (NEW
concept: the standalone scalar linguistic-complexity unit SEQ-COMPLEX-001 = vocabulary-usage ratio
`Σ Vᵢ / Σ Vmax,i`, previously only surfaced per-window inside [[windowed-sequence-complexity-profile]]);
RnaStructure/Turner_McCaskill_Partition_Function → [[turner-mccaskill-partition-function]] resolved
2026-07-17, **closing the RnaStructure domain** (last pending doc) (NEW concept: the full Turner-2004
nearest-neighbour McCaskill partition function is a genuinely-distinct engine — not the base-pair-counting
teaching model of [[rna-partition-function-mccaskill]] (RNA-PARTITION-001), which explicitly listed a
*full Turner-parameter partition function* as **not implemented**. This unit reuses the MFE folder's
Turner-2004 loop energies under Boltzmann-weighted inside recurrences (`Vexp`/`WMexp`/`Wexp`, `Z=Wexp(0,n−1)`),
exposes distinct entry points `CalculateUnpairedProbabilities`/`CalculateRegionUnpairedProbability` and
distinct outputs (`p_unpaired(i)=1−Σ_jP(i,j)`, `ΔG_ensemble=−RT·ln Z`, RNAplfold-style region accessibility
`Z_open/Z`), computes marginals by **constrained re-folds** (`Z_forbid(i)/Z`, `Z_require(i,j)/Z` — no outside
recursion, `O(n⁵)` bpp / `O(n³)` Z), and drives the TargetScan context++ **SA** feature consumed by
[[mirna-target-site-prediction]]; the simplified page was cross-linked to it and its `source_commit` bumped;
no `wiki/sources/` page — a spec, not an Evidence/Validation report);
RnaStructure/RNA_Free_Energy → [[rna-free-energy-turner-model]] resolved 2026-07-16 (REUSE: the
RNA-ENERGY-001 aggregate free-energy spec reconciles onto the existing Turner-model energy-terms
concept — titled for RNA-ENERGY-001 and already listing its Evidence in `sources:` — enriched only
with the aggregate spec's genuinely-distinct implementation surface: the Simplified status of the
whole energy layer, partial Turner-table coverage (1×1 internal loops exact; int21/int22 NOT
implemented; 2×3+ generic mismatch fallback), the multibranch free-base cost fixed at 0.0, the
`MAXLOOP = 30` bound, the retained internal `CalculateMinimumFreeEnergyClassic` baseline, and the
three-helper O(b)/O(l)/O(n³) complexity table; the whole-structure `CalculateMinimumFreeEnergy`
summation contract stays on [[rna-minimum-free-energy-folding]], not duplicated);
Repeat_Analysis/Tandem_Repeat_Detection → [[repetitive-element-detection]] resolved 2026-07-16,
**closing the Repeat_Analysis domain** (last pending doc) (REUSE: the REP-TANDEM-001 spec reconciles
against the existing repeats/tandem family anchor, which already synthesizes `GenomicAnalyzer.FindTandemRepeats`
as the consolidated GENOMIC-TANDEM-001 duplicate — the exact head-to-tail detector, its period-ambiguous
non-canonicalizing behaviour, and the 1–6 bp microsatellite class; the spec adds no new detector, so no
contradiction; enriched only with the spec's genuinely-distinct implementation surface: the `RepeatFinder.GetTandemRepeatSummary(DnaSequence, int minRepeats=3)`
aggregation helper that returns a single `TandemRepeatSummary` record by delegating to `FindMicrosatellites(sequence, 1, 6, minRepeats)`
(inheriting its null-throw and `minRepeats ≥ 2` floor), rolling 1–6 bp microsatellites into total count / total bases /
percent coverage / longest repeat / most-frequent unit / per-class counts, and the two scope caveats — the summary
sees only 1–6 bp units (minisatellite/macrosatellite tandems excluded) and its named per-class fields stop at
tetranucleotide (penta/hexa feed totals but get no dedicated field); no `wiki/sources/` page — a spec, not an
Evidence/Validation report);
Repeat_Analysis/Microsatellite_Detection → [[repetitive-element-detection]] resolved 2026-07-16
(REUSE: the REP-STR-001 spec reconciles against the existing repeats/tandem family anchor, which
already synthesizes the REP-STR-001 Evidence — the Benson TRF approximate detector, the
`ComputeBernoulliStatistics` PM/PI layer, and the 1–6 bp microsatellite class; enriched only with
the spec's genuinely-distinct implementation surface for the *perfect* default detector: the four
`RepeatFinder.FindMicrosatellites` overloads (DnaSequence/string × cancellable+progress), the
`minUnitLength`/`maxUnitLength`/`minRepeats` = 1/6/3 defaults and validation floors, the
`MicrosatelliteResult`/`RepeatType` mono–hexa output, the `IsRedundantUnit` primitive-unit filter
and the contained-interval suppression that is narrower-than-non-overlap (§5.4 Deviation 1), the
`SequenceAligner.GlobalAlign` reuse, and the `O(n·U·R)` perfect / `O(n²·P·L²)` approximate cost;
no `wiki/sources/` page — a spec, not an Evidence/Validation report);
Repeat_Analysis/Direct_Repeat_Detection → [[direct-repeat-detection]] resolved 2026-07-16 (NEW
concept: no Evidence-derived concept existed for this unit — `RepeatFinder.FindDirectRepeats`,
REP-DIRECT-001 — and it is a genuinely-distinct same-orientation dispersed-pair operation, not the
head-to-tail tandem / reverse-complement inverted sub-problems of the [[repetitive-element-detection]]
anchor nor the position-list [[longest-repeated-substring]] enumerator; created with an inbound link
from the family anchor);
Repeat_Analysis/Inverted_Repeat_Detection → [[inverted-repeat-detection]] resolved 2026-07-16 (NEW
concept: the reverse-complement sibling of [[direct-repeat-detection]] — `RepeatFinder.FindInvertedRepeats`,
REP-INV-001 — a *third*, genuinely-distinct DNA-only exact-arm implementation, separate from the
imperfect IUPACpal `W·G·W̄ᴿ` annotation model and the RNA-INVERT-001 stem model both synthesized in the
[[repetitive-element-detection]] anchor; adds the `CanFormHairpin` (loop ≥ 3) flag and the palindrome =
zero-loop special case; created with inbound links from the family anchor and the direct-repeat sibling);
Repeat_Analysis/Palindrome_Detection → [[palindrome-detection]] resolved 2026-07-16 (NEW concept:
the biological DNA palindrome — even-length windows where `S = ReverseComplement(S)` — is exposed
via a genuinely-distinct dedicated entry point `RepeatFinder.FindPalindromes` (+ the lighter
`GenomicAnalyzer.FindPalindromes`), REP-PALIN-001, with its own `PalindromeResult`/`PalindromeInfo`
return types and even-length step-by-2 scan, separate from the loop-bearing `FindInvertedRepeats`
of [[inverted-repeat-detection]] which will not report the zero-loop case under default
`minLoopLength=3`; created with inbound links from the family anchor [[repetitive-element-detection]]
and the [[inverted-repeat-detection]] sibling, cross-links [[restriction-site-detection]]);
Quality/Quality_Statistics → [[fastq-quality-statistics]] resolved 2026-07-16, closing the Quality
domain (last pending doc) (REUSE: reconciled against the existing Evidence-derived concept, which
already synthesized the mean/median/min/max, population variance/σ (÷N), inclusive Q20/Q30 thresholds,
the Q30 NGS benchmark, the arithmetic-vs-probability mean distinction, and the `5?I`/`5II?`/single-`I`
oracles; enriched only with the spec's genuinely-distinct implementation content: the three
`QualityScoreAnalyzer` entry points, the multi-read `PerPositionMeanQuality` delegate variant, the
O(n log n) sort-dominated / O(n) Q30 complexity, and the sibling `CalculateExpectedErrors` /
`PhredToErrorProbability` error-probability summaries);
Quality/Phred_Score_Handling → [[phred-quality-encoding]] resolved 2026-07-16 (REUSE: reconciled
against the existing Evidence-derived concept, which already synthesized the Q=−10log₁₀p definition,
the Phred+33/Phred+64 offsets and ranges, auto-detection, the score-preserving ±31 re-offset and its
overflow rules; enriched only with the spec's genuinely-distinct implementation content — the three
`QualityScoreAnalyzer` canonical methods (`ParseQualityString`/`ToQualityString`/`ConvertEncoding`) in
Seqeron.Genomics.IO, the O(n) single-pass cost, the `Auto`→`DetectEncoding` (parse) vs `Auto`→Phred+33
(encode) resolution, the INV-03 round-trip invariant, the ArgumentNullException/ArgumentOutOfRangeException
contract, the suffix-tree-not-applicable note, and the legacy `QualityStringToPhred`/`PhredToQualityString`
helpers that lack range validation));
ProteinPred/MoRF_Prediction → [[morf-prediction-dip-in-disorder]] resolved 2026-07-16, closing the
ProteinPred domain (last pending doc) (REUSE: reconciled against the existing Evidence-derived
concept, which already documents the dip-in-disorder criterion, the 10–70 Mohan band, the TOP-IDP
disorder score, α/β/ι sub-types and the paywalled-flank assumption; enriched only with the spec's
genuinely-distinct implementation content: the `PredictMoRFs(string, int, int)` entry point in
`DisorderPredictor.cs`, the O(n·w)/O(n) cost with the 21-residue window, the 0-based-inclusive /
non-overlapping / start-ordered output contract, the explicit `(0.5 − mean d)/0.5` clamped score
formula, the 0.5-MoRF-vs-0.542-TOP-IDP dual-threshold distinction, the suffix-tree-not-applicable
decision and the P/L/I worked oracle (Start 29–End 50, Score ≈ 0.2759 → 0.3996));
ProteinMotif/Signal_Peptide_Prediction → [[protein-domain-and-signal-peptide-prediction]] resolved
2026-07-16 (REUSE: the dedicated von Heijne 1986 signal-peptide spec, unit PROTMOTIF-SP-001, is already
fully owned by this Evidence-derived concept — which documents the current EMBOSS `sigcleave` weight-matrix
method with the P17644 oracle; the spec describes the SAME correct method (no superseded tripartite model —
that was the sibling `Domain_Prediction.md`), so no contradiction; enriched only with the spec's
genuinely-distinct implementation content: the `PredictSignalPeptide`/`BuildWeightMatrix` entry points with
the log-odds matrix built once at static init, the O(n·15)=O(n) time / O(1) space cost, the fixed-width-PWM
suffix-tree-not-applicable decision, and the min-15-residue-window accepted deviation);
ProteinMotif/Profile_HMM_Domain_Detection → [[protein-domain-and-signal-peptide-prediction]] resolved
2026-07-16 (REUSE: the per-algorithm profile-HMM spec is the SAME test unit — PROTMOTIF-DOMAIN-001 — as
`Domain_Prediction`, already owned by this Evidence-derived concept, which had a dedicated `FindDomainsByHmm`
Plan7 section covering Viterbi/Forward, hmmsearch-parity `pre_score`, null2, the Gumbel/exponential E-value
layer, multi-domain envelope decomposition and stochastic-traceback clustering; enriched only with the
spec's genuinely-distinct implementation content: the HMMER3/f parser layout (−ln-to-5-decimals storage,
`*`→−∞, COMPO/BEGIN, the 7-transition node order), the two-row O(n·M)/O(M) Viterbi/Forward DP shape, the
glocal-default-vs-local-multihit-opt-in distinction (spec Deviation #1) with `minBitScore`/`Z` defaults and
the `FindDomainHitsByHmm`/`ScoreDomainHmmEValue`/`FindDomainEnvelopes` method surface, and the INV-HMM-01/03/05 invariants);
ProteinMotif/Pattern_Matching_Methods → [[protein-motif-pattern-search]] resolved 2026-07-16 (REUSE:
this Evidence-derived concept already cites the PROTMOTIF-PATTERN-001 Evidence and is ahead of the spec;
added the spec to `sources:` and enriched only its genuinely-distinct implementation content — `FindDomains`
as the fourth pattern-matching primitive delegating to the same lookahead+IC engine and wrapping hits as
`ProteinDomain`, plus the suffix-tree-evaluated-not-used matcher decision);
ProteinMotif/PROSITE_Pattern_Matching → [[protein-motif-pattern-search]] resolved 2026-07-16 (REUSE:
the dedicated PROSITE spec, unit PROTMOTIF-PROSITE-001, is already fully covered by this Evidence-derived
concept — which is ahead of the spec; enriched only with the spec's explicit scope boundary: syntax
converter + regex-search wrapper with no PROSITE profile/matrix scanning, no external-catalog lookup, no
ScanProsite result metadata, and repository-defined `Score`/`EValue`);
ProteinMotif/Motif_Search → [[protein-motif-pattern-search]] resolved 2026-07-16 (REUSE: the primary
per-algorithm spec is exactly this Evidence-derived concept, unit PROTMOTIF-FIND-001; enriched only with
the `ParseRegexAllowedCounts` regex-walk helper that supplies the per-position allowed-residue count
feeding the information-content score);
ProteinMotif/Domain_Prediction → [[protein-domain-and-signal-peptide-prediction]] resolved 2026-07-16
(REUSE: the primary per-algorithm spec reconciles against the existing Evidence-derived concept, which is
already ahead of the spec — it documents the current von Heijne 1986 weight-matrix signal-peptide method
and the Plan7 profile-HMM engine, whereas the spec still describes the superseded/fabricated tripartite
n/h/c model; enriched with the `FindMotifByPattern` score-delegation / information-content-score provenance
and the O(n·d), d=3 domain-scan cost);
ProteinMotif/Common_Motif_Finding → [[common-protein-motifs]] resolved 2026-07-16 (REUSE: the primary
per-algorithm spec is exactly this Evidence-derived concept; enriched with the O(p·n) `FindCommonMotifs`
complexity, the `FindMotifByPattern` delegation + entry-point decomposition, the suffix-tree-not-applicable
design decision, and the `FindAllKnownMotifs` registry-naming deviation);
ProteinMotif/Coiled_Coil_Prediction → [[coiled-coil-prediction]] resolved 2026-07-16 (REUSE: the primary
per-algorithm spec is exactly this Evidence-derived concept; enriched with the `BestHeptadOccupancy`/`BuildRegion`
helper decomposition, the single-pass array precompute + lazy-yield shape, and the suffix-tree-not-applicable note);
Population_Genetics/Runs_Of_Homozygosity → [[runs-of-homozygosity-inbreeding]] resolved 2026-07-16,
closing the Population_Genetics domain (REUSE: reconciled against the existing Evidence-derived concept;
enriched with the O(n log n)/O(n) `FindROH` cost, the O(m)/O(1) F_ROH cost, the inclusive-run vs
half-open-segment coordinate distinction, and the eager argument-validation / `genomeLength ≤ 0 → 0`
contract);
K-mer, Metagenomics, MolTools, Oncology, PanGenome and Pattern_Matching domains now fully covered;
Phylogenetics/Tree_Construction → [[distance-based-tree-construction]] resolved 2026-07-15 (REUSE: the primary UPGMA/NJ spec is exactly this Evidence-derived hinge concept; enriched with the `distanceMethod=JukesCantor` default, O(n²) space, the cluster-index dictionary detail, and the §5.4 accepted-deviation / §6.2 builder limitations);
Phylogenetics/Newick_Format → [[distance-based-tree-construction]] resolved 2026-07-15 (REUSE: the Newick I/O layer serializes this hinge concept's `PhyloNode` output — a format serializer, not a separate algorithm, so no dedicated concept, per [[phylo-newick-001-evidence]]);
Phylogenetics/Distance_Matrix → [[evolutionary-distance-matrix]] resolved 2026-07-15 (reconciled against the existing Evidence-derived concept);
Phylogenetics/Tree_Comparison → [[tree-comparison-metrics]] resolved 2026-07-15 (reconciled against the existing Evidence-derived concept; enriched with the O(n² log n) impl cost, helper methods, and test-pinned oracles);
Phylogenetics/Tree_Statistics → [[tree-statistics]] resolved 2026-07-15, closing the Phylogenetics domain (REUSE: reconciled against the existing Evidence-derived concept; enriched with the `PhylogeneticAnalyzer` method set, the O(n)/O(h) recurrences, the `yield` pre-order leaf order, the `EmptyTreeHeight` constant, and the topological-vs-Biopython-`depths()` distinction);
Pattern_Matching/Suffix_Tree → [[suffix-tree]] resolved 2026-07-15, closing the Pattern_Matching domain (last pending doc);
PanGenome/Phylogenetic_Marker_Selection → [[phylogenetic-marker-selection]] resolved 2026-07-15, closing the PanGenome domain;
Oncology/Variant_Allele_Frequency → [[variant-allele-frequency-and-binomial-ci]] resolved 2026-07-15, closing the Oncology domain;
Oncology/Tumor_Purity_Estimation → [[tumor-purity-from-mutation-vaf]] resolved 2026-07-15;
Oncology/Tumor_Ploidy_Estimation → [[tumor-ploidy-estimation-and-whole-genome-doubling]] resolved 2026-07-15;
Oncology/Tumor_Phylogeny_Reconstruction → [[tumor-phylogeny-clonal-tree-reconstruction]] resolved 2026-07-15;
Oncology/Tumor_Mutational_Burden → [[tumor-mutational-burden]] resolved 2026-07-15;
Oncology/Tumor_Heterogeneity_Analysis → [[intratumor-heterogeneity-metrics]] resolved 2026-07-15;
Oncology/Tumor_Gene_Expression_Outlier → [[expression-outlier-zscore-signature-score]] resolved 2026-07-15;
Oncology/Sequencing_Artifact_Detection → [[sequencing-artifact-detection]] resolved 2026-07-15;
Oncology/SBS96_Trinucleotide_Context_Catalog → [[sbs96-mutational-signature-catalog]] resolved 2026-07-14;
Oncology/Neoantigen_Peptide_Generation → [[neoantigen-peptide-generation]] resolved 2026-07-14;
Oncology/Mutational_Signature_Fitting → [[mutational-signature-fitting-and-extraction]] resolved 2026-07-14;
Oncology/Mutational_Signature_Extraction_NMF → [[mutational-signature-fitting-and-extraction]] resolved 2026-07-14;
Oncology/Mutational_Signature_Exposure_Bootstrap → [[signature-exposure-bootstrap-confidence-intervals]] resolved 2026-07-14;
Oncology/Mutational_Process_Classification → [[mutational-process-classification]] resolved 2026-07-14;
Oncology/Microsatellite_Instability_Detection → [[microsatellite-instability-detection]] resolved 2026-07-14;
Oncology/MRD_Detection → [[tumor-informed-mrd-detection]] resolved 2026-07-14;
Oncology/MHC_Peptide_Binding_Classification → [[mhc-peptide-binding-prediction]] resolved 2026-07-14;
Oncology/Loss_Of_Heterozygosity → [[loss-of-heterozygosity-detection]] resolved 2026-07-14;
Oncology/Known_Fusion_Database_Lookup → [[gene-fusion-nomenclature-known-fusion-lookup]] resolved 2026-07-14;
Oncology/Immune_Infiltration_Estimation → [[immune-infiltration-deconvolution]] resolved 2026-07-14;
Oncology/Homozygous_Deletion_Detection → [[homozygous-deletion-detection]] resolved 2026-07-14;
Oncology/HRD_Score → [[homologous-recombination-deficiency-score]] resolved 2026-07-14;
Oncology/HLA_Nomenclature_And_Allele_Specific_LOH → [[hla-nomenclature-and-allele-specific-loh]] resolved 2026-07-14;
Oncology/Fusion_Gene_Detection → [[gene-fusion-detection-read-evidence]] resolved 2026-07-14;
Oncology/Fusion_Breakpoint_Analysis → [[fusion-breakpoint-frame-and-protein-prediction]] resolved 2026-07-14;
Oncology/Focal_Amplification_Detection → [[focal-amplification-detection]] resolved 2026-07-14;
Oncology/Driver_Mutation_Detection → [[driver-gene-classification-20-20-rule]] resolved 2026-07-14;
Oncology/CtDNA_Analysis → [[ctdna-detection-and-tumor-fraction]] resolved 2026-07-14;
Oncology/Copy_Number_Alteration_Classification → [[copy-number-alteration-classification]] resolved 2026-07-14;
Oncology/Complex_Rearrangement_Classification → [[chromothripsis-inference]] resolved 2026-07-14;
Oncology/Clonal_Hematopoiesis_Filtering → [[clonal-hematopoiesis-cfdna-filtering]] resolved 2026-07-14;
Oncology/Clinical_Actionability_Assessment → [[clinical-actionability-oncokb-levels]] resolved 2026-07-14;
Oncology/Cancer_Variant_Annotation → [[cancer-variant-tier-classification-amp-asco-cap]] resolved 2026-07-14;
Oncology/Cancer_Cell_Fraction_Estimation → [[cancer-cell-fraction-clonal-clustering]] resolved 2026-07-14;
Oncology/Clonal_Subclonal_Classification → [[clonal-subclonal-classification-ccf-posterior]] resolved 2026-07-14; MolTools/Restriction_Site_Detection →
[[restriction-site-detection]] resolved 2026-07-14, closing the MolTools domain; K-mer_Search and PanGenome_Core_Accessory
resolved 2026-07-13; DNA_Dimer_Tm, DNA_Hairpin_Folding_Tm, DNA_Hairpin_Special_Loop_Bonus,
LNA_Adjusted_Nearest_Neighbor_Tm and NearestNeighbor_Salt_Corrected_Tm →
[[primer-dimer-thermodynamics-tm]] resolved 2026-07-13;
Hybridization_Probe_Design → [[hybridization-probe-design]] resolved 2026-07-13;
MolTools/Melting_Temperature → [[melting-temperature]] resolved 2026-07-13;
Primer3_Penalty_Objective → [[primer3-weighted-penalty-objective]] resolved 2026-07-13;
Primer_Design → [[primer-design]] resolved 2026-07-13;
Primer_Structure_Analysis → [[primer-structure-qc-screens]] resolved 2026-07-13;
Probe_Validation → [[probe-offtarget-specificity-scan]] resolved 2026-07-13).

## Covered via concept (done) — moved

The **Covered via concept (done)** table (227 algorithm docs → concept pages) now
lives in **[[backlog-covered-table]]** (split out 2026-07-20 to keep this
append-only history page under the size cap).
