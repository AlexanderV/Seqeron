# Seqeron.Genomics — Validated Limitations & Operating Envelope

**Library:** Seqeron.Genomics (mission-critical)   **Last reviewed:** 2026-06-24

This file lists the library's **current, sourced limitations** — every place where an algorithm is a
faithful but *simplified / subset* realisation of a fuller published method, consumes a caller-supplied
input rather than computing it upstream, is a deliberate from-first-principles **heuristic**, exposes an
**opt-in** alternative convention, or has a documented **scope / API** boundary. Each is `BY-DESIGN` and
**validated correct for its stated contract** (every unit below is `✅ CLEAN`; these are `PASS-WITH-NOTES`
description notes, never defects). The purpose is to make the honest operating envelope visible in one place.

Every row traces to its per-unit report under `docs/Validation/reports/{UNIT}.md`.

---

## 1. Algorithmic simplification (a fuller published method exists; the implemented one is correct for its scope)

| Unit | Not implemented / simplified | Note |
|------|------------------------------|------|
| RNA-STRUCT-001 | Pseudoknot classes outside the pknotsRG canonical csr-PK grammar (kissing hairpins, triple-crossing / chained "complex" helix interactions, non-canonical bulged or unequal-length helices); tertiary-stabilised knots as the MFE structure | Classes Reeder & Giegerich (2004) explicitly exclude from csr-PK; tertiary-stabilised knots (e.g. BWYV / PDB 437D) are not recoverable by *any* nearest-neighbour thermodynamic model — an energy-model floor, not an algorithm gap. |
| PRIMER-TM-001 | Mismatch / dangling-end / secondary-structure NN terms | An opt-in **nearest-neighbour salt-corrected design Tm** (`CalculateMeltingTemperatureNN`) is now provided: SantaLucia (1998) unified ΔH°/ΔS° + bimolecular Tm, with SantaLucia Eq. 5 and Owczarzy (2004) monovalent and Owczarzy (2008) divalent Mg²⁺ corrections; the published 35.8 °C worked example reproduces exactly. The default Wallace/Marmur-Doty Tm is unchanged. **Residual:** the Watson-Crick NN model does not include mismatch / dangling-end terms or hairpin/secondary-structure Tm (use Biopython `Tm_NN`, MELTING 5, or UNAFold for those). |
| REP-STR-001 | Full TRF probabilistic k-tuple seeding + statistical-significance scoring | Approximate (imperfect / interrupted) tandem-repeat detection is now provided opt-in (`FindApproximateTandemRepeats`): each candidate period is aligned against tandem copies of the majority-rule consensus with the Tandem Repeats Finder recommended scoring (match +2, mismatch −7, indel −7; report threshold Minscore = 50), reporting period size, copy number, percent matches, percent indels, consensus, and alignment score — Benson (1999). The default perfect-repeat `FindMicrosatellites` (1–6 bp exact STRs) is unchanged. **Residual:** TRF's probabilistic k-tuple distance-list **seeding** and the **sum-of-Bernoulli statistical-significance** scoring of detected repeats are not bundled — candidate discovery is a deterministic exhaustive (start, period) scan, not whole-genome scale (use the reference TRF tool for that). |
| SPLICE-ACCEPTOR-001 | Full MaxEntScan maximum-entropy 3' model | Explicit branch-point detection is now provided opt-in (`FindAcceptorBranchPoint`: human `yUnAy` consensus, branch A at position 0, scanned 18–40 nt upstream of the AG, reporting branch-point position + conservation-weighted score + polypyrimidine-tract fraction) — Gao et al. (2008) / Mercer et al. (2015). The default `FindAcceptorSites` PWM-consensus + polypyrimidine-tract scoring is unchanged. **Residual:** the Yeo & Burge (2004) MaxEntScan 23-nt maximum-entropy 3' acceptor model is **not bundled** — it needs the Burge-lab precomputed score tables (`me2x3acc*`), which are large trained data files with no clean redistribution licence (use the MaxEntScan reference tool for that). |
| ONCO-IMMUNE-001 | Full LM22 deconvolution + cohort normalisation | Uses a 5-marker immune/stromal signature matrix (vs CIBERSORT LM22) and single-sample un-normalised ssGSEA, so the ESTIMATE-style purity is a **relative**, not clinically-absolute, value. |
| ONCO-MRD-001 | In-pipeline (automatic) background re-estimation | INVAR2 formulas reproduced verbatim. A KDE-smoothed fragment-size weighting is now provided opt-in (`FragmentSizeProfile.FromKernelDensity`: Gaussian KDE, Silverman-rule bandwidth + `adjust`, analytic Gaussian-bin integral — Silverman 1986 / INVAR2 `estimate_real_length_probability`); the default discrete `COUNT/TOTAL` profile is unchanged. **Residual:** the per-locus background `e` can be caller-supplied; `EstimateLocusBackground` re-estimates it from control plasma as a separate (not auto-wired) step, matching INVAR2's own parse/detection separation. |
| META-BIN-001 | Markov-corrected TETRA + single-copy-marker QC | Tetranucleotide signature uses raw-frequency Pearson distance (vs z-score normalisation); bin completeness/contamination use length-ratio & GC-stddev proxies instead of single-copy marker genes (CheckM). |

## 2. Threshold / aggregation layer (classifies a caller-supplied input; does not predict it upstream)

| Unit | Not modelled / caller-supplied |
|------|--------------------------------|
| ONCO-MHC-001 | The peptide–MHC affinity / %Rank prediction (NetMHCpan-style learned model) is caller-supplied: the library classifies a supplied IC50 / %Rank into strong/weak binder but does not bundle a trained predictor (no redistributable, cross-verifiable model available). |

## 3. Declared from-first-principles heuristics (a deliberate, sourced heuristic — not a trained predictor / reference tool)

| Unit | Heuristic in place of |
|------|-----------------------|
| CHROM-CENT-001 | Alpha-satellite-specific detection is now provided opt-in (`DetectAlphaSatellite`: ~171-bp tandem periodicity + AT-richness + CENP-B box `YTTCGTTGGAARCGGGA` count; `FindCenpBBoxes`) — Willard (1985) / Waye & Willard (1987) / Masumoto et al. (1989). The default `AnalyzeCentromere` (and its `AlphaSatelliteContent` generic repeat score) is unchanged, and Levan (1964) arm-ratio classification is exact. **Residual:** higher-order repeat (HOR) structure / suprachromosomal-family classification is out of scope — detection is monomer-level only. |
| PROBE-DESIGN-001 | TaqMan (5'-nuclease hydrolysis probe) rules are now provided opt-in (`EvaluateTaqManProbe`, `SelectTaqManStrand`): no G at the 5' end (a 5' G quenches the reporter even after cleavage), more Cs than Gs, no run of ≥4 Gs, GC 30–80%, length 18–22 nt, and probe Tm ≥ primer Tm + 10 °C — Applied Biosystems / Thermo Fisher / PREMIER Biosoft. The default generic designer (per-application length/Tm/GC windows + symmetric 5'/3' G·C heuristic penalty) is unchanged. **Residual:** MGB (minor-groove binder), LNA, and dual-quencher probe chemistries are out of scope — the TaqMan rules target standard single reporter/quencher hydrolysis probes. |
| PROBE-VALID-001 | A **gapped** (Smith–Waterman) off-target scan is now provided opt-in (`ScanOffTargetsGapped`): it reuses the library's validated local aligner (`SequenceAligner.LocalAlign`, `BlastDna` scoring) to find off-targets reachable through insertions/deletions that the ungapped Hamming scan misses — Smith & Waterman (1981) / Altschul et al. (1990) — reports per-hit alignment identity and coverage, and **separates the intended on-target match from genuine off-target hits** (`OnTargetHits` vs `OffTargetHits`), with the off-target identity threshold (default 0.75 over the probe length) from Kane et al. (2000). The default ungapped `ValidateProbe` Hamming scan (whose `OffTargetHits` still pools on- plus off-target matches) and `CheckSpecificity` are unchanged. **Residual:** `ScanOffTargetsGapped` is an exhaustive sliding Smith–Waterman scan, not a seeded BLAST k-mer index over a whole genome (no genome-scale seed-and-extend), and it has no duplex-Tm / E-value model. |
| DISORDER-REGION-001 | `RegionType` / `Confidence` labelling heuristics (0.25 enrichment, priority order, confidence rescaling) affect the **labels only**, never the region boundaries (which follow the validated TOP-IDP threshold). |
| MIRNA-TARGET-001 | The TargetScan **context++** score (Agarwal et al. 2015, eLife 4:e05005) is now provided opt-in (`ScoreTargetSiteContextPlusPlus`): a per-site-type multiple-linear-regression sum using the **verbatim fitted coefficients** from the TargetScan distribution `Agarwal_2015_parameters.txt`, with the feature computation/scaling ported exactly from `targetscan_70_context_scores.pl`. It realises the **locally-computable feature subset** — site-type Intercept, Local_AU, sRNA position-1/8 nucleotide identity, and (7mer-A1/6mer) target site-position-8 identity. The default `Score` (Grimson-2007-proportional) and the simplified ΔG model are unchanged; exact seed-match site typing (8mer/7mer-m8/7mer-A1/6mer) is unchanged. **Residual:** the full-transcript context++ features are NOT computed and are reported in `OmittedFeatures` — `3P_score` (3' supplementary pairing), `SPS`, `TA_3UTR` (target-site abundance), `Min_dist`, `SA` (structural accessibility), `PCT`, `Len_3UTR`, `Len_ORF`, `ORF8m`, `Off6m` — so `ContextScorePartial` is a partial context++ (an upper bound, since most omitted coefficients are negative), not the published headline CS. |
| MIRNA-PRECURSOR-001 | A full **MFE-structure-folding** caller is now provided opt-in (`FindPreMiRnaHairpinsByMfe` / `AssessHairpinByMfe`): each candidate is folded by the validated Zuker–Stiegler engine (RNA-STRUCT-001, Turner 2004) and the hairpin is read from the real MFE dot-bracket — single dominant hairpin, stem bp ≥ 16 (Ambros 2003), loop 3–25 nt (Bartel 2004), MFEI ≥ 0.85 (Zhang 2006). It detects natural miRBase precursors (hsa-mir-21, hsa-let-7a-1) the default rejects. The default `FindPreMiRnaHairpins` consecutive-pairing heuristic (stem ≥18 bp, loop 3–25 nt, Turner-2004 MFE score) is unchanged. **Residual:** Drosha/Dicer cleavage-site (mature/star excision-coordinate) prediction and a trained natural-vs-background classifier (miRDeep2) are out of scope. |

## 4. Opt-in compatibility / convention modes (default behaviour unchanged; an alternative convention is available on request)

| Unit | Default | Opt-in alternative |
|------|---------|--------------------|
| SEQ-GC-001, SEQ-GCSKEW-001, SEQ-STATS-001 | GC as a **percentage** | `fraction:true` returns [0,1]; `CalculateGcFraction(GcAmbiguityMode {Remove,Ignore,Weighted})` reproduces Biopython `gc_fraction` (incl. ambiguity-code weights) byte-for-byte. |
| PARSE-VCF-001 | 0-based internal `Position` | `Variant.VcfPosition = Position + 1` exposes the VCF 1-based POS; additive and non-breaking. POS is held as a 32-bit int. |

## 5. Documented scope & API notes (correct contract; recorded for transparency)

| Unit | Note |
|------|------|
| CRISPR-PAM-001 | On the reverse strand, `PamSite.TargetStart` is an index into the reverse-complement string; `Position` is the forward-strand coordinate and is correct. |
| ANNOT-ORF-001 | An optional `requireStartCodon:false` mode seeds ORFs from a run's start (non-canonical); the canonical ATG / alt-start (GTG/TTG) 6-frame path is standard and exact. |
| CODON-USAGE-001 | Returns per-codon **counts** + a total-variation-distance similarity only; frequency-per-1000 / RSCU / within-family fraction live in CODON-RSCU-001 / CODON-STATS-001 / SEQ-CODON-FREQ-001. |
| TRANS-PROT-001 | Ambiguous codon → `'X'` is reachable only via the `string` overload; typed `DnaSequence`/`RnaSequence` reject IUPAC ambiguity at construction. |
| PARSE-FASTA-001 | DNA-only parser (rejects IUPAC / RNA / protein FASTA); on multi-space headers the `Description` keeps a single leading space. |
| PARSE-FASTQ-001 | `DetectEncoding` is a range heuristic that can mis-flag very-high-quality Phred+33 bases as Phred+64; decode with an explicit offset to avoid the ambiguity. |
| PARSE-GENBANK-001 | Feature-table parsing covers the common INSDC location grammar (range / complement / join / `<`,`>` partials). (The multi-line qualifier-wrap defect found during validation was fixed, commit 30d4f7b1.) |
| PARSE-EMBL-001 | A remote-reference location nested **inside** `complement(...)` / `join(...)` is not captured (anchored regex); matches the shared GenBank path and is pre-existing. |
| PROTMOTIF-DOMAIN-001 | WD40 / SH3 / PDZ use simplified consensus regexes (Pfam PF00400 / PF00018 / PF00595 are HMM profiles); the two PROSITE patterns — zinc-finger PS00028 and Walker-A PS00017 — are exact. |

---

## How to read this

- For **research / pipeline** use, every row above is a normal scope boundary of a from-first-principles library — the algorithm is validated correct for its stated contract.
- For **clinical / decision-grade** use, treat §2 (and the clinically-flavoured §1 oncology rows) as layers that require an external validated predictor and clinical sign-off — the library computes the rule, not the trained model behind it.
- §3 heuristics are honest first-principles approximations; swap in a reference tool's output where decision-grade accuracy is required.
- §4 are pure conventions: the default is unchanged, the alternative is opt-in for cross-tool parity.
- §5 are correct-contract notes, surfaced so callers are not surprised by an edge case.
- Each row traces to its per-unit validation report under `docs/Validation/reports/` (all `✅ CLEAN`, `PASS-WITH-NOTES`).
