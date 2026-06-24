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
| PRIMER-TM-001 | Full salt/Mg²⁺-corrected Tm | Wallace rule omits the −7 °C term; Marmur-Doty is the simplified GC form; SantaLucia 1998 nearest-neighbour ΔG is used for 3'-end stability, not for the design Tm. Primer3 `p_obj_fn` penalty reproduced to 1e-10. |
| REP-STR-001 | Imperfect / interrupted tandem repeats | Detects **perfect** STR/microsatellites only (1–6 bp motifs); the Tandem Repeats Finder approximate-repeat model is out of scope. |
| SPLICE-ACCEPTOR-001 | Branch-point / full MaxEntScan model | 3' acceptor scored by a position-weight-matrix consensus + polypyrimidine-tract fraction; no explicit branch-point or maximum-entropy splice model. |
| ONCO-IMMUNE-001 | Full LM22 deconvolution + cohort normalisation | Uses a 5-marker immune/stromal signature matrix (vs CIBERSORT LM22) and single-sample un-normalised ssGSEA, so the ESTIMATE-style purity is a **relative**, not clinically-absolute, value. |
| ONCO-MRD-001 | Smoothed fragment-size model; in-pipeline background | INVAR2 formulas reproduced verbatim, but the fragment-size weighting uses the empirical discrete proportion (not a smoothed KDE) and the per-locus background can be caller-supplied rather than re-estimated from controls. |
| META-BIN-001 | Markov-corrected TETRA + single-copy-marker QC | Tetranucleotide signature uses raw-frequency Pearson distance (vs z-score normalisation); bin completeness/contamination use length-ratio & GC-stddev proxies instead of single-copy marker genes (CheckM). |

## 2. Threshold / aggregation layer (classifies a caller-supplied input; does not predict it upstream)

| Unit | Not modelled / caller-supplied |
|------|--------------------------------|
| ONCO-MHC-001 | The peptide–MHC affinity / %Rank prediction (NetMHCpan-style learned model) is caller-supplied: the library classifies a supplied IC50 / %Rank into strong/weak binder but does not bundle a trained predictor (no redistributable, cross-verifiable model available). |

## 3. Declared from-first-principles heuristics (a deliberate, sourced heuristic — not a trained predictor / reference tool)

| Unit | Heuristic in place of |
|------|-----------------------|
| CHROM-CENT-001 | Centromere region detection is a **generic tandem-repeat-density** heuristic (k-mer recurrence × GC-uniformity), not alpha-satellite / 171-bp-monomer-specific; the result field `AlphaSatelliteContent` carries this generic repeat score. Levan (1964) arm-ratio classification is exact. |
| PROBE-DESIGN-001 | Generic hybridization-probe designer; the symmetric 5'/3' G·C position penalty is a declared heuristic. Chemistry-specific rules (TaqMan no-5'-G / more-C-than-G / probe-Tm+10) are out of this unit's scope. |
| PROBE-VALID-001 | Off-target specificity via an ungapped Hamming-distance scan (not BLAST-grade); `OffTargetHits` pools on- plus off-target matches. |
| DISORDER-REGION-001 | `RegionType` / `Confidence` labelling heuristics (0.25 enrichment, priority order, confidence rescaling) affect the **labels only**, never the region boundaries (which follow the validated TOP-IDP threshold). |
| MIRNA-TARGET-001 | Site efficacy / ΔG scoring is a Grimson-2007-proportional simplification, not the full TargetScan context++ model. Seed-match site typing (8mer/7mer-m8/7mer-A1/6mer) is exact. |
| MIRNA-PRECURSOR-001 | Pre-miRNA hairpin = consecutive-complementary-pairing heuristic (stem ≥18 bp, loop 3–25 nt) with a Turner-2004 MFE score; not a full MFE-structure-folding caller. |

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
| ANNOT-GFF-001 | GFF3 **export** emits `source`/`score` as "." and CDS `phase` as "0" (spec-permitted minimal output); parsing preserves all 9 columns. |
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
