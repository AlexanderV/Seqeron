# Hybridization Probe Design

| Field | Value |
|-------|-------|
| Algorithm Group | MolTools |
| Test Unit ID | PROBE-DESIGN-001 |
| Related Projects | N/A |
| Implementation Status | Complete (Primer3 probe picker + heuristic ranking designer) |
| Last Reviewed | 2026-10-01 |

## 1. Overview

Hybridization probe design generates oligonucleotide probes that detect specific nucleic acid sequences through complementary binding. In this repository, probe design supports multiple applications including FISH, DNA microarrays, Northern blots, qPCR, and Southern blots. The documented implementation is a heuristic candidate-generation and ranking workflow driven by GC content, melting temperature, self-complementarity, secondary structure, simple repeats, and application-specific parameter ranges.

A **Primer3 hybridization-probe picker** (`DesignProbesPrimer3`) reproduces Primer3's `PRIMER_TASK=pick_hyb_probe_only` (internal-oligo) selection exactly (§4.4).

An **opt-in TaqMan (5'-nuclease hydrolysis probe) rule set** is also provided (`EvaluateTaqManProbe`, `SelectTaqManStrand`) for the qPCR/TaqMan chemistry. It is separate from and does not alter the generic designer, which remains the default.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Hybridization probes are used in sequence-detection assays such as FISH, DNA microarrays, Northern blots, and qPCR. Probe quality depends on a balance among target complementarity, melting temperature, GC content, and avoidance of self-structure or repetitive motifs that can reduce specificity. The original document ties these considerations to nucleic-acid thermodynamics and to application-specific probe lengths and Tm targets. Sources: Wikipedia (Hybridization probe, Fluorescence in situ hybridization, DNA microarray, Nucleic acid thermodynamics), SantaLucia (1998), Breslauer et al. (1986).

### 2.2 Core Model

**Ranking (`ProbeParameters.Ranking`, audit round 3, A3-11).** Two orders are available; the candidate set (windows with additive score > 0) is the same for both:

- `ProbeRanking.AdditiveScore` (default, unchanged) — the additive score below, score descending (ties by length, then start). **This score is a library heuristic**: the penalty values (0.3, 0.3, 0.2, 0.2, 0.15, 0.1, 0.02) have no published source; the score ranks candidates but is not a hybridization probability.
- `ProbeRanking.Primer3Penalty` — **the sourced ranking**: Primer3's internal-oligo objective `p_obj_fn` (`libprimer3.cc`, `OT_INTL` branch; PRIMER_INTERNAL_n_PENALTY) with Primer3's default internal-oligo weights (PRIMER_INTERNAL_WT_TM_GT/_LT = PRIMER_INTERNAL_WT_SIZE_GT/_LT = 1, all other weights 0), i.e. $|T_m - OptTm| + |N - OptLength|$ with the probe Tm at the parameters' conditions, computed by the canonical `PrimerDesigner.CalculatePrimer3Penalty` (the same call `DesignProbesPrimer3` makes) and reported in `Probe.Primer3Penalty`; order = Primer3 `primer_rec_comp` (penalty ↑, start ↓, length ↑); a probe whose Tm is not computable (non-ACGT base, which Primer3 rejects) ranks last with `Primer3Penalty = null`. `ProbeParameters.OptTm` / `OptLength` default to Primer3's PRIMER_INTERNAL_OPT_TM = 60 °C / PRIMER_INTERNAL_OPT_SIZE = 20 and must lie within [MinTm, MaxTm] / [MinLength, MaxLength] (`ArgumentOutOfRangeException`, Primer3 `_pr_data_control` "Optimum internal oligo Tm lower than minimum or higher than maximum" / "PRIMER_INTERNAL_{OPT,DEFAULT}_SIZE > MAX_SIZE / < MIN_SIZE"; checked eagerly) — so a preset such as Microarray (82–90 °C, 50–60 nt) needs its own optima. `Score` stays the additive score. Cross-check: primer3-py 2.3.1 `design_primers` PRIMER_TASK = pick_hyb_probe_only (PRIMER_PICK_INTERNAL_OLIGO, internal size 18–27, Tm 57–63, GC 20–80 %, poly-X 5, optima (60, 20), (61.5, 22), (59, 19), (60, 24)) on four random 70-nt templates: every one of the 48 Primer3 probes is a `DesignProbes` candidate with the same PRIMER_INTERNAL_n_PENALTY and Tm (|Δ| ≤ 1e-9) and the same relative order, including the start-descending ties (`ProbeDesigner_Primer3Ranking_Tests`). `DesignProbesPrimer3` remains Primer3's complete picker (its acceptance limits and five-prime-problem enumeration as well as its ranking).

The additive candidate score starts at `1.0` and is reduced by a fixed set of penalties:

| Factor | Penalty |
|--------|---------|
| GC content outside the configured range | `-0.3` |
| Tm outside the configured range | `-0.3` |
| Homopolymer run above the configured maximum | `-0.2` |
| Self-complementarity above the configured maximum | `-0.2` |
| Secondary-structure potential | `-0.15` |
| Simple repeats | `-0.1` |
| Starts with `G/C` | `-0.02` |
| Ends with `G/C` | `-0.02` |

**Probe Tm.** Every probe Tm (designer, tiling, molecular beacon, TaqMan check, `AnalyzeOligo`) is Primer3's `seqtm` (`PrimerDesigner.CalculateMeltingTemperaturePrimer3`, identical to primer3-py `calc_tm(..., max_nn_length=ProbeParameters.MaxNearestNeighborLength)`, default 36 = the MAX_NN_TM_LENGTH libprimer3 uses when picking oligos; the Microarray preset uses 60): SantaLucia (1998) nearest-neighbour ΔH/ΔS with the SantaLucia entropy salt correction and the von Ahsen Mg²⁺→Na⁺-equivalent for ≤ 36 nt, and Primer3 `long_seq_tm` = 81.5 + 16.6·log10([Na⁺]eq) + 0.41·%GC − 600/N above 36 nt, at the hybridization conditions in `ProbeParameters` (`ProbeParameters` default = Primer3 internal-oligo defaults: 50 nM probe, 50 mM monovalent, 0 Mg²⁺, 0 dNTP; the Microarray preset — the `DesignProbes` / tiling default — uses OligoArray's 1 M / 1 µM, see §4.2). A probe with a non-ACGT base has no computable Tm (reported 0 with a warning and the Tm penalty; `AnalyzeOligo`/TaqMan report NaN).

**Self-structure.** With `StructureScreen = Thermodynamic` (default) a probe of ≤ `ThermodynamicScreenMaxLength` (default 60) nt of A/C/G/T is screened like a Primer3 hybridization probe: ntthal self-dimer (ANY), 3′ self-dimer (END1) and hairpin Tm (`PrimerDesigner.CalculatePrimer3OligoStructure`, bit-exact to primer3-py) must not exceed `MaxStructureTm` (PRIMER_INTERNAL_MAX_SELF_ANY_TH / _SELF_END_TH / _HAIRPIN_TH = 47 °C); the self-dimer limit drives the "self-complementarity" penalty and the hairpin limit the "secondary structure" penalty. Longer probes (thal.c THAL_MAX_ALIGN, default 60 — see the opt-in below), non-ACGT probes and `StructureScreen = Heuristic` use the fallback screens: the self-dimer criterion is Primer3's alignment-mode internal-oligo screen (PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT = 0, `oligo_compl`): dpal `self_any` (`PrimerDesigner.CalculatePrimerSelfAnyComplementarity`, DPAL_LOCAL) and `self_end` (`CalculatePrimerSelfEndComplementarity`, DPAL_GLOBAL_END), which have no length limit, must not exceed `MaxSelfAny` / `MaxSelfEnd` (PRIMER_INTERNAL_MAX_SELF_ANY / _SELF_END = 12.00); the hairpin criterion stays the sequence-only inverted-repeat stem ≥ 4 bp with a 3-nt loop and ≥ 80 % matches (with the default THAL_MAX_ALIGN = 60; for a thermodynamic hairpin of longer probes use the opt-in below). Note: Primer3 itself picks internal oligos of ≤ 36 nt, so 12.00 is calibrated for short oligos; random 200–500-nt probes exceed it by chance in roughly half of the cases (a −0.2 ranking penalty in `DesignProbes`) — raise `MaxSelfAny`/`MaxSelfEnd` for long FISH/blot probes. `MaxSelfComplementarity` (fold-back fraction limit) is no longer used by any screen (kept for source compatibility).

**Simple repeats** are di-/trinucleotide microsatellites of ≥ 4 complete copies (canonical `RepeatFinder.FindMicrosatellites`, MISA convention: primitive units, A/C/G/T only).

#### 2.2.1 TaqMan (5'-nuclease hydrolysis probe) rules — opt-in

For TaqMan qPCR chemistry, the published Applied Biosystems / Thermo Fisher guidelines [7][8][9] add chemistry-specific constraints that the generic designer does not enforce. `EvaluateTaqManProbe` checks each as a boolean; `SelectTaqManStrand` chooses the better strand.

| Rule | Threshold | Source |
|------|-----------|--------|
| No G at the 5' end | first base ≠ `G` (a 5' G adjacent to the reporter dye quenches reporter fluorescence even after cleavage) | [7][9] |
| More Cs than Gs | `count(C) > count(G)` | [7] |
| No run of ≥4 consecutive Gs | max G-run `< 4` | [7] |
| G+C content | `30%–80%` | [7] |
| Probe length | `18–22 nt` (default; configurable) | [7] |
| Probe Tm above primer Tm | probe Tm `≥ primerTm + 10 °C` | [7][8] |

If a 5' G cannot be avoided on the sense strand, the probe is designed on the complement (antisense) strand [8]; `SelectTaqManStrand` returns the reverse-complement strand when it better satisfies the rules.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | The ranking pass retains only raw-score-positive candidates before optional specificity rescaling | `EnumerateRankedProbes(...)` rejects candidates when `score <= 0`; the suffix-tree overload then rescales (or filters) every retained candidate |
| INV-02 | Probe GC content is a fraction: G+C over the valid (A/C/G/T/U) bases, N excluded | canonical `CalculateGcFractionFast` (the `DesignProbes` prefix sums count with the same `CountGcAndValidNucleotides`); = Primer3 `gc_and_n_content` (100·num_gc/num_gcat) |
| INV-03 | `CheckSpecificity(...)` returns `0` for no hits, `1` for a unique hit, and `1 / hits` otherwise | That mapping is explicit in source |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `targetSequence` | `string` | required | Sequence from which probe candidates are generated | Uppercased before processing |
| `parameters` | `ProbeParameters?` | `Defaults.Microarray` | Application-specific probe design limits | Includes length, Tm, GC, homopolymer, and self-complementarity thresholds |
| `maxProbes` | `int` | `10` | Maximum number of returned probes | Applied after ranking (and, in the genome-index overload, after the specificity filter / re-rank); `<= 0` returns none |
| `parameters.Ranking` | `ProbeRanking` | `AdditiveScore` | `AdditiveScore` (library heuristic) or `Primer3Penalty` (Primer3 internal-oligo `p_obj_fn`, sourced) | Enum value must be defined |
| `parameters.OptTm` / `OptLength` | `double` / `int` | `60` / `20` | PRIMER_INTERNAL_OPT_TM / _OPT_SIZE of the Primer3 ranking | Within [MinTm, MaxTm] / [MinLength, MaxLength] when `Ranking = Primer3Penalty` |
| `genomeIndex` | `ISuffixTree` | required for specificity overload | Pre-built suffix tree for genome-wide uniqueness filtering | Used only by the overload with specificity checking |
| `requireUnique` | `bool` | `true` | Whether non-unique probes are excluded when `genomeIndex` is provided | Filters candidates with specificity `< 1.0` |
| `bothStrands` | `bool` | `false` | Genome-index overload: also count the probe's reverse-complement occurrences (`CheckSpecificity(…, bothStrands)`) — a probe binds the other strand of a double-stranded genome wherever its reverse complement occurs, as blastn `-strand both` searches | Optional (audit round 4, A4-2); reverse-palindromic probes counted once; default = indexed strand only |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Sequence` | `string` | Probe sequence |
| `Start` / `End` | `int` | Probe coordinates in the source sequence |
| `Tm` | `double` | Melting-temperature estimate |
| `GcContent` | `double` | GC fraction |
| `Score` | `double` | Additive library-heuristic quality score (unsourced penalty values) |
| `Primer3Penalty` | `double?` | Primer3 PRIMER_INTERNAL_n_PENALTY (`p_obj_fn`) with `Ranking = Primer3Penalty`; `null` otherwise or when the Tm is not computable |
| `Type` | `ProbeType` | Probe category such as `Standard`, `Tiling`, `Antisense`, `LNA`, or `MolecularBeacon` |
| `Warnings` | `IReadOnlyList<string>` | Quality warnings recorded during evaluation |

### 3.3 Preconditions and Validation

`DesignProbes(...)` returns no probes when the target sequence is null, empty, shorter than the configured minimum length, or when `maxProbes <= 0`. All sequences are converted to uppercase before processing. The suffix-tree overload walks **every** candidate in the ranking order (lazily, stopping once `maxProbes` probes have been produced) and applies the specificity value returned by `CheckSpecificity(...)` (with `bothStrands` on both strands of the indexed genome — occurrences of the probe and of its reverse complement, a reverse-palindromic probe once; default the indexed strand only): with `requireUnique` a candidate occurring more than once is dropped (the surviving probes keep their score and order), otherwise the score is multiplied by the specificity and the probes are re-ranked on the scaled score (stable — equal scores keep the documented `(length, start)` tie order; the `Primer3Penalty` order is unchanged because the specificity does not enter the penalty). With `requireUnique = false` a probe whose specificity is `0` (not present in the index) therefore still appears, with final score `0`, last.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalize the input sequence to uppercase.
2. Precompute GC prefix sums for the target sequence.
3. Enumerate all candidate windows within the configured length range.
4. Evaluate each candidate for GC content, Tm, homopolymers, self-complementarity, secondary structure, repeats, and terminal G/C penalties.
5. Keep only candidates with positive raw scores, sort them by score (default) or by Primer3 internal-oligo penalty (`Ranking = Primer3Penalty`), and return the top results.
6. In the genome-index overload, walk the whole ranked candidate list lazily and apply suffix-tree specificity (uniqueness filter, or score scaling followed by a stable re-rank) before yielding up to `maxProbes` probes.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Application-specific defaults from `ProbeDesigner.Defaults`:

| Application | Length (bp) | Tm (°C) | GC Content | Tm conditions / scale | Notes |
|-------------|-------------|---------|------------|-----------------------|-------|
| Microarray | 50-60 | 82-90 | 0.40-0.60 | OligoArray 2.0: [Na⁺] 1 M, 1 µM oligo, nearest-neighbour Tm over the whole oligo (`MaxNearestNeighborLength` 60) | Default parameter set; Tm window and conditions from Rouillard, Zuker & Gulari 2003 |
| FISH | 200-500 | 70-90 | 0.35-0.65 | Primer3 probe conditions (50 nM / 50 mM), `long_seq_tm` | Library convention; attainable 71.25–85.35 |
| Northern Blot | 100-300 | 65-80 | 0.40-0.60 | Primer3 probe conditions, `long_seq_tm` | Library convention; attainable 70.30–82.50 |
| qPCR | 20-30 | 68-70 | 0.30-0.80 | Primer3 probe conditions, NN ≤ 36 nt | TaqMan probe Tm/GC per Applied Biosystems Primer Express guidelines [10] |
| Southern Blot | 150-500 | 65-75 | 0.35-0.65 | Primer3 probe conditions, `long_seq_tm` | Library convention; attainable 70.32–85.35 |

The implementation uses prefix sums for O(1) GC lookup and a suffix-tree overload for O(m) uniqueness checking of each candidate probe. The thermodynamic self-structure screens only lower a score, so they are evaluated lazily in descending base-score order with a branch-and-bound stop (identical result to an exhaustive evaluation; ties by length, then start).

Every preset's Tm window is reachable on its own Tm scale for its length × G+C window (audit round 3, A3-10; primer3-py `calc_tm` witnesses and exact `long_seq_tm` ranges locked in `ProbeDesigner_Primer3Probe_Tests`). The Microarray preset follows OligoArray 2.0 (Rouillard, Zuker & Gulari 2003): nearest-neighbour Tm at [Na⁺] = 1 M and 1 µM oligo, Tm window 82–90 °C — computed as Primer3 `oligotm` (SantaLucia 1998, C/4 for non-self-complementary oligos; identical to Biopython `Tm_NN(DNA_NN3)` at the same conditions) with `ProbeParameters.MaxNearestNeighborLength` = 60 (Primer3 `seqtm` `nn_max_len`); attainable Tm of random 50–60-mers with 40–60 % G+C ≈ 83–97 °C. Its former window 75–85 °C at the 50 mM Primer3 conditions was unreachable (`long_seq_tm` maximum 74.50 °C for a 60-mer with 60 % G+C). The Microarray structure screen runs at the same 1 M / 1 µM conditions. `ValidateProbe`, `AssessCrossHybridization`, `DesignMolecularBeacon` and `AnalyzeOligo` keep the Primer3 probe conditions (50 nM / 50 mM, NN ≤ 36 nt) as their default. The FISH, Northern and Southern Tm windows are library conventions (no published Tm window for these long probes was found) — set `MinTm`/`MaxTm` and the hybridization conditions for your protocol. Lengths: microarray long oligos 50–60 nt (Kane et al. 2000 50-mers; Agilent 60-mers).

### 4.4 Primer3 hybridization-probe picker (`DesignProbesPrimer3`)

Exact port of Primer3 `make_internal_oligo_list` → `pick_primer_range` → `calc_and_check_oligo_features` → `p_obj_fn` (internal-oligo branch) in the default thermodynamic mode, with `Primer3ProbeSettings` = the `PRIMER_INTERNAL_*` defaults (size 18/20/27, Tm 57/60/63 °C, GC 20–80 %, max poly-X 5, Ns 0, self-any/self-end/hairpin Tm ≤ 47 °C, 50 mM monovalent, 0 Mg²⁺, 0 dNTP, 50 nM). Accepted windows are ranked by penalty = |Tm − OptTm| + |len − OptSize| (+ PRIMER_INTERNAL_WT_GC_PERCENT_GT/_LT × the GC deviation from PRIMER_INTERNAL_OPT_GC_PERCENT when `WeightGcPercentGt/Lt` ≠ 0, optimum `OptGcPercent`, null = 50 %; illegal conditions — salt or DNA concentration ≤ 0, negative divalent/dNTP — throw as Primer3 `_pr_data_control`; audit round 3, A3-1, primer3-py parity in `Primer_Design.md` §5.3) and ordered as `primer_rec_comp` (penalty ↑, start ↓, length ↑); `numReturn` = PRIMER_NUM_RETURN (5; < 1 throws `ArgumentOutOfRangeException`, Primer3 `_pr_data_control` "PRIMER_NUM_RETURN < 1" — audit round 3, A3-15). Verified against primer3-py 2.3.1 `design_primers` on 950 random templates/settings (positions, Tm, penalty, SELF_ANY_TH, SELF_END_TH, HAIRPIN_TH identical, max |Δ| = 0).

Mishybridization library (PRIMER_INTERNAL_MISHYB_LIBRARY; audit round 3, A3-3 part 2): `Primer3ProbeSettings.MishybLibrary` (a `PrimerMisprimingLibrary`, primer3-py `mishyb_lib`: name → sequence, optional `*weight`, IUPAC codes, Primer3 adds each entry's reverse complement), `MaxLibraryMishyb` (PRIMER_INTERNAL_MAX_LIBRARY_MISHYB, 12), `WeightLibraryMishyb` (PRIMER_INTERNAL_WT_LIBRARY_MISHYB, 0), `LibraryAmbiguityCodesConsensus` (PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS, 0). Every window passing the structure limits is scored by `libprimer3.c` `oligo_repeat_library_mispriming` (OT_INTL): $w_i = weight_i \cdot align(probe, y_i)$, dpal unanchored local alignment (`local` / `local_ambig`: +1 match, −1 mismatch, −0.25 vs N, −2 per single-base gap, max gap 1; an entry < 3 nt scores its length; `PrimerDesigner.CalculateLibraryMishyb`, the same dpal port as the primer mispriming library); a window with any $w_i$ > (short) `MaxLibraryMishyb` is rejected and — Primer3 OP_HIGH_SIM_TO_NON_TEMPLATE_SEQ is a "five-prime problem" — ends the 5′ extension of that 3′ end; the penalty gains `WeightLibraryMishyb` × score (`p_obj_fn`). `Primer3Probe.LibraryMishyb` / `LibraryMishybName` report Primer3's `repeat_sim.max` score and entry (PRIMER_INTERNAL_n_LIBRARY_MISHYB; primer3-py key `…_LIBRARY_MISPRIMING`). A weight without a library and, in alignment mode, a limit > 32767 throw (`_pr_data_control`). The same core serves `PrimerDesigner.DesignPrimers` with PRIMER_PICK_INTERNAL_OLIGO (`Primer_Design.md` §2.2 item 7). Cross-check: primer3-py 2.3.1 `design_primers(mishyb_lib=…)` with PRIMER_TASK = pick_hyb_probe_only: 814/814 random 40–220-nt templates (7875 probes; random 1–8-entry libraries — template fragments, mutations, IUPAC codes, 1–2-nt and N-rich entries, `*w` weights — PRIMER_INTERNAL_MAX_LIBRARY_MISHYB ∈ {6, 8, 10, 12, 12.9, 15, 20}, PRIMER_INTERNAL_WT_LIBRARY_MISHYB ∈ {0, 0.1, 0.5, 1}, both PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT modes, PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS ∈ {0, 1}, PRIMER_NUM_RETURN 5–50; probes rejected by the library in 699 templates) identical on position, PRIMER_INTERNAL_k_PENALTY (|Δ| ≤ 1e-9) and the library score + entry name.

Fraction bound (PRIMER_ANNEALING_TEMP; audit round 3, A3-5 part 1): `Primer3ProbeSettings.AnnealingTemperature` (default −10 = off; ≤ 100), `MinBound` / `MaxBound` / `OptBound` (PRIMER_INTERNAL_MIN/MAX/OPT_BOUND, −10 / 110 / 97 %), `WeightBoundGt` / `WeightBoundLt` (PRIMER_INTERNAL_WT_BOUND_GT/_LT, 0). With $T_a$ > 0 every window gets Primer3's `oligotm` fraction bound at $T_a$ (`PrimerDesigner.CalculateFractionBoundPrimer3` at the probe conditions; `Primer3Probe.Bound` = PRIMER_INTERNAL_n_BOUND) and is rejected outside [MinBound, MaxBound] right after the Tm check (not a five-prime problem). The bound terms of the internal-oligo `p_obj_fn` are not gated by $T_a$: without it (or for a window > 36 nt) the bound is Primer3's OLIGOTM_ERROR −999999.9999, so `WeightBoundLt` w adds w·(OptBound + 999999.9999), reproduced. OptBound outside [MinBound, MaxBound] or $T_a$ > 100 → `ArgumentOutOfRangeException` (`_pr_data_control`). For the internal oligo of a primer pair the annealing temperature is `PrimerParameters.AnnealingTemperature` (one global Primer3 setting). Cross-check: 85 pick_hyb_probe_only lists (432 probes) within the F44 harness identical to primer3-py 2.3.1 (position, penalty, bound).

Sequence quality (SEQUENCE_QUALITY; audit round 3, A3-5 part 2a): `DesignProbesPrimer3(template, settings, numReturn, sequenceQuality)` takes one integer quality per template base; every window gets `Primer3Probe.MinSequenceQuality` (PRIMER_INTERNAL_n_MIN_SEQ_QUALITY = min(PRIMER_QUALITY_RANGE_MAX, qualities over the window), `PrimerDesigner.CalculateSequenceQualityPrimer3`), is rejected when it is below `Primer3ProbeSettings.MinQuality` (PRIMER_INTERNAL_MIN_QUALITY, checked after the GC check; a five-prime problem) and gets `WeightSequenceQuality` × (`QualityRangeMax` − min quality) (PRIMER_INTERNAL_WT_SEQ_QUAL, the last internal-oligo `p_obj_fn` term). PRIMER_INTERNAL_WT_END_QUAL (`WeightEndQuality`) has no effect (never read by Primer3 2.3.1). `_pr_data_control`: quality length ≠ template, MinQuality ≠ 0 without quality or outside [`QualityRangeMin`, `QualityRangeMax`] (0 / 100), a value outside the range, a weight without quality → `ArgumentException`. For the internal oligo of a primer pair the quality comes from `PrimerPairOptions.SequenceQuality` and the range from `PrimerParameters.QualityRangeMin/Max`. Cross-check: the F45 harness (primer3-py 2.3.1, probe-only and PRIMER_PICK_INTERNAL_OLIGO cases, position / penalty / MIN_SEQ_QUALITY identical).

GC optimum and lower-case masking (audit round 3, A3-25 + A3-26): `Primer3ProbeSettings.OptGcPercent` (PRIMER_INTERNAL_OPT_GC_PERCENT) is undefined (null) by default, as in Primer3's code (`DEFAULT_OPT_GC_PERCENT = PR_UNDEFINED_INT_OPT`); a non-zero `WeightGcPercentGt/Lt` without it throws `ArgumentException` with Primer3's `_pr_data_control` message "Hyb probe GC content is part of objective function while optimum gc_content is not defined". `Primer3ProbeSettings.LowercaseMasking` (PRIMER_LOWERCASE_MASKING, default false) rejects every window whose 3′ base (its rightmost template base) is a lower-case a/c/g/t of the template as given (`is_lowercase_masked` on `trimmed_orig_seq`, before every other check); lower case elsewhere is accepted. Verified against primer3-py 2.3.1 pick_hyb_probe_only on random mixed-case templates (F46).

### 4.5 Oligo properties

- Molecular weight: canonical `SequenceStatistics.CalculateNucleotideMolecularWeight` (= Biopython `molecular_weight`, single-stranded; `CalculateMolecularWeight(seq)` infers RNA when the oligo has U and no T).
- ε260: `CalculateExtinctionCoefficient` = sum of mononucleotide values (no hypochromicity); `CalculateExtinctionCoefficientNearestNeighbor` = nearest-neighbour model ε = Σ ε(NᵢNᵢ₊₁) − Σ ε(internal Nᵢ) with the Cantor, Warshaw & Shapiro (1970) DNA / Warshaw & Tinoco (1966) RNA tables [11].
- Molecular beacon: loop Tm window [T+7, T+10] °C and stem-loop (ntthal hairpin) Tm ≥ T+7 °C for a detection temperature T (Tyagi & Kramer 1996 design rules [12]); without T the window is 55–65 °C. The stem-loop Tm is reported for beacons ≤ `maxAlignLength` nt (default 60 = THAL_MAX_ALIGN; opt-in up to 10 000).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `DesignProbes` | `O(n × m)` | `O(k)` | `n` is sequence length, `m` is the scanned length range, and `k` is the number of retained candidates |
| `CheckSpecificity` | `O(m)` | `O(1)` | Uses suffix-tree lookups per probe |
| `DesignTilingProbes` | `O(n)` over fixed-length windows | `O(k)` | Produces overlapping probes for coverage |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ProbeDesigner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs), [ThermoConstants.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Infrastructure/ThermoConstants.cs)

- `ProbeDesigner.DesignProbes(string, ProbeParameters?, int)`: Main probe-generation and ranking routine.
- `ProbeDesigner.DesignProbes(string, ISuffixTree, ProbeParameters?, int, bool, bool)`: Uniqueness-aware overload using a suffix tree (lazy walk of all candidates; `EnumerateRankedProbes(...)` is the shared ranking stream; opt-in `bothStrands` counts reverse-complement occurrences too).
- `ProbeDesigner.DesignTilingProbes(...)`: Generates overlapping tiling probes for coverage.
- `ProbeDesigner.CheckSpecificity(string, ISuffixTree)`: Maps suffix-tree hit counts to a specificity score.
- `ProbeDesigner.EvaluateTaqManProbe(string, double?, int, int)`: Opt-in TaqMan rule check; returns a `TaqManProbeEvaluation` with one boolean per rule and a `PassesAll` conjunction.
- `ProbeDesigner.SelectTaqManStrand(string, double?)`: Chooses the sense strand or its reverse complement, whichever better satisfies the TaqMan rules (no 5'-G, more C than G first).
- `ProbeDesigner.DesignProbesPrimer3(string, Primer3ProbeSettings?, int, IReadOnlyList<int>?)`: Primer3 `pick_hyb_probe_only` picker (§4.4), incl. the mishybridization library.
- `PrimerDesigner.CalculateLibraryMishyb(string, PrimerMisprimingLibrary, bool)`: Primer3 internal-oligo library mishybridization score and entry.

### 5.2 Current Behavior

The implementation evaluates candidates with prefix-sum GC optimization and keeps the raw-score-positive ones. Probe sequences are uppercased before evaluation. `EnumerateRankedProbes(...)` yields them lazily in the ranking order — the self-structure screens (and the specificity scaling) can only lower a score, so finished candidates are released through a priority queue as soon as no unfinished candidate can outrank them, and the enumeration equals an exhaustive evaluation followed by a stable sort. The suffix-tree overload consumes that stream, so it considers **every** candidate (`maxProbes * 5` shortlist removed, audit round 4, A4-1, F59) while costing only the candidates a prefix of the output needs: it enforces uniqueness (`CheckSpecificity(...) < 1` dropped) or scales the score by `1 / hitCount` and re-ranks on the scaled score; with `requireUnique = false` a probe with specificity `0` stays in the output with score `0` (ranked last). `DesignTilingProbes(...)` includes suboptimal probes when needed for coverage and reports coverage, mean Tm, and Tm range. The source also defines probe types `Standard`, `Tiling`, `Antisense`, `LNA`, and `MolecularBeacon`.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Application-specific probe-length, GC, and Tm windows.
- Heuristic penalties for GC, Tm, self-complementarity, secondary structure, and repeats.
- Genome-index-based uniqueness checking through a suffix tree, over the **whole** ranked candidate list (`requireUnique` returns the best `maxProbes` unique probes; without it the score is scaled by `1 / hits` and the probes are re-ranked) — audit round 4, A4-1, F59; opt-in on both strands of a double-stranded genome (`bothStrands`, as blastn `-strand both`) — audit round 4, A4-2, F60.
- Opt-in TaqMan rules (no 5'-G, more C than G, no ≥4-G run, GC 30–80%, length 18–22 nt, probe Tm ≥ primer Tm + 10 °C) and strand selection. [7][8][9]

**Intentionally simplified:**

- `DesignProbes` ranks by default with a fixed additive penalty score — a library heuristic without a published source; the sourced Primer3 internal-oligo objective is the opt-in `ProbeParameters.Ranking = ProbeRanking.Primer3Penalty` (and Primer3's complete picker is `DesignProbesPrimer3`); **consequence:** additive scores rank candidates but are not hybridization probabilities.
- Probes > 60 nt (Northern/Southern/FISH presets): by default ntthal does not align a self-structure whose two strands are both > 60 nt (thal.h `THAL_MAX_ALIGN` = 60; primer3-py raises "At least one sequence must be equal to or shorter than 60bp for thermodynamic calculations"), so their self-dimer criterion is Primer3's alignment-mode internal-oligo screen (dpal `self_any` / `self_end` ≤ `MaxSelfAny` / `MaxSelfEnd` = 12.00, no length limit; F37) and the hairpin criterion is the sequence-only inverted-repeat stem screen (§2.2).
- **Opt-in thermodynamic screen for > 60 nt (audit round 3, A3-9, F55).** `THAL_MAX_ALIGN` is a compile-time guard of thal.c (`#ifndef THAL_MAX_ALIGN` in thal.h; used only by `thal_check_errors` — the DP tables are allocated from the actual lengths), not a limit of the recursions. `ProbeParameters.ThermodynamicScreenMaxLength` (default 60 = unchanged behaviour; 60–10 000 = THAL_MAX_SEQ, else `ArgumentOutOfRangeException`) raises it: ACGT probes up to that length get the full ntthal self-dimer (ANY), 3′ self-dimer (END1) and hairpin Tm ≤ `MaxStructureTm` screen in `DesignProbes` / `DesignTilingProbes` / `ValidateProbe`. The engines are exposed as `PrimerDesigner.CalculateHairpinThermodynamicsNtthal` / `CalculateHairpinStructureNtthal` / `CalculateDimerThermodynamicsNtthal` / `CalculateDimerStructureNtthal` / `CalculatePrimer3OligoStructure` overloads with a trailing `maxAlignLength` (constants `PrimerDesigner.NtthalMaxAlignLength` = 60, `NtthalMaxSequenceLength` = 10 000). Verified against thal.c + thal_parameters.c of primer3-py 2.3.1 compiled with `-DTHAL_MAX_ALIGN=10000` (driven through `thal()` exactly as primer3-py's `calc_hairpin` / `calc_homodimer` / `calc_end_stability`; the same build reproduces primer3-py on 600 ≤ 60-nt cases): 560 oligos of 61–120 nt (random, designed stem-loops with mismatches, palindromes, GC-rich) × hairpin/ANY/END1 at four condition sets (50/0/0/50 nM 37 °C; 50/1.5/0.6; 100/3/0.8/250 nM 55 °C max_loop 20; 10/0/0/1 µM 25 °C max_loop 10) → 0 mismatches (|ΔTm| = 0, |ΔG| ≤ 1.5e-11 cal/mol). Caveats: (1) Primer3 chose 60 as "the maximum reasonable length for nearest neighbor models … only two states of melting" (thal.h) — beyond it the same single-structure two-state model is applied, which over-simplifies long-probe melting; (2) cost is O(n²·30²) per probe and screen (the three ntthal runs together ≈ 0.08 s for a 120-mer, 0.3 s for a 250-mer, 1.4 s for a 500-mer; Release build), multiplied by the number of windows `DesignProbes` evaluates — keep the value to the probe lengths in use. The same opt-in reaches `AssessCrossHybridization` (the conditions' `ThermodynamicScreenMaxLength` is the THAL_MAX_ALIGN of the site duplex; the duplex Tm is computed when the probe **or** the site is within it, as thal.c `thal_check_errors` — so a > 60-nt probe on a ≤ 60-nt site gets primer3-py's `calc_heterodimer` value by default) and `DesignMolecularBeacon(…, maxAlignLength)` (stem-loop Tm of beacons up to that length); MCP `design_probes` / `design_tiling_probes` / `design_antisense_probes` / `validate_probe` (`thermodynamic_screen_max_length`) and `design_molecular_beacon` (`max_align_length`) expose it (audit round 3, A3-27 + A3-28, F56; 600/600 random probe/site strands of 40–120-nt probes and 60/60 beacons of 61–112 nt = thal.c `-DTHAL_MAX_ALIGN=10000`, the default equal to primer3-py on all 563 strands it accepts).

**Not implemented:**

- Database-style alignment or experimentally calibrated hybridization prediction; **users should rely on:** external probe-validation workflows when those are required.
- Dual-quencher probe chemistries and the quantitative MGB ΔTm (vendor `MGB_dds` parameters unpublished — BLOCKED, B07 F28). The TaqMan rules implemented here target standard (single reporter/quencher) hydrolysis probes.

**Implemented elsewhere:** LNA-modified probe Tm — `PrimerDesigner.CalculateMeltingTemperatureNNLna` / `CalculateNearestNeighborThermodynamicsLna` (Owczarzy 2011 default, McTigue 2004 option; MELTING 5 parity; [LNA_Adjusted_Nearest_Neighbor_Tm.md](LNA_Adjusted_Nearest_Neighbor_Tm.md), PROBE-LNATM-001); the citable 3′-MGB design rules (Kutyavin 2000: 3′ attachment, 12–20-mer) — `ProbeDesigner.EvaluateMgbProbeDesign`. `DesignProbes` itself designs unmodified DNA probes (`ProbeType.LNA` is a label only).

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Probe Tm = Primer3 `seqtm` at Primer3 internal-oligo conditions | Assumption | Primer Express' own Tm scale (on which the 68–70 °C TaqMan window is defined) is proprietary; Primer3 is the open reference | accepted | Conditions are configurable in `ProbeParameters` / TaqMan arguments |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null or empty target sequence | Returns no probes | Explicit early return in source |
| Sequence shorter than `MinLength` | Returns no probes | No valid candidate window exists |
| Candidate with score `<= 0` | Rejected | The evaluator returns `null` for non-positive scores |
| Specificity check with no genome hits | Returns `0` | The probe does not match the indexed genome |
| `DesignProbesPrimer3` with `numReturn` < 1 | `ArgumentOutOfRangeException` | Primer3 `_pr_data_control` "PRIMER_NUM_RETURN < 1" |
| `DesignProbesPrimer3` with `WeightGcPercentGt/Lt` ≠ 0 and no `OptGcPercent` | `ArgumentException` | Primer3 "Hyb probe GC content is part of objective function while optimum gc_content is not defined" |
| `DesignProbesPrimer3` with `WeightLibraryMishyb` ≠ 0 and no mishyb library | `ArgumentException` | Primer3 "Internal oligo mispriming score is part of objective function while mishyb library is not defined" |

### 6.2 Limitations

`DesignProbes` ranks with heuristic additive penalties over Primer3-exact measurements (seqtm Tm; ntthal self-dimer / hairpin Tm for ACGT probes ≤ `ThermodynamicScreenMaxLength` (default 60 = Primer3 THAL_MAX_ALIGN; opt-in up to 10 000, F55/F56), dpal self_any / self_end and a sequence-only hairpin stem otherwise); specificity is exact-hit uniqueness through the suffix tree over **all** candidates (mismatch-aware off-target assessment is `ValidateProbe` / `ScanOffTargetsGapped`, PROBE-VALID-001). A unique-probe request can therefore cost a suffix-tree lookup per candidate window when few candidates are unique (the lookup is O(probe length); the walk stops at `maxProbes` probes). It is suitable for fast candidate generation and filtering, but not for high-confidence experimental validation by itself.

## 8. References

1. Wikipedia. "Nucleic acid thermodynamics." https://en.wikipedia.org/wiki/Nucleic_acid_thermodynamics
2. Wikipedia. "Hybridization probe." https://en.wikipedia.org/wiki/Hybridization_probe
3. Wikipedia. "Fluorescence in situ hybridization." https://en.wikipedia.org/wiki/Fluorescence_in_situ_hybridization
4. Wikipedia. "DNA microarray." https://en.wikipedia.org/wiki/DNA_microarray
5. SantaLucia, J. (1998). "A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics." PNAS 95(4):1460-5.
6. Breslauer, K.J. et al. (1986). "Predicting DNA Duplex Stability from the Base Sequence." PNAS 83:3746-3750.
7. PREMIER Biosoft. "TaqMan® probe design tips." http://www.premierbiosoft.com/tech_notes/TaqMan.html (accessed 2026-06-24).
8. Applied Biosystems / Thermo Fisher Scientific. "Designing a TaqMan Gene Expression Assay." https://www.thermofisher.com/us/en/home/life-science/pcr/real-time-pcr/real-time-pcr-learning-center/gene-expression-analysis-real-time-pcr-information/designing-taqman-gene-expression-assay.html (accessed 2026-06-24).
9. ScienceDirect Topics. "TaqMan — an overview." https://www.sciencedirect.com/topics/biochemistry-genetics-and-molecular-biology/taqman (accessed 2026-06-24).
10. Applied Biosystems. Primer Express / "Designing TaqMan® MGB Probe and Primer Sets" guidelines: probe Tm 68–70 °C, no 5′ G, more C than G, avoid ≥ 4 G runs, G+C 30–80 %, probe ≥ 13 nt (via WebSearch extracts, 2026-10-01).
11. Cantor CR, Warshaw MM, Shapiro H (1970) Biopolymers 9:1059–1077; Warshaw MM, Tinoco I (1966) J Mol Biol 20:29–38 (ε260 nearest-neighbour tables as tabulated by OligoCalc / ATDBio; values confirmed via WebSearch extracts).
12. Tyagi S, Kramer FR (1996) Nat Biotechnol 14:303–308; Marras SAE et al. molecular-beacon design rules (stem 5–7 bp; probe and stem Tm 7–10 °C above the detection temperature).
13. Rozen S, Skaletsky H (2000) Primer3 on the WWW; Untergasser A et al. (2012) Nucleic Acids Res 40:e115; primer3 `libprimer3.cc` (raw.githubusercontent.com/primer3-org/primer3).
