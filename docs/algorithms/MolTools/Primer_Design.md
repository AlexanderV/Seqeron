# Primer Pair Design

| Field | Value |
|-------|-------|
| Algorithm Group | MolTools |
| Test Unit ID | PRIMER-DESIGN-001 |
| Related Projects | N/A |
| Implementation Status | Implemented (Primer3 pair search: product-size ranges, PRIMER_NUM_RETURN, pair objective, internal oligo; Primer3 thermodynamic structure screen by default) |
| Last Reviewed | 2026-10-01 |

## 1. Overview

Primer pair design selects forward and reverse oligonucleotides that can amplify a target DNA region by PCR. In this repository, primer design enumerates every candidate on either side of the target inside the search region (Primer3 SEQUENCE_TARGET / SEQUENCE_INCLUDED_REGION / PRIMER_PRODUCT_SIZE_RANGE), filters each by per-primer constraints, and then — exactly as Primer3 (`libprimer3.cc` `choose_pair_or_triple`) — returns the pair(s) with the lowest Primer3 pair penalty among all pairs meeting the pair constraints (product size and Tm, Tm agreement, secondary structure and primer-dimer avoidance, optionally an internal hybridization oligo). The Tm used everywhere in design is Primer3's default primer Tm, the scale on which the Primer3 Tm window 57–63 °C is defined.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

PCR primer design balances primer length, GC content, melting temperature, repetitive sequence content, and 3' end properties so that both primers bind specifically and amplify the desired product. The original document cites Primer3 and Addgene for standard design ranges and also notes that 3' terminal stability affects extension efficiency. Sources: Primer3 Manual, Addgene primer-design guidance, Wikipedia (Primer (molecular biology)), SantaLucia (1998).

### 2.2 Core Model

1. **Tm (per primer).** Primer3 `seqtm`/`oligotm` with PRIMER_TM_FORMULA = SantaLucia 1998 and
   PRIMER_SALT_CORRECTIONS = SantaLucia 1998 (`CalculateMeltingTemperaturePrimer3`):
   $[Mon]_{eq} = [Mon] + 120\sqrt{[Mg^{2+}] - [dNTP]}$ (mM; von Ahsen 2001),
   $\Delta S = \Delta S^\circ_{1M} + 0.368(N-1)\ln([Mon]_{eq}/1000)$,
   $T_m = \Delta H / (\Delta S + R\ln(C/4)) - 273.15$ ($C/1$ if self-complementary, $R = 1.987$),
   with the SantaLucia (1998) Table 2 NN and terminal-initiation terms as tabulated in `oligotm.c`;
   for $N > 36$ the `long_seq_tm` formula $81.5 + 16.6\log_{10}([Mon]_{eq}/1000) + 41\,GC/N - 600/N$.
   Conditions are `PrimerParameters.DnaConcentrationNanomolar` / `MonovalentMillimolar` / `DivalentMillimolar` /
   `DntpMillimolar` (PRIMER_DNA_CONC, PRIMER_SALT_MONOVALENT, PRIMER_SALT_DIVALENT, PRIMER_DNTP_CONC; null = Primer3
   defaults 50 nM oligo, 50 mM monovalent, 1.5 mM Mg²⁺, 0.6 mM dNTP). As in Primer3, the same primer conditions
   (`p_args`) drive the ntthal self-any / self-end / hairpin values (`create_thal_arg_holder(p_args)`), the pair
   compl-any / compl-end values and the product Tm (`long_seq_tm`); the internal oligo uses its own PRIMER_INTERNAL_*
   conditions (`Primer3ProbeSettings`, defaults 50 nM / 50 mM / 0 / 0). Primer3 `_pr_data_control` limits apply
   (salt and DNA concentration > 0, divalent ≥ 0; the library also rejects a negative dNTP and NaN/∞ —
   `ArgumentOutOfRangeException`). PRIMER_TM_FORMULA / PRIMER_SALT_CORRECTIONS stay fixed at SantaLucia 1998 (the
   Breslauer / Schildkraut / Owczarzy alternatives of `oligotm.c` are not ported).
2. **Per-primer penalty** (Primer3 `p_obj_fn`, default weights):
   $penalty = |T_m - OptimalTm| + |length - OptimalLength|$ (`CalculatePrimer3Penalty`); with non-zero
   PRIMER_WT_GC_PERCENT_GT/_LT (`PenaltyWeights.GcGt/GcLt`) the GC terms are taken around
   `PrimerParameters.OptimalGcPercent` (PRIMER_OPT_GC_PERCENT; null = 50 %, the manual's default — Primer3's code keeps
   it undefined and rejects GC weights without it). Internal oligo: `Primer3ProbeSettings.OptGcPercent` /
   `WeightGcPercentGt` / `WeightGcPercentLt` (PRIMER_INTERNAL_OPT_GC_PERCENT / _WT_GC_PERCENT_GT / _LT).
3. **Pair constraints** (`characterize_pair` order): product size in the current PRIMER_PRODUCT_SIZE_RANGE
   range (default 100–300 bp; ranges are tried in order, the next only when no further pair fits), product
   Tm within PRIMER_PRODUCT_MIN_TM/MAX_TM (if set), $|T_{m,f} - T_{m,r}| \le$ PRIMER_PAIR_MAX_DIFF_TM
   (library default 5 °C; Primer3 100), each primer's structure screen, pair compl-any/compl-end ntthal Tm
   ≤ 47 °C, the pair library mispriming score (below) and, with PRIMER_PICK_INTERNAL_OLIGO, an internal oligo
   strictly between the primers.
4. **Pair objective** (`obj_fn`): $W_{pr}(penalty_f + penalty_r) + W_{io}\,penalty_{io} + W_{\Delta Tm}|T_{m,f} - T_{m,r}|$
   $+ W_{any}\,g(compl\_any) + W_{end}\,g(compl\_end)$ $+ W_{tm<}(T_{opt} - T_{prod})^+ + W_{tm>}(T_{prod} - T_{opt})^+$
   $+ W_{size<}(S_{opt} - S)^+ + W_{size>}(S - S_{opt})^+$, with $g(x) = x - (T_{low} - 5 - 1)$ if $T_{low} - 5 \le x$,
   else $1/(T_{low} - 5 + 1 - x)$ ($T_{low}$ = lower primer Tm, `temp_cutoff` 5); defaults $W_{pr} = 1$, all
   others 0 (`Primer3PairWeights`). Product Tm (`long_seq_tm`):
   $81.5 + 16.6\log_{10}([Mon]_{eq}/1000) + 41\,GC/N - 600/N$ at the primer conditions.
5. **Selection:** the minimum-objective pair; ties within $10^{-6}$ broken as `compare_primer_pair` (left
   primer further 3′, right primer 5′ end further left, shorter left, shorter right). `DesignPrimerPairs`
   repeats the search PRIMER_NUM_RETURN times, removing each selected pair. By default primers may be reused
   (PRIMER_MIN_LEFT/RIGHT_THREE_PRIME_DISTANCE = −1, `PrimerPairOptions.MinLeftThreePrimeDistance` /
   `MinRightThreePrimeDistance`, shorthand `MinThreePrimeDistance` = PRIMER_MIN_THREE_PRIME_DISTANCE); with
   $d \ge 0$ a selected pair excludes from all later pairs every left (right) primer whose 3′ end lies fewer than $d$
   bases from the selected left (right) primer's 3′ end ($d = 0$: only the identical primer) — Primer3
   `choose_pair_or_triple` + `left/right_oligo_in_pair_overlaps_used_oligo` (left 3′ end = start + length − 1,
   right 3′ end = its leftmost top-strand base). PRIMER_INTERNAL_MIN_THREE_PRIME_DISTANCE acts only with
   SEQUENCE_INTERNAL_OVERLAP_JUNCTION_LIST (not modelled).
6. **Mispriming library** (PRIMER_MISPRIMING_LIBRARY, `PrimerParameters.MisprimingLibrary` =
   `PrimerMisprimingLibrary`, primer3-py `misprime_lib`; `libprimer3.c` `oligo_repeat_library_mispriming`,
   `pair_repeat_sim`, `p3_seq_lib.c`): entries name → sequence, an optional `*w` weight (0–100, `parse_seq_name`) in
   the name; sequences upper-cased, whitespace removed, IUPAC kept, other characters → N; Primer3 appends each
   entry's reverse complement as "reverse <name>". Primer score $w_i = weight_i \cdot align(primer, y_i)$ with dpal
   (+100 match, −100 mismatch, −25 N, −200 gap, max gap 1, score / 100, floored at 0; an entry shorter than 3 nt scores
   its length) anchored at the primer 3′ end (`DPAL_LOCAL_END`, port of `_dpal_long_nopath_maxgap1_local_end`); a left
   primer is aligned with entry $i$, a right primer (5′→3′) with its reverse complement. PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS
   (`LibraryAmbiguityCodesConsensus`, Primer3 default 0): 0 — IUPAC codes never align (`set_dpal_args`: INT_MIN) and,
   Primer3 quirk, a right primer uses unanchored `DPAL_LOCAL`; 1 — `dpal_set_ambiguity_code_matrix` (a code matches every
   base it represents). A primer fails when any $w_i$ > (short) PRIMER_MAX_LIBRARY_MISPRIMING (default 12); the reported
   score is Primer3's `repeat_sim.max` entry (first entry whose $w$ exceeds the integer part of the running maximum).
   Pair score = $\max_i \lfloor w_i^{left} + w_i^{right} \rfloor$ (≥ 0), failing above PRIMER_PAIR_MAX_LIBRARY_MISPRIMING
   (`PrimerPairOptions.MaxLibraryMispriming`, 24). Weights: PRIMER_WT_LIBRARY_MISPRIMING
   (`Primer3PenaltyWeights.LibraryMispriming`, per-primer `p_obj_fn`; primers are then scored at pick time) and
   PRIMER_PAIR_WT_LIBRARY_MISPRIMING (`Primer3PairWeights.LibraryMispriming`, `obj_fn` after the product-size terms);
   otherwise primers are scored lazily in `characterize_pair` (after the structure checks). Outputs:
   `PrimerCandidate.LibraryMispriming` / `LibraryMisprimingName`, `PrimerPairResult.LibraryMispriming` /
   `LibraryMisprimingName` (PRIMER_LEFT/RIGHT/PAIR_k_LIBRARY_MISPRIMING).
7. **Internal oligo** (`choose_internal_oligo`): the lowest-penalty oligo of the Primer3 internal-oligo list
   (`ProbeDesigner.DesignProbesPrimer3` rules, PRIMER_INTERNAL_* defaults; its self-any/self-end/hairpin checks
   postponed until the oligo is considered) with start > left primer 3′ end and end < right primer 5′ base.

`PrimerCandidate.Score` (100 − 2|len − opt| − 2|Tm − opt| − 0.5|GC − 50| − 5·homopolymer + 5 GC-clamp
bonus) is reported for information only and does not drive selection.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `DesignPrimers(...)` returns `IsValid = false` when either side has no valid candidates | The source returns an invalid `PrimerPairResult` when either best candidate is missing |
| INV-02 | Pair validity requires the product size in a PRIMER_PRODUCT_SIZE_RANGE range, `|Tm_f - Tm_r| <= MaxTmDifference` (unrounded Tm) and no primer-dimer (default: Primer3 `compl_any_th`/`compl_end_th` ≤ 47 °C, `CalculatePrimer3PairComplementarity`; `Heuristic` screen: `!HasPrimerDimer(...)`); if any such pair exists among the valid candidates, `IsValid = true` | Exhaustive pair search |
| INV-04 | The returned valid pair minimises the pair objective (`PairPenalty`; default `Forward.Penalty + Reverse.Penalty`) over all compatible pairs of the first product-size range that has one | Primer3 `choose_pair_or_triple` |
| INV-05 | `DesignPrimers(...)` = rank 0 of `DesignPrimerPairs(...)`; ranks are in non-decreasing `PairPenalty` order within a range and pairwise distinct | Primer3 `choose_pair_or_triple` while loop |
| INV-03 | `ProductSize = reverse.Position + reverse.Sequence.Length - forward.Position` | The source computes product size directly from the chosen candidates |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `template` | `DnaSequence` | required | Template DNA sequence |
| `targetStart` | `int` | required | Start of the target region | Must satisfy `targetStart >= 0` |
| `targetEnd` | `int` | required | Exclusive end of the target region (target = `[targetStart, targetEnd)`, Primer3 SEQUENCE_TARGET) | Must satisfy `targetEnd < template.Length` and `targetStart < targetEnd` |
| `parameters` | `PrimerParameters?` | `PrimerDesigner.DefaultParameters` | Primer design thresholds | Defaults are `18-25` bp length, `40-60%` GC, `57-63°C` Tm, `OptimalLength = 20`, `OptimalTm = 60`, `MaxHomopolymer = 4`, `MaxDinucleotideRepeats = 4`, `Avoid3PrimeGC = false`, and `Check3PrimeStability = true` (deprecated, no effect); `GcClamp` = 0, `MaxEndGc` = null (5) and `MaxEndStability` = null (100) are Primer3's defaults in both parameter sets; `PrimerDesigner.Primer3DefaultParameters` = Primer3's (18/20/27 nt, GC 20–80 %, poly-X 5, no dinucleotide limit) |
| `pairOptions` | `PrimerPairOptions?` | `PrimerPairOptions.Default` | Pair options: `ProductSizeRanges` (PRIMER_PRODUCT_SIZE_RANGE, default 100–300), `MaxTmDifference` (PRIMER_PAIR_MAX_DIFF_TM, default 5; `Primer3Defaults`: 100), `NumReturn` (PRIMER_NUM_RETURN, 5), `IncludedRegion` (SEQUENCE_INCLUDED_REGION), `ProductOptSize`/`ProductOptTm`/`ProductMinTm`/`ProductMaxTm`, `Weights` (PRIMER_PAIR_WT_*), `PickInternalOligo` + `InternalOligo` (PRIMER_PICK_INTERNAL_OLIGO, PRIMER_INTERNAL_*) | Primer3 `_pr_data_control` errors throw `ArgumentException` (weight without optimum, PRIMER_MAX_SIZE or PRIMER_INTERNAL_MAX_SIZE > min product size, NUM_RETURN < 1, target outside the included region) |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Forward` | `PrimerCandidate?` | Selected forward primer or `null` |
| `Reverse` | `PrimerCandidate?` | Selected reverse primer or `null` |
| `IsValid` | `bool` | Pair validity flag |
| `Message` | `string` | Result explanation |
| `ProductSize` | `int` | Predicted amplicon size |
| `PairPenalty` | `double?` | PRIMER_PAIR_k_PENALTY (valid pair) |
| `ProductTm` | `double?` | PRIMER_PAIR_k_PRODUCT_TM (`long_seq_tm`) |
| `ComplAnyTh` / `ComplEndTh` | `double?` | PRIMER_PAIR_k_COMPL_ANY_TH / _COMPL_END_TH (thermodynamic screen) |
| `InternalOligo` | `ProbeDesigner.Primer3Probe?` | PRIMER_INTERNAL_k_* (with `PickInternalOligo`) |

`DesignPrimerPairs(...)` returns up to `NumReturn` such valid results, best first.

### 3.3 Preconditions and Validation

`DesignPrimers(...)` throws `ArgumentNullException` for a null template and `ArgumentException` when the requested target region or the options are invalid. Forward candidates end at or before `targetStart` and reverse candidates start at or after `targetEnd`, both inside the included region (default: whole template); only positions that can form a product inside some product-size range are enumerated (Primer3 `pick_primer_range` additionally drops left primers starting past `n − min product` and right primers ending before `min product`; the further restriction is result-equivalent). Reverse-primer candidates are reverse-complemented before evaluation so that they are scored in primer orientation.

## 4. Algorithm

### 4.1 High-Level Steps

1. Search region: included region, primers outside the target, product within the product-size ranges.
2. Enumerate all candidate primers within the configured length range (reverse candidates reverse-complemented).
3. Evaluate each candidate for GC content, Tm, homopolymers, dinucleotide repeats and 3′ stability (structure screen postponed).
4. Sort each side by Primer3 penalty; scan reverse × forward candidates with Primer3's pruning (stop once `W_pr·(penalty_f + penalty_r)` exceeds the best pair found; stop at objective 0); characterize each pair once (cached): product size, product Tm, ΔTm, lazy primer structure, pair complementarity, internal oligo, objective.
5. Keep the best pair; for `DesignPrimerPairs` mark it used and repeat until `NumReturn` pairs or no pair is left (then the next product-size range). If none exists, `DesignPrimers` returns the individually best primers with `IsValid = false` and a message naming the violated constraint.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Parameter ranges documented in the original file and current source:

| Parameter | Min | Optimal | Max | Notes |
|-----------|-----|---------|-----|-------|
| Length (bp) | 18 | 20 | 25 | Matches current source defaults |
| GC Content (%) | 40 | 50 | 60 | Matches current source defaults |
| Melting Temp (°C) | 55 | 60 | 65 | Original document summary; current source defaults narrow this to `57-63` |
| Homopolymer Run | N/A | N/A | 4 | Current source default |
| Dinucleotide Repeats | N/A | N/A | 4 | Current source default |
| `GcClamp` (PRIMER_GC_CLAMP) | 0 | N/A | `MinLength` | Number of consecutive G/C required at the 3′ end (Primer3 default 0; > PRIMER_MIN_SIZE throws) |
| `MaxEndGc` (PRIMER_MAX_END_GC) | 0 | N/A | 5 | Max G/C among the five 3′-most bases; checked only when < 5 (default 5) |
| `MaxEndStability` (PRIMER_MAX_END_STABILITY) | 0 | N/A | 100 | Fails when `end_stability` = −`Calculate3PrimeStability` > limit (default 100; an ACGT pentamer reaches at most 6.86, GCGCG/CGCGC) |
| `Avoid3PrimeGC` | N/A | N/A | `false` | Deprecated library rule (not Primer3): despite the name it *requires* at least one `G`/`C` in the last two bases; kept for source compatibility — use `GcClamp` / `MaxEndGc` |
| `MisprimingLibrary` (PRIMER_MISPRIMING_LIBRARY) | N/A | N/A | null | name → sequence entries (`PrimerMisprimingLibrary`); null = no library check |
| `MaxLibraryMispriming` (PRIMER_MAX_LIBRARY_MISPRIMING) | N/A | N/A | 12 | Truncated to a C `short`; > 32767 rejected in alignment mode (`_pr_data_control`) |
| `PrimerPairOptions.MaxLibraryMispriming` (PRIMER_PAIR_MAX_LIBRARY_MISPRIMING) | N/A | N/A | 24 | Pair score = max over entries of ⌊left + right⌋ |
| `LibraryAmbiguityCodesConsensus` (PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS) | N/A | N/A | false (0) | Primer3 v2 default (`pr_set_default_global_args_2`); 1 makes IUPAC codes match their bases |
| `Check3PrimeStability` | N/A | N/A | `true` | Deprecated, no effect: it gated ΔG(3′ pentamer) < −9 kcal/mol, which no ACGT pentamer reaches (minimum −6.86), so it never fired; the Primer3 limit is `MaxEndStability` |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `DesignPrimers` / `DesignPrimerPairs` | `O(c log c + p)` | `O(c + p)` | `c` candidates (positions × lengths in the product window), `p` characterized pairs (each up to 5 ntthal alignments; Primer3's pruning bounds `p`) |
| `EvaluatePrimer` | `O(n)` | `O(1)` | Per-primer scan and helper calculations |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [PrimerDesigner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs)

- `PrimerDesigner.DesignPrimers(DnaSequence, int, int, PrimerParameters?, PrimerPairOptions?)`: Designs and validates the best primer pair around a target region.
- `PrimerDesigner.DesignPrimerPairs(DnaSequence, int, int, PrimerParameters?, PrimerPairOptions?)`: PRIMER_NUM_RETURN ranked pairs.
- `PrimerDesigner.CalculateProductMeltingTemperaturePrimer3(string, ...)`: Primer3 `long_seq_tm` product Tm.
- `PrimerDesigner.EvaluatePrimer(string, int, bool, PrimerParameters?)`: Scores a single primer candidate.
- `PrimerDesigner.CalculateLibraryMispriming(string, bool, PrimerMisprimingLibrary, bool)`: Primer3 library mispriming score of one primer ([PrimerDesigner.MisprimingLibrary.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.MisprimingLibrary.cs)).
- `PrimerDesigner.CalculateMeltingTemperaturePrimer3(string, ...)`: Primer3-default primer Tm used by design.
- `PrimerDesigner.CalculatePrimer3Penalty(...)`: Primer3 per-primer penalty used for ranking.
- `PrimerDesigner.CalculatePrimerScore(...)` (private): informational heuristic score.

### 5.2 Current Behavior

Forward primers are taken directly from the template; reverse primers are reverse-complemented before evaluation, and their `Position` is the leftmost template coordinate of the binding site. Per-primer hard constraints: length, GC%, Primer3 3′-end checks (PRIMER_GC_CLAMP / PRIMER_MAX_END_GC / PRIMER_MAX_END_STABILITY, `calc_and_check_oligo_features`), Primer3-default Tm window, homopolymer, dinucleotide repeat, secondary structure (default `PrimerStructureScreen.Primer3Thermodynamic`: Primer3 ntthal self-any / self-end / hairpin Tm ≤ `MaxStructureTm` = 47 °C, evaluated lazily in the pair loop like Primer3's `characterize_pair`; `Heuristic`: `HasHairpinPotential`), the deprecated library rule `Avoid3PrimeGC` (≥ 1 G/C in the last two bases), and no non-ACGT base (Primer3 PRIMER_MAX_NS_ACCEPTED = 0). Pair selection is the Primer3 pair search described in §2.2; the internal oligo reuses the `ProbeDesigner` Primer3 internal-oligo list (`EnumeratePrimer3InternalOligos`, shared with `DesignProbesPrimer3`).

### 5.3 Conformance to Theory / Spec

**Implemented (verified against primer3-py 2.3.1):**

- Primer3 default Tm: bit-identical to `primer3.calc_tm` (max |Δ| = 0 over 3 000 random 2–45-mers, incl. self-complementary and > 36 nt).
- Pair search with product-size ranges, PRIMER_NUM_RETURN, pair objective and internal oligo (audit round 1, 2026-10-01): with `Primer3DefaultParameters` + `PrimerPairOptions.Primer3Defaults`, `DesignPrimerPairs` reproduced `primer3.design_primers` (only SEQUENCE_TEMPLATE/SEQUENCE_TARGET set) for 1000/1000 random templates (150–700 bp, 4990 pairs; left/right position+length, PRIMER_PAIR_k_PENALTY, _PRODUCT_TM, _COMPL_ANY_TH, _COMPL_END_TH, _PRODUCT_SIZE, ranks 0–4, |Δ| ≤ 1e-9); with random non-default PRIMER_PAIR_WT_* / PRODUCT_OPT_* / PRODUCT_MIN/MAX_TM / PAIR_MAX_DIFF_TM / multi-range PRIMER_PRODUCT_SIZE_RANGE / SEQUENCE_INCLUDED_REGION settings and with PRIMER_PICK_INTERNAL_OLIGO = 1 (+ PRIMER_PAIR_WT_IO_PENALTY; PRIMER_INTERNAL_k position, penalty, Tm, self-any/self-end/hairpin Tm) the counts are in the PRIMER-DESIGN-001 F-entry of `docs/Validation/review-2026-09/B07.md`.
- Primer3 alignment mode (audit round 2, A1, 2026-10-02): with `PrimerStructureScreen.Primer3Alignment` (PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT = 0: dpal self_any/self_end ≤ PRIMER_MAX_SELF_ANY/_END 8/3, pair compl_any/compl_end ≤ PRIMER_PAIR_MAX_COMPL_ANY/_END 8/3, internal oligo ≤ 12/12, PRIMER_WT_SELF_ANY/_END, PRIMER_PAIR_WT_COMPL_ANY/_END) `DesignPrimerPairs` reproduced `primer3.design_primers` for 2000/2000 random templates (9135 pairs: positions, PRIMER_PAIR_k_PENALTY, _COMPL_ANY, _COMPL_END, PRIMER_LEFT/RIGHT_k_SELF_ANY/_SELF_END/_PENALTY, product Tm, PRIMER_INTERNAL_k_* ; 1000 templates with random non-default limits/weights/internal oligos); per-primer thermodynamic weights PRIMER_WT_SELF_ANY_TH/_SELF_END_TH/_HAIRPIN_TH (`PrimerParameters.PenaltyWeights`) 60/60 templates (280 pairs) identical.
- Reaction conditions and GC optimum (audit round 3, A3-1, 2026-10-02): random PRIMER_SALT_MONOVALENT 10–200 mM, PRIMER_SALT_DIVALENT {0, 0–6} mM, PRIMER_DNTP_CONC {0, 0–2} mM, PRIMER_DNA_CONC 10–500 nM, PRIMER_OPT_GC_PERCENT 25–75 with PRIMER_WT_GC_PERCENT_GT/_LT ∈ {0, 0.25–1}, the PRIMER_INTERNAL_* counterparts (+ PRIMER_INTERNAL_OPT_GC_PERCENT / _WT_GC_PERCENT_GT/_LT) with PRIMER_PICK_INTERNAL_OLIGO ∈ {0, 1}, both PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT modes and an optional product-Tm objective: 2500/2500 templates (12 440 pairs) identical on positions, PRIMER_PAIR/LEFT/RIGHT_k_PENALTY, _TM, PRIMER_PAIR_k_PRODUCT_TM, the self-any/self-end/hairpin (or dpal) and pair compl values and PRIMER_INTERNAL_k position/penalty/Tm (relative |Δ| ≤ 1e-9); 1/2500 before (fixed 50 mM / 1.5 mM / 0.6 mM / 50 nM, GC optimum 50).
- Primer3 per-primer penalty and pair search (before the audit, ±200 bp flanks): `DesignPrimers` returned exactly Primer3's `PRIMER_LEFT_0`/`PRIMER_RIGHT_0` in 553/553 random templates (3 seeds × 300) where Primer3's best pair also passes this library's extra screens (settings mirroring `DefaultParameters`, thermodynamic structure limits disabled).

**Deviations from Primer3 defaults (documented; `Primer3DefaultParameters` / `PrimerPairOptions.Primer3Defaults` restore them):** length 18–25 (Primer3 18–27), GC 40–60 % (20–80 %), poly-X 4 (5), pair ΔTm ≤ 5 °C (100), dinucleotide-repeat limit (no Primer3 equivalent). The product-size range is Primer3's (default 100–300 bp; before audit round 1 the ±200 bp flanks bounded the product instead). Structure limits are Primer3's default thermodynamic ones (PRIMER-STRUCT-001): with them `DesignPrimers` returned primer3-py's pair (same settings) in 574/600 random templates; the 26 differences all trace to ntthal engine values (PRIMER-DIMER-001 / PRIMER-HAIRPIN-001). After PRIMER-DIMER-001 (dimer engine bit-exact) a re-run on 1800 random templates (seeds 1–9 × 200) gave 1733 identical; all 67 differences trace to `calc_hairpin` values only (PRIMER-HAIRPIN-001). After PRIMER-HAIRPIN-001 (hairpin engine bit-exact) the same 1800 templates are 1800/1800 identical (primers and pair penalty).

- 3′-end checks and 3′ distance (audit round 3, A3-6 + A3-7, 2026-10-02): random 150–600-bp templates (GC bias 35–60 %), PRIMER_GC_CLAMP ∈ {0–3}, PRIMER_MAX_END_GC ∈ {0–5}, PRIMER_MAX_END_STABILITY ∈ {4–9, 100}, PRIMER_MIN_THREE_PRIME_DISTANCE or PRIMER_MIN_LEFT/RIGHT_THREE_PRIME_DISTANCE ∈ {−1, 0, 1, 2, 3, 5, 10, 20}, PRIMER_NUM_RETURN 5–10, PRIMER_WT_END_STABILITY ∈ {0, 0.5}, both alignment modes, PRIMER_PICK_INTERNAL_OLIGO ∈ {0, 1}: 2200/2200 templates (9590 pairs) identical on left/right start + length, PRIMER_PAIR/LEFT/RIGHT_k_PENALTY, _END_STABILITY and PRIMER_INTERNAL_k position (|Δ| ≤ 1e-9); without the new options 36/2000.

- Mispriming library (audit round 3, A3-3 part 1, 2026-10-02): random 150–500-bp templates with 1–8-entry random libraries (55 % template fragments, forward or reverse-complemented, with point mutations and IUPAC codes; random IUPAC-containing sequences; 1–2-nt entries; N-rich entries; `*w` weights 0–10), PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS ∈ {0, 1}, PRIMER_MAX_LIBRARY_MISPRIMING ∈ {6, 8, 10, 12, 12.9, 15, 20}, PRIMER_PAIR_MAX_LIBRARY_MISPRIMING ∈ {12–30}, PRIMER_WT_LIBRARY_MISPRIMING ∈ {0, 0.1, 0.5, 1}, PRIMER_PAIR_WT_LIBRARY_MISPRIMING ∈ {0, 0.2, 1}, both alignment modes, PRIMER_PICK_INTERNAL_OLIGO ∈ {0, 1}: 1240/1240 templates (4693 pairs; primers rejected by the library in 1005 templates, pairs in 202) identical to primer3-py 2.3.1 `design_primers(misprime_lib=…)` on left/right start + length, PRIMER_PAIR/LEFT/RIGHT_k_PENALTY (|Δ| ≤ 1e-9), PRIMER_LEFT/RIGHT/PAIR_k_LIBRARY_MISPRIMING score and entry name, and PRIMER_INTERNAL_k position.

**Not implemented:** internal-oligo mishybridization library (PRIMER_INTERNAL_MISHYB_LIBRARY, A3-3 part 2) / template mispriming, position penalties (PRIMER_INSIDE/OUTSIDE_PENALTY), sequence quality, PRIMER_INTERNAL_MIN_THREE_PRIME_DISTANCE (needs SEQUENCE_INTERNAL_OVERLAP_JUNCTION_LIST), genome-wide specificity.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Invalid target region | Throws `ArgumentException` | Explicit source guard |
| No valid forward or reverse candidates | Returns an invalid `PrimerPairResult` with null candidates | Explicit fallback in source |
| No pair within `MaxTmDifference` (default 5 °C) | Returns the individually best primers with `IsValid = false` | PRIMER_PAIR_MAX_DIFF_TM |
| No pair with a product size in any range | `IsValid = false`, message "No primer pair with a product size in … bp (PRIMER_PRODUCT_SIZE_RANGE)." | Primer3 "unacceptable product size" |
| Included region (default: template) shorter than the smallest product size | `IsValid = false`, message "SEQUENCE_INCLUDED_REGION length < min PRIMER_PRODUCT_SIZE_RANGE" | Primer3 `_pr_data_control` per-sequence error (primer3-py raises it) |
| `PickInternalOligo` and no acceptable oligo between the primers | Pair fails (Primer3 "no internal oligo") | `choose_internal_oligo` |
| Primer-dimer detected for every pair | Returns `IsValid = false` | Pair compatibility requires no dimer signal |
| Non-ACGT base in a candidate | Candidate invalid (Tm 0, issue "Tm not computable") | Primer3 PRIMER_MAX_NS_ACCEPTED = 0 |
| PRIMER_WT_LIBRARY_MISPRIMING / PRIMER_PAIR_WT_LIBRARY_MISPRIMING ≠ 0 without a library | `ArgumentException` | Primer3 `_pr_data_control` "Mispriming score is part of objective function, but mispriming library is not defined" |
| Library entry with an empty sequence or an illegal `*weight` (missing, < 0, > 100) | `ArgumentException` | `add_seq_to_seq_lib` / `parse_seq_name` (primer3-py raises OSError) |
| Library entry with a non-IUPAC character | Character becomes N, `PrimerMisprimingLibrary.Warnings` | `upcase_and_check_char` (primer3-py 2.3.1 aborts here: it passes a NULL `errfrag` to the warning) |

### 6.2 Limitations

There is no template-mispriming check and no internal-oligo mishybridization library (the primer mispriming library
is Primer3's, §2.2 item 6). The secondary-structure screen is Primer3's thermodynamic one by default (the sequence-only screen is available as `PrimerStructureScreen.Heuristic`).

## 7. Examples and Related Material

### 7.2 Applications and Use Cases (Optional)

Related material called out in the original document:

- `PRIMER-TM-001`: Melting temperature calculation (prerequisite).
- `PRIMER-STRUCT-001`: Hairpin and dimer detection (used in evaluation).

## 8. References

1. [Primer (molecular biology)](https://en.wikipedia.org/wiki/Primer_(molecular_biology)) - Standard primer design criteria.
2. [How to Design a Primer](https://www.addgene.org/protocols/primer-design/) - Addgene protocol guidance.
3. [primer3.org/manual.html](https://primer3.org/manual.html) - Primer3 manual.
4. SantaLucia JR (1998). "A unified view of polymer, dumbbell and oligonucleotide DNA nearest-neighbor thermodynamics", PNAS 95:1460-65.
5. Untergasser A et al. (2012). "Primer3 — new capabilities and interfaces", NAR 40(15):e115.
6. Primer3 source (primer3-org/primer3, `src/oligotm.c`: `oligotm`, `seqtm`, `long_seq_tm`, `divalent_to_monovalent`; `src/libprimer3.cc`: `make_detection_primer_lists`, `pick_primer_range`, `choose_pair_or_triple`, `characterize_pair`, `obj_fn`, `choose_internal_oligo`, `primer_rec_comp`, `compare_primer_pair`, `p_obj_fn`, `_pr_data_control`, `pr_set_default_global_args_1`).
7. von Ahsen N, Wittwer CT, Schütz E (2001). Clin Chem 47:1956-61 (divalent→monovalent equivalence).
