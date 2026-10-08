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
   `PrimerParameters.OptimalGcPercent` (PRIMER_OPT_GC_PERCENT; null = undefined, as in Primer3's code —
   `DEFAULT_OPT_GC_PERCENT = PR_UNDEFINED_INT_OPT`; the manual's "50" is not applied). A non-zero GC weight without the
   optimum is Primer3's `_pr_data_control` error ("Primer GC content is part of objective function while optimum
   gc_content is not defined", `ArgumentException`). Internal oligo: `Primer3ProbeSettings.OptGcPercent` /
   `WeightGcPercentGt` / `WeightGcPercentLt` (PRIMER_INTERNAL_OPT_GC_PERCENT / _WT_GC_PERCENT_GT / _LT; "Hyb probe GC
   content is part of objective function while optimum gc_content is not defined" — checked whether or not an internal
   oligo is picked).
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
   right 3′ end = its leftmost top-strand base). PRIMER_INTERNAL_MIN_THREE_PRIME_DISTANCE (which Primer3's
   shorthand also sets, `read_boulder.c`) is read by `choose_pair_or_triple` only when
   SEQUENCE_INTERNAL_OVERLAP_JUNCTION_LIST is non-empty (`pick_internal_oligo && intl_overlap_junctions_count > 0`);
   this unit takes no junction list (§5.3, outside its scope), so it is inert here exactly as in Primer3 without
   that list — no separate option exists (primer3-py 2.3.1, PRIMER_PICK_INTERNAL_OLIGO = 1: 1000/1000 designs on 200
   random 200–500-bp templates, 1319 pairs, identical for PRIMER_INTERNAL_MIN_THREE_PRIME_DISTANCE ∈ {−1, 0, 1, 5, 20}
   and unset; audit round 5, A5-3).
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
   postponed until the oligo is considered) with start > left primer 3′ end and end < right primer 5′ base. With a
   mishybridization library (PRIMER_INTERNAL_MISHYB_LIBRARY, `Primer3ProbeSettings.MishybLibrary`, primer3-py
   `mishyb_lib`; same `PrimerMisprimingLibrary` type and dpal port as item 6) the oligo is scored by
   `oligo_repeat_library_mispriming` OT_INTL: $w_i = weight_i \cdot align(oligo, y_i)$ with the **unanchored**
   `DPAL_LOCAL` (`local` / `local_ambig` — PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS is global, so the pair search uses
   `PrimerParameters.LibraryAmbiguityCodesConsensus`); an oligo with any $w_i$ > (short) PRIMER_INTERNAL_MAX_LIBRARY_MISHYB
   (`MaxLibraryMishyb`, 12) is rejected. With PRIMER_INTERNAL_WT_LIBRARY_MISHYB (`WeightLibraryMishyb`, default 0) ≠ 0
   oligos are scored while the list is built (a rejection is a "five-prime problem" ending the 5′ extension; the
   weighted score enters the internal-oligo `p_obj_fn` and, via PRIMER_PAIR_WT_IO_PENALTY, the pair objective),
   otherwise in `choose_internal_oligo` after the postponed structure checks. Output: `Primer3Probe.LibraryMishyb` /
   `LibraryMishybName` (PRIMER_INTERNAL_k_LIBRARY_MISHYB; primer3-py key `…_LIBRARY_MISPRIMING`).
8. **Template mispriming** (`oligo_template_mispriming`; `PrimerDesigner.CalculateTemplateMispriming`,
   PrimerDesigner.TemplateMispriming.cs): each primer (5′→3′ as synthesised) is aligned with the template strand it
   anneals to, once 5′ and once 3′ of its own site ($T$ = the larger), and with the whole opposite strand ($T_r$).
   PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT = 0 (`PrimerParameters.ThermodynamicTemplateAlignment` false, default):
   dpal `DPAL_LOCAL_END` (the item-6 port; a segment < 3 nt scores its length); = 1: the ntthal THAL_END1 Tm
   (`use_end_for_th_template_mispriming`, `NtthalDimer`, primer reaction conditions; empty segment / no structure /
   negative Tm → 0; Primer3 swaps the strand roles of the two thermodynamic calls, reproduced; template ≤ 10000 nt).
   PRIMER_LEFT/RIGHT_k_TEMPLATE_MISPRIMING[_TH] = max($T$, $T_r$) (`PrimerCandidate.TemplateMispriming`); a primer fails
   when PRIMER_MAX_TEMPLATE_MISPRIMING[_TH] (`MaxTemplateMispriming[Th]`, default −100) ≥ 0 and is exceeded. Scored
   at pick time when PRIMER_WT_TEMPLATE_MISPRIMING[_TH] ≠ 0 (the weighted term enters `p_obj_fn`: linear, or the 5 °C
   `temp_cutoff` rule in thermodynamic mode), otherwise in `characterize_pair` — only when the primer's library scores
   were not computed at pick time (no library or PRIMER_WT_LIBRARY_MISPRIMING = 0), exactly as Primer3. Pair:
   PRIMER_PAIR_k_TEMPLATE_MISPRIMING[_TH] = max($T^L + T_r^R$, $T_r^L + T^R$) (`PrimerPairResult.TemplateMispriming`),
   computed when PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING[_TH] ≥ 0 or PRIMER_PAIR_WT_TEMPLATE_MISPRIMING[_TH] > 0
   (`PrimerPairOptions.MaxTemplateMispriming[Th]`, `Primer3PairWeights.TemplateMispriming[Th]`); alignment mode fails
   when the limit is ≥ 0 and exceeded, thermodynamic mode when the limit is **non-zero** and exceeded (Primer3's
   test — with the default −100 and a pair weight every pair fails; 0 = no limit); the weighted value enters `obj_fn`
   (thermodynamic: `temp_cutoff` rule relative to the lower primer Tm). Only the weights/limits of the active mode apply.
9. **Fraction bound** (`oligotm` with PRIMER_ANNEALING_TEMP; `PrimerDesigner.CalculateFractionBoundPrimer3`,
   PrimerDesigner.BoundAndPosition.cs): with `PrimerParameters.AnnealingTemperature` $T_a$ > 0 every primer (and the
   internal oligo — one global Primer3 setting) gets, from the same SantaLucia 1998 ΔH/ΔS (SantaLucia salt correction)
   as its Tm, $\Delta G = \Delta H - (T_a + 273.15)\Delta S$, $K = e^{-\Delta G / (1.987 (T_a + 273.15))}$,
   bound $= 100 / (1 + \sqrt{1/((C/x) K)})$ % ($x = 4\cdot10^9$, $10^9$ for a self-complementary oligo; C in nM;
   an oligo > 36 nt or with a non-ACGT base has none — Primer3's OLIGOTM_ERROR −999999.9999). Checked after the Tm
   (`MinBound` / `MaxBound`, PRIMER_MIN/MAX_BOUND −10 / 110 %; not a five-prime problem) and penalised in `p_obj_fn`
   by PRIMER_WT_BOUND_GT / _LT (`Primer3PenaltyWeights.BoundGt` / `BoundLt`) around PRIMER_OPT_BOUND (`OptBound`, 97 %).
   Primer terms only when $T_a$ > 0; the internal-oligo terms (`Primer3ProbeSettings.WeightBoundGt` / `WeightBoundLt`,
   `MinBound` / `MaxBound` / `OptBound`) are not gated, so without $T_a$ PRIMER_INTERNAL_WT_BOUND_LT adds
   w·(opt + 999999.9999), exactly as Primer3. Output: `PrimerCandidate.Bound`, `Primer3Probe.Bound`
   (PRIMER_LEFT/RIGHT/INTERNAL_k_BOUND, when $T_a$ > 0).
10. **Position penalty** (`compute_position_penalty`; `PrimerDesigner.CalculatePositionPenaltyPrimer3`): with
   `PrimerPairOptions.InsidePenalty` / `OutsidePenalty` other than Primer3's −1 / 0 (PRIMER_INSIDE/OUTSIDE_PENALTY)
   primers are no longer kept off the target: `make_detection_primer_lists` searches the whole included region, a
   left primer is allowed while its 3′ end ≤ the target's last base (a right primer while its 3′ end — leftmost
   top-strand base — ≥ the target's first base), a 3′ end inside the target costs `InsidePenalty` × (distance from the
   near target end + 1), one outside `OutsidePenalty` × (bases between it and the target), × PRIMER_WT_POS_PENALTY
   (`Primer3PenaltyWeights.PositionPenalty`, default 1) in `p_obj_fn`, and a pair must still span the target
   (`pair_spans_target`: left 3′ end < right 3′ end). Output `PrimerCandidate.PositionPenalty`
   (PRIMER_LEFT/RIGHT_k_POSITION_PENALTY). Primer3 multiplies with the value as given, so changing only the outside
   penalty leaves the inside multiplier at −1 (negative penalties); a negative pair penalty makes Primer3 abort
   (`obj_fn` `PR_ASSERT(sum >= 0.0)`) — the library throws `InvalidOperationException`.
11. **Sequence quality** (`sequence_quality_is_ok`; `PrimerDesigner.CalculateSequenceQualityPrimer3`,
   PrimerDesigner.SequenceQuality.cs): with SEQUENCE_QUALITY (`PrimerPairOptions.SequenceQuality`, one integer per
   template base; `DesignProbesPrimer3(…, sequenceQuality)`) every oligo gets `seq_quality` = min(PRIMER_QUALITY_RANGE_MAX,
   qualities over the oligo) and `seq_end_quality` = the same over its five 3′-most bases (left primer / internal oligo:
   the last five template positions; right primer: the first five). Checked after the end-GC check: `seq_quality` <
   PRIMER_MIN_QUALITY (`PrimerParameters.MinQuality`; internal oligo `Primer3ProbeSettings.MinQuality`) or, for primers
   only, `seq_end_quality` < PRIMER_MIN_END_QUALITY (`MinEndQuality`) rejects the oligo (both five-prime problems; a
   5′ extension keeps the low base, so the result equals Primer3's extension break). `p_obj_fn` adds
   PRIMER_WT_SEQ_QUAL × (PRIMER_QUALITY_RANGE_MAX − `seq_quality`) (`Primer3PenaltyWeights.SequenceQuality`; internal
   oligo `WeightSequenceQuality`, the last internal-oligo term); PRIMER_[INTERNAL_]WT_END_QUAL is parsed by Primer3 but
   never read by `p_obj_fn` — accepted with no effect (`EndQuality` / `WeightEndQuality`). Output
   `PrimerCandidate.MinSequenceQuality` / `Primer3Probe.MinSequenceQuality` (PRIMER_LEFT/RIGHT/INTERNAL_k_MIN_SEQ_QUALITY).
   `_pr_data_control`: SEQUENCE_QUALITY length ≠ template length ("Error in sequence quality data"), a non-zero
   PRIMER_[INTERNAL_]MIN_QUALITY without quality ("Sequence quality data missing"), with quality a non-zero minimum
   outside [PRIMER_QUALITY_RANGE_MIN, PRIMER_QUALITY_RANGE_MAX] (`QualityRangeMin` / `QualityRangeMax`, 0 / 100, global —
   also used for the internal oligo) or a quality value outside it ("Sequence quality score out of range"), a non-zero
   PRIMER_[INTERNAL_]WT_SEQ_QUAL without quality — all `ArgumentException`; the PRIMER_INTERNAL_* values are checked even
   without PRIMER_PICK_INTERNAL_OLIGO, as in Primer3. `EvaluatePrimer` has no template, hence no quality: a non-zero
   `MinQuality` or quality weight throws there. With PRIMER_PICK_INTERNAL_OLIGO and no acceptable internal oligo at all,
   Primer3's `make_internal_oligo_list` fails and no pair is examined (reproduced: no pairs). A non-zero
   PRIMER_PAIR_WT_IO_PENALTY without PRIMER_PICK_INTERNAL_OLIGO throws ("Internal oligo quality is part of objective
   function while internal oligo choice is not required", audit A3-23).
12. **Template masking** (Primer3 masker, Kõressaar et al. 2018; masker.c; `PrimerDesigner.MaskTemplatePrimer3`,
   `CalculateMaskFailureRatePrimer3`, PrimerDesigner.Masking.cs): with `PrimerPairOptions.MaskTemplate`
   (PRIMER_MASK_TEMPLATE) and `MaskKmerLists` (`PrimerMaskingKmerLists`: caller-supplied 11-mer / 16-mer genome counts,
   or the GenomeTester4 files `{PRIMER_MASK_KMERLIST_PATH}{PRIMER_MASK_KMERLIST_PREFIX}_11.list` / `_16.list` via
   `FromGenomeTester4Files`), the predicted PCR failure rate of a primer ending with a 16-nt window is
   $F = e^{s-4.336}/(1+e^{s-4.336})$, $s = 0.1772\ln n_{11} + 0.239\ln n_{16}$ ($n_k$ = count of the 3′-terminal k-mer,
   else of its reverse complement, else 1; $F = 0$ when $s = 0$). The included region is masked once
   (`read_and_mask_sequence`, both strands separately, soft masking): every window whose forward score exceeds
   PRIMER_MASK_FAILURE_RATE (`MaskFailureRate`, 0.1; 0 = none) lower-cases its last base and the
   PRIMER_MASK_5P_DIRECTION − 1 (`MaskFivePrimeDirection`, 1) bases before it and the PRIMER_MASK_3P_DIRECTION
   (`MaskThreePrimeDirection`, 0) bases after it on the forward copy; a reverse score above the rate masks the window's
   first base, the 5P − 1 bases after it and the 3P bases before it on the reverse copy (masker.c's 5000-character ring
   buffer is reproduced). A left primer whose 3′ base is masked on the forward copy, or a right primer whose 3′ base (its
   leftmost template base) is masked on the reverse copy, is rejected (`is_lowercase_masked`, "3' end overlaps masked
   sequence"); every primer ≥ 16 nt gets `failure_rate` from its last 16 nt (`PrimerCandidate.MaskFailureRate`), weighted by
   PRIMER_WT_MASK_FAILURE_RATE (`Primer3PenaltyWeights.MaskFailureRate`, default 0) after the size terms of `p_obj_fn`.
   Internal oligos are not masked (only the lower-case check below applies to them); without PRIMER_MASK_TEMPLATE the
   failure rate is 0.
   PRIMER_MASK_TEMPLATE without lists throws (primer3-py: "masking template chosen, but path to
   PRIMER_MASK_KMERLIST_PATH not specified"); negative directions or PRIMER_MASK_3P_DIRECTION > 4984 (undefined in
   masker.c's unsigned counters / buffer) throw `ArgumentOutOfRangeException`.
13. **Lower-case masking** (PRIMER_LOWERCASE_MASKING, `PrimerPairOptions.LowercaseMasking`; `calc_and_check_oligo_features` →
   `is_lowercase_masked`): with a case-preserving template — the `string` overloads `DesignPrimers(string, …)` /
   `DesignPrimerPairs(string, …)` (validated as `DnaSequence`, A/C/G/T only) — a left primer or internal oligo whose 3′ base
   (rightmost template base), or a right primer whose 3′ base (leftmost template base), is a lower-case a/c/g/t of the
   template as given (`trimmed_orig_seq`) is rejected ("3' end overlaps lower-case masked sequence"); lower case elsewhere
   in the oligo is accepted. PRIMER_MASK_TEMPLATE implies it (primer3-py sets `lowercase_masking = mask_template`): primers
   are then checked on the masked copies, which keep the input's lower case (masker.c COND0), internal oligos on the
   case-preserving template. `DesignProbesPrimer3` has `Primer3ProbeSettings.LowercaseMasking`. A `DnaSequence` template
   is upper case, so the option has no effect there.

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
| `pairOptions` | `PrimerPairOptions?` | `PrimerPairOptions.Default` | Pair options: `InsidePenalty` / `OutsidePenalty` (PRIMER_INSIDE/OUTSIDE_PENALTY, −1 / 0, §2.2 item 10), `ProductSizeRanges` (PRIMER_PRODUCT_SIZE_RANGE, default 100–300), `MaxTmDifference` (PRIMER_PAIR_MAX_DIFF_TM, default 5; `Primer3Defaults`: 100), `NumReturn` (PRIMER_NUM_RETURN, 5), `IncludedRegion` (SEQUENCE_INCLUDED_REGION), `ProductOptSize`/`ProductOptTm`/`ProductMinTm`/`ProductMaxTm`, `Weights` (PRIMER_PAIR_WT_*), `PickInternalOligo` + `InternalOligo` (PRIMER_PICK_INTERNAL_OLIGO, PRIMER_INTERNAL_*) | Primer3 `_pr_data_control` errors throw `ArgumentException` (weight without optimum, PRIMER_MAX_SIZE or PRIMER_INTERNAL_MAX_SIZE > min product size, NUM_RETURN < 1, target outside the included region) |

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
| Dinucleotide Repeats | N/A | N/A | 4 | Library default of an unsourced library screen (`FindLongestDinucleotideRepeat`; no Primer3 equivalent, `Primer3DefaultParameters` disables it — A3-8) |
| `GcClamp` (PRIMER_GC_CLAMP) | 0 | N/A | `MinLength` | Number of consecutive G/C required at the 3′ end (Primer3 default 0; > PRIMER_MIN_SIZE throws) |
| `MaxEndGc` (PRIMER_MAX_END_GC) | 0 | N/A | 5 | Max G/C among the five 3′-most bases; checked only when < 5 (default 5) |
| `MaxEndStability` (PRIMER_MAX_END_STABILITY) | 0 | N/A | 100 | Fails when `end_stability` = −`Calculate3PrimeStability` > limit (default 100; an ACGT pentamer reaches at most 6.86, GCGCG/CGCGC) |
| `Avoid3PrimeGC` | N/A | N/A | `false` | Deprecated library rule (not Primer3): despite the name it *requires* at least one `G`/`C` in the last two bases; kept for source compatibility — use `GcClamp` / `MaxEndGc` |
| `MisprimingLibrary` (PRIMER_MISPRIMING_LIBRARY) | N/A | N/A | null | name → sequence entries (`PrimerMisprimingLibrary`); null = no library check |
| `MaxLibraryMispriming` (PRIMER_MAX_LIBRARY_MISPRIMING) | N/A | N/A | 12 | Truncated to a C `short`; > 32767 rejected in alignment mode (`_pr_data_control`) |
| `PrimerPairOptions.MaxLibraryMispriming` (PRIMER_PAIR_MAX_LIBRARY_MISPRIMING) | N/A | N/A | 24 | Pair score = max over entries of ⌊left + right⌋ |
| `LibraryAmbiguityCodesConsensus` (PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS) | N/A | N/A | false (0) | Primer3 v2 default (`pr_set_default_global_args_2`); 1 makes IUPAC codes match their bases |
| `ThermodynamicTemplateAlignment` (PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT) | N/A | N/A | false (0) | true = ntthal END1 Tm instead of dpal score (§2.2 item 8) |
| `MaxTemplateMispriming` / `MaxTemplateMisprimingTh` (PRIMER_MAX_TEMPLATE_MISPRIMING[_TH]) | N/A | N/A | −100 | Negative = not checked; only the active mode's limit applies |
| `PrimerPairOptions.MaxTemplateMispriming` / `MaxTemplateMisprimingTh` (PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING[_TH]) | N/A | N/A | −100 | Alignment: ≥ 0 checked; thermodynamic: non-zero checked (Primer3) |
| `Primer3PenaltyWeights.TemplateMispriming[Th]`, `Primer3PairWeights.TemplateMispriming[Th]` (PRIMER_[PAIR_]WT_TEMPLATE_MISPRIMING[_TH]) | N/A | N/A | 0 | Must be ≥ 0 |
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
- `PrimerDesigner.DesignPrimers(string, …)` / `DesignPrimerPairs(string, …)`: the same on a case-preserving template string (PRIMER_LOWERCASE_MASKING, `PrimerPairOptions.LowercaseMasking`).
- `PrimerDesigner.CalculateProductMeltingTemperaturePrimer3(string, ...)`: Primer3 `long_seq_tm` product Tm.
- `PrimerDesigner.EvaluatePrimer(string, int, bool, PrimerParameters?)`: Scores a single primer candidate.
- `PrimerDesigner.CalculateLibraryMispriming(string, bool, PrimerMisprimingLibrary, bool)`: Primer3 library mispriming score of one primer ([PrimerDesigner.MisprimingLibrary.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.MisprimingLibrary.cs)).
- `PrimerDesigner.CalculateTemplateMispriming(DnaSequence, int, int, bool, bool, …)`: Primer3 template mispriming scores (`TemplateMisprimingScore` same-strand / other-strand / max) of one primer site ([PrimerDesigner.TemplateMispriming.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.TemplateMispriming.cs)).
- `PrimerDesigner.CalculateMeltingTemperaturePrimer3(string, ...)`: Primer3-default primer Tm used by design.
- `PrimerDesigner.CalculateFractionBoundPrimer3(string, double, ...)` / `CalculatePositionPenaltyPrimer3(...)`: Primer3 fraction bound at the annealing temperature and position penalty relative to the target ([PrimerDesigner.BoundAndPosition.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.BoundAndPosition.cs)); `PrimerParameters.AnnealingTemperature` / `MinBound` / `MaxBound` / `OptBound`.
- `PrimerDesigner.CalculateSequenceQualityPrimer3(IReadOnlyList<int>, int, int, bool, int)`: Primer3 `seq_quality` / `seq_end_quality` of one oligo ([PrimerDesigner.SequenceQuality.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.SequenceQuality.cs)); `PrimerPairOptions.SequenceQuality`, `PrimerParameters.MinQuality` / `MinEndQuality` / `QualityRangeMin` / `QualityRangeMax`.
- `PrimerDesigner.MaskTemplatePrimer3(string, PrimerMaskingKmerLists, double, int, int)` / `CalculateMaskFailureRatePrimer3(string, PrimerMaskingKmerLists)` / `PrimerMaskingKmerLists` (+ `FromGenomeTester4Files`): Primer3's k-mer masker ([PrimerDesigner.Masking.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.Masking.cs)); `PrimerPairOptions.MaskTemplate` / `MaskKmerLists` / `MaskFailureRate` / `MaskFivePrimeDirection` / `MaskThreePrimeDirection`.
- `PrimerDesigner.CalculatePrimer3Penalty(...)`: Primer3 per-primer penalty used for ranking.
- `PrimerDesigner.CalculatePrimerScore(...)` (private): informational heuristic score.

### 5.2 Current Behavior

Forward primers are taken directly from the template; reverse primers are reverse-complemented before evaluation, and their `Position` is the leftmost template coordinate of the binding site. Per-primer hard constraints: length, GC%, Primer3 3′-end checks (PRIMER_GC_CLAMP / PRIMER_MAX_END_GC / PRIMER_MAX_END_STABILITY, `calc_and_check_oligo_features`), Primer3-default Tm window, homopolymer, dinucleotide repeat, secondary structure (default `PrimerStructureScreen.Primer3Thermodynamic`: Primer3 ntthal self-any / self-end / hairpin Tm ≤ `MaxStructureTm` = 47 °C, evaluated lazily in the pair loop like Primer3's `characterize_pair`; `Heuristic`: `HasHairpinPotential`), the deprecated library rule `Avoid3PrimeGC` (≥ 1 G/C in the last two bases), and no non-ACGT base (Primer3 PRIMER_MAX_NS_ACCEPTED = 0). Pair selection is the Primer3 pair search described in §2.2; the internal oligo reuses the `ProbeDesigner` Primer3 internal-oligo list (`EnumeratePrimer3InternalOligos`, shared with `DesignProbesPrimer3`).

### 5.3 Conformance to Theory / Spec

**Implemented (verified against primer3-py 2.3.1):**

- Primer3 default Tm: bit-identical to `primer3.calc_tm(..., max_nn_length=36)` — libprimer3's MAX_NN_TM_LENGTH used when picking oligos (max |Δ| = 0 over 3 000 random 2–45-mers, incl. self-complementary and > 36 nt). primer3-py's `calc_tm` itself defaults to `max_nn_length=60`, so for 37–60 nt the five-argument `CalculateMeltingTemperaturePrimer3` differs from a bare `calc_tm(seq)` (39-mer `CTGTATAGTCCCACCTGGTGATCCTATGCTTGTGAGTAC`: 73.0492 vs 70.4249 °C); the `maxNearestNeighborLength` overload with 60 reproduces the `calc_tm` default (audit round 5, A5-4).
- Pair search with product-size ranges, PRIMER_NUM_RETURN, pair objective and internal oligo (audit round 1, 2026-10-01): with `Primer3DefaultParameters` + `PrimerPairOptions.Primer3Defaults`, `DesignPrimerPairs` reproduced `primer3.design_primers` (only SEQUENCE_TEMPLATE/SEQUENCE_TARGET set) for 1000/1000 random templates (150–700 bp, 4990 pairs; left/right position+length, PRIMER_PAIR_k_PENALTY, _PRODUCT_TM, _COMPL_ANY_TH, _COMPL_END_TH, _PRODUCT_SIZE, ranks 0–4, |Δ| ≤ 1e-9); with random non-default PRIMER_PAIR_WT_* / PRODUCT_OPT_* / PRODUCT_MIN/MAX_TM / PAIR_MAX_DIFF_TM / multi-range PRIMER_PRODUCT_SIZE_RANGE / SEQUENCE_INCLUDED_REGION settings and with PRIMER_PICK_INTERNAL_OLIGO = 1 (+ PRIMER_PAIR_WT_IO_PENALTY; PRIMER_INTERNAL_k position, penalty, Tm, self-any/self-end/hairpin Tm) the counts are in the PRIMER-DESIGN-001 F-entry of `docs/Validation/review-2026-09/B07.md`.
- Primer3 alignment mode (audit round 2, A1, 2026-10-02): with `PrimerStructureScreen.Primer3Alignment` (PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT = 0: dpal self_any/self_end ≤ PRIMER_MAX_SELF_ANY/_END 8/3, pair compl_any/compl_end ≤ PRIMER_PAIR_MAX_COMPL_ANY/_END 8/3, internal oligo ≤ 12/12, PRIMER_WT_SELF_ANY/_END, PRIMER_PAIR_WT_COMPL_ANY/_END) `DesignPrimerPairs` reproduced `primer3.design_primers` for 2000/2000 random templates (9135 pairs: positions, PRIMER_PAIR_k_PENALTY, _COMPL_ANY, _COMPL_END, PRIMER_LEFT/RIGHT_k_SELF_ANY/_SELF_END/_PENALTY, product Tm, PRIMER_INTERNAL_k_* ; 1000 templates with random non-default limits/weights/internal oligos); per-primer thermodynamic weights PRIMER_WT_SELF_ANY_TH/_SELF_END_TH/_HAIRPIN_TH (`PrimerParameters.PenaltyWeights`) 60/60 templates (280 pairs) identical.
- Reaction conditions and GC optimum (audit round 3, A3-1, 2026-10-02): random PRIMER_SALT_MONOVALENT 10–200 mM, PRIMER_SALT_DIVALENT {0, 0–6} mM, PRIMER_DNTP_CONC {0, 0–2} mM, PRIMER_DNA_CONC 10–500 nM, PRIMER_OPT_GC_PERCENT 25–75 with PRIMER_WT_GC_PERCENT_GT/_LT ∈ {0, 0.25–1}, the PRIMER_INTERNAL_* counterparts (+ PRIMER_INTERNAL_OPT_GC_PERCENT / _WT_GC_PERCENT_GT/_LT) with PRIMER_PICK_INTERNAL_OLIGO ∈ {0, 1}, both PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT modes and an optional product-Tm objective: 2500/2500 templates (12 440 pairs) identical on positions, PRIMER_PAIR/LEFT/RIGHT_k_PENALTY, _TM, PRIMER_PAIR_k_PRODUCT_TM, the self-any/self-end/hairpin (or dpal) and pair compl values and PRIMER_INTERNAL_k position/penalty/Tm (relative |Δ| ≤ 1e-9); 1/2500 before (fixed 50 mM / 1.5 mM / 0.6 mM / 50 nM, GC optimum 50).
- Primer3 per-primer penalty and pair search (before the audit, ±200 bp flanks): `DesignPrimers` returned exactly Primer3's `PRIMER_LEFT_0`/`PRIMER_RIGHT_0` in 553/553 random templates (3 seeds × 300) where Primer3's best pair also passes this library's extra screens (settings mirroring `DefaultParameters`, thermodynamic structure limits disabled).

**Deviations from Primer3 defaults (documented; `Primer3DefaultParameters` / `PrimerPairOptions.Primer3Defaults` restore them):** length 18–25 (Primer3 18–27), GC 40–60 % (20–80 %), poly-X 4 (5), pair ΔTm ≤ 5 °C (100), dinucleotide-repeat limit (unsourced library screen, no Primer3 equivalent; A3-8). The product-size range is Primer3's (default 100–300 bp; before audit round 1 the ±200 bp flanks bounded the product instead). Structure limits are Primer3's default thermodynamic ones (PRIMER-STRUCT-001): with them `DesignPrimers` returned primer3-py's pair (same settings) in 574/600 random templates; the 26 differences all trace to ntthal engine values (PRIMER-DIMER-001 / PRIMER-HAIRPIN-001). After PRIMER-DIMER-001 (dimer engine bit-exact) a re-run on 1800 random templates (seeds 1–9 × 200) gave 1733 identical; all 67 differences trace to `calc_hairpin` values only (PRIMER-HAIRPIN-001). After PRIMER-HAIRPIN-001 (hairpin engine bit-exact) the same 1800 templates are 1800/1800 identical (primers and pair penalty).

- 3′-end checks and 3′ distance (audit round 3, A3-6 + A3-7, 2026-10-02): random 150–600-bp templates (GC bias 35–60 %), PRIMER_GC_CLAMP ∈ {0–3}, PRIMER_MAX_END_GC ∈ {0–5}, PRIMER_MAX_END_STABILITY ∈ {4–9, 100}, PRIMER_MIN_THREE_PRIME_DISTANCE or PRIMER_MIN_LEFT/RIGHT_THREE_PRIME_DISTANCE ∈ {−1, 0, 1, 2, 3, 5, 10, 20}, PRIMER_NUM_RETURN 5–10, PRIMER_WT_END_STABILITY ∈ {0, 0.5}, both alignment modes, PRIMER_PICK_INTERNAL_OLIGO ∈ {0, 1}: 2200/2200 templates (9590 pairs) identical on left/right start + length, PRIMER_PAIR/LEFT/RIGHT_k_PENALTY, _END_STABILITY and PRIMER_INTERNAL_k position (|Δ| ≤ 1e-9); without the new options 36/2000.

- Mispriming library (audit round 3, A3-3 part 1, 2026-10-02): random 150–500-bp templates with 1–8-entry random libraries (55 % template fragments, forward or reverse-complemented, with point mutations and IUPAC codes; random IUPAC-containing sequences; 1–2-nt entries; N-rich entries; `*w` weights 0–10), PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS ∈ {0, 1}, PRIMER_MAX_LIBRARY_MISPRIMING ∈ {6, 8, 10, 12, 12.9, 15, 20}, PRIMER_PAIR_MAX_LIBRARY_MISPRIMING ∈ {12–30}, PRIMER_WT_LIBRARY_MISPRIMING ∈ {0, 0.1, 0.5, 1}, PRIMER_PAIR_WT_LIBRARY_MISPRIMING ∈ {0, 0.2, 1}, both alignment modes, PRIMER_PICK_INTERNAL_OLIGO ∈ {0, 1}: 1240/1240 templates (4693 pairs; primers rejected by the library in 1005 templates, pairs in 202) identical to primer3-py 2.3.1 `design_primers(misprime_lib=…)` on left/right start + length, PRIMER_PAIR/LEFT/RIGHT_k_PENALTY (|Δ| ≤ 1e-9), PRIMER_LEFT/RIGHT/PAIR_k_LIBRARY_MISPRIMING score and entry name, and PRIMER_INTERNAL_k position.

- Internal-oligo mishybridization library (audit round 3, A3-3 part 2, 2026-10-02): see §2.2 item 7 and the PROBE-DESIGN-001 / PRIMER-DESIGN-001 F42 entry of `docs/Validation/review-2026-09/B07.md` — primer3-py 2.3.1 `design_primers(mishyb_lib=…)` with PRIMER_PICK_INTERNAL_OLIGO = 1: 386/386 random 150–450-bp templates (1325 pairs; random 1–8-entry libraries, PRIMER_INTERNAL_MAX_LIBRARY_MISHYB ∈ {6–20}, PRIMER_INTERNAL_WT_LIBRARY_MISHYB ∈ {0, 0.1, 0.5, 1}, PRIMER_PAIR_WT_IO_PENALTY ∈ {0, 1}, both alignment modes and consensus settings, 137 with a primer mispriming library too; internal oligos rejected by the library in 315) identical on left/right start + length, PRIMER_PAIR_k_PENALTY, PRIMER_INTERNAL_k position / penalty (|Δ| ≤ 1e-9) and PRIMER_INTERNAL_k_LIBRARY_MISHYB score + entry name.

- Template mispriming (audit round 3, A3-4, 2026-10-02): random 150–500-bp templates (thermodynamic template mode 150–260 bp), 70 % made repetitive (1–6 copied 12–40-nt fragments, forward or reverse-complemented, 0–3 mutations; 30 % with a tandem repeat), PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT ∈ {0, 1}, PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT ∈ {0, 1}, PRIMER_MAX_TEMPLATE_MISPRIMING ∈ {−100, 6–15}, _TH ∈ {−100, 15–50}, PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING ∈ {−100, 14–24}, _TH ∈ {−100, 0, 30–80}, PRIMER_WT_TEMPLATE_MISPRIMING[_TH] ∈ {0, 0.1, 0.5, 1}, PRIMER_PAIR_WT_TEMPLATE_MISPRIMING[_TH] ∈ {0, 0.2, 1}, 14 % with a mispriming library, 15 % PRIMER_PICK_INTERNAL_OLIGO: 352/352 templates (1590 pairs; 183 thermodynamic; primers rejected by the template limit in 124 templates, pairs in 67) identical to primer3-py 2.3.1 `design_primers` on left/right start + length, PRIMER_PAIR/LEFT/RIGHT_k_PENALTY (|Δ| ≤ 1e-9), PRIMER_LEFT/RIGHT_k_TEMPLATE_MISPRIMING (alignment mode; primer3-py omits the per-primer _TH keys in thermodynamic mode), PRIMER_PAIR_k_TEMPLATE_MISPRIMING[_TH] and PRIMER_INTERNAL_k position; per-strand END1 values against `primer3.calc_end_stability`.

- Fraction bound and position penalty (audit round 3, A3-5 part 1, 2026-10-02): random 150–450-bp templates with a random 5–50-nt target, PRIMER_ANNEALING_TEMP ∈ {unset, −10, 0, 40–70}, random PRIMER_MIN/MAX/OPT_BOUND and PRIMER_INTERNAL_MIN/MAX/OPT_BOUND, PRIMER_[INTERNAL_]WT_BOUND_GT/_LT ∈ {0, 1e-7–0.5}, random reaction conditions, PRIMER_INSIDE_PENALTY ∈ {−1, 0, 0.1, 0.5, 1, 2} / PRIMER_OUTSIDE_PENALTY ∈ {0, 0.05–1} (60 % non-default), PRIMER_WT_POS_PENALTY ∈ {0, 0.5, 1, 2}, optional product-size ranges, PRIMER_PICK_INTERNAL_OLIGO (+ PRIMER_PAIR_WT_IO_PENALTY), both alignment modes, plus 85 pick_hyb_probe_only probe lists: 400/400 cases identical to primer3-py 2.3.1 `design_primers` (1304 pairs + 432 probes; left/right start + length, PRIMER_PAIR/LEFT/RIGHT/INTERNAL_k_PENALTY |Δ| ≤ 1e-9, PRIMER_LEFT/RIGHT/INTERNAL_k_BOUND, PRIMER_LEFT/RIGHT_k_POSITION_PENALTY, PRIMER_INTERNAL_k position); the 34 settings where primer3-py aborts on a negative pair penalty throw `InvalidOperationException`. Details: F44 in `docs/Validation/review-2026-09/B07.md`.
- Sequence quality and PRIMER_PAIR_WT_IO_PENALTY consistency (audit round 3, A3-5 part 2a + A3-23, 2026-10-02): 400 random cases (150–400-bp pair templates with a random 5–50-nt target, 20 % pick_hyb_probe_only lists of 40–160 nt; SEQUENCE_QUALITY uniform or high with low dips in 85 % of cases, occasionally one base short or one value out of range; PRIMER_QUALITY_RANGE_MIN/MAX ∈ {0–100, 0–60, 10–40, 5–93}; random PRIMER_[INTERNAL_]MIN_QUALITY, PRIMER_MIN_END_QUALITY, PRIMER_[INTERNAL_]WT_SEQ_QUAL / _WT_END_QUAL; PRIMER_PICK_INTERNAL_OLIGO 35 %, PRIMER_PAIR_WT_IO_PENALTY ∈ {0, 0.5, 1}; both alignment modes): 394/394 comparable cases identical to primer3-py 2.3.1 `design_primers` (1104 pairs + probes on start + length, PRIMER_PAIR/LEFT/RIGHT/INTERNAL_k_PENALTY |Δ| ≤ 1e-9 and PRIMER_LEFT/RIGHT/INTERNAL_k_MIN_SEQ_QUALITY; 110 `_pr_data_control` error cases raise `ArgumentException` with Primer3's first message); 6 probe-only cases that set primer-side (p_args) quality settings have no `DesignProbesPrimer3` counterpart. Details: F45 in `docs/Validation/review-2026-09/B07.md`.
- Template masking (audit round 3, A3-5 part 2b, 2026-10-02): masker.c compiled unchanged (primer3-py 2.3.1 sources) with a driver setting libprimer3's masker parameters: 300/300 random templates (1–12000 nt, i.e. across the 5000-character ring buffer; random lists incl. zero counts, k-mers listed in both orientations; PRIMER_MASK_FAILURE_RATE ∈ {0–0.5}, 5P ∈ {0–20}, 3P ∈ {0–4}) give identical forward / reverse masked copies (535 189 masked bases) and identical failure rates for 6000 primers (2025 non-zero, |Δ| ≤ 1e-15), with the lists given as dictionaries and read back from GenomeTester4 files. primer3-py 2.3.1 `design_primers` with PRIMER_MASK_TEMPLATE on 400 random 150–450-bp templates (per-case GenomeTester4 lists built from template k-mers, 10 % with PRIMER_MASK_TEMPLATE 0, random PRIMER_MASK_FAILURE_RATE / 5P / 3P / PRIMER_WT_MASK_FAILURE_RATE ∈ {0, 0.5–20}, included regions, internal oligos, product ranges, both alignment modes; primers rejected for a masked 3′ end in 227 templates): 400/400 identical — 1370 pairs on left/right start + length, PRIMER_PAIR/LEFT/RIGHT/INTERNAL_k_PENALTY (|Δ| ≤ 1e-9). Details: F45 (part 2b) in `docs/Validation/review-2026-09/B07.md`.

- Undefined GC optimum and lower-case masking (audit round 3, A3-25 + A3-26, 2026-10-02): 600 random cases (150–450-bp templates, 85 % with random lower-case runs and isolated lower-case bases; PRIMER_LOWERCASE_MASKING ∈ {unset, 0, 1}; 20 % of the pair cases with PRIMER_MASK_TEMPLATE and per-case k-mer lists; 25 % pick_hyb_probe_only; PRIMER_PICK_INTERNAL_OLIGO 35 %; random PRIMER_[INTERNAL_]OPT_GC_PERCENT set or unset with PRIMER_[INTERNAL_]WT_GC_PERCENT_GT/_LT ∈ {unset, 0, 0.25–1}; both alignment modes): 600/600 identical to primer3-py 2.3.1 `design_primers` (1316 pairs + probes on start + length and PRIMER_PAIR/LEFT/RIGHT/INTERNAL_k_PENALTY |Δ| ≤ 1e-9; 266 `_pr_data_control` GC-optimum errors reproduced with Primer3's message). Ignoring the template case: 471/600. Details: F46 in `docs/Validation/review-2026-09/B07.md`.

**Outside this unit's scope (Primer3 tags not modelled; not in the PRIMER-DESIGN-001 test spec):**
SEQUENCE_OVERLAP_JUNCTION_LIST / SEQUENCE_INTERNAL_OVERLAP_JUNCTION_LIST with PRIMER_MIN_5/3_PRIME_OVERLAP_OF_JUNCTION
and PRIMER_INTERNAL_MIN_5/3_PRIME_OVERLAP_OF_JUNCTION — Primer3's junction-spanning mode (e.g. exon–exon junctions of a
cDNA template), a separate feature that the API has no input for. Without a junction list Primer3 itself never reads
PRIMER_INTERNAL_MIN_THREE_PRIME_DISTANCE, so the library equals Primer3 for every input it accepts (§2.2 item 5;
audit round 5, A5-3).

**BLOCKED — genome-wide specificity:** needs a BLAST search of the primers against a genome / nucleotide database
(NCBI Primer-BLAST [8] is a web service that couples Primer3 with BLAST over NCBI databases); Primer3 itself has no
genome search and the API receives only the template. The in-scope Primer3 specificity checks are the mispriming
library (§2.2 item 6) and template mispriming (item 8); approximate probe off-target search is PROBE-VALID-001.

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
| SEQUENCE_QUALITY length ≠ template, PRIMER_[INTERNAL_]MIN_QUALITY ≠ 0 without quality or outside the quality range, a quality value outside [PRIMER_QUALITY_RANGE_MIN, MAX], PRIMER_[INTERNAL_]WT_SEQ_QUAL ≠ 0 without quality | `ArgumentException` | Primer3 `_pr_data_control` ("Error in sequence quality data", "Sequence quality data missing", "PRIMER_[INTERNAL_]MIN_QUALITY < / > PRIMER_QUALITY_RANGE_MIN / MAX", "Sequence quality score out of range", "Sequence quality is part of objective function but sequence quality is not defined") |
| PRIMER_PAIR_WT_IO_PENALTY ≠ 0 without `PickInternalOligo` | `ArgumentException` | Primer3 `_pr_data_control` "Internal oligo quality is part of objective function while internal oligo choice is not required" |
| `PickInternalOligo` and no acceptable internal oligo anywhere in the included region | No pairs (`DesignPrimers`: invalid result naming the internal oligo) | `make_internal_oligo_list` fails, Primer3 returns before the pair search |
| PRIMER_[INTERNAL_]WT_GC_PERCENT_GT/_LT ≠ 0 without PRIMER_[INTERNAL_]OPT_GC_PERCENT (also the internal-oligo values without `PickInternalOligo`) | `ArgumentException` | Primer3 `_pr_data_control` "Primer GC content …" / "Hyb probe GC content is part of objective function while optimum gc_content is not defined" (optimum undefined by default) |
| Lower-case template base at an oligo's 3′ end with `LowercaseMasking` (or `MaskTemplate`) and the `string` template overload | Oligo rejected | `is_lowercase_masked` |
| PRIMER_WT_END_QUAL / PRIMER_INTERNAL_WT_END_QUAL ≠ 0 | No effect | Primer3 2.3.1 `p_obj_fn` never reads `weights.end_quality` |
| `MaskTemplate` without `MaskKmerLists`; negative mask directions, PRIMER_MASK_3P_DIRECTION > 4984, non-finite failure rate | `ArgumentException` / `ArgumentOutOfRangeException` | primer3-py "masking template chosen, but path to PRIMER_MASK_KMERLIST_PATH not specified"; masker.c undefined for the others |
| PRIMER_WT_MASK_FAILURE_RATE without `MaskTemplate` | No effect (failure rate 0) | `calc_and_check_oligo_features` computes `failure_rate` only with `mask_template` |
| PRIMER_INTERNAL_WT_LIBRARY_MISHYB ≠ 0 without a mishyb library (with `PickInternalOligo`) | `ArgumentException` | Primer3 `_pr_data_control` "Internal oligo mispriming score is part of objective function while mishyb library is not defined" |
| PRIMER_MAX_TEMPLATE_MISPRIMING / PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING > 32767 (alignment mode), negative template weight | `ArgumentOutOfRangeException` / `ArgumentException` | Primer3 `_pr_data_control`; a negative weight would trip Primer3's `p_obj_fn` assertion |
| Pair template limit/weight with PRIMER_WT_LIBRARY_MISPRIMING ≠ 0 but no per-primer template weight | `InvalidOperationException` | Primer3 never scores the primers' template mispriming then and aborts (`PR_ASSERT` in `characterize_pair`) |
| PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT = 1 with a template > 10000 nt | `ArgumentException` | thal THAL_MAX_SEQ |
| PRIMER_OPT_BOUND / PRIMER_INTERNAL_OPT_BOUND outside [MIN, MAX]; PRIMER_ANNEALING_TEMP > 100 | `ArgumentOutOfRangeException` | Primer3 `_pr_data_control` ("Optimum primer fraction binding lower than minimum or higher than maximum", "Annealing temperature higher than 100 C") |
| Bound weights without PRIMER_ANNEALING_TEMP | primer terms inert; internal-oligo terms use bound = −999999.9999 | `p_obj_fn` gates only the primer branch |
| Non-default inside/outside penalty giving a negative pair penalty (e.g. only PRIMER_OUTSIDE_PENALTY changed, inside stays −1) | `InvalidOperationException` | Primer3 aborts (`obj_fn` `PR_ASSERT(sum >= 0.0)`) |
| PRIMER_WT_LIBRARY_MISPRIMING / PRIMER_PAIR_WT_LIBRARY_MISPRIMING ≠ 0 without a library | `ArgumentException` | Primer3 `_pr_data_control` "Mispriming score is part of objective function, but mispriming library is not defined" |
| Library entry with an empty sequence or an illegal `*weight` (missing, < 0, > 100) | `ArgumentException` | `add_seq_to_seq_lib` / `parse_seq_name` (primer3-py raises OSError) |
| Library entry with a non-IUPAC character | Character becomes N, `PrimerMisprimingLibrary.Warnings` | `upcase_and_check_char` (primer3-py 2.3.1 aborts here: it passes a NULL `errfrag` to the warning) |

### 6.2 Limitations

Template mispriming, the primer mispriming library and the internal-oligo mishybridization library are Primer3's
(§2.2 items 6–8); `EvaluatePrimer` has no template, so it never applies the template terms (use
`CalculateTemplateMispriming` + `CalculatePrimer3Penalty`). The secondary-structure screen is Primer3's thermodynamic one by default (the sequence-only screen is available as `PrimerStructureScreen.Heuristic`).

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
8. Ye J, Coulouris G, Zaretskaya I, Cutcutache I, Rozen S, Madden TL (2012). "Primer-BLAST: a tool to design target-specific primers for polymerase chain reaction", BMC Bioinformatics 13:134.
