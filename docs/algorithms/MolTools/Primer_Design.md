# Primer Pair Design

| Field | Value |
|-------|-------|
| Algorithm Group | MolTools |
| Test Unit ID | PRIMER-DESIGN-001 |
| Related Projects | N/A |
| Implementation Status | Implemented (Primer3 pair selection; heuristic structure screens) |
| Last Reviewed | 2026-09-28 |

## 1. Overview

Primer pair design selects forward and reverse oligonucleotides that can amplify a target DNA region by PCR. In this repository, primer design enumerates every candidate in the flanking regions, filters each by per-primer constraints, and then — exactly as Primer3 (`libprimer3.cc` `choose_pair_or_triple`) — returns the pair with the lowest Primer3 pair penalty among all pairs meeting the pair constraints (Tm agreement, primer-dimer avoidance). The Tm used everywhere in design is Primer3's default primer Tm, the scale on which the Primer3 Tm window 57–63 °C is defined.

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
   Defaults 50 nM oligo, 50 mM monovalent, 1.5 mM Mg²⁺, 0.6 mM dNTP.
2. **Per-primer penalty** (Primer3 `p_obj_fn`, default weights):
   $penalty = |T_m - OptimalTm| + |length - OptimalLength|$ (`CalculatePrimer3Penalty`).
3. **Pair selection:** minimise $penalty_f + penalty_r$ (PRIMER_PAIR_WT_PR_PENALTY = 1, other pair
   weights 0) over all pairs with $|T_{m,f} - T_{m,r}| \le 5$ °C and no primer-dimer; ties within
   $10^{-6}$ broken as `compare_primer_pair` (left primer further 3′, right primer 5′ end further
   left, shorter left, shorter right).

`PrimerCandidate.Score` (100 − 2|len − opt| − 2|Tm − opt| − 0.5|GC − 50| − 5·homopolymer + 5 GC-clamp
bonus) is reported for information only and does not drive selection.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `DesignPrimers(...)` returns `IsValid = false` when either side has no valid candidates | The source returns an invalid `PrimerPairResult` when either best candidate is missing |
| INV-02 | Pair validity requires both `|Tm_f - Tm_r| <= 5` (unrounded Tm) and `!HasPrimerDimer(...)`; if any such pair exists among the valid candidates, `IsValid = true` | Exhaustive pair search |
| INV-04 | The returned valid pair minimises `Forward.Penalty + Reverse.Penalty` over all compatible pairs | Primer3 `choose_pair_or_triple` |
| INV-03 | `ProductSize = reverse.Position + reverse.Sequence.Length - forward.Position` | The source computes product size directly from the chosen candidates |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `template` | `DnaSequence` | required | Template DNA sequence |
| `targetStart` | `int` | required | Start of the target region | Must satisfy `targetStart >= 0` |
| `targetEnd` | `int` | required | Exclusive end of the target region (target = `[targetStart, targetEnd)`, Primer3 SEQUENCE_TARGET) | Must satisfy `targetEnd < template.Length` and `targetStart < targetEnd` |
| `parameters` | `PrimerParameters?` | `PrimerDesigner.DefaultParameters` | Primer design thresholds | Defaults are `18-25` bp length, `40-60%` GC, `57-63°C` Tm, `OptimalLength = 20`, `OptimalTm = 60`, `MaxHomopolymer = 4`, `MaxDinucleotideRepeats = 4`, `Avoid3PrimeGC = false`, and `Check3PrimeStability = true` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Forward` | `PrimerCandidate?` | Selected forward primer or `null` |
| `Reverse` | `PrimerCandidate?` | Selected reverse primer or `null` |
| `IsValid` | `bool` | Pair validity flag |
| `Message` | `string` | Result explanation |
| `ProductSize` | `int` | Predicted amplicon size |

### 3.3 Preconditions and Validation

`DesignPrimers(...)` throws `ArgumentException` when the requested target region is invalid. Forward candidates are searched up to 200 bp upstream of `targetStart`; reverse candidates are searched up to 200 bp downstream of `targetEnd`. Reverse-primer candidates are reverse-complemented before evaluation so that they are scored in primer orientation.

## 4. Algorithm

### 4.1 High-Level Steps

1. Define a forward search region up to 200 bp upstream of the target start.
2. Define a reverse search region up to 200 bp downstream of the target end.
3. Enumerate all candidate primers within the configured length range.
4. Evaluate each candidate for GC content, Tm, homopolymers, dinucleotide repeats, hairpin potential, and 3' stability.
5. Sort each side by Primer3 penalty; scan reverse × forward candidates with Primer3's pruning (stop a row once `penalty_f + penalty_r` exceeds the best pair found).
6. Keep the lowest-penalty pair with Tm difference ≤ `5°C` and no primer-dimer; if none exists, return the individually best primers with `IsValid = false` and a message naming the violated constraint.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Parameter ranges documented in the original file and current source:

| Parameter | Min | Optimal | Max | Notes |
|-----------|-----|---------|-----|-------|
| Length (bp) | 18 | 20 | 25 | Matches current source defaults |
| GC Content (%) | 40 | 50 | 60 | Matches current source defaults |
| Melting Temp (°C) | 55 | 60 | 65 | Original document summary; current source defaults narrow this to `57-63` |
| Homopolymer Run | N/A | N/A | 4 | Current source default |
| Dinucleotide Repeats | N/A | N/A | 4 | Current source default |
| `Avoid3PrimeGC` | N/A | N/A | `false` | When enabled, the current check requires at least one `G`/`C` in the last two bases |
| `Check3PrimeStability` | N/A | N/A | `true` | Gates whether `EvaluatePrimer(...)` records the `ΔG < -9` issue |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `DesignPrimers` | `O(n²)` | `O(k)` | Candidate enumeration over positions and lengths |
| `EvaluatePrimer` | `O(n)` | `O(1)` | Per-primer scan and helper calculations |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [PrimerDesigner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs)

- `PrimerDesigner.DesignPrimers(DnaSequence, int, int, PrimerParameters?)`: Designs and validates a primer pair around a target region.
- `PrimerDesigner.EvaluatePrimer(string, int, bool, PrimerParameters?)`: Scores a single primer candidate.
- `PrimerDesigner.CalculateMeltingTemperaturePrimer3(string, ...)`: Primer3-default primer Tm used by design.
- `PrimerDesigner.CalculatePrimer3Penalty(...)`: Primer3 per-primer penalty used for ranking.
- `PrimerDesigner.CalculatePrimerScore(...)` (private): informational heuristic score.

### 5.2 Current Behavior

Forward primers are taken directly from the template; reverse primers are reverse-complemented before evaluation, and their `Position` is the leftmost template coordinate of the binding site. Per-primer hard constraints: length, GC%, Primer3-default Tm window, homopolymer, dinucleotide repeat, heuristic hairpin (`HasHairpinPotential`), 3′ ΔG (`< −9` kcal/mol flagged; note that the SantaLucia 5-mer ΔG never goes below −6.86, so this gate never fires — consistent with Primer3's default PRIMER_MAX_END_STABILITY = 100), optional GC clamp, and no non-ACGT base (Primer3 PRIMER_MAX_NS_ACCEPTED = 0). Pair selection is the exhaustive Primer3 pair search described in §2.2.

### 5.3 Conformance to Theory / Spec

**Implemented (verified against primer3-py 2.3.1):**

- Primer3 default Tm: bit-identical to `primer3.calc_tm` (max |Δ| = 0 over 3 000 random 2–45-mers, incl. self-complementary and > 36 nt).
- Primer3 per-primer penalty and pair search: `DesignPrimers` returned exactly Primer3's `PRIMER_LEFT_0`/`PRIMER_RIGHT_0` in 553/553 random templates (3 seeds × 300) where Primer3's best pair also passes this library's extra screens (settings mirroring `DefaultParameters`, thermodynamic structure limits disabled).

**Deviations from Primer3 defaults (documented):** length 18–25 (Primer3 18–27), GC 40–60 % (20–80 %), poly-X 4 (5), pair ΔTm ≤ 5 °C (100), dinucleotide-repeat limit (no Primer3 equivalent), heuristic hairpin / primer-dimer screens instead of ntthal PRIMER_MAX_HAIRPIN_TH / PRIMER_PAIR_MAX_COMPL_*_TH (owned by PRIMER-STRUCT-001), no PRIMER_PRODUCT_SIZE_RANGE (the ±200 bp flanks bound the product). In ~30 % of random templates Primer3's best pair is rejected by the heuristic hairpin screen, so the returned pair differs from Primer3's.

**Not implemented:** mispriming libraries, internal oligos, multiple returned pairs, genome-wide specificity.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Invalid target region | Throws `ArgumentException` | Explicit source guard |
| No valid forward or reverse candidates | Returns an invalid `PrimerPairResult` with null candidates | Explicit fallback in source |
| No pair within 5 °C | Returns the individually best primers with `IsValid = false` | Pair compatibility requires `<= 5°C` |
| Primer-dimer detected for every pair | Returns `IsValid = false` | Pair compatibility requires no dimer signal |
| Non-ACGT base in a candidate | Candidate invalid (Tm 0, issue "Tm not computable") | Primer3 PRIMER_MAX_NS_ACCEPTED = 0 |

### 6.2 Limitations

Structure screens (hairpin, primer-dimer) are the heuristic ones of PRIMER-STRUCT-001 rather than Primer3's ntthal Tm limits, and there is no product-size range or mispriming check.

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
6. Primer3 source (primer3-org/primer3, `src/oligotm.c`: `oligotm`, `seqtm`, `long_seq_tm`, `divalent_to_monovalent`; `src/libprimer3.cc`: `choose_pair_or_triple`, `primer_rec_comp`, `compare_primer_pair`, `p_obj_fn`).
7. von Ahsen N, Wittwer CT, Schütz E (2001). Clin Chem 47:1956-61 (divalent→monovalent equivalence).
