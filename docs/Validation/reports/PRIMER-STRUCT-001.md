# Validation Report: PRIMER-STRUCT-001 — Primer Structure Analysis

- **Validated:** 2026-06-12   **Area:** Molecular Tools (Seqeron.Genomics.MolTools)
- **Canonical method(s):** `PrimerDesigner.HasHairpinPotential(seq, minStemLength, minLoopLength)`,
  `PrimerDesigner.HasPrimerDimer(primer1, primer2, minComplementarity)`,
  `PrimerDesigner.Calculate3PrimeStability(seq)`,
  `PrimerDesigner.FindLongestHomopolymer(seq)`, `PrimerDesigner.FindLongestDinucleotideRepeat(seq)`
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs`
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/PrimerDesigner_PrimerStructure_Tests.cs`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN

---

## Stage A — Description

### Sources opened & what they confirm

- **Wikipedia — Stem-loop** (https://en.wikipedia.org/wiki/Stem-loop): A stem-loop/hairpin
  forms when two complementary (self-complementary) regions of a single strand base-pair into
  a duplex "stem" with unpaired bases as a "loop". Confirms the minimum-loop constraint:
  *"loops that are fewer than three bases long are sterically impossible and thus do not form"*;
  optimal loop 4–8 nt. Stability rises with G·C content (3 H-bonds vs 2 for A·T). This validates
  hairpin = internal inverted repeat (stem) + loop, with min loop ≈ 3 nt.
- **Wikipedia — Primer dimer** (https://en.wikipedia.org/wiki/Primer_dimer): a PD is "two primer
  molecules that have attached (hybridized) to each other because of strings of complementary
  bases in the primers." Formation: (1) primers anneal **at their respective 3' ends**, with
  high GC and overlap length stabilising; (2) polymerase extends; (3) amplification. Distinguishes
  self-complementarity (a primer binding itself = self-dimer) from complementarity between
  different primers (cross-dimer). Validates the 3'-end emphasis and Watson-Crick complementarity.
- **Primer3 manual / Primer3Web 4.1.0 help** (primer3.org, primer3.ut.ee): `PRIMER_MAX_SELF_ANY`
  limits self-complementarity over all alignments; `PRIMER_MAX_SELF_END` limits 3'-terminal
  self-complementarity. `PRIMER_MAX_END_STABILITY` (default 100) = max ΔG (kcal/mol) for duplex
  disruption of the **five 3' bases** computed from nearest-neighbour parameters; bigger = more
  stable 3' end. Validates the 5-mer 3'-stability ΔG concept.
- **SantaLucia (1998)** PNAS 95:1460-65, Table 1 unified NN ΔG°37 parameters. Published stability
  trend: GC > CG > GG > GA≈GT≈CA > CT > AA > AT > TA, which matches the implementation's
  dictionary exactly. Initiation: +0.98 (terminal G·C), +1.03 (terminal A·T).

### Detection logic validated

- **Self-dimer / cross-dimer**: a primer's region aligns antiparallel to the reverse complement
  of the other primer; Watson-Crick pairing A-T / G-C; deleterious pairing is at the **3' ends**.
  Score = number of complementary base pairs in the 3'-terminal window; "problematic" when count
  ≥ a configurable threshold. (Self-dimer is the special case primer1 == primer2.)
- **Hairpin**: an internal inverted repeat forms a stem of ≥ minStemLength bp separated by a
  loop ≥ minLoopLength (≥3 nt steric minimum). Watson-Crick antiparallel pairing across the stem.
- **3' stability**: ΔG°37 of the last 5 bases via SantaLucia NN sums + two initiation terms.

This is a **base-pair-counting / NN-ΔG heuristic**, which is the standard, honestly-labelled
approach (Primer3 itself reports an alignment/secondary-structure score, not a hard ΔG for SELF_ANY/END).

### Worked examples (hand computation)

**3' stability ΔG°37** (NN sums + init +0.98 G·C / +1.03 A·T):
- `GCGCG` = GC+CG+GC+CG = (−2.24−2.17−2.24−2.17) + 0.98 + 0.98 = **−6.86** (most stable 5-mer ✓ Primer3)
- `TATAT` = TA+AT+TA+AT = (−0.58−0.88−0.58−0.88) + 1.03 + 1.03 = **−0.86** (most labile 5-mer ✓ Primer3)
- `TACGT` = TA+AC+CG+GT = (−0.58−1.44−2.17−1.44) + 1.03 + 1.03 = **−3.57** ✓ (S6)
- `GGGGG` = GG×4 + 0.98×2 = −7.36 + 1.96 = **−5.40** ✓ (C5 problematic primer)

**3' self-dimer (positive case)** `AAAAAAAA` vs `AAAAAAAA`: revcomp(primer2)=`TTTTTTTT`; the 3'
window of primer1 (`AAAAAAAA`) pairs A-T at all 8 positions → 8 complementary bp ≥ 4 ⇒ dimer ✓.

**No 3' dimer (negative case)** `AAAACCCCCCCC` vs `GGGGGGGGTTTT`: revcomp(primer2)=`AAAACCCCCCCC`;
primer1 3' window `CCCCCCCC` vs window-start `AAAACCCC` → C pairs only with G, 0 complementary bp
< 4 ⇒ no dimer ✓.

**Hairpin (positive)** `ACGTACGTACGT` (stem=4, loop=3): 5' stem `ACGT` (pos 0–3) is antiparallel-
complementary to 3' stem `ACGT` (pos 8–11) with a 4-nt gap (≥ loop 3) ⇒ hairpin ✓.

**Hairpin (negative)** `AAAACCCCAAAA`: no `TTTT`/`GGGG` to pair the A/C stems ⇒ no hairpin ✓.

### Edge-case semantics

Empty/null → no structure (0 / false). Hairpin requires length ≥ 2·stem + loop. Stability requires
≥5 bp. Dinucleotide-repeat requires ≥4 bp. Case-insensitive (universal DNA convention). All defined and sourced.

### Findings / divergences

None. Description is biologically correct and matches authoritative sources.

---

## Stage B — Implementation

### Code path reviewed

`PrimerDesigner.cs`: `FindLongestHomopolymer` (246), `FindLongestDinucleotideRepeat` (273),
`HasHairpinPotential` (307) → `HasHairpinPotentialSimple` (328) / `HasHairpinPotentialWithSuffixTree`
(357), `HasPrimerDimer` (398), `Calculate3PrimeStability` (427), helpers `AreComplementary` (536),
`IsComplementary` (547), `Reverse` (529).

### Logic realised correctly? (evidence)

- **Complementary pairing** (`IsComplementary`, line 547): exactly Watson-Crick A-T / T-A / G-C / C-G.
- **Hairpin antiparallel stem** (line 338): checks `AreComplementary(fragment, Reverse(target))`.
  For 5' stem `fragment = seq[i..i+k]` and 3' stem `target = seq[j..j+k]`, antiparallel pairing
  requires `fragment[m]` to pair `target[k-1-m]`; `Reverse(target)[m] = target[k-1-m]`, so this is
  precisely correct antiparallel Watson-Crick pairing. Inner loop starts at
  `j = i + minStemLength + minLoopLength` (line 335) so the unpaired gap ≥ minLoopLength (≥3 default),
  honoring the steric loop minimum. Length guard `len ≥ 2·stem + loop` at line 309.
- **3'-end dimer emphasis** (lines 404–418): `seq2 = revcomp(primer2)`; inspects only the last
  `checkLength = min(8, …)` bases of primer1 vs the 3'-terminal window of primer2, counts
  complementary positions, returns count ≥ `minComplementarity` (default 4). Matches the validated
  3'-3' antiparallel self/cross-dimer model.
- **3' stability** (lines 437–472): SantaLucia (1998) NN ΔG°37 dictionary matches Table 1 values
  and the published stability trend exactly; adds the two initiation terms (+0.98 G·C / +1.03 A·T).
- **Case handling**: `ToUpperInvariant` applied in every method.

### Cross-verification table recomputed vs code (tests run)

| Case | Expected (source) | Code result |
|------|-------------------|-------------|
| GCGCG ΔG (M16) | −6.86 (Primer3/SantaLucia) | −6.86 ✓ |
| TATAT ΔG (M17) | −0.86 (Primer3/SantaLucia) | −0.86 ✓ |
| TACGT ΔG (S6) | −3.57 (SantaLucia) | −3.57 ✓ |
| GGGGG×4 ΔG (C5) | −5.40 (SantaLucia) | −5.40 ✓ |
| A₈/A₈ self-dimer (M13) | true | true ✓ |
| non-comp 3' (M12) | false | false ✓ |
| ACGTACGTACGT hairpin (M10) | true | true ✓ |
| AAAACCCCAAAA (M9) | false | false ✓ |
| AAAAA poly-A homopolymer / ACACACAC dinuc | 6 / 4 | 6 / 4 ✓ |

All values reproduced by the passing test suite.

### Variant/delegate consistency

Simple O(n²) path and suffix-tree path (≥100 bp) both validated (C3 tests). Suffix-tree builds on
`seq`, searches `revComp` patterns, and applies the same `j ≥ i + stem + loop` loop constraint
(both stem orientations), agreeing with the simple path on the tested long sequences.

### Test quality audit

`PrimerDesigner_PrimerStructure_Tests.cs` asserts **exact** sourced values (ΔG to ±0.01, exact
counts, exact booleans), not "no-throw" tautologies; covers all Stage-A edge cases (null/empty/short,
custom stem/loop/complementarity params, mixed case, long suffix-tree path). 42 test instances.

### Findings / defects

None. Implementation faithfully realises the validated description. It is a configurable
bp-counting (dimer/hairpin) + NN-ΔG (3' stability) heuristic, honestly labelled as such — CLEAN.

#### Minor non-blocking note (not a defect)
The dimer/hairpin detectors count complementary base pairs (with a 3'-window for dimers) rather than
computing a full alignment ΔG; this matches the spec's documented heuristic and Primer3's
score-based concept. The cross-spec note in the TestSpec (PRIMER-DESIGN-001 `stability3Prime < -9`
being unreachable since max 5-mer = −6.86) concerns PRIMER-DESIGN-001 and is out of scope here.

---

## Verdict & follow-ups

- **Stage A:** PASS — detection logic, loop minimum, 3'-end emphasis, Watson-Crick pairing, and
  SantaLucia/Primer3 ΔG values all confirmed against authoritative sources.
- **Stage B:** PASS — code matches the validated description; all worked examples reproduced.
- **End-state:** ✅ CLEAN — no defect; no code changes required.
- **Tests:** filter `~PrimerStructure` 42/42 passed; full suite `Seqeron.Genomics.Tests` 4461 passed, 0 failed.
