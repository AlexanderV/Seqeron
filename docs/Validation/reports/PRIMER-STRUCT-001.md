# Validation Report: PRIMER-STRUCT-001 — Primer Structure Analysis

- **Validated:** 2026-06-24   **Area:** Molecular Tools (Seqeron.Genomics.MolTools)
- **Canonical method(s):**
  `PrimerDesigner.HasHairpinPotential(seq, minStemLength, minLoopLength)`,
  `PrimerDesigner.HasPrimerDimer(primer1, primer2, minComplementarity)`,
  `PrimerDesigner.Calculate3PrimeStability(seq)`,
  `PrimerDesigner.FindLongestHomopolymer(seq)`, `PrimerDesigner.FindLongestDinucleotideRepeat(seq)`
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs`
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/PrimerDesigner_PrimerStructure_Tests.cs`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN

This is an independent re-validation (fresh context). It confirms the prior report
(`git show cb113ce:docs/Validation/reports/PRIMER-STRUCT-001.md`) and additionally re-derives
every SantaLucia NN ΔG value from the primary paper.

---

## Stage A — Description

### Sources opened & what they confirm

- **SantaLucia (1998)** PNAS 95:1460-65, Table 1 unified NN parameters (1 M NaCl), fetched from
  PMC (https://pmc.ncbi.nlm.nih.gov/articles/PMC19045/). Extracted ΔG°37 (kcal/mol):
  AA/TT −1.00, AT −0.88, TA −0.58, CA/TG −1.45, GT/AC −1.44, CT/AG −1.28, GA/TC −1.30,
  CG −2.17, GC −2.24, GG/CC −1.84; initiation terminal G·C **+0.98**, terminal A·T **+1.03**.
  **Every value matches the implementation dictionary exactly** (including the symmetric
  TT=AA, TG=CA, AC=GT, AG=CT, TC=GA, CC=GG entries).
- **Primer3 manual / PRIMER_MAX_END_STABILITY**: max ΔG (kcal/mol) for disrupting the duplex of
  the **five 3' bases**, NN-computed; published anchors "most stable 5-mer = 6.86 (GCGCG)",
  "most labile 5-mer = 0.86 (TATAT)". Confirms the 5-mer 3'-stability ΔG concept and the two
  reference values used by the tests.
- **Wikipedia — Stem-loop**: hairpin = self-complementary inverted repeat (stem) + loop;
  *"loops fewer than three bases long are sterically impossible"* → min loop ≈ 3 nt. Validates the
  `minLoopLength` default and the steric minimum.
- **Wikipedia — Primer dimer**: a primer-dimer forms when two primers anneal **at their 3' ends**
  via complementary bases; the 3' emphasis is the deleterious case (polymerase extends). Validates
  the 3'-window dimer model; self-dimer is the special case primer1 == primer2.

### Detection logic validated

- **3' stability**: ΔG°37 of the last 5 bases = sum of 4 NN terms + 2 initiation terms.
- **Self/cross-dimer**: align primer1's 3' window antiparallel against revcomp(primer2);
  count Watson-Crick pairs in the 3'-terminal window; "dimer" when count ≥ threshold.
- **Hairpin**: internal inverted repeat — 5' stem of ≥minStemLength pairs antiparallel with a 3'
  stem, separated by an unpaired gap ≥ minLoopLength (≥3).
- **Homopolymer / dinucleotide repeat**: longest run / longest 2-mer microsatellite (Primer3
  PRIMER_MAX_POLY_X / repeat avoidance).

This is a configurable base-pair-counting (dimer/hairpin) + NN-ΔG (3' stability) heuristic,
honestly labelled — the same family of metric Primer3 reports (a score, not a hard ΔG, for
SELF_ANY/SELF_END).

### Worked examples (hand computation, against the re-derived Table 1)

- `GCGCG` = (−2.24−2.17−2.24−2.17) + 0.98 + 0.98 = **−6.86** (Primer3 most-stable 5-mer ✓)
- `TATAT` = (−0.58−0.88−0.58−0.88) + 1.03 + 1.03 = **−0.86** (Primer3 most-labile 5-mer ✓)
- `TACGT` = (−0.58−1.44−2.17−1.44) + 1.03 + 1.03 = **−3.57** ✓ (S6)
- `GGGGG` = (−1.84×4) + 0.98×2 = **−5.40** ✓ (C5 problematic primer)
- Self-dimer `AAAAAAAA`/`AAAAAAAA`: revcomp(p2)=`TTTTTTTT`; 8 A-T pairs ≥4 ⇒ dimer ✓.
- Non-dimer `AAAACCCCCCCC`/`GGGGGGGGTTTT`: revcomp(p2)=`AAAACCCCCCCC`; 3' `CCCCCCCC` vs `AAAACCCC` → 0 pairs ⇒ no dimer ✓.
- Hairpin `ACGTACGTACGT` (stem 4, loop 4): `ACGT`(0–3) vs Reverse(`ACGT`@8–11)=`TGCA` → A-T,C-G,G-C,T-A ⇒ hairpin ✓.

### Edge-case semantics

Empty/null → 0/false. Hairpin needs len ≥ 2·stem+loop; stability needs ≥5 bp; dinuc needs ≥4 bp.
Case-insensitive everywhere. All defined and sourced.

### Findings / divergences

None. Description is biologically/thermodynamically correct.

---

## Stage B — Implementation

### Code path reviewed (file:line)

`PrimerDesigner.cs`: `FindLongestHomopolymer` (246), `FindLongestDinucleotideRepeat` (273),
`HasHairpinPotential` (307) → `HasHairpinPotentialSimple` (328) / `HasHairpinPotentialWithSuffixTree`
(357), `HasPrimerDimer` (398), `Calculate3PrimeStability` (427), helpers `Reverse` (621),
`AreComplementary` (628), `IsComplementary` (639).

### Logic realised correctly? (evidence)

- **NN dictionary** (437–455): every key/value matches the re-derived SantaLucia Table 1.
  Initiation (469–470): `+0.98` if terminal base is G/C else `+1.03`, applied to both 5'/3' termini.
- **3' window** (433): `last5 = seq[^5..]`, 4 NN sums (458–463) + 2 inits → matches worked values.
- **Dimer** (404–418): `seq2 = revcomp(primer2)`, `checkLength = min(8, …)`, counts WC pairs of
  primer1's 3'-terminal window vs primer2's; `count ≥ minComplementarity` (default 4). Realises the
  validated 3'-3' model.
- **Hairpin (simple)**: inner loop starts at `j = i + minStemLength + minLoopLength` (335) so gap ≥
  loop; `AreComplementary(fragment, Reverse(target))` (338) = exact antiparallel Watson-Crick.
- **`IsComplementary`** (639): exactly A-T/T-A/G-C/C-G.
- **Dinucleotide** boundary `while (j+1 < len && …)` is correct (reading `Substring(j,2)` needs
  `j+1 < len`); `ACACACAC` → 4 traced by hand. Overlapping start index `i` ensures any frame is found.
- Case: `ToUpperInvariant` / `char.ToUpperInvariant` in every method.

### Cross-verification table recomputed vs code (tests run)

| Case | Expected (source) | Code result |
|------|-------------------|-------------|
| GCGCG ΔG (M16) | −6.86 (Primer3/SantaLucia) | −6.86 ✓ |
| TATAT ΔG (M17) | −0.86 (Primer3/SantaLucia) | −0.86 ✓ |
| TACGT ΔG (S6)  | −3.57 (SantaLucia) | −3.57 ✓ |
| GGGGG×4 ΔG (C5) | −5.40 (SantaLucia) | −5.40 ✓ |
| A₈/A₈ self-dimer (M13) | true | true ✓ |
| non-comp 3' (M12) | false | false ✓ |
| ACGTACGTACGT hairpin (M10) | true | true ✓ |
| AAAACCCCAAAA (M9) | false | false ✓ |
| AAAAAA / ACACACAC | homopolymer 6 / dinuc 4 | 6 / 4 ✓ |

All reproduced by the passing suite.

### Variant/delegate consistency

Simple O(n²) and suffix-tree (≥100 bp) paths both pass the C3 long-sequence tests (positive +
negative). Suffix-tree maps revComp position `p` to `j = n − p − minStemLength` and applies the same
`j ≥ i + stem + loop` constraint in both stem orientations, agreeing with the simple path on tests.

### Test quality audit

`PrimerDesigner_PrimerStructure_Tests.cs` (42 instances) asserts **exact** sourced values
(ΔG ±0.01, exact counts, exact booleans), not no-throw tautologies; covers all Stage-A edge cases
(null/empty/short, custom stem/loop/complementarity, mixed case, long suffix-tree path).

### Findings / defects

None. Implementation faithfully realises the validated description.

---

## Verdict & follow-ups

- **Stage A:** PASS — NN ΔG values re-derived from SantaLucia (1998) Table 1 match exactly;
  initiation terms, Primer3 5-mer anchors (6.86/0.86), 3' dimer emphasis and ≥3 nt loop all confirmed.
- **Stage B:** PASS — code matches the validated description; all worked examples reproduced.
- **End-state:** ✅ CLEAN — no defect; no code changes required.
- **Tests:** filter `~PrimerDesigner_PrimerStructure_Tests` = 42/42 passed.
- **Non-blocking note (out of scope):** the TestSpec's cross-spec remark that PRIMER-DESIGN-001's
  `stability3Prime < -9` threshold is unreachable (max 5-mer = −6.86) concerns PRIMER-DESIGN-001.
