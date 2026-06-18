# MiRNA-Target Pairing Analysis

| Field | Value |
|-------|-------|
| Algorithm Group | MiRNA |
| Test Unit ID | MIRNA-PAIR-001 |
| Related Projects | Seqeron.Genomics.Annotation |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-13 |

## 1. Overview

This unit computes the base-pairing structure of a microRNA (miRNA) bound antiparallel to a target mRNA segment. It reports which positions form Watson-Crick pairs (A-U, G-C), which form G:U wobble pairs, and which are mismatches, and produces a human-readable alignment string and an estimated duplex folding free energy. The base-pairing classification is exact and specification-driven (Watson-Crick + Crick's 1966 wobble rule); the free-energy value is a simplified nearest-neighbor stacking estimate. It is used as the pairing primitive underlying miRNA target-site prediction.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A mature miRNA (~22 nt) guides the silencing complex to mRNA targets primarily through Watson-Crick pairing between the miRNA 5′ "seed" (nucleotides 2–8) and complementary sites in the target 3′UTR [1][3]. The target site is the reverse complement of the seed read antiparallel [5]. Beyond standard Watson-Crick pairing (A-U, G-C), RNA tolerates the G-U "wobble" pair [4], which is weaker and counted separately from canonical pairs [3].

### 2.2 Core Model

For miRNA `m` (5′→3′) and target `t` (5′→3′) the antiparallel duplex pairs miRNA index `i` with target index `len(t)−1−i`. A position is classified by the pair (m[i], t[len(t)−1−i]):

- **Watson-Crick match** if {A-U, U-A, G-C, C-G} [3].
- **G:U wobble** if {G-U, U-G} [4].
- **Mismatch** otherwise [3] (A pairs only with U; C only with G [3]).

Duplex folding free energy is estimated with the Turner 2004 nearest-neighbor model [6]: ΔG ≈ Σ stacking(5′-m[i]m[i+1]-3′ / 3′-t[i]t[i+1]-5′) over consecutive paired positions, using the NNDB Turner 2004 stacking table. All paired-stack stacking free energies at 37 °C are negative (stabilising) [6].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Fixed 37 °C, standard buffer (Turner 2004 conditions) | Free-energy magnitude shifts; sign generally preserved |
| ASM-02 | Ungapped antiparallel duplex over the overlap region | Bulges/loops not modelled; ΔG of looped duplexes overestimated |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `CanPair` true ⇔ {A-U,U-A,G-C,C-G,G-U,U-G} | Watson-Crick [3] + wobble [4] |
| INV-02 | `IsWobblePair` true ⇔ {G-U,U-G}; Watson-Crick pairs are not wobble | G:U distinguished from WC [3][4] |
| INV-03 | Wobble pairs ⊆ pairable pairs | G-U is both a wobble and a valid pair [4] |
| INV-04 | `GetReverseComplement` reverses + complements (A↔U, G↔C, T→A), length-preserving | A-U/G-C complementarity [3]; reverse-complement seed [5] |
| INV-05 | Matches + Mismatches + GUWobbles = min(len(miRNA),len(target)); no gaps | every overlap position is classified into exactly one class |
| INV-06 | Fully Watson-Crick-paired duplex has FreeEnergy ≤ 0; all-mismatch duplex has FreeEnergy ≥ 0 | Turner 2004 paired-stack energies are negative; no paired stacks ⇒ sum 0 [6] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| miRnaSequence | string | required | miRNA, 5′→3′ | DNA or RNA; T normalised to U; case-insensitive |
| targetSequence | string | required | target mRNA segment, 5′→3′ | DNA or RNA; T→U; case-insensitive |
| base1, base2 (CanPair/IsWobblePair) | char | required | two nucleotides | case-insensitive; T treated via U where relevant |
| rnaSequence (GetReverseComplement) | string | required | nucleotide sequence | DNA/RNA; T→A complement, output uses U |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| MiRnaDuplex.MiRnaSequence / TargetSequence | string | normalised (uppercase, RNA) inputs |
| MiRnaDuplex.AlignmentString | string | per miRNA position: `\|`=Watson-Crick, `:`=G:U wobble, space=mismatch |
| MiRnaDuplex.Matches / Mismatches / GUWobbles | int | counts over the overlap |
| MiRnaDuplex.Gaps | int | always 0 (ungapped) |
| MiRnaDuplex.FreeEnergy | double | estimated ΔG (kcal/mol, 37 °C), nearest-neighbor stacking sum |
| CanPair / IsWobblePair | bool | pairing / wobble predicate |
| GetReverseComplement | string | antiparallel RNA reverse complement |

### 3.3 Preconditions and Validation

Indexing is 0-based. Inputs are uppercased and `T`→`U` normalised (RNA alphabet). Pairing is antiparallel: miRNA index `i` pairs target index `len−1−i`. Empty/null `miRnaSequence` or `targetSequence` to `AlignMiRnaToTarget` returns an empty `MiRnaDuplex` (all counts 0, empty strings). Empty/null input to `GetReverseComplement` returns `""`. Unrecognised bases complement to `N` and never pair. No exceptions are thrown for these cases.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalise both sequences (uppercase, T→U).
2. For each overlap index `i` in `0..min(len)−1`, pair miRNA[i] with target[len−1−i].
3. Classify the pair as Watson-Crick match, G:U wobble, or mismatch; emit the alignment symbol and increment the matching counter.
4. Sum Turner 2004 nearest-neighbor stacking energies over consecutive paired positions to estimate ΔG.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Pairing table: Watson-Crick {A-U, G-C} [3] + G:U wobble [4].
- Free-energy table: NNDB Turner 2004 nearest-neighbor stacking energies (`StackingEnergies`), the same set used by the pre-miRNA hairpin model; origin Xia et al. (1998) / Turner 2004 [6].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| AlignMiRnaToTarget | O(min(n,m)) | O(min(n,m)) | single antiparallel pass + alignment buffer |
| GetReverseComplement | O(n) | O(n) | |
| CanPair / IsWobblePair | O(1) | O(1) | |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [MiRnaAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs)

- `MiRnaAnalyzer.AlignMiRnaToTarget(string, string)`: antiparallel duplex base-pairing alignment + ΔG estimate.
- `MiRnaAnalyzer.GetReverseComplement(string)`: RNA reverse complement (T→U).
- `MiRnaAnalyzer.CanPair(char, char)`: Watson-Crick + G:U pairing predicate.
- `MiRnaAnalyzer.IsWobblePair(char, char)`: G:U wobble predicate.

### 5.2 Current Behavior

The aligner is ungapped and antiparallel over the overlap region (the shorter of the two sequences). `CalculateDuplexEnergy` sums Turner 2004 nearest-neighbor stacking only over runs of consecutive paired positions; mismatches break stacking but incur no explicit penalty term, and loop/bulge/initiation/terminal terms are omitted. This is a substring/positional scan, not a search: there is no occurrence enumeration of a pattern within a larger text, so the repository suffix tree is **not** applicable here (suffix-tree reuse evaluated; N/A — single positional alignment, not exact-match search).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Watson-Crick pairing A-U, G-C [3] and G:U wobble [4], with wobble counted separately from canonical matches [3].
- Antiparallel reverse-complement orientation of miRNA-to-target pairing [5].
- Nearest-neighbor stacking energies from the NNDB Turner 2004 table [6].

**Intentionally simplified:**

- Free energy: only nearest-neighbor stacking over consecutive paired positions is summed; **consequence:** ΔG magnitude is approximate (no loop, bulge, initiation, terminal-mismatch, or AU/GU-end terms), so only its sign and relative ordering are reliable, not its absolute kcal/mol value.
- Ungapped duplex; **consequence:** bulges and internal loops between paired regions are not represented.

**Not implemented:**

- Seed-context efficacy scoring and site-type classification: out of scope here; **users should rely on:** `MiRnaAnalyzer.FindTargetSites` (MIRNA-TARGET-001) and [Target_Site_Prediction](Target_Site_Prediction.md).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Turner stacking numeric values not re-retrieved this session | Assumption | Free-energy magnitude unverified this session; sign tested only | accepted | NNDB pages 404 to fetcher; values reused from prior committed unit; see Evidence ASSUMPTION 1 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty/null miRNA or target | empty `MiRnaDuplex`, all counts 0 | defensive contract |
| Empty/null to GetReverseComplement | `""` | defensive contract |
| Unequal lengths | align over the shorter overlap | INV-05 |
| DNA input (T) | normalised to U before pairing | T↔U RNA convention |
| G:U pair | counted as wobble, not match | Crick (1966) [4]; G:U ≠ WC [3] |
| Unknown base | complements to N; never pairs | no valid pairing defined |

### 6.2 Limitations

The free-energy magnitude is approximate (stacking-only); do not use it as a calibrated thermodynamic ΔG. The aligner is ungapped and does not model bulges, internal loops, or 3′ supplementary/compensatory pairing geometry. The A-opposite-position-1 preference is an Argonaute recognition event, not base pairing, and is therefore not represented here.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var duplex = MiRnaAnalyzer.AlignMiRnaToTarget("AAAA", "UUUU");
// duplex.Matches == 4, duplex.GUWobbles == 0, duplex.Mismatches == 0
// duplex.AlignmentString == "||||", duplex.FreeEnergy <= 0
```

**Numerical / biological walk-through:**

`GGGG` vs `UUUU`: each miRNA G (i) pairs antiparallel with target U → G:U wobble at all 4 positions → AlignmentString `::::`, GUWobbles 4, Matches 0, Mismatches 0. `GetReverseComplement("GAGGUAG")` = `CUACCUC` (let-7a-5p seed reverse complement).

### 7.3 Related Tests, Evidence, or Documents

- Tests: [MiRnaAnalyzer_AlignMiRnaToTarget_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/MiRnaAnalyzer_AlignMiRnaToTarget_Tests.cs) — covers `INV-01`..`INV-06`
- Evidence: [MIRNA-PAIR-001-Evidence.md](../../../docs/Evidence/MIRNA-PAIR-001-Evidence.md)
- Related algorithms: [Target_Site_Prediction](Target_Site_Prediction.md), [Seed_Sequence_Analysis](Seed_Sequence_Analysis.md)

## 8. References

1. Bartel DP. 2009. MicroRNAs: Target Recognition and Regulatory Functions. Cell 136(2):215–233. https://www.cell.com/fulltext/S0092-8674(09)00008-7
2. Agarwal V, Bell GW, Nam JW, Bartel DP. 2015. Predicting effective microRNA target sites in mammalian mRNAs. eLife 4:e05005. https://pmc.ncbi.nlm.nih.gov/articles/PMC4532895/
3. Transcriptome-Wide Non-Canonical Interactions (review). PMC4870184. https://pmc.ncbi.nlm.nih.gov/articles/PMC4870184/
4. Crick FHC. 1966. Codon–anticodon pairing: the wobble hypothesis. J. Mol. Biol. 19(2):548–555 — via Wikipedia "Wobble base pair". https://en.wikipedia.org/wiki/Wobble_base_pair
5. Lewis BP, Burge CB, Bartel DP. 2005. Conserved seed pairing, often flanked by adenosines, indicates that thousands of human genes are microRNA targets. Cell 120(1):15–20. https://pubmed.ncbi.nlm.nih.gov/15652477/
6. Turner DH, Mathews DH. 2010. NNDB: the nearest neighbor parameter database. Nucleic Acids Res. 38:D280–D282. Turner 2004 parameters: https://rna.urmc.rochester.edu/NNDB/turner04/index.html
