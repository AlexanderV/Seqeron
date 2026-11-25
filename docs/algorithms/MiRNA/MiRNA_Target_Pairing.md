# MiRNA-Target Pairing Analysis

| Field | Value |
|-------|-------|
| Algorithm Group | MiRNA |
| Test Unit ID | MIRNA-PAIR-001 |
| Related Projects | Seqeron.Genomics.Annotation |
| Implementation Status | Implemented |
| Last Reviewed | 2026-11-25 |

## 1. Overview

This unit computes the base-pairing structure of a microRNA (miRNA) bound antiparallel to a target mRNA segment. It reports which positions form Watson-Crick pairs (A-U, G-C), which form G:U wobble pairs, and which are mismatches, and produces a human-readable alignment string and the Turner 2004 nearest-neighbour free energy of that duplex. The base-pairing classification is exact and specification-driven (Watson-Crick + Crick's 1966 wobble rule); the free energy is the full NNDB Turner 2004 duplex model (initiation, stacks, loops, terminal AU/GU, dangling ends / terminal mismatches) evaluated by the canonical `RnaSecondaryStructure` kernels. It is used as the pairing primitive underlying miRNA target-site prediction. `FoldMiRnaDuplex` additionally computes the minimum-free-energy gapped duplex (bulges and internal loops allowed), the ViennaRNA `RNAduplex` model [7].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A mature miRNA (~22 nt) guides the silencing complex to mRNA targets primarily through Watson-Crick pairing between the miRNA 5′ "seed" (nucleotides 2–8) and complementary sites in the target 3′UTR [1][3]. The target site is the reverse complement of the seed read antiparallel [5]. Beyond standard Watson-Crick pairing (A-U, G-C), RNA tolerates the G-U "wobble" pair [4], which is weaker and counted separately from canonical pairs [3].

### 2.2 Core Model

For miRNA `m` (5′→3′) and target `t` (5′→3′) the antiparallel duplex pairs miRNA index `i` with target index `len(t)−1−i`. A position is classified by the pair (m[i], t[len(t)−1−i]):

- **Watson-Crick match** if {A-U, U-A, G-C, C-G} [3].
- **G:U wobble** if {G-U, U-G} [4].
- **Mismatch** otherwise [3] (A pairs only with U; C only with G [3]).

Duplex free energy follows the NNDB Turner 2004 nearest-neighbour model [6] for the structure in which every in-register Watson-Crick / G·U position is paired:

ΔG°37 = ΔG(intermolecular initiation, +4.09) [+ ΔG(symmetry, +0.43) for a duplex of two identical strands] + Σ ΔG(stack) (incl. the 5′GGUC/3′CUGG special stack) + Σ ΔG(internal loop / bulge) between helices + ΔG(terminal AU/GU, +0.45) per exterior helix end + ΔG(exterior stacking) at each duplex end — a terminal mismatch when the end is extended on both strands, otherwise the single 5′ or 3′ dangling end (NNDB: "if a helix end is extended on both the 5′ and 3′ strands, then a terminal mismatch exists (not the sum of 5′ and 3′ dangling ends)").

This is what ViennaRNA `eval_structure` computes for the same `miRNA&target` structure [7] (Vienna uses 0.1-kcal-rounded parameters, init 4.10, AU 0.50, and no symmetry term). With no base pair there is no duplex and ΔG = 0. Note that a short duplex can be destabilising: AAAA·UUUU = 4.09 − 3·0.93 + 2·0.45 = +2.20 kcal/mol.

`FoldMiRnaDuplex` minimises the same energy over all intermolecular structures (ViennaRNA `duplexfold` recursion: C(i,j) = init + exterior(i,j) or min over the next-outer pair (k,l) of C(k,l) + loop(k,l;i,j), loops ≤ 30 nt; ΔG = min C(i,j) + exterior on the inner end) [7].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Fixed 37 °C, 1 M NaCl (Turner 2004 conditions) | Free-energy magnitude shifts |
| ASM-02 | `AlignMiRnaToTarget`: ungapped fixed register (miRNA i · target L−1−i); every pairable position is paired | A bulged site is mis-registered — use `FoldMiRnaDuplex` for the MFE gapped duplex |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `CanPair` true ⇔ {A-U,U-A,G-C,C-G,G-U,U-G} | Watson-Crick [3] + wobble [4] |
| INV-02 | `IsWobblePair` true ⇔ {G-U,U-G}; Watson-Crick pairs are not wobble | G:U distinguished from WC [3][4] |
| INV-03 | Wobble pairs ⊆ pairable pairs | G-U is both a wobble and a valid pair [4] |
| INV-04 | `GetReverseComplement` reverses + complements (A↔U, G↔C, T→A, IUPAC R↔Y K↔M B↔V D↔H, S/W/N fixed; other symbols pass through), length-preserving | A-U/G-C complementarity [3]; reverse-complement seed [5]; Biopython `reverse_complement_rna` |
| INV-05 | Matches + Mismatches + GUWobbles = min(len(miRNA),len(target)); no gaps | every overlap position is classified into exactly one class |
| INV-06 | No base pair ⇒ FreeEnergy = 0; otherwise FreeEnergy = NNDB Turner 2004 duplex energy (e.g. NNDB worked examples AGCGCU·AGCGCU −7.94, GCACG·CGUGC −6.04) | NNDB helix / dangle / terminal-mismatch pages [6] |
| INV-07 | `FoldMiRnaDuplex(m,t).FreeEnergy` ≤ `AlignMiRnaToTarget(m,t).FreeEnergy` whenever the fixed-register structure has loops ≤ 30 nt (non-identical strands) | the fixed-register structure is one of the minimised structures [7] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| miRnaSequence | string | required | miRNA, 5′→3′ | DNA or RNA; T normalised to U; case-insensitive |
| targetSequence | string | required | target mRNA segment, 5′→3′ | DNA or RNA; T→U; case-insensitive |
| base1, base2 (CanPair/IsWobblePair) | char | required | two nucleotides | case-insensitive; T treated via U where relevant |
| rnaSequence (GetReverseComplement) | string | required | nucleotide sequence | DNA/RNA; T→A complement, output uses U; IUPAC complemented |
| maxLoopSize (FoldMiRnaDuplex) | int | 30 | largest bulge / internal loop (total unpaired nt) | ≥ 0 (ViennaRNA MAXLOOP = 30) |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| MiRnaDuplex.MiRnaSequence / TargetSequence | string | normalised (uppercase, RNA) inputs |
| MiRnaDuplex.AlignmentString | string | per miRNA position: `\|`=Watson-Crick, `:`=G:U wobble, space=mismatch |
| MiRnaDuplex.Matches / Mismatches / GUWobbles | int | counts over the overlap |
| MiRnaDuplex.Gaps | int | always 0 (ungapped) |
| MiRnaDuplex.FreeEnergy | double | Turner 2004 ΔG°37 (kcal/mol) of the reported duplex (§2.2); 0 if nothing pairs |
| MiRnaDuplexStructure | record | `FoldMiRnaDuplex` result: normalised sequences, 0-based duplex span (MiRnaStart..MiRnaEnd, TargetStart..TargetEnd), RNAduplex-notation Structure (`((..((&))..))`), BasePairs (miRNA, target) outermost first, FreeEnergy; null if no pair is possible |
| CanPair / IsWobblePair | bool | pairing / wobble predicate |
| GetReverseComplement | string | antiparallel RNA reverse complement |

### 3.3 Preconditions and Validation

Indexing is 0-based. Inputs are uppercased and `T`→`U` normalised (RNA alphabet). Pairing is antiparallel: miRNA index `i` pairs target index `len−1−i`. Empty/null `miRnaSequence` or `targetSequence` to `AlignMiRnaToTarget` returns an empty `MiRnaDuplex` (all counts 0, empty strings). Empty/null input to `GetReverseComplement` returns `""`. IUPAC codes are complemented by the canonical `SequenceExtensions.GetRnaComplementBase`; other symbols pass through. Non-ACGU symbols never pair and are not used as dangling / mismatch neighbours. No exceptions are thrown for these cases. `FoldMiRnaDuplex` returns null for null/empty input or when no pair is possible and throws `ArgumentOutOfRangeException` for a negative `maxLoopSize`.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalise both sequences (uppercase, T→U).
2. For each overlap index `i` in `0..min(len)−1`, pair miRNA[i] with target[len−1−i].
3. Classify the pair as Watson-Crick match, G:U wobble, or mismatch; emit the alignment symbol and increment the matching counter.
4. Evaluate the Turner 2004 energy of the paired positions (§2.2): helices = runs of consecutive pairs (`RnaSecondaryStructure.CalculateStemEnergy`, AU/GU end penalty kept only on exterior ends), loops between helices (`CalculateInternalLoopEnergy` / `CalculateBulgeLoopEnergy`), dangling ends / terminal mismatches (`GetDanglingEndEnergy` / `GetTerminalMismatchEnergy`), plus initiation (and symmetry for identical strands).
5. `FoldMiRnaDuplex`: fill C(i,j) over miRNA i ascending, target j descending; backtrack from the minimum; re-evaluate the structure with step 4.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Pairing table: Watson-Crick {A-U, G-C} [3] + G:U wobble [4].
- Energy terms: NNDB Turner 2004 (`Turner2004Parameters` + `RnaSecondaryStructure` kernels): intermolecular initiation 4.09 and symmetry 0.43 (NNDB helix page), stacks (Xia et al. 1998), terminal AU/GU 0.45, dangling-end and terminal-mismatch tables, internal-loop (incl. int11) and bulge models [6].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| AlignMiRnaToTarget | O(min(n,m)) | O(min(n,m)) | single antiparallel pass + alignment buffer |
| FoldMiRnaDuplex | O(n·m·L²) | O(n·m) | L = maxLoopSize (30) |
| GetReverseComplement | O(n) | O(n) | |
| CanPair / IsWobblePair | O(1) | O(1) | |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [MiRnaAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs)

- `MiRnaAnalyzer.AlignMiRnaToTarget(string, string)`: antiparallel ungapped duplex + Turner 2004 ΔG.
- `MiRnaAnalyzer.FoldMiRnaDuplex(string, string, int)`: MFE gapped duplex (RNAduplex model).
- `MiRnaAnalyzer.GetReverseComplement(string)`: RNA reverse complement via `SequenceExtensions.GetRnaComplementBase`.
- `MiRnaAnalyzer.CanPair(char, char)`: delegates to `RnaSecondaryStructure.CanPair`.
- `MiRnaAnalyzer.IsWobblePair(char, char)`: `RnaSecondaryStructure.GetBasePairType(a, b) == Wobble`.

### 5.2 Current Behavior

The aligner is ungapped and antiparallel over the overlap region (the shorter of the two sequences). `EvaluateDuplexEnergy` scores the paired positions with the full Turner 2004 duplex model through the canonical `RnaSecondaryStructure` kernels (§2.2, §4.1). Cross-check against ViennaRNA 2.8.0 `eval_structure` on 51 fixtures: max |Δ| 0.31 kcal/mol where no 2×2 internal loop occurs (Vienna's 0.1-kcal rounding); with 2×2 loops up to 2.2 kcal/mol because the canonical internal-loop kernel has no int22 table (RNA-STRUCT-001 residual, B12). `FoldMiRnaDuplex` vs `RNA.duplexfold`: identical base-pair set on 48/50 fixtures (the other two are near-degenerate alternatives within 0.3 kcal/mol); energies differ by the rounding residual, and by up to 1.4 kcal/mol where NNDB's 1-nt-bulge special-C bonus and degeneracy term apply (ViennaRNA omits both). This is a substring/positional scan, not a search: there is no occurrence enumeration of a pattern within a larger text, so the repository suffix tree is **not** applicable here (suffix-tree reuse evaluated; N/A — single positional alignment, not exact-match search).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Watson-Crick pairing A-U, G-C [3] and G:U wobble [4], with wobble counted separately from canonical matches [3].
- Antiparallel reverse-complement orientation of miRNA-to-target pairing [5].
- NNDB Turner 2004 duplex energy (initiation, symmetry, stacks, loops, AU/GU ends, dangling ends / terminal mismatches) [6]; NNDB worked examples reproduced exactly (−7.94, −6.04).
- MFE intermolecular duplex (ViennaRNA `duplexfold` recursion) as `FoldMiRnaDuplex` [7].

**By design (contract, not an approximation):**

- `AlignMiRnaToTarget` keeps the fixed ungapped register used by `FindTargetSites`; the gapped MFE duplex is `FoldMiRnaDuplex`.

**Not implemented:**

- Seed-context efficacy scoring and site-type classification: out of scope here; **users should rely on:** `MiRnaAnalyzer.FindTargetSites` (MIRNA-TARGET-001) and [Target_Site_Prediction](Target_Site_Prediction.md).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | 2×2 internal loops scored by the generic internal-loop formula (no int22 table in the canonical kernel) | Deviation (shared kernel) | up to 2.2 kcal/mol vs ViennaRNA on duplexes with 2×2 loops | open in B12 (RNA-STRUCT-001) | fixed centrally when the canonical kernel gains int22/int21 |
| 2 | Exterior dangling ends in the canonical MFE recursion use the swapped NNDB key | Deviation (shared kernel, not used here) | MFE only; this unit passes the keys in NNDB orientation (confirmed vs Vienna) | reported to B12 | see B13 report F13 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty/null miRNA or target | empty `MiRnaDuplex`, all counts 0 | defensive contract |
| Empty/null to GetReverseComplement | `""` | defensive contract |
| Unequal lengths | align over the shorter overlap | INV-05 |
| DNA input (T) | normalised to U before pairing | T↔U RNA convention |
| G:U pair | counted as wobble, not match | Crick (1966) [4]; G:U ≠ WC [3] |
| IUPAC / unknown base | complemented per IUPAC (R↔Y …) or passed through; never pairs | Biopython `reverse_complement_rna` |
| Bulged target site | `AlignMiRnaToTarget` loses register; `FoldMiRnaDuplex` finds the bulge | §2.2 |

### 6.2 Limitations

`AlignMiRnaToTarget` is a fixed-register ungapped duplex (use `FoldMiRnaDuplex` for bulges / asymmetric loops). Neither method models intramolecular structure of the target (use `CalculateSiteAccessibility`) or salt / temperature corrections. The A-opposite-position-1 preference is an Argonaute recognition event, not base pairing, and is therefore not represented here.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var duplex = MiRnaAnalyzer.AlignMiRnaToTarget("AAAA", "UUUU");
// duplex.Matches == 4, duplex.GUWobbles == 0, duplex.Mismatches == 0
// duplex.AlignmentString == "||||", duplex.FreeEnergy == +2.20 (4.09 − 3×0.93 + 2×0.45)
```

**Numerical / biological walk-through:**

`GGGG` vs `UUUU`: each miRNA G (i) pairs antiparallel with target U → G:U wobble at all 4 positions → AlignmentString `::::`, GUWobbles 4, Matches 0, Mismatches 0. `GetReverseComplement("GAGGUAG")` = `CUACCUC` (let-7a-5p seed reverse complement).

NNDB worked example (self-complementary duplex): `AlignMiRnaToTarget("AGCGCU", "AGCGCU")` → 4.09 + 2×0.45 + 0.43 − 2.08 − 3.42 − 2.36 − 3.42 − 2.08 = −7.94 kcal/mol. let-7a against its site with one extra A (`AACUAUACAACACUACUACCUCA`): `AlignMiRnaToTarget` −10.81 (register lost), `FoldMiRnaDuplex` −32.10 with a 1-nt bulge (ViennaRNA `duplexfold` same pairs, −31.8).

### 7.3 Related Tests, Evidence, or Documents

- Tests: [MiRnaAnalyzer_AlignMiRnaToTarget_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Annotation/MiRnaAnalyzer_AlignMiRnaToTarget_Tests.cs) — covers `INV-01`..`INV-06`
- Evidence: [MIRNA-PAIR-001-Evidence.md](../../../docs/Evidence/MIRNA-PAIR-001-Evidence.md)
- Related algorithms: [Target_Site_Prediction](Target_Site_Prediction.md), [Seed_Sequence_Analysis](Seed_Sequence_Analysis.md)

## 8. References

1. Bartel DP. 2009. MicroRNAs: Target Recognition and Regulatory Functions. Cell 136(2):215–233. https://www.cell.com/fulltext/S0092-8674(09)00008-7
2. Agarwal V, Bell GW, Nam JW, Bartel DP. 2015. Predicting effective microRNA target sites in mammalian mRNAs. eLife 4:e05005. https://pmc.ncbi.nlm.nih.gov/articles/PMC4532895/
3. Transcriptome-Wide Non-Canonical Interactions (review). PMC4870184. https://pmc.ncbi.nlm.nih.gov/articles/PMC4870184/
4. Crick FHC. 1966. Codon–anticodon pairing: the wobble hypothesis. J. Mol. Biol. 19(2):548–555 — via Wikipedia "Wobble base pair". https://en.wikipedia.org/wiki/Wobble_base_pair
5. Lewis BP, Burge CB, Bartel DP. 2005. Conserved seed pairing, often flanked by adenosines, indicates that thousands of human genes are microRNA targets. Cell 120(1):15–20. https://pubmed.ncbi.nlm.nih.gov/15652477/
6. Turner DH, Mathews DH. 2010. NNDB: the nearest neighbor parameter database. Nucleic Acids Res. 38:D280–D282. Turner 2004 pages opened 2026-11-25: https://rna.urmc.rochester.edu/NNDB/rna_2004/rna_2004_watson_crick_helices.html, …_dangling_ends.html, …_terminal_mismatches.html, …_exterior_loops.html, …_internal_loops.html; table `parameter_tables/turner_2004_rnastructure/turner_2004_dangle_dg.txt`.
7. Lorenz R et al. 2011. ViennaRNA Package 2.0. Algorithms Mol. Biol. 6:26. Reference implementation ViennaRNA 2.8.0 (`RNA.fold_compound(m+'&'+t).eval_structure`, `RNA.duplexfold`).
