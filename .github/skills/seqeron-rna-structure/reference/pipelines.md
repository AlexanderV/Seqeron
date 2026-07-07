# seqeron-rna-structure — fuller pipelines & parameter guidance

Dual-mode recipes for the RNA-secondary-structure family (analysis server, `RnaSecondaryStructure.*`).
Rigor (tool-only, provenance, cross-check, alpha caveat) is delegated to
[`bio-rigor`](../../bio-rigor/SKILL.md). Schemas: `docs/mcp/tools/analysis/<tool>.md`. Tool index:
[`tool-map.md`](tool-map.md).

## Parameter defaults (verify per tool doc)

- `minStemLength` = 3 · `minLoopSize` = 3 (values < 3 clamped) · `maxLoopSize` = 10 · `allowWobble` = true.
- `find_rna_inverted_repeats`: `minLength` = 4 arm, `minSpacing` = 3, `maxSpacing` = 100 (loop window).
- Energies are **kcal/mol at 37 °C**, Turner 2004 nearest-neighbour. `minimum_free_energy` returns
  **0** when no structure can form (not an error).
- Base-pair positions are **0-based**; dot-bracket pairs are `(openingIndex, closingIndex)`.

## 1. Fold an RNA (main entry) + independent MFE cross-check

**Goal.** Turn a sequence into a structure and a trustworthy energy.

1. **[MCP]** `predict_rna_structure`(rnaSequence, minStemLength=3, minLoopSize=3, maxLoopSize=10)
   → `dotBracket`, `basePairs`, `stemLoops`, `pseudoknots`, `mfe`.
2. **[MCP]** `minimum_free_energy`(rnaSequence, minLoopSize=3) → DP MFE (Zuker-style, Turner 2004).
3. Cross-check: the two energies should have consistent sign/magnitude. They can differ — prediction
   is **greedy nested-only stem-loops**, so its total is an *upper bound* on the DP MFE, not equal to it.

- **[C# API]** `RnaSecondaryStructure.PredictStructure(...)` → `.CalculateMinimumFreeEnergy(...)`.

```
Provenance
1) predict_rna_structure(rnaSequence,minStemLength=3,minLoopSize=3,maxLoopSize=10) → dotBracket,basePairs,stemLoops,pseudoknots,mfe
2) minimum_free_energy(rnaSequence,minLoopSize=3) → mfe (independent DP)
Cross-check: DP mfe ≤ greedy total (greedy is nested-only, non-optimal).
Envelope: RNA-STRUCT-001 documented-only — NN-model result, not a guaranteed native fold.
Caveat: alpha — validate before use.
```

## 2. Stem-loop / hairpin enumeration + arm corroboration

**Goal.** Find candidate hairpins and independently confirm the stem arms.

1. **[MCP]** `find_stem_loops`(rnaSequence, minStemLength=3, minLoopSize=3, maxLoopSize=10, allowWobble=true)
   → per candidate: stem 5'/3' bounds, base pairs, stem energy; loop type/bounds/size/sequence;
   dot-bracket; Turner total energy.
2. **[MCP]** `find_rna_inverted_repeats`(sequence, minLength=4, minSpacing=3, maxSpacing=100)
   → left-arm / right-arm spans + arm length (antiparallel complementary — RNA rules).
3. Cross-check: a genuine hairpin from step 1 should have its stem arms overlapped by an inverted
   repeat from step 2 (two independent detectors of the same complementary duplex).

- **[C# API]** `RnaSecondaryStructure.FindStemLoops(...)` · `.FindInvertedRepeats(...)`.

```
Provenance
1) find_stem_loops(seq,minStemLength=3,minLoopSize=3,maxLoopSize=10,allowWobble=true) → stems,loops,energies,dotBracket
2) find_rna_inverted_repeats(seq,minLength=4,minSpacing=3,maxSpacing=100) → arm spans
Cross-check: step-1 stems overlap step-2 inverted-repeat arms.
Envelope: RNA-STRUCT-001 (model floor). Caveat: alpha.
```

## 3. Pseudoknot detection (from any base-pair set)

**Goal.** Decide whether a structure contains crossing pairs.

1. Obtain base pairs: either **[MCP]** `parse_dot_bracket`(dotBracket) → `(open,close)` tuples
   (map to `{position1,position2,base1,base2,type}` using the sequence), or take `basePairs` straight
   from `predict_rna_structure`.
2. **[MCP]** `detect_pseudoknots`(basePairs) → crossing pair-of-pairs where `i<k<j<l`; nested
   (`i<k<l<j`) and disjoint (`j<k`) excluded. Each crossing pair reported once.

- **[C# API]** `RnaSecondaryStructure.ParseDotBracket(...)` → `.DetectPseudoknots(pairs)`.
- **STOP (RNA-STRUCT-001):** `detect_pseudoknots` is exact only for the **csr-PK grammar**. Kissing
  hairpins / loop-loop, triple-crossing / chained, and non-canonical bulged-helix pseudoknots are
  **not** faithfully detectable (unmeasured loop-loop energy). If asked to find those classes, report
  the limitation ([`../../../../docs/Validation/LIMITATIONS.md`](../../../../docs/Validation/LIMITATIONS.md))
  rather than presenting the output as complete.

```
Provenance
1) parse_dot_bracket(dotBracket) → (open,close) tuples  [or basePairs from predict_rna_structure]
2) detect_pseudoknots(basePairs) → crossing pairs (i<k<j<l)
Envelope: RNA-STRUCT-001 — only csr-PK crossings in scope; non-csr-PK classes not detectable.
Caveat: alpha.
```

## 4. Dot-bracket QC (validate → parse)

**Goal.** Trust an externally supplied structure string before using it.

1. **[MCP]** `validate_dot_bracket`(dotBracket) → balanced & family-matched? A mispairing like `(]`
   is rejected even if bracket counts balance. **Do not parse an invalid string.**
2. **[MCP]** `parse_dot_bracket`(dotBracket) → `(openingIndex, closingIndex)` tuples; multi-family
   (`()[]{}<>` + upper/lower letter pairs), unpaired `.,-:_~` skipped.

- **[C# API]** `RnaSecondaryStructure.ValidateDotBracket(...)` → `.ParseDotBracket(...)`.

```
Provenance
1) validate_dot_bracket(dotBracket) → true
2) parse_dot_bracket(dotBracket) → base-pair index tuples
Caveat: alpha. (No guarded unit here — dot-bracket handling is exact.)
```

## 5. Manual base-pairing & energy-term inspection

**Goal.** Probe individual pairings or motif energies (e.g. building/verifying a hand model).

- **Pairing:** **[MCP]** `can_pair`(base1, base2) → bool · `base_pair_type`(base1, base2) →
  WatsonCrick / Wobble / null · `rna_complement_base`(base) → complement.
- **Energy terms (Turner 2004, kcal/mol):** **[MCP]** `stem_energy` · `hairpin_loop_energy` ·
  `internal_loop_energy` · `bulge_loop_energy` · `multibranch_loop_energy` · `dangling_end_energy` ·
  `terminal_mismatch_energy` (see [`tool-map.md`](tool-map.md) for each signature/notes).
- **[C# API]** `RnaSecondaryStructure.CanPair` · `.GetBasePairType` · `.GetComplement` ·
  `.CalculateStemEnergy` · `.CalculateHairpinLoopEnergy` · `.CalculateInternalLoopEnergy` ·
  `.CalculateBulgeLoopEnergy` · `.CalculateMultibranchLoopEnergy` · `.GetDanglingEndEnergy` ·
  `.GetTerminalMismatchEnergy`.

```
Provenance
e.g. base_pair_type(base1=G,base2=U) → Wobble ; hairpin_loop_energy(...) → kcal/mol
Caveat: alpha; energies are Turner 2004 NN at 37 °C — not experimentally re-measured here.
```

## Scope reminders

- **No base-pair-probability / partition-function MCP tool.** McCaskill partition function is a
  documented C#-only algorithm (`docs/algorithms/RnaStructure/RNA_Partition_Function.md`,
  `Turner_McCaskill_Partition_Function.md`); there is no exposed tool — do not fabricate one.
- **Pre-miRNA hairpins, miRNA seed/target pairing, RNA-seq, splicing** → [`bio-annotation`](../../bio-annotation/SKILL.md).
- **DNA inverted repeats / palindromes** are a different tool (`find_inverted_repeats`, DNA rules) →
  `bio-assembly` / `bio-annotation`.
