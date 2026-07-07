# Canonical Algorithm Map

This map defines canonical algorithm IDs and documentation aliases to avoid duplicate units and split ownership.

## Canonical ID Aliases

| Alias ID | Canonical ID | Status | Note |
|---|---|---|---|
| `SEQ-COMPOSITION-001` | `SEQ-STATS-001` | Consolidated | Same canonical fixture and behavior |
| `SEQ-TM-001` | `SEQ-THERMO-001` | Consolidated | Same canonical fixture and behavior |
| `GENOMIC-TANDEM-001` | `REP-TANDEM-001` | Consolidated | Same method/class, unified under repeats |

## Taxonomy Aliases (Folder Naming)

| Alias / Variant | Canonical Bucket | Policy |
|---|---|---|
| `Molecular_Tools` | `MolTools` | Keep one owner bucket, cross-link from old location |
| `PopGen` | `Population_Genetics` | Keep short prefix in IDs only |
| `K-mer_Analysis` | `K-mer` | Merge docs under one k-mer bucket |
| `RNA_Secondary_Structure`, `RNA_Structure` | `RnaStructure` | Keep one RNA structure bucket |

## Canonicalization Rules

1. One concept must have exactly one canonical ID.
2. Alias IDs remain searchable but must point to canonical IDs.
3. One canonical behavior document per concept; domain pages should link rather than copy.
4. Legacy-baseline methods are retained and explicitly labeled `legacy-baseline`.

## Legacy / Baseline Methods

Retained on purpose as standard baselines/comparators — not garbage, but not the default choice.
Each doc carries a short baseline note linking back here.

| Method | Canonical Doc | Why keep | Why not default |
|---|---|---|---|
| UPGMA | [Tree_Construction.md](Phylogenetics/Tree_Construction.md) | Standard educational baseline for distance trees | Assumes a molecular clock; can bias topology/branch lengths (prefer Neighbor-Joining) |
| Jukes-Cantor / Kimura-2P | [Distance_Matrix.md](Phylogenetics/Distance_Matrix.md) | Canonical substitution corrections for distance matrices | Simplifying assumptions; not always best for real datasets |
| Chi-square Hardy-Weinberg | [Hardy_Weinberg_Test.md](Population_Genetics/Hardy_Weinberg_Test.md) | Fast, deterministic QC screen | Exact tests can be preferable for sparse counts |
| Nussinov-style classic RNA baseline | [Minimum_Free_Energy.md](RnaStructure/Minimum_Free_Energy.md) | Useful comparator for performance and behavior | Simplified energy model, not full thermodynamic fidelity |
| OLC assembly | [Overlap_Layout_Consensus.md](Assembly/Overlap_Layout_Consensus.md) | Important classical assembler paradigm | Often less practical than de Bruijn for many read regimes |

## Status

Folder normalization for the four canonical buckets above is **complete** — the alias folders
(`K-mer_Analysis`, `Molecular_Tools`, `PopGen`, `RNA_Structure`, `RNA_Secondary_Structure`) have
been merged into their canonical buckets and every inbound reference repointed. Remaining candidates
(`Motif_Analysis`, `Sequence_Comparison`, `Genomic_Analysis`, `Extended_Annotation`,
`Extended_Assembly`) are listed in the folder index and not yet consolidated.
