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

## Next Additions

- Add concept-level alias rows as cleanup continues.
- Add file-level redirects after folder normalization.
