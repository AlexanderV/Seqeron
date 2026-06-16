# Algorithm Taxonomy Cleanup Matrix

Date: 2026-06-16
Source: `ALGORITHMS_CHECKLIST_V2.md`, `docs/algorithms/*`

## Goal

Normalize algorithm naming and folder taxonomy, remove duplicate registry units, and keep legacy methods clearly marked as baseline/reference.

## Confirmed Checklist Duplicates

| Duplicate Unit | Canonical Unit | Evidence | Recommended Action |
|---|---|---|---|
| `SEQ-COMPOSITION-001` | `SEQ-STATS-001` | Checklist status says duplicate | Keep only canonical in summary tables; keep duplicate ID as alias for backward compatibility |
| `SEQ-TM-001` | `SEQ-THERMO-001` | Checklist status says duplicate | Same as above |
| `GENOMIC-TANDEM-001` | `REP-TANDEM-001` | Consolidated note + duplicate status | Same as above |

## Taxonomy Aliases (Folder Level)

| Current Folders | Canonical Folder | Issue | Recommended Action |
|---|---|---|---|
| `MolTools`, `Molecular_Tools` | `MolTools` | Split domain with overlapping intent | Keep `MolTools` as canonical; move docs from `Molecular_Tools`; leave redirect notes |
| `PopGen`, `Population_Genetics` | `Population_Genetics` | Same domain split by naming style | Keep `Population_Genetics` as canonical; keep `PopGen` only as short label in IDs |
| `K-mer`, `K-mer_Analysis` | `K-mer` | Operational and analytical k-mer docs split across two folders | Merge into `K-mer`; organize by subheadings in docs |
| `RNA_Secondary_Structure`, `RNA_Structure`, `RnaStructure` | `RnaStructure` | Triple naming variant for same area | Keep one folder only; add alias map in index |

## Concept-Level Duplicates in Docs

| Concept | Duplicate Paths | Canonical Doc |
|---|---|---|
| Relative synonymous codon usage | `Annotation/Relative_Synonymous_Codon_Usage.md`, `Codon/Relative_Synonymous_Codon_Usage.md` | `Codon/Relative_Synonymous_Codon_Usage.md` |
| Tandem repeat detection | `Repeat_Analysis/Tandem_Repeat_Detection.md`, `Genomic_Analysis/Tandem_Repeat_Detection.md` | `Repeat_Analysis/Tandem_Repeat_Detection.md` |
| Melting temperature | `Molecular_Tools/Melting_Temperature.md`, `Statistics/Melting_Temperature.md` | `MolTools` or `Statistics` (choose one owner and cross-link) |
| Low complexity region detection (protein) | `ProteinMotif/Low_Complexity_Region_Detection.md`, `ProteinPred/Low_Complexity_Region_Detection.md` | Keep one canonical behavior doc + one domain-specific usage note |

## Legacy/Baseline Methods To Keep (Not Remove)

These are not garbage; they should be explicitly tagged as baseline/reference in docs and checklist metadata.

| Method Family | Why Keep | Why Not Default |
|---|---|---|
| UPGMA | Standard educational baseline for distance trees | Assumes molecular clock; can bias topology/branch lengths |
| Jukes-Cantor / K2P | Canonical substitution corrections for distance matrices | Simplifying assumptions; not always best for real datasets |
| Chi-square Hardy-Weinberg | Fast, deterministic QC screen | Exact tests can be preferable for sparse counts |
| Nussinov-style classic RNA baseline | Useful comparator for performance and behavior | Simplified energy model, not full thermodynamic fidelity |
| OLC assembly | Important classical assembler paradigm | Often less practical than de Bruijn approaches for many read regimes |

## Proposed Canonicalization Rules

1. One algorithm concept = one canonical unit ID.
2. Legacy IDs stay as aliases with explicit mapping table.
3. One concept doc per method; domain pages should link to canonical docs rather than duplicate content.
4. Folder names use one convention only (PascalCase or snake_case, not mixed).
5. Add metadata tags in each doc frontmatter:
   - `status: canonical|alias|legacy-baseline`
   - `canonical_id: <ID>`
   - `owned_by: <domain>`

## Suggested Execution Order

1. Add `docs/algorithms/CANONICAL_MAP.md` with ID alias mapping.
2. Normalize folder taxonomy (no content rewrite yet).
3. Merge duplicated concept docs into canonical docs and replace old files with short redirect stubs.
4. Update `docs/algorithms/README.md` to reference only canonical sections.
5. Add legacy/baseline badges to UPGMA, JC/K2P, HWE-chi-square, Nussinov-classic, OLC docs.

## Minimal Canonical ID Alias Map (Initial)

| Alias ID | Canonical ID |
|---|---|
| `SEQ-COMPOSITION-001` | `SEQ-STATS-001` |
| `SEQ-TM-001` | `SEQ-THERMO-001` |
| `GENOMIC-TANDEM-001` | `REP-TANDEM-001` |

