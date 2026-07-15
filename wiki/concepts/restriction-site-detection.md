---
type: concept
title: "Restriction site detection (both-strand IUPAC recognition scan + cut-position computation)"
tags: [moldesign, restriction]
mcp_tools:
  - find_restriction_sites
  - find_all_restriction_sites
  - get_enzyme
sources:
  - docs/algorithms/MolTools/Restriction_Site_Detection.md
source_commit: f49a7c3f60f774ce1c9bdb64647502119d973aaa
created: 2026-07-14
updated: 2026-07-14
---

# Restriction site detection

The **RESTR** family's **location** primitive — given a DNA sequence and an enzyme,
find *where* the enzyme's recognition sequence occurs and *where* it cuts. This is the
operation that both siblings depend on but neither owns: [[restriction-enzyme-filtering]]
selects enzymes from a library by cut properties (no sequence), and
[[restriction-digest-simulation]] consumes the forward-strand cut positions detection
produces to partition a sequence into fragments. Detection sits between them —
sequence-in, per-site records out. Implemented on `RestrictionAnalyzer`
(`Seqeron.Genomics.MolTools`); status **Simplified**. A primary per-algorithm spec
(Test Unit ID N/A), synthesized against the
[[algorithm-validation-evidence|validation-evidence]] pattern used across MolTools.

## Core model — both-strand IUPAC scan

For each enzyme, slide a window of the recognition-sequence length `m` across the input
and test each window character-by-character under **IUPAC ambiguity rules**
(`IupacHelper.MatchesIupac`), so degenerate motifs (e.g. `SfiI`, `HincII`) match their
ambiguity codes. A **forward-strand** match at start `i` has cut position:

```
cutPosition = i + CutPositionForward
```

The **reverse strand** is searched by scanning the reverse complement, then mapping each
reverse-frame start `i` back to forward coordinates:

```
forwardPos   = sequenceLength - i - patternLength
cutPosition  = forwardPos + CutPositionReverse
```

Both strands are reported independently, so a **palindromic** site (which reads the same
on both strands) surfaces **twice at the same position** with different `IsForwardStrand`
flags. This is the double-report that [[restriction-digest-simulation]] deliberately
collapses (it keeps forward-strand cuts only) so fragment counts are not doubled.

## Contract and output

Each detected site is a `RestrictionSite`:

| Field | Type | Meaning |
|-------|------|---------|
| `Position` | `int` | Start of the recognition sequence (forward coordinates) |
| `Enzyme` | `RestrictionEnzyme` | Enzyme definition that matched |
| `IsForwardStrand` | `bool` | `true` forward match, `false` reverse-strand match |
| `CutPosition` | `int` | Computed cut position (per formulas above) |
| `RecognizedSequence` | `string` | The actual matched segment |

Entry points: `FindSites(DnaSequence, string)` (one enzyme), `FindSites(DnaSequence,
params string[])` (union of per-enzyme results), `FindSites(string, string)` (raw-string
overload — uppercases input, empty/null yields no sites), `FindAllSites(DnaSequence)`
(scans the full built-in catalog), and `GetEnzyme(string)` (case-insensitive lookup).

## Overhang classification

The spec derives overhang type from the two cut positions — the same
blunt/sticky distinction [[restriction-enzyme-filtering]] filters on, here read off the
enzyme's own cut offsets:

| Type | Condition | Example |
|------|-----------|---------|
| 5' overhang | `CutPositionForward < CutPositionReverse` | EcoRI |
| 3' overhang | `CutPositionForward > CutPositionReverse` | PstI |
| Blunt end | `CutPositionForward == CutPositionReverse` | EcoRV |

Recognition-length categories: 4-cutters (AluI, MspI, TaqI), 6-cutters (EcoRI, BamHI,
HindIII), 8-cutters (NotI, PacI, AscI).

## Key invariants and contract

- **In-bounds position:** `0 ≤ Position ≤ sequence.Length − recognitionLength` for every
  reported site (yielded only from valid window starts).
- **Exact slice:** `RecognizedSequence.Length == Enzyme.RecognitionSequence.Length`.
- **Case-insensitive lookup:** the built-in enzyme dictionary uses
  `StringComparer.OrdinalIgnoreCase`.
- **Empty raw-string input → no sites** (raw-string overload short-circuits on null/empty).
- **Unknown enzyme name → `ArgumentException`**; null `DnaSequence` → throws; null custom
  `enzyme` → `ArgumentNullException`.
- **Palindrome double-report:** both strands searched independently, so a palindromic
  site appears on both strands at the same coordinate — must be de-duplicated downstream
  when unique cut counts are required.
- **Complexity:** single-enzyme `O(n·m)` streaming `O(1)`; multi-enzyme /
  `FindAllSites` `O(n·k·m)` for `k` enzymes.

## Scope and limitations

Sequence-based matching against a **built-in enzyme catalog** only — enzymes absent from
the catalog are unavailable by name (custom-enzyme overload aside). It does **not** model
methylation sensitivity, assay conditions, or digestion efficiency; those require external
enzymology references. Double-reporting of palindromic sites must be interpreted downstream.

## Sources

`docs/algorithms/MolTools/Restriction_Site_Detection.md` (spec), which cites Wikipedia
*Restriction enzyme / Restriction site / EcoRI*, Roberts (1976) *Restriction
endonucleases*, REBASE, and IUPAC-IUB (1970) nucleic-acid notation. See
[[restriction-enzyme-filtering]] (enzyme selection) and
[[restriction-digest-simulation]] (fragmentation) for the sibling RESTR units.
