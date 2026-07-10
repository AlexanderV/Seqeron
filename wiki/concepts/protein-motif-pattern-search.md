---
type: concept
title: "Protein motif search by pattern (general PROSITE→regex engine + IC scoring)"
tags: [analysis, algorithm, protein, motif]
mcp_tools:
  - find_motif_by_pattern
  - find_motif_by_prosite
  - prosite_to_regex
sources:
  - docs/Evidence/PROTMOTIF-FIND-001-Evidence.md
  - docs/Evidence/PROTMOTIF-PATTERN-001-Evidence.md
  - docs/Evidence/PROTMOTIF-PROSITE-001-Evidence.md
source_commit: 0908c5f04255fcb3d51c7706d74a689eee481faa
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: protmotif-find-001-evidence
      evidence: "Test Unit ID: PROTMOTIF-FIND-001 ... Algorithm: Protein Motif Search (Pattern-based)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:common-protein-motifs
      source: protmotif-find-001-evidence
      evidence: "This unit is the general FindMotifByPattern engine (arbitrary caller-supplied PROSITE pattern → regex via ConvertPrositeToRegex); FindCommonMotifs is the fixed-dictionary application over the same CommonMotifs catalog documented in the same Evidence file"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: protmotif-pattern-001-evidence
      evidence: "Test Unit ID: PROTMOTIF-PATTERN-001 ... Algorithm: Protein Pattern Matching Methods (FindMotifByPattern, FindMotifByProsite, ConvertPrositeToRegex, FindDomains)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:protein-domain-and-signal-peptide-prediction
      source: protmotif-pattern-001-evidence
      evidence: "PROTMOTIF-PATTERN-001 groups FindDomains with the pattern-matching methods; FindDomains is the deterministic PROSITE-pattern domain scan covered by protein-domain-and-signal-peptide-prediction"
      confidence: medium
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: protmotif-prosite-001-evidence
      evidence: "Test Unit ID: PROTMOTIF-PROSITE-001 ... Algorithm: PROSITE Pattern Matching (ConvertPrositeToRegex, FindMotifByProsite)"
      confidence: high
      status: current
---

# Protein motif search by pattern (general PROSITE→regex engine + IC scoring)

**Pattern-based protein motif search** is the **general, caller-supplies-the-pattern**
protein motif engine: the user hands in an arbitrary **PROSITE pattern** (or an equivalent
regex) and it reports every occurrence in an amino-acid sequence. Seqeron exposes it as
`ProteinMotifFinder.FindMotifByPattern`. The fixed-dictionary scanner
[[common-protein-motifs|`FindCommonMotifs`]] is just a curated application of this engine
over a built-in catalog — this page is the **engine underneath it**. Validated under test
unit **PROTMOTIF-FIND-001** ("Protein Motif Search (Pattern-based)") and revalidated by
**PROTMOTIF-PATTERN-001** ("Protein Pattern Matching Methods") and **PROTMOTIF-PROSITE-001**
("PROSITE Pattern Matching"); the validation records are [[protmotif-find-001-evidence]],
[[protmotif-pattern-001-evidence]] and [[protmotif-prosite-001-evidence]], and
[[test-unit-registry]] tracks the units. See [[algorithm-validation-evidence]] for the
artifact pattern.

Within the ProteinMotif family it is the general engine, sibling of the fixed-catalog
[[common-protein-motifs]], the windowed [[coiled-coil-prediction]], and the longer-signature
[[protein-domain-and-signal-peptide-prediction]] — all on `ProteinMotifFinder`. It is the
**protein** analogue of the caller-supplies-the-query DNA matchers (contrast the DNA fixed
catalog [[regulatory-element-detection]]).

## The three primitives

| Primitive | Role |
|-----------|------|
| `ConvertPrositeToRegex` | translate a PROSITE pattern into a .NET regex |
| `FindMotifByPattern` | run the regex against a sequence and emit all `MotifMatch` hits (incl. overlaps) |
| `FindMotifByProsite` | end-to-end convenience path: PROSITE string → `ConvertPrositeToRegex` → `FindMotifByPattern` (PROTMOTIF-PATTERN-001 pins it on PS00001 and PS00016) |
| `CalculateMotifScore` / `CalculateEValue` | information-content score and expected random-match count for a motif |

## PROSITE → regex conversion

`ConvertPrositeToRegex` maps each PROSITE syntax element to its regex equivalent (PROSITE
User Manual `PA` line rules):

| PROSITE element | Meaning | Regex |
|-----------------|---------|-------|
| `x` | any amino acid | `.` |
| `[ABC]` | ambiguity — one of the set | `[ABC]` |
| `{ABC}` | exclusion — any except the set | `[^ABC]` |
| `x(n)` | fixed repetition | `.{n}` |
| `x(n,m)` | range repetition | `.{n,m}` |
| `A(n)` | fixed count on a residue letter | `A{n}` |
| `<` | N-terminus anchor | `^` |
| `>` | C-terminus anchor | `$` |
| `[G>]` | residue **or** end-of-sequence (rare; PS00267/PS00539) | `(?:G\|$)` |
| `-` | element separator | (removed) |
| `.` | pattern terminator (trailing or mid-pattern) | (ends the PA line) |

### PA-line grammar corner cases

PROTMOTIF-PATTERN-001 sharpens three grammar rules the converter must honour (ScanProsite
doc / PROSITE User Manual §IV.E):

- **Ranges only on `x`.** `x(2,4)` is valid; a range on a residue letter (`A(2,4)`) is **not**
  a valid PROSITE element, though a fixed count `A(3)` is.
- **Trailing period terminates the pattern** — characters after the `.` are ignored.
- **Reject the `*` Kleene star.** `<{C}*>` uses the ScanProsite *query* extension, not the
  standard PA-line grammar; the converter must raise `FormatException` rather than silently
  treat `*` as a residue ("reject, don't silently drop").
- **`>` inside brackets = "residue OR end-of-sequence".** PROTMOTIF-PROSITE-001 pins this rare
  case (only PS00267 and PS00539 use it): `[G>]` means Gly-or-C-terminus and converts to the
  regex alternation `(?:G|$)`, **not** a plain class. Verified: PS00267 `F-[IVFY]-G-[LM]-M-[G>]`
  → `F[IVFY]G[LM]M(?:G|$)`, PS00539 `F-[GSTV]-P-R-L-[G>]` → `F[GSTV]PRL(?:G|$)`. Matching must
  succeed via **both** the `G` branch and the C-terminus branch, and fail mid-sequence without
  a `G`.
- **Period terminates parsing (mid-pattern too).** Beyond a trailing `.`, a period anywhere
  ends the pattern — `R-G-D.A-B-C` → `RGD` (content after `.` is dropped), per User Manual §IV.E.

## Overlapping-match discovery

`FindMotifByPattern` wraps the pattern in a **zero-width lookahead** `(?=(pattern))` so that
the regex engine advances one position at a time and reports **all** matches — including
overlapping occurrences — mirroring ScanProsite's default behaviour (De Castro et al. 2006).
Plain (non-lookahead) regex matching consumes each match and would miss overlaps.

## API contract and invariants

| Aspect | Behaviour |
|--------|-----------|
| Coordinates | **0-based** `MotifMatch.Start`/`End` (repository convention; ScanProsite reports 1-based) |
| Case | matching is **case-insensitive** (protein sequences may be lower/upper; PROSITE codes are uppercase) |
| Empty / null input | returns no results |
| Invalid regex | handled gracefully (a malformed caller pattern does not throw an unhandled exception) |
| Overlaps | all occurrences reported, including overlapping ones (lookahead) |

## Information-content scoring

Scoring follows the sequence-logo information measure of Schneider & Stephens (1990) —
**not** BLAST statistics (an earlier BLAST citation was explicitly removed from the Evidence
because the E-value here is a direct combinatorial probability, not a Karlin–Altschul
extreme-value fit):

- **`CalculateMotifScore`** — information content `IC = Σ log₂(20 / allowed_count)` summed
  over each **constrained** position. An unconstrained `x` position contributes 0 bits; a
  fully fixed single-residue position contributes `log₂ 20 ≈ 4.32` bits; an `[ABC]`
  three-residue set contributes `log₂(20/3)` bits.
- **`CalculateEValue`** — `E = (N − L + 1) · 2^(−IC)`, where `N` = sequence length,
  `L` = motif length, `IC` = total information content. This is the expected number of random
  matches under a uniform amino-acid background.

## Verified reference patterns and two corrected bugs

The Evidence pins the engine's PROSITE-pattern catalog against the official database and, in
the process, **fixed two implementation patterns** that had diverged from PROSITE:

- **PS00007 (TYR_PHOSPHO_SITE_1)** — official `[RK]-x(2)-[DE]-x(3)-Y` → correct regex
  `[RK].{2}[DE].{3}Y`. The implementation had used a loosened `[RK].{2,3}[DE].{2,3}Y`
  (wrong repeat counts) — corrected to exact counts.
- **PS00018 (EF_HAND_1)** — official
  `D-{W}-[DNS]-{ILVFYW}-[DENSTG]-[DNQGHRK]-{GP}-[LIVMC]-[DENQSTAGC]-x(2)-[DE]-[LIVMFYW]`.
  The implementation had (1) used `x` (any) at position 2 instead of `{W}` (any except Trp)
  and (2) dropped the trailing `[LIVMFYW]` element — both corrected.

The full verified accession set (patterns confirmed against ExPASy PROSITE) includes PS00001,
PS00004–PS00009, PS00016–PS00018, PS00028 (C2H2 zinc finger) and PS00029 (leucine zipper).
Several (PS00001/00004–00009/00016/00017/00029) carry PROSITE's **SKIP-FLAG=TRUE** — very
common sites that generate many false positives, which callers may choose to skip.
PROTMOTIF-PROSITE-001 adds a **real-protein positive control**: Human Transferrin (UniProt
P02787, TRFE_HUMAN) scanned against PS00001 `N-{P}-[ST]-{P}` yields exactly 2 N-glycosylation
sites at ScanProsite 1-based 432–435 and 630–633 (0-based 431–434 and 629–632).

## Non-PROSITE literature-verified patterns

Beyond PROSITE, the `CommonMotifs` dictionary carries five short linear motifs, each verified
against its primary publication (not heuristics):

| ID | Pattern | Motif | Reference |
|----|---------|-------|-----------|
| NLS1 | `K-[KR]-x-[KR]` | monopartite nuclear localization signal | Dingwall & Laskey 1991 |
| NES1 | `[LIVFM]-x(2,3)-[LIVFM]-x(2,3)-[LIVFM]-x-[LIVFM]` | nuclear export signal | la Cour et al. 2004 |
| SIM1 | `[VIL]-x-[VIL]-[VIL]` | SUMO-interacting motif (type b) | Hecker et al. 2006 |
| WW1 | `P-P-x-Y` | WW-domain PPxY ligand | Chen & Sudol 1995 |
| SH3_1 | `[RK]-x(2)-P-x(2)-P` | SH3 Class I ligand (+xxPxxP) | Mayer 2001 |

## Deviations and assumptions

No algorithm deviations remain open — the Evidence's change history records that all
assumptions were eliminated: the two mispatterned PROSITE entries were corrected, the five
non-PROSITE patterns were re-derived from literature, heuristic scoring was replaced by
information-content scoring (Schneider & Stephens 1990), and overlapping-match discovery was
implemented via regex lookahead. The only standing API-shape convention is 0-based vs
ScanProsite's 1-based coordinates (no correctness effect; shared with the sibling ProteinMotif
units).

## References

ExPASy PROSITE database and User Manual (accessed 2026-02-12); Hulo et al. 2007 (*Nucleic
Acids Res.* 36:D245–9, "20 years of PROSITE"); De Castro et al. 2006 (ScanProsite);
Schneider & Stephens 1990 (*Nucleic Acids Res.* 18:6097–6100, sequence logos / information
content); Dingwall & Laskey 1991; la Cour et al. 2004; Hecker et al. 2006; Chen & Sudol 1995;
Mayer 2001. Full citations in [[protmotif-find-001-evidence]],
[[protmotif-pattern-001-evidence]] and [[protmotif-prosite-001-evidence]] (do not duplicate
here).
