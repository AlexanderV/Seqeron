# Protein Motif Search

## Documented Theory

### Purpose

Protein motif searching identifies occurrences of short, conserved sequence patterns
within protein sequences. These patterns correspond to biologically functional sites
such as post-translational modification sites, binding domains, and structural motifs.
The PROSITE database (Hulo et al., 2007) is the authoritative collection of such
patterns, expressed in a formal pattern notation that can be converted to regular
expressions for computational scanning.

### Core Mechanism

The algorithm converts PROSITE-format patterns to regular expressions and uses
regex matching to find all non-overlapping occurrences of each pattern in a
protein sequence. The PROSITE pattern syntax is formally defined as follows
(PROSITE User Manual, https://prosite.expasy.org/prosuser.html):

- Single amino acid letters match themselves (IUPAC one-letter codes)
- `x` matches any amino acid
- `[ABC]` matches A, B, or C (alternative residues)
- `{ABC}` matches any amino acid except A, B, or C (exclusion)
- `x(n)` matches exactly n copies of any amino acid
- `x(n,m)` matches between n and m copies of any amino acid
- `<` anchors to the N-terminus; `>` anchors to the C-terminus
- Elements are separated by `-`

The conversion to regex is:
- `x` → `.`
- `[ABC]` → `[ABC]`
- `{ABC}` → `[^ABC]`
- `x(n)` → `.{n}`
- `x(n,m)` → `.{n,m}`
- `<` → `^`
- `>` → `$`
- `-` separators are dropped

### Properties

- **Deterministic:** The same sequence and pattern always yield the same matches.
- **Complete for non-overlapping matches:** All non-overlapping regex matches are returned.
- **Case-insensitive:** Input sequences are normalized to uppercase before matching.

### Complexity

| Aspect | Value | Source |
|--------|-------|--------|
| Time | O(n × m) where n = sequence length, m = pattern count | String matching with finite automata |
| Space | O(n) for match storage | Standard regex matching |

---

## Implementation Notes

**Implementation location:** [ProteinMotifFinder.cs](src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs)

- `FindMotifByPattern(proteinSequence, regexPattern, motifName, patternId)`: Scans a protein sequence for all non-overlapping matches of a regex pattern. Returns `IEnumerable<MotifMatch>` with position (0-based Start/End), matched sequence, name, score, and E-value.
- `FindCommonMotifs(proteinSequence)`: Scans a protein sequence against all entries in the `CommonMotifs` dictionary. Delegates to `FindMotifByPattern` for each motif.
- `FindMotifByProsite(proteinSequence, prositePattern, motifName)`: Converts a PROSITE-format pattern to regex via `ConvertPrositeToRegex`, then delegates to `FindMotifByPattern`. (Scope: PROTMOTIF-PROSITE-001)
- `ConvertPrositeToRegex(prositePattern)`: Converts PROSITE pattern notation to .NET regex. (Scope: PROTMOTIF-PROSITE-001)
- `CommonMotifs`: Read-only dictionary mapping PROSITE accession numbers to `PrositePattern` records containing both PROSITE and regex representations.

---

## Deviations and Assumptions

| # | Item | Type | Status | Resolution |
|---|------|------|--------|------------|
| 1 | **PS00007 wrong repeat ranges** | Deviation | 🔧 FIXED | Implementation had `[RK]-x(2,3)-[DE]-x(2,3)-Y`, PROSITE defines `[RK]-x(2)-[DE]-x(3)-Y`. Fixed pattern and regex. |
| 2 | **PS00018 wrong position 2 and missing trailing element** | Deviation | 🔧 FIXED | Implementation had `D-x-[DNS]-...` and omitted `[LIVMFYW]`. PROSITE defines `D-{W}-[DNS]-...-[DE]-[LIVMFYW]`. Fixed both pattern and regex. |
| 3 | **Non-PROSITE motifs in CommonMotifs** | Assumption | ⚠ ASSUMPTION | Entries NLS1, NES1, SIM1, WW1, SH3_1 are not from PROSITE. Their patterns are based on literature but not PROSITE-verified. Documented. |
| 4 | **Scoring heuristic** | Assumption | ⚠ ASSUMPTION | `CalculateMotifScore` and `CalculateEValue` are implementation-specific heuristics not from any authoritative source. Scores are not tested for exact values. |
| 5 | **Non-overlapping matching** | Assumption | ⚠ ASSUMPTION | .NET Regex.Matches returns non-overlapping matches. PROSITE scanning tools may differ. This is a known scope limitation. |

---

## Sources

- PROSITE Database. SIB Swiss Institute of Bioinformatics. https://prosite.expasy.org/ (Accessed 2026-02-12)
- PROSITE User Manual. https://prosite.expasy.org/prosuser.html (Accessed 2026-02-12)
- Hulo N, Bairoch A, Bulliard V, et al. (2007). "The 20 years of PROSITE." Nucleic Acids Res. 36(Database issue):D245-9. https://doi.org/10.1093/nar/gkm977
- De Castro E, Sigrist CJA, Gattiker A, et al. (2006). "ScanProsite." Nucleic Acids Res. 34:W362-365. https://doi.org/10.1093/nar/gkl124
- Wikipedia. "Sequence motif." https://en.wikipedia.org/wiki/Sequence_motif (Accessed 2026-02-12)
