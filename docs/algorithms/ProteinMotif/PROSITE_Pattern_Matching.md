# PROSITE Pattern Matching

## Documented Theory

### Purpose

PROSITE pattern matching converts protein motif patterns expressed in the PROSITE
notation (Hulo et al. 2007) into regular expressions, then scans protein sequences
for all occurrences of the resulting pattern. This enables identification of
biologically functional sites (post-translational modifications, binding domains,
structural motifs) using the formal PROSITE syntax standard.

### Core Mechanism

The PROSITE pattern notation is formally defined in the PROSITE User Manual
(https://prosite.expasy.org/prosuser.html, Section IV.E — The PA line):

| Syntax Element | Meaning | Regex Equivalent |
|---------------|---------|------------------|
| Uppercase letter (e.g., `A`) | Specific amino acid | Same letter |
| `x` | Any amino acid | `.` |
| `[ABC]` | Any of listed amino acids | `[ABC]` |
| `{ABC}` | Any amino acid except listed | `[^ABC]` |
| `-` | Element separator | (dropped) |
| `(n)` | Repeat preceding element n times | `{n}` |
| `(n,m)` | Repeat `x` between n and m times | `{n,m}` |
| `<` | N-terminus anchor | `^` |
| `>` | C-terminus anchor | `$` |
| `.` (period) | Pattern terminator (in data files) | (dropped) |

**Constraints (from ScanProsite documentation):**
- Range `(n,m)` is only valid with `x`, not with specific amino acids.
- Fixed repetition `(n)` is valid with any element including character classes.
- In rare cases (PS00267, PS00539), `>` can appear inside brackets: `[G>]` means
  G or end-of-sequence.

After conversion, the regex is used with lookahead-based matching to detect
overlapping occurrences, consistent with ScanProsite's overlap-enabled default
mode (De Castro et al. 2006).

### Properties

- **Deterministic:** Same input always produces same output.
- **Case-insensitive:** Sequences are normalized to uppercase before matching.
- **Complete:** All occurrences (including overlapping) are returned via regex lookahead.

### Complexity

| Aspect | Value | Source |
|--------|-------|--------|
| Time | O(n × m) where n = sequence length, m = pattern length | Regex matching with NFA |
| Space | O(k) where k = number of matches | Match storage |

---

## Implementation Notes

**Implementation location:** [ProteinMotifFinder.cs](src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs)

- `ConvertPrositeToRegex(prositePattern)`: Converts PROSITE pattern notation to
  .NET-compatible regular expression string. Handles all standard syntax elements.
- `FindMotifByProsite(proteinSequence, prositePattern, motifName)`: End-to-end
  PROSITE matching — converts pattern via `ConvertPrositeToRegex`, then delegates
  to `FindMotifByPattern` for regex-based scanning with overlapping match detection.

---

## Deviations and Assumptions

None. All PROSITE syntax elements defined in the PROSITE User Manual §IV.E are fully
implemented, including `[G>]` C-terminus inside brackets (PS00267, PS00539).

---

## Sources

- PROSITE User Manual. SIB Swiss Institute of Bioinformatics. https://prosite.expasy.org/prosuser.html (Accessed 2026-02-12)
- ScanProsite Documentation. https://prosite.expasy.org/scanprosite/scanprosite_doc.html (Accessed 2026-02-12)
- Hulo N, Bairoch A, Bulliard V, et al. (2007). "The 20 years of PROSITE." Nucleic Acids Res. 36(Database issue):D245-9. https://doi.org/10.1093/nar/gkm977
- De Castro E, Sigrist CJA, Gattiker A, et al. (2006). "ScanProsite." Nucleic Acids Res. 34:W362-365. https://doi.org/10.1093/nar/gkl124
