# Genomic analysis: one argument contract (GenomicAnalyzer)

> Audit R14-22 (B09 F295). The methods of the genomic-analysis unit used to treat the same kind of argument differently (`FindRepeats` / `FindCommonRegions` took `minLength <= 0` as 1, `FindTandemRepeats` threw for `minUnitLength = 0`, `FindOpenReadingFrames` accepted 0 and -5 silently, the MCP tools refused `minLength < 1` in two tools and tolerated it in the third). The contract below is what every method does now; `GenomicAnalyzer_ArgumentContract_Tests` and the MCP tests lock it.

## The two kinds of numeric argument

| Kind | Examples | Library | MCP tool |
|------|----------|---------|----------|
| **A minimum LENGTH of a reported substring / region / ORF** (`minLength`) | `FindRepeats`, `FindMaximalRepeats`, `FindRepeatsBothStrands`, `FindCommonRegions` (every mode), `FindCommonRegionSpans`, `FindOpenReadingFrames` | a value **<= 0 means "no lower bound"**: accepted, and the result is exactly that of `minLength = 1` (a repeat or a region cannot be the empty string; every ORF is at least one codon, so 1 and 3 admit the same ORFs). No exception for any `int`. | **refused**: `ArgumentOutOfRangeException` ("Minimum repeat length / common-region length / ORF length must be at least 1") for `minLength < 1` in `find_repeats`, `find_common_regions`, `find_open_reading_frames` (the last one also with `orfFinderMinLength`, which applies the NCBI ORFfinder rule instead: below 30 is changed to 30) |
| **A STRUCTURAL parameter** (a count or a size that defines the search) | `FindTandemRepeats(minUnitLength, minRepetitions, maxUnitLength)`: `minUnitLength >= 1`, `minRepetitions >= 2`, `maxUnitLength >= minUnitLength` | **throws** `ArgumentOutOfRangeException` naming the parameter (a zero-length unit never terminates the scan; one copy is not a repeat) | the same refusal (`find_tandem_repeats`, `find_microsatellites`, ...) |

Other rules of the unit that are NOT "minimum lengths" and stay as documented on their methods: `FindPalindromes` (`minLength` even and >= 4, `maxLength >= minLength`), the output limit of the repeat finders (`RepeatOutputLimits`, F292), the per-program ranges of the EMBOSS / MUMmer / T-Coffee ports (the program's own range or refusal).

## Why "accept and clamp" for the library and "refuse" for the tools

A library caller computes a length from data (a fraction of a read, a user slider) and "no minimum" is a legitimate answer; the result is well defined (every repeat / region / ORF). A tool call comes from a client that may not know the convention, and a `0` or `-5` there is far more often a mistake than a request for "everything": the tools therefore fail fast with the same message shape in the three tools. Where the library cannot give a meaningful answer (a zero-length tandem unit, one repetition) it throws in both.

## Verification

- `GenomicAnalyzer_ArgumentContract_Tests`: for `0, -1, -5, int.MinValue` the repeat finders (DnaSequence, string, maximal, both strands), the common-region finders (DnaSequence, string, spans, MEM mode) and `FindOpenReadingFrames` (DnaSequence, string) return exactly the `minLength = 1` result; `FindTandemRepeats` throws for `minUnitLength <= 0`, `minRepetitions < 2` and `maxUnitLength < minUnitLength`.
- MCP `FindCommonRegionsTests.MinLength_BelowOne_IsRefused_LikeFindRepeatsAndFindOpenReadingFrames`: the three tools refuse 0, -1, -5.
