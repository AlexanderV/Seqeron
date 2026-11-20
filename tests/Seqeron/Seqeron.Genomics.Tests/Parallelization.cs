// Tier 1: Enable parallel test execution for stateless test fixtures.
// All static algorithm methods (MotifFinder, DisorderPredictor, SequenceAligner, etc.)
// have no shared mutable state, so tests can safely run in parallel.
// I/O-bound fixtures (FastaParser, etc.) should be marked [NonParallelizable].
//
// Memory (audit A3-15, B09 F144): "no shared mutable state" is not "no shared memory": eight workers run fixtures side by side, and a few
// fixtures allocate hundreds of MB to GB (5 600 x 5 600 affine matrices, suffix trees of 3 M bases, mega-character oracle inputs). The peak working set of the
// heavy Alignment / Analysis fixtures is measured in docs/Validation/review-2026-09/LOCAL_RUN.md (section 5); every fixture above ~0.5 GB is [NonParallelizable]
// (SequenceAligner_Cancellation_Tests 2.5 GB, SequenceAligner_LinearSpace_Tests 2.2 GB, GenomicAnalyzer_FindOpenReadingFrames_Oracle_Tests, SequenceAligner_TopNLinearSpace_Tests,
// GenomicAnalyzer_FindCommonRegion_Tests). A NEW fixture that allocates more than that must carry the attribute too.
[assembly: Parallelizable(ParallelScope.Children)]
[assembly: LevelOfParallelism(8)]
