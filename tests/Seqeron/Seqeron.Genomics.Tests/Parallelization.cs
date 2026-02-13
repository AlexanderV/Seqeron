using NUnit.Framework;

// Tier 1: Enable parallel test execution for stateless test fixtures.
// All static algorithm methods (MotifFinder, DisorderPredictor, SequenceAligner, etc.)
// have no shared mutable state, so tests can safely run in parallel.
// I/O-bound fixtures (FastaParser, etc.) should be marked [NonParallelizable].
[assembly: Parallelizable(ParallelScope.Children)]
[assembly: LevelOfParallelism(8)]
