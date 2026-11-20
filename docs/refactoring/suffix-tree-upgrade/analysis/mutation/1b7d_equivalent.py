# Mutants of the audit of 1b.7d that no test kills, in two groups (id -> argument). Read by summarize_1b7d.py and quoted in the ADR subsection.
# EQUIVALENT: the mutant behaves like the original in every state that a caller or a test can reach (the argument says why); DEFENCE: a second line of defence behind a check that is tested at another level, or a
# failure path that no injectable fault reaches (an exception in the middle of a scan, after the creation of a file inside a call): not equivalent, not killed, reported as such.
EQUIVALENT = {
    'B13': 'InMemory differs only on a disposed run (both fields are null after Dispose); no caller asks a disposed run',
    'D03': 'Verify is used by tests only and differs only for a residue at or above the period, which Create never produces (the constructor is private)',
    'D05': 'IsPerfect is used by tests only; the counter of the difference 0 is not read (the check loops from d = 1)',
    'D21': 'IsPrime(1) is never called: both callers loop from q = 2',
    'L01': 'the recursion ends by the shrink rule (a reduced string that is not shorter than its text goes to the base case) long before 64 levels; a limit of 6,400 is never reached',
    'M05': 'the limit of the common prefix of two suffixes is clamped inside MappedTextReader.Lce, which every caller uses (the same argument as the PlcpScan limit mutants of 1b.7c)',
    'M09': 'a leaf of no records contributes an empty run to the merge, which adds nothing; the splitter ranges hold at least one record',
    'M15': 'no worker thread is started: the calling thread does all the work; the output and the shape do not depend on the thread count (tested), only the time does',
    'M21': 'merging an empty sample: the merge of zero records writes nothing and returns, as the early return does',
    'M31': 'a leaf with one group of names handled as an inner node of one range assigns the same name to the same records; the output and the counted work are the same',
    'M51': 'a reservoir of two is refused at the entry of BuildRankFile (ReservoirSize >= 3, tested); the check in Split repeats it',
    'M78': 'the order of two equal elements in the merge: two different positions never compare equal under the total order of the suffixes',
    'F03': 'a bucket of 25 sorted by insertion instead of by radix: both give the sorted order of the keys, equal keys keep their order in both',
    'F11': 'the binary search over splitters with a tie: a splitter and a position are different suffixes and never compare equal',
    'F13': 'the merge takes the right element on a tie: the same argument (no two elements compare equal)',
    'H21': 'no extra thread: the calling thread does all the work; the output and the shape do not depend on the thread count (tested), only the time does',
    'H42': 'the fan-out is also bounded by 2 m / K + 1 <= m for K >= 2 (validated), so the bound m is never the smaller one',
    'H54': 'the degenerate split cannot happen: the splitters are suffixes of the sample, so no child gets all m positions once the sample is the whole bucket; the retry is dead code (no test reaches even one retry)',
    'H56': 'the retry counter of the dead retry (see H54)',
    'H68': 'the number of attempts of the dead retry (see H54)',
    'H72': 'a missing directory at the cleanup: the exception is a DirectoryNotFoundException, an IOException like the one that is being rethrown from the failure',
    'Q43': 'a negative arena cap: the constructor of ArenaBudget, called before any file is created, refuses it with the same exception type',
    'Q99': 'fsync of the suffix array is not observable in the bytes or in the counters (the option is passed through to RandomAccess.FlushToDisk)',
    'Q101': 'the length of the suffix array file: it is created at exactly 4 (n + 1) bytes and never resized, so the check cannot fail',
    'Q121': 'the scratch subdirectories of the phases are created by the phases themselves; the creation in the builder is redundant',
    'A10': 'a chunk of zero slots is refused at the entry of BuildRankFile (ChunkSlots >= 1, tested); Math.Max repeats it',
    'X07': 'the width of a reduced reader: every caller passes ReducedBits (32); no path gives another width',
}
DEFENCE = {
    'E13': 'a text above twice the limit at the entry of BuildRankFile: the same limit is refused at the entry of SaBuilder (tested); a direct call with such a text needs a file of 2 GB and a run of PH1 on it',
    'H07': 'the same limit at the entry of BucketSorter.Run, behind the one of SaBuilder',
    'L15': 'the tag of the files of a nested sort: the nested sorter removes its own files when it fails (tested in BucketSorterTests); the tag lets the cleanup of the level find them, which only matters for a failure that the nested sorter does not handle',
    'L17': 'Complete of the suffix-array file of a level: it checks that no writer of the file is open, which the sorter guarantees (every writer disposed: tested)',
    'Q97': 'Complete of the suffix-array file: the same check, no writer is open (BucketSorter disposes all of them: tested)',
    'M19': 'the check that the scatter placed every record of the sample: it fires only on a defect of the scatter, which no injectable fault produces',
    'L29': 'the flag that the rank file of level 0 exists, which the cleanup uses after a failure inside the inversion (after the creation of the file): no injectable fault fails inside the inversion after the creation',
    'A08': 'the same flag in the base case: nothing fails between the creation of the file and the end of the base case',
    'H64': 'disposal of the writers of a partition when the scan fails in the middle (an I/O error of a spill while other writers hold files): the injectable spill failures fail at the creation of the first file of every writer',
    'H66': 'disposal of the children when the hand-over to the pool fails: nothing between the split and the push can fail',
}

# Work package 1b.9c (2026-11-20): the mutants of the lists that the audit did not run, and the ten "second line of defence" mutants above, closed by an argued equivalence (a reviewer can check each argument against the
# code at the head named in it). The mutants of the second group that a new test kills are in 1b9c_closed.py; summarize_1b9c.py reads both. Every one of these survived the whole default set of SuffixTree.Esa.Tests (stage 2 of the runner).
EQUIVALENT_1B9C = {
    'Q18': 'Resolve: every reference that passes through Scale when the run is not scaled (2^20, 16 x 2^20, 4 x 2^20, 2^20 and 64 x 2^10) is above its floor (3, 4, 16, 3, 0), and the factor is exactly 1.0, so Math.Max(floor, Round(reference x 1.0)) is the reference',
    'Q20': 'Resolve: every floor given to Scale is repeated by the Math.Max of the same value at the call (3, 4, 16, 3) or is 0 for the run spill, which is not negative (the factor is not), so dropping the inner floor changes nothing',
    'M19': 'SampleNamer.Merge: the check cannot fire: the loop before it throws unless every window received exactly its own number of records, and Scattered of a reader grows by the number that Drain returns (BucketRun.cs: both are i - _pos), so the sum over the readers is the sum of the windows, the number of samples',
    'F12': 'LeafSorter.MergeSort: the merge of two sorted halves is skipped when the halves are in order; CompareTied of two different positions is never 0 (different suffixes of one text are different), so "<= 0" and "< 0" agree (the same argument as F03, F11 and F13)',
    'L17': 'ExternalDcs.BuildLevel: Complete(fsync: false) of the file of the suffix array of a level only asserts that no writer of the file is open and that it is not disposed (it flushes nothing), and the sorter has disposed every writer (using, finally) when Run returned; the assertion itself is tested in isolation (MappedOutputFileTests.Complete_NeedsEveryWriterToBeDisposed_AndMayBeRepeated_AndSyncsOrNot)',
    'Q97': 'SaBuilder.BuildCore: the same as L17 for the file of the suffix array of the build (Complete(options.SyncOutput) with the option false is that assertion; with it true it is the fsync, which is not observable, see Q99)',
    'H33': 'BucketSorter.RunTasks: Partition returns the number of degenerate retries, and the retry is dead code (H54: the splitters are suffixes of the bucket, so no child gets all the positions), so the value is 0 on every path and `retries += 0` is the statement without it',
    'H41': 'BucketSorter.RunTasks: VerifyCanaries of the worker arena checks the guard bands of a test arena (a production arena has none), which are damaged only if the stage writes outside a buffer, and the stage does not; the check is tested in isolation (Safety/BufAndArenaTests: a damaged canary is reported, a production arena has nothing to verify) and runs at the end of every other stage that has an arena',
    'L16': 'ExternalDcs.BuildLevel: the nested sort without the ledger of the level makes its own ScratchLedger over the same IoCounters (BucketSorter.Run: options.Ledger ?? new ScratchLedger(options.Counters)), so the live and the peak scratch bytes of the counters are the same, and both ledgers give their bytes back at a failure (the sorter releases its own, the level releases the rest); the ledger itself is not visible outside the build',
}
