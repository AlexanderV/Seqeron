#!/usr/bin/env python3
"""Work package 1b.9e: the mutants that survived every test of the default set (and the tests written for the survivors) and are argued equivalent, one by one. Each argument was also checked by the whole default
set of SuffixTree.Esa.Tests (stage 2 of run_1b9e_audit.py): an equivalent mutant must survive it. summarize_1b9e.py reads EQUIVALENT."""
EQUIVALENT = {
    'O07': ('MappedOutputFile.OpenExisting: a negative cut is not refused by the explicit check',
            'The call that follows, RandomAccess.SetLength(handle, cut), refuses a negative length with the same exception type (ArgumentOutOfRangeException), and GetLength(handle) > cut is then true for every file. '
            'The test asserts the type, and the mutant gives it; only the text of the message differs. The explicit check is kept for its clearer message.'),
    'L04': ('SesaBuildLock.Dispose: Interlocked.Exchange replaced by a plain read',
            'The exchange makes a second Dispose a no-op even from two threads at once. With a plain read a second call (same thread) finds the stream still set, deletes the already deleted file (no exception: File.Delete of a missing file is '
            'silent) and disposes the stream again (no exception). The difference needs two threads disposing one lock at the same instant, which nothing in the build does (the lock is a using variable of one method); not observable by a deterministic test.'),
    'L06': ('SesaBuildLock.Acquire: FileMode.Create instead of OpenOrCreate',
            'The lock is the exclusive handle, not the content of the file: nobody reads the lock file, so whether Create truncates it, and whether the truncation comes before or after the lock is taken, changes nothing that a build or a test can see. '
            'Both modes refuse a second holder with the same IOException (FileShare.None) and both take a stale file; the tests that take the lock twice and that find a stale file pass with either.'),
    'B36': ('SesaBuilder.Build: the scratch directory is not removed by the finally of a published build',
            'The body of the build removes the scratch directory itself after the rename (DeleteDirectoryQuietly(scratch) three lines after File.Move), so for a published build the finally has nothing left to remove; the mutant only drops a second, empty-handed call.'),
    'B55': ('the adapter of the checkpoint reports zero seconds of recording',
            'The seconds that the recording took are subtracted from the seconds that PH1 and PH2 / PH3 report (SaCheckpoint tests: a checkpoint that reports 1000 s and 2000 s moves the two figures by exactly that, mutants S08 and S09 are killed by them). '
            'The adapter hands the real figure; with zero the two phase times of a resumed build include the 1-10 milliseconds of the recording. No output, no counter and no decision depends on a phase time.'),
    'B56': ('the adapter of the checkpoint hands no cover (the build makes its own)',
            'Both ways give the same cover (the cover is a function of the period and the exactness); the checkpoint only saves making it twice, which costs 0.27 s per build at X = 993 (the marker `created` of `sesa resume build` comes 0.28 s after `scanned` with the cover made twice and 0.003 s once; '
            'the resume build was 12-17 % slower at 10^7 before the change). SaCheckpoint1b9eTests.TheCoverThatTheCheckpointHolds kills the same change in PH1-PH3 itself (mutant S11), so the contract is tested where it is used; the value that the adapter hands cannot change a byte.'),
    'C31': ('SesaResumeRun.Cover returns null',
            'The same as B56 from the other side: the run still computes the cover (and uses it to check the rank file); the build then makes its own for PH1. Same value, 0.27 s more.'),
    'C11': ('Commit: the new stage is added to all the stages of the state, not to the first _kept',
            'At every Commit the state holds exactly _kept stages: Load sets _kept to the number of verified stages and Prepare (always called for a resumed build, tested by B09 and V70/V71) cuts the state to them before anything runs, and a fresh build starts with none; '
            'Commit sets _kept to the new count. Truncate(_kept) is a no-op there. It is kept as a guard against a future caller.'),
    'C13': ('WriteState: FileStream.Flush(false) instead of Flush(Fsync)',
            'Whether the state file reaches the device before it replaces the old one is durability after a power loss; the process, the page cache, kill -9, the counters (the flush is counted separately, tested) and every byte that a later run reads are the same. '
            'Nothing in a test or in the kill campaign can lose power. The cost of the real flush is in the cost tables (two fsyncs per stage).'),
    'C19': ('PrimaryRow: a kept BWT record is looked up before the PSI record',
            'A state never holds both: SesaResumeState.With replaces the BWT stage when the PSI stage is recorded and the parser refuses a state that has both in the wrong order; and the two records carry the same primary row. The order of the lookup cannot matter.'),
    'C26': ('CoverDigest: the digest does not include the period',
            'The period is recorded and compared on its own (V23, killed), and a cover is a function of its period; two covers that differ only in the period have different residues, so the digest of the residues alone still differs. The digest of the period adds nothing.'),
    'V04': ('HeaderMismatch: another sigma is accepted',
            'sigma = the number of letters, and the state is refused unless its list of letters equals the letters of this input (V06, killed) and the parser refuses a state whose letters are not as many as its sigma (F16, killed): a state with another sigma has other letters.'),
    'V52': ('Expected: an unknown stage name is accepted',
            'The parser refuses a state with a stage name that is not one of the six (IndexOf = -1 is never above the previous stage): the default branch is unreachable from a parsed state. It stays as the safe answer to a state that is built in memory.'),
    'V72': ('Prepare: a build whose state was not accepted also prepares',
            'The caller calls Prepare only for a resumed build (if (resumed) run.Prepare()); the guard at the top of Prepare is a second line of defence and the mutant removes it only for a call that never happens.'),
}
