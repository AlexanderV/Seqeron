# Mutants of the audit of 1b.7a that are argued equivalent (id -> argument). Read by make_1b7a_round2.py (they are reported, not re-run) and quoted in the ADR subsection.
EQUIVALENT = {
    'A18': 'BudgetedArena.Dispose without marking the arena disposed (replacement of the first A18, which did not compile): every later call reaches NativeArena, which throws the same ObjectDisposedException, and a second Dispose releases 0 bytes because the live total was zeroed',
    'A19': 'Dispose keeps the live total: a second Dispose returns at the _disposed guard before it reads the total, so the total is released once either way',
    'A20': 'Mark without its own disposal check: the inner NativeArena.Mark throws the same ObjectDisposedException',
    'A23': 'ReleaseTo without the upper-bound check of the mark: the loop does not run for mark > count, freed is 0, and NativeArena.ReleaseTo throws the same ArgumentOutOfRangeException before the finally releases 0 bytes',
    'C04': 'SafeRead maps an empty read to -1 instead of 0: both callers (ProcIoSnapshot.TryRead, RssReader.Read) treat got <= 0 as unavailable',
    'C06': 'the key at the very end of the text is not searched (< instead of <=): a key with nothing after it never parses a number, the result is false either way',
    'E06': 'CountedFile.Flush (fsync) does nothing: the only difference is durability after a crash or the block-layer write time, which no assertion of a test on a correct file system can observe',
    'F05': 'the range test without offset < 0: a negative offset reaches OpenWindow, whose CreateViewAccessor throws the same ArgumentOutOfRangeException',
    'F14': 'SetLength(Math.Max(length, fs.Length)): FileMode.Create has truncated the file already, so fs.Length is 0 and the maximum is the length',
    'G18': 'the block name "dist.tuples" appears only in the message of a checked-build violation',
    'H10': 'Phi accepts an entry equal to n after the first: the key range test of Distribute (key >= slots with slots = n) refuses it with the same InvalidDataException',
    'H35': 'Dist without the final VerifyCanaries: on a run whose memory accesses are correct the guard bands are intact, the call is defence in depth that turns a memory bug into an exception (all Dist tests run with guarded arenas)',
    'I28': 'SpillStack.Dispose without the _disposed guard: a second Dispose finds _file null and returns',
    'I04': 'a pop from the last entry of a block retreats only when height > 2 instead of > 1: with blocks of at least two entries a top block that holds one entry and has a block below it means a height of at least 3',
}
