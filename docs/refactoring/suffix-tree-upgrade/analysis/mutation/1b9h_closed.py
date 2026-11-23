# CLOSED = {mutant id: ('test' | 'equivalent', text)}: how every survivor of the audit of work package 1b.9h was closed (read by summarize_1b9h.py). 'test': the test that kills it at the head
# that holds it (the survivors were run again, results_1b9h_closed_*.tsv); 'equivalent': the argument.
CLOSED = {
    # the registry of the mapped files (checked build)
    '1b9h_G01': ('test', 'ARelativeNameAndTheAbsoluteNameOfAFile_AreTheSameFile'),
    '1b9h_G08': ('test', 'ARelativeNameAndTheAbsoluteNameOfAFile_AreTheSameFile'),
    '1b9h_G07': ('test', 'ANameThatWasNeverListed_IsTakenOffWithoutAnError_AndLeavesTheListAsItWas'),
    '1b9h_G13': ('test', 'TwoMappingsOfOneFile_AreOneFileInTheList'),
    '1b9h_G14': ('test', 'TwoNamesThatDifferInCaseOnly_AreTwoFiles_WhereTheFileSystemSaysSo'),
    # the bounds checks (checked build)
    '1b9h_T06': ('test', 'TheGuardOfTheRawReads_RefusesANegativeOffset_AndACountThatLeavesTheRegion (every caller refuses a negative index first, so the guard is reached only by a call of its own)'),
    '1b9h_T16': ('test', 'ABulkRead_ThatWouldLeaveTheRegion_IsAMemoryViolation_AtTheExactByte now with 5 and 37 symbols (the bytes of a range are rounded up)'),
    '1b9h_T22': ('test', 'ADisposeThatIsRepeated_TakesOffOnlyItsOwnMapping (two readers of one file)'),
    '1b9h_R21': ('test', 'ADisposeThatIsRepeated_TakesOffOnlyItsOwnMapping (two tables of one file)'),
    '1b9h_R07': ('equivalent', 'the slot is block * classes + index with block = p / period <= n / period for p <= n (checked first) and 0 <= index < classes, so slot <= (n / period) * classes + classes - 1 = Slots - 1: the check of the slot cannot fire once the position and the residue are checked (defence in depth for the memory-safety condition itself)'),
    '1b9h_R15': ('equivalent', 'ReadSlot is called by Verify for slot = 0 .. Slots - 1 only; the check of the slot cannot fire'),
    '1b9h_R16': ('equivalent', 'Verify refuses a disposed table at its entry (the test of R18); a Dispose from another thread during a Verify is a race that the contract excludes'),
    '1b9h_R17': ('equivalent', 'the checked reads of Verify are ReadSlot (slot in [0, Slots), not disposed), and the loop and the entry of Verify give both; the raw read of a slot in range is the same read'),
    # the allocation meter
    '1b9h_A04': ('test', 'AReadingAboveTheLimit_WithACollection_IsInconclusive_NotAFailure (the first run of the audit had no result: a host that the kernel killed for lack of memory)'),
    '1b9h_A13': ('test', 'ARegionThatTheCodeEnds_IsReported_AlthoughNoCollectionRan'),
    '1b9h_A14': ('test', 'ACollectionThatIsFollowedByANewRegion_IsReported_AlthoughTheModeIsTheRegionAgain'),
    '1b9h_A09': ('equivalent', 'the counter read before the start of the region against after: the start is a collection made by this thread, which retires the allocation context and subtracts its unused rest (the probe, helper 0, window l: 2,429 collections made by the thread inside windows, 0 readings off), so the reading is the same'),
    '1b9h_A12': ('equivalent', 'without a region opened by the meter or by the caller the mode is not NoGCRegion at the end (a region opened by another thread inside the window is a race that the use excludes)'),
    '1b9h_A17': ('equivalent', 'the runtime does not end the region for large objects beyond the reserve of the large-object heap while the total budget holds (8 MB against a reserve of 1 KB: the region is intact); the reserve is the documented intent, not an observable'),
    '1b9h_A18': ('equivalent', 'reached only if another thread opens a region between the check of the mode and the start (a race): the catch makes the window measured without the promise'),
    '1b9h_A19': ('equivalent', 'EndNoGCRegion throws when the region was left, and the catch of InvalidOperationException takes it: the check of the mode only avoids the exception'),
    '1b9h_A20': ('equivalent', 'the check of the mode already avoids the call when the region was left; the exception is reached only in a race between the check and the call'),
    # the discovery limit
    '1b9h_D04': ('test', 'the tests of the limit now have [HangLimit(30)]: the mutant no longer stops the host, it fails AComputationThatDoesNotReturn_... at its limit'),
    '1b9h_D09': ('test', 'the same'),
    '1b9h_D07': ('equivalent', 'at the limit the computation has not returned, so its result is still the default value'),
    # the hang breaker
    '1b9h_H12': ('test', 'TheAttribute_FailsATestThatBeginsWhileAnotherIsStillInside_AndNamesIt (the attribute takes a breaker of its own for the tests, so that the breaker of the host is not opened by them)'),
    '1b9h_H04': ('equivalent', 'after the first runaway is recorded every Begin returns it whatever the mark inside says, and the mark is read only to find a runaway when there is none: a refused test that marks itself is not observable through Begin, End or Runaway'),
    '1b9h_H11': ('equivalent', 'the constructor of the attribute computes the same product with checked (the mutants H13 and H14 are on that expression) and throws first, so Seconds * 1000 in ApplyToContext cannot overflow'),
    '1b9h_H03': ('equivalent', 'once the first runaway is recorded no test is marked as inside (a refused test does not mark itself), so the mark stays that of the runaway or null and an exchange writes the same name'),
    # the difference cover table
    '1b9h_C15': ('test', 'TheTableIsCheckedWhenItIsUsed_AndADamagedTableIsRefused now with the first two residues swapped and the whole set reversed'),
    '1b9h_C01': ('equivalent', 'for a period other than 993 the table (32 residues) fails the length (q + 1) or the check of the cover, FromTable returns null and the search runs: the cover is the same, and so is the number of searches'),
    '1b9h_C08': ('equivalent', 'a set of residues in Z_X in which every non-zero difference occurs exactly once (IsPerfect) has k (k - 1) = X - 1 = q^2 + q ordered differences, so k = q + 1: the length is implied by the perfect-set check'),
    '1b9h_C09': ('equivalent', 'the same count: a perfect set of 32 residues exists only for X = 993, so a period that is not 993 never passes IsPerfect with 32 residues'),
    '1b9h_C13': ('equivalent', 'a repeated residue leaves at most q distinct residues, whose q (q - 1) ordered differences cannot cover the X - 1 = q^2 + q non-zero values: Verify and IsPerfect refuse it'),
    '1b9h_C16': ('equivalent', 'with q + 1 distinct residues in range, "every difference occurs" (Verify) and "every non-zero difference occurs exactly once" (IsPerfect) are the same statement: (q + 1) q = X - 1 ordered differences cover X - 1 values iff each occurs once'),
    '1b9h_C17': ('equivalent', 'the same: Verify is implied by IsPerfect, and IsPerfect by Verify, for q + 1 distinct residues'),
    '1b9h_C18': ('equivalent', 'the same'),
}
