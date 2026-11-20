# Work package 1b.9c: the fresh mutants (F01-F62, make_1b9c_fresh.py) that the new tests did not kill when the audit first ran them (pass 1: the new test classes alone), and what closed each: a test (red against the mutant, green on the head: run again
# against the final tests, results_1b9c_f_*.tsv) or an argued equivalence. A fresh mutant that pass 1 killed is not listed.
FRESH_FATE = {
    'F60': 'test: SaSortStatsGap1b9cTests.ThePartitionOfTheRoot_CountsEveryPositionThatItScans',
    'F58': 'test: SesaGap1b9cTests.TheBwtPass_RefusesAWindowOfZeroBytes_BeforeItOpensTheSuffixArray (a read of three bytes is refused at the entry, four pass it)',
    'F53': 'equivalent: the period 256 is the only one whose number of rounds changes (12 instead of 40), and its cover after 12 rounds is the cover after 40 (the same 23 residues: pruned_cover_1b9c.py)',
    'F54': 'test: SaCoverGap1b9cTests.ThePrunedCover_IsThePrototypesCover_... (the periods 1,169, 1,282 and 1,288, whose cover is smaller after a seventh round)',
    'F05': 'test: SesaGap1b9cTests.TheFlushCounter_CountsOneFlushAtATime_AndIsResetWithTheOtherCounters',
    'F45': 'test: SaBuildGap1b9cTests.AWaitInOnePhase_ShowsInTheSecondsOfThatPhaseAlone',
    'F46': 'test: SaGap1b9cTests.ATextAboveTheLimitOfThePositions_... and ATextAtTheLimitOfThePositions_PassesTheEntry... (the build itself, not only PH1 and PH3)',
    'F49': 'test: SaBuildGap1b9cTests.ABuildThatSucceeded_LeavesNoScratchByteInTheAccount_AndPh3ReportsItsOwnPeakWhichHoldsTheRankFile',
    'F62': 'test: SaBuildGap1b9cTests.ABuildThatSucceeded_LeavesNoScratchByteInTheAccount_AndPh3ReportsItsOwnPeakWhichHoldsTheRankFile',
    'F32': 'test: SaGap1b9cTests.ASampleThatExactlyFitsTheBucketCapacity_IsOneLeaf_AndOneMoreRecordIsSplit (the scans of the root)',
    'F24': 'killed by the default set in pass 2 (the test of 1b.9b, SesaGap1b9bTests.ABuildInfoThatWouldOverrunItsRoom_IsRefused_...)',
    'F35': 'test: SesaGap1b9cTests.WithoutIncludePsi_NoPsiSectionIsWritten_AtAnyThreadCount_AndWithItTheSectionIsThere (the PSI section is compared with the oracle at every thread count)',
}
