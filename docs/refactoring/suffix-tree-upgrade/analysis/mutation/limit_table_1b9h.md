# work package 1b.9h, item 1: the limits of the time limit per test: the classes (the slowest quiet run is the baseline default set, the slowest loaded run the default classes with four busy processes beside the host), then the tests that have a limit of their own (10 x the larger of the two runs, at least 30 s)

| class | guarded tests | limit of the class | slowest quiet run (s) | slowest loaded run (s) | longest limit of a test of the class (s) |
|---|---|---|---|---|---|
| `SesaBuilderTests` | 16 | 30 s | 32.19 | 48.95 | 490 |
| `SaBuilderTests` | 36 | 30 s | 32.07 | 44.43 | 450 |
| `ExternalDcsTests` | 29 | 30 s | 18.67 | 28.72 | 290 |
| `SesaGap1b9cTests` | 17 | 30 s | 1.66 | 27.72 | 280 |
| `SesaFailureTests` | 17 | 30 s | 15.53 | 22.17 | 230 |
| `SesaGapTests` | 8 | 30 s | 10.52 | 12.44 | 130 |
| `Ph6Tests` | 16 | 30 s | 10.93 | 11.19 | 120 |
| `BucketSorterTests` | 19 | 30 s | 1.43 | 9.82 | 100 |
| `SesaSectionTests` | 21 | 30 s | 2.29 | 6.75 | 70 |
| `FullDisk1b9dTests` | 10 | 30 s | 17.66 | 6.30 | 180 |
| `SesaOverflowSectionTests` | 1 | 30 s | 4.21 | 5.84 | 60 |
| `SaGap1b9cTests` | 25 | 30 s | 3.61 | 5.62 | 60 |
| `SaPinnedTests` | 5 | 30 s | 2.73 | 4.01 | 50 |
| `SaGapTests` | 15 | 30 s | 2.60 | 3.87 | 40 |
| `SesaAuditTests` | 22 | 30 s | 2.86 | 2.95 | 30 |
| `SaPlumbingTests` | 5 | 30 s | 1.38 | 2.52 | 30 |
| `SaBuildGap1b9cTests` | 12 | 30 s | 0.99 | 1.23 | 30 |
| `SaDcsGap1b9cTests` | 6 | 30 s | 0.28 | 0.49 | 30 |
| `SaPositionGap1b9cTests` | 10 | 30 s | 0.02 | 0.04 | 30 |
| `SaSortStatsGap1b9cTests` | 4 | 30 s | 0.01 | 0.02 | 30 |
| `SaNamesGap1b9cTests` | 1 | 30 s | 0.00 | 0.02 | 30 |

tests whose limit is below 10 x the slowest run: []
tests with their own limit:
| `SesaBuilderTests.TheBytes_DoNotDependOnThreadsSizesSeedsOrArenaCaps` | 32.19 | 48.95 | 490 |
| `SaBuilderTests.TheOutput_IsIdenticalForEveryThreadCount_BucketCapacity_AndSeed` | 32.07 | 44.43 | 450 |
| `SesaBuilderTests.EveryStringOfSmallAlphabets_GivesTheInMemoryIndex` | 5.17 | 39.85 | 400 |
| `ExternalDcsTests.TheOutput_IsIdentical_ForEveryThreadCount_BucketCapacity_AndSeed` | 18.67 | 28.72 | 290 |
| `SesaGap1b9cTests.WithoutIncludePsi_NoPsiSectionIsWritten_AtAnyThreadCount_AndWithItTheSectionIsThere` | 0.54 | 27.72 | 280 |
| `SesaFailureTests.AFailureAtEveryHookPoint_LeavesNothing_AndKeepsTheTypeOfTheException` | 15.53 | 22.17 | 230 |
| `FullDisk1b9dTests.AVolumeOfAnySize_GivesASuccessOrAnIoException_NeverASignal_WhateverThePreflightKnows` | 17.66 | 6.30 | 180 |
| `FullDisk1b9dTests.AFullDisk_AtEveryBoundaryOfTheBuild_IsASuccessOrAnIoException_NeverASignal_AndNothingIsLeft` | 15.44 | 4.03 | 160 |
| `ExternalDcsTests.EveryString_OverSigma1To3_OfTheDefaultTier_IsTheOracle` | 5.97 | 14.15 | 150 |
| `SesaFailureTests.AFailureInsideAStage_LeavesNothing_AtEveryTracePoint` | 7.27 | 12.23 | 130 |
| `SesaGapTests.Q52_ThePh3Cap_IsEnforcedOnItsOwn_WhenPh1Fits` | 10.52 | 12.44 | 130 |
| `SesaBuilderTests.TheFamilies_GiveTheInMemoryIndex_AtTheDefaultAndATinyGeometry` | 10.56 | 11.35 | 120 |
| `Ph6Tests.TheDefaultWindows_OfARealRun_AreNotSmallerThanOneMiB` | 10.93 | 11.19 | 120 |
| `FullDisk1b9dTests.ADiskThatFillsInTheReservationOfTheTail_IsAnIoException_AndLeavesNothing` | 9.45 | 3.47 | 100 |
| `BucketSorterTests.EveryString_OverSigma1To3_OfTheDefaultTier_IsTheSuffixArray` | 1.25 | 9.82 | 100 |
| `SaBuilderTests.Families_WithTheDefaultParameters_AndWithTheScaledGeometry_AreTheSuffixArray` | 7.20 | 3.67 | 80 |
| `Ph6Tests.EveryStringOfSmallAlphabets_GivesTheDefinitionsBwtAndPsi` | 6.77 | 7.13 | 80 |
| `FullDisk1b9dTests.AFullDisk_AtEveryPointThatThePhasesReport_IsASuccessOrAnIoException_NeverASignal_AndNothingIsLeft` | 6.04 | 0.99 | 70 |
| `SesaSectionTests.TheSectionStage_EveryStringOfSmallAlphabets_EqualsTheOracle` | 2.29 | 6.75 | 70 |
| `SesaFailureTests.TheLatePreflight_AfterB2_CountsTheGrowthAndTheSpillBoundOfTheStack_ExactlyToTheByte` | 3.57 | 5.06 | 60 |
| `SesaOverflowSectionTests.AnOverflowListOfAMebibyte_IsAlignedToTheLargeBoundary_AndTheWholeFileConforms` | 4.21 | 5.84 | 60 |
| `SaGap1b9cTests.TheSortOfALevel_SpillsAtLeast64BytesWhateverTheIoBufferIs_SoTinyBuffersAndTheFloorGiveTheSameFiles` | 3.61 | 5.62 | 60 |
| `SaBuilderTests.Families_AtEveryTinyVariant_AreTheSuffixArray` | 3.83 | 5.85 | 60 |
| `SaPinnedTests.TheBuild_HasThePinnedShape_AndCounts` | 2.73 | 4.01 | 50 |
| `ExternalDcsTests.ALeafThatOutgrowsItsCapacity_IsSplitAgain_AndTheRanksAreStillRight` | 2.30 | 3.63 | 40 |
| `SaGapTests.APh1ThatFailsAtAnyCallOfTheHook_OrUnderAnyArenaCap_LeavesNothing_AndGivesTheAccountBack` | 2.60 | 3.87 | 40 |
| `FullDisk1b9dTests.TheVolumeFillsRightAfterTheGrowth_TheLcpIsStillWritten_AndTheBuildEndsWithoutASignal` | 3.50 | 1.71 | 40 |
| `SaBuilderTests.SmallTexts_AreTheSuffixArray_WithRotatingTinyParameters` | 0.65 | 3.12 | 40 |
| `BucketSorterTests.TheOutput_IsIdenticalForEveryThreadCount_LeafCapacity_FanOut_AndSeed` | 1.43 | 3.79 | 40 |
