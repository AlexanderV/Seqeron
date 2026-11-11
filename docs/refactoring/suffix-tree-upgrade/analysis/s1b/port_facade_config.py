"""
Configuration of port_facade_tests.py (work package 1b.6): which existing test files are derived into backend-parameterised copies,
which methods are removed from a derived copy (with the reason), and which files are out of scope (with the reason).
Every removal and every out-of-scope file carries a concrete reason; the traceability table is generated from this file.
"""

UK = 'Ukkonen'
TRYBUILD = 'SuffixTree.TryBuild is a static member of the Ukkonen class that returns false instead of throwing; the other backends have no TryBuild (the factories throw ArgumentNullException / ArgumentException, asserted in FacadeInt32LimitTests)'
EMPTYSINGLETON = 'SuffixTree.Empty is the shared empty-tree singleton of the Ukkonen class; no factory counterpart'
PE = 'Persistent'
MC = 'McpCore'

# origin -> {relative path: {'exclude': {method: reason}, 'subs': [(regex, replacement)], 'note': text, 'no_fixture': bool}}
PORTED = {
    UK: {
        'Algorithms/AllDistinctLongestCommonSubstringsTests.cs': {'subs': [(r'\(SuffixTree inner\)', '(ISuffixTree inner)')]},
        'Algorithms/AllLongestRepeatedSubstringsTests.cs': {},
        'Algorithms/BruteForceVerificationTests.cs': {},
        'Algorithms/DistinctSubstringCountTests.cs': {},
        'Algorithms/ExactMatchAnchorTests.cs': {},
        'Algorithms/KCommonSubstringTests.cs': {},
        'Algorithms/LiteratureExamplesTests.cs': {},
        'Algorithms/LongestCommonSubstringOracleTests.cs': {},
        'Algorithms/LongestCommonSubstringTests.cs': {},
        'Algorithms/LongestRepeatedSubstringTests.cs': {},
        'Algorithms/MaximalMatchTests.cs': {'subs': [(r'\(SuffixTree inner\)', '(ISuffixTree inner)'), (r'private static \(int, int, int\)\[\] Mem', 'private (int, int, int)[] Mem'), (r'private static \(int, int, int\)\[\] Mum', 'private (int, int, int)[] Mum')]},
        'Algorithms/MaximalRepeatedPairsTests.cs': {'subs': [(r'private static \(int, int, int\)\[\] Pairs', 'private (int, int, int)[] Pairs')]},
        'Algorithms/MetamorphicTests.cs': {},
        'Algorithms/SuffixEnumerationTests.cs': {},
        'Algorithms/SuffixTreeCoreContractTests.cs': {},
        'Algorithms/SuffixTreePropertyTests.cs': {},
        'Compatibility/BinaryDataTests.cs': {},
        'Compatibility/UnicodeTests.cs': {},
        'Core/BuildTests.cs': {'subs': [(r'SuffixTree\.Empty\.IsEmpty', 'Tree("").IsEmpty')], 'exclude': {
            'TryBuild_NullString_ReturnsFalse': TRYBUILD,
            'TryBuild_ValidString_ReturnsTrueWithTree': TRYBUILD,
            'TryBuild_EmptyString_ReturnsTrueWithEmptyTree': TRYBUILD,
            'TryBuild_NullTextSource_ReturnsFalse': TRYBUILD,
            'TryBuild_ValidTextSource_ReturnsTrue': TRYBUILD,
            'Empty_StaticProperty_ReturnsEmptyTree': EMPTYSINGLETON,
            'Empty_ReturnsSameInstance': EMPTYSINGLETON,
        }, 'note': 'verbatim assertions through the backend factory; the shared SuffixTree.Empty instance in IsEmpty_CorrelatesWithTextLength is replaced by a freshly built empty tree; Build(ITextSource) cases go through TryBuild of the test base (the Ukkonen null / length semantics on top of the backend factory)'},
        'Core/DiagnosticsTests.cs': {'exclude': {
            'ToString_EmptyTree_ContainsEmptyIndicator': 'Object.ToString() override of the Ukkonen class (format "SuffixTree (Nodes: n, Leaves: l, Text: ...)"); not a member of ISuffixTree',
            'ToString_ShortText_ContainsFullText': 'Object.ToString() override of the Ukkonen class; not a member of ISuffixTree',
            'ToString_LongText_TruncatesAt50Chars': 'Object.ToString() override of the Ukkonen class; not a member of ISuffixTree',
            'EstimatePrintTreeCapacity_IsPositiveAndBounded': 'internal static helper of the Ukkonen PrintTree (initial StringBuilder capacity arithmetic); the ESA PrintTree has no such estimate',
            'EstimatePrintTreeCapacity_IsMonotonic_AcrossTheInt32WrapPoint': 'internal static helper of the Ukkonen PrintTree',
            'EstimatePrintTreeCapacity_SmallTexts_KeepTheEstimateOf100CharsPerTextChar': 'internal static helper of the Ukkonen PrintTree',
        }},
        'Core/InvariantTests.cs': {},
        'Properties/PropertyBasedTests.cs': {},
        'Regression/RegressionTests.cs': {},
        'Robustness/EdgeCaseTests.cs': {},
        'Robustness/Int32LimitGuardTests.cs': {'exclude': {'TryBuild_TextBeyondLimit_ReturnsFalse_DoesNotThrow': TRYBUILD}, 'backend_skip': {
            'Build_TextOneBeyondLimit_ThrowsArgumentException_BeforeReadingContent': {'persistent-heap': 'ParamName differs by design: the Ukkonen constructor names its argument value, the v6 and ESA factories text (asserted in FacadeInt32LimitTests)', 'esa': 'ParamName differs by design: the Ukkonen constructor names its argument value, the v6 and ESA factories text (asserted in FacadeInt32LimitTests)'},
            'Build_TextExactlyAtLimit_PassesTheLengthGuard': {'persistent-heap': 'the v6 heap mode rejects a text of 1073741823 characters earlier, for its estimated storage size (ArgumentException from ThrowIfTooLargeForHeap)'},
        }},
        'Robustness/StressTests.cs': {},
        'Robustness/ThreadSafetyTests.cs': {'subs': [(r'ConcurrentBag<SuffixTree>', 'ConcurrentBag<ISuffixTree>')]},
        'Search/ContainsTests.cs': {},
        'Search/CountOccurrencesTests.cs': {},
        'Search/FindAllOccurrencesTests.cs': {},
    },
}

# origin -> {relative path: {'reason': text, 'methods': {method: own reason}, 'disposition': text, 'target': text}}
NOT_PORTED = {
    UK: {
        'Core/StringTextSourceTests.cs': {'reason': 'class under test is StringTextSource (a text source), not a suffix tree; no ISuffixTree member is asserted'},
        'Algorithms/PerformanceGuardTests.cs': {'reason': 'wall-clock scaling guards of the Ukkonen builder (Build/Contains/FindAll/LCS linearity); the ESA gates are measured under the 1b.9 protocol, a timing assert in the default run would be flaky'},
        'Performance/ComplexityGuardTests.cs': {'reason': 'wall-clock complexity guards of the Ukkonen builder and its search; the ESA gates are measured under the 1b.9 protocol'},
    },
}

PORTED[PE] = {
    'Core/SuffixTreeTestBase.cs': {
        'note': 'the abstract base of the persistent heap/mmf/large-format fixtures; here a concrete fixture per backend',
        'subs': [
            (r'public abstract class SuffixTreeTestBase', '[TestFixture]\n    [Category("Core")]\n    public class SuffixTreeTestBase'),
            (r'protected abstract ISuffixTree CreateTree\(string text\);', 'protected ISuffixTree CreateTree(string text) => Tree(text);'),
        ]},
    'Validation/CrossImplementationContractTests.cs': {},
    'Validation/DiagnosticsContractTests.cs': {},
    'Validation/MathematicalInvariantsTests.cs': {},
    'Validation/EmptyPatternContractTests.cs': {
        'subs': [(r'private PersistentSuffixTree BuildTree\(string text\)\s*\{.*?\n    \}', 'private ISuffixTree BuildTree(string text) => Tree(text);')]},
    'Validation/BoundaryAndEdgeCaseTests.cs': {'exclude': {
        'Mmf_FrequentResizing_MaintainsIntegrity': 'memory-mapped storage with a 64-byte initial capacity (MappedFileStorageProvider, PersistentSuffixTreeBuilder): v6 storage mechanics',
    }},
    'Validation/TreeContractTests.cs': {
        'subs': [(r'private PersistentSuffixTree BuildTree\(string text\)\s*\{.*?\n    \}', 'private ISuffixTree BuildTree(string text) => Tree(text);')]},
    'Validation/SpanOverloadContractTests.cs': {
        'subs': [(r'\[TestCase\(true, TestName = "[^"]*"\)\]\n', ''),
                 (r'private ISuffixTree CreateTree\(string text, bool useMappedStorage\)\s*\{.*?\n        \}', 'private ISuffixTree CreateTree(string text, bool useMappedStorage) => Tree(text);')]},
    'Validation/LcsAndEnumerationContractTests.cs': {
        'exclude': {
            'LcsInfoAndFindAll_MatchReference': 'calls LongestCommonSubstringInfo(ReadOnlySpan<char>), which exists on the concrete PersistentSuffixTree class only and is not a member of ISuffixTree; the string-overload assertions are covered by FacadeChecks.CheckMatching, the span overload of the facade by FacadeSurfaceTests.LongestCommonSubstringInfo_SpanOverload_EqualsTheStringOverload',
        },
        'subs': [(r'\[TestCase\(true, TestName = "[^"]*"\)\]\n', ''),
                 (r'private PersistentSuffixTree CreateTree\(string text, bool useMappedStorage\)\s*\{.*?\n    \}', 'private ISuffixTree CreateTree(string text, bool useMappedStorage) => Tree(text);')]},
    'Parity/ParityTests.cs': {},
    'Parity/TopologyParityTests.cs': {},
    'Safety/NullGuardTests.cs': {'exclude': {
        'Builder_Build_NullText_ThrowsArgumentNullException': 'PersistentSuffixTreeBuilder is the v6 builder, not a suffix-tree member',
        'Load_NullStorage_ThrowsArgumentNullException': 'PersistentSuffixTree.Load takes a v6 storage provider',
    }, 'subs': [(r'Tree\(null!\)', 'Tree((ITextSource)null!)')]},
}

HYBRID_STMT = (r'^[ \t]*using \(var hybrid = PersistentSuffixTreeFactory\.CreateCore\([^\n]*\)\)\n[^\n]*\n', '')
LIFETIME = 'ObjectDisposedException after Dispose is the lifetime contract, asserted for every disposable backend in FacadeLifetimeTests (the Ukkonen tree is not disposable)'
MMF = 'memory-mapped file storage (CreatePersistent with a path / LoadPersistent): v6 storage, not an ISuffixTree behaviour'

PORTED[PE].update({
    'Algorithms/AllLongestRepeatedSubstringsParityTests.cs': {'subs': [HYBRID_STMT], 'exclude': {
        'MemoryMappedTree_EqualsInMemory': MMF, 'DisposedTree_Throws': LIFETIME}},
    'Algorithms/DistinctSubstringCountParityTests.cs': {'exclude': {'DisposedTree_Throws': LIFETIME}},
    'Algorithms/MaximalMatchParityTests.cs': {'subs': [HYBRID_STMT], 'exclude': {
        'MemoryMappedTree_EqualsInMemory': MMF, 'DisposedTree_Throws': LIFETIME}},
    'Algorithms/RepeatedPairsAndCommonSubstringsParityTests.cs': {'subs': [
        (r'^[ \t]*using \(var hybrid = [^\n]*\)\n[ \t]*\{\n.*?\n[ \t]*\}\n', ''),], 'exclude': {
        'MemoryMappedTree_EqualsInMemory': MMF, 'DisposedTree_Throws': LIFETIME}},
})

PORTED[PE].update({
    'Parity/B05SuffixTreeAuditPropertyTests.cs': {'subs': [(r'private static IEnumerable<ISuffixTree> Trees\(', 'private IEnumerable<ISuffixTree> Trees(')],
        'note': 'the reference tree (Ukkonen) stays; the persistent tree of the original is the backend under test'},
    'Parity/B05SuffixTreePropertyParityTests.cs': {'exclude': {
        'RandomUtf16Text_RoundTripsVerbatim_HeapMmfAndSerializer': 'persistence round trips (memory-mapped file, SuffixTreeSerializer export/import)',
        'NodeWithMoreThan256Children_KeepsAllChildren_HeapMmfHybrid': 'checks the v6 child-array layout across heap, memory-mapped and hybrid storage; the ISuffixTree-level behaviour (alphabets of 257, 300, 1000 symbols) is exercised by FacadeFuzzTests and the fuzz families with sigma above 256',
    }},
    'Safety/ConcurrencyTests.cs': {'subs': [(r'var tree = SuffixTree\.Build\(text\);', 'var tree = Tree(text);')], 'exclude': {
        'MmfSuffixTree_ParallelReads_ShouldBeThreadSafe': MMF,
        'ConcurrentDispose_HeapProvider_NoException': 'HeapStorageProvider (v6 storage provider) disposal',
        'ConcurrentDispose_MappedProvider_NoException': 'MappedFileStorageProvider (v6 storage provider) disposal',
        'ConcurrentDispose_MemoryMappedTextSource_NoDoubleRelease': 'MemoryMappedTextSource (v6 text source) disposal',
    }, 'note': 'the in-memory parallel-read test, here for every backend; concurrent Dispose of the tree itself is FacadeLifetimeTests.Dispose_FromSeveralThreads_DisposesOnce_WithoutException'},
})

FL = 'Facade/FacadeLifetimeTests.cs :: FacadeLifetimeTests.'
FI = 'Facade/FacadeInt32LimitTests.cs :: FacadeInt32LimitTests.'
FS = 'Facade/FacadeSurfaceTests.cs :: FacadeSurfaceTests.'
FK = 'Facade/FacadeKnownAnswerTests.cs :: FacadeKnownAnswerTests.'
FD = 'Facade/FacadeDifferentialTests.cs :: FacadeDifferentialTests.'
COVERED = 'covered by new test'
V6STORE = 'v6 storage / loader / builder mechanics (storage providers, header validation, node layout), no ISuffixTree behaviour'

NOT_PORTED[PE] = {
    'Algorithms/SuffixTreeFuzzTests.cs': {'reason': 'corrupts the header and body of a v6 file and expects the v6 loader to throw (PersistentSuffixTreeFactory.Load)'},
    'Core/BuilderGuardTests.cs': {'reason': 'PersistentSuffixTreeBuilder guards (second Build call, deepest-node offset, in-memory collections); ' + V6STORE},
    'Core/HeapSuffixTreeTests.cs': {'reason': 'no test methods of its own: the heap fixture of SuffixTreeTestBase, which is ported (SuffixTreeTestBase.cs, persistent-heap fixture)'},
    'Core/LargeFormatSuffixTreeTests.cs': {'reason': 'no test methods of its own: the large-format fixture of SuffixTreeTestBase (v6 Large layout); the base is ported'},
    'Core/MmfSuffixTreeTests.cs': {'reason': 'memory-mapped fixtures of the v6 tree (file creation, reload, trim); the base assertions are ported in SuffixTreeTestBase.cs'},
    'Core/PersistentSuffixTreeFactoryLifecycleTests.cs': {'reason': 'concrete PersistentSuffixTree type and v6 file loading', 'methods': {
        'CreatePersistent_NullText_ThrowsArgumentNullException': (COVERED, FI + 'Create_NullText_ThrowsArgumentNullException', 'the factory of the ESA refuses a null text with ArgumentNullException'),
    }},
    'Format/HybridTransitionZoneTests.cs': {'reason': 'v6 binary format: hybrid compact/large transition zones'},
    'Format/MagicNumberTests.cs': {'reason': 'v6 binary format: magic number and version'},
    'Format/StorageFormatTests.cs': {'reason': 'v6 binary format: header and node layout'},
    'Parity/SafeVsUnsafeDifferentialTests.cs': {'reason': 'safe against unsafe storage providers of the v6 tree'},
    'Parity/UnsafeParityTests.cs': {'reason': 'unsafe storage provider of the v6 tree against the in-memory tree; ' + V6STORE},
    'Safety/AsciiMemoryMappedTextSourceTests.cs': {'reason': 'AsciiMemoryMappedTextSource (v6 text source), not a suffix tree'},
    'Safety/MemoryMappedTextSourceTests.cs': {'reason': 'MemoryMappedTextSource (v6 text source), not a suffix tree'},
    'Safety/NodeApiSafetyTests.cs': {'reason': 'PersistentSuffixTreeNode / navigator API of the v6 tree (internal node handles)'},
    'Safety/OverflowAndValidationTests.cs': {'reason': 'v6 text sources and storage providers (substring overflow, zero-capacity heap provider)'},
    'Safety/ProviderSafetyTests.cs': {'reason': 'v6 storage providers'},
    'Safety/StorageBoundsOverflowTests.cs': {'reason': 'v6 storage bounds arithmetic'},
    'Safety/DisposeBehaviorTests.cs': {'reason': V6STORE, 'methods': {
        'Dispose_HeapTree_AllMethodsThrowObjectDisposedException': (COVERED, FL + 'UseAfterDispose_EverySearchAnalysisAndDiagnosticsMember_ThrowsObjectDisposedException', 'same member list plus the navigation and suffix-array members, for every disposable backend'),
        'Dispose_DoubleDispose_DoesNotThrow': (COVERED, FL + 'Dispose_IsIdempotent', 'idempotent Dispose for every disposable backend'),
        'FindAllLCS_ReturnsAllPositions': (COVERED, FS + 'FindAllLongestCommonSubstrings_ReturnsAllPositions_AndNoMatchIsEmpty', 'the two FindAllLongestCommonSubstrings cases'),
        'FindAllLCS_NoMatch_ReturnsEmpty': (COVERED, FS + 'FindAllLongestCommonSubstrings_ReturnsAllPositions_AndNoMatchIsEmpty', 'the two FindAllLongestCommonSubstrings cases'),
        'ReadOnlyProvider_WriteThrowsInvalidOperation': 'read-only MappedFileStorageProvider write rejection',
        'TrimToSize_ReducesFileSize': 'memory-mapped file trimming',
        'Dispose_TextSourceThrows_StorageStillDisposed': 'order of disposal of the v6 text source and storage',
        'Dispose_StorageThrows_CompletesWithoutCorruption': 'disposal of a throwing v6 storage provider',
        'HeapStorageProvider_DoubleDispose_DoesNotThrow': 'HeapStorageProvider (v6 storage provider)',
    }},
    'Safety/Int32LimitGuardTests.cs': {'reason': V6STORE, 'methods': {
        'Builder_Build_TextOneBeyondLimit_ThrowsArgumentException_BeforeReadingContent': (COVERED, FI + 'Create_TextOneBeyondTheLimit_ThrowsArgumentException_BeforeReadingContent', 'through the ESA factory (the original calls PersistentSuffixTreeBuilder.Build)'),
        'Builder_Build_TextOfIntMaxValue_ThrowsArgumentException': (COVERED, FI + 'Create_TextOfIntMaxValue_ThrowsArgumentException', 'through the ESA factory'),
        'Factory_Create_HeapMode_TextBeyondLimit_ThrowsArgumentException_WithoutAllocatingStorage': (COVERED, FI + 'Create_TextOneBeyondTheLimit_ThrowsArgumentException_BeforeReadingContent', 'ParamName text, content not read'),
        'MaxTextLength_IsTheLargestLengthWhoseWorstCaseNodeCountFitsInt32': (COVERED, FI + 'MaxTextLength_IsTheLargestLengthWhoseWorstCaseNodeCountFitsInt32', 'same assertion'),
        'Builder_Build_TextExactlyAtLimit_PassesTheLengthGuard': (COVERED, FI + 'Create_TextExactlyAtTheLimit_PassesTheLengthGuard_AndAllocatesNoGigabytes', 'through the ESA factory'),
        'Builder_Build_RejectedLength_DoesNotConsumeTheBuilder': 'PersistentSuffixTreeBuilder state after a rejected length (the factory holds no state)',
        'Builder_NodeCountLimit_ExactlyEnough_Builds_OneLess_ThrowsInvalidOperationException': 'NodeCountLimit of the v6 builder; the ESA counterpart is the NodeCount guard of the facade, FacadeInt32LimitTests.NodeCount_AboveInt32_IsRefusedByTheFacade_AndEverythingElseKeepsWorking',
        'Builder_NodeCountLimit_DefaultsToInt32MaxValue': 'NodeCountLimit of the v6 builder',
        'Factory_Create_HeapMode_TextTooLargeForManagedArray_ThrowsClearArgumentException': 'v6 heap storage estimate of 72 bytes per character; the ESA has no such estimate',
        'HeapProvider_EnsureCapacity_AtIntMaxValue_ThrowsInvalidOperationException': 'HeapStorageProvider',
        'HeapProvider_EnsureCapacity_JustAboveArrayMaxLength_ThrowsInvalidOperationException': 'HeapStorageProvider',
        'HeapProvider_Constructor_InitialCapacityAboveArrayMaxLength_ThrowsArgumentOutOfRange': 'HeapStorageProvider',
        'HeapProvider_FailedGrowth_LeavesProviderUsable': 'HeapStorageProvider',
        'Load_NonMmfStorage_AsciiTextLongerThanByteArrayLimit_ThrowsInvalidOperationException': 'v6 loader',
        'Load_NonMmfStorage_Utf16TextLongerThanByteArrayLimit_ThrowsInvalidOperationException': 'v6 loader',
        'Load_NonMmfStorage_Utf16TextBeyondInt32Bytes_ThrowsInvalidOperationException': 'v6 loader',
    }},
    'Serialization/ImportTruncationTests.cs': {'reason': 'SuffixTreeSerializer import of a truncated stream (logical persistence format)'},
    'Serialization/LogicalPersistenceTests.cs': {'reason': 'SuffixTreeSerializer export/import and file round trips (logical persistence format)'},
    'Serialization/SerializerHashTests.cs': {'reason': 'builds the v6 trees with PersistentSuffixTreeBuilder', 'methods': {
        'CalculateLogicalHash_IsDeterministic': (COVERED, FD + 'LogicalHash_EqualsThatOfTheOtherBackends', 'the logical hash of the facade equals that of the Ukkonen and v6 trees for the same text, so it is deterministic and layout independent'),
        'CalculateLogicalHash_UsesLittleEndianByteOrder': 'asserts the hash of a v6 tree is 32 bytes and repeatable; the byte order is an implementation detail of the serializer',
    }},
    'Serialization/SetSizeAndExportTests.cs': {'reason': 'SetSize and Export of the v6 serializer / storage'},
    'Validation/LoadValidationTests.cs': {'reason': 'v6 header validation on load'},
    'Validation/StrictLoadValidationTests.cs': {'reason': 'v6 structure validation on load'},
    'Validation/PersistentTextFidelityAndParityTests.cs': {'reason': 'v6 text persistence (heap, memory-mapped, hybrid, reload, export)', 'methods': {
        'Parity_Heap_Mmf_Hybrid_Reload': (COVERED, FD + 'NamedTexts_EveryMember', 'the heap part (node and leaf counts, Traverse stream, LRS, suffixes, LCS, anchors, occurrences against the Ukkonen tree) for the ESA; the memory-mapped, hybrid and reload parts are v6 storage'),
    }},
}

# Mcp.Core: the MCP wrappers build their own Ukkonen tree; the tests assert argument validation and the DTO mapping, plus literal answers
MCPW = 'asserts the argument validation / DTO mapping of the MCP wrapper (SuffixTreeCoreTools), which constructs its own tree internally'
NOT_PORTED[MC] = {
    'Contracts/SizeLimitTests.cs': {'reason': 'size limits of the MCP tool inputs'},
    'Contracts/SuffixTreeCoreWrapperParityTests.cs': {'reason': 'wrapper against library parity for the MCP layer'},
    'Contracts/ToolRegistrationContractTests.cs': {'reason': 'MCP tool registration metadata'},
    'Genomics/CalculateSimilarityTests.cs': {'reason': 'genomics helper (similarity), not a suffix-tree member'},
    'Genomics/CountApproximateOccurrencesTests.cs': {'reason': 'genomics helper (approximate matching), not a suffix-tree member'},
    'Genomics/EditDistanceTests.cs': {'reason': 'genomics helper (edit distance), not a suffix-tree member'},
    'Genomics/FindLongestCommonRegionTests.cs': {'reason': 'genomics wrapper around the LCS; ' + MCPW},
    'Genomics/FindLongestRepeatTests.cs': {'reason': 'genomics wrapper around the LRS; ' + MCPW},
    'Genomics/HammingDistanceTests.cs': {'reason': 'genomics helper (Hamming distance), not a suffix-tree member'},
    'Core/SuffixTreeContainsTests.cs': {'reason': MCPW, 'methods': {'SuffixTreeContains_ReturnsExpectedResult': (COVERED, FK + 'ContainsCountFindAll_LiteralAnswers', 'the six (text, pattern, expected) rows')}},
    'Core/SuffixTreeCountTests.cs': {'reason': MCPW, 'methods': {'SuffixTreeCount_ReturnsExpectedCount': (COVERED, FK + 'ContainsCountFindAll_LiteralAnswers', 'the six (text, pattern, count) rows')}},
    'Core/SuffixTreeFindAllTests.cs': {'reason': MCPW, 'methods': {
        'SuffixTreeFindAll_ReturnsExactPositions': (COVERED, FK + 'ContainsCountFindAll_LiteralAnswers', 'literal position lists'),
        'SuffixTreeFindAll_ReturnsPositionsInAscendingOrder': (COVERED, FK + 'ContainsCountFindAll_LiteralAnswers', 'literal position lists (sorted: the order of FindAllOccurrences is engine order by contract, the MCP wrapper sorts)'),
    }},
    'Core/SuffixTreeLcsTests.cs': {'reason': MCPW, 'methods': {
        'SuffixTreeLcs_GusfieldExampleAndTieBreak_MatchBruteForce': (COVERED, FK + 'LongestCommonSubstring_LiteralAnswers', 'Gusfield example and the tie rule'),
        'SuffixTreeLcs_ReturnsExpectedForRepresentativeInputs': (COVERED, FK + 'LongestCommonSubstring_LiteralAnswers', 'banana/panama, abcdef/xyzcdeww, no match'),
        'SuffixTreeLcs_ResultSatisfiesCoreContract': (COVERED, FK + 'LongestCommonSubstring_LiteralAnswers', 'length equals the dynamic-programming length, the result occurs in both texts'),
    }},
    'Core/SuffixTreeLrsTests.cs': {'reason': MCPW, 'methods': {
        'SuffixTreeLrs_ClassicStrings_MatchBruteForce': (COVERED, FK + 'LongestRepeatedSubstring_LiteralAnswers', 'mississippi, abracadabra, ATATAT'),
        'SuffixTreeLrs_ReturnsExpectedForRepresentativeInputs': (COVERED, FK + 'LongestRepeatedSubstring_LiteralAnswers', 'banana, aaaaa, abcdef'),
    }},
    'Core/SuffixTreeAllLcsTests.cs': {'reason': MCPW, 'methods': {
        'SuffixTreeAllLcs_Ties_AllReported': (COVERED, FK + 'AllDistinctLongestCommonSubstrings_LiteralAnswers', 'abcxyzabc / xyzqabcpxyz'),
        'SuffixTreeAllLcs_NothingShared_Empty': (COVERED, FK + 'AllDistinctLongestCommonSubstrings_LiteralAnswers', 'AAAA / CCC'),
        'SuffixTreeAllLcs_EqualsLibrary': (COVERED, FK + 'AllDistinctLongestCommonSubstrings_LiteralAnswers', 'GATTACA / TAGACCA'),
    }},
    'Core/SuffixTreeAllLrsTests.cs': {'reason': MCPW, 'methods': {
        'SuffixTreeAllLrs_Ties_AllReportedInFirstOccurrenceOrder': (COVERED, FK + 'AllLongestRepeatedSubstrings_LiteralAnswers', 'abcxbcaxcab'),
        'SuffixTreeAllLrs_SingleRepeat_MatchesLockedCases': (COVERED, FK + 'AllLongestRepeatedSubstrings_LiteralAnswers', 'banana, mississippi, aaaa'),
        'SuffixTreeAllLrs_NoRepeat_Empty_AndEqualsLibrary': (COVERED, FK + 'AllLongestRepeatedSubstrings_LiteralAnswers', 'abcd, xyzzyxxz'),
    }},
    'Core/SuffixTreeFindMemsTests.cs': {'reason': MCPW, 'methods': {
        'SuffixTreeFindMems_GattacaRepeat_MatchesMummer': (COVERED, FK + 'MaximalMatches_LiteralAnswers', 'GATTACAGATTACA / ATTACAGGATTACAT'),
        'SuffixTreeFindMems_ThreeSetsDiffer_MatchesMummer_AndEqualsLibrary': (COVERED, FK + 'MaximalMatches_LiteralAnswers', 'MUMmer 3.23 three sets'),
    }},
    'Core/SuffixTreeFindMumsTests.cs': {'reason': MCPW, 'methods': {
        'SuffixTreeFindMums_Both_MatchesMummerMum': (COVERED, FK + 'MaximalMatches_LiteralAnswers', 'MUMmer -mum'),
        'SuffixTreeFindMums_Reference_MatchesMummerMumreference_AndEqualsLibrary': (COVERED, FK + 'MaximalMatches_LiteralAnswers', 'MUMmer -mumreference'),
        'SuffixTreeFindMums_LengthOne_KeepsMumAtReferenceStart': (COVERED, FK + 'MaximalMatches_LiteralAnswers', 'ACGT / TTAG, length 1'),
    }},
    'Core/SuffixTreeKCommonSubstringsTests.cs': {'reason': MCPW, 'methods': {
        'SuffixTreeKCommonSubstrings_RosalindLcsmSample': (COVERED, FK + 'KCommonSubstrings_LiteralAnswers', 'Rosalind LCSM sample'),
        'SuffixTreeKCommonSubstrings_MinSupportTwo': (COVERED, FK + 'KCommonSubstrings_LiteralAnswers', 'minimum support 2'),
        'SuffixTreeKCommonSubstrings_NothingShared_Empty': (COVERED, FK + 'KCommonSubstrings_LiteralAnswers', 'AAA / CCC'),
    }},
    'Core/SuffixTreeMaximalRepeatsTests.cs': {'reason': MCPW, 'methods': {
        'SuffixTreeMaximalRepeats_GusfieldExample_MatchesRepeatMatch': (COVERED, FK + 'MaximalRepeatedPairs_LiteralAnswers', 'xabcyabcwabcyz'),
        'SuffixTreeMaximalRepeats_UniqueSymbols_AndNonAcgt': (COVERED, FK + 'MaximalRepeatedPairs_LiteralAnswers', 'ACGTNACGTNACGT with and without the unique symbol N'),
    }},
    'Core/SuffixTreeStatsTests.cs': {'reason': MCPW, 'methods': {
        'SuffixTreeStats_MatchesFundamentalInvariants': (COVERED, FK + 'Statistics_FundamentalInvariants', 'leaf count, depth and node count bounds'),
        'SuffixTreeStats_IsConsistentWithOtherTools': (COVERED, FK + 'Statistics_FundamentalInvariants', 'MaxDepth at least the LRS length, CountOccurrences("") equals LeafCount'),
    }},
    'Core/SuffixTreeCoreConsistencyTests.cs': {'reason': MCPW, 'methods': {
        'SuffixTreeCoreTools_AreMutuallyConsistent_AndMatchBruteForceOracle': (COVERED, FS + 'SearchMembers_AgainstTheDefinition_OnTheCorpus', 'Contains, count and positions against the brute-force oracle on a pattern set per text'),
    }},
}
