# K0 gaps (cells <= 2) and hard limits. Paths are repo-relative, audited 2026-11-02.

## Cells scored <= 2
- C2 (1): no wildcard/class matching in tree (ISuffixTreeSearch.cs:39-71 take literal string/span). Only MotifFinder.cs:145 ScanDegenerate, an O(n*m) scan; DnaSequence.cs:100-110 rejects N/IUPAC in the sequence.
- C3 (1): no mismatch/edit descent over nodes. ApproximateMatcher.cs:42 (Hamming window scan) and :137 (Sellers DP per start, O(n*m*(m+k))) never touch SuffixTree (grep: no tree use in that file).
- C5 (2): ms[] never returned; AdvanceMatchingStatistics (SuffixTreeAlgorithms.cs:392) is private. Public FindExactMatchAnchors (ISuffixTreeAnalysis.cs:85) gives one first-peak anchor per run only.
- C6 (2): KCommon.cs:38,90 build a throw-away concatenated tree per call (RunColorCounting :181 disposes it); no document listing/counting, no reusable multi-text index. Separators come from unused UTF-16 units (:123-155), so k is capped by unused code units; total length > int.MaxValue throws (:136).
- C7 (1): nothing strand-aware in src/SuffixTree. MEM/MUM docs say forward strand only (ISuffixTreeAnalysis.cs:96,114). RC repeats (RepeatFinder.cs:4753) use an internal prefix-doubling SA (SequenceComplexity.cs:1909, internal, O(n log^2 n)), not the tree.
- C8 (2): ISuffixTreeNavigator is public (ISuffixTreeNavigator.cs:13) but both implementations are internal (SuffixTree.Navigator.cs:13, PersistentSuffixTreeNavigator.cs:11); SuffixTreeNode is internal (SuffixTreeNode.cs:14); PersistentSuffixTreeNode.TryGetChild/SuffixLink internal (:71,:43). No public suffix link, parent, LCA.
- C10 (1): sequential Ukkonen; measured 10 Mbp = 23.3 s, ~1.7 GB (K-mer_Positions.md:86). Persistent build = Ukkonen + FinalizePipeline + SequentialFinalize passes over MMF.
- C11 (1): SuffixTreeNode = 4 inline keys + 4 child refs + Dictionary ref + terminator ref + link (SuffixTreeNode.cs:34-40); the "20 bytes per character" claim (SuffixTree.cs:76) is contradicted by the repo's own 1.7 GB / 10 Mbp. Persistent: 24 B node + 8 B child entry per node, 72 B/char reserved (PersistentSuffixTreeFactory.cs:156-161) plus temp 40 B/char children + 10 B/char depth (:69,:73).
- C13 (2): text length int32 (ITextSource.cs:17, PersistentConstants.cs:57), NODE_COUNT int32 (:63), `int _nodeCount` (PersistentSuffixTreeBuilder.cs:35, `_nodeCount++` Core.cs:180, unchecked).
- C14 (2): see C13; also managed-heap Dictionary<long,long> _crossZoneSuffixLinks and List _deferredSuffixLinks (PersistentSuffixTreeBuilder.cs:76,81) grow with cross-zone links.
- C15 (0): no Parallel/Task/Thread in src/SuffixTree (grep).
- C17 (0): ExtendTree is private (SuffixTree.Construction.cs:6); ctor consumes whole text (SuffixTree.cs:125-139).
- C18 (0): no compressed index in repo (grep FM-index/BWT/wavelet: no hits in src).
- C19 (2): keys are char/uint (SuffixTreeNode.cs:46,86; PersistentSuffixTreeNode.cs:71); no int alphabet, no 2-bit path.
- C22 (1): no PackageId in src/SuffixTree/**/*.csproj; apps/SuffixTree.* are harnesses.
- C24 (2): git tag list empty; Version defaults to 1.0.0 (Directory.Build.props).

## Hard caps / limits noticed
- .NET string max ~1.07 Gchar caps StringTextSource inputs; custom ITextSource goes through virtual indexer / IEnumerable<char> (SuffixTree.cs:128-139, Builder.cs:300-306), much slower.
- 2n+1 nodes overflow int32 NodeCount above ~1.07 Gchar with a custom ITextSource; no guard.
- HeapStorageProvider hard cap 2 GB (HeapStorageProvider.cs:55-56): heap-backed persistent tree ~29 M chars; Create(filePath:null) clamps to int.MaxValue (PersistentSuffixTreeFactory.cs:176).
- PersistentSuffixTree.LoadStringInternal throws above int.MaxValue bytes (PersistentSuffixTree.cs:80-83) when storage is not MMF.
- Compact zone ends at 0xFFFFFFFE; hybrid jump table reserve is 25% of compact zone, 16 MB floor (Builder.cs:62-66); a text needing more cross-zone entries than the reserve is not guarded here.
- Query hot loops read raw pointers with Debug.Assert-only bounds (HybridLayout.cs:136-160, MappedFileStorageProvider.Unchecked.cs:21-34); Load validates header only (PersistentStorageHeaderValidation.cs), so a corrupt body of a valid-header file is not detected.
- Positions/counts are int (ISuffixTreeSearch.cs:55,71); FindAllOccurrences returns tree order, not sorted (:49-54). MEM/MUM/pairs return fully materialised List (no streaming); query must be a string.
- GetAllSuffixes/EnumerateSuffixes O(n^2) (Suffix_Tree.md:43-44); GenomicAnalyzer.FindRepeats uses GetAllSuffixes (GenomicAnalyzer.cs:60) so it is toy-sized only.
- PrintTree capacity `estimatedNodes * 50` is int arithmetic that can wrap for n > ~21 M (SuffixTree.Diagnostics.cs:81-82).
- MCP tools rebuild the tree on every call from a JSON string (SuffixTreeCoreTools.cs:18 etc.); no session/persisted index, no size limit coded.
- Recursion: none found in tree code; traversals use explicit stacks (SuffixTree.Diagnostics.cs:13-62, README "Iterative Algorithms").
