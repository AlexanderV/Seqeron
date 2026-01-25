# MCP Tools Traceability Matrix

| ToolName | Server | MethodId | DocRef | Tests | mcp.json | md | Status |
|----------|--------|----------|--------|-------|----------|-----|--------|
| suffix_tree_contains | Core | SuffixTree.Contains | SuffixTree.Search.cs#L22 | 2/2 | ✓ | ✓ | Ready |
| suffix_tree_count | Core | SuffixTree.CountOccurrences | SuffixTree.Search.cs#L132 | 2/2 | ✓ | ✓ | Ready |
| suffix_tree_find_all | Core | SuffixTree.FindAllOccurrences | SuffixTree.Search.cs#L115 | 2/2 | ✓ | ✓ | Ready |
| suffix_tree_lrs | Core | SuffixTree.LongestRepeatedSubstring | SuffixTree.LRS.cs#L20 | 2/2 | ✓ | ✓ | Ready |
| suffix_tree_lcs | Core | SuffixTree.LongestCommonSubstring | SuffixTree.LCS.cs#L20 | 2/2 | ✓ | ✓ | Ready |
| suffix_tree_stats | Core | SuffixTree.NodeCount/LeafCount/MaxDepth | SuffixTree.cs | 2/2 | ✓ | ✓ | Ready |
| find_longest_repeat | Core | GenomicAnalyzer.FindLongestRepeat | GenomicAnalyzer.cs#L20 | 2/2 | ✓ | ✓ | Ready |
| find_longest_common_region | Core | GenomicAnalyzer.FindLongestCommonRegion | GenomicAnalyzer.cs#L178 | 2/2 | ✓ | ✓ | Ready |
| calculate_similarity | Core | GenomicAnalyzer.CalculateSimilarity | GenomicAnalyzer.cs#L238 | 2/2 | ✓ | ✓ | Ready |
| hamming_distance | Core | ApproximateMatcher.HammingDistance | ApproximateMatcher.cs#L163 | 2/2 | ✓ | ✓ | Ready |
| edit_distance | Core | ApproximateMatcher.EditDistance | ApproximateMatcher.cs#L186 | 2/2 | ✓ | ✓ | Ready |
| count_approximate_occurrences | Core | ApproximateMatcher.CountApproximateOccurrences | ApproximateMatcher.cs#L283 | 2/2 | ✓ | ✓ | Ready |
| dna_validate | Sequence | DnaSequence.TryCreate | DnaSequence.cs#L129 | 2/2 | ✓ | ✓ | Ready |
| dna_reverse_complement | Sequence | DnaSequence.GetReverseComplementString | DnaSequence.cs#L149 | 2/2 | ✓ | ✓ | Ready |
| rna_validate | Sequence | RnaSequence.TryCreate | RnaSequence.cs#L176 | 2/2 | ✓ | ✓ | Ready |
| rna_from_dna | Sequence | RnaSequence.FromDna | RnaSequence.cs#L147 | 2/2 | ✓ | ✓ | Ready |
| protein_validate | Sequence | ProteinSequence.TryCreate | ProteinSequence.cs#L357 | 2/2 | ✓ | ✓ | Ready |
| nucleotide_composition | Sequence | SequenceStatistics.CalculateNucleotideComposition | SequenceStatistics.cs#L48 | 2/2 | ✓ | ✓ | Ready |
| amino_acid_composition | Sequence | SequenceStatistics.CalculateAminoAcidComposition | SequenceStatistics.cs#L98 | 2/2 | ✓ | ✓ | Ready |
| molecular_weight_protein | Sequence | SequenceStatistics.CalculateMolecularWeight | SequenceStatistics.cs#L159 | 2/2 | ✓ | ✓ | Ready |
| molecular_weight_nucleotide | Sequence | SequenceStatistics.CalculateNucleotideMolecularWeight | SequenceStatistics.cs#L180 | 2/2 | ✓ | ✓ | Ready |
| isoelectric_point | Sequence | SequenceStatistics.CalculateIsoelectricPoint | SequenceStatistics.cs#L228 | 2/2 | ✓ | ✓ | Ready |
