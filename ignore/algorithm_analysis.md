# Analysis of Algorithms and Parameters for Stream/Path Support

This report identifies all MCP tools and underlying library algorithms that currently rely on large string parameters and should be updated to support streams or file paths to prevent LLM context overflow.

## Summary of Findings

The current implementation of MCP tools typically takes `string content` as an argument. This forces the LLM to read and include the entire file content in the conversation context, which is inefficient for large genomic files (FASTA, FASTQ, VCF, etc.) and core string data.

Most library parsers already support `TextReader` or file paths, making it easy to add `_file` or `_path` versions of these tools.

## Proposed Tools for Update

### 1. Genomic Parsers (Seqeron.Mcp.Parsers)

These tools should have a corresponding version that accepts a `filePath` instead of `content`.

| Tool Name | Large Parameter(s) | Underlying Algorithm | Potential New Tool Name |
| :--- | :--- | :--- | :--- |
| `fasta_parse` | `content` | `FastaParser.Parse` | `fasta_parse_file` |
| `fastq_parse` | `content` | `FastqParser.Parse` | `fastq_parse_file` |
| `fastq_statistics` | `content` | `FastqParser.CalculateStatistics` | `fastq_statistics_file` |
| `fastq_filter` | `content` | `FastqParser.FilterByQuality` | `fastq_filter_file` |
| `fastq_trim_quality`| `content` | `FastqParser.TrimByQuality` | `fastq_trim_quality_file` |
| `fastq_trim_adapter`| `content` | `FastqParser.TrimAdapter` | `fastq_trim_adapter_file` |
| `bed_parse` | `content` | `BedParser.Parse` | `bed_parse_file` |
| `bed_filter` | `content` | `BedParser.Filter...` | `bed_filter_file` |
| `bed_merge` | `content` | `BedParser.MergeOverlapping` | `bed_merge_file` |
| `bed_intersect` | `contentA`, `contentB` | `BedParser.Intersect` | `bed_intersect_files` |
| `vcf_parse` | `content` | `VcfParser.Parse` | `vcf_parse_file` |
| `vcf_statistics` | `content` | `VcfParser.CalculateStatistics" | `vcf_statistics_file` |
| `vcf_filter` | `content` | `VcfParser.Filter...` | `vcf_filter_file` |
| `gff_parse` | `content` | `GffParser.Parse` | `gff_parse_file` |
| `genbank_parse` | `content` | `GenBankParser.Parse` | `genbank_parse_file` |
| `embl_parse` | `content` | `EmblParser.Parse` | `embl_parse_file` |

### 2. Suffix Tree Operations (SuffixTree.Mcp.Core)

These tools involve building a suffix tree, which is a memory-intensive operation for large texts.

| Tool Name | Large Parameter(s) | Potential New Tool Name |
| :--- | :--- | :--- |
| `suffix_tree_contains` | `text` | `suffix_tree_contains_file` |
| `suffix_tree_count` | `text` | `suffix_tree_count_file` |
| `suffix_tree_find_all` | `text` | `suffix_tree_find_all_file` |
| `suffix_tree_lrs` | `text` | `suffix_tree_lrs_file` |
| `suffix_tree_lcs` | `text1`, `text2` | `suffix_tree_lcs_files` |
| `suffix_tree_stats` | `text` | `suffix_tree_stats_file` |

### 3. Sequence Analyzers (SuffixTree.Mcp.Core / Seqeron.Genomics)

These tools analyze sequences which are often stored in files.

| Tool Name | Large Parameter(s) | Potential New Tool Name |
| :--- | :--- | :--- |
| `find_longest_repeat` | `sequence` | `find_longest_repeat_file` |
| `find_longest_common_region`| `sequence1`, `sequence2` | `find_longest_common_region_files` |
| `calculate_similarity` | `sequence1`, `sequence2` | `calculate_similarity_files` |
| `hamming_distance` | `sequence1`, `sequence2" | `hamming_distance_files` |
| `edit_distance" | `sequence1`, `sequence2` | `edit_distance_files` |
| `count_approximate_occurrences`| `sequence` | `count_approximate_occurrences_file` |

## Implementation Recommendation

1.  **Library Level:** Add `BuildFromFile` methods to `SuffixTree` and ensure all parsers in `Seqeron.Genomics` have robust `Stream` support if not already present.
2.  **MCP Level:** Add new tool methods that accept `filePath` (or multiple paths). These methods will:
    *   Validate the path.
    *   Read the file (or use a `Stream`).
    *   Call the existing algorithm.
    *   Return the same results as the content-based tools.
