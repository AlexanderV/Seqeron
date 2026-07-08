# SuffixTree.Mcp.Core

MCP server — **Suffix-tree search, edit/Hamming distance, k-mer similarity.**

Exposes **12 tools** — the same validated `Seqeron.Genomics` algorithms as the C# API, callable over
MCP. Every tool carries an explicit JSON input/output schema and a Schema+Binding test, with a
per-tool doc under [`docs/mcp/tools/core/`](../../../../docs/mcp/tools/core). Rollout status:
[`docs/mcp/MCP_STATUS.md`](../../../../docs/mcp/MCP_STATUS.md).

## Run

```bash
dotnet run --project SuffixTree.Mcp.Core
```

Register it in any MCP client as a stdio server (`command: dotnet`, `args: ["run","--project","SuffixTree.Mcp.Core"]`). New to MCP? The [hub guide](../../../../docs/mcp/README.md) lists all 11 servers and how to wire them up.

## Tools (12)

| Tool | Description |
|------|-------------|
| `calculate_similarity` | Calculate similarity between two DNA sequences using k-mer Jaccard index (0-100 percentage scale). |
| `count_approximate_occurrences` | Count approximate occurrences of a pattern in a sequence, allowing up to maxMismatches substitutions. |
| `edit_distance` | Calculate edit distance (Levenshtein distance) between two sequences. |
| `find_longest_common_region` | Find the longest common region between two DNA sequences. |
| `find_longest_repeat` | Find the longest repeated region in a DNA sequence. |
| `hamming_distance` | Calculate Hamming distance between two sequences of equal length. |
| `suffix_tree_contains` | Check if a pattern exists in text using suffix tree. |
| `suffix_tree_count` | Count the number of occurrences of a pattern in text using suffix tree. |
| `suffix_tree_find_all` | Find all positions where a pattern occurs in text using suffix tree. |
| `suffix_tree_lcs` | Find the longest common substring between two texts using suffix tree. |
| `suffix_tree_lrs` | Find the longest repeated substring in text using suffix tree. |
| `suffix_tree_stats` | Get statistics about a suffix tree: node count, leaf count, max depth, and text length. |
