# SuffixTree.Mcp.Core

MCP server — **Suffix-tree search, edit/Hamming distance, k-mer similarity.**

Exposes **18 tools** — the same validated `Seqeron.Genomics` algorithms as the C# API, callable over
MCP. Every tool carries an explicit JSON input/output schema and a Schema+Binding test, with a
per-tool doc under [`docs/mcp/tools/core/`](../../../../docs/mcp/tools/core). Rollout status:
[`docs/mcp/MCP_STATUS.md`](../../../../docs/mcp/MCP_STATUS.md).

## Run

```bash
dotnet run --project SuffixTree.Mcp.Core
```

Register it in any MCP client as a stdio server (`command: dotnet`, `args: ["run","--project","SuffixTree.Mcp.Core"]`). New to MCP? The [hub guide](../../../../docs/mcp/README.md) lists all 11 servers and how to wire them up.

## Size limit

The `suffix_tree_*` tools, `find_longest_repeat` and `find_longest_common_region` rebuild an in-memory suffix tree from their
string arguments on every call (measured: 169 bytes of heap per character for random DNA, up to 193 for repetitive text, about
15 s per 10 Mbp). A call whose text arguments total more than **5 000 000 characters** (about 1 GiB at 200 bytes per character) is
rejected with an `ArgumentException` ("Input too large ...") before anything is built. The limit counts all text arguments of the call
together (text + pattern, text + query, the sum of the k texts plus one separator each).

Change it when starting the server (a value that is not a positive integer stops the server with a clear message):

```bash
dotnet run --project SuffixTree.Mcp.Core -- --SuffixTree:MaxTextChars=20000000
SuffixTree__MaxTextChars=20000000 dotnet run --project SuffixTree.Mcp.Core
```

`calculate_similarity`, `hamming_distance`, `edit_distance` and `count_approximate_occurrences` build no tree and are not covered.
Persistent sessions that keep an index across calls (so large inputs need not be re-sent and rebuilt) are planned in stage S8 of the
suffix-tree upgrade plan. See `SuffixTreeMcpLimits` for how the default was derived.

## Tools (18)

| Tool | Description |
|------|-------------|
| `calculate_similarity` | Calculate similarity between two DNA sequences using k-mer Jaccard index (0-100 percentage scale). |
| `count_approximate_occurrences` | Count approximate occurrences of a pattern in a sequence, allowing up to maxMismatches substitutions. |
| `edit_distance` | Calculate edit distance (Levenshtein distance) between two sequences; optional insertion/deletion/substitution costs give the weighted distance. |
| `find_longest_common_region` | Find the longest common region between two DNA sequences. |
| `find_longest_repeat` | Find the longest repeated region in a DNA sequence. |
| `hamming_distance` | Calculate Hamming distance between two sequences of equal length. |
| `suffix_tree_all_lcs` | Find every distinct longest common substring of two texts (all length ties, unlike suffix_tree_lcs which returns one), each with all 0-based… |
| `suffix_tree_all_lrs` | Find every distinct longest repeated substring of text (all length ties, unlike suffix_tree_lrs which returns one representative), each wit… |
| `suffix_tree_contains` | Check if a pattern exists in text using suffix tree. |
| `suffix_tree_count` | Count the number of occurrences of a pattern in text using suffix tree. |
| `suffix_tree_find_all` | Find all 0-based start positions (including overlapping occurrences) where a pattern occurs in text using suffix tree, in ascending order. |
| `suffix_tree_find_mems` | Find all maximal exact matches (MEMs) of length >= minLength between a reference text and a query: every left- and right-maximal match with… |
| `suffix_tree_find_mums` | Find maximal unique matches (MUMs) of length >= minLength: MEMs whose string occurs exactly once in the reference and, with uniqueness 'bot… |
| `suffix_tree_k_common_substrings` | Find every distinct longest substring that occurs in at least minSupport of the given texts (default: all of them = longest common substring… |
| `suffix_tree_lcs` | Find the longest common substring between two texts using suffix tree. |
| `suffix_tree_lrs` | Find the longest repeated substring in text using suffix tree. |
| `suffix_tree_maximal_repeats` | Find every maximal repeated pair (firstPosition < secondPosition, length >= minLength) of a text: two equal substrings that cannot be extend… |
| `suffix_tree_stats` | Get statistics about a suffix tree: node count, leaf count, max depth, and text length. |
