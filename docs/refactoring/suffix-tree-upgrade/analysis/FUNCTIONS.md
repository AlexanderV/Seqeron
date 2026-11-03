# Functionality matrix — FUNCTIONS (v1, 2026-11-02)

Same project, rules and privacy rule as SPEC.md in this folder (read SPEC.md "Hard rules" first). This file replaces the
FUNCTIONAL part of the 0-4 matrix with a fine-grained function list. No weights: result = % coverage.

## Scale (per function, per tool)
0 absent: not offered by the tool's API / CLI (state the API index page you checked)
1 partial / workaround: only via low-level primitives plus user code; or limited (one strand only, tiny n, one mode,
  CLI-only when an API is expected, experimental/undocumented)
2 full: documented, directly callable (API method or CLI command/option), works as the function name says
Evidence: for 1 and 2 give the documentation URL (anchor to the page that documents it); for 0 give the URL of the API
index / manual you checked. Flag V = you read the API reference / manual; U = you only saw an overview or paper.
Best-of rule for families (K8, K9, K11): score the best member per function and name it in the evidence.

## Functions
A Construction & structure variants
F01 Build index from one in-memory text
F02 Build one index over several texts (collection / generalized)
F03 Build from a file or stream without loading the whole text first
F04 Append characters to an existing index (online)
F05 Delete / replace text in an existing index (dynamic)
F06 Merge two indexes
F07 Expose suffix array SA[i] / inverse SA
F08 Expose LCP array / PLCP
F09 Expose BWT / LF-mapping / psi
F10 Expose child table / LCP-interval tree (enhanced suffix array tree)
F11 Sampled / compressed locate (SA sampling trade-off)
B Exact search
F12 Existence test
F13 Count occurrences
F14 Locate all occurrences
F15 Locate with limit / lazy iterator / top-k
F16 Occurrences returned in text order
F17 Reusable handle for a pattern (SA interval / node) for later steps
F18 Extend a match by one character to the left (backward step)
F19 Extend a match by one character to the right
F20 Bidirectional extension (either side, same search state)
F21 Extract substring / display text from the index
F22 Map text position -> document id + offset (collections)
F23 Document listing / document counting for a pattern
F24 Top-k / ranked document retrieval
C Pattern types
F25 IUPAC / degenerate symbols in the pattern
F26 Wildcard / don't-care positions
F27 Character classes / regex-like patterns
F28 Case-insensitive search
F29 Search both strands (reverse complement) in one call
F30 Palindromes / inverted repeats finder
D Approximate search
F31 k mismatches (Hamming) on the index
F32 k edits (Levenshtein, indels) on the index
F33 Search schemes / bidirectional pigeonhole for approximate search
F34 Approximate search tolerant of N / ambiguity in the text
F35 Seed generation for read mapping (MEM / SMEM seeds)
E Matching & comparison between texts
F36 Matching statistics array of a query against the index
F37 Maximal exact matches (MEM) query vs index, all occurrences
F38 Maximal unique matches (MUM; both / reference variants)
F39 Super-maximal exact matches (SMEM)
F40 Chaining / clustering of matches into alignments (nucmer-style)
F41 Longest common substring of two texts
F42 LCS over k texts / k-common substrings
F43 Shortest unique substring / minimal unique substrings
F44 Minimal absent words
F F Repeats & text statistics
F45 Longest repeated substring (single / all ties)
F46 Maximal repeated pairs
F47 Supermaximal repeats
F48 Tandem repeats (exact)
F49 Distinct substring count (total / by length)
F50 k-mer counting / spectrum from the index
F51 Linguistic / substring complexity profile
F52 LZ factorization / longest previous factor
F53 Mappability / uniqueness map
G Tree / interval navigation
F54 Child by character
F55 Parent
F56 Suffix link
F57 Weiner link (inverse suffix link)
F58 String depth of a node
F59 Lowest common ancestor / LCP(i,j) of two suffixes
F60 Level ancestor / ancestor at a given depth
F61 Leaf rank / leaf label (suffix start)
F62 Subtree leaf count / size
F63 Enumerate children / sibling traversal
F64 DFS / visitor traversal over the whole tree
F65 Edge / path label access
H Persistence, interoperability, operations
F66 Serialize / deserialize
F67 Memory-mapped load (no full read)
F68 Versioned on-disk format with documented layout
F69 Integrity check / hash of the index
F70 Command-line tool
F71 Server / RPC / MCP-style tool surface
F72 Index statistics / inspection (nodes, depth, bytes)
F73 Visual / graph export (DOT etc.)
F74 Documented thread-safe shared read-only use
F75 Batch query API (many patterns in one call, parallel)
I Alphabets
F76 ASCII / byte text
F77 Unicode text
F78 Integer alphabet (arbitrary symbols)
F79 2-bit packed DNA
F80 Protein alphabet
F81 Policy for N / unknown symbols / separators

## Output (UTF-8 TSV, tab-separated, no tabs/newlines inside cells)
analysis/func/<agent-id>.tsv
header: competitor_id	func_id	score(0|1|2)	flag(V|U)	evidence(<=160 chars, English)	url	accessed
Exactly one row per (assigned competitor x F01..F81) = 81 rows per competitor.
Second file func/<agent-id>.extra.tsv: competitor_id	function_name_as_the_tool_calls_it	one_line_description	url
= every documented index / search / analysis function of the tool that does NOT fit any of F01..F81 (completeness check;
ignore alignment-scoring, I/O parsing and other unrelated modules).

## Competitor ids (same as SPEC.md)
K1 sdsl-lite (use the active xxsds v3 fork for current features; note upstream differences)
K2 SeqAn3, K3 GenomeTools, K4 MUMmer4, K5 libsais+libdivsufsort, K6 Go index/suffixarray, K7 rust-bio + fm-index crate,
K8 Python (suffix-trees + pydivsufsort), K9 .NET (strongest package: BelNytheraSeiche.WaveletMatrix 1.0.6; use NetBio.Core
for MEM/MUM functions and say so), K10 BWA / Bowtie2, K11 r-index family, K12 ropebwt3
