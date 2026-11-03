# Suffix-tree competitive matrix — SPEC (v1, 2026-11-02)

Goal: a competitor × capability score matrix for a .NET suffix-tree library ("Seqeron.SuffixTree", genomic use),
later used to plan upgrades. YOU score ONLY the competitors assigned to you. Do NOT score Seqeron.

## Hard rules
1. Every cell needs a primary-source URL (official docs, README/source file at a pinned tag, paper/arXiv, release page).
   No URL = no verified score. Never invent a URL; open it (WebFetch) before citing it. A search snippet is not a source.
2. If after real effort a cell cannot be verified, score it OPTIMISTICALLY (the level you think it most plausibly has,
   rounded up) and flag it `U`. Verified cells get flag `V`. (Optimistic-for-competitor keeps the final claim conservative.)
3. Capability absent from the tool = 0 (that is what a user experiences). "Needs your own code on top of primitives" = 1.
4. Privacy rule (binding for every subagent): NEVER put any e-mail, name, account id or other identity of the owner into any
   header, User-Agent, URL, query or payload when fetching anything. Use the tool defaults / a generic UA only.
5. Date: today is 2026-11-02. Record the version / release date you scored. Prefer the latest stable release.
6. Do not run anything that needs credentials. Do not clone huge repos; fetch single files (raw) when possible.
7. Honest reporting: if a site is blocked, say which and what you scored instead. Do not paper over gaps.

## Scale (same for every cell)
0 absent / impossible without writing the algorithm yourself
1 partial or workaround: reachable via low-level primitives with substantial user code, OR works only on toy sizes,
  OR a documented major limitation
2 supported, basic: public API provides it, works as documented, some limits (size, alphabet, single mode)
3 full: complete, documented, usable at genome scale (≥100 Mbp) or the equivalent for that capability, no major limitation
4 best-in-class: the reference/state-of-the-art way to do this among the compared tools (the published best algorithm,
  complete, and documented as such), or — for measurable capabilities — measured/published within 1.5× of the best tool

Measurable capabilities (C10, C11, C12, C13, C15): give the number you found (units!) in the evidence column.
Where a head-to-head number is not published, score by the published/theoretical class and flag `U`.
(A separate benchmark pass will measure some tools on one machine and may override these cells.)

## Capabilities (id, weight for a genomic .NET library, what level 3 means)
C1  w3  Exact search: contains / count / locate all occurrences, text access (extract). L3: all three, any pattern length, on genome-size text.
C2  w2  Pattern classes: IUPAC / wildcards / character classes in the pattern. L3: first-class API, not a post-filter.
C3  w3  Approximate search on the index (k mismatches and/or k edits, with backtracking/pruning over the index). L3: documented API, both modes.
C4  w3  MEM / MUM / maximal repeats between texts (MUMmer-style), incl. maximal repeated pairs. L3: complete sets, documented, large inputs.
C5  w2  Matching statistics / streaming match anchors of a query against the index. L3: public API returning ms[] or equivalent.
C6  w2  Generalized index over several texts (k-common substrings, LCS of many, document listing/counting). L3: public API.
C7  w3  Both strands: reverse-complement aware search/indexing for DNA. L3: single call searches/indexes both strands.
C8  w2  Tree/interval navigation API: suffix links / Weiner links, string depth, parent/child, LCA, node-level visitors. L3: public, documented.
C9  w2  Analytics: LRS, distinct-substring counts, k-mer/substring statistics, supermaximal repeats. L3: ≥4 of these available.
C10 w3  Construction speed (bp/s) on DNA. L4: within 1.5× of the fastest compared tool for the same index kind.
C11 w3  Index memory (bytes per input char, in RAM, after build and during build peak). L4: within 1.5× of the smallest for same functionality.
C12 w3  Query speed (exact count/locate latency per pattern; throughput). L4: within 1.5× of the fastest.
C13 w3  Scale: > 2^32 characters / 64-bit offsets (human-genome 3 Gbp and beyond) without redesign. L3: documented working at ≥ 3 Gbp.
C14 w2  External / larger-than-RAM operation (mmap, disk-resident, semi-external build). L3: build AND query without holding index in RAM.
C15 w2  Parallel construction (multi-thread/GPU). L3: documented near-linear speed-up on ≥ 8 threads.
C16 w2  Concurrent queries: thread-safe read-only index shared by many threads. L3: documented + tested.
C17 w1  Online / dynamic: append characters (online), insert/delete texts. L3: documented dynamic updates.
C18 w2  Compression: compressed/succinct index (CSA, FM-index, r-index, run-length, pangenome-scale). L3: published ratio ≥ 4× smaller than plain suffix array on DNA.
C19 w2  Alphabet generality + DNA optimisation: arbitrary/Unicode/integer alphabets AND 2-bit/σ=4 optimised path. L3: both.
C20 w2  Persistence: documented on-disk format, instant (mmap) load, versioning, integrity check/hash. L3: all four.
C21 w2  API quality & integration for .NET-style use: typed/generic API, spans/zero-copy, error contracts; native .NET = L3+. For non-.NET tools score the best realistic .NET route (P/Invoke / CLI / port): native lib with C ABI = 2, C++-only = 1, CLI only = 1, none = 0.
C22 w1  Distribution & exposure: package-manager release (NuGet/PyPI/crates/conda), CLI, server/MCP tool surface. L3: ≥ 3 channels.
C23 w3  Verification: automated tests, oracle/reference parity, fuzzing, documented validation of algorithms. L3: tests in CI + reference comparison documented.
C24 w1  Maintenance & licence: releases within the last 24 months, issue activity, OSI licence. L3: all.

## Output (one TSV per agent, UTF-8, header row, tab-separated, no tabs/newlines inside cells)
File: analysis/results/<agent-id>.tsv
Columns: competitor_id	cap_id	score	flag(V|U)	evidence(<=220 chars, English)	url	accessed(YYYY-MM-DD)
Exactly one row per (assigned competitor × C1..C24) = 24 rows per competitor.
Second file results/<agent-id>.meta.tsv: competitor_id	name	version	release_date	licence	repo_or_home_url	scope_note(<=200 chars)

## Competitor ids
K1 sdsl-lite (CST/CSA/FM-index/wavelet trees, C++; also note sdsl-lite v3 / successors if any)
K2 SeqAn3 (FM-index, bidirectional FM, search schemes; mention SeqAn2 Index<ESA> as note only)
K3 GenomeTools (suffixerator / enhanced suffix arrays, Vmatch lineage; Vmatch as note)
K4 MUMmer4 (nucmer/mummer; MEM/MUM engine; mention sparseMEM/essaMEM/E-MEM as notes)
K5 libsais + libdivsufsort (SA / BWT / LCP / PLCP construction libraries)
K6 Go standard library index/suffixarray
K7 rust-bio (suffix array, BWT, FM-index, bidirectional) + the `fm-index` crate
K8 Python: `suffix-trees` (PyPI) + pydivsufsort
K9 .NET ecosystem survey: best suffix-tree / suffix-array / FM-index packages on NuGet and GitHub (score the strongest one; list the others in scope_note)
K10 BWA / Bowtie2 (FM-index read mappers; score their index as a reusable index)
K11 r-index family: r-index, Big-BWT, MONI, PHONI, SPUMONI (pangenome-scale run-length indexes)
K12 ropebwt3 (dynamic FM-index / BWT for large collections, Heng Li)
