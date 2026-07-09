# Wiki Log

Append-only chronological record of operations on the wiki. Each entry begins with `## [YYYY-MM-DD] <op> | <description>` so it's parseable with `grep "^## \[" log.md | tail -N`.

Operations:
- `ingest` — a source was processed into the wiki.
- `query` — a question was answered against the wiki (typically only logged when the answer was filed back as synthesis).
- `lint` — a health check was run.
- `schema` — the schema was modified.
- `shard` — an index was sharded.

---

## [2026-07-09] ingest | README.md → readme (source) + 5 concepts + 1 gotcha
   Established hub pages from the project front page: three-front-doors, skill-layer,
   layered-architecture, scientific-rigor, validation-and-testing, research-grade-limitations.
   graph: +7 nodes, +4 typed edges

## [2026-07-09] ingest | ALGORITHMS_CHECKLIST_V2.md → algorithms-checklist-v2 (source) + 2 concepts
   Created test-unit-registry, definition-of-done. Updated validation-and-testing and
   research-grade-limitations (campaign-added pending re-validation, unverified complexity claims).
   graph: +3 nodes, +2 typed edges

## [2026-07-09] ingest | docs/ADVANCED_TESTING_CHECKLIST.md → advanced-testing-checklist (source)
   New source page (technique effectiveness/gap analysis, 10 techniques, P0–P3, 2026-03-19 baseline).
   Updated validation-and-testing (coverage-gap paragraph + typed edge). Flagged internal 79-vs-86
   unit discrepancy and that only architecture testing is complete. No new concepts (elaborates existing).
   graph: +1 node, +1 typed edge

## [2026-07-09] ingest | docs/MCP-Methods-Audit.md → mcp-methods-audit (source)
   New source page: 277 public static methods / 54 classes census of the C# API surface, with
   counting rules (overloads, extensions, SAM-flag one-liners, StatisticsHelper). Updated
   three-front-doors (quantified the C# API door). Flagged 277-methods-vs-427-MCP-tools and
   54-vs-57-class-denominator discrepancies as data points, not contradictions.
   graph: no typed edges (flat inventory; mentions auto-derived); +1 node

## [2026-07-09] ingest | docs/mcp-checklist.md → mcp-checklist (source)
   New source page for the MCP Implementation Checklist v4. Flagged prominently as SUPERSEDED
   (2026-07-01): a 12-server/241-tool plan never built; live status is docs/mcp/MCP_STATUS.md (11
   real servers). Captured DoD gates, 1000-5999 error catalog, two-tests-per-tool (Schema+Binding,
   no business asserts), MethodId/HasDocs/DocRef traceability, G1-G5 gates, and the SuffixTree.Mcp→
   Seqeron.Mcp rename. Added 241-vs-277-vs-427 tool-count reconciliation; cross-linked from
   mcp-methods-audit. No typed edges (supersession target docs/mcp/MCP_STATUS.md is not yet a wiki
   node; mentions auto-derived).

## [2026-07-09] ingest | docs/mcp-plan.md → mcp-plan (source)
   New source page for the MCP Implementation Plan v4 — the sibling *design doc* to mcp-checklist.
   Flagged SUPERSEDED (2026-07-01): 12-server/241-tool design never built (real: 11 servers/427
   tools, docs/mcp/MCP_STATUS.md). Distinctive content vs the checklist: v3→v4 context-budget
   rationale (5-7%/24%/2-8%), full 241-tool inventory across all 12 servers, and sections 6-8
   standards the banner declares still-valid (naming, JSON-Schema 2020-12, error ranges 1000-5999,
   SemVer, 2-tests-per-tool, .mcp.json+.md doc contract). Reused the 241/277/427 reconciliation.
   Cross-linked bidirectionally with mcp-checklist. No typed edges (superseded near-sibling source;
   ontology has no fitting source→source predicate, and count deltas are reconciled not contradictory;
   mentions auto-derived).

## [2026-07-09] ingest | docs/mcp-prompt.md → mcp-prompt (source)
   New source page for the CURRENT one-tool-per-session MCP-completion subagent prompt (the live
   successor to the superseded mcp-plan/mcp-checklist). Captured: docs/mcp/MCP_STATUS.md as authoritative
   B/T/D ledger; the shipped 11-server decomposition named concretely (server→project→tools file) —
   first source to enumerate it (no Variants/Assembly/Epigenetics/Structure servers; Analysis+Annotation+
   MolTools consolidate them; Core still under SuffixTree.Mcp.Core); Sequence/Parsers/Core gold standard;
   the 3-part tool-wrapper DoD (binding attribute+record return, ≥2 NUnit tests, .mcp.json+.md docs);
   execution flow + full-green-gate + stale-bin/obj caveat. FLAGGED CONTRADICTION: the prompt's DoD
   requires evidence-based Binding tests asserting exact documented values, reversing the
   no-business-asserts policy in mcp-checklist/mcp-plan (annotated both pages surgically). Cross-linked
   both ways with mcp-checklist and mcp-plan.
   graph: +1 node, +2 typed edges (contradicts → mcp-checklist, mcp-plan)

## [2026-07-09] ingest | docs/sonar-gate-plan.md → sonar-gate-plan (source) + 1 concept
   New source page for the Sonar gate ratchet tracker (66/66 SonarAnalyzer rules → blocking or
   silenced-with-justification; green under TreatWarningsAsErrors; 14 assemblies / 20,266 core
   tests). Created the build-quality-gate concept (static-analysis gate + warnings-as-errors,
   fix-vs-silence ratchet, review-not-blind-fix on S1244/S125). Linked it from validation-and-testing
   (added docs/sonar-gate-plan.md to that page's sources). Flagged the doc's internal staleness: a
   "remaining 31 rules" planning section survives alongside the 66/66 completion banner (Log rows are
   ground truth). Captured the S4456 fail-fast behaviour change and the pre-existing flaky FsCheck
   properties.
   graph: +2 nodes, +1 typed edge (build-quality-gate relates_to validation-and-testing)

## [2026-07-09] ingest | docs/Evidence/ALIGN-GLOBAL-001-Evidence.md → align-global-001-evidence (source) + 2 concepts
   First per-algorithm Evidence file (of ~213). Created the shared hub concept
   algorithm-validation-evidence (templated 5-part structure: header/online-sources/dataset/
   deviations/references) so future evidence ingests link in rather than duplicate. Created the
   genuinely-distinct algorithm concept global-alignment-needleman-wunsch (linear-gap recurrence,
   O(nm), traceback, GapExtend=d / GapOpen-unused, affine-as-extension). Concise source page for
   the ALIGN-GLOBAL-001 artifact (Wikipedia sources, GCATGCG/GATTACA example, score 0). Linked
   the evidence hub from test-unit-registry. No contradictions; deviations = None.
   graph: +3 nodes, +3 typed edges

## [2026-07-09] ingest | docs/Evidence/ALIGN-MULTI-001-Evidence.md → align-multi-001-evidence (source) + 1 concept
   Second per-algorithm Evidence file. Created the genuinely-distinct concept
   multiple-sequence-alignment (star `MultipleAlign`, iterative-refinement `MultipleAlignIterative`
   = MUSCLE Stage 3, consistency `MultipleAlignConsistency` = T-Coffee; SP-score objective, majority
   consensus, invariants, NP-completeness). Concise source page capturing the unusually rich file
   (main doc + MUSCLE and T-Coffee addenda, 2026-06-23) and its design-choice deviations (cosine-
   similarity center pick, gap-gap=0, gap-votes/tie-to-nucleotide consensus, UPGMA vs NJ, gap-0
   progressive DP). Linked the new source + concept into the algorithm-validation-evidence hub and
   added the ALIGN-MULTI source to that hub's frontmatter. No contradictions; all deviations flagged
   as API/structure choices preserving the verified properties.
   graph: +2 nodes, +2 typed edges (relates_to test-unit-registry, alternative_to global-alignment-needleman-wunsch)

## [2026-07-09] ingest | docs/Evidence/ALIGN-SEMI-001-Evidence.md → align-semi-001-evidence (source) + 1 concept
   Third per-algorithm Evidence file. Created the genuinely-distinct concept
   semi-global-alignment-fitting (ends-free "glocal" hybrid; fitting/query-in-reference variant
   = Rosalind SIMS; NW recurrence with no zero floor, first row = 0 free reference start gaps,
   first column = d·i, traceback from max of last row; overlap/OAP and full-semiglobal/SMGB
   noted as sibling variants; INV-1..5). Concise source page for the ALIGN-SEMI-001 artifact
   (Wikipedia + Rosalind SIMS/SMGB + Brudno 2003 glocal sources, corner cases, fitting-variant
   design choice). Linked new source + concept into the algorithm-validation-evidence hub and
   added ALIGN-SEMI to that hub's frontmatter; added a reciprocal nav link from
   global-alignment-needleman-wunsch. No contradictions; deviation = deliberate fitting-variant
   selection + standard .NET null contract.
   graph: +2 nodes, +2 typed edges (relates_to test-unit-registry, alternative_to global-alignment-needleman-wunsch)

## [2026-07-09] ingest | docs/Evidence/ALIGN-STATS-001-Evidence.md → align-stats-001-evidence (source) + 1 concept
   Fourth per-algorithm Evidence file. Created the genuinely-distinct concept
   alignment-statistics (post-alignment metric layer, not an aligner): percent
   identity/similarity/gaps under the EMBOSS/BLAST convention (count / Length × 100,
   denominator includes gap columns; Similarity = identical OR positively-scoring columns, so
   Similarity ≥ Identity; "positive substitution score ⇒ similar"); DNA SimpleDna ⇒ Similarity
   = Identity vs Mismatch=+1 ⇒ Similarity > Identity; srspair three-line markup (|/:/space, the
   graded `.` tier unreachable → rendering-only). Concise source page for the ALIGN-STATS-001
   artifact (EMBOSS needle/AlignFormats + BLAST NBK1734 + pseqsid sources, the 149-column
   HBA/HBB worked example as a formula cross-check, two hand-built DNA datasets, empty/null/
   lineWidth contracts). Linked new source + concept into the algorithm-validation-evidence hub
   and added ALIGN-STATS to that hub's frontmatter; added a reciprocal nav link from
   global-alignment-needleman-wunsch. No contradictions; only deviation is the rendering-only
   `.`-tier collapse, non-correctness-affecting.
   graph: +2 nodes, +2 typed edges (relates_to test-unit-registry, relates_to global-alignment-needleman-wunsch)

## [2026-07-09] ingest | docs/Evidence/ANNOT-CODING-001-Evidence.md → annot-coding-001-evidence (source) + 1 concept
   Fifth per-algorithm Evidence file and the first from the annotation domain. Created the
   genuinely-distinct concept coding-potential-hexamer-score (CPAT hexamer usage-bias coding
   potential, Wang et al. 2013): score = mean of ln(coding[k]/noncoding[k]) over in-frame
   hexamers (wordSize 6, stepSize 3, natural log, codon-boundary words via word_generator);
   sign convention positive=coding / negative=noncoding; pseudo-score branches (coding-only
   +1, noncoding-only -1), both-zero-in-both hexamer skipped-not-counted (verified verbatim
   vs canonical liguowang/cpat + WGLab/lncScore 2026-06-15), short-seq and missing-key skips;
   unit-agnostic tables (counts vs proportions differ by additive constant). Concise source
   page for the ANNOT-CODING-001 artifact (CPAT paper + FrameKmer.py + EMBOSS tcode/Fickett
   sources, worked example 0.34657359027997264, two assumption records). Linked new source +
   concept into the algorithm-validation-evidence hub and added ANNOT-CODING to that hub's
   frontmatter. Contradictions: none between sources; only deviation is the C# port's -1->0
   no-scorable-hexamer sentinel (affects meaningless empty scores only). Follow-up: Fickett
   TESTCODE recorded as a related not-implemented alternative if it is ever ingested.
   graph: +2 nodes, +1 typed edges (relates_to test-unit-registry)

## [2026-07-09] ingest | docs/Evidence/ANNOT-CODONUSAGE-001-Evidence.md → annot-codonusage-001-evidence (source) + 1 concept
   Sixth per-algorithm Evidence file; first of the large codon-usage family (CODON-CAI/ENC/OPT/
   RARE/RSCU/STATS/USAGE, SEQ-CODON-FREQ, TRANS-CODON still to come). Created the
   genuinely-distinct concept relative-synonymous-codon-usage (RSCU, Sharp & Li 1986): per-codon
   codon-usage-bias normalization RSCU = n_i·x_{i,j}/Σx over a synonymous family; 1.0=no bias,
   >1 preferred / <1 under-represented, bounded [0,n_i], Σ-over-family=n_i invariant; counts
   pooled across all reference sequences, sense codons only (forward_table), single-codon Met/Trp
   always 1.0, unobserved family → 0.0; Standard NCBI table 1 default. Deliberately positioned as
   the base anchor of the codon family so future codon ingests link in; distinguished from CAI's
   0.5 pseudocount (Sharp & Li 1987, CAI-only, NOT applied to plain RSCU). Concise source page for
   the ANNOT-CODONUSAGE-001 artifact (LIRMM formula page + PMC2528880 + Sharp & Li 1986 primary +
   CodonU internal_comp.py::rscu + NCBI table 1 sources; Leu CTTCTTCTGTTA → 3/1.5/1.5/0/0/0,
   uniform Phe → 1.0, Met → 1.0 datasets; two API-default assumptions). Linked new source + concept
   into the algorithm-validation-evidence hub and added ANNOT-CODONUSAGE to that hub's frontmatter.
   Contradictions: none — LIRMM formula, PMC2528880 definition, and CodonU code are algebraically
   identical. Follow-up: when CODON-RSCU-001 (an apparent RSCU duplicate) and the rest of the codon
   family are ingested, share this concept rather than duplicating; CAI/ENC/etc. may each warrant
   their own concept.
   graph: +2 nodes, +1 typed edges (relates_to test-unit-registry)

## [2026-07-09] ingest | docs/Evidence/ANNOT-REPEAT-001-Evidence.md → annot-repeat-001-evidence (source) + 1 concept
   Seventh per-algorithm Evidence file. Created the genuinely-distinct concept
   repetitive-element-detection, deliberately scoped as the shared anchor for the whole
   repeats/tandem family (GENOMIC-REPEAT, GENOMIC-TANDEM, microsatellite/STR, low-complexity)
   so future repeat ingests link in rather than re-deriving definitions. Covers the three
   sub-problems: tandem repeats (head-to-tail, ≥2 copies, STR 1-6bp / minisatellite 10-60bp,
   primitive-shortest-period rule), inverted repeats (IUPACpal grammar W W̄ᴿ / W G W̄ᴿ, imperfect
   δ_H ≤ k, zero-gap = even-length palindrome), and RepeatMasker-class assignment (SINE/LINE/LTR/
   DNA/Satellite/Simple_repeat/Low_complexity/Small RNA/Unknown). Concise source page for the
   ANNOT-REPEAT-001 artifact (Wikipedia Tandem/Inverted + IUPACpal Hampson 2021 PMC7866733 +
   RepeatMasker sources; ATTCGATTCGATTCG/GAATTC/TTACGAAAAAACGTAA datasets; six MUST tests).
   Captured the one assumption: ClassifyRepeat matches by exact-substring containment (element ⊆
   query, longest match, one-directional) with motif-size Simple_repeat fallback, NOT Smith-
   Waterman-Gotoh homology against a curated Repbase library — a Framework/Simplified limitation,
   vocabulary source-backed. Linked new source + concept into the algorithm-validation-evidence hub
   and added ANNOT-REPEAT to that hub's frontmatter. Contradictions: none (Wikipedia & IUPACpal IR
   definitions are the same grammar; RepeatMasker class list is shared vocabulary). Follow-up: when
   GENOMIC-REPEAT/GENOMIC-TANDEM and other repeat-family units are ingested, share
   repetitive-element-detection rather than duplicating.
   graph: +2 nodes, +1 typed edges (relates_to test-unit-registry)

## [2026-07-09] ingest | docs/Evidence/ASSEMBLY-CONSENSUS-001-Evidence.md → assembly-consensus-001-evidence (source) + 1 concept
   Eighth per-algorithm Evidence file; first of the large Assembly family (DBG/OLC/SCAFFOLD/COVER/
   STATS/TRIM/CONSENSUS/CORRECT/MERGE still to come). Rejected a broad "genome-assembly" hub as too
   vague/heterogeneous for a useful anchor; instead created the genuinely-distinct concept
   consensus-sequence (column-wise majority/threshold consensus — the C of Overlap-Layout-Consensus
   and the same operation as the MSA consensus step). Decision rule traced verbatim to Biopython
   dumb_consensus: tally non-gap residues only, emit iff unique max AND max_size/num_atoms >= threshold
   (strict >=) else ambiguous; tie→ambiguous (not arbitrary pick); all-gap column→ambiguous with no
   div-by-zero (short-circuit and); consensus length = longest read (ragged handled). EMBOSS cons
   plurality + Wikipedia definition/IUPAC corroborate. Two parameterized presentation-only default
   assumptions: ambiguous symbol N-not-X (DNA/RNA IUPAC), default threshold 0.5-not-Biopython-0.7
   (simple-majority); neither alters the source-backed rule. Concise source page for the artifact
   (three sources, four datasets, nine MUST tests). Linked new source + concept into the
   algorithm-validation-evidence hub and added ASSEMBLY-CONSENSUS to that hub's frontmatter.
   Contradictions: none. Follow-up: future Assembly-family and MSA-consensus units share
   consensus-sequence rather than duplicating; graph/OLC/scaffold units likely warrant their own
   distinct concepts.
   graph: +2 nodes, +2 typed edges (relates_to test-unit-registry, relates_to multiple-sequence-alignment)

## [2026-07-09] ingest | docs/Evidence/ASSEMBLY-CORRECT-001-Evidence.md → assembly-correct-001-evidence (source) + 1 concept
   Ninth per-algorithm Evidence file; second of the Assembly family (after CONSENSUS). Created the
   genuinely-distinct concept kmer-spectrum-error-correction — the anchor for the assembly CORRECT
   family (distinct from consensus-sequence). Two-sided k-mer-spectrum read error correction traced
   verbatim to Musket (Liu 2013) + Quake (Kelley 2010): trusted k-mer = multiplicity > coverage
   cut-off, base trusted if covered by any trusted k-mer; find the unique alternative base making
   ALL k-mers covering position i trusted (evaluate leftmost AND rightmost covering k-mer), ≤1
   substitution per k-mer; >1 valid alternative → unchanged (ambiguity), no correcting set →
   unchanged; substitution-only so read count + per-read length preserved. Song & Florea 2018
   corroborate (solid/weak k-mers). Concise source page for the artifact (three rank-1 sources, the
   k=3/cut-off=2 single-substitution ACGTACGT worked example, the k=1 ambiguity example, six MUST
   tests). Linked new source + concept into the algorithm-validation-evidence hub and added
   ASSEMBLY-CORRECT to that hub's frontmatter. One assumption: kmerSize=15/minKmerFrequency=2 fixed
   defaults (vs reference auto-cut-off-from-histogram-valley) are non-behavioral — every behavioral
   test passes k and cut-off explicitly. Contradictions: none (all three sources describe the same
   trusted/untrusted two-sided model). Follow-up: remaining Assembly-family units (DBG/OLC/SCAFFOLD/
   COVER/STATS/TRIM/MERGE) likely warrant their own distinct concepts.
   graph: +2 nodes, +1 typed edge (relates_to test-unit-registry)

## [2026-07-09] ingest | docs/Evidence/ASSEMBLY-COVER-001-Evidence.md → assembly-cover-001-evidence (source) + 1 concept
   Tenth per-algorithm Evidence file; third of the Assembly family (after CONSENSUS, CORRECT).
   Created the genuinely-distinct concept coverage-depth-calculation — the anchor for the assembly
   COVER family. Per-base sequencing depth = count of placed reads spanning each reference position
   (exact, model-free); average depth = Σdepth/G = Lander-Waterman C=LN/G; breadth = (#depth≥1)/G =
   1−e^−c. Boundary-clip at reference end + all-zero/empty-input rules; hand-built ACGTTGCAAT oracle
   (depth [1,1,1,2,2,2,2,2,1,1], avg 1.5, breadth 1.0); Lander-Waterman Poisson (P(uncovered)=e^−c,
   1×→0.37, 5×→0.0067) captured explicitly as a property/derivation check only — the per-base array
   is exact regardless of uniformity. Concise source page for the artifact (Illumina rank-2 + Daniel
   Cook + Metagenomics Wiki rank-3 + Daley PMC7398442 rank-1 + Lander-Waterman 1988 primary sources,
   two datasets, seven MUST/SHOULD/COULD tests). Linked new source + concept into the
   algorithm-validation-evidence hub and added ASSEMBLY-COVER to that hub's frontmatter. One
   assumption: read-placement model (ungapped minOverlap best-match FindBestAlignment) is
   implementation-level and out of scope — tests use unambiguous exact-match reads to isolate the
   source-defined counting rule. Contradictions: none (all sources give the same depth/average/breadth
   definitions). Follow-up: remaining Assembly-family units (DBG/OLC/SCAFFOLD/STATS/TRIM/MERGE) likely
   warrant their own distinct concepts.
   graph: +2 nodes, +1 typed edge (relates_to test-unit-registry)

## [2026-07-09] ingest | docs/Evidence/ASSEMBLY-DBG-001-Evidence.md → assembly-dbg-001-evidence (source) + 1 concept
   Eleventh per-algorithm Evidence file; fourth of the Assembly family (after CONSENSUS, CORRECT,
   COVER). Created the genuinely-distinct concept de-bruijn-graph-assembly — the anchor for the
   assembly DBG family (BuildDeBruijnGraph + AssembleDeBruijn). Graph construction traced verbatim to
   Langmead's JHU DBG notes: distinct (k-1)-mers are nodes, each k-mer is one directed prefix→suffix
   edge, repeated k-mers make a directed multigraph; chop bound range(0,len-(k-1)) ⇒ reads < k yield
   no k-mers. Reconstruction as an Eulerian walk under Jones & Pevzner Theorems 8.1 (cycle iff all
   balanced) / 8.2 (path iff ≤2 semi-balanced), O(|E|) Hierholzer, spelled as path[0] + last char of
   each subsequent node; Compeau-Pevzner-Tesler 2011 supply the assembly application (Eulerian-path
   tractable vs NP-complete Hamiltonian/overlap). Unique-walk oracles (AAABBBA k=3 with full node/edge
   set, a_long_long_long_time k=5, to_every… k=4-correct/k=3-wrong turn-repeat, ATGGCGTGCA k=4) plus
   the AAABBBBA multiedge case; failure modes (repeat≥k-1 → multiple walks, gap → disconnected/multi-
   contig, extra copy/error → non-Eulerian, Superwalk NP-hard). Concise source page for the artifact.
   Linked new source + concept into the algorithm-validation-evidence hub and added ASSEMBLY-DBG to
   that hub's frontmatter. Three assumptions: walk-selection unspecified (exact asserts on unique-walk
   inputs only; non-unique checked on invariants/branch structure), empty/null → empty AssemblyResult
   (mirrors OLC), reads < k contribute no k-mers. Contradictions: none — Langmead cites the same J&P
   Euler theorems Compeau builds on; Compeau 2011 PDFs are image-only so cited for metadata only.
   Follow-up: OLC (the alternative fragment-assembly formulation) warrants its own concept when
   ingested; remaining Assembly units (SCAFFOLD/STATS/TRIM/MERGE) likely warrant distinct concepts.
   graph: +2 nodes, +1 typed edge (relates_to test-unit-registry)

## [2026-07-09] ingest | docs/Evidence/ASSEMBLY-MERGE-001-Evidence.md → assembly-merge-001-evidence (source) + 1 concept
   Twelfth per-algorithm Evidence file; fifth of the Assembly family (after CONSENSUS, CORRECT,
   COVER, DBG). Created the genuinely-distinct concept contig-merge-overlap-collapse — the anchor
   for the assembly MERGE family: the suffix–prefix overlap collapse primitive
   MergeContigs(contig1, contig2, overlapLength) behind greedy shortest-common-superstring and the
   OLC layout step. Overlap traced verbatim to Langmead's JHU SCS/OLC notes + MIT 7.91J Lecture 6:
   overlap = length-l suffix of X exactly matching a length-l prefix of Y (l ≤ min(|X|,|Y|)),
   suffixPrefixMatch returns the longest such match else 0, collapse keeps one copy so
   |merge| = |c1|+|c2|−l; overlap 0 → plain concatenation X+Y. Published oracles BAA+AAB(ov2)→BAAB,
   {AAA,AAB,ABB,BBB,BBA} chain→AAABBBA (len 7), BAA+AAB(ov0)→BAAAAB. Two API-contract assumptions
   (caller-supplied overlap length trusted not re-verified — verification is FindOverlap's job;
   out-of-range overlap ≤0 or >min → concatenation), both derived directly from the source facts,
   neither a correctness/scoring parameter. Concise source page for the artifact (Langmead SCS +
   Langmead OLC + MIT 7.91J rank-1 sources, Compeau 2011 background-only, three oracles, MUST/SHOULD/
   COULD tests). Linked new source + concept into the algorithm-validation-evidence hub and added
   ASSEMBLY-MERGE to that hub's frontmatter; added a reciprocal nav link from de-bruijn-graph-assembly
   (MERGE is the overlap-based sibling of the DBG k-mer/Eulerian formulation). Contradictions: none —
   the three sources give the identical suffix-of-X/prefix-of-Y overlap definition and corroborate one
   another. Follow-up: an end-to-end OLC concept (and remaining Assembly units SCAFFOLD/STATS/TRIM)
   warrant their own pages when ingested; FindOverlap/FindAllOverlaps (the overlap-discovery side)
   would share contig-merge-overlap-collapse.
   graph: +2 nodes, +2 typed edges (relates_to test-unit-registry, relates_to de-bruijn-graph-assembly)

## [2026-07-09] ingest | docs/Evidence/ASSEMBLY-OLC-001-Evidence.md → assembly-olc-001-evidence (source) + 1 concept
   Thirteenth per-algorithm Evidence file; sixth of the Assembly family (after CONSENSUS, CORRECT,
   COVER, DBG, MERGE). Created the genuinely-distinct concept overlap-layout-consensus-assembly —
   the anchor for the assembly OLC family and the second of the two canonical fragment-assembly
   paradigms (FindAllOverlaps + AssembleOLC). Three stages traced verbatim to Compeau, Pevzner &
   Tesler 2011 + Langmead OLC/SCS notes: Overlap (read=node overlap graph, directed edge A→B on the
   longest suffix-of-A/prefix-of-B ≥ threshold, report only longest per pair), Layout (exact = a
   Hamiltonian path = NP-complete → heuristic transitive reduction + non-branching-stretch contigs),
   Consensus (majority vote per column). Complexity suffix-tree O(N+a) vs all-pairs DP O(N²).
   Published oracles: GTACGTACGAT 6-mers minOverlap4 → exactly 12 directed edges (lengths 4/5,
   re-derived), 5-overlap tiling → single AAAAACCCCCGGGGGTTTTT, CTCTAGGCC/TAGGCCCTC l=3 → overlap 6.
   Failure modes: NP-complete layout, repeats>read-length split contigs, error dead-end subgraphs,
   greedy-SCS suboptimal, sub-resolution repeats collapse. Two assumptions: exact-match identity 1.0
   for canonical numeric cases (minIdentity generalizes; separate threshold test 0.875 accepted@0.85/
   rejected@0.95), empty read set → empty AssemblyResult. Concise source page for the artifact. Linked
   new source + concept into the algorithm-validation-evidence hub and added ASSEMBLY-OLC to that hub's
   frontmatter; added reciprocal nav links from de-bruijn-graph-assembly (fulfilling its flagged OLC
   follow-up) and contig-merge-overlap-collapse. Modeled OLC as alternative_to de-bruijn-graph-assembly
   (Hamiltonian/overlap-graph vs Eulerian/k-mer, the contrast Compeau 2011 draws explicitly).
   Contradictions: none — Compeau 2011 and both Langmead notes give the identical overlap-graph/
   Hamiltonian-path/three-stage account; re-derived numeric oracles match the source slides. Follow-up:
   remaining Assembly units (SCAFFOLD/STATS/TRIM) warrant their own pages; FindOverlap/FindAllOverlaps
   overlap-discovery shares both this concept and contig-merge-overlap-collapse.
   graph: +2 nodes, +2 typed edges (relates_to test-unit-registry, alternative_to de-bruijn-graph-assembly)

## [2026-07-09] ingest | docs/Evidence/ASSEMBLY-SCAFFOLD-001-Evidence.md → assembly-scaffold-001-evidence (source) + 1 concept
   Fourteenth per-algorithm Evidence file; seventh of the Assembly family (after CONSENSUS, CORRECT,
   COVER, DBG, MERGE, OLC). Created the genuinely-distinct concept scaffolding — the anchor for the
   assembly SCAFFOLD family, deliberately positioned as a *downstream* step orthogonal to the
   overlap-vs-k-mer DBG/OLC contrast (it lays finished ordered contigs onto a coordinate frame with
   sized gaps rather than reconstructing sequence). Construction rule traced verbatim to Jackman et
   al. ABySS 2.0 (Genome Research 2017): scaffold = ordered path contigs concatenated interspersed
   with runs of `N` whose length = the (upstream ML-estimated) inter-contig distance; positive gap
   g → exactly g fill chars, scaffold length = Σ|contig|+Σgap; each contig in ≤1 scaffold, unlinked →
   singleton; fill char parameterized (source fixes default `N`). Non-positive (zero/negative) gap →
   AGP unknown-size default 100 N: NCBI AGP v2.1 ("gap lengths must be positive ... use U and 100
   ... GenBank/EMBL/DDBJ standard for unknown-size gaps") supplies the source-backed 100 constant,
   Sahlin et al. 2012 confirm the negative-gap case is frequent (de Bruijn one-k-mer overlap), ABySS
   says negative = overlap → merge if found. Oracles ACGTNNNTTGGNNCCAA (len 17, 1 scaffold) and
   AAAA+100N+TTTT (len 108). One assumption, a scoping decision not an invented value: unresolved-
   overlap placeholder falls back to the AGP unknown-gap length 100 (this unit does no overlap
   resolution) — the 100 is source-backed, only the fall-back-rather-than-resolve choice is assumed.
   Concise source page for the artifact (ABySS 2.0 + AGP v2.1 + Sahlin 2012 + Bambus sources, two
   oracles, MUST/SHOULD/COULD coverage). Linked new source + concept into the algorithm-validation-
   evidence hub and added ASSEMBLY-SCAFFOLD to that hub's frontmatter; added a reciprocal nav link
   from contig-merge-overlap-collapse (scaffolding hands off to the suffix–prefix merge primitive on
   the negative-gap = overlap case). Contradictions: none — ABySS/AGP/Sahlin/Bambus give the same
   ordered-contigs + sized-`N`-gap model; the AGP 100-N default and the ABySS negative-gap=overlap
   rule are complementary. Follow-up: remaining Assembly units (STATS/TRIM) warrant their own pages;
   an overlap-resolving scaffolder (if ever built) would compose scaffolding with contig-merge.
   graph: +2 nodes, +2 typed edges (relates_to test-unit-registry, relates_to contig-merge-overlap-collapse)

## [2026-07-09] ingest | docs/Evidence/ASSEMBLY-STATS-001-Evidence.md → assembly-stats-001-evidence (source) + 1 concept
   Fifteenth per-algorithm Evidence file; eighth of the Assembly family (after CONSENSUS, CORRECT,
   COVER, DBG, MERGE, OLC, SCAFFOLD). Created the genuinely-distinct concept assembly-statistics —
   the anchor for the assembly STATS family: the read-only contiguity/QC summary metrics over a set
   of contig lengths (N50/L50/Nx/Lx/N90/L90/auN + totals/largest/smallest/GC/N-run gaps), downstream
   of and orthogonal to the DBG/OLC/merge build steps and the scaffolding layout step. Definitions
   traced verbatim to Miller, Koren & Sutton 2010 (N50 = smallest of the fewest largest contigs whose
   combined length is "at least 50%") + Wikipedia worked examples + QUAST `N50.py` (`NG50_and_LG50`
   stop test `s <= limit` = inclusive cumulative ≥ threshold; `au_metric` = ΣL²/ΣL) + Heng Li 2020
   (Nx covers x%, auN = area under Nx curve = ΣLᵢ²/ΣLⱼ). Nx is a length, Lx a count; N90 ≤ N50 /
   L90 ≥ L50 monotonicity; boundary inclusive (≥). Published oracles: Assembly A {80,70,50,40,30,20}
   total 290 → N50 70/L50 2/N90 30/L90 5/auN≈57.586, Assembly B (+10,+5) total 305 → N50 50/L50 3,
   auN check {100,80,60,40,20} → 22000/300 = 73.333…/N50 80/L50 2. Two assumptions outside the cited
   contract: empty input → all-zero AssemblyStatistics (QUAST returns None) — an API-shape choice
   changing no defined value (§6.1); and CalculateStatistics.MedianLength reports the upper median
   (lengths[count/2]), an auxiliary field the N50 literature doesn't define, tested-as-implemented and
   flagged not-source-derived. Concise source page for the artifact (Miller 2010 rank-1 + Wikipedia
   rank-4 + QUAST rank-3 + Heng Li rank-3 sources, three datasets, MUST/SHOULD/COULD coverage). Linked
   new source + concept into the algorithm-validation-evidence hub and added ASSEMBLY-STATS to that
   hub's frontmatter. Contradictions: none — Miller/Wikipedia/QUAST/Heng Li give identical
   largest-first inclusive-≥ definitions and QUAST au_metric matches Heng Li's ΣL²/ΣL exactly.
   Follow-up: remaining Assembly unit (TRIM) warrants its own page when ingested.
   graph: +2 nodes, +1 typed edge (relates_to test-unit-registry)
