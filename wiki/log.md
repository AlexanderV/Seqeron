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

## [2026-07-09] ingest | docs/Evidence/ASSEMBLY-TRIM-001-Evidence.md → assembly-trim-001-evidence (source) + 1 concept
   Sixteenth per-algorithm Evidence file; ninth and last of the Assembly family (after CONSENSUS,
   CORRECT, COVER, DBG, MERGE, OLC, SCAFFOLD, STATS). Created the genuinely-distinct concept
   quality-trimming-running-sum — the anchor for the assembly TRIM family, deliberately positioned as
   a read-QC *preprocessing* step (operates on one read's Phred quality string, reconstructs no
   sequence) upstream of error-correction and the DBG/OLC build steps. Running-sum core traced
   verbatim to cutadapt algorithm docs (which state the algorithm "is the same as the one used by
   BWA"): subtract the cutoff from every quality, compute partial sums from each index to the 3' end,
   cut at the argmin; "repeat for the other end" on the 5' pass. BWA `bwa_trim_read` (bwaseqio.c, Heng
   Li) gives the algebraically-equivalent argmax form (accumulate threshold−(q−33) from the 3' end,
   track argmax max_l) plus two BWA-specifics: `s<0` early break and the `BWA_MIN_RDLEN=35` hard floor
   (bwtaln.h). Phred+33 decode `q = ASCII−33` from Cock et al. 2010 (NAR, rank 1). Published oracle:
   qualities 42,40,26,27,8,7,11,4,2,3 (`KI;<)(,%#$`) @ threshold 10 → partial sums min −25 at index 4
   → first 4 bases kept (with full ASCII derivation). Failure modes: threshold<1 disables (BWA guard /
   cutoff 0 → nothing trimmed), all-high unchanged, all-low fully removed, good-base-among-bad retained
   only if no new minimum reached (cutadapt "refinement"). Two assumptions, both outside the running-sum
   optimum: both-end pass order (3'-then-5' on the surviving window — not numerically significant since
   passes act on disjoint ends), and the `minLength` post-trim filter (cutadapt `--minimum-length`,
   drops trimmed length < minLength — a documented downstream filter, not the core). Concise source page
   for the artifact. Linked new source + concept into the algorithm-validation-evidence hub and added
   ASSEMBLY-TRIM to that hub's frontmatter. Contradictions: none — cutadapt explicitly identifies its
   algorithm with BWA's; the BWA argmax of accumulated (threshold−q) is the argmin of cutadapt's partial
   sums of (q−threshold); Cock supplies the shared Phred+33 encoding. Follow-up: the Assembly family is
   now fully ingested (9/9); FASTQ-quality-parsing units, if ingested, would relate to this concept.
   graph: +2 nodes, +1 typed edge (relates_to test-unit-registry)

## [2026-07-09] ingest | docs/Evidence/CHROM-ANEU-001-Evidence.md → chrom-aneu-001-evidence (source) + 1 concept
   Seventeenth per-algorithm Evidence file; first of the Chromosome-analysis family. Created the
   genuinely-distinct concept aneuploidy-detection — the anchor for the chromosome copy-number/ploidy
   family (karyotype/centromere/arm-ratio/synteny units will get their own concepts). Two-stage
   algorithm: per-bin copy number from read depth (logRatio = log2(observedDepth/medianDepth), CN =
   round(2^logRatio × 2) clamped [0,10], ×2 rescales the ratio onto the diploid baseline so ratio
   1.0 → CN 2) then whole-chromosome classification requiring a dominant CN across ≥ minFraction
   (default 80%) of bins, returning only CN ≠ 2 (nullisomy/monosomy/trisomy/tetrasomy/pentasomy/
   "Copy number = N"); the ≥80% gate is also the mosaicism tolerance. Confidence = 1 − min(1,
   |expected − observed|) with expected = CN/2, observed = 2^logRatio; = 1.0 at every integer-CN
   ratio (S1 boundary test 0.0/0.5/1.0/1.5/2.0). Concise source page for the artifact (Wikipedia
   Aneuploidy + CNV + Griffiths 2000 + Santaguida-Amon 2015 + McCarroll-Altshuler 2007 sources;
   Down/Edwards/Patau/Turner/Klinefelter clinical oracles). Linked new source + concept into the
   algorithm-validation-evidence hub and added CHROM-ANEU to that hub's frontmatter. Two documented
   limitations (artifact §7): sex chromosomes not special-cased (X/Y scored vs CN=2 baseline, normal
   male single-X would flag monosomic — research-grade simplification, autosome-focused) and partial
   aneuploidy detected per-bin but not whole-chromosome (needs consistent CN ≥80% bins). Contradictions:
   none — Wikipedia supplies the definition + CN terminology ladder, the depth→CN model and confidence
   formula are implementation definitions the sources don't contradict. Follow-up: remaining
   Chromosome-analysis units (karyotype, centromere/telomere, arm-ratio, synteny, GC-skew) warrant their
   own concepts when ingested.
   graph: +2 nodes, +1 typed edge (relates_to test-unit-registry)

## [2026-07-09] ingest | docs/Evidence/CHROM-CENT-001-Evidence.md → chrom-cent-001-evidence (source) + 1 concept
   Eighteenth per-algorithm Evidence file; second of the Chromosome-analysis family (after ANEU).
   Created the genuinely-distinct concept centromere-analysis — the anchor for the chromosome
   centromere / alpha-satellite family, sibling of aneuploidy-detection. Unusually this artifact is a
   layered multi-session record: base `AnalyzeCentromere` (generic tandem-repeat-density heuristic,
   sliding-window k-mer + low GC-variability + k=15 repeat content; AlphaSatelliteContent is a repeat
   score NOT alpha-satellite-specific) + Levan 1964 q/p arm-ratio classification (exact thresholds
   1.7/3.0/7.0/∞ → Metacentric/Submetacentric/Subtelocentric/Acrocentric/Telocentric) + four opt-in
   additive detectors: DetectAlphaSatellite/FindCenpBBoxes (171-bp tandem period ±5, ≥0.50
   self-similarity, AT>0.50, 17-bp CENP-B box IUPAC `YTTCGTTGGAARCGGGA` — no embedded monomer string),
   DetectHigherOrderRepeat (split into 171-bp monomers, GlobalAlign+CalculateStatistics, HOR period =
   smallest k with monomers k apart ≥95% identical / <5% divergence across ≥90% of array; intra-HOR
   50–70% vs inter-HOR 97–100%; period 1 = homogeneous 1-mer not multi-monomer HOR), and
   AssignSuprachromosomalFamily (bundled CC0 Dfam ALR/ALRa=A, ALRb=B via CENP-B box; ≥60% gate; SF3
   pentameric period%5==0 / SF4 monomeric A-type / {SF1,SF2} dimeric A→B / SF5 irregular). Sources:
   Wikipedia Centromere/Karyotype/Chromosome + Levan 1964 + Hartley/O'Neill 2019 & McNulty/Sullivan
   2018 (PMC6121732) + Masumoto 1989 (PMC4843215) + Rosandić 2024 (PMC11050224) & Alkan 2007/ColorHOR +
   Shepelev 2009 + Dfam (CC0) + T2T/CHM13 (CC0). Concise source page for the artifact. Linked new
   source + concept into the algorithm-validation-evidence hub and added CHROM-CENT to that hub's
   frontmatter; added a reciprocal sibling nav link from aneuploidy-detection. Two flagged ASSUMPTION
   parameters (≥60% alpha-satellite gate, SF3⇔period%5==0 pentameric proxy). Contradictions: none — the
   encyclopedic + alphoid-DNA literature + Dfam/T2T reference agree (171-bp monomer, 17-bp box, <5%
   inter-HOR recur). Residual data-blocked limitation: SF1-vs-SF2 not separated and diverged-pentamer
   SF3 (e.g. DXZ1 period 12) not tagged — needs an SF-resolved consensus monomer library that is
   non-redistributable (no LICENSE / non-machine-retrievable supplements). Follow-up: remaining
   Chromosome-analysis units (telomere, arm-ratio, synteny, GC-skew) warrant their own concepts.
   graph: +2 nodes, +2 typed edges (relates_to test-unit-registry, relates_to aneuploidy-detection)

## [2026-07-09] ingest | docs/Evidence/CHROM-TELO-001-Evidence.md → chrom-telo-001-evidence (source) + 1 concept
   Twenty-first per-algorithm Evidence file; fifth of the Chromosome-analysis family (after ANEU, CENT,
   KARYO, SYNT). Created the genuinely-distinct concept telomere-analysis — the anchor for the chromosome
   telomere family, sibling of aneuploidy-detection. Two parts: repeat detection at each end (3′ forward
   TTAGGG / 5′ reverse-complement CCCTAA, 6-bp vertebrate unit) with a configurable 70%-per-window purity
   threshold (5/6 match for 6-bp, 5/7 for Arabidopsis TTTAGGG; divergent TTAGGA → purity 5/6≈0.833),
   purity∈[0,1] tracked with length; and length estimation direct from the run or via the qPCR T/S ratio
   (Cawthon 2002 linearity EstimatedLength=referenceLength×(tsRatio/referenceRatio)). Invariants: length≥0,
   purity∈[0,1], threshold consistency (has⇒len≥minTelomereLength), IsCriticallyShort=(hasTelomere&&
   len<criticalLength)OR empty, orientation (5′=revcomp / 3′=forward). Four configurable parameters flagged
   as implementation defaults NOT biological constants: criticalLength 3000, minTelomereLength 500,
   searchLength (truncates reported length), referenceLength 7000. Sources: Wikipedia Telomere + Meyne 1989
   (TTAGGG conserved across vertebrates, PMID 2780561) + Cawthon 2002 (T/S ∝ length, r²=0.677, PMID
   12000852) + Blackburn-Gall 1978. Oracles: 200×TTAGGG on 1000 A's → len 1200/purity 1.0, both-ends,
   no-telomere/empty→critically-short, TTAGGA×200→0.833, TTAGGG×2000→12000, searchLen 600→truncate to 600;
   T/S table {1.0,1.5,0.5,2.0}@7000→{7000,10500,3500,14000}, refRatio 2.0→3500, 0.0→0. Concise source page
   for the artifact. Linked new source + concept into the algorithm-validation-evidence hub and added
   CHROM-TELO to that hub's frontmatter; added a reciprocal sibling nav link from aneuploidy-detection.
   Contradictions: none — Deviations and Assumptions is None; Wikipedia repeat table, Meyne 1989 repeat
   conservation, and Cawthon 2002 T/S proportionality agree. Follow-up: remaining Chromosome-analysis units
   (arm-ratio, GC-skew) warrant their own concepts; non-vertebrate telomere repeats (documented in the
   species table) need their own repeat unit + per-window match count.
   graph: +2 nodes, +2 typed edges (relates_to test-unit-registry, relates_to aneuploidy-detection)

## [2026-07-09] ingest | docs/Evidence/CHROM-KARYO-001-Evidence.md → chrom-karyo-001-evidence (source) + 1 concept
   Nineteenth per-algorithm Evidence file; third of the Chromosome-analysis family (after ANEU, CENT).
   Created the genuinely-distinct concept karyotype-analysis — the anchor for the chromosome
   karyotyping / ploidy-detection family, sibling of aneuploidy-detection and centromere-analysis. Two
   independent algorithms: `AnalyzeKaryotype` karyotypes chromosome *descriptors* (split sex/autosome,
   group autosomes by base name, count copies, compare to expected ploidy, label nullisomy..pentasomy
   by absolute copy count — same cytogenetic ladder as ANEU but keyed on descriptor counts not depth
   log-ratios; TotalChromosomes/TotalGenomeSize/MeanChromosomeLength invariants); and `DetectPloidy`
   estimates whole-genome ploidy from read depth (true median of sorted depths, ratio=medianDepth/
   expectedDiploidDepth, ploidy=round(ratio×2) clamped [1,8] — note the [1,8] whole-genome clamp and
   direct ratio vs ANEU's [0,10] CN clamp via 2^logRatio; confidence=1−|ratio×2−ploidy|×2, empty→(2,0)).
   Sources all Wikipedia (Karyotype/Ploidy/Aneuploidy, verified 2026-03-08). Oracles: normal diploid
   human 46/no-aneuploidy, Down (3×chr21→Trisomy), Turner (45,X→Monosomy), disomy-in-tetraploid,
   tetrasomy/pentasomy, diploid/tetraploid/haploid depth ratios. Five design decisions captured (DD1
   empty-karyotype, DD2 empty-depth→(2,0), DD3 [1,8] clamp, DD4 nullisomy unreachable via GroupBy —
   absent chromosomes form no group, term mapped for completeness only, DD5 disomy is aneuploidy only
   in non-diploid/ISCN contexts). Concise source page for the artifact. Linked new source + concept
   into the algorithm-validation-evidence hub and added CHROM-KARYO to that hub's frontmatter; added a
   reciprocal sibling nav link from aneuploidy-detection. Contradictions: none — the artifact's
   Deviations and Assumptions section is None; the three Wikipedia sources agree and DD4/DD5 are
   architecture/nomenclature notes not departures. Follow-up: remaining Chromosome-analysis units
   (telomere, arm-ratio, synteny, GC-skew) warrant their own concepts when ingested.
   graph: +2 nodes, +2 typed edges (relates_to test-unit-registry, relates_to aneuploidy-detection)

## [2026-07-09] ingest | docs/Evidence/CHROM-SYNT-001-Evidence.md → chrom-synt-001-evidence (source) + 1 concept
   Twentieth per-algorithm Evidence file; fourth of the Chromosome-analysis family (after ANEU, CENT,
   KARYO). Created the genuinely-distinct concept synteny-and-rearrangement-detection — deliberately
   named as the SHARED synteny anchor so the upcoming comparative-genomics COMPGEN-SYNTENY-001 unit
   reuses it rather than re-deriving syntenic-block definitions. Two algorithms: `FindSyntenyBlocks`
   (group ortholog pairs by chromosome pair, sort by reference position, identify collinear runs, merge
   consecutive segments under maxGap, emit blocks ≥ minGenes; each block carries strand '+'/'-',
   GeneCount, and SequenceIdentity=NaN — not computable from coordinate-only input per MCScanX; I1–I5)
   and `DetectRearrangements` (sort blocks by ref chr/position, compare adjacent pairs: different target
   chr → Translocation, same target chr + different strand → Inversion, gap asymmetry → Deletion,
   overlapping source coords + different targets → Duplication; Type recognized-value + Position1 non-null
   + Chromosome2-differs invariants). Sources: Wikipedia Synteny/Comparative-genomics/Chromosomal-
   rearrangement + MCScanX (Wang 2012), SyRI (Goel 2019), Liu 2018, MUMmer. Oracles: collinear-forward
   (4 genes chr1→chrA → 1 block '+', 1000–8000), inverted block ('-'), translocation (chrA→chrB @ 50000),
   inversion (positions 50000/60000, size 10000). Captured artifact §7 coverage-strengthen (8 weak→exact,
   2 duplicate removed, 1 missing M16 maxGap-split implemented). Concise source page for the artifact.
   Linked new source + concept into the algorithm-validation-evidence hub and added CHROM-SYNT to that
   hub's frontmatter; added a reciprocal sibling nav link from aneuploidy-detection. Contradictions: none
   — the artifact's Deviations section is None; Wikipedia synteny/rearrangement definitions and MCScanX/
   SyRI tool descriptions agree; SequenceIdentity=NaN is MCScanX-backed. Follow-up: remaining Chromosome-
   analysis units (telomere, arm-ratio, GC-skew) warrant their own concepts; COMPGEN-SYNTENY-001 shares
   synteny-and-rearrangement-detection when ingested.
   graph: +2 nodes, +2 typed edges (relates_to test-unit-registry, relates_to aneuploidy-detection)

## [2026-07-09] ingest | docs/Evidence/CODON-CAI-001-Evidence.md → codon-cai-001-evidence (source) + 1 concept
   Twenty-second per-algorithm Evidence file; second of the codon-usage family (after ANNOT-CODONUSAGE/
   RSCU). Created the genuinely-distinct concept codon-adaptation-index (CAI, Sharp & Li 1987): a
   whole-gene directional codon-bias / expression-proxy score in [0,1] = geometric mean of relative
   adaptiveness w_i = f_i/max_synonym_f (family-MAX normalization, one level above RSCU's family
   normalization), equivalently exp((1/L)Σ ln w); the geometric mean makes it low-value-sensitive (one
   rare codon drags CAI down). Captured: stop codons excluded; single-codon Met/Trp w≡1 → canonical
   Sharp & Li 1987 / Jansen 2003 EXCLUSION rule (quoted verbatim in the artifact), exposed as the opt-in
   excludeSingleCodonAminoAcids flag (default includes them, historical); exclude can yield L=0→CAI 0;
   E. coli K12 (Kazusa 316407) oracles AUG→1.0, CUG-CCG-ACC→1.0, CUA-CCA-ACA→0.1980, plus the four
   exclusion-mode cases. Concise source page for the artifact (Wikipedia + Sharp & Li 1987 + Jansen 2003
   PMC2684136 + Kazusa sources). Linked new source + concept into the algorithm-validation-evidence hub
   and added CODON-CAI to that hub's frontmatter; cross-linked bidirectionally with
   relative-synonymous-codon-usage (CAI reuses RSCU-style weights). One deviation: the Seqeron
   implementation clamps zero-frequency codons (freq=0 but family maxFreq>0) to w=1e-6 (incomplete-table
   protection) rather than strict w=0/log(0); unknown-AA/maxFreq=0 → NaN skipped; empty → 0.
   FLAGGED cross-page nuance (not a source contradiction): the RSCU page described CAI's log(0) guard as
   Sharp & Li's "0.5 pseudocount" (a reference-table-build convention) whereas this implementation uses
   a 1e-6 score-time clamp — reconciled the RSCU page wording and noted it on both pages. Sources agree
   internally (Wikipedia formulae = Sharp & Li = Jansen exclusion quote). Follow-up: remaining codon-usage
   units (CODON-ENC/OPT/RARE/STATS/USAGE, SEQ-CODON-FREQ, TRANS-CODON) — ENC/optimization likely warrant
   their own concepts, raw frequency/usage tables may share existing concepts.
   graph: +2 nodes, +2 typed edges (relates_to test-unit-registry, relates_to relative-synonymous-codon-usage)

## [2026-07-09] ingest | docs/Evidence/CODON-ENC-001-Evidence.md → codon-enc-001-evidence (source) + 1 concept
   Twenty-third per-algorithm Evidence file; third of the codon-usage family (after ANNOT-CODONUSAGE/
   RSCU, CODON-CAI). Created the genuinely-distinct concept effective-number-of-codons (ENC / Nc,
   Wright 1990): a reference-free whole-gene codon-bias measure in [20,61] = reciprocal of codon
   homozygosity F̂=(n·Σp_i²−1)/(n−1) aggregated by degeneracy class as Nc=2+9/F̂₂+1/F̂₃+5/F̂₄+3/F̂₆
   (constants = the standard-code partition: 2 singlets Met/Trp, 9 doublets, 1 triplet Ile, 5 quartets,
   3 sextets; stops excluded). 20 = max bias (one codon per aa), 61 = uniform usage; sampling
   WITHOUT replacement (Fuglsang 2006 superior estimator). Corner cases: n≤1 → F̂ undefined, drop and
   within-class-average (Eq. 4); Ile-absent 3-fold → F̂₃=(F̂₂+F̂₄)/2 fallback (Eq. 5a); Eq. 3 overshoot
   > 61 → re-adjust to 61; per-aa overshoot on small n. Oracles: unbiased→61, max-bias→20, Fuglsang
   no-bias-discrepancy sim→40.5, Phe TTT×3/TTC×1→F̂=0.5/Nc=2 (even split→3). One assumption: the
   Math.Max(20,…) lower clamp is a defensive bound NOT a Wright instruction (only re-adjust-DOWN-to-61
   is source-prescribed; 20 is the structural minimum). Concise source page for the artifact (Fuglsang
   2004 rank-1 verbatim Wright equations + Fuglsang 2006 Genetics rank-1 + NCBI degeneracy partition
   rank-2 + Wright 1990 Gene primary). Linked new source + concept into the algorithm-validation-evidence
   hub and added CODON-ENC to that hub's frontmatter; cross-linked bidirectionally with
   codon-adaptation-index (ENC = reference-free counterpart, modeled alternative_to) and
   relative-synonymous-codon-usage (F̂ built from the same p_i). Contradictions: none — Fuglsang 2004 &
   2006 reproduce identical Wright equations, NCBI partition matches Eq. 3 constants. Follow-up:
   remaining codon-usage units (CODON-OPT/RARE/STATS/USAGE, SEQ-CODON-FREQ, TRANS-CODON) — optimization
   likely warrants its own concept, raw frequency/usage tables may share existing concepts.
   graph: +2 nodes, +3 typed edges (relates_to test-unit-registry, alternative_to codon-adaptation-index, relates_to relative-synonymous-codon-usage)

## [2026-07-09] ingest | docs/Evidence/CODON-OPT-001-Evidence.md → codon-opt-001-evidence (source) + 1 concept
   Twenty-fourth per-algorithm Evidence file; fourth of the codon-usage family (after ANNOT-CODONUSAGE/
   RSCU, CODON-CAI, CODON-ENC). Created the genuinely-distinct concept codon-optimization — the family's
   sole *rewriting* operation (`OptimizeSequence`): synonymous-codon substitution to improve heterologous
   host expression, deliberately positioned as the ACTUATOR to RSCU/CAI/ENC's measurement. Five strategies
   each traced to a source point: MaximizeCAI (most-frequent codon, Sharp & Li 1987 CAI), BalancedOptimization
   (CAI vs 40-60% GC, rebuilds Changes list after GC balancing), HarmonizeExpression (match host distribution,
   Mignon 2018), AvoidRareCodons (replace only sub-threshold codons), MinimizeSecondary (delegates to
   BalancedOptimization for selection + dedicated ReduceSecondaryStructure). Invariants: protein preservation
   across all strategies (synonymous only), Met/AUG & Trp/UGG fixed points (single-codon families), stop
   preserved, CAI∈(0,1]; RNA notation with T→U, trim-to-complete-codon, case-insensitive, 1e-6 zero-freq CAI
   clamp (same guard as codon-adaptation-index). Organism fixtures: E. coli K12 (Kazusa 316407), S. cerevisiae
   (4932), H. sapiens (9606) preferred-codon tables. Concise source page for the artifact (Wikipedia
   codon-usage-bias + CAI + Plotkin-Kudla 2011 + Mignon 2018 + Kazusa sources). Linked new source + concept
   into the algorithm-validation-evidence hub and added CODON-OPT to that hub's frontmatter; cross-linked the
   new concept from codon-adaptation-index (MaximizeCAI drives CAI→1), relative-synonymous-codon-usage, and
   effective-number-of-codons (all three "sibling still in docs/Evidence" lines now resolve to the ingested
   page). Also removed stray `</content></invoke>` tags left at the tail of effective-number-of-codons.md
   (pre-existing Write artifact). Contradictions: none — Wikipedia strategy catalogue, Sharp & Li 1987,
   Plotkin-Kudla 2011 review and Mignon 2018 harmonization agree; behaviours recorded "from theory/
   implementation" so the correctness anchor is the protein-preservation invariant + CAI formula, both
   source-backed. Follow-up: remaining codon-usage units (CODON-RARE/STATS/USAGE, SEQ-CODON-FREQ, TRANS-CODON)
   — rare-codon analysis may warrant its own concept, raw frequency/usage tables may share existing concepts.
   graph: +2 nodes, +2 typed edges (relates_to test-unit-registry, relates_to codon-adaptation-index)

## [2026-07-09] ingest | docs/Evidence/CODON-RARE-001-Evidence.md → codon-rare-001-evidence (source) + 1 concept
   Twenty-fifth per-algorithm Evidence file; fifth of the codon-usage family (after ANNOT-CODONUSAGE/RSCU,
   CODON-CAI, CODON-ENC, CODON-OPT). Judged genuinely distinct — beyond a thresholded-frequency view, the
   2026-06-24 addendum adds two published cluster-detection algorithms with their own formulas and oracles —
   so created the concept rare-codon-analysis (the family's localization/diagnostic unit vs RSCU/CAI/ENC's
   whole-gene summaries; feeds codon-optimization's AvoidRareCodons). Base `FindRareCodons`: codon rare when
   per-family frequency < threshold (default 0.15), reports 0-indexed position(×3)/AA/actual-freq; E. coli
   K12 rare set AGA 0.04/AGG 0.02/CGA 0.06/CUA 0.04 (Kazusa MG1655); invariants (pos multiples of 3 in
   [0,len−3], freq∈[0,1], reported<threshold, deterministic), edges (empty→empty, non-÷3 trailing ignored,
   T→U, threshold 0/1 extremes, unknown codon→freq 0 always reported). Addendum cluster methods: Clarke &
   Clark 2008 %MinMax (per-AA synonymous Xij/Xmax/Xmin/Xavg, signed %Min/%Max ∈[−100,+100], default 18-codon
   window, rare clusters = negative %Min peaks) + Chartier/Gaudreault/Najmanovich 2012 Sherlocc (7-codon
   window, ≥4 of 7 "slow"=freq≤threshold positions → rare-codon cluster). Corner cases: single-codon AA
   (Met/Trp) contributes 0 to %MinMax num+denom → no divide-by-zero/NaN; window>seq → none; overlapping
   qualifying windows merged into one maximal cluster (flagged implementation choice, Sherlocc reports
   regions). Arg-family oracles AGA³→−86.36% / CGC³→+100% / CUG·AGA→+36.47% and Sherlocc 7×AGA→1 cluster /
   4+3→cluster / 3+4→none. Sources: Wikipedia codon-usage-bias + GenScript GenRCA (Fan 2024 BMC
   Bioinformatics) + Kazusa + Shu 2006 (5×CUA ~3-fold inhibition, PMC6032470) + Sharp & Li 1987 + Clarke &
   Clark 2008 (PLoS ONE) + Rodriguez 2018 (%MinMax) + Chartier 2012 (Bioinformatics, DOI bts149) +
   mtthchrtr/sherlocc README. Concise source page for the artifact. Linked new source + concept into the
   algorithm-validation-evidence hub and added CODON-RARE to that hub's frontmatter; added reciprocal nav
   links from codon-optimization (AvoidRareCodons actuator) and codon-adaptation-index (localizes the low-w
   codons that pull CAI down). Contradictions: none — Deviations and Assumptions is None; the base
   threshold-frequency approach and the two complementary cluster methods each cite peer-reviewed sources
   plus a reference implementation; the overlapping-window merge is an explicitly flagged choice. Follow-up:
   remaining codon-usage units (CODON-STATS/USAGE, SEQ-CODON-FREQ, TRANS-CODON) — raw frequency/usage tables
   may share existing concepts.
   graph: +2 nodes, +2 typed edges (relates_to test-unit-registry, relates_to codon-optimization)

## [2026-07-09] ingest | docs/Evidence/CODON-RSCU-001-Evidence.md → codon-rscu-001-evidence (source)
   Twenty-sixth per-algorithm Evidence file; sixth of the codon-usage family and the SECOND RSCU
   unit (after ANNOT-CODONUSAGE-001). REUSED the existing relative-synonymous-codon-usage concept
   rather than duplicating — this file validates the same n_i·x/Σx measure but adds the supporting
   `CountCodons` counting operation and a broader reference panel. No new concept created. Concise
   source page for the CODON-RSCU-001 artifact: Sharp/Tuohy/Mosurski 1986 (NAR 14(13):5125-5143, the
   RSCU-introducing paper) + LIRMM RSCU RS + GenomicSig (CRAN) + seqinr `uco` + cubar `est_rscu` +
   PMC2528880 sources; Leu CTGCTGCTGCTA→4.5/1.5/0 (Σ=6), Phe TTTTTTTTC→4/3,2/3 (Σ=2), unbiased
   TTTTTC→1.0, Met ATGATG→1.0, CountCodons frame/exclusion (ATGAAATGA/ATGAA/ATGNNNAAA) datasets; two
   assumptions (absent-family 0/0→0 vs cubar pseudocount default 1; stop codons as a degeneracy-3
   family). Surgically updated the RSCU concept (added CODON-RSCU source + a 2nd relates_to
   test-unit-registry edge; documented the CountCodons counting contract and the primary-attribution
   note) and the algorithm-validation-evidence hub (frontmatter + source-list). FLAGGED two cross-page
   nuances, neither a source contradiction: (1) stop-codon handling — the concept/ANNOT-CODONUSAGE say
   stops are EXCLUDED (Biopython forward_table) whereas CODON-RSCU says the repo treats the 3 stops as
   a degeneracy-3 synonymous family; both agree it never changes an amino-acid codon's RSCU; (2)
   primary attribution — CODON-RSCU + seqinr cite Sharp, Tuohy & Mosurski 1986 (the RSCU-introducing
   paper) whereas the concept/begomovirus restatement wrote "Sharp & Li 1986"; noted both on the
   concept. Sources internally consistent (LIRMM/GenomicSig/seqinr/Sharp-Tuohy-Mosurski algebraically
   identical; cubar pseudocount an explicit zero-division convention). Follow-up: remaining codon-usage
   units (CODON-STATS/USAGE, SEQ-CODON-FREQ, TRANS-CODON) — raw frequency/usage tables may share
   existing concepts.
   graph: +1 node, +1 typed edge (relates_to test-unit-registry)

## [2026-07-09] ingest | docs/Evidence/CODON-STATS-001-Evidence.md → codon-stats-001-evidence (source)
   Twenty-seventh per-algorithm Evidence file; seventh of the codon-usage family. NO new concept
   created — CODON-STATS-001 is the family's aggregation/reporting view (`GetStatistics` bundles codon
   counts + RSCU + ENC + CAI + positional GC + total codons; `CalculateCai` = the same CAI validated as
   CODON-CAI-001), so it REUSES the existing family concepts (relative-synonymous-codon-usage,
   codon-adaptation-index, effective-number-of-codons; cross-refs codon-optimization, rare-codon-analysis)
   rather than duplicating them — consistent with the "aggregation view reuses concepts" rule. The one
   piece not covered by an existing concept — positional GC composition (GC1/GC2/GC3, OverallGc, and
   GC3s = GC of *synonymous* third positions, excluding Met/Trp/stop, Peden 1999 §1.8.2.1.3 / EMBOSS cusp
   / PMC7596632 "59 synonymous codons") — is documented inline on the source page, not promoted to a
   concept (small facet of a stats bundle; a dedicated positional-GC unit can mint one later if needed).
   Sources: Sharp & Li 1987 (+ Biopython SharpEcoliIndex `w` reproduction) + Wikipedia + seqinr + CodonW/
   Peden thesis + EMBOSS cusp + Kazusa H. sapiens. Oracles: ATGGCA→GC3s 0 vs GC3 50 (shows the Met/Trp/
   stop exclusion), GCTGCC→CAI √0.122=0.34928…, ATGTGGTAA→CAI 0 (no scorable codon), CTGGTTAAA→GC1/GC2/
   GC3 66.67/0/33.33; EColiOptimalCodons reproduces Sharp&Li w, HumanOptimalCodons reproduces Kazusa RSCU
   (CTG≈2.3713). Two documented deviations: GC3s reported as a percentage (×100, EMBOSS-style, vs CodonW
   fraction — labeling only, subset exactly per Peden) and zero-`w` codons skipped rather than floored to
   0.01 (Bulmer 1988) — real-CDS CAI unaffected (no reference synonymous codon has w=0). Surgically updated
   the algorithm-validation-evidence hub (frontmatter sources + source-link list) and the RSCU base-anchor
   concept (noted CODON-STATS as the aggregation view + the GC3s definition). Contradictions: none — all
   sources agree on the formulae and the synonymous-codon exclusion set; the CAI zero-handling wording here
   (skip-zero-w) and CODON-CAI-001's (1e-6 clamp) describe the same guard from different angles, flagged on
   the source page. Follow-up: remaining codon-usage units (CODON-USAGE, SEQ-CODON-FREQ, TRANS-CODON) — raw
   frequency/usage tables likely share existing concepts.
   graph: no typed edges (aggregation source reusing existing concepts; source pages can't be relates_to
   subjects per the ontology, so no new node warrants an edge; mentions auto-derived from wikilinks)

## [2026-07-09] ingest | docs/Evidence/CODON-USAGE-001-Evidence.md → codon-usage-001-evidence (source) + 1 concept
   Twenty-eighth per-algorithm Evidence file; eighth of the codon-usage family. Created ONE new concept
   codon-usage-comparison — the *raw* end of the family (`CalculateCodonUsage` + `CompareCodonUsage`).
   Judged genuinely distinct despite the "raw table likely reuses RSCU" hint: the raw counting IS the
   RSCU primitive (documented as reuse + a relates_to edge to relative-synonymous-codon-usage), but
   `CompareCodonUsage`'s Total Variation Distance similarity — Similarity = 1 − Σ|f₁−f₂|/2 ∈ [0,1] between
   two codon-frequency distributions, with proven identity 1.0 / symmetry / range / disjoint→0 /
   partial-overlap-exact properties — is a distribution-comparison operation no existing bias concept
   (RSCU/CAI/ENC/optimization/rare/stats) provides, so it warrants its own page. `CalculateCodonUsage`:
   non-overlapping triplets from offset 0, T→U internally, uppercase, drop partial trailing codon, returns
   unnormalized Dictionary<codon,int> (Σcounts=total invariant). Oracles: AUGGCUGCU→{AUG:1,GCU:2}, all-64
   codons→64 keys count 1; TVD sims identical→1.0, disjoint UUU/GGG→0, 2/3-shared→2/3, 0.5/0.75/0.75/0.25
   cases (all analytically derivable). Sources: Wikipedia codon-usage-bias (degeneracy) + Kazusa CUTG
   format + Sharp & Li 1987 (per-AA normalization) + Plotkin-Kudla 2011 + Athey 2017; predefined E. coli
   K12 / S. cerevisiae / H. sapiens tables Kazusa-verified (all 64 relative fractions, March 2026). Two
   deviations, both deliberate/benign: TVD-not-cosine metric choice (Wikipedia lists cosine + correlation;
   every test value derivable from the TVD formula and the 4 proven properties follow from TVD theory) and
   empty→similarity 0 (no data → 0, not NaN/exception). Concise source page for the artifact. Linked new
   source + concept into the algorithm-validation-evidence hub (frontmatter sources + both link lists) and
   cross-linked from the RSCU base-anchor concept (raw-table sibling). Contradictions: none — Wikipedia,
   Kazusa, and Sharp & Li agree on the codon-usage biology; the TVD similarity is an implementation metric
   choice the sources don't contradict. Follow-up: remaining codon-usage units (SEQ-CODON-FREQ, TRANS-CODON)
   — raw frequency/usage tables may share codon-usage-comparison or RSCU rather than minting new concepts.
   graph: +2 nodes, +2 typed edges (relates_to test-unit-registry, relates_to relative-synonymous-codon-usage)
