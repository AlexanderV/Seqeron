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

## [2026-07-09] ingest | docs/Evidence/COMPGEN-ANI-001-Evidence.md → compgen-ani-001-evidence (source) + 1 concept
   Twenty-ninth per-algorithm Evidence file; FIRST of the Comparative-genomics (COMPGEN) family. Created
   the genuinely-distinct concept average-nucleotide-identity — the anchor for the COMPGEN ANI
   genome-similarity family, sibling of the already-existing shared synteny anchor
   synteny-and-rearrangement-detection (ANI = how nucleotide-identical two genomes are; synteny = whether
   their gene order is conserved). ANIb definition traced verbatim to Goris et al. 2007 (IJSEM): fragment
   the query into consecutive 1020 nt pieces (mirrors ~1 kb DDH shearing), best-match place each against
   the reference, keep only fragments with >30% identity over ≥70% alignable length — BOTH recalculated
   over the full query-fragment length (pyani ani_pid=ani_alnids/qlen, ani_coverage=ani_alnlen/qlen, not
   the local sub-region) — and average the qualifying per-fragment identities; non-conserved fragments are
   discarded, not zero-scored. Species boundary ANI≈95%↔70% DDH (Goris; Konstantinidis & Tiedje 2005 ≈94%).
   2026-06-23 refresh resolved the old ungapped assumption: gapped Smith-Waterman placement (pyani BLASTN
   -xdrop_gap_final 150, ani_alnlen=blast_alnlen-blast_gaps) recovers indels (AAAACCCC/AAAATCCCC 0.875→1.0),
   and CalculateReciprocalAni implements the reciprocal/symmetric value = mean of both directions
   (order-independent, since only the query is fragmented → single-direction is asymmetric). Exact-arithmetic
   oracles (identical→1.0, one mismatch→0.9375, AATT→0.875, CGTC-excluded→1.0, ref<frag→0, query<frag→0);
   null/empty→0, non-positive fragmentLength→ArgumentOutOfRangeException. Concise source page for the
   artifact (Goris 2007 + Konstantinidis & Tiedje 2005 rank-1 + pyani rank-3 sources). Linked new source +
   concept into the algorithm-validation-evidence hub (frontmatter sources + both link lists) and added a
   reciprocal sibling nav link from synteny-and-rearrangement-detection. One documented DECISION (not a
   correctness gap, not a deviation): the gapped path uses SequenceAligner.LocalAlign (full-DP Smith-Waterman,
   BLAST DNA scoring) rather than the NCBI BLASTN engine — more sensitive than BLAST's heuristic seeding, same
   recalculated-over-fragment identity/coverage; numeric ANI may differ slightly from NCBI-BLASTN pipelines,
   indel handling correct (algorithm doc §5.3). Contradictions: none — Goris, Konstantinidis & Tiedje, and
   pyani agree on fragmentation, cut-offs, averaging, gapped placement, and reciprocal computation. Follow-up:
   remaining COMPGEN units (orthologs/RBH, COMPGEN-SYNTENY-001 which reuses synteny-and-rearrangement-detection,
   reversal distance, dot-plot, conserved gene clusters) warrant their own concepts or share existing ones when
   ingested.
   graph: +2 nodes, +2 typed edges (concept relates_to test-unit-registry, relates_to synteny-and-rearrangement-detection)

## [2026-07-09] ingest | docs/Evidence/COMPGEN-CLUSTER-001-Evidence.md → compgen-cluster-001-evidence (source) + 1 concept
   Thirtieth per-algorithm Evidence file; second of the Comparative-genomics (COMPGEN) family (after ANI).
   Created the genuinely-distinct concept conserved-gene-clusters-common-intervals — a conserved gene
   cluster = a gene-label SET that is contiguous (an interval) in EVERY genome, order- and strand-free
   inside the window: the common-interval model. Distinct from both COMPGEN siblings — ANI measures
   nucleotide identity, synteny requires a COLLINEAR ORDERED block, a common interval only requires the
   same gene SET contiguous in each genome. Definitions traced verbatim to Uno & Yagiura 2000
   (Algorithmica, originating common-interval model, O(n²) LHP + output-sensitive O(n+K) RC), Heber &
   Stoye 2001 (CPM, k-permutation generalisation, optimal O(kn+z)/O(n)), Bui-Xuan/Habib/Paul 2013
   (MinMax-Profiles arXiv:1304.5140, unifying view: interval [i,j] defined only for i<j → size ≥2,
   singletons excluded, whole set (1..n) always common; golden-vector Example 1), and Didier et al. 2013
   (arXiv:1310.4290, extension permutations→sequences with duplicates: paralogs handled, a set is common
   iff SOME contiguous window/location in each genome has exactly that label set). Seqeron does the simple
   strict O(n²·K_genomes) check with a minClusterSize filter; K≥2 required (single genome → every interval
   trivially common). Oracles: golden vector Id₇ vs (7 2 1 3 6 4 5) → non-trivial {1,2}/{1,2,3}/{3,4,5,6}/
   {4,5}/{4,5,6}/{1..6} + trivial {1..7} (brute-force reproduced); split-negative {2,3} (positions 2,4 in
   genome 2 → non-adjacent); sequence-with-duplicates T/S → {1,2} not common, {1,2,3,4} common. Concise
   source page for the artifact (Uno & Yagiura 2000 + Heber & Stoye 2001 + Bui-Xuan/Habib/Paul 2013 +
   Didier et al. 2013, all rank-1). Linked new source + concept into the algorithm-validation-evidence hub
   (frontmatter sources + both link lists) and added reciprocal sibling nav links from
   average-nucleotide-identity and synteny-and-rearrangement-detection. ONE documented ASSUMPTION (API-shape,
   not a correctness gap): the public method keeps a maxGap parameter but the validated/tested behaviour is
   the STRICT gap-free common-interval model — maxGap does NOT relax it, and the gene-teams gapped extension
   (Bergeron, Corteel & Raffinot 2002) is NOT implemented (source not retrievable). Contradictions: none —
   the four sources agree on the interval definition, the contiguous-in-every-genome cluster rule, the
   size-≥2 constraint, and the sequence-with-duplicates generalisation. Follow-up: remaining COMPGEN units
   (orthologs/RBH, COMPGEN-SYNTENY-001 reusing synteny-and-rearrangement-detection, reversal distance,
   dot-plot, genome-comparison pipeline) warrant their own concepts or share existing ones when ingested.
   graph: +2 nodes, +2 typed edges (concept relates_to test-unit-registry, relates_to synteny-and-rearrangement-detection)

## [2026-07-09] ingest | docs/Evidence/COMPGEN-COMPARE-001-Evidence.md → compgen-compare-001-evidence (source) + 1 concept
   Thirty-first per-algorithm Evidence file; third of the Comparative-genomics (COMPGEN) family (after ANI,
   CLUSTER). This one is an ORCHESTRATING PIPELINE, not a single-metric unit: `CompareGenomes` performs the
   end-to-end two-genome comparison, partitioning each genome's genes into a CORE (conserved) set and a
   DISPENSABLE (genome-specific) set and reporting an OverallSynteny fraction. Created the genuinely-distinct
   concept genome-comparison-core-dispensable — distinct because it composes sub-units into the pairwise
   pan-genome model (Tettelin et al. 2005 PNAS, the paper that coined pan-genome/core/dispensable): core =
   the reciprocal-best-hit ortholog pairs (Moreno-Hagelsieb & Latimer 2008 + Tatusov 1997, COMPGEN-RBH-001),
   dispensable = the rest of each genome; outputs Conserved/Specific1/Specific2 + OverallSynteny = (genes in
   MCScanX syntenic blocks)/min(|g1|,|g2|) clamped ≤1 (fraction-of-syntenic-genes metric; blocks from MCScanX
   Wang 2012 = COMPGEN-SYNTENY-001, min 5 collinear anchors). Oracles: one-shared-one-unique → 1/1/1,
   disjoint → 0/2/2, identical-5-collinear+1-unique → Conserved 5, Specific 1/1, OverallSynteny 5/6=0.8333,
   0 rearrangements; symmetric partition (swap g1/g2 swaps Specific1/Specific2); empty genomes → all 0.
   Concise source page for the artifact (Tettelin 2005 + Moreno-Hagelsieb 2008/Tatusov 1997 + ScienceDirect/
   Wikipedia synteny overview + MCScanX Wang 2012 sources). Linked new source + concept into the
   algorithm-validation-evidence hub (frontmatter sources + both link lists) and added reciprocal sibling nav
   links from average-nucleotide-identity, conserved-gene-clusters-common-intervals, and
   synteny-and-rearrangement-detection. TWO documented ASSUMPTIONs, both source-backed, neither a
   partition-logic gap: (1) alignment-free 5-mer Jaccard RBH gate (identity ≥0.3, coverage ≥0.5) replaces the
   Tettelin 50%/50% alignment gate, inherited verbatim from COMPGEN-RBH-001 — partition logic unchanged;
   (2) MCScanX ≥5-collinear-anchor block threshold means OverallSynteny can be 0 even when conserved orthologs
   exist. Contradictions: none — Tettelin (core/dispensable), Moreno-Hagelsieb/Tatusov (RBH), and the synteny
   sources each govern a distinct pipeline output and are mutually consistent; Deviations = None. Follow-up:
   remaining COMPGEN units (RBH orthologs, reversal distance, dot-plot) warrant their own concepts when ingested.
   graph: +2 nodes, +3 typed edges (concept relates_to test-unit-registry, relates_to synteny-and-rearrangement-detection, relates_to average-nucleotide-identity)

## [2026-07-09] ingest | docs/Evidence/COMPGEN-DOTPLOT-001-Evidence.md → compgen-dotplot-001-evidence (source) + 1 concept
   Thirty-second per-algorithm Evidence file; fourth of the Comparative-genomics (COMPGEN) family (after ANI,
   CLUSTER, COMPARE). Created the genuinely-distinct concept dot-plot-word-match — genuinely distinct from the
   metric (ANI), ordered-block (synteny), gene-set (conserved-clusters), and pipeline (genome-comparison)
   siblings because it is the VISUAL word-match / k-tuple dot matrix that keeps the whole match relation as a
   2-D plot. Algorithm traced verbatim: dot at (x,y) iff the length-`wordSize` word starting at x in sequence1
   exactly matches the word at y in sequence2 (EMBOSS `dottup` exact word match — NOT scored `dotmatcher`),
   case-insensitive (both upper-cased), x=seq1 / y=seq2, all overlapping occurrences via suffix tree; `wordSize`
   default 10 (EMBOSS) is the noise-vs-sensitivity trade-off (longer=less noise/faster/less sensitive), `stepSize`
   subsamples x. Diagonals = similarity, main diagonal = self-comparison, repeats = extra diagonals, indels break
   the diagonal (Wikipedia). Oracles: Huttley `AGCGT`/`AT` k=1 → exactly {(0,0),(4,1)}; `ACGTACGT` self wordSize4
   → {(0,0),(0,4),(1,1),(2,2),(3,3),(4,0),(4,4)} (all overlapping word starts); `ACGT` self main diagonal.
   Corner cases: word>sequence / null / empty / disjoint-alphabet → no dots; non-positive wordSize/stepSize →
   ArgumentOutOfRangeException. Sources: Gibbs & McIntyre 1970 (Eur.J.Biochem. 16:1–11, rank 1, paywalled → method
   via secondaries, only citation/DOI attributed) + EMBOSS `dottup` manual+manpage (rank 3, default wordsize 10) +
   Wikipedia Dot plot (rank 4) + Huttley TIB Dotplot (rank 4, k=1 worked example). Concise source page for the
   artifact. Linked new source + concept into the algorithm-validation-evidence hub (frontmatter sources + both
   link lists) and added a reciprocal sibling nav link from average-nucleotide-identity. TWO ASSUMPTIONs, both
   explicitly non-correctness-affecting: (1) coordinate orientation x=seq1/y=seq2 (a presentation convention;
   transposing mirrors the plot but not the match set as a relation); (2) case-insensitive comparison (dottup/Gibbs
   do not mandate case folding; impl upper-cases both). Contradictions: none — Gibbs & McIntyre (via secondaries),
   Wikipedia, EMBOSS dottup, and Huttley agree on the exact-word match rule, diagonals-as-similarity, and the
   wordSize noise/sensitivity trade-off; Deviations = None. Follow-up: remaining COMPGEN units (RBH orthologs,
   COMPGEN-SYNTENY-001 reusing synteny-and-rearrangement-detection, reversal distance) warrant their own concepts
   or share existing ones when ingested.
   graph: +2 nodes, +3 typed edges (concept relates_to test-unit-registry, relates_to average-nucleotide-identity, relates_to synteny-and-rearrangement-detection)

## [2026-07-09] ingest | docs/Evidence/COMPGEN-ORTHO-001-Evidence.md → compgen-ortho-001-evidence (source) + 1 concept
   Thirty-third per-algorithm Evidence file; fifth of the Comparative-genomics (COMPGEN) family (after ANI,
   CLUSTER, COMPARE, DOTPLOT). Created the genuinely-distinct concept ortholog-detection-reciprocal-best-hits
   — the homology-classification unit and the shared RBH/ortholog anchor deliberately scoped so the future
   COMPGEN-RBH-001 unit reuses it, and the already-ingested genome-comparison-core-dispensable pipeline's
   conserved/core set IS these RBH pairs. Two rules traced verbatim: (1) ORTHOLOGS by Reciprocal Best Hits
   (Moreno-Hagelsieb & Latimer 2008: two genes in two genomes are orthologs iff each is the other's best hit;
   Tatusov 1997 COG symmetrical best hits; Fitch 1970 orthology=speciation / paralogy=duplication) — best hit
   = max-similarity candidate with deterministic tie-break (descending bit-score then ascending E-value),
   RECIPROCITY MANDATORY so a one-directional best hit (A→B, B→C≠A) is NOT an ortholog (the guarded defect
   class), ≥50% coverage gate + max E-value 1e-6 significance gate; (2) recent (IN-)PARALOGS by within-genome
   mutual best hits (Remm/Storm/Sonnhammer 2001 InParanoid in-paralog rule; out-paralogs pre-speciation
   excluded). Partial-matching output, determinism, empty-sequence/null contracts. Oracles: reciprocity
   {a1↔b1,a2↔b2}, non-reciprocity (b2=a1's superstring shares all 5-mers but a1↛b2 → RBH count 1), in-paralog
   {p1↔p2} with unrelated q1 excluded, empty→no orthologs / single-gene→no paralogs. Sources: Fitch 1970
   (Syst.Zool. 19:99-113, via Koonin 2011 PMC3178060 verbatim quote) + Tatusov 1997 (Science 278:631-637,
   full text 403/404-blocked, method via search summary+scirp, DOI confirmed) + Moreno-Hagelsieb 2008
   (Bioinformatics 24:319-324) + Remm 2001 (JMB 314:1041-1052, PMC5674930 corroboration) + Li 2003 OrthoMCL +
   Ondov 2016 Mash (alignment-free basis). Concise source page for the artifact. Linked new source + concept
   into the algorithm-validation-evidence hub (frontmatter sources + both link lists); cross-linked the
   already-ingested genome-comparison-core-dispensable concept (replaced its bare "COMPGEN-RBH-001" reference
   with [[ortholog-detection-reciprocal-best-hits]] in the intro + core/conserved bullet) and added a
   navigation link in the compgen-compare-001-evidence source page. ONE ASSUMPTION, source-backed and
   non-correctness-affecting: alignment-free 5-mer Jaccard replaces the BLAST bit-score ranking (the
   ComparativeGenomics class does not reference the Alignment project; cf. Mash) — affects only which
   near-identical pair wins ties; the RBH reciprocity rule, coverage gate (→ shared k-mers ≥50% of smaller
   set), and min-similarity gate are source-backed. Contradictions: none — Fitch/Tatusov/Moreno-Hagelsieb/Remm
   are mutually consistent, each governing a distinct part of the rule; Deviations = None beyond the metric
   substitution. Follow-up: COMPGEN-RBH-001 (apparent RBH duplicate) shares this concept rather than
   duplicating; remaining COMPGEN units (COMPGEN-SYNTENY-001 reusing synteny-and-rearrangement-detection,
   reversal distance) warrant their own concepts or share existing ones.
   graph: +2 nodes, +3 typed edges (concept relates_to test-unit-registry, genome-comparison-core-dispensable, average-nucleotide-identity)

## [2026-07-09] ingest | docs/Evidence/COMPGEN-RBH-001-Evidence.md → compgen-rbh-001-evidence (source)
   Thirty-fourth per-algorithm Evidence file; sixth of the Comparative-genomics (COMPGEN) family (after
   ANI, CLUSTER, COMPARE, DOTPLOT, ORTHO). NO new concept — this file is the RBH-only slice of the
   already-ingested COMPGEN-ORTHO-001, and the ortholog-detection-reciprocal-best-hits concept was
   deliberately scoped (during the ORTHO ingest) as the shared RBH anchor COMPGEN-RBH-001 reuses.
   Reused that concept: added COMPGEN-RBH-001 to its sources frontmatter and rewrote its intro to cite
   BOTH validation records (COMPGEN-ORTHO-001 = RBH + within-genome in-paralog; COMPGEN-RBH-001 = the
   between-genome ortholog slice, no in-paralog rule). Distinctive content vs ORTHO: only TWO sources
   (Moreno-Hagelsieb & Latimer 2008 operational RBH + thresholds; Tatusov 1997 COG genome-specific BeTs/
   mutually-consistent-best-hit triangles) — no Fitch, no Remm, no in-paralog dataset; verbatim quotes
   from the fetched OUP article (best hit = descending bit-score then ascending E-value; ≥50% coverage of
   "any of the protein sequences"; max E-value 1e-6); and an in-file DATA-QUALITY note where a
   search-engine summary claiming a 60% coverage threshold was rejected in favor of the article body's
   50% (recorded as a resolved discrepancy, not a source-vs-source contradiction). Datasets: reciprocity
   {a1↔b1,a2↔b2}, non-reciprocity (b2 = a1's superstring shares all 5-mers but a1↛b2 → RBH count 1),
   coverage/min-identity gate. Concise source page written; linked into the algorithm-validation-evidence
   hub (frontmatter sources + link list). ONE ASSUMPTION, source-backed: alignment-free 5-mer Jaccard
   replaces the BLAST bit-score ranking (cf. Mash) — affects only near-identical tie-breaks; reciprocity/
   coverage/threshold semantics unchanged. Contradictions: none between sources (Tatusov symmetrical BeTs
   and Moreno-Hagelsieb operational RBH are consistent, the latter the pairwise operationalization of the
   former); Deviations = None beyond the metric substitution. No new typed graph edges (reused existing
   concept, no new concept/relationship; mentions auto-derived).

## [2026-07-09] ingest | docs/Evidence/COMPGEN-REARR-001-Evidence.md → compgen-rearr-001-evidence (source) + 1 concept
   Thirty-fifth per-algorithm Evidence file; seventh of the Comparative-genomics (COMPGEN) family
   (after ANI, CLUSTER, COMPARE, DOTPLOT, ORTHO, RBH). Created a NEW concept
   genome-rearrangement-breakpoint-distance — the signed-permutation / breakpoint formulation of
   rearrangement detection (Hannenhalli–Pevzner / Bafna–Pevzner), genuinely distinct from the existing
   block-signal synteny-and-rearrangement-detection (CHROM-SYNT-001), which classifies from adjacent
   synteny-block coordinates. Modeled the two as alternative_to: this unit counts breakpoints b(α) on a
   signed permutation (extended (0,…,n+1), breakpoint = consecutive pair where neither (x,y) nor (−y,−x)
   survives in β), reports the breakpoint distance d_BP=n−sim(common adjacencies) and the reversal-distance
   lower bound d≥b/2, and ClassifyRearrangement returns Inversion (sign-flip reversal) vs Transposition
   (orientation-preserving block move) — Translocation/Deletion/Insertion/Duplication are a documented
   "Not implemented" (a single in-order permutation can't express them). Sources all rank 1: Hunter
   College Lecture 16 (verbatim signed-permutation/reversal/breakpoint/lower-bound), Tannier–Zheng–Sankoff
   PMC3887456 (adjacency vocabulary + d=n−sim + telomeres), Bafna–Pevzner 1998 (transposition vs inversion).
   Oracles: Hunter α=(−2,−3,+1,+6,−5,−4)→b=6/d≥3 with (−5,−4) excluded via (−y,−x); identity→0;
   single reversed block (+1,−4,−3,−2,+5)→b=2. Three source-backed ASSUMPTIONS (orthologMap-supplied
   anchors delegating anchor generation to the ORTHO/synteny units; strand '+'/'-'=sign; only
   Inversion/Transposition classified). Concise source page written; linked into the
   algorithm-validation-evidence hub (frontmatter sources + source-list + concept-list); added a reciprocal
   "two formulations" note on synteny-and-rearrangement-detection; index updated (source + concept lines).
   Contradictions: none among sources; Deviations = None beyond the three scoping assumptions.
   graph: +2 nodes, +2 typed edges (relates_to test-unit-registry, alternative_to synteny-and-rearrangement-detection)

## [2026-07-09] ingest | docs/Evidence/COMPGEN-REVERSAL-001-Evidence.md → compgen-reversal-001-evidence (source)
   Thirty-sixth per-algorithm Evidence file; eighth of the Comparative-genomics (COMPGEN) family
   (after ANI, CLUSTER, COMPARE, DOTPLOT, ORTHO, RBH, REARR). NO new concept — this file is NOT a
   distinct algorithm: CalculateReversalDistance computes the SAME reversal-distance lower bound
   d≥b/2 already documented in genome-rearrangement-breakpoint-distance (COMPGEN-REARR-001), only on
   UNSIGNED gene-order indices (breakpoint = |π_{i+1}−π_i|≠1 on the extended perm, the magnitude
   specialization of the signed criterion per Bafna–Pevzner §2 / Hübotter 2020) and returning the
   tightest integer ⌈b/2⌉=(b+1)/2 rather than a raw breakpoint count. It is a lower bound, NOT the
   exact distance (no Hannenhalli–Pevzner cycle/hurdle refinement). Reused that concept: added a new
   "Unsigned reversal distance ⌈b/2⌉" subsection with the unsigned breakpoint rule + oracle table,
   added COMPGEN-REVERSAL to its frontmatter sources (source_commit bumped to c6c3b01), and added a
   relates_to test-unit-registry typed edge for the new source. Sources: Bafna–Pevzner 1998 (rank 1,
   breakpoint def + lower-bound construction, identity=only 0-breakpoint perm), Hunter Lecture 16
   (rank 2, reversal removes ≤2, d≥b/2, symmetry d_β(α)=d_α(β)), Hübotter 2020 survey (rank 4,
   unsigned |Δ|≠1 specialization, corroboration only), Bergeron–Mixtacki–Stoye 2009 (rank 1,
   adjacency-vs-breakpoint). Oracles (unsigned): [2,3,1,6,5,4]→b=4→2, [4,3,2,1]→b=2→1, identity→0.
   Two source-backed ASSUMPTIONS: integer ⌈b/2⌉ rounding (tightest integer the theorem guarantees),
   unequal-length inputs throw ArgumentException (distance defined only within one marker set).
   Concise source page written; linked into the algorithm-validation-evidence hub (frontmatter sources
   + source-list) and the genome-rearrangement-breakpoint-distance concept; index updated (source line
   + concept-entry now marks COMPGEN-REARR-001 + COMPGEN-REVERSAL-001). Contradictions: none —
   signed/unsigned criteria are the same rule under |Δ|; Deviations = None beyond the two assumptions.
   graph: +1 node, +1 typed edge (relates_to test-unit-registry)

## [2026-07-09] ingest | docs/Evidence/COMPGEN-SYNTENY-001-Evidence.md → compgen-synteny-001-evidence (source)
   Thirty-seventh per-algorithm Evidence file; ninth of the Comparative-genomics (COMPGEN) family
   (after ANI, CLUSTER, COMPARE, DOTPLOT, ORTHO, RBH, REARR, REVERSAL). NO new concept — REUSED the
   existing shared synteny anchor synteny-and-rearrangement-detection (created for CHROM-SYNT-001,
   which was pre-named as the anchor COMPGEN-SYNTENY would reuse). This file is the comparative-genomics
   whole-genome counterpart and supplies the concrete MCScanX collinearity DP scoring parameters behind
   that anchor's FindSyntenyBlocks: chain DP Score(v)=max(MatchScore(v), max_u[Score(u)+MatchScore(v)+
   GapPenalty×NumberofGaps(u,v)]), MatchScore 50 / GapPenalty −1 / MAX_GAPS 25, report non-overlapping
   chains scoring over 250 (≥5 collinear anchor pairs), both transcriptional directions → forward +
   inverted (IsInverted) blocks, anchors from BLASTP E≤1e-5 with <5-gene collapse (generation delegated
   to COMPGEN-ORTHO-001). Enriched the concept with a new "MCScanX collinearity DP model" section + the
   scoring block, and updated its intro to record dual validation (CHROM-SYNT-001 chromosome-scale +
   COMPGEN-SYNTENY-001 whole-genome). Sources: MCScanX (Wang 2012, PMC3326336, rank 1, verbatim
   recurrence/params) + MCScanX Oxford HTML (synteny-vs-collinearity, anchors=homologs) + Wikipedia
   Synteny (rank 4, definitions). Oracles: 5 forward anchors→score 250→forward block; reversed order→
   inverted block; 4 anchors→score 200→no block; ≥25-gene gap→chain breaks; empty→no blocks. Two
   source-backed ASSUMPTIONS: report rule ≥250 AND ≥5 pairs (resolves the "over 250" vs "≥5 pairs"
   wording tension in favour of the explicit 5-pair minimum); anchors supplied as an orthologMap.
   Concise source page written; linked into the algorithm-validation-evidence hub (frontmatter sources +
   source-list + concept-list); index updated (source line). Contradictions: none — the two MCScanX
   renderings and Wikipedia agree; Deviations = None beyond the two scoping assumptions.
   graph: +1 node, +1 typed edge (relates_to test-unit-registry)

## [2026-07-09] ingest | docs/Evidence/DISORDER-LC-001-Evidence.md → disorder-lc-001-evidence (source) + 1 concept
   Thirty-eighth per-algorithm Evidence file; FIRST of the protein disorder / features family
   (DISORDER-LC / MORF / PRED / PROPENSITY / REGION). Created the genuinely-distinct concept
   protein-low-complexity-seg — the anchor for the protein-disorder/features family: the SEG algorithm
   (Wootton & Federhen 1993/1996) partitioning a protein into low- and high-complexity segments.
   Complexity = Shannon entropy H=−Σpᵢ·log₂pᵢ in bits/residue (max log₂20≈4.322), matching NCBI
   `blast_seg.c` `s_Entropy`; two-stage scan with three parameters W=12 (trigger window) / K1=2.2
   (trigger/locut cutoff) / K2=2.5 (extension/hicut cutoff), all verbatim NCBI/GCG defaults: stage-1
   triggers windows with entropy ≤ K1, stage-2 extends while ≤ K2. Judged genuinely distinct — SEG
   low-complexity is a different algorithm from intrinsic-disorder prediction (TOP-IDP) / MoRFs, so
   PRED/PROPENSITY/MORF/REGION are expected to warrant their own concept(s); and it is the PROTEIN
   counterpart of the genomic-DNA low-complexity under repetitive-element-detection (different alphabet
   + complexity measure), so I did NOT fold it into that repeats anchor. Hand-derived oracle window
   entropies (L=12): QQ..→0.0 triggers, AAAAAALLLLLL→1.0 triggers, AAABBBCCCDDD→2.0 triggers@K1=2.2
   (not strict 0.5), ACDEFGHIKLMN 12-distinct→3.584963>K2 no segment; corner cases seq<W→empty,
   homopolymer≥W→one full-span segment, all-distinct→none. Concise source page for the artifact (NCBI
   `blast_seg.c` rank-3 reference impl + GCG/Weizmann SEG help & `ncbi-seg` manpage rank-3 + Wootton &
   Federhen 1993 C&C 17(2):149-163 / 1996 Meth.Enzymol. 266:554-571 rank-1 primary). Linked new source
   + concept into the algorithm-validation-evidence hub (frontmatter sources + both link lists) and
   updated the index (source + concept lines). TWO documented ASSUMPTIONs, both flagged as deviations
   from Wootton & Federhen but neither moving segment boundaries on the canonical cases: (1) region-type
   label string "X-rich"/"X/Y-rich" (dominant-residue >50% presentation extension — SEG defines only
   segment location, not a label); (2) greedy single-residue extension (grow contig one residue at a
   time while whole-segment entropy ≤ K2 vs the reference merge of length-W extension windows — identical
   boundaries on homopolymer/dipeptide oracles). Contradictions: none — the NCBI reference impl, the
   GCG/manpage program docs, and the Wootton & Federhen primary literature agree on W=12/K1=2.2/K2=2.5,
   the Shannon-entropy bits/residue measure, and the two-stage trigger/extend scan. Follow-up: remaining
   protein-disorder units (MORF/PRED/PROPENSITY/REGION) warrant their own concept(s) — likely a shared
   intrinsic-disorder (TOP-IDP) anchor distinct from this low-complexity one — when ingested.
   graph: +2 nodes, +1 typed edge (relates_to test-unit-registry)

## [2026-07-09] ingest | docs/Evidence/DISORDER-MORF-001-Evidence.md → disorder-morf-001-evidence (source) + 1 concept
   Thirty-ninth per-algorithm Evidence file; SECOND of the protein disorder / features family
   (after DISORDER-LC-001 / SEG low-complexity). Created the genuinely-distinct concept
   morf-prediction-dip-in-disorder — MoRF (Molecular Recognition Feature) prediction by the "dip
   within disorder" heuristic. A MoRF = a short ordered segment embedded in a longer intrinsically
   disordered region that undergoes a disorder-to-order transition on partner binding. Criterion
   traced verbatim: Seqeron reports a MoRF where an ordered run (per-residue disorder score < 0.5,
   the PMC2570644 threshold) of TOTAL length within the Mohan 2006 10–70 residue band is flanked on
   BOTH sides by a disordered residue (score ≥ 0.5) inside a disordered region; terminal dips (not
   flanked both sides) excluded. Per-residue score from `PredictDisorder` = TOP-IDP scale (Campen
   2008) normalized `(raw+0.884)/1.871` to [0,1]; window averaging smooths boundaries. Mohan 2006
   α/β/ι bound-state sub-types recorded; MoRF score∈[0,1] rising with dip depth (bounded
   normalization = documented derivation, 0.5 threshold source-backed). Judged distinct from the SEG
   [[protein-low-complexity-seg]] sibling — SEG partitions by compositional complexity, MoRF reads a
   per-residue disorder profile for an ordered dip (the DISORDER-LC ingest had pre-flagged MORF as
   warranting its own concept). Oracle: synthetic ordered-L-dip in long P/E disordered flanks → one
   MoRF; corner cases fully-ordered/fully-disordered/out-of-10–70-band/terminal-dip → none. Concise
   source page for the artifact (Mohan 2006 J Mol Biol PMID 16935303 rank-1 + Cheng/Oldfield
   PMC2570644 rank-1 "dip" operational def + Oldfield 2005 Biochemistry PMID 16156658 rank-1 +
   Wikipedia rank-4; Campen 2008 TOP-IDP for the underlying score). Linked new source + concept into
   the algorithm-validation-evidence hub (frontmatter + both link lists) and updated the index
   (source + concept lines); added a reciprocal sibling nav link from protein-low-complexity-seg.
   ONE documented ASSUMPTION, scoped to the flank-length detail only: Oldfield 2005's exact numeric
   dip parameters (flank lengths, ordered-run window) are PAYWALLED and unretrievable, so the unit
   implements the fully-retrievable qualitative criterion — the 0.5 threshold, the 10–70 band, and
   the order-within-disorder shape are all source-traceable and NOT assumptions. Contradictions:
   none — Mohan/Cheng-Oldfield/Oldfield/Wikipedia agree on the 10–70 length, the short-order-within-
   longer-disorder shape, and the disorder-to-order transition. Follow-up: remaining protein-disorder
   units (PRED/PROPENSITY/REGION) warrant their own concept(s) — a shared intrinsic-disorder (TOP-IDP)
   `PredictDisorder` anchor is the likely next distinct concept.
   graph: +2 nodes, +2 typed edges (relates_to test-unit-registry, relates_to protein-low-complexity-seg)

## [2026-07-09] lint | structural + staleness + coverage + graph (89 pages)
   Structural: 1 orphan (readme) — fixed by adding a reciprocal [[readme]] content link
   from three-front-doors (README is that concept's source). Staleness: clean. Graph
   lint: clean (no typed-edge issues). Coverage: 1135 uncovered — the step-7 "report too
   long" signal; dominated by generated/reference material, not real gaps. Triaged with
   the user in batch:
   - Added a Coverage exclude policy to SCHEMA.md for generated subtrees (docs/mcp/tools/**
     427 tool docs, docs/Validation/reports/** per-run reports, docs/refactoring/**,
     docs/skills/_generated/**, docs/templates/**) → residual coverage drops to 405.
   - docs/algorithms/** kept IN scope. Reconciliation found step-1 already done: every
     existing algorithm concept already lists its docs/algorithms doc in sources: (commit
     9ce49ba, staleness-clean) → no frontmatter edits/source_commit bumps needed. 37 docs
     covered-via-concept, 208 pending (+2 index/map docs README/CANONICAL_MAP noted).
   - Created wiki/backlog.md (type: index): covered→concept table, pending grouped by 34
     domains with expected slugs (Oncology 37 / MolTools 17 / RnaStructure 13 / … the
     largest gaps), and the three approved queued source batches (10 testing checklists,
     4 Validation ledgers, 3 MCP top-level docs incl. MCP_STATUS.md). Pending list folds
     into the main per-algorithm ingest campaign, not a separate one. Registered under a
     new Meta section in index.md.
   - Known residual lint noise, accepted by design: backlog.md reads as an orphan (index-
     only inbound link, same quirk as readme pre-fix — it is an index-discoverable meta
     page) and trips the 400-line soft cap (488 lines; a flat reference table, well under
     the 800 hard cap, not worth sharding).
   Cadence note: coverage should always be run with the SCHEMA exclude globs from now on;
   without them every lint re-surfaces the ~693-doc generated long-tail.

## [2026-07-09] ingest | docs/Evidence/DISORDER-PRED-001-Evidence.md → disorder-pred-001-evidence (source) + 1 concept
   Fortieth per-algorithm Evidence file; THIRD of the protein disorder / features family (after
   DISORDER-LC-001 / SEG and DISORDER-MORF-001 / MoRF). Created the genuinely-distinct SHARED
   anchor concept intrinsic-disorder-prediction-top-idp — the TOP-IDP `PredictDisorder`
   sliding-window intrinsic-disorder predictor that MoRF prediction + disordered-region detection
   sit on (the "likely next distinct concept" flagged by the MORF ingest follow-up; NOT previously
   represented — only referenced inline on the morf page). Per-residue Sᵢ = mean over a window
   (default 21, edge-clipped) of min-max-normalized TOP-IDP `(p−(−0.884))/1.871 ∈ [0,1]`, disordered
   when Sᵢ ≥ 0.542 (Campen 2008 maximum-likelihood cutoff); Dunker 2001 disorder{A,R,G,Q,S,P,E,K} /
   order{W,C,F,I,Y,V,L,N} / ambiguous{D,H,M,T} sets; Kyte-Doolittle `CalculateHydropathy` utility;
   W30→0.0 / P30→1.0 / E30→0.866 homopolymer oracles; null/empty→empty, non-canonical residues
   skipped (all-unknown window→0.0), case-insensitive, O(n·w). Sources: Campen 2008 (PMC2676888
   rank-1 primary, TOP-IDP Table 2 + cutoff) + Dunker 2001 + Kyte-Doolittle 1982 + Wikipedia rank-3;
   concept also lists docs/algorithms/ProteinPred/Disorder_Prediction.md as a second source (SEG
   precedent; aids backlog reconciliation). Cross-linked: MoRF concept now `depends_on` this anchor
   (replaced its inline PredictDisorder prose with a wikilink) and SEG concept + index updated to
   point at the now-existing concept; linked into the algorithm-validation-evidence hub (frontmatter
   + both link lists) and index (source + concept lines). Assumptions: None in the evidence file —
   every parameter (TOP-IDP scale, 0.542 cutoff, Dunker sets, hydropathy) is source-traced;
   implementation-side note only = explicitly SIMPLIFIED single-feature TOP-IDP heuristic, not
   competitive with IUPred2A/MobiDB-lite, non-canonical residues skipped, edge windows clipped.
   Contradictions: none — the PRED unit's 0.542 max-likelihood cutoff vs the MoRF unit's 0.5
   order/disorder threshold (PMC2570644) are different published thresholds for different purposes,
   flagged as such, NOT a contradiction. Follow-up: remaining protein-disorder units (PROPENSITY /
   REGION — Disorder_Propensity.md, Disordered_Region_Detection.md) can reuse this anchor;
   Disordered_Region_Detection is the aggregation layer over PredictDisorder's per-residue calls.
   graph: +2 nodes, +3 typed edges (concept relates_to test-unit-registry + relates_to protein-low-complexity-seg; morf depends_on intrinsic-disorder-prediction-top-idp)

## [2026-07-09] ingest | docs/Evidence/DISORDER-PROPENSITY-001-Evidence.md → disorder-propensity-001-evidence (source)
   Fourth protein-disorder-family unit. NOT a new algorithm — the RAW per-residue propensity
   primitive layer beneath PredictDisorder, so REUSED/enriched the existing anchor concept
   [[intrinsic-disorder-prediction-top-idp]] rather than creating a new one (as the DISORDER-PRED
   log entry anticipated). Four in-scope methods: GetDisorderPropensity (returns the RAW
   un-normalized TOP-IDP Table-2 value W−0.884…P+0.987 — explicitly NOT the [0,1] normalized Sᵢ the
   windowed predictor uses; unknown residue→0.0 via GetValueOrDefault; case-folded), IsDisorderPromoting,
   and the two public sets DisorderPromotingAminoAcids={A,E,G,K,P,Q,R,S} / OrderPromotingAminoAcids=
   {C,F,I,L,N,V,W,Y} (with ambiguous {D,H,M,T} in neither; three sets disjoint + cover all 20 = 8+8+4).
   Sources: Campen 2008 (PMC2676888 rank-1, Table 2 raw values + anchors) + Wikipedia IDP (rank-4, for
   the Dunker 2001 classification it cites) + Dunker 2001 PubMed (rank-1 citation locator). Enriched
   the concept with a new "Per-residue propensity primitives" section (raw-vs-normalized value-space
   caveat), added the new evidence path to its frontmatter sources + bumped source_commit, and
   cross-linked from the intro. Updated index (new source line + concept line now names both units).
   Assumptions: two, both implementation-side, not source contradictions — (1) unknown-residue
   propensity 0.0 = GetValueOrDefault contract not a source value; (2) ranking-string vs Table-2-value
   S/K discrepancy (rendered "…Q,K,S,E,P" but S=0.341<K=0.586 → by value "…Q,S,K,E,P"; numeric values
   authoritative, no correctness impact). The 0.542 cutoff is recorded in-source for completeness but
   governs PredictDisorder, NOT this unit's four methods. Contradictions: none. Follow-up: only the
   REGION unit (Disordered_Region_Detection.md, the contiguous-run aggregation layer over
   PredictDisorder) remains in the protein-disorder family.
   graph: +1 node (source), +0 typed edges (concept↔registry / ↔SEG edges already exist from the
   DISORDER-PRED slug; source-page mentions auto-derived)

## [2026-07-09] ingest | docs/Evidence/DISORDER-REGION-001-Evidence.md → disorder-region-001-evidence (source)
   Fifth (and final anticipated) protein-disorder-family unit. NOT a new algorithm — the contiguous-run
   AGGREGATION + region-classification layer over the per-residue PredictDisorder TOP-IDP profile, so
   REUSED/enriched the existing anchor concept [[intrinsic-disorder-prediction-top-idp]] rather than
   creating a new one (exactly the growth the DISORDER-PRED/PROPENSITY log entries anticipated). Added a
   new "Disordered-region detection (DISORDER-REGION-001)" section to the anchor: contiguous run ≥
   minRegionLength(5) with MeanScore + boundary oracles (empty/all-ordered→no regions, all-disordered→one
   region, isolated<minLen→none, trailing region no off-by-one), and a two-scheme classification table.
   Default RegionType heuristic = dominant residue fraction > 0.25 → Proline-rich/Acidic/Basic/Ser-Thr-rich
   else Long IDR(len>30, van der Lee)/Standard IDR; opt-in ClassifyRegionFlavorMobiDbLite (Necci 2020 v3
   source-exact) = charge classes PA/PPE/NPE/WeaklyCharged via Das&Pappu FCR/NCPR at 0.35, then composition
   Cys→Pro→Gly→SEG→Polar{S,T,N,Q} at ≥0.32 inclusive, 9-residue window, sub-region≥9; boundaries + rescaled
   Confidence∈[0,1] unchanged when flavor scheme used (MobiDB-lite defines no per-residue confidence).
   Sources: Campen 2008 (rank-1, scale/cutoff/region idea) + Dunker 2001 (rank-1, long-IDR>30 significance)
   + van der Lee 2014 (rank-1, proline/acidic/basic/Ser-Thr subtypes + short-vs-long split) + Necci 2020
   MobiDB-lite 3.0 (rank-1 paper + rank-3 version-pinned v3 impl) + Wikipedia (rank-4). Added the new
   evidence path to the anchor frontmatter sources + bumped source_commit to 98b44f1a; cross-linked the
   intro; updated index (new source line + anchor concept line now names all three PRED/PROPENSITY/REGION
   units + region-detection layer). CONTRADICTION flagged (in-source ref #6): the default RegionType 0.25
   composition cutoff is an internal ~5×-random heuristic, NOT Das&Pappu 2013's 0.25 — that value is NCPR,
   a globule/coil conformational-state threshold, not a compositional-enrichment threshold; MobiDB-lite's
   own 0.35/0.32 thresholds ARE source-exact. Follow-up: none — the protein-disorder family (LC/MoRF/PRED/
   PROPENSITY/REGION) is now fully ingested.
   graph: +1 node (source), +0 typed edges (concept↔registry / ↔SEG edges already exist on the anchor;
   source-page [[wikilinks]] mentions auto-derived)

## [2026-07-09] ingest | docs/Evidence/EPIGEN-AGE-001-Evidence.md → epigen-age-001-evidence (source) + 1 concept
   FIRST unit of the NEW Epigenetics family (EPIGEN-AGE-001) — epigenetic ("DNAm") age via the Horvath
   2013 multi-tissue DNA-methylation clock. Surveyed wiki/concepts: NO pre-existing epigenetics /
   methylation / CpG concept to reuse (grep hits in log/backlog/registry were incidental), so created a
   genuinely-distinct anchor concept [[epigenetic-age-horvath-clock]] — the two-stage scorer: linear
   predictor Y = intercept + Σ coef_i·β_i over 353 elastic-net-selected clock CpGs (transformed-age
   units), then the two-branch Horvath inverse calibration `anti.trafo` (adult.age=20): 21·exp(Y)−1 for
   Y<0, 21·Y+20 for Y≥0, continuous at (0,20)=age 20. Embedded intercept 0.695507258 + 353
   `CoefficientTraining` weights (Additional file 3), cross-verified byte-identical vs an independent
   GitHub mirror (all 353 pairs). Invariants Y=0→20 / CpGs-absent-ignored / empty-map→F⁻¹(intercept);
   null map|coeffs→ArgumentNullException, empty coeffs→ArgumentException. Oracles: anti.trafo 0→20 / 1→41
   / −1→6.7254682646002895 / −2.5→0.7237849711018749, linear-assembly Y=0.684247258→34.369192418.
   Sources: Horvath 2013 paper (Genome Biology 14:R115, rank-1) + aldringsvitenskap reference R
   `horvath2013.R`/`StepwiseAnalysis.R` (rank-3, trafo/anti.trafo + predictedAge assembly) + Additional
   file 3 Springer supplement (rank-1) + its byte-identical GitHub mirror (rank-3, the cross-check) +
   perishky/meffonym (rank-3). Concept also lists docs/algorithms/Epigenetics/Epigenetic_Age_Estimation.md
   as a second source (backlog reconciliation). Reconciled backlog: moved Epigenetic_Age_Estimation.md
   from pending Epigenetics (6→5) to covered-via-concept (37→38 covered / 208→207 pending). Linked into
   the algorithm-validation-evidence hub (frontmatter sources + both link lists) and index (source +
   concept lines). Scope note: the Evidence unit covers ONLY the multi-tissue 353-CpG clock; the
   algorithm doc (reviewed 2026-06-23, one day later) has since extended the impl with Horvath-2018
   skin&blood (391 CpG, same anti.trafo) + Levine-2018 PhenoAge (513 CpG, NO transform) — recorded as an
   evolution, NOT a contradiction. Assumptions: table-embed assumption RESOLVED 2026-06-22
   (caller-supplied → embedded + cross-verified); no correctness gap remains. Contradictions: none.
   Follow-up: the rest of the Epigenetics family (CpG detection, methylation calling, DMRs, chromatin
   state, bisulfite) remains pending — those are distinct algorithms and will get their own concepts.
   graph: +2 nodes (source + concept), +1 typed edge (concept relates_to test-unit-registry;
   source-page [[wikilinks]] mentions auto-derived)

## [2026-07-09] ingest | docs/Evidence/EPIGEN-BISULF-001-Evidence.md → epigen-bisulf-001-evidence (source) + 1 concept
   Created bisulfite-methylation-calling concept (second Epigenetics unit: SimulateBisulfiteConversion
   Frommer 1992 + CalculateMethylationFromBisulfite Bismark meth/(meth+unmeth) + GenerateMethylationProfile
   Schultz 2012 weighted level). Cross-linked sibling epigenetic-age-horvath-clock (produces the β-values it
   consumes). Updated algorithm-validation-evidence (source list + mention), index, backlog (Bisulfite_Sequencing_Analysis
   moved to covered: 39/206). No contradictions.
   graph: +2 nodes, +1 typed edge

## [2026-07-09] ingest | docs/Evidence/EPIGEN-CHROM-001-Evidence.md → epigen-chrom-001-evidence (source) + 1 concept
   Third unit of the Epigenetics family (EPIGEN-CHROM-001) — ChromHMM-style chromatin state prediction
   from histone modification marks. Surveyed wiki/concepts: NO pre-existing chromatin / histone /
   ChromHMM concept (grep hits in epigenetic-age/centromere/log were incidental), and this is a
   genuinely-distinct algorithm from the two methylation siblings (histone ChIP-seq axis, not DNA
   methylation), so created a new anchor concept [[chromatin-state-prediction]]. Core = the ChromHMM
   binary present/absent mark model (Ernst & Kellis 2012 multivariate HMM; BinarizeBed/BinarizeBam →
   LearnModel operates on 1/0 calls; state = function of the SET of present marks, magnitude beyond the
   call irrelevant = binary invariance). `PredictChromatinState` takes the six Roadmap 18-state marks
   {H3K4me3, H3K4me1, H3K27ac, H3K36me3, H3K27me3, H3K9me3} as [0,1] signals, present > threshold
   (default 0.5), maps the pattern → ActivePromoter(TssA)/ActiveEnhancer/WeakEnhancer/Transcribed(Tx)/
   Repressed(ReprPC)/Heterochromatin(Het)/BivalentPromoter(TssBiv)/BivalentEnhancer(EnhBiv)/
   LowSignal(Quies). Two combinatorial rules captured: bivalency (H3K4me3+H3K27me3) is a state not a
   contradiction, and promoter(H3K4me3) dominates enhancer(H3K4me1) at one locus. Companions
   `AnnotateHistoneModifications` (per-mark region label) + `FindAccessibleRegions` (ATAC-like
   contiguous-above-threshold merge + minWidth exclusion). Sources: Ernst & Kellis 2012 (Nat Methods,
   rank-1) + ChromHMM manual (rank-3, binarization) + Roadmap Epigenomics 15/18-state definitions
   (rank-2) + six per-mark Wikipedia primaries (Liang 2004 H3K4me3 / Rada-Iglesias 2018 H3K4me1 /
   Creyghton 2010 H3K27ac / Ferrari 2014 H3K27me3 / Nicetto 2019 H3K9me3 / Kimura 2013 H3K36me3).
   Concept also lists docs/algorithms/Epigenetics/Chromatin_State_Prediction.md as a second source.
   Reconciled backlog: moved Chromatin_State_Prediction.md from pending Epigenetics (4→3) to
   covered-via-concept (39→40 covered / 206→205 pending). Cross-linked both methylation siblings
   ([[epigenetic-age-horvath-clock]], [[bisulfite-methylation-calling]]) bidirectionally; linked into
   the algorithm-validation-evidence hub (frontmatter sources + both link lists) and index (source +
   concept lines). Two assumptions: presence-threshold value 0.5 (ChromHMM uses a Poisson background
   from raw counts, not a fixed [0,1] cut — tests choose unambiguous magnitudes) and single-locus
   promoter-over-enhancer precedence (Roadmap derives it from spatial HMM context). Research-grade: the
   state-assignment logic is fully source-backed but this is NOT a trained HMM (no LearnModel / Poisson
   binarization / spatial context). Contradictions: none.
   graph: +2 nodes, +1 typed edge (concept relates_to test-unit-registry; source-page mentions auto-derived)

## [2026-07-09] ingest | docs/Evidence/EPIGEN-CPG-001-Evidence.md → epigen-cpg-001-evidence (source) + 1 concept
   FOURTH unit of the Epigenetics family (EPIGEN-CPG-001) — CpG site detection + canonical CpG
   observed/expected ratio + CpG-island detection. Surveyed wiki/concepts: NO pre-existing CpG /
   CpG-island concept (the sibling [[bisulfite-methylation-calling]] only *reuses* `FindCpGSites` inline),
   so created a genuinely-distinct anchor concept [[cpg-island-detection]] — a SEQUENCE-ONLY algorithm
   touching no methylation state, distinct from all three siblings (bisulfite measures state, Horvath clock
   scores age from β-values, chromatin works on histone marks). Three ops on `EpigeneticsAnalyzer`:
   `FindCpGSites` (linear O(n) scan, 0-based C position, adjacent `CGCG`→2 distinct sites; the same call
   [[bisulfite-methylation-calling]] reuses); `CalculateCpGObservedExpected` (Gardiner-Garden & Frommer
   O/E = CpG_count/((C·G)/L), UCSC-standard, div-by-zero guard → 0 when no C/G or L<2); `FindCpGIslands`
   (sliding-window merge, default ≥200 bp / GC ≥0.5 / O-E ≥0.6 INCLUSIVE, 0-based inclusive-Start /
   exclusive-End tuple, O(n·w) rescan). Oracles: CGCG×20→O/E 2.0, ACGTCGACG→3.0, ACGT→4.0, AT-only→0.0,
   400 bp CGCG→1 island. Corner cases: GpC≠CpG, case-insensitive uppercase-normalize, length-1→0 sites,
   zero C/G→O/E 0. Sources: Gardiner-Garden & Frommer 1987 (J Mol Biol, rank-1, canonical criteria +
   formula) + Takai & Jones 2002 (PNAS, rank-1, stricter ≥500/55%/0.65 + confirms the ≥ inclusive
   operators) + Saxonov 2006 (PNAS, rank-1, alt expected ((C+G)/2)²/L) + Wikipedia CpG site (rank-4).
   Concept also lists docs/algorithms/Epigenetics/CpG_Site_Detection.md as a second source (backlog
   reconciliation). Reconciled backlog: moved CpG_Site_Detection.md from pending Epigenetics (3→2) to
   covered-via-concept (40→41 covered / 205→204 pending); updated the index backlog-summary counts.
   Cross-linked all three Epigenetics siblings bidirectionally ([[bisulfite-methylation-calling]] tightest
   — shares `FindCpGSites`; [[epigenetic-age-horvath-clock]]; [[chromatin-state-prediction]]); linked into
   the algorithm-validation-evidence hub (frontmatter sources + source-list + distinct-concept list) and
   index (source + concept lines). Takai-Jones stricter + Saxonov alt-formula recorded as NOT-preset
   (supplied via custom args), an impl scoping decision not a deviation. Assumptions: None (all behaviour
   source-defined per the Evidence file). Contradictions: none.
   graph: +2 nodes (source + concept), +1 typed edge (concept relates_to test-unit-registry;
   source-page [[wikilinks]] mentions auto-derived)

## [2026-07-09] ingest | docs/Evidence/EPIGEN-DMR-001-Evidence.md → epigen-dmr-001-evidence (source) + 1 concept
   FIFTH unit of the Epigenetics family (EPIGEN-DMR-001) — differentially methylated region (DMR)
   detection by the methylKit tiling-window + Fisher's-exact-test model. Surveyed wiki/concepts: NO
   pre-existing DMR / methylation-comparison concept, and this is a genuinely-distinct algorithm —
   it COMPARES methylation between two samples (vs the single-sample siblings), consuming the per-CpG
   C/(C+T) β-values [[bisulfite-methylation-calling]] produces over the CpGs [[cpg-island-detection]]
   locates. Created anchor concept [[differentially-methylated-regions]]. Core: fixed tiling windows
   (`win.size=1000/step.size=1000/cov.bases=0`, tileMethylCounts), meth.diff=group2%−group1%
   (percentage points), per-window pooled 2×2 methylated/unmethylated counts → Fisher's exact test
   (one sample per group; replicates→logistic regression out of scope), hypergeometric single-table p
   `(a+b)!(c+d)!(a+c)!(b+d)!/(a!b!c!d!n!)` + two-sided sum-of-≤-observed; getMethylDiff reports
   q<0.01 AND |meth.diff|>25 STRICT (hyper meth.diff>+25 / hypo <−25). Corner cases: empty→no DMRs,
   zero-coverage group / zero marginal→degenerate 2×2→p=1.0 (not reported), |meth.diff|=25→excluded.
   Oracles: Fisher single-table a=1/b=9/c=11/d=3,n=24→≈0.001346076 (Wikipedia studying-by-gender);
   hyper window g1 level 0.0/cov 20 + g2 level 1.0/cov 20 × 3 sites → pooled meth {0,60}/unmeth {60,0}
   → meth.diff +100, p≈0 (complete separation), Hypermethylated. Sources: Akalin 2012 methylKit
   (Genome Biology 13:R87, PMC3491415, rank-1/3) + tileMethylCounts/calculateDiffMeth man pages +
   get.methylDiff diffMeth.R source (rank-3) + Fisher's exact test Wikipedia citing Fisher 1922/1935
   (rank-4, primary math). Concept also lists docs/algorithms/Epigenetics/Differentially_Methylated_Regions.md
   as a second source. Reconciled backlog: moved Differentially_Methylated_Regions.md from pending
   Epigenetics (2→1) to covered-via-concept (41→42 covered / 204→203 pending); updated the index
   backlog-summary counts. Cross-linked all four Epigenetics siblings ([[bisulfite-methylation-calling]]
   tightest — DMR consumes its β-values, added reciprocal link; [[cpg-island-detection]];
   [[epigenetic-age-horvath-clock]]; [[chromatin-state-prediction]]); linked into the
   algorithm-validation-evidence hub (frontmatter sources + source-list + distinct-concept list) and
   index (source + concept lines). Two evidence-backed assumptions (per-window pooling into one 2×2 =
   tileMethylCounts→Fisher documented pipeline; numC/numT=round(level×coverage) count reconstruction
   from the fractional-level MethylationSite). Research-grade: no logistic-regression replicate path /
   SLIM q-value model / overlapping-window tiling. Contradictions: none — Akalin 2012, the methylKit
   reference, and the Fisher hypergeometric definition are mutually consistent.
   graph: +2 nodes (source + concept), +1 typed edge (concept relates_to test-unit-registry;
   source-page [[wikilinks]] mentions auto-derived)

## [2026-07-09] ingest | docs/Evidence/EPIGEN-METHYL-001-Evidence.md → epigen-methyl-001-evidence (source) + methylation-context-classification (concept)
   Sixth and FINAL unit of the Epigenetics family — completes it. New source page epigen-methyl-001-evidence
   + new concept methylation-context-classification (CpG/CHG/CHH sequence-context classifier). Kept as its
   own concept (not folded into bisulfite): the distinct, wiki-worthy method is the IUPAC H="not G"
   trinucleotide classification of NON-CpG cytosines (CHG/CHH), which [[bisulfite-methylation-calling]]
   explicitly does not call from reads. The shared GenerateMethylationProfile (Schultz 2012 weighted level)
   is documented on bisulfite and only referenced here (no duplication). Sources: Cornish-Bowden 1985 IUPAC
   (H=A/C/T), Krueger-Andrews 2011 Bismark (CpG/CHG/CHH discrimination, CpG/CHG symmetric vs CHH asymmetric),
   Lister 2009 (non-CG mCHG/mCHH prevalence, IMR90 99.98% CG vs H1 ES ~25% non-CG), Schultz 2012 (weighted
   level). Oracles CGACAGCAA→CpG@0/CHG@3/CHH@6 + weighted (8+2)/(10+10)=0.5. Two assumptions (sequence-only
   level=0 placeholder, 0.5 MethylatedCpG count cutoff descriptive-only). Cross-linked all Epigenetics
   siblings (bisulfite tightest — added reciprocal cross-links on bisulfite's intro + not-implemented note);
   wired into algorithm-validation-evidence hub (frontmatter sources + source-list + distinct-concept list)
   and index (source + concept lines). Reconciled backlog: moved Methylation_Analysis.md from pending
   Epigenetics (1→0, section removed) to covered-via-concept (42→43 covered / 203→202 pending, 34→33 domains);
   Epigenetics family now 6/6 covered. Contradictions: none — IUPAC, Bismark, Lister 2009, Schultz 2012 mutually consistent.
   graph: +2 nodes (source + concept), +1 typed edge (concept relates_to test-unit-registry; source-page [[wikilinks]] mentions auto-derived)

## [2026-07-09] ingest | GENOMIC-COMMON-001-Evidence.md → genomic-common-001-evidence (source) + longest-common-substring (concept)
   Ingested the Longest Common Substring / Common Region Detection evidence unit (generalized suffix
   tree). Created source page genomic-common-001-evidence and NEW concept longest-common-substring —
   genuinely distinct (a string/suffix-tree algorithm, no prior LCS concept). Synthesized: LCS = a
   longest *contiguous* substring of both strings (contiguity is THE distinction from the gapped
   longest common *subsequence*); generalized-suffix-tree solution (LCS = path label to the deepest
   internal node whose subtree holds leaves from both strings, Θ(n+m)/O(M+N) build+DFS, Gusfield 1997).
   API contract: FindLongestCommonRegion (0-based positions in both) + FindCommonRegions(minLength);
   CommonRegion.None = empty/len 0/positions −1, identical→whole seq at 0/0. One assumption =
   documented deterministic first-found-in-`other` tie-break (Wikipedia reports all, GeeksforGeeks
   returns one); DNA oracle CACAGAG/TACATAGAT ties ACA vs AGA → selects ACA. Sources Wikipedia
   "Longest common substring" + GeeksforGeeks "Suffix Tree Application 5". Cross-linked
   [[dot-plot-word-match]] (shared generalized-suffix-tree exact-match engine, reciprocal link added).
   Wired into algorithm-validation-evidence hub (frontmatter sources + source-list + distinct-concept
   list) and index (source + concept lines). Reconciled backlog: moved Common_Region_Detection.md from
   pending Sequence_Comparison (1→0, section removed) to covered-via-concept (43→44 covered / 202→201
   pending, 33→32 domains). Contradictions: none — Wikipedia and GeeksforGeeks agree on contiguity + GST
   mechanism, differ only on reporting ties (all vs one), resolved deterministically by the repo.
   graph: +2 nodes (source + concept), +2 typed edges (concept relates_to test-unit-registry + relates_to dot-plot-word-match; source-page [[wikilinks]] mentions auto-derived)

## [2026-07-09] ingest | docs/Evidence/GENOMIC-MOTIFS-001-Evidence.md → genomic-motifs-001-evidence (source) + known-motif-search (concept)
   New source page [[genomic-motifs-001-evidence]] + new concept [[known-motif-search]] — the
   Motif-Analysis "Known Motif Search" unit: multi-pattern EXACT substring matching of a set of
   known query motifs (`GenomicAnalyzer.FindMotif`), the exact-equality baseline distinct from motif
   discovery and degenerate (IUPAC/PROSITE/PWM) matching. THE correctness rule = all OVERLAPPING
   occurrences reported (AAA in AAAAA→{0,1,2}, per Gusfield/Tufts + Biopython `count_overlap` 3-not-2);
   API contract 0-based sorted positions, per-motif position lists (absent motifs omitted),
   upper-cased result keys, empty/whitespace motif→no entry (`Array.Empty<int>()`). Oracles: EcoRI
   GAATTC in GAATTCAAAGAATTC→{0,9}, `{ACGT,AA,TTT}` in ACGTACGTAA→{0,4}/{8}/omitted. Sources
   Tufts COMP 150GEN (Gusfield) + Biopython `Seq.search`/`count_overlap` + Wikipedia "Restriction
   site" (EcoRI). Two API-shape assumptions (empty-motif→no entry, upper-case keys), deviations None.
   Wired into algorithm-validation-evidence hub (frontmatter sources + source-list + distinct-concept
   list) and index (source + concept lines). Reconciled backlog: moved Known_Motif_Search.md from
   pending Motif_Analysis (1→0, section removed) to covered-via-concept (44→45 covered / 201→200
   pending, 32→31 domains). Contradictions: none — Gusfield and Biopython agree all occurrences
   including overlaps are reported.
   graph: +2 nodes (source + concept), +1 typed edge (concept relates_to test-unit-registry;
   source-page [[wikilinks]] mentions auto-derived)

- 2026-07-09 — ingest `docs/Evidence/GENOMIC-ORF-001-Evidence.md` (test unit GENOMIC-ORF-001,
  Open Reading Frame detection). Created source [[genomic-orf-001-evidence]] + NEW concept
  [[open-reading-frame-detection]] (no prior ORF/gene-prediction concept existed). Algorithm =
  `GenomicAnalyzer.FindOpenReadingFrames`: six-frame ATG→first-in-frame-stop enumeration, standard
  code (start ATG / stops TAA-TAG-TGA), reported span INCLUDES the stop (Length%3==0) while the
  translated protein EXCLUDES it; THE correctness rule = every in-frame ATG reaching a stop is
  reported so nested ORFs sharing a stop are both returned (Rosalind MGMTPRLGLESLLE/MTPRLGLESLLE),
  ATG-with-no-stop→none; minLength in NUCLEOTIDES inclusive (default 100), 0-based Position / Frame
  1–3 / IsReverseComplement, INV-01..05, O(n²) worst / O(n) typical. Oracles: Rosalind_99→4 distinct
  proteins, ATGAAAAAATAA→MKK (pos 0 / frame 1). Sources Rosalind + Wikipedia + NCBI ORFfinder +
  NCBI transl_table=1. Three source-anchored assumptions (stop-inclusive span / nt minLength /
  ATG-only), deviations None (one fixed pre-existing greedy bug). Scope-disambiguated from the
  annotation-layer `GenomeAnnotator.FindOrfs` (ANNOT-ORF-001, ATG/GTG/TTG + aa-length + strand/start
  flags) and `Translator.FindOrfs` (genetic-code-parameterized) — deliberately NOT contract-equivalent.
  Wired into algorithm-validation-evidence hub (frontmatter sources + source-list + distinct-concept
  list) and index (source + concept lines). Reconciled backlog: moved Open_Reading_Frame_Detection.md
  from pending Analysis (2→1) to covered-via-concept (45→46 covered / 200→199 pending, 31 domains).
  Contradictions: none — Rosalind/Wikipedia/NCBI agree on six-frame ATG→stop, distinct-protein return.
  graph: +2 nodes (source + concept), +1 typed edge (concept relates_to test-unit-registry;
  source-page [[wikilinks]] mentions auto-derived)

## [2026-07-09] ingest | docs/Evidence/GENOMIC-REPEAT-001-Evidence.md → genomic-repeat-001-evidence (source) + 1 concept
  Per-algorithm Evidence file. Created the genuinely-distinct concept longest-repeated-substring —
  the Repeat-Analysis LRS + all-repeats-enumeration unit (`GenomicAnalyzer.FindLongestRepeat` /
  `FindRepeats`): LRS = deepest internal node with ≥2 leaves in a *single-string* suffix tree
  (string depth = repeat length, CMU 15-451 §2.1 verbatim / Wikipedia / GeeksforGeeks App-3 /
  Gusfield 5.4 via JHU); FindRepeats enumerates every substring occurring ≥2× via sorted-suffix
  adjacent-LCP *every-prefix* expansion (O(n²), the FINDINGS_REGISTER short-prefix fix). Positioned
  as the one-string sibling of longest-common-substring (added a comparison table + reciprocal nav
  link on both pages) and explicitly distinguished from the tandem/inverted repetitive-element-detection
  anchor (§2.5 LRS-vs-FindTandemRepeats contrast). Oracles ATCGATCGA→ATCGA{0,4}, AAAAAAAAAA→AAAAAAAAA{0,1}
  overlap, ATATATA→ATATA{0,2}, ACGT/empty→None, ACGTACGTTTTTACGT@3→8-substring set. Corner cases:
  overlaps counted, minLength≤0→max(1,minLength), ACGT-only, no reverse-complement / maximal-repeat
  classification. Concise source page for the GENOMIC-REPEAT-001 artifact (four sources, LRS + brute-force
  enumeration datasets, tie-break + ascending-positions assumptions). Wired into algorithm-validation-evidence
  hub (frontmatter sources + source-list + distinct-concept list) and index (source + concept lines).
  Reconciled backlog: moved Repeat_Analysis/Repeat_Detection.md from pending Repeat_Analysis (6→5) to
  covered-via-concept (46→47 covered / 199→198 pending, 31 domains). Contradictions: none — all four
  sources agree on the deepest-internal-node characterisation and overlap allowance; deviations None.
  graph: +2 nodes (source + concept), +3 typed edges (concept relates_to test-unit-registry,
  relates_to longest-common-substring, relates_to repetitive-element-detection; source-page
  [[wikilinks]] mentions auto-derived)

## [2026-07-09] ingest | docs/Evidence/GENOMIC-SIMILARITY-001-Evidence.md → genomic-similarity-001-evidence (source) + kmer-jaccard-similarity (concept)
  Per-algorithm Evidence file. Created the genuinely-distinct concept kmer-jaccard-similarity — the
  Analysis family's alignment-free pairwise-similarity unit (`GenomicAnalyzer.CalculateSimilarity`):
  Jaccard index `|A∩B|/|A∪B|` over the two sequences' *distinct* k-mer sets (HashSet, within-sequence
  repeats collapse), exact (no MinHash sketch) J×100 in [0,100], O(n+m). Sources: Jaccard 1901 (index
  definition, [0,1] range, non-empty-set scope, distance 1−J) + Ondov 2016 *Mash* (k-mer-set Jaccard =
  fraction of shared k-mers, sketch estimate |A_s∩B_s|/s) + Mash distance docs. INV symmetry / identical→100 /
  disjoint→0 / distinct-set; k=3 oracles 80.0 / 100⁄3 / 100 / 0 / (AAAAAA vs AAAA→100). Three source-backed
  assumptions (empty-union→0.0 ASM-1 convention, ×100 scaling, default k=5) + suffix-tree-evaluated-not-used
  note. Positioned as `alternative_to` the positional [[alignment-statistics]] (§2.5 set-resemblance vs
  residue-by-residue), and cross-linked (body/mentions) to the 5-mer-Jaccard metric behind
  ortholog-detection-reciprocal-best-hits and the exact-set basis Mash sketches for average-nucleotide-identity.
  Concise source page for the artifact. Wired into index (source + concept lines). Reconciled backlog: moved
  Analysis/Sequence_Similarity.md from pending Analysis (section emptied, 1→0) to covered-via-concept
  (47→48 covered / 198→197 pending, 31 domains). Contradictions: none — Jaccard's set definition and Mash's
  k-mer-set application are consistent; deviations None.
  graph: +2 nodes (source + concept), +2 typed edges (concept relates_to test-unit-registry,
  alternative_to alignment-statistics; source-page [[wikilinks]] mentions auto-derived)

## [2026-07-09] ingest | docs/Evidence/GENOMIC-TANDEM-001-Evidence.md → genomic-tandem-001-evidence (source) + repetitive-element-detection (concept, enriched)
  Per-algorithm Evidence file validating `GenomicAnalyzer.FindTandemRepeats` (exact tandem-repeat
  detection). REUSED the existing repeats/tandem anchor [[repetitive-element-detection]] rather than
  creating a new page — GENOMIC-TANDEM-001 is a consolidated duplicate of REP-TANDEM-001 (same method,
  same brute-force scan, canonical fixture, no new tests), and tandem detection is already sub-problem #1
  of that concept. Enriched the concept's tandem section with the two entry points over the same
  exact-copy model: `GenomicAnalyzer.FindTandemRepeats` (reports EVERY unit-length/period interpretation,
  no primitive-unit canonicalization — `AAAA` → period 1×4 AND period 2×2) vs the annotation
  `RepeatAnalyzer` path (primitive-unit rule); both exact-only, neither reports Benson TRF's approximate
  copies (Framework/Simplified limitation). Sources: Benson 1999 (Tandem Repeats Finder, period/copy-number/
  k≥2 definition, approximate-vs-exact) + Wikipedia "Tandem repeat" (ATTCG×3 worked example, STR/mini/
  macrosatellite classes, ~8% genome / >50 diseases). Oracles ATTCGATTCGATTCG→ATTCG/period5/3copies/len15
  and ATGATGATG→ATG/3. Concise source page written. Wired into index (source line). Reconciled backlog:
  moved Genomic_Analysis/Tandem_Repeat_Detection.md from pending Genomic_Analysis (section emptied, 1→0)
  to covered-via-concept under [[repetitive-element-detection]] (48→49 covered / 197→196 pending, 30
  domains); Repeat_Analysis/Tandem_Repeat_Detection.md (REP-TANDEM-001, a separate unit) left pending.
  Contradictions: none among sources; the two entry points' period-handling divergence documented, not a
  source conflict.
  graph: +1 node (source page), +1 typed edge (concept relates_to test-unit-registry from
  genomic-tandem-001-evidence; source-page [[wikilinks]] mentions auto-derived)

## [2026-07-09] ingest | KMER-ASYNC-001-Evidence.md → kmer-async-001-evidence (source) + asynchronous-kmer-counting (concept)
   First K-mer family unit. KMER-ASYNC-001 validates the asynchronous k-mer count `KmerAnalyzer.CountKmersAsync`
   — the cooperatively cancelable, progress-reporting `Task.Run` wrapper over synchronous `CountKmers`
   (KMER-COUNT-001). Determined this is NOT a distinct counting algorithm: the numeric result is fixed by
   the k-mer formula L−k+1 and is identical to the sync reference; the uniquely validated content is the
   .NET cooperative-cancellation + progress contract (ThrowIfCancellationRequested→OperationCanceledException
   + Canceled state, pre-start Task.Run(func,token) cancellation, awaiting a canceled task throws,
   IProgress 0→1.0). No existing k-mer-counting concept to enrich (sync KMER-COUNT-001 not yet ingested), so
   created a focused concept `asynchronous-kmer-counting` (folds the execution contract + inherited count
   into one wiki-worthy page rather than a thin gotcha). Sources: Wikipedia K-mer (L−k+1 / nᵏ, ATGG→ATG+TGG,
   GTAGAGCTGT k=2/3/4 total 9/8/7 distinct 7/8/7) + Microsoft Learn Task Cancellation / Task.Run. One
   assumption = numeric contract identical to sync (non-correctness-affecting); not parallelized; suffix-tree
   evaluated-not-used. Wired into index (source + concept lines). Reconciled backlog: moved
   K-mer/Asynchronous_K-mer_Counting.md from pending K-mer (10→9) to covered-via-concept under
   [[asynchronous-kmer-counting]] (49→50 covered / 196→195 pending, 30 domains); the other 9 K-mer docs
   (incl. K-mer_Counting.md / KMER-COUNT-001) left pending. Contradictions: none (count definition and .NET
   cancellation contracts are orthogonal and mutually consistent).
   graph: +2 nodes (concept + source page), +1 typed edge (concept relates_to test-unit-registry from
   asynchronous-kmer-counting; source/concept [[wikilinks]] mentions auto-derived)

## [2026-07-09] ingest | KMER-BOTH-001-Evidence.md → kmer-both-001-evidence (source) + both-strand-kmer-counting (concept)
   Second K-mer family unit. KMER-BOTH-001 validates `KmerAnalyzer.CountKmersBothStrands` — additive
   strand-aware counting. Determined this IS a genuinely distinct method (not a thin wrapper like
   KMER-ASYNC): the ADDITIVE / kPAL-"balance" convention `count[w]=forward[w]+forward[RC(w)]` (count
   k-mers of S and of RC(S), sum per key), which keeps a key per observed k-mer — explicitly NOT the
   canonical-collapsing convention (lexicographically-smaller of {w,RC(w)} as one key) of Jellyfish `-C`
   / Mash, which the algorithm does not implement. Created concept `both-strand-kmer-counting` (core
   model, INV-01..05 incl. grand-total 2·(L−k+1) / strand-symmetry / palindrome-doubling, additive-vs-
   canonical table, contract, three oracles, complexity, deviations). Sources: kPAL Methodology + Anvar
   2014 (Genome Biology 15:555, balance = sum of k-mer and its RC) + Shporer 2016 (inversion symmetry,
   grounds INV-01) + Marçais-Kingsford 2011 Jellyfish (single-strand primitive + `-C` contrast) + Mash
   issue #45 (canonical def) + Clavijo 2018 (strand rationale). Oracles ATGGC k=2→{AT:2,TG:1,GG:1,GC:2,
   CC:1,CA:1}, palindromic ACGT→{AC:2,CG:2,GT:2}, AAA→{AA:2,TT:2}. Two API-shape assumptions (empty/k>L→
   empty dict, k≤0→ArgumentOutOfRangeException inherited from CountKmers); Deviations = None. Cross-linked
   with sibling [[asynchronous-kmer-counting]] (shared sync CountKmers primitive; linked its both-strand
   mention). Concept covers the Both_Strand algorithm doc too. Wired into index (source + concept lines)
   and the algorithm-validation-evidence hub (frontmatter + body source list + own-concept list; also
   back-filled the missing kmer-async-001-evidence hub-body link). Reconciled backlog: moved
   K-mer/Both_Strand_Kmer_Counting.md from pending K-mer (9→8) to covered-via-concept (50→51 covered /
   195→194 pending, 30 domains); the other 8 K-mer docs (incl. K-mer_Counting.md / KMER-COUNT-001) left
   pending. Contradictions: none — kPAL balance and inversion symmetry give identical additive semantics;
   canonical wording cited only to contrast the not-implemented collapsing mode.
   graph: +2 nodes (concept + source page), +2 typed edges (concept relates_to test-unit-registry +
   relates_to asynchronous-kmer-counting; source/concept [[wikilinks]] mentions auto-derived)

## [2026-07-09] ingest | docs/Evidence/KMER-DIST-001-Evidence.md → kmer-dist-001-evidence (source) + k-mer-euclidean-distance (concept)
   K-mer Euclidean distance (`KmerAnalyzer.KmerDistance`): alignment-free L2 distance over normalized
   k-mer FREQUENCY vectors f_s(w)=count/(L−k+1), summed over the union of observed k-mers (absent word→0
   component). Genuinely distinct from the presence/absence set measure — created a dedicated concept
   rather than enriching [[kmer-jaccard-similarity]] (Euclidean captures k-mer abundance; Jaccard does
   not), wired as `alternative_to` it. Sources Zielezinski 2017 (word-vector model, Fig.1 x=ATGTGTG/
   y=CATGTG k=3) + Lau 2022 (frequency normalization + Euclidean metric) + Vinga-Almeida 2003 (4^k vector)
   + Boden 2014 (relative-frequency Euclidean). Oracles √0.11≈0.3316624790, AAAA/AAAT k=1 √0.125,
   identical→0, disjoint-single-kmer→√2. Two assumptions (ASM-01 case-fold, ASM-02 empty/L<k→zero-vector);
   count-based/Manhattan/Canberra/Chebyshev/cosine/D2/spaced-word not implemented; Deviations = None.
   Concept lists both the Evidence file and K-mer_Euclidean_Distance.md algorithm doc in sources. Wired
   into index (source + concept lines) + the algorithm-validation-evidence hub (frontmatter sources +
   body evidence-link + own-concept list). Reconciled backlog: moved K-mer/K-mer_Euclidean_Distance.md
   from pending K-mer (8→7) to covered-via-concept (51→52 covered / 194→193 pending, 30 domains).
   Contradictions: none.
   graph: +2 nodes (concept + source page), +2 typed edges (concept relates_to test-unit-registry +
   alternative_to kmer-jaccard-similarity; source/concept [[wikilinks]] mentions auto-derived)

## [2026-07-09] ingest | KMER-GENERATE-001-Evidence.md → kmer-generate-001-evidence (source) + k-mer-generation (concept)
   Fourth K-mer family Evidence. Test unit KMER-GENERATE-001, KmerAnalyzer.GenerateAllKmers — exhaustive
   enumeration of ALL possible k-mers of length k over an alphabet (the complete n^k universe Σ^k, 4^k for
   DNA), sequence-independent. Judged GENUINELY DISTINCT from the counting siblings → new dedicated concept
   [[k-mer-generation]] (generation = full n^k word set / frequency-array address space; counting = observed
   substrings + counts). Sources Wikipedia K-mer (n^k / 4^k, AGAT example) + BioInfoLogics 4^k (per-position
   Cartesian product) + Python itertools.product (k-fold product, odometer/lexicographic emission on sorted
   alphabet). Model = k-fold Cartesian product, lazy recursive prefix-extension, INV-01..04 (n^k count /
   all-distinct-set / length-k / sorted→lexicographic). Oracles k=1→{A,C,G,T}, k=2→16 AA..TT, k=3→64
   (AAA..TTT), protein 20^2=400, single-letter 1^4=1. Edge cases k≤0→ArgumentOutOfRangeException,
   empty alphabet→ArgumentException, unsorted→positional order, no dedup. One assumption (default "ACGT"
   sorted, documented property); Deviations = None. Cross-linked to sibling K-mer concepts; enriched
   both-strand-kmer-counting with an inbound [[k-mer-generation]] wikilink. Concept lists both the Evidence
   file and K-mer_Generation.md algorithm doc in sources. Reconciled backlog: moved K-mer/K-mer_Generation.md
   from pending K-mer (7→6) to covered-via-concept (52→53 covered / 193→192 pending, 30 domains).
   Contradictions: none.
   graph: +2 nodes (concept + source page), +2 typed edges (concept relates_to test-unit-registry +
   relates_to both-strand-kmer-counting; source/concept [[wikilinks]] mentions auto-derived)

## [2026-07-09] ingest | KMER-POSITIONS-001-Evidence.md → kmer-positions-001-evidence (source) + k-mer-positions (concept)
   Fifth K-mer family Evidence. Test unit KMER-POSITIONS-001, `KmerAnalyzer.FindKmerPositions(sequence, kmer)`
   — the ascending 0-based positions where a given k-mer occurs in a sequence (a *position / occurrence
   index*: *where*, not *how many*), solving the exact Pattern Matching Problem `Occ(P,T)={i∈[0,L−k]:
   T[i..i+k)=P}` with all overlapping starts reported. Judged GENUINELY DISTINCT from the counting units
   (positions = an ordered IEnumerable<int> of offsets for one k-mer, the inverse index to the
   Dictionary<string,int> count table) → new dedicated concept [[k-mer-positions]] rather than enriching a
   counting concept. Positioned as the single-pattern K-mer-family sibling of the multi-pattern exact
   matcher [[known-motif-search]] (`GenomicAnalyzer.FindMotif`) — same 0-based ascending all-overlapping
   semantics, one list vs a per-motif map. Sources: Rosalind BA1D (binding 0-based worked example ATAT/
   GATATATGCATATACTT→1 3 9, overlapping all reported) + Wikipedia k-mer (L−k+1 candidates, AGAT 2-mers) +
   Compeau & Pevzner Pattern Matching Problem (textbook 1-based prose deferred to BA1D's machine-checked
   0-based). INV-01..04 (match predicate / ascending / count=overlap-occurrence-count / range [0,L−k], empty
   when k>L). Oracles ATAT→[1,3,9], AA/AAAA→[0,1,2] self-overlap, whole-seq→[0], absent/longer/null-empty→
   empty (no throw). Impl = O(L·k) naive span scan (ReadOnlySpan.SequenceEqual, lazy yield), suffix-tree
   evaluated-and-rejected (unordered leaves + no single-query amortization). Three API-shape / repo-interop
   assumptions (0-based, case-insensitive upper-casing per sibling CountKmers, null/empty→empty); Deviations
   = None. Cross-linked reciprocally with [[known-motif-search]] (single-pattern counterpart note) and
   [[both-strand-kmer-counting]] (inverse-index note). Concept lists both the Evidence file and
   K-mer_Positions.md algorithm doc in sources. Wired into index (source + concept lines) + the
   algorithm-validation-evidence hub (frontmatter sources + body evidence-link + own-concept list).
   Reconciled backlog: moved K-mer/K-mer_Positions.md from pending K-mer (6→5) to covered-via-concept
   (53→54 covered / 192→191 pending, 30 domains). Contradictions: none.
   graph: +2 nodes (concept + source page), +2 typed edges (concept relates_to test-unit-registry +
   relates_to both-strand-kmer-counting; source/concept [[wikilinks]] mentions auto-derived)

## [2026-07-09] ingest | docs/Evidence/KMER-STATS-001-Evidence.md → kmer-stats-001-evidence (source) + 1 concept
   Sixth K-mer family Evidence file (after ASYNC, BOTH, DIST, GENERATE, POSITIONS). Judged
   `KmerAnalyzer.AnalyzeKmers` a GENUINELY DISTINCT companion summary layer over the shared CountKmers
   multiset — it reduces the count profile to a KmerStatistics bundle {TotalKmers, UniqueKmers, MaxCount,
   MinCount, AverageCount, Entropy} and adds the one formula unique to this unit, the Shannon **k-entropy**
   `E_k=−Σ p(α)log₂p(α)`, `p(α)=mult/(L−k+1)` (Manca 2021 arXiv:2106.15351 + Entropy–Rank Ratio
   arXiv:2511.05300), so it warrants its own concept rather than enriching a counting concept. Created
   concept [[k-mer-statistics]]. Captured the naming GOTCHA: `UniqueKmers` holds the **distinct** count
   (each different k-mer once), NOT the count==1 singletons (that is the separate KMER-UNIQUE-001 /
   Unique_And_MinCount_Kmers.md unit). Count facts (TotalKmers=L−k+1, distinct) from Wikipedia + BioInfoLogics
   count tables; AverageCount=total/distinct. Oracles GTAGAGCTGT k=1 (10/4/max4(G)/min1(C)/avg2.5/H1.846439…) +
   k=3 (8/8/1/1/H log₂8=3.0) + ATCGATCAC k=3 (7/6/2(ATC)/1/avg1.17/H2.521640…) + AAAA k=2 homopolymer
   (3/1/3/3/H0); corner cases homopolymer→H0/max=min=total, all-distinct→H log₂D/max=min=1, k>L·empty→all-zero,
   k≤0→ArgumentOutOfRangeException, case-insensitive. Two presentation-only assumptions (AverageCount rounded
   2dp via Math.Round; Entropy unrounded bits, tests within 1e-10), neither correctness-affecting. Concise
   source page for the artifact; concept lists both the Evidence file and K-mer_Statistics.md algorithm doc in
   sources. Wired into index (source + concept lines) + the algorithm-validation-evidence hub (frontmatter
   source + body evidence-link + own-concept list). Reconciled backlog: moved K-mer/K-mer_Statistics.md from
   pending K-mer (5→4) to covered-via-concept (54→55 covered / 191→190 pending, 30 domains). Contradictions:
   none — count tables and both k-entropy sources are mutually consistent.
   graph: +2 nodes (concept + source page), +2 typed edges (concept relates_to test-unit-registry +
   relates_to asynchronous-kmer-counting; source/concept [[wikilinks]] mentions auto-derived)

## [2026-07-09] ingest | docs/Evidence/KMER-UNIQUE-001-Evidence.md → kmer-unique-001-evidence (source) + unique-and-mincount-kmers (concept)
   Seventh K-mer family Evidence file (after ASYNC, BOTH, DIST, GENERATE, POSITIONS, STATS). Judged
   `KmerAnalyzer.FindUniqueKmers` + `FindKmersWithMinCount` a GENUINELY DISTINCT frequency-filtering
   unit — confirmed by the prior KMER-STATS-001 flag that "unique" (count==1 singletons) is THIS unit,
   separate from k-mer-statistics' `UniqueKmers`=distinct-count field. Created concept
   [[unique-and-mincount-kmers]] (expected backlog slug). Two operations filter the shared CountKmers
   multiset by per-k-mer Count at opposite ends of the distribution: FindUniqueKmers = Count==1
   singletons, FindKmersWithMinCount = Count≥minCount recurrent k-mers ordered by count desc. Captured
   the total/distinct/unique terminology (BioInfoLogics: unique="appear only once") and the reciprocal
   GOTCHA against [[k-mer-statistics]] (ATCGATCAC k=3 → 5 unique singletons vs 6 distinct; ATC=2
   excluded), surgically enriching the k-mer-statistics gotcha to link the new concept. Sources:
   Wikipedia K-mer (L−k+1 total, AGAT) + BioInfoLogics (distinct/unique, ATCGATCAC 7/6/5) + Compeau &
   Pevzner (`Count(Text,Pattern)`, most-frequent / Count≥t recurrent). Oracles ATCGATCAC k=3→{TCG,CGA,
   GAT,TCA,CAC}, AGAT k=2→{AG,GA,AT}, ACGTACGT k=4 (ACGT=2) FindKmersWithMinCount(…,2)→{(ACGT,2)} /
   (…,1)→all-4-count-desc / FindUniqueKmers→{CGTA,GTAC,TACG}, AAAAA k=3→∅; corner cases empty/k>L→empty,
   k≤0→ArgumentOutOfRangeException, case-insensitive. Two source-consistent assumptions (minCount≤1 ⇒
   Count≥minCount holds for all ⇒ returns all distinct count-desc; upper-casing per sibling methods),
   neither correctness-affecting. Concise source page; concept lists both the Evidence file and
   Unique_And_MinCount_Kmers.md algorithm doc in sources. Wired into index (source + concept lines) +
   the algorithm-validation-evidence hub (frontmatter source + body evidence-link + own-concept list).
   Reconciled backlog: moved K-mer/Unique_And_MinCount_Kmers.md from pending K-mer (4→3) to
   covered-via-concept (55→56 covered / 190→189 pending, 30 domains). Contradictions: none — Wikipedia,
   BioInfoLogics, and Compeau & Pevzner are mutually consistent; deviations None.
   graph: +2 nodes (concept + source page), +2 typed edges (concept relates_to test-unit-registry +
   relates_to k-mer-statistics; source/concept [[wikilinks]] mentions auto-derived)

## [2026-07-09] ingest | docs/Evidence/META-ALPHA-001-Evidence.md → meta-alpha-001-evidence (source) + 1 concept
   First per-algorithm Evidence file from the Metagenomics domain (new topic area — confirmed no
   existing metagenomics/diversity concept). Created the genuinely-distinct concept alpha-diversity —
   the anchor for the Metagenomics diversity family: within-sample diversity indices from one
   taxon→abundance map via `MetagenomicsAnalyzer.CalculateAlphaDiversity` → `AlphaDiversity` record with
   six fields. Formulas traced verbatim to the primary literature: observed richness S_obs=|{pᵢ>0}|,
   Shannon H=−Σpᵢln(pᵢ) using Math.Log (nats) per Shannon 1948, Simpson concentration λ=Σpᵢ² per Simpson
   1949, inverse Simpson 1/λ = Hill order-2 effective species per Hill 1973, Pielou evenness J=H/ln(S)
   for S>1 else 0 (standard ecological convention, ln(1)=0) per Pielou 1966, Chao1 S_obs+f₁²/(2f₂) with
   the f₂=0 bias-corrected branch S_obs+f₁(f₁−1)/2 per Chao 1984; Whittaker 1960 α/β/γ framing.
   Counts-or-proportions accepted (positive values internally normalized to sum 1), non-positive
   filtered (ln(0) undefined), O(n). INV-01..05 + empty/null→all-0 + single-species H0/λ1/J0 corner
   cases; oracles single→H0/λ1, (0.5,0.5)→ln2/0.5/2/J1, 4-equal→ln4/0.25/4/J1, (0.9,0.1)→H0.325/J0.469.
   FLAGGED NUANCE (not a contradiction): the Evidence file says "Deviations: None — all formulas match
   exactly", while the algorithm doc §5.4 records one accepted deviation — Chao1 falls back to
   ObservedSpecies for non-integer/proportional abundance input (data-type gate, not a formula change);
   captured on both the source and concept pages as consistent. Concise source page (Wikipedia
   Diversity-index/Alpha-diversity/Species-richness/Species-evenness + Shannon/Simpson/Hill/Chao/Pielou
   primaries). Wired into index (source + concept lines) + the algorithm-validation-evidence hub
   (frontmatter source + body evidence-link + own-concept list). Reconciled backlog: moved
   Metagenomics/Alpha_Diversity.md from pending Metagenomics (10→9) to covered-via-concept (56→57
   covered / 189→188 pending, 30 domains). Cross-linked [[beta-diversity]] as a not-yet-created future
   sibling (no stub). Contradictions: none.
   graph: +2 nodes (concept + source page), +1 typed edge (concept relates_to test-unit-registry;
   source/concept [[wikilinks]] mentions auto-derived)

## [2026-07-09] ingest | META-BETA-001-Evidence.md → meta-beta-001-evidence (source) + beta-diversity (concept)
   Second Metagenomics-family unit: between-sample dissimilarity CalculateBetaDiversity → Bray-Curtis
   (abundance, not a true metric) + Jaccard distance (presence/absence, true metric), from
   Whittaker 1960 (α/β/γ) + Bray & Curtis 1957 + Jaccard 1901 + Wikipedia primaries. Created a dedicated
   [[beta-diversity]] concept (the [[alpha-diversity]] page already referenced it as its expected
   sibling) + [[meta-beta-001-evidence]] source page. Wired into index (source + concept lines) + the
   algorithm-validation-evidence hub (frontmatter source + body evidence-link + own-concept list).
   Cross-linked reciprocally with [[alpha-diversity]] (within- vs between-sample halves of Whittaker's
   framework; added a relates_to edge on each). Noted ecological-Jaccard shares the index math but not
   the domain with sequence [[kmer-jaccard-similarity]] (prose mention only, no typed edge). Reconciled
   backlog: moved Metagenomics/Beta_Diversity.md from pending Metagenomics (9→8) to covered-via-concept
   (57→58 covered / 188→187 pending, 30 domains). Contradictions: none.
   graph: +2 nodes (concept + source page), +3 typed edges (beta relates_to test-unit-registry;
   beta relates_to alpha-diversity; alpha relates_to beta-diversity — reciprocal); body [[wikilinks]]
   mentions auto-derived

## [2026-07-09] ingest | docs/Evidence/META-BIN-001-Evidence.md → meta-bin-001-evidence (source) + metagenomic-binning (concept)
   Third Metagenomics-family unit. Created source page [[meta-bin-001-evidence]] and new concept
   [[metagenomic-binning]] (MetagenomicsAnalyzer.BinContigs — k-means over composite distance
   |ΔGC|+|Δcoverage|+TNF-Pearson-distance → MAGs; completeness/contamination are length-ratio/GC-variance
   PROXIES, not CheckM marker calls; opt-in TETRA z-score signature CalculateTetranucleotideZScores/
   TetranucleotideZScoreCorrelation, z(ACGT)=√5 oracle). Flagged the CheckM marker-gene QC as an
   explicit honest residual (not implemented) and the three now-resolved prior assumptions (deviations
   None). Cross-linked as Metagenomics-family sibling of [[alpha-diversity]]/[[beta-diversity]] in prose
   (different question — genome reconstruction vs community diversity — so no typed sibling edge, source
   does not assert one). Reconciled backlog: moved Metagenomics/Genome_Binning.md from pending
   Metagenomics (8→7) to covered-via-concept (58→59 covered / 187→186 pending, 30 domains).
   Contradictions: none.
   graph: +2 nodes (concept + source page), +1 typed edge (metagenomic-binning relates_to
   test-unit-registry); body [[wikilinks]] mentions auto-derived

## [2026-07-09] ingest | docs/Evidence/META-BIN-001-MarkerQC-Evidence.md → meta-bin-001-markerqc-evidence (source); enriched metagenomic-binning (concept)
   ADDENDUM to META-BIN-001 — validates the CheckM-style single-copy marker-gene completeness/
   contamination now built on top of the TNF/coverage binning. Created source page
   [[meta-bin-001-markerqc-evidence]]. ENRICHED the existing [[metagenomic-binning]] concept rather
   than creating a new one (the marker QC is the quality-metric layer of binning, not a separate
   wiki-worthy algorithm): rewrote the proxy-vs-CheckM GOTCHA (the residual is now BUILT but exposed
   through a distinct opt-in API `EstimateBinQualityFromMarkerCounts`/`EstimateBinQualityFromMarkers`/
   `DetectMarkers`, NOT wired into `BinContigs`, whose fields stay length-ratio/GC-variance proxies);
   added a Marker-gene QC section (CheckM Eqs. 1–2 over collocated sets `M`, multi-copy counts once
   toward present + N−1 toward contamination; bundled CC0 Pfam sets = 9 ribosomal + bac120 6 + ar122
   35 as singleton sets, TIGRFAM CC BY-SA NOT bundled/caller-supplied; glocal Plan7 Viterbi ≥ Pfam GA1
   gate vs HMMER local+null2 engine diff; oracles 250/3%≈83.333 comp / 100/9%≈11.111 cont, uS8→PF00410
   +176 bits); refreshed the scope/limitations paragraph. Added the new source + HEAD source_commit to
   the concept frontmatter. Added a forward-pointer on the base [[meta-bin-001-evidence]] source page
   (its "honest residual" note now flags the addendum built it). Wired into index (new source line +
   refreshed the metagenomic-binning concept summary). Hub [[algorithm-validation-evidence]] frontmatter
   NOT edited — its per-file list drifted (base META-BIN/META-BETA absent too); the source page links
   the hub in prose (mention edge). Backlog: no change — base Genome_Binning.md already covered-via-
   concept, the addendum has no separate docs/algorithms doc. Contradictions: none (the addendum
   supersedes the base file's "not implemented" residual; recorded as an evolution, not a conflict).
   graph: +1 node (source page); no new typed edges (concept already relates_to test-unit-registry);
   body [[wikilinks]] mentions auto-derived

## [2026-07-09] ingest | docs/Evidence/META-CLASS-001-Evidence.md → meta-class-001-evidence (source) + taxonomic-classification (concept)
   Fourth Metagenomics-family Evidence unit. Source page [[meta-class-001-evidence]] + new concept
   [[taxonomic-classification]] (faithful Kraken k-mer/LCA/RTL per-read classifier: canonical-k-mer→
   LCA-of-owning-taxa database, classification-tree max-scoring root-to-leaf path, tie→LCA of leaves,
   Confidence=C/Q, no-hit→Unclassified root). Genuinely distinct concept — per-read/LCA assignment, not
   diversity or binning; deliberately scoped to classification (abundance profiling Taxonomic_Profile
   left as a separate future unit). Cross-linked to siblings [[metagenomic-binning]] +
   [[alpha-diversity]]/[[beta-diversity]]. index.md: +1 source +1 concept. Backlog: moved
   Metagenomics/Taxonomic_Classification.md pending→covered (59→60 covered / 186→185 pending; §Metagenomics 7→6).
   Contradictions: none (pre-C1 flat best-hit wording superseded by the LCA/RTL enhancement — recorded
   as evolution, not conflict; Evidence file lists no open questions / no deviations).
   graph: +2 nodes (source + concept), +3 typed edges (concept relates_to test-unit-registry /
   metagenomic-binning / alpha-diversity); body [[wikilinks]] mentions auto-derived.

## [2026-07-09] ingest | docs/Evidence/META-FUNC-001-Evidence.md → meta-func-001-evidence (source) + functional-prediction (concept)
   Fifth Metagenomics-family Evidence unit. Source page [[meta-func-001-evidence]] + new concept
   [[functional-prediction]] (PICRUSt/KO-style functional prediction in two exact-numeric pieces:
   (A) homology-based annotation transfer `PredictFunctions` — exact-signature `string.Contains` hit
   scored by BLOSUM62 self-score, BLAST bit `S'=(λS−lnK)/ln2` + E-value `E=K·m·n·e^(−λS)=m·n·2^(−S')`
   with ungapped BLOSUM62 λ=0.3176/K=0.134 (Altschul tutorial + NCBI blast_stat.c + BLOSUM62 diagonals),
   best hit = lowest E-value; (B) hypergeometric pathway ORA `FindPathwayEnrichment` — right-tail
   P(X≥x) in log-Gamma space, x/M/n=0→p=1, sorted ascending). Genuinely distinct concept (functional
   capability, not who-is-there / diversity). Cross-linked to all four siblings [[taxonomic-classification]]
   (added a reciprocal "who is there vs what can they do" nav link there) / [[metagenomic-binning]] /
   [[alpha-diversity]] / [[beta-diversity]], and to [[alignment-statistics]] (BLAST significance is a
   different layer from percent-id). One assumption ASM-01 = ungapped exact-match model (affects which
   hits found, not the bit-score/E-value formulas); Evidence lists no contradictions. Oracles WWW→
   S'18.0202932787533/E 3.3852730346546e−5 (both forms agree) + ORA N8000/M400/n100/x20→7.88e−8.
   Hub [[algorithm-validation-evidence]]: added META-FUNC to frontmatter sources (bumped source_commit to
   HEAD) + source-list + concept-list. index.md: +1 source +1 concept. Backlog: moved
   Metagenomics/Functional_Prediction.md pending→covered (60→61 covered / 185→184 pending; §Metagenomics
   6→5). SCOPE NOTE: the shared ORA half (`FindPathwayEnrichment`/`HypergeometricUpperTail`) is its OWN
   unit META-PATHWAY-001 (Pathway_Enrichment_ORA.md, separate META-PATHWAY-001-Evidence.md, not yet
   ingested) — META-FUNC-001 validates only Functional_Prediction.md; flagged on both the source and
   concept pages so META-PATHWAY-001 can share this material later.
   graph: +2 nodes (source + concept), +1 typed edge (concept relates_to test-unit-registry); body
   [[wikilinks]] mentions auto-derived.

## [2026-07-09] ingest | docs/Evidence/META-PATHWAY-001-Evidence.md → meta-pathway-001-evidence (source) + pathway-enrichment-ora (concept)
   Sixth Metagenomics-family Evidence unit and the DEDICATED unit for the ORA / hypergeometric machinery
   that META-FUNC-001 exercised as component B. DECISION: created a focused new concept
   [[pathway-enrichment-ora]] that OWNS the method (rather than only enriching [[functional-prediction]]) —
   the evidence is substantial enough to stand alone: its own GO::TermFinder (Boyle 2004) + PNNL ORA §8.2
   sources, the M↔n symmetry invariant, and exact hand-derived rational oracles. The concept synthesizes
   the hypergeometric right-tail `P(X≥x)=1−Σ_{i=0}^{x−1}C(M,i)C(N−M,n−i)/C(N,n)` (`phyper(x−1,M,N−M,n,
   lower.tail=FALSE)`, N=background/M=pathway/n=query/x=overlap, upper-tail/without-replacement), log-Gamma
   summation to N=8000, p∈[0,1], sorted ascending; p=1 when x/M/n=0; background = explicit else
   union-of-pathway-members default (query unioned in, members intersected); NO BH/Bonferroni FDR. Oracles
   PNNL N8000/M100/n400/x20→7.88e−8 + exact 1/252 / 5/6 / 1 / 251/252. Created source page
   [[meta-pathway-001-evidence]]. Reciprocally cross-linked with [[functional-prediction]]: rewrote its
   component-B blockquote to defer ownership here (was "not yet ingested / may get its own page") and added
   a typed edge functional-prediction relates_to pathway-enrichment-ora; also updated the
   [[meta-func-001-evidence]] source page's scope note ("now ingested"). Hub
   [[algorithm-validation-evidence]]: added META-PATHWAY to frontmatter sources (bumped source_commit to
   HEAD 14005a6) + source-list + concept-list. index.md: +1 source +1 concept, refreshed the
   functional-prediction / meta-func lines. Backlog: moved Metagenomics/Pathway_Enrichment_ORA.md
   pending→covered (61→62 covered / 184→183 pending; §Metagenomics 5→4). Contradictions: none — Boyle 2004
   and PNNL §8.2 give the identical right-tail formula; the background-defaulting assumption is
   formula-preserving and caller-overridable. Note: the ORA statistic is generic (GO/proteomics sources)
   though registered under metagenomics via `FindPathwayEnrichment`.
   graph: +2 nodes (source + concept), +3 typed edges (pathway-enrichment-ora relates_to
   test-unit-registry + relates_to functional-prediction; reciprocal functional-prediction relates_to
   pathway-enrichment-ora); body [[wikilinks]] mentions auto-derived.

## [2026-07-09] ingest | docs/Evidence/META-PROF-001-Evidence.md → meta-prof-001-evidence (source) + 1 concept
   Seventh ingested Metagenomics-family Evidence file (META-PROF-001). Decision: created the
   genuinely-distinct concept [[taxonomic-profile]] rather than enriching [[taxonomic-classification]] —
   profiling is the aggregation/estimation step the classification unit explicitly deferred, with its own
   method `MetagenomicsAnalyzer.GenerateTaxonomicProfile(IEnumerable<TaxonomicClassification>)` producing a
   `TaxonomicProfile` (relative-abundance maps at four ranks kingdom/phylum/genus/species = count(taxon)/
   Σcount(classified), inline species-level Shannon H=−Σpᵢln(pᵢ) nats + Simpson concentration λ=Σpᵢ², and
   TotalReads/ClassifiedReads). Counting rules: Unclassified excluded from denominators, empty rank strings
   filtered, per-rank Σ≈1.0. Invariants ClassifiedReads≤TotalReads / =Σ(counts any rank) / Shannon≥0 /
   0≤Simpson≤1; oracles Shannon=ln(3) (3 uniform), Simpson=0.375 ([2,1,1]), TotalReads3/ClassifiedReads2;
   empty→0/0/empty & single taxon→1.0/H0/λ1 vs empty→λ0 (empty-sum convention). Sources Wikipedia
   Metagenomics + Relative-abundance-distribution + MetaPhlAn docs + Segata 2012 (Nature Methods). Created
   source page [[meta-prof-001-evidence]]. Cross-linked reciprocally: rewrote taxonomic-classification's
   deferred-profiling sentence to point at [[taxonomic-profile]] ("not yet ingested"→link + input-shape).
   Hub [[algorithm-validation-evidence]]: added META-PROF to frontmatter sources (bumped source_commit to
   HEAD 02f28f4) + source-list + concept-list. index.md: +1 source +1 concept. Backlog: moved
   Metagenomics/Taxonomic_Profile.md pending→covered (62→63 covered / 183→182 pending; §Metagenomics 4→3).
   Contradictions: none — the verified design decisions (nats log, concentration-index λ, empty→0) are
   mathematical facts, no literature deviations. Scope note: count-based tally, NOT MetaPhlAn marker-gene
   coverage estimation; no genome-size/copy-number correction; inherits upstream classifier accuracy.

## [2026-07-09] ingest | docs/Evidence/META-TAXA-001-Evidence.md → meta-taxa-001-evidence (source) + significant-taxa-detection (concept)
   Eighth ingested Metagenomics-family Evidence file (META-TAXA-001). Decision: created the
   genuinely-distinct concept [[significant-taxa-detection]] rather than folding into an existing unit —
   community **differential abundance** via the per-taxon two-group **Mann–Whitney U / Wilcoxon rank-sum**
   test is a distinct *statistical test* from the hypergeometric [[pathway-enrichment-ora]] and the
   Fisher's-exact [[differentially-methylated-regions]] (cross-linked as alternatives-by-test). Two methods
   `MetagenomicsAnalyzer.MannWhitneyU(group1,group2,useContinuityCorrection=true)` (core, U1/U2/z/p) +
   `FindSignificantTaxa(profiles,groups,pThreshold=0.05,useContinuityCorrection=true)` (per-taxon →
   SignificantTaxon ascending by p). Model: pool→midranks (Σ(t³−t)) → U1=R1−n1(n1+1)/2, U2=n1·n2−U1,
   m_U=n1·n2/2, tie-corrected σ_U, z=(|U−m_U|−cc)/σ_U on max(U1,U2), two-tailed p=2·(1−Φ(z)) via shared
   `StatisticsHelper.NormalCDF` (A&S 7.1.26 erf, ≈1e−6). INV-01..06 incl. all-tied→σ0→p1 and group-swap
   symmetry; oracles SciPy x[19,22,16,29,24]/y[20,11,17,12]→U1=17/U2=3/σ=sqrt(200/12)/z_cc=1.5922→p≈0.11135
   & z_nocc=1.7146→p≈0.08641, tortoise/hare U_T=11/U_H=25/sum=36. Sources Wikipedia Mann–Whitney U (Mann &
   Whitney 1947) + SciPy mannwhitneyu + Xia & Sun 2017 (PMC6128532, microbiome domain) + A&S 7.1.26.
   Created source page [[meta-taxa-001-evidence]]. Cross-linked: [[significant-taxa-detection]] depends_on
   [[taxonomic-profile]] (consumes its per-sample abundance vectors) + reciprocal mention added to
   taxonomic-profile's scope paragraph. Hub [[algorithm-validation-evidence]]: added META-TAXA to
   frontmatter sources (bumped source_commit→HEAD b8447d68) + source-list + concept-list. index.md: +1
   source +1 concept. Backlog: moved Metagenomics/Significant_Taxa_Detection.md pending→covered (64→65
   covered / 181→180 pending; §Metagenomics 2→1). Contradictions: none — three source-backed assumptions
   (continuity-correction-on default = SciPy, two-tailed, two-label/absence=0); only simplifications are
   asymptotic-not-exact p and A&S-7.1.26 Φ numerics. Scope: two-group only, no FDR (caller applies BH),
   rank test ignores compositionality.
   graph: +2 nodes, +3 typed edges
   graph: +2 nodes (source + concept), +3 typed edges (taxonomic-profile relates_to test-unit-registry +
   depends_on taxonomic-classification + relates_to alpha-diversity); body [[wikilinks]] mentions auto-derived.
- 2026-07-09 — ingest `docs/Evidence/META-RESIST-001-Evidence.md` (test unit META-RESIST-001,
  Antibiotic-Resistance Gene Detection; seventh Metagenomics-family unit). Created source
  [[meta-resist-001-evidence]] + NEW concept [[antibiotic-resistance-gene-detection]] (genuinely
  distinct method — no prior AMR/resistance concept). Algorithm =
  `MetagenomicsAnalyzer.FindAntibioticResistanceGenes(contigs, referenceGenes, id=0.90, cov=0.60)`:
  ResFinder-style screen of assembled contigs vs a CALLER-SUPPLIED resistance-gene reference DB
  (curated CARD/ResFinder tables not embedded). Private `BestUngappedMatch` slides each reference
  across the contig at every offset −(m−1)..n−1 (overhanging both ends so contig-edge/split genes
  score against the reference length), keeps the max-match window (tie→shorter=higher identity),
  then identity=matches/w (BLAST gapless denominator, Heng Li 2018) & coverage=w/m (fraction of
  REFERENCE length); reports the reference iff identity≥idThreshold AND coverage≥covThreshold;
  best-matching gene per contig = max identity, tie→max coverage (Zankari 2012 "best-matching
  gene"; CARD RGI best-hit by bit score). INV-01..05; defaults 0.90 ID / 0.60 cov named constants;
  oracles CGTACGT@AAACGTACGT→1.0/1.0, CGTTCGT vs CGTACGT→6/7≈0.857/1.0, contig-edge CGTA→1.0 /
  4⁄7≈0.571. Sources: Zankari 2012 (original ResFinder) + ResFinder GitHub (-t 0.80/-l 0.60) + Sci
  Rep 2023 + JAC 2016 (98% ID/60% cov, edge/split rationale) + Heng Li 2018 (identity formula) +
  CARD RGI. One assumption ASM-01 = gapless ungapped model (indel-requiring matches under-scored vs
  gapped BLAST; substitution divergence + contig-edge truncation scored exactly). Cross-linked
  [[functional-prediction]] as the sibling BLAST-style homology screen (shared machinery; AMR scores
  nucleotide identity/coverage, PredictFunctions a BLOSUM62 protein bit-score/E-value) — comparison
  table on the concept. Hub [[algorithm-validation-evidence]]: added META-RESIST to frontmatter
  sources (bumped source_commit to HEAD c81ef58a) + source-list + concept-list. index.md: +1 source
  +1 concept. Backlog: moved Metagenomics/Antibiotic_Resistance_Detection.md pending→covered (63→64
  covered / 182→181 pending; §Metagenomics 3→2). Contradiction flagged (non-blocking): the evidence
  file's extracted ResFinder README default is 0.80 ID (and the study SELECTED 0.98), while the
  implementation ships 0.90 ID as the default — recorded as a threshold-provenance note on the source
  page; the 0.90 constant is user-selectable so it does not change the algorithm, only the operating
  point. graph: +2 nodes (source + concept), +2 typed edges (antibiotic-resistance-gene-detection
  relates_to test-unit-registry + relates_to functional-prediction); body [[wikilinks]] mentions
  auto-derived.
- 2026-07-09 — ingest `docs/Evidence/MIRNA-PAIR-001-Evidence.md` (test unit MIRNA-PAIR-001,
  MiRNA-Target Pairing Analysis; FIRST MiRNA-family unit — NEW topic area, no prior RNA
  base-pairing / miRNA concept existed). Created source [[mirna-pair-001-evidence]] + NEW concept
  [[rna-base-pairing]] ("RNA base pairing (Watson-Crick + G-U wobble) and the miRNA-target
  duplex"). Algorithm = `MiRnaAnalyzer.AlignMiRnaToTarget` + `CanPair`/`IsWobblePair`/
  `GetReverseComplement`: `CanPair`⟺{A-U,U-A,G-C,C-G,G-U,U-G} = Watson-Crick {A-U,G-C}
  (Agarwal 2015 / PMC4532895) + the single standard **G-U wobble** (Crick 1966), `IsWobblePair`⟺
  {G-U,U-G} (wobble⊆pairable, counted separately from matches per PMC4870184); `GetReverseComplement`
  = antiparallel RNA reverse complement for seed→target motif (Lewis 2005; let-7a `GAGGUAG`→
  `CUACCUC`); `AlignMiRnaToTarget` pairs miRNA[i]↔target[len−1−i] over the shorter overlap, ungapped,
  `|`(WC)/`:`(wobble)/space(mismatch), counts sum to min(len)/Gaps=0, ΔG = simplified Turner-2004
  stacking sum over consecutive paired runs (sign reliable — fully-WC ≤0, all-mismatch ≥0 —
  magnitude NOT). Oracles AAAA/UUUU→4 matches, GGGG/UUUU→4 wobbles, AAAA/AAAA→4 mismatches. Made the
  concept the **shared base-pairing primitive** anchor (per ingest brief: Watson-Crick/G-U wobble is
  a primitive both RNA-structure and miRNA use) — documented so a future RnaStructure
  `RNA_Base_Pairing.md` ingest can reference/enrich the same page rather than duplicate the rule.
  Hub [[algorithm-validation-evidence]]: added MIRNA-PAIR to frontmatter sources (bumped
  source_commit to HEAD da06ef55) + source-list + concept-list. index.md: +1 source +1 concept.
  Backlog: moved MiRNA/MiRNA_Target_Pairing.md pending→covered (65→66 covered / 180→179 pending;
  §MiRNA 4→3). One ASSUMPTION recorded (Turner stacking numerics not re-retrieved this session →
  tests assert base-pairing structure + ΔG sign, not kcal/mol magnitude); A-opposite-position-1 is
  Argonaute recognition not base pairing (out of scope). No contradictions. graph: +2 nodes
  (source + concept), +1 typed edge (rna-base-pairing relates_to test-unit-registry); body
  [[wikilinks]] mentions auto-derived.
- 2026-07-09 — ingest `docs/Evidence/MIRNA-PRECURSOR-001-Evidence.md` (test unit MIRNA-PRECURSOR-001,
  Pre-miRNA Hairpin Detection; SECOND MiRNA-family unit). Created source
  [[mirna-precursor-001-evidence]] + NEW concept [[pre-mirna-hairpin-detection]] (genuinely distinct
  method — precursor stem-loop hairpin detection, not the miRNA-target duplex of [[rna-base-pairing]]).
  Algorithm = `MiRnaAnalyzer`: DEFAULT heuristic `FindPreMiRnaHairpins` counts uninterrupted
  complementary pairs ({A-U,G-C}+G-U wobble — the [[rna-base-pairing]] primitive) from both ends
  inward → accept iff stem ≥18 bp (Krol 2004) + loop 3-25 nt (Bartel 2004); extracts mature(5' arm)/
  star(3' arm), balanced dot-bracket, Turner-2004 ΔG (stacking+loop+terminal-mismatch+0.45 AU/GU).
  DOCUMENTED LIMITATION (accepted, not a bug): consecutive-pairing is stricter than real structure →
  rejects natural miRBase precursors (hsa-mir-21 16 end-pairs, let-7a-1 5, tests M18/M19). Three
  OPT-IN production paths (default unchanged): (1) `AssessHairpinByMfe`/`FindPreMiRnaHairpinsByMfe`
  fold via the RNA-STRUCT-001 Zuker–Stiegler engine and read the hairpin from the real MFE structure
  (single dominant hairpin/no multibranch + stem bp ≥16 (Ambros 2003) + loop 3-25 + MFEI ≥0.85 (Zhang
  2006, AMFE=100·|ΔG°|/n, MFEI=AMFE/GC%)) → detects hsa-mir-21 (ΔG° −35.13/32 bp/MFEI 1.0037) &
  let-7a-1 (ΔG° −34.31/MFEI 1.0091) the heuristic rejects; a 120-nt multibranch 5S-rRNA-like fold is
  REJECTED on STRUCTURE (multibranch, not a single dominant hairpin) despite a strongly negative
  ΔG° −47.04 — proving acceptance rests on topology, not merely a weak ΔG°. (2)
  `PredictDroshaDicerCleavage` = published measuring ruler only — Drosha
  +11 bp from basal junction (Han 2006), Dicer 22-nt 5'-counting (Park 2011), RNase III 2-nt 3'
  overhang (Lee 2003), optional CNNC 16-18 nt confidence flag (Auyeung 2013); hsa-miR-21-5p
  cross-check reproduces `UAGCUUAUCAGACUGAUGUUGA` (22 nt) exactly. (3) `ClassifyPreMiRna` = trained
  logistic regression over [ΔG,AMFE,MFEI,GC,%paired], 13 public-domain miRBase positives vs
  Altschul-Erickson 1985 di-shuffle negatives (Bonnet 2004 convention), held-out accuracy=AUC=1.0 —
  NO GPL miRDeep2 code/weights. Sources: Bartel 2004/2009 + Ambros 2003 + Krol 2004 + miRBase +
  Wikipedia + Bonnet 2004 + Zhang 2006 + Meyers 2008 + Han 2006 + Park 2011 + Lee 2003 + Auyeung 2013
  + Altschul-Erickson 1985 + Turner 2004. Two accepted assumptions (ASM-03 5'-arm mature extraction;
  ASM-01 uninterrupted-stem strictness — both mitigated by the opt-in MFE fold); residual read-stacking
  miRDeep2 signal data-blocked (needs caller's reads). No contradictions. Hub
  [[algorithm-validation-evidence]]: added MIRNA-PRECURSOR to frontmatter sources (bumped source_commit
  to HEAD e0541d58) + source-list + concept-list. rna-base-pairing: added reciprocal sibling nav link.
  index.md: +1 source +1 concept. Backlog: moved MiRNA/Pre_miRNA_Detection.md pending→covered (66→67
  covered / 179→178 pending; §MiRNA 3→2). graph: +2 nodes (source + concept), +2 typed edges
  (pre-mirna-hairpin-detection relates_to test-unit-registry + depends_on rna-base-pairing); body
  [[wikilinks]] mentions auto-derived.
- 2026-07-09 — ingest `docs/Evidence/MIRNA-SEED-001-Evidence.md` (test unit MIRNA-SEED-001, Seed
  Sequence Analysis; THIRD MiRNA-family unit). Created source [[mirna-seed-001-evidence]] + NEW concept
  [[seed-sequence-analysis]] (genuinely distinct — string-level seed extraction / family equality, not
  the base-pairing predicate/duplex of [[rna-base-pairing]] nor the precursor hairpins of
  [[pre-mirna-hairpin-detection]]). Algorithm = `MiRnaAnalyzer`: `GetSeedSequence` returns positions
  **2-8** (7-nt extended seed) via `Substring(1,7)` uppercase — casing only, **no** T→U (that is
  `CreateMiRna`), `<8 nt`/null/empty → `""`; `CreateMiRna(name, sequence)` normalises
  `ToUpperInvariant()`+`T→U`, extracts the seed from the normalised sequence, stores `SeedSequence` +
  fixed zero-based `SeedStart=1`/`SeedEnd=7`; `CompareSeedRegions` = Hamming over the 7-nt seed
  (`Matches`+`Mismatches`=7, mismatches also count length diff), `IsSameFamily` ⟺ exact seed equality,
  empty seed → zeroed. **miRNA family = identical 2-8 seed** (let-7a/-7b/-7c-5p all `GAGGUAG` → same
  family; miR-21-5p `AGCUUAU` differs; self→0 mismatches). Sources: Wikipedia MicroRNA + TargetScan
  FAQ/7mer + Lewis 2005 + Bartel 2009 + Agarwal 2015 + Grimson/Friedman + miRBase. Domain context: site
  ladder 8mer/7mer-m8/7mer-A1/6mer over the 2-7 (6-nt canonical) vs 2-8 (7-nt extended) distinction, but
  matching-to-target + site-class assignment DEFERRED to target-site prediction (MIRNA-TARGET-001,
  future); seed→target reverse complement owned by [[rna-base-pairing]] (`GetReverseComplement`).
  Intentionally simplified: exact-7-mer family equality (no isomiR/offset/noncanonical seeds, not a
  curated taxonomy). Terminology nuance FLAGGED (2-7-vs-2-8 collapse), no source contradictions. Hub
  [[algorithm-validation-evidence]]: added MIRNA-SEED to frontmatter sources (bumped source_commit to
  HEAD 989c8a14) + source-list + concept-list. rna-base-pairing: added reciprocal seed-extraction nav
  link (§2 seed→target). index.md: +1 source +1 concept. Backlog: moved MiRNA/Seed_Sequence_Analysis.md
  pending→covered (67→68 covered / 178→177 pending; §MiRNA 2→1). graph: +2 nodes (source + concept),
  +2 typed edges (seed-sequence-analysis relates_to test-unit-registry + relates_to rna-base-pairing);
  body [[wikilinks]] mentions auto-derived.

## [2026-07-09] ingest | docs/Evidence/MIRNA-TARGET-001-Evidence.md → mirna-target-001-evidence (source) + 1 concept
   miRNA target-site prediction — the FOURTH and FINAL MiRNA-family unit (COMPLETES the family).
   Created concept [[mirna-target-site-prediction]]: two-pass antiparallel seed-RC scan classifying the
   Bartel/TargetScan hierarchy (8mer=2-8+A1 / 7mer-m8=2-8 / 7mer-A1=2-7+A1 / 6mer=2-7 / offset-6mer=3-8,
   higher classes suppress overlapping offset-6mer), heuristic score (base 1.0/0.52/0.32/0.15/0.10,
   +0.05 >10 matches, −0.01/mismatch, clamp [0,1]) + heuristic ΔG; opt-in TargetScan context++ scorer
   (per-site-type MLR, min-max-scaled continuous + raw indicators; computed Local_AU/3P_score/Min_dist/
   Len_3UTR/Off6m + ComputeTa3Utr TA=log10 N + McCaskill-partition SA + Friedman-Bls PCT; SPS/Len_ORF/
   ORF8m/PCT-sigmoid caller-supplied → partial CS + OmittedFeatures). Sources: Bartel 2009 + Lewis 2005 +
   Grimson 2007 + Agarwal 2015 + Garcia 2011 + Friedman 2009 + McCaskill/ViennaRNA + TargetScan 8 +
   miRBase. let-7a GAGGUAG→CUACCUC site oracles; 8mer partial CS −0.7561913315126536; TA=log10(5)=0.69897.
   No source contradictions (heuristic-score + partial-CS + unemitted Centered/Supplementary enum are
   intentional simplifications). Hub [[algorithm-validation-evidence]]: added MIRNA-TARGET to frontmatter
   sources (bumped source_commit to HEAD aa11631f) + source-list + concept-list. Reciprocal nav links added
   on [[seed-sequence-analysis]] (target now depends_on it) and [[rna-base-pairing]] (finder depends on
   GetReverseComplement + AlignMiRnaToTarget). index.md: +1 source +1 concept. Backlog: moved
   MiRNA/Target_Site_Prediction.md pending→covered (68→69 covered / 177→176 pending; §MiRNA now 0, 30→29
   domains). graph: +2 nodes (source + concept), +3 typed edges (mirna-target-site-prediction relates_to
   test-unit-registry + depends_on seed-sequence-analysis + depends_on rna-base-pairing); body [[wikilinks]]
   mentions auto-derived.

## [2026-07-09] ingest | docs/Evidence/MOTIF-CONS-001-Evidence.md → motif-cons-001-evidence (source) + 1 concept
   Consensus from a multiple alignment (MotifFinder.CreateConsensusFromAlignment) — a Motif-Analysis unit
   distinct from the assembly [[consensus-sequence]] (ASSEMBLY-CONSENSUS-001). Created concept
   [[consensus-from-alignment]]: PURE most-frequent (plurality) column consensus over equal-length aligned
   strings, deterministic ALPHABETICAL tie-break (A<C<G<T), NO threshold (always emits — no n/x
   no-consensus output). Sources: Wikipedia "Consensus sequence" (Schneider & Stephens 1990) + Rosalind
   CONS (profile matrix + equal-length precondition + ties→multiple valid) + EMBOSS cons (the plurality-
   threshold alternative NOT adopted) + Geneious/LANL (alphabetical tie-break). Oracles: Rosalind 7×8
   sample → profile A=`5 1 0 0 5 5 0 0`/C/G/T → consensus ATGCAACT; tie-break AT+GT→AT; identical→unchanged;
   single→unchanged. Two documented assumptions (alphabetical tie-break, no-threshold scope — the area's
   IUPAC-degenerate GenerateConsensus + PWM CreatePwm are separate methods, not stubbed). Contract:
   equal-length→ArgumentException, non-ACGT→ArgumentException, null→ArgumentNullException, empty→"".
   No source contradictions. Reciprocal nav cross-link added on [[consensus-sequence]] (kept its own
   frontmatter sources per precedent). Hub [[algorithm-validation-evidence]]: added MOTIF-CONS to
   frontmatter sources (bumped source_commit to HEAD de59ece4) + source-list + concept-list. index.md:
   +1 source +1 concept. Backlog: moved Pattern_Matching/Consensus_From_Alignment.md pending→covered
   (69→70 covered / 176→175 pending; §Pattern_Matching 9→8, domains still 29).
   graph: +2 nodes (source + concept), +2 typed edges (consensus-from-alignment relates_to
   test-unit-registry + alternative_to consensus-sequence); body [[wikilinks]] mentions auto-derived.

## [2026-07-09] ingest | docs/Evidence/MOTIF-DISCOVER-001-Evidence.md → motif-discover-001-evidence (source) + 1 concept
   Per-algorithm Evidence file; second Motif-family unit (after MOTIF-CONS / GENOMIC-MOTIFS).
   Created the genuinely-distinct concept overrepresented-kmer-discovery — the de novo
   motif-discovery method (`MotifFinder.DiscoverMotifs`): enumerate every length-k k-mer of ONE
   DNA sequence, count overlapping occurrences, rank by observed/expected enrichment
   `Count / ((N−k+1)/4^k)` under a zero-order i.i.d. uniform background (Compeau & Pevzner);
   deterministic exact single-pass hash-map, 0-based positions. Distinct from the sibling motif
   concepts by *question asked*: it finds UNKNOWN over-represented words (motif = output) whereas
   [[known-motif-search]] matches a supplied set of KNOWN motifs (motif = input) — modelled
   `alternative_to` it — and [[consensus-from-alignment]] collapses an already-aligned instance
   set. Oracles: `ATGC` in `ATGCATGCATGC` k=4 → Count 3 @ {0,4,8}, E=9/256, enrichment 768/9≈85.333;
   `AAA` in `AAAAAAAAAA` k=3 → Count 8, E=0.125, enrichment 64.0. Corner cases k>N→empty, null→
   ArgumentNullException, k<1→ArgumentOutOfRangeException. One assumption: minCount (default 2) is a
   presentation threshold, not correctness-affecting (O/E defined for every k-mer). Intentional
   simplifications (not deviations): zero-order uniform background only (no higher-order Markov, so
   O/E can over/under-state on biased sequences), no closed-form p-value/E-value (self-overlap
   approximation affects only the probability statistic, not Count/E). Sources: Compeau & Pevzner
   *Bioinformatics Algorithms* Ch.2 (wikiselev wiki, rank 1) + monaLisa `getKmerFreq`/PeerJ O/E-ratio
   corroboration (rank 3). Linked new source + concept into the algorithm-validation-evidence hub
   (added MOTIF-DISCOVER to frontmatter sources + source-list + concept-list) and cross-linked both
   sibling motif concepts (known-motif-search, consensus-from-alignment) to it. index.md: +1 source
   +1 concept. Backlog: moved Motif_Discovery/Overrepresented_Kmer_Discovery.md pending→covered
   (70→71 covered / 175→174 pending; §Motif_Discovery 3→2, domains still 29). Contradictions: none.
   Follow-up: remaining Motif_Discovery units (Regulatory_Elements, Shared_Motifs / FindSharedMotifs)
   warrant their own pages when ingested; other de novo families (greedy/median-string/Gibbs) not
   implemented.
   graph: +2 nodes (source + concept), +2 typed edges (overrepresented-kmer-discovery relates_to
   test-unit-registry + alternative_to known-motif-search); body [[wikilinks]] mentions auto-derived.

## [2026-07-09] ingest | MOTIF-GENERATE-001-Evidence.md → motif-generate-001-evidence (source) + 1 concept
   IUPAC-Degenerate Consensus Generation (MotifFinder.GenerateConsensus): per-column keep every
   base with count > 0.25·n (strict >) → NC-IUB 1984 IUPAC symbol for that base set; no-pass
   fallback → most-frequent (alphabetical tie). Created concept iupac-degenerate-consensus and
   cross-linked the plurality [[consensus-from-alignment]] (str_replace: named GenerateConsensus
   as MOTIF-GENERATE-001 with a wikilink) and the exact [[known-motif-search]] / de-novo
   [[overrepresented-kmer-discovery]] siblings. index.md: +1 source +1 concept. Backlog: moved
   Pattern_Matching/IUPAC_Degenerate_Consensus.md pending→covered (71→72 covered / 174→173
   pending; §Pattern_Matching 8→7, domains still 29). Contradictions: none.
   Follow-up: sibling degenerate units IUPAC_Degenerate_Matching (scanning direction) +
   Position_Weight_Matrix (CreatePwm) still pending; warrant their own pages when ingested.
   graph: +2 nodes (source + concept), +2 typed edges (iupac-degenerate-consensus relates_to
   test-unit-registry + alternative_to consensus-from-alignment); body [[wikilinks]] mentions auto-derived.

## [2026-07-09] ingest | MOTIF-REGULATORY-001-Evidence.md → motif-regulatory-001-evidence (source) + regulatory-element-detection (concept)
   Ingested the Regulatory-Elements evidence unit: scanning a DNA sequence against a curated
   `KnownMotifs` catalog of 12 canonical regulatory consensus strings (TATA/−10/−35/CAAT/GC
   promoter boxes, Kozak + Shine-Dalgarno translation signals, poly(A), E-box/AP-1/NF-κB/CREB
   TF sites), each source-anchored to its primary literature; reports Name/Pattern/Sequence per
   occurrence at 0-based start, mixes exact + one IUPAC-degenerate (E-box `CANNTG`) match.
   Decision: created a DEDICATED concept [[regulatory-element-detection]] rather than enriching
   [[known-motif-search]] — the unit's correctness is its *cited catalog of named biological
   elements* (incl. the AP-1 `TGAGTCA`→`TGACTCA` corrected-defect regression), a fixed-catalog
   specialization of the generic caller-supplied exact scan. Cross-linked as the canonical-catalog
   sibling of [[known-motif-search]] and the matching-in-practice counterpart of the generation
   [[iupac-degenerate-consensus]] (both concept pages updated with back-links). index.md: +1
   source +1 concept. Backlog: moved Motif_Discovery/Regulatory_Elements.md pending→covered
   (72→73 covered / 173→172 pending; §Motif_Discovery 2→1, domains still 29). Contradictions: none
   (two source-backed representative-site assumptions: NF-κB strong site `GGGACTTTCC`, Kozak exact
   `GCCGCCACCATGG`). Follow-up: distinct promoter-detection unit (Annotation/Promoter_Detection.md)
   + Motif_Discovery/Shared_Motifs.md still pending.
   graph: +2 nodes (source + concept), +2 typed edges (regulatory-element-detection relates_to
   test-unit-registry + relates_to known-motif-search); body [[wikilinks]] mentions auto-derived.

## [2026-07-09] ingest | MOTIF-SHARED-001-Evidence.md → motif-shared-001-evidence (source) + shared-motifs (concept)
   Shared motifs across a sequence set (`FindSharedMotifs`): the van Helden / RSAT oligo-analysis
   **"matching sequences"** quorum — enumerate every fixed-`k` exact word across a *set* of sequences
   and report each word present in ≥ `minSequences` of them, keyed by presence/absence per sequence
   (a within-sequence repeat contributes 1, not its occurrence multiplicity), each carrying its
   `SequenceIndices` set + `Prevalence`=matching/total. Decision: created a DEDICATED concept
   [[shared-motifs]] rather than enriching [[longest-common-substring]] — the source explicitly
   contrasts this fixed-k + quorum + ALL-qualifying-words method against the ROSALIND LCSM framing
   (variable-length single longest substring present in *all*, via generalized suffix tree), which it
   does NOT implement. Modeled as `alternative_to` [[longest-common-substring]] (the k-string quorum
   vs single-longest-in-all pair) and `relates_to` [[overrepresented-kmer-discovery]] (same van Helden
   word-enumeration family: cross-sequence quorum vs single-sequence O/E enrichment). Cross-linked both
   ways: enriched overrepresented-kmer-discovery (wikilinked its FindSharedMotifs mention) and
   longest-common-substring (added a many-string-relative nav paragraph). Oracle: S0=`ATGATG`/
   S1=`ATGCCC`/S2=`CCCGGG`, k=3 minSeq=2 → `ATG`{0,1}(2/3)/`CCC`{1,2}; Rosalind GATTACA/TAGACCA/ATACA
   contrast (all-2-mers-in-all vs single LCSM `AC`). Corner cases: within-seq repeat→1, below-quorum
   excluded, k>shortest→no words, empty→none, k<1→throws; exact-word only (Das & Dai "no variations").
   Sources: RSAT oligo-analysis manual (rank 3, reference impl — verbatim matching-sequences/occurrence
   defs) + Das & Dai 2007 (rank 1, word-enumeration family) + van Helden/André/Collado-Vides 1998 (rank
   1 primary, HTTP 403) + Rosalind LCSM (rank 4, contrast-only). Linked new source + concept into the
   algorithm-validation-evidence hub (added MOTIF-SHARED to frontmatter sources + source-list +
   concept-list). index.md: +1 source +1 concept. Backlog: moved Motif_Discovery/Shared_Motifs.md
   pending→covered (73→74 covered / 172→171 pending; §Motif_Discovery removed, domains 29→28).
   Contradictions: none; deviations None — two presentation/API assumptions (default k=6/minSeq=2,
   Prevalence as fraction). Follow-up: ProteinMotif/Common_Motif_Finding + Motif_Search still pending
   (protein-side motif family).
   graph: +2 nodes (source + concept), +3 typed edges (shared-motifs relates_to test-unit-registry +
   alternative_to longest-common-substring + relates_to overrepresented-kmer-discovery); body
   [[wikilinks]] mentions auto-derived.

## [2026-07-09] lint | structural + graph + semantic pass (172 pages)
Structural: 2 orphans → 1 fixed (mutation-testing now linked from [[validation-and-testing]] and
[[mutation-testing-analysis]]); backlog orphan is intentional (index). Broken wikilink fixed
(methylation-context-classification anchored link to bisulfite-methylation-calling — the lint resolver
does not support `#anchor` syntax, so dropped to a plain link + prose section reference). Oversize:
backlog.md 453 lines (soft cap only, working coverage tracker — left as-is). Stale: none.
Graph: 1 broken source ref fixed by creating the missing companion source page
[[mutation-testing-analysis]] for docs/Evidence/MUTATION-TESTING-ANALYSIS.md (the mutation-testing
concept had been ingested without it); edge source: mutation-testing-analysis now resolves.
Graph re-extracted: +2 nodes, +12 edges; graph lint clean.
Coverage: 324 uncovered under docs/** (172 algorithms tracked in backlog + 132 Evidence = active
per-unit campaign + 10 checklists + others) — not triaged item-by-item (that many is the "lint report
too long" signal; see recommendations). No source contradictions found in the semantic pass over the
recently-updated motif/epigenetics/testing pages.

## [2026-07-09] ingest | docs/checklists/*.md → 10 testing-methodology checklists (9 concepts + 10 sources)
Ingested the full 10-doc testing-methodology family as a coherent batch. New concepts:
[[property-based-testing]], [[metamorphic-testing]], [[fuzzing]], [[snapshot-testing]],
[[algebraic-testing]], [[architecture-testing]], [[differential-testing]],
[[combinatorial-testing]], [[characterization-testing]] (mutation already had [[mutation-testing]],
now enriched with the checklist end-state). New source pages: one per checklist
(*-checklist slugs). Wired all nine concepts + mutation into the [[validation-and-testing]] hub
bullet list (each now links its concept + P0–P3 priority) and refreshed the coverage paragraph:
per-checklist end-state (property/metamorphic/fuzzing 258/258, architecture 22/22, combinatorial
193, mutation all-files-≥80% by 2026-06-30, algebraic 89+169-N/A, differential 107) supersedes the
older "only architecture complete" 2026-03-19 baseline — a temporal progression, no contradiction.
Real remaining gap: snapshot 37/255 + on-demand characterization. Semantic note recorded on
[[mutation-testing]]: the [[mutation-testing-analysis]] 60.6% baseline (2026-02-14) and the
[[mutation-testing-checklist]] ≥80% end-state (2026-06-30) are two points in time, not a conflict.
graph: +21 nodes, +126 edges (10 typed relates_to edges to validation-and-testing/property-based/
layered-architecture/snapshot; rest are body-wikilink mentions); graph lint clean. index updated
(10 sources + 9 concepts). Wiki now 191 pages.

## [2026-07-09] ingest | docs/Evidence/ONCO-ACTION-001-Evidence.md → clinical actionability (OncoKB levels)
First Oncology-family unit. New source page [[onco-action-001-evidence]] and new anchor concept
[[clinical-actionability-oncokb-levels]] (Clinical Actionability Assessment by the OncoKB Therapeutic
Levels of Evidence). The algorithm is a pure level-ranking of caller-supplied leveled drug associations
under the fixed combined order R1 > 1 > 2 > 3A > 3B > 4 > R2 (sensitivity axis 1 > 2 > 3A > 3B > 4,
resistance axis R1 > R2), reporting the max per axis + combined, or NotActionable when a variant carries
no leveled association. Genuinely distinct from all existing concepts (no oncology page existed) →
warranted its own concept, wired into the [[algorithm-validation-evidence]] hub (frontmatter source +
evidence link + anchor bullet). Sources: Chakravarty 2017 OncoKB (JCO PO, DOI paywalled) + OncoKB
Levels-of-Evidence PDF V2 + OncoKB Curation SOP v3 + oncokb-annotator README — all mutually consistent
(SOP explicitly consistent with AMP/ASCO/CAP Li 2017). Two assumptions: NotActionable is the library's
name for OncoKB's empty-HIGHEST_LEVEL observable; the knowledgebase is a caller input (library ranks,
does not embed the OncoKB DB). No contradictions. index updated (1 source + 1 concept).
graph: +2 nodes, +1 typed edge (relates_to → test-unit-registry on the concept); graph lint clean.

## [2026-07-09] ingest | docs/Evidence/ONCO-ANNOT-001-Evidence.md → onco-annot-001-evidence (source) + 1 concept
   Second Oncology unit: Cancer-Specific Variant Annotation by the AMP/ASCO/CAP 2017 four-tier
   clinical-significance classification (AnnotateCancerVariants + GetCOSMICAnnotation). Created concept
   cancer-variant-tier-classification-amp-asco-cap (decision rule: evidence level A/B→Tier I, C/D→Tier II,
   no-level+MAF≥1%-or-no-assoc→Tier IV, no-level+rare+assoc→Tier III; 1% primary benign cutoff inclusive;
   evidence level dominates frequency; GetCOSMICAnnotation = null-on-miss caller-supplied catalog lookup).
   Sources: Li MM et al. 2017 (J Mol Diagn, four-tier consensus, Figure 2 / Tables 3-7) + Tate JG et al.
   2019 (COSMIC external DB) — mutually consistent. Cross-linked as the sibling of, and consistent with,
   clinical-actionability-oncokb-levels (OncoKB levels). Two assumptions (caller-supplied evidence inputs;
   III/IV discriminator = direct Figure 2/Table 6-7 reading). No contradictions. index + hub updated
   (1 source + 1 concept).
   graph: +2 nodes, +2 typed edges (relates_to → test-unit-registry, relates_to → clinical-actionability-oncokb-levels on the concept); graph lint clean.

## [2026-07-09] ingest | docs/Evidence/ONCO-ARTIFACT-001-Evidence.md → onco-artifact-001-evidence (source) + 1 concept
   Third Oncology unit: Sequencing Artifact Detection (FilterArtifacts) — OxoG / FFPE deamination
   substitution classification + strand-orientation bias. Created concept sequencing-artifact-detection,
   deliberately framed as the QC sibling of the two clinical-significance ONCO units (it removes
   false-positive somatic calls from DNA damage / mapping bias BEFORE clinical interpretation, rather than
   judging significance). Three disjoint signals: (1) substitution-class — OxoG oxidation G>T(R1)/C>A(R2)
   [Chen 2017] vs FFPE cytosine-deamination C>T/G>A [Do & Dobrovic 2015], else not-an-artifact; (2) GIV
   (Global Imbalance Value) = per-substitution R1/R2 count ratio (GIV_G_T = count(G>T in R1)/count(G>T in
   R2)), neutral 1 / damaged > 1.5 [Chen 2017 + Ettwiller Damage-estimator]; (3) FisherStrand FS =
   -10*log10(two-sided Fisher-exact p) on the [ref_fwd,ref_rev,alt_fwd,alt_rev] 2x2 table, MIN_PVALUE
   1e-320 [GATK]. Oracles: GIV 200/100->2.0 & balanced->1.0; FS [10,10,10,10]->0.0 & [20,0,0,20]->large;
   class table G>T/C>A->OxoG, C>T/G>A->FFPE, A>G->neither. Result subset of input. Two assumptions:
   no BAM parser (per-strand/read-mate evidence passed on the variant record, API-shape only); GIV 1/1.5
   thresholds verbatim from the Nature Methods summary of Chen 2017. No source contradictions — the four
   sources each cover a disjoint signal and are mutually consistent. Wired into algorithm-validation-evidence
   hub (frontmatter source + evidence link + anchor bullet); index updated (1 source + 1 concept).
   graph: +2 nodes, +1 typed edge (relates_to → test-unit-registry on the concept)

## [2026-07-09] ingest | docs/Evidence/ONCO-ASCAT-001-Evidence.md → onco-ascat-001-evidence (source) + 1 concept
   Fourth Oncology unit: allele-specific copy number + joint tumor purity/ploidy fit — the upstream
   copy-number layer beneath the three clinical-interpretation ONCO units. New concept
   allele-specific-copy-number-ascat spanning FOUR disjoint algorithm stages with disjoint primary
   literature: (1) ASCAT core (Van Loo 2010 PNAS + ascat.runAscat.R) — nA/nB inversion from per-locus
   logR r + BAF b, joint (ρ,ψ) grid search minimising length-weighted squared minor-allele distance to
   non-negative integers (BAF=0.5 down-weighted x0.05), GoF=(1−d/TheoretMaxdist)*100, round+clamp-0,
   major=larger, γ=1 for sequencing (0.55 arrays only); (2) ASPCF segmentation (Nilsen 2012 PCF
   `Σ(y−ȳ)²+γ|S|` O(n²) DP e_k=min_j(d_jk+e_{j−1}+γ) + Ross 2021 joint common-breakpoint separate-means
   + BAF mirroring); (3) subclonal two-state Battenberg (Nik-Zainal 2012) n_obs=f·n₁+(1−f)·n₂ over
   bracketing integers ⌊⌋/⌈⌉, integer→single clonal state; (4) multiplicity/CCF (McGranahan 2016 /
   PICTograph VAF=(m·CCF·p)/(c·p+2(1−p)) / DeCiFering c=(F·v)/(ρ·M)), clamp m to [1,major-CN]. Planted
   oracles invert the forward model: ρ₀=0.80, ψ₀∈{2,3}, segments 1+1/2+0(CN-LOH)/2+1, clonal CCF≈1.0;
   ASPCF two-level track γ=0.5→1 breakpoint; subclonal 1.4/0.6→states (2,0)/(1,1) f≈0.4. Four
   synthesis-only/scope assumptions (het-SNP BAF forward model + avg-ploidy logR normalisation used only
   to synthesise inputs; γ exposed not hard-coded; two-state uses bracketing integers, ≥3 populations
   out of scope). Genuinely distinct from the total-CN chromosome-arm aneuploidy-detection (no allelic
   contrast/purity) — cross-linked as its allele-specific counterpart. No source contradictions (four
   disjoint stages). Wired into algorithm-validation-evidence hub (frontmatter source + evidence link +
   anchor bullet); index updated (1 source + 1 concept).
   graph: +2 nodes, +1 typed edge (relates_to → test-unit-registry on the concept)

## [2026-07-09] ingest | docs/Evidence/ONCO-CCF-001-Evidence.md → onco-ccf-001-evidence (source) + 1 concept
   Fifth Oncology unit: cancer cell fraction (CCF) estimation + 1D clonal/subclonal clustering — the
   downstream clonal-structure layer above the ASCAT copy-number substrate. The CCF point formula is
   already carried by allele-specific-copy-number-ascat §4, so this reuses/cross-links ASCAT for it; the
   genuinely distinct, wiki-worthy content is the standalone EstimateCCF with the reported-value [0,1]
   cap (exposing uncapped raw — CNAqc 1.06 noise case) and ClusterCCFValues, a deterministic 1D Lloyd
   k-means (quantile seeding, no RNG) that deconvolutes the CCF vector into clones/subclones with the
   highest-centroid = clonal rule (Tarabichi 2021). New concept cancer-cell-fraction-clonal-clustering.
   Sources corroborate the CCF closed form three ways (Tarabichi 2021 Nat. Methods / Zheng 2022
   PICTograph / McGranahan 2016 Science) + CNAqc (CCF>1 from noise) + Lloyd 1982. Two source-consistent
   assumptions ([0,1] cap via invariant + McGranahan clonal definition; Lloyd k-means as the concrete 1D
   method — sources name clustering only broadly). No contradictions. Enriched ASCAT §4 with a forward
   cross-link to the clustering concept; wired into algorithm-validation-evidence hub (frontmatter
   source + evidence link + anchor bullet); index updated (1 source + 1 concept).
   graph: +2 nodes, +2 typed edges (relates_to → test-unit-registry, depends_on → allele-specific-copy-number-ascat, on the concept)

## [2026-07-09] ingest | ONCO-CHIP-001-Evidence.md → onco-chip-001-evidence (source) + clonal-hematopoiesis-cfdna-filtering (concept)
   Sixth Oncology unit: clonal-hematopoiesis (CHIP) filtering for cfDNA liquid biopsy — the pre-interpretation
   biological-origin filter. Sources Steensma 2015 (CHIP def: VAF ≥ 2% + driver gene + no malignancy) + Genovese
   2014 (recurrent CH genes) + Razavi 2019 (CH = dominant cfDNA confounder 81.6%/53.2%, matched-WBC = definitive
   origin test) + Arango-Argoty 2025 (gold standard) + Bolton 2020 (strict origin: WBC VAF ≥ 2% AND ≥ 10 reads
   AND ≥ φ× tumour VAF, φ=2.0 / 1.5 lymph node). Three methods IdentifyCHIPVariants / FilterCHIP (matched-WBC +
   conservative gene+VAF fallback) / CallVariantOrigin. Two source-consistent assumptions (canonical default gene
   set, ≥1-alt-read WBC presence test); no contradictions. New concept cross-linked as the biological-origin
   sibling of sequencing-artifact-detection (reciprocal body link added there); wired into the
   algorithm-validation-evidence hub (frontmatter source + evidence link + anchor bullet); index updated (1 source
   + 1 concept).
   graph: +2 nodes, +2 typed edges (relates_to → test-unit-registry, relates_to → sequencing-artifact-detection, on the concept)

## [2026-07-09] ingest | ONCO-CLONAL-001-Evidence.md → onco-clonal-001-evidence (source) + clonal-subclonal-classification-ccf-posterior (concept)
   Seventh Oncology unit: clonal vs subclonal mutation classification via a Bayesian CCF posterior — the
   probabilistic clonal-structure classifier. Sources Landau 2013 Cell (ABSOLUTE-style expected allele fraction
   f(c)=αc/(2(1−α)+αq), posterior P(c)∝Binom(a|N,f(c)) uniform prior on 100-point grid c∈[0.01,1], rule clonal iff
   P(CCF>0.95)>0.5) + Satas 2021 Cell Systems DeCiFering (multiplicity-general f(c)=αMc/(2(1−α)+αq), Eq. 1). Grid
   oracles A1/B2/E clonal, C1/D subclonal, E the M=2 multiplicity lift; point-estimate IdentifyClonalMutations
   strict CCF>0.95 → indices {0,2,4}. One API-shape assumption (per-variant local copy number q over a genome-wide
   ploidy scalar), no source contradictions. Judged genuinely DISTINCT from ONCO-CCF-001 (point estimate + Lloyd
   k-means clustering): new concept created and cross-linked alternative_to cancer-cell-fraction-clonal-clustering
   (reciprocal body link added there); wired into the algorithm-validation-evidence hub (frontmatter source +
   evidence link + anchor bullet); index updated (1 source + 1 concept).
   graph: +2 nodes, +3 typed edges (relates_to → test-unit-registry, alternative_to → cancer-cell-fraction-clonal-clustering, depends_on → allele-specific-copy-number-ascat, on the concept)

## [2026-07-09] ingest | ONCO-CNA-001-Evidence.md → onco-cna-001-evidence (source) + copy-number-alteration-classification (concept)
   Eighth Oncology unit: copy-number alteration classification — a single log2 copy ratio → absolute integer
   CN (n=2·2^log2, CNVkit `_log2_ratio_to_absolute_pure`, diploid ref_copies=2) → discrete CNA state via
   CNVkit `absolute_threshold` hard-threshold caller (default −1.1/−0.25/0.2/0.7 → DeepDeletion/Loss/Neutral/
   Gain/Amplification; first `log2<=thresh` boundary-inclusive→lower bin; above-last→ceil(2·2^log2); NaN→
   neutral CN). Sources CNVkit call.py + docs (germline −0.4/0.3 vs tumor −0.25/0.2, purity≥30% caveat) +
   GISTIC2 Mermel 2011 (±0.1 noise band + +0.848/−0.737 high-amplitude cutoffs) + GISTIC2 -ta/-td docs +
   SV-CNV-001 in-repo overlap check. Judged genuinely DISTINCT from allele-specific ONCO-ASCAT-001 (no allelic
   contrast / purity fit) and from whole-chromosome CHROM-ANEU-001 (per-segment 5-state oncology call vs
   ≥80%-bin chromosome vote) — both share only the n=2·2^log2 conversion; and from SV-CNV-001's round-based
   integer CN (no state classification). New concept created and cross-linked; wired into the
   algorithm-validation-evidence hub (frontmatter source + evidence link + anchor bullet); index updated
   (1 source + 1 concept). One diploid-ploidy=2 assumption, no source contradictions.
   graph: +2 nodes, +1 typed edge (relates_to → test-unit-registry, on the concept)

## [2026-07-09] ingest | ONCO-CNA-002-Evidence.md → onco-cna-002-evidence (source) + focal-amplification-detection (concept)
   Ninth Oncology unit: focal amplification detection — a two-part predicate `DetectFocalAmplifications`
   keeps segments both amplified (log2 gain > GISTIC2 t_amp 0.1) AND focal (SegLen/ArmLength <
   broad_len_cutoff 0.98 — Mermel 2011's length-based focal/arm-level split; strict < 0.98, exactly 0.98
   → arm-level), then `IdentifyAmplifiedOncogenes` maps each focal amp's arm prefix to a built-in oncogene
   panel (17q→ERBB2, 8q→MYC, 7p→EGFR, 11q→CCND1, 12q→MDM2 AND CDK4, NCBI Gene cytobands). Sources GISTIC2
   Mermel 2011 (Genome Biology, length rule) + GISTIC2 docs (broad_len_cutoff 0.98 / t_amp 0.1) + CNVkit
   (single-copy gain log2(3/2)=0.585 > 0.1 → 0.1 admits all gains) + NCBI Gene oncogene cytobands. Judged
   genuinely DISTINCT from ONCO-CNA-001 (log2→5-state classification): it asks the orthogonal LENGTH
   question and maps to oncogenes, sharing only the GISTIC2 t_amp=0.1 amplitude gate — new concept created
   and cross-linked (relates_to copy-number-alteration-classification, reciprocal body link added there).
   Worked oracles A 17q 0.50/log2 1.0→ERBB2, B 8q 0.99→arm-level, C 7p log2 0.05→not amplified, D 11q 0.98
   boundary→no. Two assumptions (amplitude+length fusion = integration choice; caller supplies arm label +
   length, no cytoband table); deletions out of scope (ONCO-CNA-003); no source contradictions. Wired into
   the algorithm-validation-evidence hub (frontmatter source + evidence link + anchor bullet); index updated
   (1 source + 1 concept).
   graph: +2 nodes, +2 typed edges (relates_to → test-unit-registry, relates_to → copy-number-alteration-classification, on the concept)

## [2026-07-09] ingest | docs/Evidence/ONCO-CNA-003-Evidence.md → homozygous / deep deletion detection (tenth Oncology unit)
   Homozygous (deep) deletion detection, the deletion mirror of ONCO-CNA-002: filter segments whose
   classified integer copy number is exactly 0 (homozygous / deep deletion), then
   `IdentifyDeletedTumorSuppressors` maps each arm prefix to a built-in tumour-suppressor panel
   (17p→TP53, 13q→RB1 AND BRCA2, 9p→CDKN2A, 10q→PTEN, 17q→BRCA1, NCBI Gene cytobands). Sources cBioPortal
   file-format + FAQ (−2 = "Deep Deletion, possibly a homozygous deletion"; −1 = shallow/heterozygous) +
   Cheng et al. 2017 Nat Commun (homozygous deletion = total copy number 0, "zero copies of both alleles",
   two hits, targets tumour suppressors) + CNVkit `absolute_threshold` (integer CN 0 ⇒ DeepDeletion —
   REUSES ONCO-CNA-001, no new threshold) + NCBI Gene tumour-suppressor cytobands. Judged genuinely
   DISTINCT and wiki-worthy — it is the loss-side counterpart of the amplification unit
   focal-amplification-detection (IdentifyDeletedTumorSuppressors mirrors IdentifyAmplifiedOncogenes) and
   a consumer of copy-number-alteration-classification's CN-0/DeepDeletion state; new concept created and
   cross-linked (relates_to test-unit-registry, copy-number-alteration-classification, and
   focal-amplification-detection; reciprocal body links + a reciprocal relates_to edge added on
   focal-amplification-detection). Oracles CN 0→homozygous→gene, CN 1 single-copy loss→not, neutral/gain/
   amp→not, boundary log2 −1.1 inclusive→CN 0. Two assumptions (CN-0 reuse of ONCO-CNA-001, caller-fixed
   tumour-suppressor panel); no source contradictions. Wired into the algorithm-validation-evidence hub
   (frontmatter source + evidence link + anchor bullet); index updated (1 source + 1 concept).
   graph: +2 nodes, +4 typed edges (concept relates_to test-unit-registry + copy-number-alteration-classification + focal-amplification-detection; focal-amplification-detection relates_to homozygous-deletion-detection)

## [2026-07-10] ingest | docs/Evidence/ONCO-CTDNA-001-Evidence.md → ctDNA detection + tumor-fraction (eleventh Oncology unit)
   ctDNA analysis (liquid-biopsy quantification / limit-of-detection layer): the Poisson detection
   probability `DetectionProbability` p = 1 − e^(−n·d·k) (n genome equivalents, d mutant allele
   fraction, k reporters) with a detectability test (caller threshold default 0.95 AND λ = n·d·k ≥ 1;
   only p is non-assumption), `CalculateTumorFraction` = 2 × mean clonal-heterozygous VAF (copy-neutral
   diploid, v = π/2), `CalculateMeanVaf` = mean altReads/totalReads across reporters, and a
   genome-equivalents helper (3.3 pg/haploid ⇒ ≈303 GE/ng). Sources Newman 2014 CAPP-Seq (detection
   0.025%–10%, 96% specificity ~0.02%, background 0.006%/0.0003%, across-reporter fraction) + Patent US
   11,085,084 restating Avanzini 2020 Sci. Adv. (Poisson λ=n·d, low-burden λ<3) + Pessoa 2023 (λ=15,000
   ×0.001=15) + Devonshire 2014 (3.3 pg/haploid) + Alcaide 2020 (303 GE/ng) + CNAqc/Antonello 2024
   (TF=2·VAF). Judged genuinely DISTINCT from the sibling clonal-hematopoiesis-cfdna-filtering (which
   FILTERS non-tumor cfDNA calls) — this QUANTIFIES the tumor signal on the same cfDNA input, so a new
   concept was created and cross-linked (relates_to test-unit-registry + clonal-hematopoiesis-cfdna-
   filtering; reciprocal body link added on the CHIP page). Oracles n=15000,d=0.001,k=1→1−e⁻¹⁵≈
   0.99999969, λ=0→p=0 not-detected, TF 0.10→0.20, 303 GE/ng, 3.3 pg→1 GE. One flagged
   detection-threshold assumption; no source contradictions (seven references cover disjoint stages).
   Wired into the algorithm-validation-evidence hub (frontmatter source + evidence link + anchor bullet);
   index updated (1 source + 1 concept).
   graph: +2 nodes, +2 typed edges (concept relates_to test-unit-registry + clonal-hematopoiesis-cfdna-filtering)

## [2026-07-10] ingest | docs/Evidence/ONCO-DRIVER-001-Evidence.md → driver-gene-classification-20-20-rule (twelfth Oncology unit)
   Driver Mutation Detection, the Vogelstein 2013 20/20 rule — a per-gene mutation-pattern heuristic
   classifying a cancer gene Oncogene (> 20% missense at recurrent positions, recurrent = same protein
   position ≥ 2×, = activating), TumorSuppressor (> 20% truncating/inactivating — nonsense, frameshift,
   splice donor/acceptor, gained/lost stop = loss of function), or Ambiguous (neither criterion, or exact
   dual-pass tie). Methods IdentifyDriverMutations (driver ⊆ somatic), MatchCancerHotspots
   (caller-supplied (gene, position) hotspot set), ScoreDriverPotential (= max of the two criterion
   fractions in [0,1]; CADD/SIFT/PolyPhen are externally trained models → caller-supplied, not
   implemented). Sources Vogelstein 2013 Science "Cancer Genome Landscapes" (originating source; PMC
   CAPTCHA + DOI 403 so wording taken verbatim from three open-access secondaries) + Tokheim & Karchin
   2020 20/20+ (verbatim rule, inactivating = nonsense/frameshift) + Schroeder 2014 OncodriveROLE
   (truncating list = frameshift / gained-or-lost stop / splice donor-acceptor; writes "≥20%") + Miller
   2017 (recurrent = ≥2×, IDH1 codon 132 R132H). Oracles IDH1 10 missense@codon132 → recurrent-missense
   1.00 → Oncogene; dispersed 8/10 truncating → 0.80 → TumorSuppressor; truncating exactly 0.20 → NOT TSG
   (strict >). Judged genuinely DISTINCT and wiki-worthy — a GENE-level driver classifier orthogonal to
   the VARIANT-level clinical classifiers cancer-variant-tier-classification-amp-asco-cap and
   clinical-actionability-oncokb-levels (body-linked as context, not typed edges); a heuristic not a
   statistical test (passenger truncations + low-recurrence drivers mislead it; 20/20+ / MutSigCV
   successors out of scope). Three assumptions: strict > 0.20 for both (Vogelstein/Tokheim ">20%" over
   OncodriveROLE's "≥20%" — the sole glyph difference, resolved to strict); dual-pass tie-break by larger
   fraction, Ambiguous on exact tie; ScoreDriverPotential = max-of-fractions proxy. No source
   contradictions. New concept + source created, wired into the algorithm-validation-evidence hub
   (frontmatter source + evidence link + anchor bullet); index updated (1 source + 1 concept).
   graph: +2 nodes, +1 typed edge (concept relates_to test-unit-registry)

## [2026-07-10] ingest | docs/Evidence/ONCO-EXPR-001-Evidence.md → onco-expr-001-evidence (source) + 1 concept
   ONCO-EXPR-001 = Tumor Gene Expression Outlier (z-score) + Signature Score, the thirteenth ingested
   Oncology unit and the wiki's first expression/transcriptome method. Per-gene outlier z = (r−μ)/σ
   against a caller-supplied reference (base) population (cBioPortal diploid or all-samples), with
   σ = sample SD divisor (n−1) — settled by the reference `NormalizeExpressionLevels.java` `std()` over
   the prose spec's silence — classified over/under-expressed under the strict ±2 default threshold
   (exactly ±2 NOT an outlier). Combined-z signature/pathway activity a = (Σzᵢ)/√k (Lee et al. 2008,
   GSVA `zscore` method, corroborated by the GSVA vignette). Zero-SD reference throws (`fatalError`),
   a behavioural deviation from the prose "z ← NA when SD = 0" — the two cBioPortal sources disagree
   and the code wins. Oracles: reference {2,2,4,6,6}→μ=4/σ=2, x=10→3.0 over / x=8→2.0 boundary-not-
   outlier / x=4→0.0 / x=−1→−2.5 under; signature {3,1,−1,1}→a=2.0, single-gene {2.5}→2.5. Corner
   cases: n≤1 SD undefined, k=0 invalid, k=1 well-defined. Two scope assumptions (caller-supplied
   cohort+signature; inputs pre-normalized / z meaningful) + one behavioural deviation (throw not NA);
   no further contradictions (z formula corroborated four ways). New concept
   [[expression-outlier-zscore-signature-score]] + source [[onco-expr-001-evidence]] created, wired
   into the algorithm-validation-evidence hub (frontmatter source + evidence link + anchor bullet);
   index updated (1 source + 1 concept).
   graph: +2 nodes, +1 typed edge (concept relates_to test-unit-registry)

## [2026-07-10] ingest | docs/Evidence/ONCO-FUSION-001-Evidence.md → onco-fusion-001-evidence (source) + 1 concept
   ONCO-FUSION-001 = Fusion Gene Detection (candidate fusion calling from breakpoint-supporting
   reads), the fourteenth ingested Oncology unit and the wiki's first gene-fusion / read-evidence
   structural-rearrangement method. Genuinely distinct from all existing ONCO concepts (copy-number,
   clonal, expression, clinical-interpretation) and from the gene-order signed-permutation
   [[genome-rearrangement-breakpoint-distance]] → new concept warranted. The STAR-Fusion / Arriba
   split-read + discordant-pair + minimum-support paradigm, corroborated across two independent tools
   + their papers (Haas 2017/2019, Uhrig 2021), no contradictions. Detection rule: DETECTED iff
   (junction ≥ MIN_JUNCTION_READS=1 AND total ≥ MIN_SUM_FRAGS=2) OR (zero junction AND discordant ≥
   MIN_SPANNING_FRAGS_ONLY=5), with total support = split_reads1+split_reads2+discordant_mates
   (Arriba) and the gene5p ≠ gene3p distinct-gene invariant; results ordered by descending support.
   Separate exon-phase in-frame check (5' coding bases − 3' start phase) mod 3 == 0 (Genomics England
   / Wikipedia Reading-frame primary cites). Oracles: EML4-ALK(3,2,4)/TMPRSS2-ERG(1,0,1)/CD74-ROS1
   (0,0,5) DETECTED, NCOA4-RET(0,0,4) span<5 / KIF5B-RET(1,0,0) sum<2 / ALK-ALK same-gene REJECTED,
   frame 300/0→in 301/0→out 301/1→in. Two scope assumptions (candidate-level counts not raw BAM —
   chimeric-read extraction is a separate FindChimericReads; phase-only in-frame, no premature-stop
   scan = ONCO-FUSION-003). New concept [[gene-fusion-detection-read-evidence]] + source
   [[onco-fusion-001-evidence]] created, wired into the algorithm-validation-evidence hub (frontmatter
   source + evidence link + anchor bullet); index updated (1 source + 1 concept).
   graph: +2 nodes, +1 typed edge (concept relates_to test-unit-registry)

## [2026-07-10] ingest | ONCO-FUSION-002-Evidence.md → onco-fusion-002-evidence (source) + gene-fusion-nomenclature-known-fusion-lookup (concept)
   Fifteenth Oncology unit. Known Fusion Database Lookup: HGNC gene-fusion designation
   (Bruford et al. 2021 — `::` double-colon separator, 5′-partner-first directional order, approved
   symbols, read-throughs keep the hyphen) `GetFusionAnnotation(5p,3p)="5p::3p"` + directional
   `MatchKnownFusions` against a caller-supplied set keyed by 5′::3′, case-insensitive. A Framework
   algorithm — format/keying source-backed, set contents caller-supplied (no bundled
   Mitelman/COSMIC/ChimerDB). BCR::ABL1 worked example; A::B ≠ B::A + hyphen ≠ :: corner cases.
   Distinct from the read-evidence caller ONCO-FUSION-001 (detection); this is the naming/annotation
   layer downstream of it (round-trips a DetectFusions FusionCall), distinct from the ONCO-FUSION-003
   premature-stop scope. New concept [[gene-fusion-nomenclature-known-fusion-lookup]] + source
   [[onco-fusion-002-evidence]] created, wired into the algorithm-validation-evidence hub (frontmatter
   source + anchor bullet) and cross-linked from [[gene-fusion-detection-read-evidence]]; index
   updated (1 source + 1 concept).
   graph: +2 nodes, +2 typed edges (concept relates_to gene-fusion-detection-read-evidence + test-unit-registry)

## [2026-07-10] ingest | docs/Evidence/ONCO-FUSION-003-Evidence.md → onco-fusion-003-evidence (source) + fusion-breakpoint-frame-and-protein-prediction (concept)
   Sixteenth Oncology unit, third member of the fusion trio (the protein-consequence layer both
   siblings explicitly deferred to). Fusion Breakpoint Analysis: junction reading-frame consequence
   + fusion protein prediction. Four-state BreakpointFrameStatus (InFrame/OutOfFrame/StopCodon/
   NotPredicted) via Arriba's two-way native-frame model (Uhrig 2021) — NOT AGFusion's three-way
   class; `in-frame (with mutation)` (contiguous ORF mult-of-3 but 3′ frameshifted) maps to
   OutOfFrame. In/out reuses ONCO-FUSION-001's exon-phase rule (5' coding bases − 3' start phase)
   mod 3 == 0; gated by breakpoint-site classification (CDS vs UTR/intron/intergenic → NotPredicted).
   PredictFusionProtein follows AGFusion model.py exactly: 5′ CDS prefix + 3′ CDS suffix →
   concatenate → translate (transl_table=1) → truncate at first stop (out-of-frame trims to whole
   codons first). Oracles ATGAAA|GATGGT→MKDG, ATGAAA|GATTAAGGT→MKD (premature stop), ATGA|AAGGT
   phase-0→OutOfFrame yet clean MKG (Arriba-vs-AGFusion divergence) / phase-1→in-frame. One
   API-shape assumption (caller supplies CDS strings + junction offsets, no bundled GTF); no
   contradictions. New concept [[fusion-breakpoint-frame-and-protein-prediction]] + source
   [[onco-fusion-003-evidence]] created, wired into the algorithm-validation-evidence hub
   (frontmatter source + anchor bullet) and cross-linked from both [[gene-fusion-detection-read-evidence]]
   (two deferral references now resolve to it) and [[gene-fusion-nomenclature-known-fusion-lookup]];
   index updated (1 source + 1 concept).
   graph: +2 nodes, +2 typed edges (concept relates_to gene-fusion-detection-read-evidence + test-unit-registry)

## [2026-07-10] ingest | docs/Evidence/ONCO-HETERO-001-Evidence.md → onco-hetero-001-evidence (source) + 1 concept
   Seventeenth Oncology unit — Tumor Heterogeneity Analysis. Created concept
   [[intratumor-heterogeneity-metrics]]: the scalar-summary ITH-metric layer — MATH score
   100·1.4826·median(|VAF−median VAF|)/median(VAF) (Mroz & Rocco 2013 / Mroz 2015 PLOS Med /
   maftools mathScore.R, three-way identical, no clustering) + Shannon clonal diversity
   H=−Σ pᵢ ln pᵢ (natural log, Liu 2017 / Shannon 1948) + subclone count (Liu richness = occupied
   CCF clusters) + subclonal fraction #(CCF<0.95)/n (Landau 2013 threshold). Oracles MATH 49.42
   (odd) / 59.304 (even), Shannon 0/ln2/ln4; zero-median-VAF→throw, MAD=0→MATH=0. Judged genuinely
   distinct — a metric/summary layer, NOT per-mutation reconstruction — so a dedicated concept vs
   reusing the CCF clustering / posterior units; it depends_on [[cancer-cell-fraction-clonal-clustering]]
   (subclone count + Shannon pᵢ consume its clusters) and reuses the 0.95 threshold of
   [[clonal-subclonal-classification-ccf-posterior]]. Two source-consistent assumptions (Shannon pᵢ =
   per-cluster mutation proportions; R even-count median), no contradictions. Wired into the
   algorithm-validation-evidence hub (frontmatter source + list link + anchor bullet), cross-linked
   from cancer-cell-fraction-clonal-clustering (downstream-summary note), index updated (1 source +
   1 concept).
   graph: +2 nodes, +3 typed edges (concept relates_to test-unit-registry, depends_on cancer-cell-fraction-clonal-clustering, relates_to clonal-subclonal-classification-ccf-posterior)

## [2026-07-10] ingest | docs/Evidence/ONCO-HLA-001-Evidence.md → onco-hla-001-evidence (source) + 1 concept
   Eighteenth Oncology unit — HLA allele nomenclature parsing/validation + allele-specific HLA LOH
   (LOHHLA), the wiki's first HLA / immuno-oncology antigen-presentation method. Created concept
   [[hla-nomenclature-and-allele-specific-loh]]: (1) WHO IPD-IMGT/HLA name parse/validate
   `HLA-[Gene]*[F1]:[F2][:F3][:F4][suffix]` (Marsh 2010 colon-delimited fields, two-field minimum /
   four-field maximum, N/L/S/C/A/Q suffixes) and (2) LOHHLA LOH call — copy number < 0.5 AND
   allelic-imbalance paired t-test p < 0.01 (both strict <, McGranahan 2017 Cell PMC5720478 +
   mskcc/lohhla LOHHLAscript.R). Oracles: HLA-A*24:02:01:02L valid / HLA-A*02 / A*02:01 / five-fields
   / ...X rejected; (1.8,0.30,0.001)→LOH allele 2 / (1.60,0.40,0.05)→no (p≥0.01 guard) /
   (1.50,0.50,0.001)→no (0.5 not <0.5) / (1.70,0.40,0.01)→no (0.01 not <0.01). Judged genuinely
   distinct (no existing HLA/MHC/neoantigen concept) → dedicated concept, cross-linked to
   [[allele-specific-copy-number-ascat]] as its HLA-locus specialization (reciprocal link added).
   One assumption (both alleles <0.5 → HomozygousLoss label, thresholds unchanged), no contradictions.
   Wired into the algorithm-validation-evidence hub (frontmatter source + list link + anchor bullet),
   index updated (1 source + 1 concept).
   graph: +2 nodes, +2 typed edges (concept relates_to test-unit-registry, relates_to allele-specific-copy-number-ascat)

## [2026-07-10] ingest | docs/Evidence/ONCO-HRD-001-Evidence.md
   Nineteenth Oncology unit — HRD composite genomic-scar score `HRD = LOH + TAI + LST`, an unweighted
   sum of three large-scale copy-number scar counts with the HRD-high cutoff >= 42 inclusive (Telli
   2016 + Stewart 2022, independently corroborated). Created concept
   [[homologous-recombination-deficiency-score]]: all three components derived per segment from the
   [[allele-specific-copy-number-ascat]] major/minor CN substrate in `DetectHRD(segments)` — HRD-LOH
   (regions >15 Mb & < whole chromosome, exclude whole-chr LOH; Abkevich 2012 / oncoscanR / scarHRD
   calc.hrd.R, no centromere table) + TAI (imbalanced major!=minor segments reaching a sub-telomere
   not crossing the centromere; Birkbak 2012 / calc.ai_new.R, sub-1 Mb dropped, single-segment
   whole-chr imbalance excluded) + LST (adjacent >=10 Mb same-arm regions <3 Mb apart after iterative
   3 Mb smoothing; Popova 2012 / calc.lst.R, autosomes only; sum via scar_score.R). TAI+LST need the
   per-chromosome centromere acen [start,end] table embedded for GRCh38/GRCh37 (UCSC cytoBand
   cross-verified vs NCBI GRC modeled centromeres — resolving the prior "centromere table
   unretrievable" blocker). Oracles (14,14,14)->42 HRD-high (boundary) / (14,13,14)->41 negative /
   (0,0,0)->0 near-diploid. Judged genuinely distinct (composite score + its own LOH/TAI/LST defs +
   cutoff + centromere tables; no existing concept covers it) -> dedicated concept, cross-linked to
   [[allele-specific-copy-number-ascat]] as the downstream genomic-scar aggregation layer (reciprocal
   link added) and distinguished from the total-CN [[aneuploidy-detection]]. One even-ploidy AI-path
   assumption (major!=minor, ASCAT ploidy column absent), no contradictions. Wired into the
   algorithm-validation-evidence hub (frontmatter source + list link + anchor bullet), index updated
   (1 source + 1 concept).
   graph: +2 nodes, +2 typed edges (concept relates_to test-unit-registry, relates_to allele-specific-copy-number-ascat)

## [2026-07-10] ingest | docs/Evidence/ONCO-IMMUNE-001-Evidence.md — Immune Infiltration Estimation (twentieth ONCO-* unit)
   Tumor immune-microenvironment quantification: CIBERSORT linear-mixture m=S·f solved by ν-SVR
   (DeconvoluteImmuneCellsNuSvr; Newman 2015 / Schölkopf 2000 — z-standardize, sweep ν∈{0.25,0.5,0.75}
   by lowest RMSE, zero-clip + normalize Σf=1; cross-checked vs scikit-learn NuSVR + planted-truth) with
   NNLS/LLSR baseline (Abbas 2009) retained, plus ESTIMATE ssGSEA immune/stromal scoring (simplified
   rank-weighted mean) + opt-in Affymetrix-only cosine tumor-purity transform (negative→NaN), and
   MCP-counter marker geometric-mean note. LM22 (547×22) caller-supplied (Stanford licence forbids
   redistribution, no exact-CIBERSORT parity); ABIS-Seq (Monaco 2019, CC BY 4.0, 1296×17) bundled via
   LoadBundledAbisSignatureMatrix. Judged genuinely distinct (deconvolution + signature scoring, no
   existing concept covers it) -> dedicated concept [[immune-infiltration-deconvolution]], cross-linked
   to [[expression-outlier-zscore-signature-score]] (shared ssGSEA signature-scoring layer, typed edge +
   reciprocal prose link) and [[hla-nomenclature-and-allele-specific-loh]] (immuno-oncology sibling,
   reciprocal prose link). Wired into the algorithm-validation-evidence hub (frontmatter source + list
   link + anchor bullet), index updated (1 source + 1 concept). Two scope assumptions (LM22
   caller-supplied, simplified ssGSEA) + Affymetrix purity domain, no contradictions.
   graph: +2 nodes, +2 typed edges (concept relates_to test-unit-registry, relates_to expression-outlier-zscore-signature-score)

## [2026-07-10] ingest | docs/Evidence/ONCO-LOH-001-Evidence.md — Loss of Heterozygosity detection (twenty-first ONCO-* unit)
   Genome-wide LOH caller: DetectLOH counts HRD-LOH regions (minor CN 0 & major CN != 0 = allelic
   loss not homozygous deletion; length strictly > 15 Mb, length = end-start; whole-chromosome
   global-LOH chromosomes excluded via chrDel; major CN capped to 1 before merge so LOH state drives
   merging) + CalculateLOHFraction (length-weighted per-chromosome LOH fraction in [0,1]). Sources
   Abkevich 2012 Br J Cancer PMID 23047548 (HRD score = number of intermediate-size LOH regions,
   BRCA1/2-correlated) + scarHRD calc.hrd.R/scar_score.R (exact R criterion + 15 Mb cut-off + chrDel,
   sizelimitLOH=15e6) + oncoscanR score_loh (independent; adjacent/overlapping LOH merged <=1 bp).
   Synthetic 7-segment dataset -> score 1 (only 20 Mb chr1; exactly-15 Mb / het-retained /
   homozygous-deletion / whole-chr-LOH excluded), LOH-fraction 0.333/0.0/1.0. Judged distinct: the
   HRD concept treated HRD-LOH only as one summed component and pointed to DetectLOH, and
   CalculateLOHFraction is unrepresented -> dedicated concept [[loss-of-heterozygosity-detection]],
   cross-linked as the LOH term of [[homologous-recombination-deficiency-score]] (reciprocal link:
   HRD-LOH bullet now names DetectLOH's page), reads segments off [[allele-specific-copy-number-ascat]]
   (reciprocal prose link added), and genome-wide counterpart of the HLA-locus
   [[hla-nomenclature-and-allele-specific-loh]]. Two assumptions (segment input shape; length-weighted
   LOH-fraction a definitional/API choice — segment criterion fully source-backed), no contradictions.
   Wired into the algorithm-validation-evidence hub (frontmatter source + list link + anchor bullet),
   index updated (1 source + 1 concept; also backfilled the missing onco-immune-001-evidence source line).
   graph: +2 nodes, +3 typed edges (concept relates_to test-unit-registry, relates_to allele-specific-copy-number-ascat, relates_to homologous-recombination-deficiency-score)
## [2026-07-10] ingest | docs/Evidence/ONCO-MHC-001-Evidence.md — MHC-Peptide Binding prediction + binder classification (twenty-second ONCO-* unit)
   The immuno-oncology affinity gate of neoantigen candidate scoring, three layers over one
   Strong/Weak/NonBinder output space. (1) Classification: IC50 tiers <50 Strong / <500 Weak /
   >=500 NonBinder (Sette 1994 ~500 nM CTL threshold + IEDB) and NetMHCpan-4.1 %Rank tiers (class I
   SB <0.5 / WB <2; class II <2 / <10, all strict <; Reynisson 2020) + peptide-length validity
   (class I 8-14, class II 13-25). (2) Matrix prediction: BIMAS product rule (T1/2 = const * prod
   coeffs, 180 coeffs + const, unlisted residue = 1.0; Parker 1994/BIMAS) and SMM additive rule
   IC50 = 50000^(1-score) (Peters & Sette 2005 / IEDB log50k; anchors 0->50000, 0.5->223.6068,
   1->1 nM). (3) Ported MHCflurry 2.0 pan-allele class-I NN (Apache-2.0; BLOSUM62
   left_pad_centered_right_pad 45x21 peptide + 37-residue allele pseudosequence 37x21, tanh Dense
   feedforward/with-skip-connections topologies, to_ic50 = 50000^(1-x), geometric-mean ensemble;
   C# port matches mhcflurry 2.1.5 single-net + full-ensemble to <0.03%). Oracles: IC50
   10/50/200/500->Strong/Weak/Weak/NonBinder, %Rank 0.4/0.5/2.0->Strong/Weak/NonBinder, BIMAS
   LMV->90 / AAA->10, MHCflurry GILGFVFTL~19 nM / SIINFEKL~11.5 uM; strict-< boundaries. Framework
   packaging boundary: trained matrix + full 80 MB ensemble caller-supplied (BIMAS CGI defunct /
   Parker paywalled / IEDB non-commercial / repo size), one member embedded for CI + LoadWeightPack,
   ~0.7 MB pseudosequence table bundled -- explicitly analogized in-file to CIBERSORT LM22 in
   ONCO-IMMUNE-001. Class I length window resolved 11->14 (full NetMHCpan-4.1 window), propagates to
   GenerateNeoantigenPeptides (ONCO-NEO-001). Judged distinct + wiki-worthy: no MHC/neoantigen concept
   existed -> dedicated concept [[mhc-peptide-binding-prediction]], cross-linked as the presentation
   affinity gate of [[hla-nomenclature-and-allele-specific-loh]] (reciprocal prose link added: HLA LOH
   removes the platform), the packaging-boundary twin of [[immune-infiltration-deconvolution]]
   (source's own analogy), and downstream of the not-yet-ingested ONCO-NEO-001. Wired into the
   algorithm-validation-evidence hub (frontmatter source + summary link + anchor bullet); index
   updated (1 source + 1 concept). Follow-up: ONCO-NEO-001 (GenerateNeoantigenPeptides) not yet ingested.
   graph: +2 nodes, +3 typed edges (concept relates_to test-unit-registry, relates_to hla-nomenclature-and-allele-specific-loh, relates_to immune-infiltration-deconvolution)

## [2026-07-10] ingest | docs/Evidence/ONCO-MRD-001-Evidence.md → onco-mrd-001-evidence (source) + 1 concept
   Twenty-third Oncology unit: tumor-informed minimal/molecular residual disease (MRD) detection.
   Judged DISTINCT from ONCO-CTDNA-001 (multi-variant MRD verdict vs single-reporter Poisson
   probability) -> dedicated concept [[tumor-informed-mrd-detection]]. Two engines: Signatera
   panel positivity DetectMRD (>=2 of 16 tracked variants = MRD-positive; Reinert 2019 / PMC9265001
   Table 1; longitudinal TrackVariantsOverTime; panel Poisson p=1-e^(-nfm) reused from ONCO-CTDNA-001)
   + INVAR GLRT (Wan 2020 + INVAR2 verbatim: per-locus mixture q=p*g+e(1-p), EM ctDNA-fraction,
   LR=logL(p̂)-logL(0), AF/SNR-weighting, IMAFv2 background-subtracted depth-weighted, fragment-size
   weighting + opt-in Gaussian-KDE size profile, repolish outlier suppression, control-derived
   background). Cross-linked ctdna-detection-and-tumor-fraction (prose link both ways via depends_on
   edge + mention). Wired into algorithm-validation-evidence hub (frontmatter source + summary link +
   anchor bullet, source_commit bumped). Index updated (1 source + 1 concept). One flagged assumption
   (per-variant "detected"=>=1 alt read, tunable; panel >=2 rule unaffected); KDE opt-in vs discrete
   default resolved; no source contradictions.
   graph: +2 nodes, +2 typed edges (concept relates_to test-unit-registry, depends_on ctdna-detection-and-tumor-fraction)

## [2026-07-10] ingest | docs/Evidence/ONCO-MSI-001-Evidence.md → onco-msi-001-evidence (source) + 1 concept
   Twenty-fourth Oncology unit: Microsatellite Instability (MSI) detection — unstable-loci fraction
   score + status classification. Judged DISTINCT and wiki-worthy (an independent immunotherapy /
   mismatch-repair biomarker not previously represented) -> dedicated concept
   [[microsatellite-instability-detection]]. Two inputs, two classifiers: continuous MSIsensor/MSIsensor2
   score = unstable/valid loci (as %) with MSI-High cutoff >=20% inclusive (Niu-lab README; per-site
   chi-square tumor-vs-normal repeat-length test at FDR 0.05, Niu 2014; tumor-only unchanged; 3.5% cohort
   separation is dataset-specific) -> binary MSI-H vs not-High; and the categorical Bethesda marker-count
   rule over the 5-marker panel (0/5->MSS, 1/5->MSI-L, >=2/5->MSI-H; Boland 1998, cross-checked vs 2004
   revised-Bethesda fraction form BJC 2014). Key modelling choice: no MSI-L band fabricated for the
   continuous score. Cross-linked homologous-recombination-deficiency-score,
   cancer-variant-tier-classification-amp-asco-cap, clinical-actionability-oncokb-levels. Wired into
   algorithm-validation-evidence hub (frontmatter source + summary link + anchor bullet, source_commit
   bumped to ea6d7a9). Index updated (1 source + 1 concept). One flagged assumption (no MSI-L band on the
   continuous score); no source contradictions.
   graph: +2 nodes, +1 typed edge (concept relates_to test-unit-registry)

## [2026-07-10] ingest | docs/Evidence/ONCO-NEO-001-Evidence.md
   Twenty-fifth Oncology unit: Neoantigen Candidate Peptide Window Generation for a somatic missense
   mutation (`GenerateNeoantigenPeptides`) — the UPSTREAM partner of the ONCO-MHC-001 affinity gate.
   Judged DISTINCT and wiki-worthy (a windowing/enumeration algorithm not previously represented; the MHC
   page even carried a "not yet ingested" placeholder for it) -> dedicated concept
   [[neoantigen-peptide-generation]]. Given a missense SNV at 1-based position p in a length-L protein,
   enumerate every length-k window of the mutant protein spanning the mutated residue (valid 0-based start
   s in [max(0,p-1-k+1), min(p-1,L-k)]; interior mutation -> exactly k windows/length), each paired with
   its matched wild-type (germline) peptide at the same coordinates + the mutation offset = the mutant/WT
   agretope (DAI, Wells 2020 TESLA). Class I lengths 8-11 (pVACtools Hundal 2020; ProGeo-neo Li 2020 21-mer
   +-10 flank; NetMHCpan 8-14/9-mer-dominant Jurtz 2017). Oracles: interior Y5C on MKTAYIAKQRSTVWLNDEFGH ->
   5 windows per k=8/9/10/11 -> 20 candidates (WT = mutant with C->Y at offset); terminal M1V k=9 -> 1
   truncated window. Two source-backed scoping assumptions (single-residue substitution only; binding
   affinity out of scope, caller-supplied via ONCO-MHC-001); no source contradictions. Cross-linked
   mhc-peptide-binding-prediction (resolved its two "not yet ingested" mentions + added a relates_to edge)
   and hla-nomenclature-and-allele-specific-loh. Wired into algorithm-validation-evidence hub (frontmatter
   source + summary link + anchor bullet, source_commit bumped to 643a974). Index updated (1 source + 1
   concept, MHC lines de-placeholdered).
   graph: +2 nodes, +4 typed edges (concept relates_to test-unit-registry / mhc-peptide-binding-prediction / hla-nomenclature-and-allele-specific-loh; mhc relates_to neoantigen-peptide-generation)

## [2026-07-10] ingest | docs/Evidence/ONCO-PHYLO-001-Evidence.md
   Twenty-sixth Oncology unit: Tumor Phylogeny Reconstruction — a clonal-evolution tree assembled from
   CCF clusters via the sum rule + lineage-precedence rule (LICHeE Popic 2015 / PICTograph Zheng 2022).
   Judged DISTINCT and wiki-worthy (a constraint-satisfaction / perfect-phylogeny tree BUILDER, the
   reconstruction step CCF clustering stops short of; no general NJ/UPGMA phylogenetics concept exists in
   the wiki to reuse) -> dedicated concept [[tumor-phylogeny-clonal-tree-reconstruction]]. Two per-sample
   CCF-ordering constraints: lineage precedence (ancestor CCF >= descendant CCF - e + presence pattern
   u.CCF=0 => v.CCF=0, LICHeE Eq. 2) and the sum rule (sum of children CCF <= parent CCF + e per node =
   pigeonhole generalization, LICHeE Eq. 5 / PICTograph sum condition -> forces branching vs nesting),
   plus trunk (CCF~1 root-path = clonal) vs branch (subclone) identification. Oracles: linear
   Normal->A->B->C (1.0/1.0/0.6/0.3, Trunk {A}); two-sample branching A->{B,C} (B 0.6/0.0, C 0.0/0.7,
   mutually non-ancestral); equal 0.6/0.6 siblings forced into a chain by the sum rule. Two
   reconstructed-tree invariants as oracles (ancestor >= descendant per edge; per-node sum rule). Two
   source-consistent assumptions (deepest-valid-ancestor deterministic tie-break for under-constrained
   private clusters; default e=0 strict inequalities); no source contradictions (LICHeE VAF form and
   PICTograph CCF form state the identical two constraints). Cross-linked cancer-cell-fraction-clonal-
   clustering (added a reciprocal "Downstream reconstruction" note + depends_on edge),
   clonal-subclonal-classification-ccf-posterior, intratumor-heterogeneity-metrics. Wired into
   algorithm-validation-evidence hub (frontmatter source + summary link + anchor bullet, source_commit
   bumped to ea992b8). Index updated (1 source + 1 concept).
   graph: +2 nodes, +3 typed edges (concept relates_to test-unit-registry; depends_on cancer-cell-fraction-clonal-clustering; relates_to clonal-subclonal-classification-ccf-posterior)

- 2026-07-10 — Ingested docs/Evidence/ONCO-PLOIDY-001-Evidence.md (Tumor Ploidy Estimation +
   Whole-Genome-Doubling detection; twenty-seventh Oncology unit). GENUINELY DISTINCT from the
   ONCO-ASCAT-001 joint grid fit: this is a post-hoc summary over already-called allele-specific
   segments, not an inference from raw logR/BAF. Two methods: (1) average ploidy `EstimatePloidy` =
   length-weighted mean per-segment total CN ψ = Σ(CN_i·L_i)/Σ(L_i) (Patchwork PMC4053982; n-scale
   Van Loo 2010, pure diploid → 2.0, >2.7n aneuploidy; CN 2/4/3 at 100/100/50 Mb → 3.0); (2)
   whole-genome doubling `DetectWholeGenomeDoubling` = facets-suite `is_genome_doubled` / Bielski 2018
   rule — WGD iff autosome-restricted fraction at major CN ≥ 2 strictly > 0.5 (mcn = tcn − lcn,
   denominator = reference chromosome-size table, ReferenceGenome {GRCh38,GRCh37}, UCSC hg38/hg19
   Ensembl-verified Σchr1–22 = 2,875,001,522 / 2,881,033,286 bp, autosomes only), legacy
   `DetectWholeGenomeDoublingFromSuppliedLength` keeps the supplied-length denominator. Created source
   onco-ploidy-001-evidence + new concept tumor-ploidy-estimation-and-whole-genome-doubling. Cross-linked
   allele-specific-copy-number-ascat (added post-hoc-summary note), loss-of-heterozygosity-detection,
   homologous-recombination-deficiency-score (shared segment substrate). Wired into
   algorithm-validation-evidence hub (frontmatter source + summary link + anchor bullet, source_commit
   bumped to 57c2be1). Index updated (1 source + 1 concept).
   graph: +2 nodes, +2 typed edges (concept relates_to test-unit-registry; relates_to allele-specific-copy-number-ascat)

## [2026-07-10] ingest | docs/Evidence/ONCO-PURITY-001-Evidence.md
   ONCO-PURITY-001 (Tumor Purity Estimation from somatic SNV VAF / allele-specific copy number),
   twenty-eighth Oncology unit. A genuinely distinct STANDALONE closed-form purity estimator, NOT an
   enrichment of ONCO-ASCAT-001: it inverts the CNAqc expected-VAF model v = m·π·c/[2(1−π)+π(n_A+n_B)]
   (Antonello 2024, verbatim) — EstimatePurityFromVAF = copy-neutral diploid het special case π = 2·VAF,
   EstimatePurity = general inversion π = 2v/(m·c + 2v − v·n_tot); FACETS 2016 confirms the denominator,
   ABSOLUTE 2012 is the inverse direction. Created source onco-purity-001-evidence + new concept
   tumor-purity-from-mutation-vaf. Cross-linked allele-specific-copy-number-ascat (alternative_to note:
   VAF inversion vs logR/BAF grid) and tumor-ploidy-estimation-and-whole-genome-doubling (reciprocal
   purity-side counterpart). Wired into algorithm-validation-evidence hub (frontmatter source + summary
   link + anchor bullet, source_commit bumped to fdf583e). Index updated (1 source + 1 concept).
   graph: +2 nodes, +2 typed edges (concept relates_to test-unit-registry; alternative_to allele-specific-copy-number-ascat)

## [2026-07-10] ingest | docs/Evidence/ONCO-SIG-001-Evidence.md
   ONCO-SIG-001 (SBS-96 single-base-substitution trinucleotide context catalog — pyrimidine-strand
   folding), twenty-ninth Oncology unit and the wiki's first mutational-signature method. Genuinely
   distinct — no existing concept covers mutational signatures. Created source onco-sig-001-evidence +
   new concept sbs96-mutational-signature-catalog: the 96-channel catalog (6 pyrimidine subtypes
   C>A/C>G/C>T/T>A/T>C/T>G × 4 5′ × 4 3′, labelled 5'[REF>ALT]3', mutated base centred; COSMIC SBS96 +
   SigProfilerMatrixGenerator Bergstrom 2019 + Alexandrov 2013, identical 6×4×4 definition) with the
   defining pyrimidine-strand folding rule (purine A/G reference reverse-complemented onto the pyrimidine
   strand — context via A↔T/C↔G + reverse, plus the substitution — before counting; C/T self-classifies).
   Seven worked folding oracles + partition invariants (exactly 96 keys, Σ counts = classifiable SBS
   variants). SCOPE NOTE: this unit is catalog/classification ONLY — the NMF/NNLS signature-exposure
   fitting against COSMIC reference signatures (mentioned in the ingest hint) is a separate downstream
   concern, not in this Evidence file. Cross-linked orthogonal ONCO biomarkers
   homologous-recombination-deficiency-score + microsatellite-instability-detection. Wired into
   algorithm-validation-evidence hub (frontmatter source + summary link + anchor bullet, source_commit
   bumped to 6fdbd84). Index updated (1 source + 1 concept). One cosmetic label-rendering assumption; no
   contradictions.
   graph: +2 nodes, +1 typed edge (concept relates_to test-unit-registry)

## [2026-07-10] ingest | docs/Evidence/ONCO-SIG-002-Evidence.md
   ONCO-SIG-002 (mutational signature fitting/refitting + de-novo NMF extraction), thirtieth Oncology
   unit and the downstream deconvolution partner of ONCO-SIG-001. Genuinely distinct from the SBS-96
   catalog — SIG-001's page explicitly deferred "NMF extraction / NNLS exposure estimation" as a separate
   concern, and no existing concept covers signature fitting. Created source onco-sig-002-evidence + new
   concept mutational-signature-fitting-and-extraction: supervised NNLS refit (min ‖Sx−d‖², x≥0;
   Lawson-Hanson active-set clamp-and-refit → S=[[1,1],[0,1]],d=[0,1]⇒[0,0.5]; cosine-≥0.95 reconstruction
   gate; raw + proportion exposures) and unsupervised de-novo NMF (Lee & Seung Frobenius+KL multiplicative
   updates, monotone non-increase, V≈WH blind-source-separation, COSMIC L1-normalized signatures), shared
   cosine metric (zero-norm→0.0 convention), Brunet-2004 cophenetic + SigProfiler silhouette rank
   selection, greedy best-cosine COSMIC reference matching; NMF non-convex→local-minimum / permutation-
   scale ambiguity / ε-guarded denominators sharp edges. Cross-linked the catalog concept
   sbs96-mutational-signature-catalog both ways (its two "separate downstream concern" mentions now point
   here; new concept depends_on it). Wired into algorithm-validation-evidence hub (frontmatter source +
   summary link + anchor bullet, source_commit bumped to 8cb9903). Index updated (1 source + 1 concept).
   Five modelling assumptions (Frobenius objective, seeded init, exposure proportions, zero-vector cosine,
   consensus/silhouette/greedy matching); no source contradictions.
   graph: +2 nodes, +2 typed edges (concept relates_to test-unit-registry; concept depends_on sbs96-mutational-signature-catalog)

## [2026-07-10] ingest | docs/Evidence/ONCO-SIG-003-Evidence.md
   ONCO-SIG-003 (signature exposure bootstrap confidence intervals), thirty-first Oncology unit and the
   uncertainty layer directly above the ONCO-SIG-002 NNLS refit. Genuinely distinct from the fitting page:
   it adds no decomposition but wraps FitSignatures in a resample→refit→percentile loop producing a
   per-signature CI. Created source onco-sig-003-evidence + new concept
   signature-exposure-bootstrap-confidence-intervals: resample the 96-channel catalog R times (default
   1000, sigminer ≥100), re-run NNLS per replicate, take [2.5%,97.5%] percentiles (Efron 1979) via the
   type-7 sample quantile (Hyndman & Fan 1996, R/NumPy default). Two resampling schemes differing only by
   whether total burden N is fixed — multinomial (sigminer fixed-N, the byte-for-byte default) vs Poisson
   (Senkin 2021 MSA variant, each channel Poisson(observedₖ), N unfixed, Poisson↔multinomial conditional
   equivalence). Discriminating corner case = single non-zero channel: multinomial collapses
   deterministically (width 0) while Poisson(λ>0) fluctuates (var=mean, positive width — the reason the
   Poisson variant was added); plus N=0→[0,0], R=1→lower=upper=mean. Type-7 oracles [0,1,2,3,4]→0.1/2.0/3.9
   and [2,4,6,8]→2.15/5.0/7.85. Cross-linked mutational-signature-fitting-and-extraction both ways (new
   concept depends_on it; its NNLS section now points here for uncertainty). Wired into
   algorithm-validation-evidence hub (frontmatter source + summary link + anchor bullet, source_commit
   bumped to 2c404cc). Index updated (1 source + 1 concept). Two source-aligned assumptions (type-7
   interpolation, fixed seed 42), multinomial the backward-compatible default; no source contradictions.
   graph: +2 nodes, +2 typed edges (concept relates_to test-unit-registry; concept depends_on mutational-signature-fitting-and-extraction)

## [2026-07-10] ingest | docs/Evidence/ONCO-SIG-004-Evidence.md
   ONCO-SIG-004 (mutational process classification), thirty-second Oncology unit and the aetiology-annotation
   layer over the ONCO-SIG-002 NNLS refit. Genuinely distinct from the fitting/bootstrap siblings: adds no
   decomposition but turns per-signature exposures into a set of active mutational processes via
   normalize → cutoff → map → aggregate. Created source onco-sig-004-evidence + new concept
   mutational-process-classification: normalize to relative contributions Wᵢ = exposureᵢ/Σexposure
   (deconstructSigs "weights between 0 and 1"); drop Wᵢ < 0.06 (verbatim signature.cutoff = 0.06 /
   weights[weights < signature.cutoff] <- 0, strict < so exactly 0.06 retained, 1.4% false-negative
   calibration); map each surviving COSMIC label to its proposed aetiology (SBS1/5→Aging, SBS2/13→APOBEC,
   SBS4→Tobacco, SBS7a–d→UV, SBS6/15/20/26→MMR; COSMIC + Alexandrov 2020); sum member contributions per
   process → active-process set + argmax dominant process. Hand-derived oracle SBS2/13/1/4 = 50/30/15/5 →
   APOBEC 0.80 (dominant) / Aging 0.15 / Tobacco 0 (SBS4 0.05 below cutoff). Corner cases: surviving mass < 1
   (rest unknown), multiple simultaneous processes, unmapped/"Unknown"-aetiology label → no process, Σ=0 →
   none. Two source-aligned assumptions (per-process summation, per-signature-cutoff-then-group); no source
   contradictions. Cross-linked mutational-signature-fitting-and-extraction (its NNLS section now points here
   for interpretation; new concept depends_on it). Wired into algorithm-validation-evidence hub (frontmatter
   source + summary link + anchor bullet, source_commit bumped to 3c5a975). Index updated (1 source + 1 concept).
   graph: +2 nodes, +2 typed edges (concept relates_to test-unit-registry; concept depends_on mutational-signature-fitting-and-extraction)

## [2026-07-10] ingest | docs/Evidence/ONCO-SOMATIC-001-Evidence.md
   ONCO-SOMATIC-001 (somatic mutation calling, tumor vs matched-normal VAF classification), thirty-third
   Oncology unit and the foundational somatic SNV caller at the head of the pipeline. Genuinely distinct
   from the downstream QC filters and copy-number/clonal units: classify each candidate variant
   Somatic/Germline/NotDetected by comparing tumor allele frequency f_t vs matched-normal f_n
   (VAF=altReads/totalReads). Created source onco-somatic-001-evidence + new concept
   somatic-variant-calling-tumor-normal. Strelka somatic state S={(f_t,f_n): f_t≠f_n} restricted to a
   ref/ref normal genotype (Saunders 2012; raw somatic prob over-calls in LOH/CN regions) / Strelka2
   continuous-VAF somatic-LOD+VAF (Kim 2018); three configurable thresholds (tumor LoD f_t≥0.05 Yan 2021
   WES 5%, normal absent ceiling f_n≤0.01 normalVafThreshold, f_n>0.01→Germline per Mutect2 germline
   filter Benjamin 2019); bounded monotone somatic score max(0,f_t−f_n)∈[0,1]; tumor-only mode
   (no matched normal → Mutect2 ℓ_n=1); FilterGermlineVariants = somatic subset of input. Oracles
   A 0.25/0.00→Somatic · B 0.48/0.50→Germline · C 0.02→NotDetected · D 0.30/0.03→Germline (CHIP-like) ·
   E tumor-only 0.20→Somatic · boundaries f_t=0.05 present / f_n=0.01 absent. Two flagged source-consistent
   assumptions (1% normal ceiling parameterized not invented; score a documented simplification not a
   caller LOD); no source contradictions. Wired into algorithm-validation-evidence hub (frontmatter source
   + summary link + anchor bullet, source_commit bumped to cd2346b7). Cross-linked upstream of the two QC
   filters sequencing-artifact-detection + clonal-hematopoiesis-cfdna-filtering (reciprocal body links
   added to both). Index updated (1 source + 1 concept).
   graph: +2 nodes, +3 typed edges (concept relates_to test-unit-registry; concept relates_to sequencing-artifact-detection; concept relates_to clonal-hematopoiesis-cfdna-filtering)

## [2026-07-10] ingest | docs/Evidence/ONCO-SV-001-Evidence.md
   ONCO-SV-001 (somatic complex-rearrangement classification / chromothripsis inference), thirty-fourth
   Oncology unit and the wiki's first complex-SV / chromothripsis method. Genuinely distinct region-level
   pattern classifier (Chromothripsis vs NotComplex) over a per-segment integer CN profile — created new
   concept chromothripsis-inference + source onco-sv-001-evidence. Korbel & Campbell 2013 six hallmark
   criteria (A clustering / B oscillating CN states / C-F heterozygosity, haplotype, randomness, derivative
   walk); computes B (oscillation = adjacent-segment CN-state reversal, ≤3 canonically 2 states) gated by
   Magrangeas-2011 ≥10 first-pass oscillation screen + Cortes-Ciriano-2020 tiers (≥7 high / 4-6 low / <4
   not-called), >60% two-state fraction, ≥6 clustered intrachromosomal SV floor; and A (breakpoint
   clustering) via exponential-null CV>1 flag. Oracles 2,1,...×11→10→Chromothripsis · 6-seg→5<10→NotComplex
   · monotone 2..7→0 oscillations/>2 states→NotComplex (progressive amp/BFB, clustering necessary-but-not-
   sufficient). Two operationalisation assumptions (oscillation=CN-state-transition count; clustering=CV>1
   vs exponential CV=1); no source contradictions. Consumes the per-segment CN states of
   copy-number-alteration-classification (ONCO-CNA-001, reciprocal body link added there); orthogonal to
   gene-fusion-detection-read-evidence + focal-amplification-detection, distinct from the gene-order
   genome-rearrangement-breakpoint-distance. Wired into algorithm-validation-evidence hub (frontmatter
   source + summary link + anchor bullet, source_commit bumped to 1d2674a9). Index updated (1 source + 1
   concept).
   graph: +2 nodes, +2 typed edges (concept relates_to test-unit-registry; concept relates_to copy-number-alteration-classification)

## [2026-07-10] ingest | docs/Evidence/ONCO-TMB-001-Evidence.md
   ONCO-TMB-001 (Tumor Mutational Burden — mutations/Mb + TMB-high classification), thirty-fifth Oncology
   unit. Genuinely distinct immunotherapy biomarker (not represented) — created new concept
   tumor-mutational-burden + source onco-tmb-001-evidence. TMB = counted somatic mutations / sequenced
   coding region in Mb (Chalmers 2017 Methods; denominator = assay coding footprint, FoundationOne
   315-gene 1.1 Mb / F1CDx 324-gene ~0.8 Mb / WES ~30-40 Mb, taken as a parameter; panel counting includes
   synonymous to reduce noise, germline/driver removed before counting). TMB-High = TMB >= 10 mut/Mb
   inclusive (FDA pembrolizumab tumor-agnostic approval 2020-06-16, F1CDx companion diagnostic, Marcus 2021;
   mut/Mb reporting + cutoff cross-confirmed by FoCR Harmonization Project, Merino 2020). Oracles 11/1.1->10.0
   High · 300/30->10.0 High · 150/10->15.0 High · 99/10->9.9 not-High · 100/10->10.0 High(boundary) · 0->not-
   High; corner cases regionMb=0->div-by-zero throws · negative rejected · <0.5 Mb computes but known-unstable
   · monotone in count/region. One flagged conflict: unsupported registry three-tier Low<6/Intermediate
   6-20/High>20 (no source for the 6/20 boundaries) resolved to the single source-backed two-tier >=10 cutoff.
   Counts the caller-supplied somatic list of somatic-variant-calling-tumor-normal (ONCO-SOMATIC-001), sibling
   immunotherapy biomarker of microsatellite-instability-detection (reciprocal body link added there),
   correlated with neoantigen-peptide-generation. Wired into algorithm-validation-evidence hub (frontmatter
   source + summary link + anchor bullet, source_commit bumped to 701e1721). Index updated (1 source + 1
   concept).
   graph: +2 nodes, +3 typed edges (concept relates_to test-unit-registry; concept relates_to somatic-variant-calling-tumor-normal; concept relates_to microsatellite-instability-detection)

## [2026-07-10] ingest | docs/Evidence/ONCO-VAF-001-Evidence.md
   ONCO-VAF-001 (Variant Allele Frequency Analysis — empirical VAF + Wilson binomial CI + purity/ploidy
   correction), thirty-sixth Oncology unit. Genuinely distinct (created new concept
   variant-allele-frequency-and-binomial-ci + source onco-vaf-001-evidence) rather than folding into the
   three neighbouring VAF units — because it owns the model-free VAF primitive and, above all, the Wilson
   score binomial confidence interval, which is not represented anywhere else in the wiki. Three quantities:
   (1) empirical VAF = altReads/totalReads (= alt AD / sum AD, GATK Mutect2 FAQ; deliberately NOT Mutect2's
   model-estimate AF FORMAT field); (2) Wilson score interval center=(p̂+z²/2n)/(1+z²/n), margin=
   (z/(1+z²/n))·√(p̂(1−p̂)/n+z²/4n²), z=1.96 for 95% (Wilson 1927 via Wikipedia), chosen over Wald for
   staying in [0,1] with non-zero width at the extremes (p̂=0→lower 0, p̂=1→upper 1); (3) AdjustVAFForPurity
   = m·CCF = VAF·(2(1−π)+π·n_tot)/π, inverting the CNAqc (Genome Biology 2024) expected-VAF model, normal
   ploidy 2. Oracles empirical 25/100→0.25 · 0/10→0.00 · 10/10→1.00; Wilson (25,100)→0.2592 [0.1755,0.3430]
   & no-overshoot (0,10)→[0,0.2775] / (10,10)→[0.7225,1]; correction (0.40,0.80,2)→1.00 / (0.20,0.50,2)→0.80
   / (0.30,0.50,4)→1.80. Corner cases VAF>1 alignment-artifact→invalid, totalReads=0→0/0 guarded, π=0→
   undefined. Cross-linked heavily: empirical VAF = the altReads/totalReads primitive somatic-variant-calling-
   tumor-normal (ONCO-SOMATIC-001) compares and ctdna-detection-and-tumor-fraction averages; AdjustVAFForPurity
   shares the CNAqc/Tarabichi model tumor-purity-from-mutation-vaf inverts for π and cancer-cell-fraction-
   clonal-clustering divides by m for CCF (its m·CCF output = the CCF-formula numerator) — reciprocal body
   links added to all three. Two source-backed assumptions (z=1.96 verbatim not 1.959964; AdjustVAFForPurity
   normal CN=2 diploid background). Contradictions: none — GATK / Wilson 1927 / CNAqc+Tarabichi cover disjoint
   facets. Wired into algorithm-validation-evidence hub (frontmatter source + summary link + anchor bullet,
   source_commit bumped to 68661290). Index updated (1 source + 1 concept). Follow-ups: none.
   graph: +2 nodes, +4 typed edges (concept relates_to test-unit-registry; concept relates_to somatic-variant-calling-tumor-normal; concept relates_to tumor-purity-from-mutation-vaf; concept relates_to cancer-cell-fraction-clonal-clustering)

## [2026-07-10] ingest | docs/Evidence/PANGEN-CLUSTER-001-Evidence.md → pan-genome gene clustering
Ingested PANGEN-CLUSTER-001 (Gene Clustering — greedy incremental homolog grouping by sequence
identity), the first pan-genome PANGEN-* unit. New concept [[pan-genome-gene-clustering]] (greedy
incremental CD-HIT-model clustering: long→short, first-match representative assignment, global
identity = identical aligned positions / alignment length; Li & Godzik 2006 + CD-HIT guide + Roary
Page 2015 95% default + EMBOSS needle). New source [[pangen-cluster-001-evidence]]. Wired into the
[[algorithm-validation-evidence]] hub (frontmatter source + source-link list + anchor bullet) and
cross-linked to the comparative-genomics siblings [[genome-comparison-core-dispensable]],
[[ortholog-detection-reciprocal-best-hits]], [[average-nucleotide-identity]]. No source
contradictions; two assumptions (alignment-free identity may underestimate with internal indels;
homolog groups only, no synteny/paralog-split). index updated (1 source + 1 concept).
NOTE: original ingest subagent died mid-run on an API error after writing both pages + the hub
frontmatter/source-link edit; the hub anchor bullet, index entries, this log line, and the graph
rebuild were completed by the orchestrator during recovery.
   graph: +2 nodes, +2 typed edges

## [2026-07-10] ingest | docs/Evidence/PANGEN-CORE-001-Evidence.md
Ingested PANGEN-CORE-001 (pan-genome partition — core/accessory/unique by occupancy + genomic
fluidity + Heaps open/closed; `ConstructPanGenome`). Genuinely distinct from siblings: it is the
N-genome occupancy-based partition (fractional Roary 99% core rule, Kislyuk fluidity, Heaps alpha),
not the clustering step [[pan-genome-gene-clustering]] (PANGEN-CLUSTER-001) nor the pairwise-RBH
[[genome-comparison-core-dispensable]] (COMPGEN-COMPARE-001). Created source [[pangen-core-001-evidence]]
+ concept [[pan-genome-core-accessory-partition]]; cross-linked both siblings; hub
[[algorithm-validation-evidence]] frontmatter/source-list/anchor + index updated. No contradictions;
two source-backed assumptions (clustering delegated to k-mer-Jaccard ClusterGenes; empty-pair fluidity term→0).
   graph: +2 nodes, +3 typed edges

## [2026-07-10] ingest | docs/Evidence/PANGEN-HEAP-001-Evidence.md
Ingested PANGEN-HEAP-001 (Pan-Genome Growth Model — Heaps'-law fit of the new-gene curve). Genuinely
distinct from the sibling PANGEN-CORE-001: it is the dedicated fitting engine (presence/absence
binarization, micropan first-appearance new-gene curve, bounded power-law least-squares
`y=K·x^(-alpha)`, permutation pooling, open ⟺ alpha<1 rule) that the occupancy partition
[[pan-genome-core-accessory-partition]] only *reports* as one open/closed output — so a dedicated
concept was warranted rather than enriching the partition. Created concept [[pan-genome-heaps-law-fit]]
+ source [[pangen-heap-001-evidence]] (micropan `heaps()` powerlaw.R + Tettelin 2008 power-law openness
+ Tettelin 2005 *S. agalactiae* anchor). Cross-linked the partition sibling (added a fitting-engine
pointer to its Heaps section) and the clustering sibling [[pan-genome-gene-clustering]]. Hub
[[algorithm-validation-evidence]] frontmatter/source-list/anchor + index (1 source + 1 concept)
updated. Exact oracles x=[2,3] y=[8,4]→alpha≈1.7095/K≈26.164/closed and constant→alpha 0/K 1/open. No
contradictions; two source-backed assumptions (optimizer method non-correctness-affecting vs L-BFGS-B,
fixed-seed permutation RNG).
   graph: +2 nodes, +3 typed edges

## [2026-07-10] ingest | docs/Evidence/PANGEN-MARKER-001-Evidence.md
Fourth PANGEN-* unit (phylogenetic marker selection). Created source summary
[[pangen-marker-001-evidence]] and a NEW concept [[phylogenetic-marker-selection]] (genuinely distinct:
single-copy core marker selection + parsimony-informative-site scoring, not covered by clustering/
partition/Heaps siblings). SelectPhylogeneticMarkers keeps single-copy core clusters (panX "all strains
exactly once" + Roary 99% core + paralog filtering) with PIS≥1, ranked by descending PIS capped at
maxMarkers; CountParsimonyInformativeSites per Zvelebil & Baum 2008 (≥2 states each in ≥2 seqs).
Cross-linked from [[pan-genome-core-accessory-partition]] and [[pan-genome-gene-clustering]]; hub
[[algorithm-validation-evidence]] frontmatter/source-list/anchor + index updated. Oracles AAAAA/AAACA/
AACCG/ACCTG→PIS 2 (cols 3,5), 3-genome selection excludes paralog/not-core/0-PIS-conserved. One
source-backed assumption (no in-repo aligner → PIS over equal-length members, unequal→PIS 0), no
contradictions.
   graph: +2 nodes, +4 typed edges

## 2026-07-10 — ingest PARSE-BED-001-Evidence
Ingested docs/Evidence/PARSE-BED-001-Evidence.md (first FileIO/PARSE-* file-parsing unit). Source
summary [[parse-bed-001-evidence]]; new concept [[bed-format-parsing]] as the PARSE-* family anchor —
UCSC BED interval parsing on the 0-based half-open coordinate model (chromStart 0-based / chromEnd
1-based-exclusive, chromStart==chromEnd = zero-length insertion), BED3→BED12 column ladder, and UCSC
validation rules (chromStart ≤ chromEnd else null, score clamp [0,1000], strand +/−/., first-line
column-count lock, BED12 block constraints). Sources UCSC FAQ + Wikipedia + BEDTools (Quinlan & Hall
2010); deviations None. Hub [[algorithm-validation-evidence]] frontmatter/source-list/anchor + index
(Sources + Concepts) updated; cross-linked to [[fuzzing]] (parsers = hottest malformed-input target).
   graph: +2 nodes, +1 typed edge

## [2026-07-10] ingest | docs/Evidence/PARSE-EMBL-001-Evidence.md → parse-embl-001-evidence (source)
   EMBL flat-file parsing (EmblParser.Parse/ParseFile): line-type records (ID/AC/DT/DE/KW/OS/OC/
   FH/FT/SQ/`//`), ID-line grammar, INSDC feature-table location descriptors (simple/^-site/</>
   partials, complement/join/order operators, remote refs accession.version:span, no nested join/
   order), data-class/division/IUPAC vocabularies, lowercase-sequence normalization. Enhancement:
   offline-first caller-supplied resolver for remote-aware location→sequence assembly (complement
   of a join reverses order; remote spans 1-based inclusive); remote-prefix per-segment strip fix
   (Location.RemoteParts); </>-verbatim-slice + missing-resolver→empty-segment assumptions.
   Sources EBI EMBL User Manual Rel. 143 + INSDC Feature Table v11.3; deviations None. No new
   concept — cross-linked to family anchor [[bed-format-parsing]] (INSDC grammar shared with the
   GenBank cousin, not yet ingested). Hub [[algorithm-validation-evidence]] frontmatter/source-list/
   per-file-link + index Sources updated; cross-linked to [[test-unit-registry]] and [[fuzzing]].
   graph: +1 node, +0 typed edges

## [2026-07-10] ingest | PARSE-FASTA-001-Evidence.md → parse-fasta-001-evidence (source)
   FileIO/PARSE-* family FASTA parsing (FastaParser Parse/ParseFile/ParseFileAsync + ToFasta/
   WriteFile): >defline (first-word Id / rest description) + sequence lines, multi-FASTA, opt-in
   SequenceAlphabet (default strict DNA; IUPAC-nucleotide/RNA/protein), round-trip w/ line-width 80.
   Sources Wikipedia FASTA + NCBI BLAST/FASTA spec + Lipman&Pearson 1985/1988 + NC-IUB 1985/IUPAC
   tables; deviations None (default-DNA + header-without-sequence-not-yielded + blank-line-skip
   assumptions). No new concept — FASTA has no coordinate/record grammar to summarize; cross-linked
   to family anchor [[bed-format-parsing]] and sibling [[parse-embl-001-evidence]]. Hub
   [[algorithm-validation-evidence]] frontmatter/source-list/per-file-link + index Sources updated;
   cross-linked to [[test-unit-registry]] and [[fuzzing]].
   graph: +1 node, +0 typed edges

## [2026-07-10] ingest | docs/Evidence/PARSE-FASTQ-001-Evidence.md → parse-fastq-001-evidence (source) + 1 concept
   FileIO PARSE-* family, fourth parsing unit (after BED/EMBL/FASTA). Source page for the FASTQ
   parsing artifact (4-line record @header/seq/+/quality, seq-len==qual-len invariant, Q20/Q30 +
   per-position stats, quality/length filter, quality+adapter trim, round-trip; Wikipedia FASTQ +
   Cock et al. 2009 + NCBI SRA sources; edge cases multi-line/@-in-quality/+-in-seq/blank-skip;
   assumptions Q93-cap-for-p<=0 + ambiguous-window->Phred+33 default; no contradictions). Unlike the
   grammar-only FASTA/EMBL siblings, FASTQ carries a genuine encoding scheme, so created the
   cross-cutting concept phred-quality-encoding (Q=-10log10p, Phred+33 vs Phred+64 ASCII offsets +
   ranges + boundary chars, deterministic per-record offset auto-detection, mis-detection = silent
   corruption) — not yet represented and shared with the Assembly trimming layer. Cross-linked the
   new concept from quality-trimming-running-sum's Phred+33 section (2 inbound links). Linked source
   into the algorithm-validation-evidence hub (frontmatter source-list + per-file-link) and both new
   pages to family anchor [[bed-format-parsing]], siblings, [[test-unit-registry]], [[fuzzing]]. Index
   Sources + Concepts updated.
   graph: +2 nodes, +1 typed edge (phred-quality-encoding relates_to test-unit-registry)

## [2026-07-10] ingest | PARSE-GENBANK-001-Evidence.md → parse-genbank-001-evidence (source) + insdc-feature-location (new concept)
   Ingested the GenBank flat-file-parsing Evidence artifact (FileIO/PARSE-* family). Wrote source
   summary parse-genbank-001-evidence (NCBI Sample Record + Wikipedia + INSDC feature-table sources;
   LOCUS/section grammar, 18 divisions, ORIGIN lowercase-normalized sequence, U49845 record,
   defensive null/empty/missing-LOCUS/missing-ORIGIN contracts, length-match/`//`/Start≤End
   invariants). Created shared concept insdc-feature-location (the DDBJ/ENA/GenBank location-descriptor
   grammar — complement/join/order/partial/remote, 1-based inclusive, operator assembly semantics +
   caller-supplied offline resolver + <>-verbatim-slice assumption; oracles join(1..3,7..9)→ACGGTA,
   complement(Y.1:1..4)→GTTT), now warranted with 2 inbound units (GenBank + EMBL). Re-pointed the
   EMBL source page's "no separate concept yet" note to the new concept. Linked GenBank source into the
   algorithm-validation-evidence hub (frontmatter + body list) and to [[bed-format-parsing]], siblings,
   [[test-unit-registry]], [[fuzzing]]. Index Sources + Concepts updated. No contradictions.
   graph: +2 nodes, +1 typed edge (insdc-feature-location relates_to test-unit-registry)

## [2026-07-10] ingest | PARSE-GFF-001-Evidence.md → parse-gff-001-evidence (source only, cross-link)
   Ingested the GFF/GTF annotation-file-parsing Evidence artifact (FileIO/PARSE-* family). Wrote
   source summary parse-gff-001-evidence (Wikipedia General Feature Format + UCSC GFF/GTF FAQ +
   Sequence Ontology GFF3 v1.26 sources; 9 tab-delimited columns, 1-based inclusive coords, phase
   0/1/2 for CDS, attribute dialects GFF3 key=value; vs GTF key "value"; vs GFF2 group, Parent
   part-of hierarchy + multi-parent + discontinuous features, RFC 3986 percent-escaping, directives,
   null/empty→empty + <8-fields→skip contracts, column/attribute/escape/hierarchy/format-detect/
   round-trip test categories; deviations None). Determination: NO new concept — GFF is a
   tab-delimited sibling of BED, cross-linked to family anchor [[bed-format-parsing]] where the
   BED-vs-GFF 0-based-vs-1-based coordinate contrast already lives; distinct GFF facts captured in
   the source page (economical per directive). Added inbound link from bed-format-parsing anchor;
   registered the source in the algorithm-validation-evidence hub (frontmatter + body list). Index
   Sources updated. No contradictions.
   graph: +1 node (parse-gff-001-evidence source), +0 typed edges (source-only page; mentions edges auto-derived)

## [2026-07-10] ingest | docs/Evidence/PARSE-VCF-001-Evidence.md → [[parse-vcf-001-evidence]]
   VCF (Variant Call Format) parsing, seventh FileIO/PARSE-* unit. Source-only cross-link (no new
   concept): VCF is a tab-delimited, 1-based sibling of BED/GFF, cross-linked to family anchor
   [[bed-format-parsing]] where the 0-based-vs-1-based coordinate contrast already lives (anchor
   names GFF/GTF/VCF as its 1-based counterparts). VCF-specific richness (##headers, 8 fixed cols
   + FORMAT/genotype samples, SNP/MNP/Ins/Del/Symbolic classification, the audit's five
   spec-compliance points incl. FILTER "." ≠ PASS and Ti/Tv-over-all-ALT) captured in the source
   page (economical per directive). Added inbound link from bed-format-parsing anchor; registered
   the source in the algorithm-validation-evidence hub (frontmatter + body list). Index Sources
   updated. No contradictions.
   graph: +1 node (parse-vcf-001-evidence source), +0 typed edges (source-only page; mentions edges auto-derived)

## [2026-07-10] ingest | docs/Evidence/PAT-APPROX-003-Evidence.md → [[pat-approx-003-evidence]] (source) + [[approximate-pattern-matching-mismatches]] (new concept)
   Approximate (Hamming-distance / k-mismatch) pattern matching — ROSALIND BA1H/BA1I/BA1N
   (Compeau & Pevzner) + go-rosalind/Rosalind-Solutions reference impls. New concept: the PAT-APPROX
   family anchor (Count_d / FindApproximateOccurrences / FindFrequentKmersWithMismatches / Neighbors
   d-neighborhood / FindBestMatch). Genuinely distinct from the exact matchers — created rather than
   folded in. Cross-linked as mismatch-tolerant alternative_to [[k-mer-positions]] and fuzzy sibling of
   [[known-motif-search]] / [[overrepresented-kmer-discovery]]; scoped against the indel-tolerant
   alignment/edit-distance family. Registered in algorithm-validation-evidence hub (frontmatter +
   source list + own-concept list); index Sources + Concepts updated. Deviations none; one FindBestMatch
   leftmost-tie-break API assumption.
   graph: +2 nodes (pat-approx-003-evidence source, approximate-pattern-matching-mismatches concept), +2 typed edges (relates_to test-unit-registry, alternative_to k-mer-positions)

## [2026-07-10] ingest | docs/Evidence/PHYLO-BOOT-001-Evidence.md → [[phylo-boot-001-evidence]] (source) + [[phylogenetic-bootstrap-support]] (new concept)
   First phylogenetics (PHYLO-*) unit. Survey confirmed NO general distance-based-phylogenetics
   concept exists yet (no Neighbor-Joining/UPGMA/distance-matrix page), so created a dedicated anchor
   concept for the family: Felsenstein's bootstrap (FBP) — resample alignment columns with replacement
   → same-length pseudo-alignments → tree per replicate (distance matrix → UPGMA/NJ) → clade support =
   fraction of replicate trees containing the identical terminal-leaf-set clade of the reference tree.
   Sources Felsenstein 1985 (Evolution) + Lemoine 2018 (Nature, PMC6030568) + Biopython
   Bio.Phylo.Consensus (get_support terminal-set matching). Oracles: two-group AAAAAAAAAA/GGGGGGGGGG
   UPGMA+JukesCantor seed 42 → support 1.0 (saturated distances), all-identical ACGTACGT → all 1.0.
   Two source-consistent assumptions (rooted-clade vs unrooted bipartition matching Biopython;
   proportion [0,1] vs percentage ×100). Explicitly distinguished from the CCF-constraint oncology
   builder [[tumor-phylogeny-clonal-tree-reconstruction]] (no distance matrix / no resampling) and
   positioned as the confidence layer over the tree [[phylogenetic-marker-selection]] feeds; added a
   reciprocal distinction link on the tumor-phylogeny page. Registered in algorithm-validation-evidence
   hub (frontmatter + source list + own-concept list); index Sources + Concepts updated. No contradictions.
   Follow-up: when NJ/UPGMA/distance-matrix or other PHYLO-* units are ingested, split tree-construction
   into its own concept and make bootstrap depend_on it.
   graph: +2 nodes (phylo-boot-001-evidence source, phylogenetic-bootstrap-support concept), +1 typed edge (relates_to test-unit-registry)

## [2026-07-10] ingest | docs/Evidence/PHYLO-COMP-001-Evidence.md → [[phylo-comp-001-evidence]] (source) + [[tree-comparison-metrics]] (new concept)
   Second phylogenetics (PHYLO-*) unit. Judged genuinely distinct from the PHYLO-BOOT-001 bootstrap
   anchor and warranted a dedicated concept: PHYLO-COMP-001 is deterministic tree COMPARISON / query
   over an already-built rooted binary PhyloNode tree, not confidence-by-resampling. Three operations —
   Robinson–Foulds distance (RobinsonFouldsDistance = symmetric difference A+B of canonical split sets;
   raw count, proper metric, even; rooted-clade = Wikipedia dummy-leaf; max RF = 2(n−2) rooted vs
   2(n−3) unrooted, reconciled by dummy-leaf equivalence), MRCA (FindMRCA, O(n) recursive leaf-name
   traversal, not-in-tree→null), patristic distance (PatristicDistance = dist(x,MRCA)+dist(y,MRCA),
   not-in-tree→NaN, zero-branches→0). Sources Wikipedia RF-metric/MRCA/Phylogenetic-tree + Robinson &
   Foulds 1981 (doi:10.1016/0025-5564(81)90043-2) + Smith 2020 (btaa614) + Day 1985. Two source-backed
   scope decisions (binary-only via PhyloNode Left/Right; rooted-only via UPGMA/NJ); no deviations, no
   contradictions. Cross-linked reciprocally with [[phylogenetic-bootstrap-support]] (agreement vs
   disagreement of the same split primitive) and distinguished from [[tumor-phylogeny-clonal-tree-reconstruction]].
   Registered in algorithm-validation-evidence hub (frontmatter source list + source_commit bump +
   source-summary list + concept list); index Sources + Concepts updated.
   graph: +2 nodes (phylo-comp-001-evidence source, tree-comparison-metrics concept), +1 typed edge (relates_to test-unit-registry)

## [2026-07-10] ingest | docs/Evidence/PHYLO-DIST-001-Evidence.md → phylo-dist-001-evidence (source) + 1 concept
   Third phylogenetics-family (PHYLO-*) Evidence file (after BOOT, COMP). Created the genuinely-distinct
   concept evolutionary-distance-matrix — the pairwise evolutionary-distance substrate the family sits on:
   CalculatePairwiseDistance / CalculateDistanceMatrix over four methods Hamming (raw diff count),
   p-distance (differences/comparableSites, uncorrected), Jukes-Cantor JC69 (-3/4·ln(1−4p/3), equal-base/
   equal-rate model), Kimura-2-parameter K2P (-1/2·ln((1−2S−V)·√(1−2V)), S=transition/V=transversion).
   Matrix invariants symmetric / zero-diagonal / non-negative / n×n, correction ordering JC69≥p & K2P≥p,
   triangle inequality expected-not-guaranteed for corrected distances; pairwise deletion (gaps + ambiguous
   IUPAC N/R/Y skipped, only A,C,G,T compared), case-insensitive; saturation JC69→+∞ at p≥3/4 and K2P→+∞ at
   V≥1/2; all-gap/empty→0 (0/n→0), unequal-length→ArgumentException, null→ArgumentNullException. Oracles
   ACGTACGT/TCGTACGT→Hamming 1/p 0.125/JC69≈0.137, pure-transition GCGT K2P≈0.34657 vs pure-transversion
   CCGT≈0.31713, mixed≈0.30679, gap case→7 comparable sites. Positioned as the UPGMA/NJ substrate that
   bootstrap wraps and re-runs per replicate and that tree-comparison's trees come out of. Sources: Wikipedia
   Models-of-DNA-evolution / Substitution-model / Distance-matrices-in-phylogeny + Jukes & Cantor 1969 +
   Kimura 1980 + Felsenstein 2004. Concise source page for the artifact. Registered in algorithm-validation-
   evidence hub (frontmatter source list + source_commit bump to 3a53115 + source-summary list + concept
   list); index Sources + Concepts updated. Updated the two prior PHYLO concepts reciprocally: fixed the now-
   stale "no distance-matrix page yet" note in phylogenetic-bootstrap-support and cross-linked evolutionary-
   distance-matrix from both bootstrap and tree-comparison-metrics. Two API-contract assumptions (empty/all-
   gap→0, pairwise deletion for gaps+ambiguity). Contradictions: none — JC69/K2P formulas, symmetric zero-
   diagonal matrix, and saturation limits are the standard textbook definitions. Follow-up: the UPGMA/NJ tree
   *construction* step itself (which consumes this matrix) remains unconcepted and would warrant its own page
   when a PHYLO tree-build unit is ingested.
   graph: +2 nodes (phylo-dist-001-evidence source, evolutionary-distance-matrix concept), +2 typed edges (relates_to test-unit-registry, relates_to phylogenetic-bootstrap-support)

## [2026-07-10] ingest | docs/Evidence/PHYLO-NEWICK-001-Evidence.md → [[phylo-newick-001-evidence]]
   Fourth phylogenetics PHYLO-* Evidence file: Newick I/O (ToNewick/ParseNewick), the tree
   serialization layer. Decision: source-only cross-link, NO new concept — Newick is a format
   serializer, not a distinct algorithm; the tree semantics it round-trips already live in the PHYLO
   concept pages. Source page synthesizes the Wikipedia/PHYLIP-Felsenstein/Olsen-1990 grammar,
   label rules, invariants N1–N9, and binary-only/no-quoted-labels/no-comments scope limits. Linked
   from the [[algorithm-validation-evidence]] hub (added to frontmatter sources list + body list) and
   index Sources; cross-linked reciprocally from [[tree-comparison-metrics]] (the PhyloNode trees it
   compares are what Newick serializes) and to [[evolutionary-distance-matrix]] /
   [[phylogenetic-bootstrap-support]]. No typed concept-to-concept edges (source-only). Contradictions:
   none — grammar and label rules are the standard Newick spec. Follow-ups: none.
   graph: +1 node (phylo-newick-001-evidence source), +0 typed edges (source-only; mentions auto-derived)

## [2026-07-10] ingest | PHYLO-STATS-001-Evidence.md → phylo-stats-001-evidence (source) + tree-statistics (concept)
   Ingested the Tree Statistics unit (PHYLO-STATS-001): GetLeaves / CalculateTreeLength / GetTreeDepth —
   whole-tree descriptive summaries (leaf count, total branch length = Σ all edges, height in edges).
   Judged genuinely distinct from [[tree-comparison-metrics]] (descriptive summaries vs compare/query
   operations) so created a modest dedicated concept [[tree-statistics]] rather than a source-only
   cross-link. Cross-linked reciprocally with [[tree-comparison-metrics]] (added a paragraph there) and
   to the PHYLO family ([[phylogenetic-bootstrap-support]], [[evolutionary-distance-matrix]],
   [[phylo-newick-001-evidence]]). Sources: Wikipedia Tree-(graph-theory)/Tree-(ADT) + Biopython
   BaseTree + DendroPy Tree.length() + Minimum-evolution. One assumption (null PhyloNode ↔ empty-tree
   height −1). Contradictions: none. Follow-ups: none.
   graph: +2 nodes (phylo-stats-001-evidence source, tree-statistics concept), +1 typed edge (tree-statistics relates_to test-unit-registry)

## [2026-07-10] ingest | PHYLO-TREE-001-Evidence.md → phylo-tree-001-evidence (source) + distance-based-tree-construction (concept)
   Ingested the Tree Construction unit (PHYLO-TREE-001): UPGMA + Neighbor-Joining — the PHYLO family's
   tree-building core that consumes a distance matrix and emits a PhyloNode tree (BuildTree /
   BuildTreeFromMatrix). Created the dedicated concept [[distance-based-tree-construction]] that prior
   PHYLO ingests explicitly flagged as MISSING/WANTED (UPGMA rooted-ultrametric height=d/2 clock;
   NJ minimum-Q additive-topology, negative branches preserved, midpoint-rooted final join). Wired it
   as the family hinge: concept depends_on [[evolutionary-distance-matrix]]; added reciprocal
   depends_on edges + wikilinks on [[phylogenetic-bootstrap-support]] (wraps it per replicate),
   [[tree-comparison-metrics]] and [[tree-statistics]] (operate on its output); linked from
   [[evolutionary-distance-matrix]]. Sources: Wikipedia UPGMA/Neighbor-joining/Phylogenetic-tree +
   Saitou & Nei 1987 + Sokal & Michener 1958 + Felsenstein 2004. Worked oracles: UPGMA 5S-rRNA
   (root 16.5, tips all 16.5) + NJ 5-taxon (Q₁(a,b)=−50, δ(a,u)=2…δ(e,w)=1). Deviations None (§8).
   Contradictions: none. Follow-ups: none.
   graph: +2 nodes (phylo-tree-001-evidence source, distance-based-tree-construction concept), +4 typed edges (distance-based-tree-construction relates_to test-unit-registry + depends_on evolutionary-distance-matrix; bootstrap/comparison/statistics each depends_on distance-based-tree-construction)

## [2026-07-10] ingest | POP-ANCESTRY-001-Evidence.md → pop-ancestry-001-evidence (source) + ancestry-estimation-admixture (concept)
   First population-genetics POP-* unit. Supervised/projection ADMIXTURE: estimate ancestry
   proportions Q by FRAPPE EM (Eq. 4) with fixed reference allele frequencies F; log-likelihood
   Eq. 2 under simplex Σ_k q_ik=1, convergence Eq. 5 ε=10⁻⁴, O(IJK²). Sources: Alexander, Novembre
   & Lange 2009 (Genome Research) + Alexander & Lange 2011 (supervised mode, Springer-gated) +
   ADMIXTURE 1.4 Manual §2.10/§2.14. Oracles: symmetric K=2/J=2 panel g=[2,0] → one EM iter
   (0.8,0.2) exactly → converges (1.0,0.0); single-SNP g=2→(0.9,0.1); identical panels uniform
   fixed point. Distinct from the PHYLO-* tree family (mixture-weight decomposition, not a tree);
   label non-identifiability pinned by fixed labelled panels. Two research-grade assumptions
   (maxIterations+ε; skip missing genotype). Contradictions: none. Follow-ups: none.
   graph: +2 nodes (pop-ancestry-001-evidence source, ancestry-estimation-admixture concept), +1 typed edge (ancestry-estimation-admixture relates_to test-unit-registry)

## [2026-07-10] ingest | POP-DIV-001-Evidence.md → pop-div-001-evidence (source) + genetic-diversity-statistics (concept)
   Second population-genetics POP-* unit, sibling of POP-ANCESTRY-001. Diversity-statistics panel:
   nucleotide diversity π = Σd_ij/(C(n,2)·L) (Nei & Li 1979), Watterson's θ_W = S/a_n with harmonic
   a_n (Watterson 1975), Tajima's D = (k̂ − S/a_1)/√(e_1·S+e_2·S(S−1)) (Tajima 1989; k̂ = pairwise-
   difference COUNT not per-site π; D<0 sweep/expansion, D>0 balancing/contraction), and Nei-1978
   unbiased gene diversity/heterozygosity (H_obs = n/(n−1)·H_exp ≡ π for haploid). Oracle: Wikipedia
   Tajima's D example n=5/L=20/S=4 → k̂=2.0/π=0.1/θ_W≈0.096/D≈0.273 (TD-C01/TD-C02). Guards n<3→D
   undefined→0, S=0/monomorphic→all 0, Var≤0→D=0. All formulae exact-match; distinct from the
   ancestry anchor (variation amount vs ancestry decomposition). Contradictions: none. Follow-ups: none.
   graph: +2 nodes (pop-div-001-evidence source, genetic-diversity-statistics concept), +2 typed edges (genetic-diversity-statistics relates_to test-unit-registry; genetic-diversity-statistics relates_to ancestry-estimation-admixture)
## [2026-07-10] ingest | POP-FREQ-001-Evidence.md → pop-freq-001-evidence (source) + allele-genotype-frequencies (concept)
   Foundational population-genetics POP-* primitive — allele/genotype frequencies, minor allele
   frequency (MAF), MAF filtering; the numeric substrate under the POP family (produces the fixed F
   for POP-ANCESTRY-001 and the per-site p_i for POP-DIV-001's heterozygosity term). Biallelic
   allele freq p=f(AA)+½f(AB)/q=f(BB)+½f(AB) via counts total=2·(n_AA+n_AB+n_BB), major=2·n_AA+n_AB,
   minor=2·n_BB+n_AB (INV p+q=1, major+minor=total; four-o'clock oracle 49/42/9→0.70/0.30); VCF/PLINK
   dosage alt_freq=Σg/(2n), MAF=min(alt_freq,1−alt_freq) (INV 0≤MAF≤0.5, symmetric alt 0.7→MAF 0.3;
   monomorphic→0; 50/50→0.5); MAF filter inclusive [minMAF,maxMAF] band, HapMap/common(>0.05)-vs-
   rare(<0.05) thresholds. Edge cases: zero samples→(0,0), empty vector→MAF 0, negative count→
   ArgumentOutOfRangeException. Scope biallelic counting/normalization only (no HWE, no multiallelic,
   no phasing). All exact-match (Wikipedia Allele/Minor-allele/Genotype frequency + Gillespie 2004 +
   NDSU). Contradictions: none. Follow-ups: none.
   graph: +2 nodes (pop-freq-001-evidence source, allele-genotype-frequencies concept), +3 typed edges (allele-genotype-frequencies relates_to test-unit-registry; relates_to ancestry-estimation-admixture; relates_to genetic-diversity-statistics)
## [2026-07-10] ingest | POP-FST-001-Evidence.md → pop-fst-001-evidence (source) + population-differentiation-fst (concept)
   Population-genetics POP-* differentiation unit — Fst (fixation index), F-statistics (Fis/Fit/Fst),
   pairwise Fst. Consumes per-population allele frequencies from POP-FREQ-001. Wright 1965 variance
   Fst = σ_S²/(pBar(1−pBar)); two-pop size-weighted pBar=(n1·p1+n2·p2)/(n1+n2),
   σ_S²=(n1(p1−pBar)²+n2(p2−pBar)²)/(n1+n2); multi-locus ratio-of-sums Σσ_S²/Σhet — computes the
   population PARAMETER from known allele freqs, explicitly NOT the Weir & Cockerham 1984 θ estimator
   (no ANOVA/finite-sample bias correction). F-statistics heterozygosity partition Fis=1−Hi/Hs,
   Fit=1−Hi/Ht, Fst=1−Hs/Ht with exact identity (1−Fit)=(1−Fis)(1−Fst). INV Fst∈[0,1] (0=panmixia,
   1=fixed differences), Fis∈[−1,1] (negative under excess heterozygotes), pairwise matrix symmetric +
   zero diagonal, NOT a metric (fails triangle inequality). Oracles: fixed p1=1/p2=0→1.0 exactly,
   pop1=(.9,.8)/pop2=(.1,.2)→1/2, unequal sizes 0.006274…, components 1/19,1/13,1/39, excess-het
   negative Fis=−2/3, pairwise cells 1/99,4/21,3/25; ref values Cavalli-Sforza 1994 + Elhaik 2012
   HapMap; Hartl-Clark interpretation bands. Edge cases denominator 0 (empty/both-fixed-same/
   monomorphic)→return 0. 25 tests (−1 dup +4 new). New concept created (distinct from POP-FREQ
   counting, POP-DIV within-sample diversity, POP-ANCESTRY decomposition). Contradictions: none.
   Follow-ups: none.
   graph: +2 nodes (pop-fst-001-evidence source, population-differentiation-fst concept), +3 typed edges (population-differentiation-fst relates_to test-unit-registry; depends_on allele-genotype-frequencies; relates_to genetic-diversity-statistics)
## [2026-07-10] ingest | POP-HW-001-Evidence.md → pop-hw-001-evidence (source) + hardy-weinberg-equilibrium-test (concept)
   Population-genetics POP-* unit — Hardy-Weinberg equilibrium (HWE) chi-square test. Consumes
   genotype counts from POP-FREQ-001 (whose scope explicitly leaves the HWE test to this unit).
   Expected genotype freqs p²/2pq/q² (Hardy 1908 / Weinberg 1908); allele freq p=(2·n_AA+n_Aa)/(2n),
   q=1−p; expected counts E={p²n,2pqn,q²n}; Pearson χ²=Σ(O−E)²/E over 3 genotype classes; df=1
   (#genotypes−#alleles=3−2); p-value via chi-square CDF (lower-incomplete-gamma approx); default
   α=0.05 critical value 3.841. TestHardyWeinberg returns InEquilibrium/ChiSquare/PValue. Oracles:
   Ford moth (1469,138,5)→p≈0.954/χ²≈0.83→in-eq, perfect (25,50,25)→χ²=0, excess-het (10,80,10)→
   χ²=36≫3.84→out-of-eq, zero samples→InEquilibrium true/PValue 1, fixed (100,0,0)→in-eq, all-het
   (0,100,0)→out-of-eq. Edge cases: zero n→PValue 1 (no evidence against H₀, hypothesis-testing
   framework not ad-hoc), expected-0 term skipped (div-by-zero guard). Scope = biallelic chi-square
   goodness-of-fit only; exact test (Wigginton 2005) and multiallelic loci noted out of scope. New
   concept created (distinct: a hypothesis test on counts vs POP-FREQ counting/normalization).
   Cross-linked allele-genotype-frequencies scope note. Contradictions: none. Follow-ups: none.
   graph: +2 nodes (pop-hw-001-evidence source, hardy-weinberg-equilibrium-test concept), +2 typed edges (hardy-weinberg-equilibrium-test relates_to test-unit-registry; depends_on allele-genotype-frequencies)
## [2026-07-10] ingest | POP-LD-001.md → pop-ld-001-evidence (source) + linkage-disequilibrium (concept)
   Population-genetics POP-* unit — linkage disequilibrium between two loci (CalculateLD: D, D', r²)
   + haplotype-block detection (FindHaplotypeBlocks). Consumes allele/haplotype frequencies from
   POP-FREQ-001. Sources: Wikipedia Linkage-disequilibrium (D=p_AB−p_A·p_B; Lewontin 1964 D'=|D|/D_max
   sign-branched clamped [0,1]; Hill & Robertson 1968 r²=D²/(p_A·q_A·p_B·q_B); diploid-frequency result
   R_AB=r_AB Wright 1933 ⇒ r² computable WITHOUT phase) + Wikipedia Haplotype-block (Gabriel 2002 /
   Patil 2001). Implementation: r² = squared Pearson correlation of 0/1/2 genotype dosage vectors
   Cov²/(Var·Var); D from diploid covariance Cov=2D ⇒ D=Cov/2; FindHaplotypeBlocks = simplified
   adjacent-pair Gabriel (consecutive r²≥threshold, default 0.7, ≥2 variants). Oracles perfect LD→r²≈1,
   no LD→r²≈0, anti-correlation→r²=1/D'=1 (sign-blind), block single→none/two-high→one/two-low→none/
   all-strong→one span/non-contiguous→multiple. INV 0≤r²≤1, 0≤|D'|≤1, empty→r²=0/D'=0, monomorphic
   (zero-variance denominator)→r²=0 guarded, distance+IDs preserved, blocks Start≤End/≥2/non-overlapping/
   ordered. Scope = two-biallelic-loci r²/D' + adjacent-pair blocks only (no full LD matrix, no
   phasing/EM, no decay-curve fit, no exact Gabriel CI). New concept created (distinct: pairwise
   inter-locus association vs per-locus counting/diversity/differentiation/HWE-test). Cross-linked
   allele-genotype-frequencies scope note. Contradictions: none. Follow-ups: none.
   graph: +2 nodes (pop-ld-001-evidence source, linkage-disequilibrium concept), +2 typed edges (linkage-disequilibrium relates_to test-unit-registry; depends_on allele-genotype-frequencies)

## [2026-07-10] ingest | docs/Evidence/POP-ROH-001-Evidence.md → pop-roh-001-evidence (source) + runs-of-homozygosity-inbreeding (concept)
   POP-ROH-001: runs of homozygosity (FindROH) + genomic inbreeding coefficient F_ROH. Window-free
   consecutive-runs scan (Marras 2015/detectRUNS): grow runs over position-sorted 0/1/2 genotypes,
   terminate on maxOppRun exceeded or gap>maxGap, retain only if minSNP AND minLengthBps pass (PLINK
   --homozyg-snp 100 AND --homozyg-kb 1000). F_ROH = ΣL_roh/L_auto (McQuillan 2008; L_auto≈2,674 Mb;
   oracle 20/100 Mb→0.20, whole-genome→1.0). Two API-encoding assumptions (0/1/2 encoding; missing
   handling out of scope). New concept created — genuinely distinct per-individual segment detection
   vs the other POP siblings (frequencies/diversity/Fst/LD/HWE). Cross-linked all POP concepts +
   ancestry family anchor. Contradictions: none. Follow-ups: none.
   graph: +2 nodes (pop-roh-001-evidence source, runs-of-homozygosity-inbreeding concept), +2 typed edges (runs-of-homozygosity-inbreeding relates_to test-unit-registry; relates_to ancestry-estimation-admixture)

## [2026-07-10] ingest | docs/Evidence/POP-SELECT-001-Evidence.md — Selection-signature detection (iHS/EHH)
   Ingested the POP-SELECT-001 evidence artifact (integrated Haplotype Score iHS + Extended Haplotype
   Homozygosity EHH scan; CalculateEhh/CalculateIHS/StandardizeIHS/ScanForSelection). Created source
   summary [[pop-select-001-evidence]] and a NEW dedicated concept [[selection-scan-ihs-ehh]] —
   genuinely distinct from the POP siblings: a haplotype-length/decay statistic, not frequencies,
   diversity, Fst, HWE, LD, or ROH. Synthesized the EHH→iHH→iHS→scan pipeline, the trapezoidal iHH
   with the 0.05 cutoff, the Voight vs selscan sign-convention pitfall (ln(iHH_A/iHH_D) vs its
   inverse), worked oracles (rehh F1205400 −1.978569274; constructed panel ln(0.25)=−1.386294361),
   invariants and edge cases. Updated hub [[algorithm-validation-evidence]] (frontmatter sources +
   body pop-* link list) and wiki/index.md (source + concept entries). Cross-linked all POP concepts.
   Contradictions: none (the sign difference is a documented convention). Follow-ups: none.
   graph: +2 nodes (pop-select-001-evidence source, selection-scan-ihs-ehh concept), +3 typed edges (selection-scan-ihs-ehh relates_to test-unit-registry; depends_on allele-genotype-frequencies; relates_to linkage-disequilibrium)

## [2026-07-10] ingest | docs/Evidence/PRIMER-TM-001-DIMER-Evidence.md → primer-tm-001-dimer-evidence (source) + primer-dimer-thermodynamics-tm (concept)
   First PCR primer-design PRIMER-* / MolTools family unit. New concept
   [[primer-dimer-thermodynamics-tm]]: self-/hetero-dimer Tm via Primer3 ntthal thermodynamic
   alignment over the SantaLucia & Hicks 2004 DNA nearest-neighbour model — 10 WC NN stacks +
   initiation + terminal A·T penalty + symmetry, bimolecular Tm with x=1(palindrome)/x=4 factor,
   [Na+] salt correction, full non-contiguous dimer DP (mismatch/loop/bulge/tstack2 overhang)
   reproducing primer3-py 2.3.0 to machine precision; poly-A/invalid → null/NaN. Cross-linked as
   the DNA counterpart of the RNA Turner-2004 folding [[rna-base-pairing]] / [[pre-mirna-hairpin-detection]].
   Updated hub [[algorithm-validation-evidence]] (frontmatter sources + body evidence-link list +
   distinct-concept-anchor list) and wiki/index.md (source + concept entries).
   Contradictions: none. Follow-ups: none.
   graph: +2 nodes (primer-tm-001-dimer-evidence source, primer-dimer-thermodynamics-tm concept), +1 typed edge (primer-dimer-thermodynamics-tm relates_to test-unit-registry)

## [2026-07-10] ingest | docs/Evidence/PRIMER-TM-001-Evidence.md → primer-tm-001-evidence (source) + primer3-weighted-penalty-objective (concept)
   Base PRIMER-TM-001 unit = Primer3 weighted per-primer penalty (objective function) `p_obj_fn` — a
   selection/scoring algorithm, NOT a Tm calc despite the unit ID. Distinct from the sibling dimer-Tm
   [[primer-dimer-thermodynamics-tm]] (same unit ID), so created a new concept
   [[primer3-weighted-penalty-objective]] and cross-linked both directions (penalty consumes Tm/self-align
   scores as terms). Default objective collapses to |Tm−60|+|len−20|. Updated hub
   [[algorithm-validation-evidence]] (frontmatter sources + body evidence-link list) and wiki/index.md
   (source + concept entries).
   Contradictions: none. Follow-ups: none.
   graph: +2 nodes (primer-tm-001-evidence source, primer3-weighted-penalty-objective concept), +2 typed edges (primer3-weighted-penalty-objective relates_to test-unit-registry; relates_to primer-dimer-thermodynamics-tm)

## [2026-07-10] ingest | docs/Evidence/PRIMER-TM-001-HAIRPIN-Evidence.md → primer-tm-001-hairpin-evidence (source)
   Hairpin/secondary-structure extension of PRIMER-TM-001 (intramolecular self-fold MFE +
   unimolecular hairpin Tm). OVERLAP resolved by ENRICHING [[primer-dimer-thermodynamics-tm]]
   (same ntthal/thal.c engine, same NnUnifiedParams stem stacks) with a new "Intramolecular
   hairpin self-folding" section (Table 4 loop increments, Eq. 11 no-concentration Tm,
   Jacobson-Stockmayer, exclusions, hairpin oracles) rather than creating a new concept.
   Updated hub [[algorithm-validation-evidence]] (frontmatter sources + body evidence-link list)
   and wiki/index.md (source + concept entries).
   Contradictions: none. Follow-ups: length-3/4 triloop/tetraloop bonus tables remain opt-in
   caller increments (supplementary material not bundled).
   graph: +1 node, +1 typed edge

## [2026-07-10] ingest | docs/Evidence/PRIMER-TM-001-NN-Evidence.md → primer-tm-001-nn-evidence (source)
   Per-oligo NN salt-corrected design Tm (opt-in) for PRIMER-TM-001. Heavy overlap with the
   existing concept: ENRICHED [[primer-dimer-thermodynamics-tm]] with a "Per-oligo design Tm and
   salt corrections" section (Eq. 3 per-primer Tm; Owczarzy 2004 monovalent-quadratic + Owczarzy
   2008 divalent-Mg²⁺/dNTP 1/Tm corrections, Biopython salt_correction methods 5–7; DNA_IMM/DNA_DE
   mismatch/dangling tables → NnInternalMismatch/NnDanglingEnd; complement-not-revcomp Tm_NN
   convention) rather than creating a new concept. Updated hub [[algorithm-validation-evidence]]
   (frontmatter sources + body evidence-link list) and wiki/index.md (source + concept entries).
   Contradictions: none (Biopython tables verified as faithful transcriptions of the primaries).
   Follow-ups: none.
   graph: +1 node, +1 typed edge
- 2026-07-10 — ingest docs/Evidence/PRIMER-TM-001-SPECIAL-LOOP-Evidence.md. Created source page
   [[primer-tm-001-special-loop-evidence]] (bundled special tri/tetraloop hairpin bonus tables:
   libprimer3 triloop.*/tetraloop.* config + thal.c calc_hairpin application + primer3-py 2.3.0
   oracles). Enriched concept [[primer-dimer-thermodynamics-tm]] rather than creating a new one —
   this unit completes the previously opt-in triloop/tetraloop increment it already flagged;
   updated the hairpin-section special-loop paragraph, failure-modes contract, intro, frontmatter
   (+source, +typed edge), and updated hub [[algorithm-validation-evidence]] (frontmatter sources +
   body evidence-link list) and wiki/index.md (source + concept entries). Contradictions: none
   (all values verbatim from libprimer3 + machine-precision verified vs primer3-py). Follow-ups: none.
   graph: +1 node, +1 typed edge

## [2026-07-10] ingest | docs/Evidence/PROBE-DESIGN-001-Evidence.md → probe-design-001-evidence (source) + 1 concept
   TaqMan 5'-nuclease hydrolysis-probe design rules (opt-in over the unchanged generic probe
   designer). Judged genuinely distinct from the primer units (probe-specific hard constraints:
   no 5'-G reporter-quench, more-C-than-G + antisense strand fallback, ≥4-G run, GC 30-80%,
   length 18-22, probe Tm ≥ primer Tm + 10) so created new PROBE-family anchor concept
   [[taqman-probe-design-rules]]; it reuses the PRIMER-TM-001-validated salt-adjusted Tm engine
   (relates_to [[primer-dimer-thermodynamics-tm]]). Updated hub [[algorithm-validation-evidence]]
   (frontmatter sources + body evidence-link list + own-concept enumeration) and wiki/index.md
   (source + concept entries). Contradictions: none (four vendor/reference sources corroborate
   point-for-point). Follow-ups: none.
   graph: +2 nodes, +2 typed edges

## [2026-07-10] ingest | docs/Evidence/PROBE-DESIGN-001-LNA-Evidence.md → probe-design-001-lna-evidence (source)
   LNA (locked nucleic acid) Tm-adjustment variant of PROBE-DESIGN-001: an LNA-adjusted
   nearest-neighbour Tm (McTigue/Peterson/Kahn 2004 — 32 LNA+DNA:DNA NN increments in cal/mol;
   internal LNA raises Tm/specificity → shorter MGB-style 13-20 nt probes). Additive-increment
   model onto the library's SantaLucia-1998-unified DNA NN engine (same as PRIMER-TM-001);
   terminal-LNA/non-ACGT/out-of-range → not-computable. Enriched existing concept
   [[taqman-probe-design-rules]] (new LNA section, LNA base-NN assumption, +typed edge relates_to
   [[primer-dimer-thermodynamics-tm]]) rather than creating a new concept — genuinely the same
   PROBE unit's modified-base Tm variant. Updated hub [[algorithm-validation-evidence]]
   (frontmatter sources + body evidence-link list) and wiki/index.md (source entry). Oracle
   CCATTGCTACC LNA@4 → Tm 63.528 °C vs all-DNA 59.692 (+3.84), MELTING mct04 63.614 to 0.086 °C.
   Contradictions: none (McTigue + MELTING + rmelting agree; 0.086 °C = documented base-NN-set diff).
   Follow-ups: none.
   graph: +1 node, +1 typed edge

## [2026-07-10] ingest | docs/Evidence/PROBE-VALID-001-Evidence.md → probe-valid-001-evidence (source) + 1 concept
   PROBE-VALID-001 = hybridization-probe off-target specificity scan via gapped Smith–Waterman
   local alignment (replaces earlier ungapped Hamming), + Kane-2000 0.75 identity threshold + opt-in
   Karlin–Altschul E-value/bit-score/λ statistics. Judged genuinely distinct from
   [[taqman-probe-design-rules]] (composition rules) → new concept [[probe-offtarget-specificity-scan]],
   cross-linked both ways as the specificity-checking sibling. Updated hub
   [[algorithm-validation-evidence]] (frontmatter sources + body evidence-link list), enriched
   [[taqman-probe-design-rules]] (sibling cross-link), wiki/index.md (source + concept entries).
   Oracles: indel copy ACGTAC-GTACGT 12/12=1.0 found by gapped / missed by ungapped; trimmed
   indel+mismatch 10/12=0.8333; λ(+1/−3)=1.3740631 (≈ blastn 1.37), bit 59.9627 / E 1.7802e−14.
   Contradictions: none (SW + BLAST gapped/ungapped + Karlin–Altschul + Kane mutually consistent).
   Follow-ups: none.
   graph: +2 nodes, +3 typed edges

## [2026-07-10] ingest | docs/Evidence/PROTMOTIF-CC-001-Evidence.md → protmotif-cc-001-evidence (source) + 1 concept
   Coiled-coil prediction (PROTMOTIF-CC-001): heptad a/d hydrophobic-core occupancy predictor
   `ProteinMotifFinder.PredictCoiledCoils` — per-window fraction of a/d core positions ∈ {I,L,V}
   maximised over 7 heptad registers, contiguous runs ≥21 residues (3 heptads) emitted with peak
   Score∈[0,1]; defaults window 28 (Lupas 1991) / threshold 0.5 / min-region 21 (Mason & Arndt).
   Judged genuinely distinct — first ingested unit of the ProteinMotif family, separate from the
   ProteinPred disorder/features family → new concept [[coiled-coil-prediction]], cross-linked from
   the [[protein-low-complexity-seg]] anchor as a sibling sequence-only protein-feature heuristic.
   Enriched with the algorithm doc docs/algorithms/ProteinMotif/Coiled_Coil_Prediction.md (INV-01..05).
   Updated hub [[algorithm-validation-evidence]] (frontmatter sources + body evidence-link + concept
   list), wiki/index.md (source + concept entries).
   Oracles: (LAALAAA)×5→(0,34,1.0); all-Gly/no-{I,L,V}→none; (LAAAAAA)→0.5 threshold boundary.
   Deviations: COILS 21×20 PSSM deliberately omitted (weights not retrievable → use COILS/Paircoil2);
   {I,L,V}-only core-set is source-verbatim. Contradictions: none. Follow-ups: none.
   graph: +2 nodes, +2 typed edges

## [2026-07-10] ingest | PROTMOTIF-COMMON-001-Evidence.md → protmotif-common-001-evidence (source) + 1 concept
   Common motif finding (`ProteinMotifFinder.FindCommonMotifs`): whole-dictionary scan of a fixed
   built-in `CommonMotifs` catalog of canonical PROSITE patterns (PS00001 N-glycosylation, PS00005
   PKC / PS00006 CK2 phospho sites, PS00016 RGD, PS00017 ATP/GTP P-loop), aggregating each hit with
   its accession/name. Judged genuinely distinct — second ingested ProteinMotif unit, a degenerate
   PROSITE-pattern dictionary scan (not the windowed a/d heuristic) → new concept
   [[common-protein-motifs]], the protein fixed-catalog analogue of the DNA [[regulatory-element-detection]]
   and sibling of [[coiled-coil-prediction]]; distinguished from the caller-supplied-set DNA exact
   [[known-motif-search]]. Updated hub [[algorithm-validation-evidence]] (frontmatter sources + body
   evidence-link + concept list) and wiki/index.md (concept entry).
   Oracles: NFTA (PS00001) / N-P-[ST] Pro-exclusion no-match / SAR / SAAE+SDED / GXXXXGKS / RGD+RGD
   overlap; 0-based-inclusive MotifMatch vs PROSITE 1-based; overlaps all reported.
   Deviations: none. Assumption: 0- vs 1-based coordinate origin (API shape, no correctness effect).
   Contradictions: none. Follow-ups: sibling units PROTMOTIF-FIND-001 / PROTMOTIF-PATTERN-001 (general
   PROSITE engine) not yet ingested.
   graph: +2 nodes, +2 typed edges

## [2026-07-10] ingest | docs/Evidence/PROTMOTIF-DOMAIN-001-Evidence.md → protmotif-domain-001-evidence (source) + 1 concept
   Third ProteinMotif-family Evidence file (after PROTMOTIF-CC coiled-coil, PROTMOTIF-COMMON common motifs).
   Created the genuinely-distinct concept protein-domain-and-signal-peptide-prediction — the ProteinMotif
   family's domain + signal-peptide unit, covering three algorithms on ProteinMotifFinder: (1) FindDomains,
   a deterministic PROSITE-PATTERN domain scan (PS00028 C2H2 zinc finger/PF00096, PS00017 P-loop-Walker A/
   PF00069, PS00678 WD40 14-element/15-residue signature/PF00400; ScanProsite→regex translation; real
   GBB1_HUMAN P62873 WD40 positive at 0-based 69/156/284) — SH3(PS50002)/PDZ(PS50106) are PROSITE PROFILEs
   with NO deterministic pattern so are excluded (honest residual; prior unsourced ad-hoc regexes removed);
   (2) the opt-in FindDomainsByHmm / Plan7ProfileHmm engine reproducing the HMMER3 pipeline over 3 bundled
   CC0 Pfam HMMs (PF00018 SH3/PF00595 PDZ/PF00400 WD40) — Viterbi/Forward log-odds (exact 1e-9 on a
   hand-built 2-symbol HMM = 0.5187937934 nats), hmmsearch-parity local-multihit pre_score scored vs
   Swiss-Prot bg->f not COMPO (SH3 68.7097/PDZ 84.8629/WD40 213.4120 bits, ~1e-5-bit parity), null2
   biased-composition correction (omega=1/256), Gumbel(MSV/Viterbi)/exponential(Forward) E-values with
   E=P·Z from STATS LOCAL, p7_domaindef multi-domain envelope decomposition (GBB1/PF00400→7 β-propeller
   blades, coords exact) + stochastic-traceback single-linkage clustering (Easel LCG seed 42) for
   closely-overlapping tandems — all cross-checked against pyhmmer 0.12.1 ground truth and a from-scratch
   Python re-derivation; (3) PredictSignalPeptide, the von Heijne tripartite n(K/R+)/h(hydrophobic
   α-helix)/c(polar) model with the −1,−3 rule {A,G,S}, score (nScore+2·hScore+cScore)/4 and evidence-based
   detection constraints (nScore>0 & hScore≥0.5, replacing the eliminated 0.4 threshold), Probability=Score.
   Oracles: C2H2 AAAACXXCXXXLXXXXXXXXHXXXHAAA→4..24, P-loop AAAAGXXXXGKSAAAA→4..11, signal
   MKRLLLLLLLLLLLLLLLLLLASAGDDDEEEFFF→detected cleavage≈25. Concise source page for the artifact. Updated
   hub [[algorithm-validation-evidence]] (frontmatter sources + body evidence-link list) and wiki/index.md
   (source entry + concept entry). Added a reciprocal ProteinMotif-family cross-link from
   [[common-protein-motifs]]. Deviations: the six previously-listed items are all RESOLVED design decisions
   (1:2:1 weights, evidence-based constraints, Probability=Score, strict {A,G,S}, PROSITE-pattern scope,
   FindDomains naming), not open assumptions. Contradictions: none — the encyclopedic + PROSITE/Pfam +
   von Heijne + HMMER/Easel/Durbin sources agree. Honest residuals: SH3/PDZ profile-only; only 3 CC0 HMMs
   bundled; MSV/bias prefilters and exact-RNG trace-ensemble bit parity not reproduced (research-grade).
   Follow-ups: PROTMOTIF-FIND-001 / PROTMOTIF-PATTERN-001 (general PROSITE engine), transmembrane-helix
   and other ProteinMotif units not yet ingested; noted the sibling protmotif-common-001-evidence source
   entry was absent from wiki/index.md (prior-ingest gap, left as-is).
   graph: +2 nodes, +2 typed edges (relates_to test-unit-registry, common-protein-motifs)
- 2026-07-10 — ingest docs/Evidence/PROTMOTIF-FIND-001-Evidence.md (Protein Motif Search, Pattern-based). Genuinely distinct from the fixed-dictionary [[common-protein-motifs]]: this is the GENERAL engine `ProteinMotifFinder.FindMotifByPattern` that takes an arbitrary caller-supplied PROSITE pattern, of which FindCommonMotifs is one application. New concept [[protein-motif-pattern-search]] covering three primitives — `ConvertPrositeToRegex` (PROSITE→regex element map), `FindMotifByPattern` (overlapping-match discovery via zero-width lookahead `(?=(pattern))` per ScanProsite; 0-based Start/End; case-insensitive; empty/null→empty; invalid regex handled gracefully), and information-content scoring `CalculateMotifScore` IC=Σlog₂(20/allowed_count) + `CalculateEValue` E=(N−L+1)·2^(−IC) (Schneider & Stephens 1990; the earlier BLAST/Altschul citation was explicitly REMOVED because the E-value is a direct combinatorial probability, not Karlin–Altschul EVD). Pins the PROSITE-pattern catalog and records two FIXED implementation bugs — PS00007 loosened `.{2,3}`→exact `[RK].{2}[DE].{3}Y`, and PS00018 EF-hand (`x`→`{W}` at pos 2 + restored dropped trailing `[LIVMFYW]`); five non-PROSITE linear motifs (NLS1/NES1/SIM1/WW1/SH3_1) re-derived from primary literature. Concise source page [[protmotif-find-001-evidence]]. Updated hub [[algorithm-validation-evidence]] (frontmatter sources + body evidence-link list + source_commit→HEAD) and wiki/index.md (source entry + concept entry). Corrected a forward-reference in [[common-protein-motifs]]: its "general engine" pointer said "PROTMOTIF-PATTERN unit" but the general engine is actually this unit, PROTMOTIF-FIND-001 (`FindMotifByPattern`) — repointed to [[protein-motif-pattern-search]]. Deviations/assumptions: all eliminated per the Evidence change history (patterns corrected, non-PROSITE patterns literature-verified, heuristic scoring replaced by IC scoring, overlapping-match lookahead implemented); only the 0-based-vs-ScanProsite-1-based coordinate convention stands (no correctness effect). Contradictions: none — PROSITE/ScanProsite + Schneider & Stephens + the five primary refs agree; flagged the stale sibling-page forward-reference as noted above. Follow-ups: PROTMOTIF-PATTERN-001 (if a distinct unit exists) and remaining ProteinMotif units (transmembrane-helix etc.) not yet ingested.
   graph: +2 nodes, +2 typed edges (relates_to test-unit-registry, common-protein-motifs)
- 2026-07-10 — ingest docs/Evidence/PROTMOTIF-LC-001-Evidence.md (Protein low-complexity region detection, SEG). SAME ALGORITHM as DISORDER-LC-001 — SEG (Wootton & Federhen 1993), identical defaults W=12/K1=2.2/K2=2.5, identical Shannon-entropy bits/residue complexity `−Σpᵢlog₂pᵢ` (max log₂20≈4.322), identical two-stage trigger(≤K1)/extend(≤K2) scan. NOT a different low-complexity method: this is the ProteinMotif-family registration of the SEG unit, a second Evidence file tracing the same method. Therefore NO new concept — enriched/cross-linked the existing anchor [[protein-low-complexity-seg]] (added PROTMOTIF-LC-001-Evidence to its frontmatter sources, added a body paragraph documenting the second unit, added a relates_to test-unit-registry graph edge sourced from protmotif-lc-001-evidence, bumped source_commit→HEAD). Wrote concise source page [[protmotif-lc-001-evidence]] recording the sources (NCBI ncbi-seg manpage + blast_seg.c constants/s_Entropy/s_LnPerm/lnfact[] + SeqComplex ce/cwf + universalmotif + Pei & Grishin 2005 + Mier et al.), the worked-window oracle table (homopolymer H=0 / 11A1B 0.413817 / 10A2B 0.650022 / 6A6B 1.0 / 12-distinct log₂12≈3.585), corner cases, and the two assumptions (Shannon bits/residue form per manpage units; short-<W→empty). Updated hub [[algorithm-validation-evidence]] (frontmatter sources + body evidence-link list + source_commit→HEAD) and wiki/index.md (Sources-section entry only). Contradictions: none — the artifact fully agrees with the DISORDER-LC-001 evidence; the only new primary vs DISORDER-LC-001 is Pei & Grishin 2005. Follow-ups: none specific; remaining ProteinMotif units (transmembrane-helix etc.) still unenumerated.
   graph: +1 node, +1 typed edge (relates_to test-unit-registry)
- 2026-07-10 — ingest docs/Evidence/PROTMOTIF-PATTERN-001-Evidence.md (Protein Pattern Matching Methods: FindMotifByPattern, FindMotifByProsite, ConvertPrositeToRegex, FindDomains). SECOND Evidence over the SAME PROSITE→regex engine as PROTMOTIF-FIND-001 — NO new concept; enriched [[protein-motif-pattern-search]] (added the end-to-end `FindMotifByProsite` primitive row; added `A(n)`→`A{n}` and trailing-`.` terminator to the conversion table; added a "PA-line grammar corner cases" subsection — ranges only on `x` (`A(2,4)` invalid, `A(3)` valid), trailing period terminates, **reject the `*` Kleene star with FormatException** since `<{C}*>` is a ScanProsite query extension not PA-line grammar; added the second source + source_commit→HEAD; noted the PROTMOTIF-PATTERN-001 revalidation in the anchor sentence and References). Wrote concise source page [[protmotif-pattern-001-evidence]] with the exact IC oracles (RGD 3·log₂20≈12.965784284662087 bits, class `[ST]`→log₂10≈3.321928094887362, wildcard→0) and the PS00001/05/16/17/29 worked-example regex table. Updated hub [[algorithm-validation-evidence]] (frontmatter sources + body evidence-link list + source_commit→HEAD) and wiki/index.md (Sources-section entry + a revalidation clause on the [[protein-motif-pattern-search]] concept entry). Two ASSUMPTIONs in the artifact (lookahead overlap enumeration is a repo contract not PROSITE-mandated; combinatorial E=(N−L+1)·2^(−IC) is a model quantity, not the ScanProsite Swiss-Prot-frequency E-value). Contradictions: none — fully consistent with PROTMOTIF-FIND-001 (same engine, IC per Schneider & Stephens 1990, pattern-is-regex per De Castro 2006). Follow-ups: FindDomains is jointly owned by [[protein-domain-and-signal-peptide-prediction]] (PROTMOTIF-DOMAIN-001) — added a medium-confidence relates_to edge rather than duplicating that unit's coverage.
   graph: +1 node, +3 typed edges (relates_to test-unit-registry ×1 from new source; relates_to protein-domain-and-signal-peptide-prediction; +existing concept node gains edges)
- 2026-07-10 — ingest docs/Evidence/PROTMOTIF-PROSITE-001-Evidence.md (PROSITE Pattern Matching: ConvertPrositeToRegex, FindMotifByProsite). THIRD Evidence over the SAME PROSITE→regex engine as PROTMOTIF-FIND-001 / PROTMOTIF-PATTERN-001 — NO new concept; enriched [[protein-motif-pattern-search]]. Distinct contributions folded in: the **`[G>]` C-terminus-inside-brackets** corner case (only PS00267 `F-[IVFY]-G-[LM]-M-[G>]`→`F[IVFY]G[LM]M(?:G|$)` and PS00539 `F-[GSTV]-P-R-L-[G>]`→`F[GSTV]PRL(?:G|$)`; residue-or-end-of-sequence → regex alternation, matched via both the G branch and the C-terminus branch, fails mid-sequence without G) added as a conversion-table row + corner-case bullet; mid-pattern period termination (`R-G-D.A-B-C`→`RGD`, §IV.E) sharpened; and a **real-protein positive control** — Human Transferrin P02787 (TRFE_HUMAN) × PS00001 `N-{P}-[ST]-{P}` → 2 N-glycosylation sites at 1-based 432–435 / 630–633 (0-based 431–434 / 629–632). Wrote concise source page [[protmotif-prosite-001-evidence]] recording the PROSITE User Manual PA-line spec + ScanProsite docs (extended syntax `-` omittable when unambiguous `MASKE`=`M-A-S-K-E`; greedy/overlap/include match modes), PS00001/00028 entries, Hulo 2007 + De Castro 2006, and the conversion/matching oracle tables. Updated hub [[algorithm-validation-evidence]] (frontmatter sources + body evidence-link list) and wiki/index.md (Sources-section entry). Bumped [[protein-motif-pattern-search]] source_commit→HEAD and added a relates_to test-unit-registry edge from protmotif-prosite-001-evidence. Contradictions: none — fully consistent with FIND-001 / PATTERN-001 (same engine, 0-based vs ScanProsite 1-based coordinate convention only). Follow-ups: none specific; remaining ProteinMotif units (transmembrane-helix etc.) still unenumerated.
   graph: +1 node, +1 typed edge (relates_to test-unit-registry)
- 2026-07-10 — ingest docs/Evidence/PROTMOTIF-SP-001-Evidence.md (Signal-peptide cleavage-site prediction, ProteinMotifFinder.PredictSignalPeptide). SAME METHOD as PROTMOTIF-DOMAIN-001 but a **REDESIGNED ALGORITHM** — the fabricated tripartite n/h/c + −1,−3 model (constants 0.95/0.825, NRegion/Probability fields, [0,1] score) was removed and REPLACED by the **von Heijne (1986) log-odds weight matrix** = EMBOSS 6.6.0 `sigcleave` (verified against current code `ProteinMotifFinder.PredictSignalPeptide` at src/…/ProteinMotifFinder.cs). Score = argmax over sites of `Σ ln(count/expect)` across positions −13..+2 (natural log; zero counts → `1.0e-10` at conserved cols −3/−1, else `1.0`); cleavage between −1/+1, CleavagePosition = 1-based mature start; IsLikelySignalPeptide ⇔ Score ≥ 3.5 (minWeight default 3.5); eukaryotic matrix (161 seqs) default, prokaryotic (36) via prokaryote:true. Worked oracle ACH2_DROME (UniProt P17644) → Score 13.739, mature start 42, window LLVLLLLCETVQA (re-derived exactly in Python). NO new concept — enriched [[protein-domain-and-signal-peptide-prediction]]: rewrote the signal-peptide section to the weight-matrix model, added a **Superseded note**, fixed the now-stale invariants/oracles + design-decisions + References + intro, added SP-001 to frontmatter sources + source_commit→HEAD. Wrote concise source page [[protmotif-sp-001-evidence]]. Updated hub [[algorithm-validation-evidence]] (frontmatter sources + body evidence-link list + source_commit→HEAD) and wiki/index.md (new Sources-section entry + a supersession clause appended to the DOMAIN-001 entry). Contradictions FLAGGED: the DOMAIN-001 evidence's tripartite signal-peptide description is now historical/superseded for this method (typed `supersedes` edge protmotif-sp-001-evidence → source:protmotif-domain-001-evidence, scoped to the signal-peptide method). Assumption: min input length = one full 15-aa window (< 15 → null). Follow-ups: none blocking; remaining ProteinMotif units (transmembrane-helix etc.) still unenumerated.
   graph: +1 node, +2 typed edges (relates_to test-unit-registry; supersedes source:protmotif-domain-001-evidence)
- 2026-07-10 — ingest docs/Evidence/PROTMOTIF-TM-001-Evidence.md (Transmembrane helix prediction, Kyte-Doolittle hydropathy sliding window; ProteinMotifFinder). GENUINELY DISTINCT ProteinMotif-family unit (hydrophobicity-based membrane-span detection) → NEW concept [[transmembrane-helix-prediction]]. Method: slide window W=19, score each window = arithmetic MEAN of per-residue Kyte-Doolittle (1982) hydropathy (`HydropathyScale` I 4.5…R −4.5, D/E/N/Q −3.5), emit contiguous runs with window mean ≥ threshold 1.6 as segments `[i₀, i₁+W−1]` with peak Score; profile length n−W+1; non-standard residues (X,B,Z,*) excluded from the mean; <W/null/empty→empty. Sources: Kyte & Doolittle 1982 (rank 1) + Davidson DGPB background (window 19 + threshold 1.6 + mean-windowing rule verbatim) + QIAGEN CLC + Davidson per-AA scores (20 scale values, matching exactly) + Biopython ProtParam protein_scale(edge=1.0)=mean / gravy + TM α-helix length ~18–21 residues / ~3–4 nm bilayer. Oracles D×10+L×20+D×10→one segment (5,34) peak 3.8 / D×40→none / L×19→(0,18,3.8). One assumption = segment End=lastPassingProfileIndex+windowSize−1 clamped (2026-06-16 off-by-one correction); no deviations. Cross-linked to [[intrinsic-disorder-prediction-top-idp]] (shares the Kyte-Doolittle scale via its CalculateHydropathy utility — no dedicated hydrophobicity concept exists in the wiki) and to the ProteinMotif siblings. Wrote source page [[protmotif-tm-001-evidence]]; created concept [[transmembrane-helix-prediction]]; updated hub [[algorithm-validation-evidence]] (frontmatter sources + body evidence-link list + source_commit→HEAD) and wiki/index.md (Sources-section + Concepts-section entries). Contradictions: none (both Davidson pages, QIAGEN, Biopython agree on scale + mean windowing). Follow-ups: no dedicated Kyte-Doolittle hydrophobicity concept yet — the scale is now referenced by both this unit and the disorder anchor; could be extracted if a third consumer appears.
   graph: +2 nodes, +2 typed edges (relates_to test-unit-registry; relates_to intrinsic-disorder-prediction-top-idp)
- 2026-07-10 — ingest docs/Evidence/QUALITY-PHRED-001-Evidence.md (Phred score handling — ParseQualityString / ToQualityString / ConvertEncoding, Phred+33 ↔ Phred+64). OVERLAPS the existing [[phred-quality-encoding]] concept (same Q=−10·log₁₀(P), same two ASCII offsets, same boundary chars) — despite the "QUALITY" family name this unit is about the ENCODING, not Q20/Q30 statistics. **NO new concept** — enriched [[phred-quality-encoding]]: added the **primary-literature anchor Cock et al. 2010** (de-facto FASTQ spec, PMC2847217) and a new "Converting between the two offsets" section (Phred score invariant → pure ±31 re-offset; **Phred+64→Phred+33 always safe** Q0–62⊆0–93; **Phred+33→Phred+64 overflows** for Q>62 → ArgumentOutOfRangeException; below-offset byte → negative Q = malformed; Solexa odds-ratio score lossy/out of scope); added QUALITY-PHRED-001 to frontmatter sources + source_commit→HEAD + cross-link to the new source. Wrote concise source page [[quality-phred-001-evidence]] with the worked oracles (`!`/`5`/`?`/`I`/`~`→0/20/30/40/93, `@h~`→0/40/62, conversion `@h~`→`!I_`, `!I`→`@h`). Updated hub [[algorithm-validation-evidence]] (frontmatter sources + body evidence-link list + source_commit→HEAD) and wiki/index.md (Sources-section entry + enriched the [[phred-quality-encoding]] Concepts-section entry with conversion + Cock anchor). Two API-shape ASSUMPTIONs (malformed byte + Q>62 overflow → ArgumentOutOfRangeException; range bounds themselves source-backed). Contradictions: none — fully consistent with PARSE-FASTQ-001 (Cock et al. 2010 vs its 2009 citation are the same paper; offsets/ranges/boundary chars agree). Follow-ups: none — Q20/Q30 / mean-quality *statistics* are a separate surface (noted in parse-fastq must-test), not covered by this unit.
   graph: +1 node, +1 typed edge (relates_to test-unit-registry)
- 2026-07-10 — ingest docs/Evidence/QUALITY-STATS-001-Evidence.md (FASTQ quality statistics — Q20/Q30 fractions, mean, median, min/max, population variance/std dev; QUALITY family, sibling of QUALITY-PHRED-001). GENUINELY DISTINCT surface from Phred encoding — the [[phred-quality-encoding]] concept explicitly flagged Q20/Q30 statistics as a SEPARATE surface, so this ingest **creates a new concept** [[fastq-quality-statistics]] that **depends on** (consumes decoded scores from) [[phred-quality-encoding]]. Statistics run over DECODED Phred scores → encoding-independent (Phred+64 same-scores → identical stats; decode is QUALITY-PHRED-001's contract, cited-only not re-tested). Contract: mean = arithmetic mean of scores (a mean over log-scaled values, NOT error-probability-averaged); median odd=middle / **even = mean of the two central order statistics**; min/max; **population** variance/σ (`(1/N)Σ(Qᵢ−μ)²`, **÷N not N−1** — quality string is the complete population); **`% ≥ Q20` / `% ≥ Q30`** with **inclusive `≥`** thresholds (Illumina: `% ≥ Q30` the NGS benchmark); `CalculateQ30Percentage` == `CalculateStatistics(...).PercentAboveQ30`. Sources: Illumina Sequencing Quality Scores (rank 2) + Newcastle Univ. ASK (population σ, rank 1) + Math is Fun (even-count median) + Wikipedia/Ewing & Green 1998 (Phred formula provenance) + Cock et al. 2010 (decode, cited-only). Oracles `5?I`→20/30/40 mean 30.0/median 30/min-max 20-40/var 200/3≈66.6667/σ≈8.16497/%≥Q20 100/%≥Q30≈66.67, even `5II?`→median (30+40)/2=35.0 mean 32.5, single `I`(Q40)→mean=median=min=max=40 σ=0 %≥Q20/Q30 100. Corner cases even/odd median branch, single-element σ=0, empty/null → zeroed `QualityStatistics` (TotalBases=0). Wrote concept [[fastq-quality-statistics]] + source page [[quality-stats-001-evidence]]; updated hub [[algorithm-validation-evidence]] (frontmatter sources + body evidence-link list) and wiki/index.md (Sources-section + Concepts-section entries); cross-linked [[phred-quality-encoding]] to the new concept. One API-shape ASSUMPTION (empty→zeroed not throw; no numeric value invented). Contradictions: none — fully consistent with QUALITY-PHRED-001 / PARSE-FASTQ-001 (same decode, same Q20/Q30 inclusive thresholds). Follow-ups: none.
   graph: +2 nodes, +2 typed edges (depends_on phred-quality-encoding; relates_to test-unit-registry)
- 2026-07-10 — ingest docs/Evidence/REP-STR-001-Evidence.md (Microsatellite / Short Tandem Repeat (STR) detection — perfect default `FindMicrosatellites` + opt-in approximate/imperfect/interrupted `FindApproximateTandemRepeats` + `ComputeBernoulliStatistics`, Benson Tandem Repeats Finder 1999 model). OVERLAPS the repeats family anchor [[repetitive-element-detection]] (tandem sub-problem, microsatellite/STR by unit length) — this unit is the concrete APPROXIMATE detector that CLOSES the "exact-copies-only" Framework/Simplified limitation the concept previously documented. **NO new concept** — enriched [[repetitive-element-detection]]: softened the "both are exact" paragraph (default paths exact; opt-in approximate path closes the gap), added an *Approximate STR detection (Benson TRF model)* subsection (seven TRF statistics, wraparound-DP alignment, majority-rule consensus so ConsensusSize==Period, weights `+2/−7/−7`, Minscore default 50; Bernoulli adjacent-copy PM/PI defaults .80/.10 distinct from consensus percent-matches, E[heads]=PM·d reproduced, R(d,k,pM)/W(d,pI) k-tuple seeding NOT reproduced = deterministic exhaustive (start,period) scan residual), added REP-STR-001-Evidence to frontmatter sources + source_commit→HEAD + a relates_to test-unit-registry graph edge sourced from rep-str-001-evidence. Wrote concise source page [[rep-str-001-evidence]] with the worked oracles (`CACACACACA`→CA×5 score 20 100%, `CAGCAGCAGTAGCAGCAG`→CAG×6 score 27 94.4% vs perfect fragments to CAG×3, 29-bp single-deletion score 51 clears gate, Bernoulli adjacent PM 13/15 / 8/10 / 0.80-on-threshold / 0.00) and the three assumptions. Updated hub [[algorithm-validation-evidence]] (frontmatter sources + body evidence-link list + source_commit→HEAD) and wiki/index.md (Sources-section entry only). Contradictions: none — Wikipedia + Benson 1999 + TRF tool docs + reference impl agree; perfect and approximate detectors are complementary (approximate is opt-in). Follow-ups: the sibling exact units GENOMIC-TANDEM-001 / ANNOT-REPEAT-001 already ingested; the k-tuple probabilistic seeding + non-redistributable simulation-table percentiles remain the documented genome-scale research-grade residual.
   graph: +1 node, +1 typed edge (relates_to test-unit-registry)
- 2026-07-10 — ingest docs/Evidence/RESTR-FILTER-001-Evidence.md (Restriction Enzyme Filtering — the FIRST RESTR-* / MolTools reagent-selection unit: `GetBluntCutters()` / `GetStickyCutters()` end-type filters + `GetEnzymesByCutLength(min,max)` recognition-length range + single-length overload). GENUINELY DISTINCT new domain (restriction enzymes) — no prior digest/enzyme concept existed in the wiki → NEW concept [[restriction-enzyme-filtering]]. Two filter axes: (1) end type is a **total, disjoint partition** — every Type II end is blunt (center cut, both strands terminate in a base pair) or sticky (staggered cut, 5'/3' overhang), no third category so blunt ∪ sticky = full library & disjoint (blunt-blunt always compatible); (2) recognition-site length over the **inclusive** `[min,max]` interval, undivided Type II sites canonically 4–8 nt. Blunt SmaI/EcoRV/AluI/HaeIII vs sticky EcoRI(5')/KpnI(3')/PstI(3')/NotI/TaqI; the **interrupted palindrome** SfiI (`GGCCNNNN^NGGCC`, 13 nt) is sticky but correctly excluded by `[4,8]` (undivided sites only). One API-shape ASSUMPTION (range bounds inclusive; recognition-length values themselves source-backed). Sources: Wikipedia *Sticky and blunt ends* + *Restriction enzyme* (Type II 4–8 nt undivided palindromes, center→blunt/staggered→sticky, EcoRI/SmaI/KpnI/PstI worked cuts) + *List of restriction enzyme cutting sites* (4/6/8-bp categories) + NEB/REBASE (KpnI 3' overhang, EcoRI 5' overhang) + PMC/REBASE (SfiI interrupted palindrome). Wrote concept [[restriction-enzyme-filtering]] + source page [[restr-filter-001-evidence]]; updated hub [[algorithm-validation-evidence]] (frontmatter sources + body evidence-link list + source_commit→HEAD) and wiki/index.md (Sources-section + Concepts-section entries). Contradictions: none — all sources agree end type is a blunt-or-overhang dichotomy and the undivided-site length range is 4–8 nt. Follow-ups: the complementary RESTR units (cut-site finding on a target, digest simulation, compatible-overhang/ligation planning) are not yet ingested.
   graph: +2 nodes, +1 typed edge (relates_to test-unit-registry)
- 2026-07-10 — ingest docs/Evidence/RNA-DOTBRACKET-001-Evidence.md (Dot-Bracket / extended WUSS notation — parse & validate a structure string; the notation/representation layer of the RNA secondary-structure family, `RnaSecondaryStructure.ParseDotBracket` / `ValidateDotBracket`). GENUINELY DISTINCT surface — no notation/dot-bracket concept existed (the base-pairing chemistry [[rna-base-pairing]] and hairpin/MFE folding [[pre-mirna-hairpin-detection]] are the neighbours, both of which *emit/consume* dot-bracket but neither *is* the notation-parse layer) → NEW concept [[rna-dot-bracket-notation]]. Core algorithm: **one balanced-bracket stack per family** — `()`/`<>`/`{}`/`[]` + uppercase(5' open)/lowercase(3' close) letter pairs, each an **independent pairing system** (ViennaRNA + Infernal `vrna_db_from_WUSS()`: any matched pair = a base pair, exact symbol has no meaning if partners match; flatten treats letter-pair pseudoknots as unpaired); a shared stack would mis-pair `([)]`. Validate ⟺ every family's stack empty at end & never underflows, closer must match a same-family opener → **crossing families (pseudoknots) valid**, **mismatched families `(]` invalid**; non-bracket WUSS symbols `-`/`,`/`:`/`.` are single-stranded (Rfam). Oracles: parse `((((....))))`→(0,11),(1,10),(2,9),(3,8), `([)]`→`(`:(0,2)+`[`:(1,3), `<<<<[[[[....>>>>]]]]`≡`((((AAAA....))))aaaa`→two crossing 4-bp helices; validate `(((...)))`/`(([[]]))`/`([)]`→true, `(((...)` / `...)` / `)(` / `(]`→false. Sources: ViennaRNA RNA-Structure-Notations + Dot-Bracket (rank 3) + WUSS/`vrna_db_from_WUSS()` + Infernal Nawrocki & Eddy 2013 (rank 3) + Rfam glossary (rank 5); all agree, no contradictions. Two API-contract assumptions (malformed → best-effort parse dropping unmatched closers, gate with `ValidateDotBracket`; empty/null → valid pair-free) — neither invents a numeric value. Wrote concept [[rna-dot-bracket-notation]] + source page [[rna-dotbracket-001-evidence]]; updated hub [[algorithm-validation-evidence]] (frontmatter sources + body evidence-link list + own-concept paragraph + source_commit→HEAD) and wiki/index.md (Sources-section + Concepts-section entries); cross-linked [[rna-base-pairing]] and [[pre-mirna-hairpin-detection]] to the new notation concept. Contradictions: none. Follow-ups: the sibling RNA-STRUCT-* folding units (MFE / stem-loop enumeration / pseudoknot detection / base-pair classification / loop-energy terms) are not yet individually ingested.
   graph: +2 nodes, +1 typed edge (relates_to test-unit-registry)
- 2026-07-10 — ingest docs/Evidence/RNA-ENERGY-001-Evidence.md (Free Energy Calculation — the thermodynamic/energy layer of the RNA secondary-structure family, `RnaSecondaryStructure.CalculateStackingEnergy`/`CalculateStemEnergy`/`CalculateHairpinLoopEnergy`/`CalculateMinimumFreeEnergy`). GENUINELY DISTINCT layer — the base-pairing chemistry [[rna-base-pairing]] and the [[rna-dot-bracket-notation]] notation are neighbours, but neither IS the Turner-2004 free-energy model → NEW concept [[rna-free-energy-turner-model]]. Core: NN model `ΔG°total = init + Σstacking + Σloops` over the NNDB Turner04 tables at 37 °C. WC stacking all negative (GC-rich most stable `GC/CG` −3.42); G-U wobble variable with **two POSITIVE** (`UG/GU` +0.30, `GU/UG` +1.29), note-a `GG/UU`=−0.5, note-b special 3-stack `5'GGUC/3'CUGG`=−4.12 (vs −1.77); hairpin-loop init positive & **non-monotonic** (3→5.4/4→5.6/6→5.4/9→6.4, Jacobson-Stockmayer beyond 9); special UNCG/GNRA tri/tetra/hexaloop total energies replace the model; all-C loop penalty (3-nt +1.5, >3-nt 0.3n+1.6); terminal mismatch (96); +0.45 per-AU/GU-end; single base pair / empty / poly-A → ΔG°=0 (stacking needs ≥2 adjacent pairs). Oracles GC 3-bp stem −5.78, NNDB hairpin example-1 `GGGAUAAAUCCC` −3.42, GGUC/CUGG −4.12. Wrote concept [[rna-free-energy-turner-model]] + source page [[rna-energy-001-evidence]]; updated hub [[algorithm-validation-evidence]] (frontmatter sources + body evidence-link list + own-concept paragraph + source_commit→HEAD) and wiki/index.md (Sources-section + Concepts-section entries); cross-linked the folding-family neighbours [[pre-mirna-hairpin-detection]] (its Turner ΔG `FreeEnergy`) and [[rna-dot-bracket-notation]] (its "does not assign energies" scope note) to the new concept. Graph: concept-to-concept edges on the new concept page — relates_to test-unit-registry + relates_to rna-base-pairing (stacking energies assigned per base-pair stack). Contradictions: none — all parameter sets are exact NNDB Turner04 matches; the three recorded items (37 °C standard state, 2-dp precision, unknown stacks→0.0) are DEFINED CONDITIONS, not assumptions. Follow-ups: the sibling RNA-STRUCT-001 MFE folder (which consumes these terms) is still not individually ingested.
   graph: +2 nodes, +2 typed edges (relates_to test-unit-registry; relates_to rna-base-pairing)
- 2026-07-10 — ingest docs/Evidence/RNA-HAIRPIN-001-Evidence.md (Hairpin Loop and Stem Free-Energy Calculation, Turner 2004 NN model — `RnaSecondaryStructure.CalculateHairpinLoopEnergy` + `CalculateStemEnergy`). CONTEXT hypothesis (generic stem-loop *enumerator*) did NOT match the actual content: this unit is the *energy* of a hairpin's loop + stem, i.e. a focused deep-dive on the same thermodynamic layer already synthesized by [[rna-free-energy-turner-model]] (RNA-ENERGY-001, whose concept page already exposes these two methods). It is NOT enumeration and NOT miRNA-specific, so it is DISTINCT from the miRNA hairpin *finder* [[pre-mirna-hairpin-detection]] (which merely sums a downstream `FreeEnergy`). **NO new concept** — REUSED/enriched [[rna-free-energy-turner-model]]: added RNA-HAIRPIN-001 to its frontmatter sources (+source_commit→HEAD), updated the intro to name both sibling records, and added the two NNDB worked hairpin oracles (Example 1 closing A-U 6-nt loop → loop +4.6 + helix −6.01 = −1.4; Example 2 5-nt loop G…G with GG first-mismatch bonus → +4.1 + −6.01 = −1.9; 3-nt loops get no first-mismatch term; P pairs → P−1 stacks). Source: NNDB Turner 2004 hairpin pages (retrieved via Wayback — live server down) + Mathews et al. 2004 PNAS 101:7287. Key params: first-mismatch bonuses UU/GA −0.9, GG −0.8; special-GU-closure −2.2 (**G-U only, not U-G** — documented asymmetry); all-C penalty 3-nt +1.5 / >3-nt 0.3n+1.6; special tri/tetra/hexaloop totals override the model; loops <3 nt prohibited; +0.45/AU-end. Wrote concise source page [[rna-hairpin-001-evidence]] (cross-linked to [[rna-free-energy-turner-model]], [[pre-mirna-hairpin-detection]], [[rna-base-pairing]], [[rna-dot-bracket-notation]]). Updated hub [[algorithm-validation-evidence]] (frontmatter sources + body evidence-link list) and wiki/index.md (Sources-section entry only, "no new concept"). Contradictions: none — all parameters exact NNDB Turner 2004 match; only recorded item is a 2-dp intermediate-rounding display choice (tests assert `.Within(1e-9)`). No new typed edges (source page reuses the existing concept's edges). Follow-ups: the RNA-STRUCT-001 MFE folder that consumes these terms is still not individually ingested.
   graph: +1 node, +0 typed edges
- 2026-07-10 — ingest docs/Evidence/RNA-INVERT-001-Evidence.md (RNA Inverted Repeats / potential stem regions — antiparallel reverse-complement arms that form a stem-loop, `RnaSecondaryStructure`, RNA secondary-structure family). CONTEXT check: the inverted-repeat model is ALREADY synthesized on the repeats-family anchor [[repetitive-element-detection]] (its Inverted-repeats section, from the SAME IUPACpal source), and the RNA antiparallel complement is [[rna-base-pairing]]; a looped IR IS a stem-loop, the object of [[pre-mirna-hairpin-detection]] / [[rna-dot-bracket-notation]]. **NO new concept** — genuinely NOT distinct, so REUSED/cross-linked those existing concepts. Sources: IUPACpal (Alamro 2021 PMC7866733 — `W G W̄ᴿ` gapped model + k-Hamming mismatch, perfect=k=0, RNA complement A⟷U/C⟷G) + Wikipedia/Ussery 2008 (`5'---TTACGnnnnnnCGTAA---3'`, zero-gap ⇒ palindrome) + EMBOSS einverted (IR = stem-loop = local alignment of a sequence vs its reverse complement). Oracles `UUACGAAAAAACGUAA` (arm `UUACG` 0–4 / loop 5–10 / arm `CGUAA` 11–15) and palindromic `GGCCAAAGGCC` (left 0–3 / right 7–10 / len 4). Scope restriction: perfect ungapped k=0 arms only (einverted scored mismatch/gap DP Not Implemented), loop bounds via minSpacing/maxSpacing, maximal-arm non-overlapping greedy reporting. Wrote source page [[rna-invert-001-evidence]] (cross-linked to the four reused concepts + hub + registry). Updated hub [[algorithm-validation-evidence]] (frontmatter sources + body evidence-link list) and [[repetitive-element-detection]] (Inverted-repeats cross-ref to the RNA sibling + frontmatter source + source_commit→HEAD + typed edge). Updated wiki/index.md (Sources-section entry only, "no new concept"). Contradictions: none — IUPACpal, Ussery/Wikipedia, and EMBOSS einverted agree on the reverse-complement-arms definition. Follow-ups: the RNA-STRUCT-001 MFE folder is still not individually ingested.
   graph: +1 node, +1 typed edge (relates_to test-unit-registry on repetitive-element-detection)
- 2026-07-10 — ingest docs/Evidence/RNA-MFE-001-Evidence.md (Minimum Free Energy (MFE) RNA secondary-structure prediction — the Zuker–Stiegler dynamic-programming folder; the folding/search layer of the RNA secondary-structure family, `RnaSecondaryStructure.CalculateMinimumFreeEnergy` / `PredictStructure`). GENUINELY DISTINCT layer — prior RNA ingests repeatedly flagged this as the not-yet-ingested MFE folder that CONSUMES the Turner terms; the energy layer [[rna-free-energy-turner-model]] only *scores* a given structure and the [[rna-dot-bracket-notation]] notation only *parses* it, but neither IS the folding/search DP → NEW concept [[rna-minimum-free-energy-folding]]. Core: Zuker–Stiegler (1981) loop decomposition (hairpin / stacking / bulge-interior / multibranch); Ward 2017 DP matrices C(i,j)=min(hairpin, interior/bulge over an inner pair, multiloop) + multiloop M/M1 + exterior F; **standard affine multiloop model → O(n³) time / O(n²) space** (logarithmic would be O(n⁴); ViennaRNA/Lorenz 2011 confirm O(n³) & the Zuker–Stiegler derivation). Oracles `CalculateMinimumFreeEnergy("CACAAAAAAAUGUG")`=−1.41 (NNDB Example 1, `PredictStructure`→`((((......))))`), `CACAGAAAGUGUG`=−1.91 (Example 2, GG first mismatch). Invariants INV-01 MFE ≤ 0 (empty open-chain always in search set), INV-02 suffix-monotone `MFE(s) ≤ MFE(prefix)`, INV-03 optimized DP == classic O(n³) baseline; empty/null / homopolymer `AAAAAAAA` / sequence `< minLoopSize+2` (`GCGC`) → 0 (hairpin loop ≥ 3 nt); intramolecular ⇒ no helix-init constant. NAMING RECONCILIATION: prior RNA pages referred to this folder generically as *RNA-STRUCT-001* (the id the pre-miRNA `AssessHairpinByMfe` path cites) — the Evidence artifact records it under its own id **RNA-MFE-001**; both denote this one MFE folder (left prior *RNA-STRUCT-001* mentions in place, added the concept link where natural). Wrote concept [[rna-minimum-free-energy-folding]] + source page [[rna-mfe-001-evidence]]; updated hub [[algorithm-validation-evidence]] (frontmatter sources + body evidence-link list + concept-list entry + source_commit→HEAD) and wiki/index.md (Sources-section + Concepts-section entries); cross-linked the family — [[rna-free-energy-turner-model]] ("MFE folder's job") and [[rna-dot-bracket-notation]] ("MFE folder reads its hairpin") now link the new folder concept. Graph: concept-to-concept edges on the new concept page — relates_to test-unit-registry + **depends_on rna-free-energy-turner-model** (consumes the Turner terms) + relates_to rna-dot-bracket-notation (produces dot-bracket output). Contradictions: none — Zuker & Stiegler 1981, Lorenz 2011, Ward 2017, NNDB Turner 2004 / Mathews 2004 agree on the DP decomposition, O(n³)/O(n²) affine complexity, and worked-example energies. Two documented simplifications (multiloop per-unpaired `c=0` with `a=9.25`/helix `c=−0.63`; 2-dp rounding, tests `.Within(1e-9)` + `Round(mfe,1)`==NNDB). Follow-ups: the RNA secondary-structure family's remaining folding surfaces (stem-loop enumeration, pseudoknot detection, base-pair classification) are still not individually ingested.
   graph: +2 nodes, +3 typed edges (relates_to test-unit-registry; depends_on rna-free-energy-turner-model; relates_to rna-dot-bracket-notation)
- 2026-07-10 — ingest docs/Evidence/RNA-PAIR-001-Evidence.md (RNA Base Pairing — `RnaSecondaryStructure.CanPair` / `GetBasePairType` / `GetComplement`, the RNA-secondary-structure family's own base-pairing primitive). OVERLAPS the shared pairing rule already synthesized on [[rna-base-pairing]] (the MIRNA-PAIR-001 `MiRnaAnalyzer` sibling) — same {A-U, G-C} + G-U wobble rule, same `T`→`U` normalisation, same case-insensitivity — so **NO new concept**: enriched [[rna-base-pairing]] with a "RNA-secondary-structure family's own copy (RNA-PAIR-001)" subsection documenting the two shape differences (the **typed `GetBasePairType` classifier** returning `WatsonCrick`/`Wobble`/`null` as a first-class value, where the miRNA surface splits WC-vs-wobble across `CanPair`+`IsWobblePair`; and the single-base `GetComplement`, base-level counterpart of `GetReverseComplement`), added RNA-PAIR-001-Evidence to its frontmatter sources + source_commit→HEAD + a second relates_to test-unit-registry graph edge sourced from rna-pair-001-evidence. Sources: Crick 1966 *J Mol Biol* 19:548 (wobble hypothesis — G-U the only standard wobble, distinct from WC) + Wikipedia Base pair (A•U 2 H-bonds / G•C 3 H-bonds, reciprocal) / Wobble base pair + IUPAC-IUB 1970 (complement table) + Biopython `complement_rna("CGAUT")`→`"GCUAA"`. Rule {A-U,U-A,G-C,C-G}=WatsonCrick / {G-U,U-G}=Wobble / else false+null; complement A→U/U→A/G→C/C→G/**T→A**(DNA T=U)/N→N/R→Y; symmetry `f(x,y)==f(y,x)`; non-alphabet→false/null no exception. One non-correctness normalization (case-insensitive upper-casing). Wrote source page [[rna-pair-001-evidence]] (cross-linked to [[rna-base-pairing]], [[rna-free-energy-turner-model]], [[rna-dot-bracket-notation]], hub + registry). Updated hub [[algorithm-validation-evidence]] (frontmatter sources + body evidence-link list + source_commit→HEAD) and wiki/index.md (Sources-section entry only, "no new concept"). Contradictions: none — Crick 1966, Wikipedia, IUPAC-IUB 1970, Biopython agree; the miRNA and RNA-structure copies of the pairing rule are identical chemistry. Follow-ups: none new — the remaining RNA-structure surfaces (stem-loop enumeration, pseudoknot detection, loop-energy terms deep-dives) already tracked by prior RNA ingests.
   graph: +1 node, +1 typed edge (relates_to test-unit-registry on rna-base-pairing)
- 2026-07-10 — ingest docs/Evidence/RNA-PARTITION-001-Evidence.md (RNA Partition Function (McCaskill) and Boltzmann Structure Probability — the probabilistic/ensemble layer of the RNA secondary-structure family, `RnaSecondaryStructure`). GENUINELY DISTINCT layer — the Boltzmann-weighted **ensemble** counterpart of the single-optimum MFE folder [[rna-minimum-free-energy-folding]]: instead of one lowest-energy fold it computes the equilibrium partition function `Z` over ALL pseudoknot-free structures + per-base-pair binding probabilities + the Boltzmann probability of a given structure → NEW concept [[rna-partition-function-mccaskill]]. Core: McCaskill 1990 O(n³) time / O(n²) space recursion; inside `Q_ij = Q_{i,j-1} + Σ_{i≤k<j-m} Q_{i,k-1}·Q^b_{kj}`, base `Q_ij=1` for i≥j−m, total `Z = Q_{1n}` (disjoint/unambiguous decomposition); Boltzmann `Pr[P|S]=Z⁻¹exp(−βE(P))`, `p(s)=e^(−βE(s))/Z` (ViennaRNA), RT=0.61626805 at 37 °C=310.15 K. KEY CORRECTION (Evidence 2026-06-16): base-pair probability requires the **outside** recursion `p_kl = Q^b_kl·O_kl/Z`, `O_kl = Q_{1,k-1}·Q_{l+1,n} + Σ w·Q_{i+1,k-1}·Q_{l+1,j-1}·O_ij`; the external-only term is WRONG for nestable pairs (`GGGAAACCC` P(2,6)=6/20 not 1/20; `GGGGCCCC` P(1,5)=3/16 not 1/16) — verified to 3.3e-16 vs Boltzmann brute force (an earlier "external suffices" claim matched a since-fixed impl bug). Oracles (E_bp=0 ⇒ Z counts admissible structures, two independent derivations) `AAAA`→1, `GC`→1, `GGGGCCCC`→16, `GGGAAACCC`→20; invariants Z ≥ 1 (empty structure weight 1), P(i,j) ∈ [0,1] + symmetric, per-base pairing sum ≤ 1 (300 random seqs max 0.983), monotone in E_bp; WC {A-U,G-C}+GU pairing only, min-loop m forbids j−k ≤ m. Wrote concept [[rna-partition-function-mccaskill]] + source page [[rna-partition-001-evidence]]; updated hub [[algorithm-validation-evidence]] (frontmatter sources + source_commit→HEAD + body evidence-link list + concept-list entry) and wiki/index.md (Sources-section + Concepts-section entries); cross-linked the RNA family — [[rna-minimum-free-energy-folding]] now names the partition function as its Boltzmann-weighted ensemble counterpart (+RNA-PARTITION source + source_commit→HEAD). Graph: concept-to-concept edges on the new concept page — relates_to test-unit-registry + **alternative_to rna-minimum-free-energy-folding** (ensemble vs single-optimum counterpart) + depends_on rna-free-energy-turner-model (Boltzmann-weights the Turner energies; medium confidence — Seqeron uses a documented simplified fixed-per-pair E_bp model). Contradictions: none — McCaskill 1990, MIT 18.417 slides, Freiburg tool, ViennaRNA agree on the recurrence, Boltzmann form, and O(n³)/O(n²) complexity. One documented assumption: simplified per-pair `E_bp` energy model vs full Turner NN (energy model only; recurrence + probabilities + invariants conformant with McCaskill 1990). Follow-ups: exact Turner-parameter ensemble energies out of scope; remaining RNA folding surfaces (pseudoknot detection, stem-loop enumeration, base-pair classification) still not individually ingested.
   graph: +2 nodes, +3 typed edges (relates_to test-unit-registry; alternative_to rna-minimum-free-energy-folding; depends_on rna-free-energy-turner-model)
- 2026-07-10 — **lint pass** (structural + staleness + coverage + graph). Started: 11 broken wikilinks (all false positives — dot-bracket `[[[[....>>>>]]]]` and matrix `[[1.0]]` inside inline code, plus 2 valid intra-page `[[#anchor]]` links in primer-dimer-thermodynamics-tm), 1 orphan (backlog, a `type: index` meta page already linked from index.md), 1 soft-cap oversize (backlog 453 ln). Staleness clean, graph lint clean, coverage = 223 uncovered = the tracked ingest campaign in [[backlog]] (not a gap). Fixes (approved): (1) hardened `scripts/wiki_lint.py` — strip fenced/inline code spans and skip `#`-anchor targets before wikilink extraction (clears all 11 false positives, no content edits), and exempt `type: index` pages from the orphan check; (2) refreshed the stale [[backlog]] count in index.md Meta (54/191 → 74/171 to match backlog.md). Post-fix: clean except the intentionally-kept backlog soft-cap. Semantic pass: RNA secondary-structure cluster (base-pairing/dot-bracket/Turner-energy/MFE/partition) internally consistent, no contradictions; recurring gap = the not-yet-ingested RNA folding surfaces (pseudoknot detection, stem-loop enumeration, base-pair classification), tracked in [[backlog]].
- 2026-07-10 — ingest docs/Evidence/RNA-PKPREDICT-001-Evidence.md (Pseudoknot Structure Prediction — canonical H-type, pknotsRG class; the **crossing-helix layer** of the RNA secondary-structure family, `RnaSecondaryStructure`). GENUINELY DISTINCT — predicts the optimal fold that may contain a single **pseudoknot** (two helices whose pairs **cross**, `i<k<j<l` Antczak 2018), the one feature the nested MFE folder [[rna-minimum-free-energy-folding]] and McCaskill ensemble [[rna-partition-function-mccaskill]] are definitionally blind to → NEW concept [[rna-pseudoknot-prediction]]. pknotsRG (Reeder & Giegerich 2004 *BMC Bioinformatics* 5:104, PMC514697): "two crossing helices with three intervening loops", grammar `a~~~u~~~b~~~v~~~a'~~~w~~~b'`, H-type 5'→3' stem1-5'→loop1→stem2-5'→loop2→stem1-3'→loop3→stem2-3', **O(n⁴)/O(n²)**, two-layer dot-bracket `((((..[[[[..))))..]]]]`. Canonization rules bound the search (equal-length bulge-free helices / maximal extent / fixed overlap boundary). Energy = Turner NN stacking on BOTH helices (same model as nested, no extra per-pair penalty — pknotsRG `Energy.lhs`) + penalties **initiation 9.0** (anti-spurious-knot gate) / **unpaired loop nt 0.3** / **base pair inside knot 0.0** kcal/mol. Oracles designed H-type `GGGGAACCCCAACCCCAAGGGG`→`HasPseudoknot==true` two crossing 4-bp helices (0,15)…(3,12)+(6,21)…(9,18); plain `GGGGAAAACCCC`→no knot = MFE `((((....))))`; BWYV `GGCGCGGCACCGUCCGCGGAACAAACGG` (PDB 437D, Su 1999) NOT recovered (tertiary-stabilized triplex/ion coordination outside NN model — documented limit of all NN-only pseudoknot predictors). Invariants MFE fallback `FreeEnergy ≤ CalculateMfeStructure().FreeEnergy`, no spurious knot, each position paired ≤1× + ≥1 genuine crossing, empty/null/too-short→empty pair-free ΔG 0. One documented scope note: PARTIAL pknotsRG coverage — single canonical H-type only; recursively-nested / over-arching / multiple knots NOT implemented; loops u/v/w fold with the existing MFE. Wrote concept [[rna-pseudoknot-prediction]] + source page [[rna-pkpredict-001-evidence]]; updated hub [[algorithm-validation-evidence]] (frontmatter sources + source_commit→HEAD + body evidence-link list) and wiki/index.md (Sources-section + Concepts-section entries); cross-linked the RNA family — [[rna-minimum-free-energy-folding]] now names the pseudoknot predictor as its crossing-helix extension, [[rna-dot-bracket-notation]] links the two-layer output. Graph: concept-to-concept edges on the new concept page — relates_to test-unit-registry + depends_on rna-free-energy-turner-model (both helices scored with Turner terms) + depends_on rna-minimum-free-energy-folding (MFE fallback baseline + folds internal loops) + relates_to rna-dot-bracket-notation (two-layer output). Contradictions: none — Reeder & Giegerich 2004, pknotsRG `Energy.lhs`, Wikipedia/Rivas & Eddy H-type geometry, PDB 437D/Su 1999, Antczak 2018 mutually consistent. Follow-ups: pseudoknot detection now ingested (closes a recurring RNA-family gap); remaining RNA surfaces (stem-loop enumeration, base-pair classification) still not individually ingested; full pknotsRG recursive/multiple-knot grammar is a documented library limitation, not a wiki gap.
   graph: +1 node, +4 typed edges (relates_to test-unit-registry; depends_on rna-free-energy-turner-model; depends_on rna-minimum-free-energy-folding; relates_to rna-dot-bracket-notation)
- 2026-07-10 — ingest docs/Evidence/RNA-PKRECURSIVE-001-Evidence.md (Recursive pknotsRG pseudoknot prediction — nested / multiple / over-arching H-type knots; the recursive-grammar extension of the single-knot RNA-PKPREDICT-001, `RnaSecondaryStructure`). NOT a new concept — it is the SAME pknotsRG class (same Turner-NN energy on both helices, same penalties 9.0/0.3/0.0, same canonization rules) that fills exactly the PARTIAL-coverage gap the single-knot unit recorded, so ENRICHED the existing [[rna-pseudoknot-prediction]] with a new §6 (recursive extension) rather than creating a page. The delta: the three loops u/v/w now fold by the SAME recursive folder (a loop may contain a further knot), the top level CHAINS multiple knots, and an enclosing helix may OVER-ARCH a knot in its loop. Sources (all re-used from RNA-PKPREDICT-001): Reeder & Giegerich 2004 *BMC Bioinformatics* 5:104 (loops "fold internally … including simple recursive pseudoknots", O(n⁴)/O(n²), canonization 8→4 boundaries vs Rivas & Eddy O(n⁶)/O(n⁴)) + Reeder, Steffen & Giegerich 2007 *NAR* 35:W320 (per-interval COMPETITION with unknotted foldings — the whole-sequence mechanism enabling multiple/nested knots) + pknotsRG `Energy.lhs` (verbatim 9.0/0.3/0.0) + Antczak 2018 (crossing i<k<j<l). Constructed fully-derivable oracles: over-arching `AAAAAAAAGGGGAACCCCAACCCCAAGGGGUUUUUUUU` (38 nt)→`((((((((((((..[[[[..))))..]]]]))))))))` ΔG −14.37 (single-knot/MFE both −13.05, no combined structure); two-knot 80-nt→two crossing knots (crossing-count 32) ΔG −28.74 (single/MFE −27.14, none); plain `GGGGAAAACCCC`→no knot = MFE −5.28; single-knot parity `GGGGAACCCCAACCCCAAGGGG`→identical −8.76. Invariants recursive ΔG ≤ MFE (0 violations on a 150-seq random sweep, seed 20260623), no spurious knots, each position paired ≤1× + ≥1 crossing, empty/null/too-short→empty pair-free. Excluded (verbatim): triple-crossing helices, kissing hairpins, bulged/unequal-length helices. Two scope notes: PARTIAL recursion (realizes the recursive CLASS via a maximal-extent helix start/end scan, not bit-identical to the reference 4-boundary ADP parser) + two-simultaneous-knot cases are ENGINEERED (isolated A·U clamps — two strong G·C knots are the genuine MFE only when the cross-region nested alternative is suppressed). Wrote source page [[rna-pkrecursive-001-evidence]]; enriched concept [[rna-pseudoknot-prediction]] (intro now names both records, §5 PARTIAL note reframed as split-across-two-units, new §6, +RNA-PKRECURSIVE source + source_commit→HEAD); updated wiki/index.md (Sources-section entry, "enriches — no new concept"). Graph: no new typed edges — the recursive source supports the same depends_on rna-free-energy-turner-model / depends_on rna-minimum-free-energy-folding / relates_to test-unit-registry / relates_to rna-dot-bracket-notation edges already declared on the concept from rna-pkpredict-001-evidence. Contradictions: none — 2004/2007 papers, `Energy.lhs`, and Antczak 2018 mutually consistent. Follow-ups: the pknotsRG PARTIAL limitation the single-knot unit flagged is now covered by this unit for the recursive class (no longer an open wiki gap); remaining RNA surfaces (stem-loop enumeration, base-pair classification) still not individually ingested.
   graph: +1 node, +0 typed edges (new source node; concept-to-concept edges already declared from the sibling single-knot evidence)
- 2026-07-10 — ingest docs/Evidence/RNA-PSEUDOKNOT-001-Evidence.md (Pseudoknot Detection — identify crossing base pairs in a GIVEN structure; the detection/analysis facet of the RNA crossing-helix family, `RnaSecondaryStructure.DetectPseudoknots`). GENUINELY DISTINCT from the existing energy-driven predictor concept [[rna-pseudoknot-prediction]] (RNA-PKPREDICT/PKRECURSIVE): this unit takes a **base-pair set** (not a sequence) and runs a pure **O(n²) combinatorial scan** for crossing pairs — no folding, no energy model — whereas prediction folds a sequence via O(n⁴) Turner-NN energy DP. Detection is exactly the crossing primitive the predictor's validity invariant leans on (`DetectPseudoknots` finds ≥1 genuine crossing when a knot is returned) → NEW concept [[rna-pseudoknot-detection]] (per the ingest note that this could be a distinct facet — confirmed distinct). Core: two pairs (i,j),(k,l) written open<close **cross** iff `i<k<j<l` (Antczak 2018 verbatim `i<j<i'<j'`); two exhaustive negatives **nested** `i<k<l<j` + **disjoint** `j<k`. Each crossing pair-of-pairs = one pseudoknot (binary relation); invariants ≥2 pairs required / endpoints normalized min-max to open<close before the test / deterministic / every reported knot satisfies `i<k<j<l`. Oracles `([)]`=(0,2)+(1,3)→one pseudoknot, nested (0,5)+(1,4)→none, disjoint (0,2)+(3,5)→none. Sources: Antczak et al. 2018 *Bioinformatics* 34(8):1304 (rank 1 — crossing/conflict + pseudoknot **order** = min base-pair-set decompositions to nested + DBL notation order 0 `()`/1 `[]`/2 `{}`/3 `<>`/4–8 letters, H-type `([)]`) + Smit, Rother, Heringa & Knight 2008 *RNA* 14(3):410 (rank 1 — presence requires crossing pairs; pseudoknot-removal / order-assignment family) + biotite.structure.pseudoknots (rank 3 — nested order 0 / knotted order 1+) + Wikipedia Pseudoknot (rank 4, cites Rivas & Eddy 1999); all agree, no contradictions. One scope note: pseudoknot-**order** grouping (DBL layering) Not Implemented — reports the binary crossing relations, not higher-order layering (documented, not an invented parameter). Wrote concept [[rna-pseudoknot-detection]] + source page [[rna-pseudoknot-001-evidence]]; enriched sibling [[rna-pseudoknot-prediction]] (validity-invariant sentence now cross-links the detection primitive + RNA-PSEUDOKNOT source + source_commit→HEAD); updated hub [[algorithm-validation-evidence]] (frontmatter sources + source_commit→HEAD + body evidence-link list + concept enumeration — also backfilled the previously-missing [[rna-pseudoknot-prediction]] entry alongside the new detection concept) and wiki/index.md (Sources-section + Concepts-section entries). Graph: concept-to-concept edges on the new concept page — relates_to test-unit-registry + relates_to rna-pseudoknot-prediction (the crossing primitive the predictor's validity invariant leans on; shared Antczak 2018 crossing condition). Contradictions: none. Follow-ups: the RNA-family follow-ups list can drop "pseudoknot detection" (now ingested); remaining RNA surfaces (stem-loop enumeration, base-pair classification) still not individually ingested.
   graph: +2 nodes, +2 typed edges (relates_to test-unit-registry; relates_to rna-pseudoknot-prediction)
- 2026-07-10 — ingest docs/Evidence/RNA-STEMLOOP-001-Evidence.md (Stem-Loop / Hairpin Detection, area `RnaStructure`, `RnaSecondaryStructure.FindStemLoops`/`FindHairpins`/`FindPseudoknots`). GENUINELY DISTINCT — this is the **enumeration layer** of the RNA secondary-structure family: it SCANS a sequence for **every** stem-loop it can form (extending an antiparallel WC/wobble stem outward around each candidate loop, **O(n²·L)**), which none of the existing pages cover. It is NOT the hairpin *energy* calculator [[rna-hairpin-001-evidence]] (which explicitly disclaims being an enumerator), NOT the single-optimal-fold DP [[rna-minimum-free-energy-folding]] (O(n³)), and NOT the miRNA-specific precursor detector [[pre-mirna-hairpin-detection]] — prior ingests/lints repeatedly flagged "stem-loop enumeration" as a not-yet-ingested RNA surface, and the pre-miRNA page named a "general RNA secondary-structure (stem-loop) concept" as a future neighbour. → NEW concept [[rna-stem-loop-enumeration]]. Distinctive surfaces: (1) the general sequence-scanning enumerator with size params (minStem/minLoop/maxLoop/allowWobble); (2) **tetraloops** — GNRA (GAAA/GCAA/GGAA/GUAA) / UNCG (UACG/UCCG/UGCG/UUCG) / CUUG, UUCG the most stable (Antao 1991), UUCG+GNRA ~70% of 16S-rRNA tetraloops (Woese 1990), ~3.0 kcal/mol bonus via Turner 2004; (3) steric loop floor ≥3 nt / optimal 4-8 nt; (4) `FindPseudoknots` sharing the `i<k<j<l` crossing primitive with [[rna-pseudoknot-detection]]. Sources: Wikipedia (Stem-loop/Tetraloop/Pseudoknot) + Woese 1990 (PNAS 87:8467) + Heus & Pardi 1991 (Science 253:191) + Antao 1991 (NAR 19:5901) + Rivas & Eddy 1999 (JMB 285:2053) + NNDB Turner 2004. Oracles `GGGAAAACCC`→`(((....)))`, `GGGCGAAAGCCC`→GNRA tetraloop, `AAAA…`/`GC`/`GCAUC`→none, base pairs (0,6)+(3,9)→crossing. Invariants non-overlap / ≤1 pair per base / contiguous loop / antiparallel stem; defaults wobble-on, minLoop 3, minStem 3. Limitations: no PK prediction from sequence, simplified energy vs ViennaRNA, no internal loops/bulges (hairpin loops only), single structure (no suboptimals). Wrote concept [[rna-stem-loop-enumeration]] + source page [[rna-stemloop-001-evidence]]; enriched [[pre-mirna-hairpin-detection]] (future-neighbour note now points to the enumerator + RNA-STEMLOOP source + source_commit→HEAD) and [[rna-hairpin-001-evidence]] (its "not a stem-loop enumerator" disclaimer now links the actual enumerator); updated hub [[algorithm-validation-evidence]] (frontmatter sources + source_commit→HEAD + body evidence-link list) and wiki/index.md (Sources-section + Concepts-section entries). Graph: concept-to-concept edges on the new concept — relates_to test-unit-registry + depends_on rna-base-pairing (stem extension checks CanPair) + relates_to rna-dot-bracket-notation (emits dot-bracket) + relates_to rna-pseudoknot-detection (shared crossing primitive) + alternative_to rna-minimum-free-energy-folding (enumerate-all vs single-optimal). Contradictions: none — Wikipedia, Woese 1990, Antao 1991, Heus & Pardi 1991, Rivas & Eddy 1999, NNDB Turner 2004 mutually consistent. Follow-ups: RNA-family follow-up list can drop "stem-loop enumeration" (now ingested); the base-pair *classification* surface (`GetBasePairType`) is already covered by [[rna-base-pairing]]; remaining not-individually-ingested RNA surface — none prominent left from the earlier gap list.
   graph: +2 nodes, +5 typed edges (relates_to test-unit-registry; depends_on rna-base-pairing; relates_to rna-dot-bracket-notation; relates_to rna-pseudoknot-detection; alternative_to rna-minimum-free-energy-folding)
- 2026-07-10 — ingest docs/Evidence/RNA-STRUCT-001-Evidence.md (Secondary Structure Prediction, area `RnaStructure`, `RnaSecondaryStructure.Predict`/`PredictWithConstraints`/`ToDotBracket`/`FromDotBracket`). GENUINELY DISTINCT — the **top-level structure-prediction umbrella** whose headline algorithm is a **Nussinov & Jacobson 1980** O(n³)/O(n²) base-pair-**maximizing** DP (weighted pair scores WC −2.0 / wobble −1.0, "relative stability, not physical energy units"), plus **constraint folding** (`PredictWithConstraints`, Mathews 2004 forced pairs) and **dot-bracket round-trip conversion** (`ToDotBracket`/`FromDotBracket`) → NEW concept [[rna-secondary-structure-prediction]]. CONTRADICTION RECONCILED: prior RNA pages (MFE concept [[rna-minimum-free-energy-folding]] + source [[rna-mfe-001-evidence]] + index + hub) claimed the generic id *RNA-STRUCT-001* was an **alias** of the physical-MFE folder RNA-MFE-001 ("the two names denote this one MFE-folding unit"); now that this artifact is ingested that is **superseded** — RNA-STRUCT-001 is a **distinct sibling** test unit (Nussinov base-pair-max + constraints + notation I/O) that merely **shares** the same `RnaSecondaryStructure` Zuker MFE machinery (its deviation D5 added `CalculateMfeStructure`/`PredictStructureMfe`, the traceback partner of RNA-MFE-001's `PredictStructure`). Fixed the alias claim in all four places + added RNA-STRUCT-001 to the MFE concept's frontmatter sources + source_commit→HEAD. Sources: Nussinov & Jacobson 1980 (*PNAS* 77(11):6309) + Zuker & Stiegler 1981 (*NAR* 9(1):133) + MIT 6.047 Lec 08 (Washietl 2012, explicit F/C/M/M¹ recurrences) + Turner 2004/NNDB + Mathews 2004 (*PNAS* 101(19):7287) + Wikipedia (Nussinov/secondary-structure). Oracles simple hairpin `GGGGAAAACCCC`→`((((....))))`, GNRA `GCGCGAAACGCGC`→GA first-mismatch −0.9, tRNA-like 72-nt cloverleaf, poly-A→no pairs/MFE 0. Invariants dot-bracket balance / MFE sign ≤ 0 / WC stacking negative / loop init positive / non-overlap / ≤1 pair per base; empty/null→MFE 0, min hairpin 3-bp stem + 3-nt loop, case-insensitive; pseudoknots **detected not predicted** (`i<k<j<l`, shared with [[rna-pseudoknot-detection]]). Deviations D1 bulge degeneracy + D2 dangling ends + D5 Zuker traceback RESOLVED; D3 int21 (2,304-entry) / D4 int22 (36,864-entry) internal-loop lookup tables BLOCKED (too large for inline static data). Wrote concept [[rna-secondary-structure-prediction]] + source page [[rna-struct-001-evidence]]; updated hub [[algorithm-validation-evidence]] (frontmatter sources + source_commit→HEAD + inline evidence-link enumeration + concept enumeration + reconciled the MFE-entry alias phrase) and wiki/index.md (Sources + Concepts entries + reconciled the MFE alias phrase); reconciled [[rna-minimum-free-energy-folding]] + [[rna-mfe-001-evidence]]. Graph: concept-to-concept edges on the new concept — relates_to test-unit-registry + **alternative_to rna-minimum-free-energy-folding** (base-pair-max vs physical energy-min, same problem) + depends_on rna-base-pairing (Nussinov maximizes WC/wobble pairs) + relates_to rna-dot-bracket-notation (ToDotBracket/FromDotBracket I/O) + relates_to rna-free-energy-turner-model (stem-loop energy model) + relates_to rna-pseudoknot-detection (shared crossing test). Contradictions: the one reconciled above (RNA-STRUCT-001 ≠ RNA-MFE-001); no source-vs-source contradictions (Nussinov 1980, Zuker 1981, MIT 6.047, Turner 2004/NNDB, Mathews 2004 mutually consistent). Follow-ups: none new — the RNA secondary-structure family's headline surfaces (base-pairing, dot-bracket, Turner energy, MFE folding, partition function, stem-loop enumeration, pseudoknot detection/prediction, and now the top-level Nussinov prediction umbrella) are all individually ingested.
   graph: +2 nodes, +6 typed edges (relates_to test-unit-registry; alternative_to rna-minimum-free-energy-folding; depends_on rna-base-pairing; relates_to rna-dot-bracket-notation; relates_to rna-free-energy-turner-model; relates_to rna-pseudoknot-detection)
- 2026-07-10 — ingest docs/Evidence/SEQ-ATSKEW-001-Evidence.md (AT skew — `(A−T)/(A+T)`, the A/T strand-asymmetry sibling of GC skew; a nucleotide-composition/skew statistic). CONTEXT check: searched wiki/concepts for gc-skew / composition / nucleotide-composition / replication-origin — NO existing skew or composition-skew concept (only passing GC-skew mentions in [[centromere-analysis]] flagging it as a future concept, and the dinucleotide CpG O/E in [[cpg-island-detection]]); the backlog lists `at-skew`/`gc-skew` slugs but neither page existed. AT skew is genuinely unrepresented → created ONE reusable FAMILY concept [[nucleotide-composition-skew]] covering BOTH AT skew and its GC-skew sibling (rather than a page per member), so a future GC-skew ingest enriches it. Formula fully sourced: Lobry 1996 *Mol Biol Evol* 13(5):660 (PMID 8676740, primary — founding intra-strand base-asymmetry observation) + Charneski et al. 2011 *PLoS Genet* 7(9):e1002283 (verbatim `(A−T)/(A+T)`; Firmicute AT skew from SELECTION not mutation) + Wikipedia "GC skew" (both formulas + range −1…+1, AT skew −1⇔A=0 / +1⇔T=0) + Biopython `Bio.SeqUtils.GC_skew` (symbol conventions: case-insensitive counting, zero-denominator ⇒ 0.0, ambiguous/non-canonical bases ignored). Hand-derived oracles `AAAA→1.0`, `TTTT→−1.0`, `ATAT→0.0`, `AAAT→0.5`, `ATTT→−0.5`, `GGCC→0.0` (no A/T), `AAATGGGCCC→0.5` (G/C ignored), `aaat→0.5` (case-insensitive). One documented ASSUMPTION: the lowercase + non-ACGT handling for the AT-skew member is inferred by analogy from the shipped `GC_skew` (Biopython ships no AT-skew line) — the formula itself is fully sourced, only the symbol-handling convention is by analogy, and it matches the repository (`ToUpperInvariant`, counts only A/T). Wrote concept [[nucleotide-composition-skew]] + source page [[seq-atskew-001-evidence]] (cross-linked to hub + registry + cpg-island/centromere cousins); updated hub [[algorithm-validation-evidence]] (frontmatter sources +SEQ-ATSKEW-001) and wiki/index.md (Sources + Concepts entries); enriched [[cpg-island-detection]] (CpG O/E now names its single-base skew cousin) and [[centromere-analysis]] (GC-skew mention now links the concept). Graph: no typed concept-to-concept edges — the source supports the AT/GC sibling relation WITHIN the one concept, but no cross-concept edge is explicitly source-backed (CpG/centromere links are wiki-navigational, not asserted by the AT-skew source). Contradictions: none — Lobry, Charneski, Wikipedia, and the Biopython convention agree on formula and range. Follow-ups: the sibling GC-skew unit (`docs/algorithms/Sequence_Composition/GC_Skew.md`, backlog slug `gc-skew`) is not yet ingested — when it is, it enriches [[nucleotide-composition-skew]] rather than creating a new page; the cumulative-skew replication-origin locator is a chromosome-scale application not yet a separate unit.
   graph: +2 nodes, +0 typed edges (new concept + source nodes; no source-backed concept-to-concept edges)
- 2026-07-10 — ingest docs/Evidence/SEQ-CODON-FREQ-001-Evidence.md (Codon Frequencies — `SequenceStatistics.CalculateCodonFrequencies(dnaSequence, readingFrame=0)`, Analysis assembly). Surveyed the existing codon-usage family (concepts codon-usage-comparison / relative-synonymous-codon-usage / codon-adaptation-index / effective-number-of-codons / rare-codon-analysis / codon-optimization; sources codon-usage-001 / codon-stats-001 / annot-codonusage-001 / codon-rscu/cai/rare/opt). This unit is the **normalized, frame-aware frequency** view — a genuinely distinct METHOD but the SAME family, so per the economy directive NO new concept was created; enriched the closest concept [[codon-usage-comparison]] instead. Distinct from the raw-count sibling `CodonOptimizer.CalculateCodonUsage` (CODON-USAGE-001) on four points: returns count/total fractions (`IReadOnlyDictionary<string,double>`) not `int` counts; adds a **reading-frame offset** (0/1/2) so the same sequence yields a different multiset (the distinctive new semantic); **excludes non-ACGT/ambiguous triplets** from count and total (Kazusa CUTG "ambiguous codons excluded"); keeps **DNA-native** keys (no T→U rewrite). Shared family behaviour: case-insensitive, incomplete trailing 1–2-nt codon dropped, frequencies sum to 1.0 (INV-02). Sources (rank in parens): Kazusa CUTG README (5, canonical count/total per-thousand convention + ambiguous exclusion) + EMBOSS `cusp` (3, its **Fraction** column is the per-AA RSCU-style metric ≠ this frequency; **Frequency**=count/1000; verbatim 386-codon sample cross-checks 22/386=56.995‰, 23/386=59.585‰) + Wikipedia codon-usage-bias (4) + Nakamura, Gojobori, Ikemura 2000 *NAR* 28(1):292 (1, the paper behind Kazusa CUTG). Oracles: `ATGATGAAA` f0→ATG 2/3,AAA 1/3 · f1→TGA 1.0; `ATGNNNAAA` f0→ATG 1/2,AAA 1/2 (NNN excluded); `ATGAA`→ATG 1.0 (trailing AA dropped); `atgaaa`→ATG 1/2,AAA 1/2. Single ASSUMPTION: empty table when total=0 (Kazusa leaves total=0 undefined; empty is the only count/total-consistent value, matches the guard). Wrote source page [[seq-codon-freq-001-evidence]]; enriched [[codon-usage-comparison]] (new "normalized, frame-aware sibling" paragraph + frontmatter sources +SEQ-CODON-FREQ-001 + source_commit→HEAD) and cross-linked the raw-count sibling source [[codon-usage-001-evidence]] (Related-units line); updated hub [[algorithm-validation-evidence]] (frontmatter sources +SEQ-CODON-FREQ-001 + source_commit→HEAD + body evidence-link list) and wiki/index.md (Sources-section entry). Graph: no typed concept-to-concept edges — the source supports the raw-count↔frequency sibling relation, but `relates_to` does not accept a source-typed subject and no NEW concept was created, so the relation is captured in prose + `mentions` edges only. Contradictions: none — Kazusa, EMBOSS cusp, Wikipedia, Nakamura 2000 agree on the count/total (per-thousand) convention + ambiguous exclusion; the one nuance recorded is cusp's Fraction≠this-frequency. Follow-ups: none — the codon-usage family's counting/frequency/RSCU/CAI/ENC/rare/optimization/stats surfaces are all individually ingested.
   graph: +1 node, +0 typed edges (new source node; no source-backed concept-to-concept edge — no new concept, source subjects invalid for relates_to)
- 2026-07-10 — ingest docs/Evidence/SEQ-COMPLEX-COMPRESS-001-Evidence.md (Lempel–Ziv (LZ76) compression-based sequence complexity `c(S)` = number of distinct phrases in a left-to-right exhaustive-history parse). CONTEXT check: searched wiki/concepts for complexity/entropy/low-complexity/dust/linguistic — the only complexity concepts were the **protein** Shannon-entropy SEG detector [[protein-low-complexity-seg]], the fixed-`k` Shannon k-entropy [[k-mer-statistics]], and the explicit-repeat [[repetitive-element-detection]]; NO nucleotide/general sequence-complexity concept and NO compression-based measure existed. LZ76 is a **genuinely distinct scalar measure** (adaptive variable-length phrase count over the whole sequence, sensitive to ordered pattern buildup — not composition entropy, not a repeat detector) → created ONE dedicated concept [[sequence-complexity-compression-lempel-ziv]], cross-linked to those three siblings as the compression member of the complexity/entropy family. Definition/parse fully sourced: Lempel & Ziv 1976 *IEEE Trans. Inf. Theory* 22(1):75–81 (primary, paywalled) + Wikipedia "Lempel–Ziv complexity" (verbatim definition + delimiter parsing rule + O(n) pointer-scan pseudocode) + Naereen/Lempel-Ziv_Complexity (set-of-seen-substrings reference impl + 4 binary doctests) + entropy/AntroPy `lziv_complexity` (normalization `LZ_n = c/(n/log_b n)`, b = distinct symbols; cites Zhang et al. 2009 *J. Math. Chem.* 46(4):1203). Oracles: raw 8/7/9/10 (Naereen doctests, component lists given); normalized 2.0 for `1001111011000010`; homopolymer `"0"×16` raw 5 (`0/00/000/0000/00000`, general `⌊(√(8n+1)−1)/2⌋`) / normalized **1.25** (b<2 clamps log base to 2 and returns the *normalized* value, NOT the raw count — an earlier raw reading was corrected to 1.25 on 2026-06-16); `ACGT`→4, `AAAA`→2, empty→0; monotone in repetitiveness. Two ASSUMPTIONs, both source-flagged: Naereen trailing-component convention (Wikipedia pseudocode adds 1, Naereen set does not — ±1 on last component), and the b<2 log-base clamp. Wrote source page [[seq-complex-compress-001-evidence]] + concept [[sequence-complexity-compression-lempel-ziv]]; updated hub [[algorithm-validation-evidence]] (frontmatter sources +SEQ-COMPLEX-COMPRESS-001 + source_commit→HEAD + body evidence-link list + concept-links list) and wiki/index.md (Sources + Concepts entries); added reciprocal cross-links from [[k-mer-statistics]] and [[protein-low-complexity-seg]] (navigation-only, frontmatter untouched — those pages are not derived from this source). Graph: one typed edge on the new concept, `relates_to concept:test-unit-registry` (the standard per-unit registration edge); no source-backed concept-to-concept edge (the source defines the measure, it does not assert a typed relation to SEG/k-mer/repeats — those are wiki-navigational siblings). Contradictions: none — Lempel-Ziv 1976, Wikipedia, Naereen, and entropy/antropy agree on the LZ76 "distinct substrings" measure and the `n/log_b n` normalization; the only internal item was the b<2 clamp correction, resolved by reading the antropy source. Follow-ups: the other complexity/entropy-family members named in the bio-qc framing (Shannon character-entropy, DUST, linguistic complexity) are not yet ingested — when they are, each is a distinct measure that can sit alongside this concept (a future broad "sequence complexity" hub could group them if the family grows).
   graph: +2 nodes, +1 typed edge (relates_to test-unit-registry)
- 2026-07-10 — ingest docs/Evidence/SEQ-COMPLEX-DUST-001-Evidence.md (DUST triplet-frequency low-complexity score — `score = ∑_t c(c−1)/2 / (L−2)` over nucleotide triplet counts; the DNA-sequence low-complexity masker, sibling of the just-ingested Lempel–Ziv compression complexity). CONTEXT check: read the complexity-family concepts ([[sequence-complexity-compression-lempel-ziv]], [[k-mer-statistics]], [[protein-low-complexity-seg]]) + searched wiki for dust/low-complexity/entropy — the protein-side low-complexity masker (SEG, Shannon entropy) and the DNA repeats anchor ([[repetitive-element-detection]]) exist, but NO DNA triplet-frequency DUST score and NO general nucleotide low-complexity masking score. DUST is a **distinct, well-known algorithm** (Morgulis et al. 2006, triplet-based, the canonical low-complexity DNA masker) — like LZ76 it got its own concept → created dedicated concept [[dust-low-complexity-score]], cross-linked as the triplet-frequency masking member of the complexity/entropy family. Formula fully sourced: Morgulis, Gertz, Schäffer & Agarwala 2006 *J Comput Biol* 13(5):1028 (PMID 16796549, primary — "triplet frequencies in 64-base windows"; the SDUST rewrite kept the scoring function and only made masking symmetric/context-insensitive) + Li 2025 longdust arXiv:2509.07357 (verbatim `∑_t c_x(t)(c_x(t)−1)/2 / (L−2)`, k=3 hardcoded, L−2 = number of triplets, default window 64 / threshold 2.0 = level 20, **HIGH score ⇒ LOW complexity**) + lh3/sdust `sdust.c` (incremental `rw += cw[t]++` proven = closed-form `c(c−1)/2` sum, `if (rw*10 > L*T)` with W=64/T=20 ⇒ score>2.0). Hand-derived oracles `ATGC`→0.0 (all-distinct triplets ⇒ max complexity) / `ACGTACGT`→0.333… / `AAAAAA`→1.5 / `ACACACAC`→1.0 / `AAAAAAAAAA`→3.5; homopolymer max `(L−3)/2`; `L<3` undefined. Two ASSUMPTIONs, both source-flagged: general `wordSize` normalization `L−w+1` (only k=3 source-backed, exact oracles asserted for k=3 only), and `L<wordSize`→0 defined-output convention. Wrote source page [[seq-complex-dust-001-evidence]] + concept [[dust-low-complexity-score]]; updated hub [[algorithm-validation-evidence]] (frontmatter sources +SEQ-COMPLEX-DUST-001 + source_commit→HEAD + body evidence-link list + concept-enumeration entry) and wiki/index.md (Sources + Concepts entries); added a reciprocal DUST sibling bullet to [[sequence-complexity-compression-lempel-ziv]] (navigation-only, frontmatter untouched — not derived from this source). Graph: one typed edge on the new concept, `relates_to concept:test-unit-registry` (the standard per-unit registration edge); no source-backed concept-to-concept edge (the source defines the DUST score, it does not assert a typed relation to LZ/SEG/k-mer/repeats — those are wiki-navigational siblings). Contradictions: none — Morgulis 2006, Li 2025, and lh3/sdust agree on the score, the w=64/T=20 (score 2.0) defaults, and the HIGH-score⇒LOW-complexity direction; the incremental accumulation is proven algebraically equal to the closed form. Follow-ups: the remaining named complexity/entropy members (Shannon character-entropy, linguistic complexity) are still not ingested — when they are, each is a distinct measure alongside this concept (a future broad "sequence complexity" hub could group DUST + LZ + Shannon + linguistic + k-entropy if the family keeps growing).
   graph: +2 nodes, +1 typed edge (relates_to test-unit-registry)
- 2026-07-10 — ingest docs/Evidence/SEQ-COMPLEX-KMER-001-Evidence.md (k-mer entropy — Shannon entropy `H = −Σ pᵢ log₂ pᵢ` of the overlapping k-mer frequency distribution, `pᵢ = nᵢ/(L−k+1)`, bits; the entropy member of the `SEQ-COMPLEX-*` sequence complexity family, sibling of the just-ingested LZ76 and DUST). CONTEXT check: read the four complexity-family concepts ([[k-mer-statistics]], [[sequence-complexity-compression-lempel-ziv]], [[dust-low-complexity-score]], [[protein-low-complexity-seg]]) and confirmed in source (`SequenceComplexity.CalculateKmerEntropy`, `SequenceComplexity.cs:136-185`) that this unit computes the **identical formula** to the Shannon **k-entropy** already written up as the `Entropy` field of [[k-mer-statistics]] (KMER-STATS-001, `KmerAnalyzer.AnalyzeKmers`) — same `H=−Σ (nᵢ/N)log₂(nᵢ/N)`, N=L−k+1 over the same overlapping multiset, just a different class (the `SEQ-COMPLEX-*` complexity-family entry point, sibling of `CalculateLinguisticComplexity`/`EstimateCompressionRatio`/`CalculateDustScore`). NOT genuinely distinct → per the ingest guidance, **enriched [[k-mer-statistics]]** with a "Second entry point" section rather than creating a redundant concept. Formula fully sourced: Li 2025 longdust arXiv:2509.07357 (overlapping `L−k+1` k-mers; `H=−Σpᵢlog₂pᵢ`; skewed⇒low / uniform⇒high) + Çakır 2025 Entropy–Rank Ratio arXiv:2511.05300 (same formula, base-2⇒bits, single-nt max log₂4=2, saturation to log λ) + Shannon 1948 via expositions (bounds `0≤H≤log(k)`, deterministic⇒0 / uniform⇒log_b n). Hand-derived oracles (bits): `ACGT` k=1→2.0 (uniform) / `ACGT` k=2→log₂3≈1.5849625 (all-distinct) / `ATATAT` k=2→0.9709505945 (binary entropy of p=0.6) / `AAAA` k=2→0.0 (homopolymer) / `AAACGT` k=2→1.9219280949 (=log₂5−0.4) / `AC` k=5→0.0 (L<k). Two ASSUMPTIONs, both API-shape only (no entropy value affected): `L<k`→0 (empty-multiset convention, matches SequenceComplexity siblings), and invalid `k<1`→ArgumentOutOfRangeException / null DnaSequence→ArgumentNullException / null-empty string→0 (sibling method guards). Wrote source page [[seq-complex-kmer-001-evidence]]; enriched [[k-mer-statistics]] (frontmatter +SEQ-COMPLEX-KMER-001 source + source_commit→HEAD + updated→2026-07-10 + new graph edge `relates_to test-unit-registry` from the new evidence slug; body "Second entry point" section + expanded the complexity-family bullet to name LZ76 + DUST siblings). Reciprocal navigation cross-links added to [[sequence-complexity-compression-lempel-ziv]] and [[dust-low-complexity-score]] (their "vs k-mer k-entropy" bullets now name the validated SEQ-COMPLEX-KMER-001 = `CalculateKmerEntropy` unit; frontmatter untouched — not derived from this source). Updated hub [[algorithm-validation-evidence]] (frontmatter sources +SEQ-COMPLEX-KMER-001 + source_commit→HEAD + body evidence-link list) and wiki/index.md (Sources entry + a note appended to the [[k-mer-statistics]] concept line). Graph: one typed edge, `relates_to concept:test-unit-registry` (source=seq-complex-kmer-001-evidence, the standard per-unit registration edge), placed on the [[k-mer-statistics]] concept page (no new concept page created); no source-backed concept-to-concept edge. Contradictions: none — Li 2025, Entropy–Rank Ratio, and the Shannon expositions agree on formula/base-2/bounds; this unit and KMER-STATS-001 compute the identical measure from two classes (a duplicate entry point, explicitly noted, not a contradiction). Follow-ups: the remaining named complexity/entropy members (Shannon character-entropy, linguistic complexity) are still un-ingested — when they arrive, a future broad "sequence complexity" hub could group DUST + LZ + k-entropy + Shannon + linguistic if the family keeps growing.
   graph: +1 node, +1 typed edge (relates_to test-unit-registry; edge on k-mer-statistics concept, no new concept)
- 2026-07-10 — ingest docs/Evidence/SEQ-COMPLEX-WINDOW-001-Evidence.md (Windowed Sequence Complexity — a sliding-window scan emitting a per-position `ComplexityPoint` *profile*, each window carrying a Shannon entropy of base composition AND a linguistic complexity value; `SequenceComplexity` area). CONTEXT check: read the four complexity-family concepts ([[sequence-complexity-compression-lempel-ziv]] LZ76, [[dust-low-complexity-score]] DUST, [[k-mer-statistics]] k-entropy, [[protein-low-complexity-seg]] SEG). GENUINELY DISTINCT → NEW concept [[windowed-sequence-complexity-profile]]: every existing SEQ-COMPLEX-* sibling reduces a whole sequence to a single **scalar**, whereas this unit is the family's **profiling / scanning driver** (per-position profile) and its **first home of linguistic complexity** (`LC = Σ Vᵢ / Σ Vmax,i`, `Vmax,i = min(4^i, N−i+1)`, `maxWordLength = min(6, w)`; the standalone SEQ-COMPLEX-001 LC unit is not yet separately ingested). Per-window Shannon is over single-base composition (n=4, uniform ⇒ log₂4 = 2.0, homopolymer ⇒ 0). Fully sourced: Troyanskaya et al. 2002 *Bioinformatics* 18(5):679 (PMID 12050064, primary — complexity **profiles** as a function of position; linguistic complexity = distinct-subword counts, suffix-tree linear-time) + Wikipedia "Linguistic sequence complexity" (Gabrielian & Bolshoy 1999 / Trifonov 1990 — `Uᵢ = Vᵢ/min(4^i,N−i+1)`, range 0<C<1) + Wikipedia "Shannon entropy" (Shannon 1948 — `H = −Σ p log₂ p` bits, uniform max log₂4=2, deterministic⇒0, `0·log0=0`). Hand-derived oracles `ACGTACGT` window → Shannon 2.0 / LC 23/29 = 0.7931034482758621, `AAAAAAAA` poly-A window → Shannon 0.0 / LC 6/29 = 0.20689655172413793, geometry L=24/w=8/s=8 → 3 windows ([0,7]c4 / [8,15]c12 / [16,23]c20). Contract: window count `floor((L−w)/s)+1` for L≥w else 0, only fully-contained windows emitted (`i+w ≤ L`, NO trailing partial ⇒ L<w gives empty profile), Position = `WindowStart + w/2` (int div, 0-based end-inclusive), bounds 0≤H≤2 / 0≤LC≤1; null DnaSequence → ArgumentNullException, windowSize<1 & stepSize<1 → ArgumentOutOfRangeException. Three ASSUMPTIONs, all non-value-affecting: center-position label, default windowSize=64/stepSize=10, maxWordLength cap 6. Wrote source page [[seq-complex-window-001-evidence]] + concept [[windowed-sequence-complexity-profile]]; updated hub [[algorithm-validation-evidence]] (frontmatter sources +SEQ-COMPLEX-WINDOW-001 + source_commit→HEAD + body evidence-link list + concept enumeration) and wiki/index.md (Sources + Concepts entries); added reciprocal navigation cross-links from [[sequence-complexity-compression-lempel-ziv]] and [[dust-low-complexity-score]] (their "Where it sits" family sections now name the windowed profiling member; frontmatter untouched — not derived from this source). Graph: one typed edge on the new concept, `relates_to concept:test-unit-registry` (the standard per-unit registration edge); no source-backed concept-to-concept edge (the source defines the profiling method + linguistic complexity, it does not assert a typed relation to the scalar siblings — those are wiki-navigational). Contradictions: none — Troyanskaya 2002, the Wikipedia linguistic-complexity `Uᵢ` definition, and Shannon 1948 are mutually consistent; the `ACGTACGT` / poly-A oracles follow directly. Follow-ups: the standalone linguistic-complexity unit SEQ-COMPLEX-001 (`docs/algorithms/Sequence_Composition/Linguistic_Complexity.md`) is still not separately ingested — when it is, it enriches [[windowed-sequence-complexity-profile]] (the per-sequence LC scalar this windows) rather than creating a new page; Shannon character-entropy remains the last un-ingested named complexity member.
   graph: +2 nodes, +1 typed edge (relates_to test-unit-registry)
- 2026-07-10 — ingest docs/Evidence/SEQ-COMPOSITION-001-Evidence.md (Sequence Composition — exact nucleotide counts/fractions A/T/G/C/U/N/Other partitioning Length, GC content `(G+C)/(A+T+G+C+U)` ∈ [0,1], and amino-acid residue counts; `Composition` area). CONTEXT check: searched wiki/concepts for composition/gc-content/base-count/nucleotide/sequence-statistics and read the sibling [[nucleotide-composition-skew]] (SEQ-ATSKEW-001). The **skew** family (GC/AT skew asymmetry) is already a concept, but NO concept covered the foundational **base/residue composition** layer — the per-symbol counts, the A/T/G/C/U/N/Other partition, GC content, or amino-acid composition. GENUINELY DISTINCT (magnitude/fraction view vs the skew asymmetry view; different outputs) → NEW concept [[base-composition]], cross-linked as the sibling of [[nucleotide-composition-skew]] (which derives from the same counts) and the base layer under [[cpg-island-detection]] (dinucleotide O/E), [[windowed-sequence-complexity-profile]] (window Shannon of base composition), and the GC constraints in [[taqman-probe-design-rules]]/[[codon-optimization]]. Fully sourced: Biopython `Bio.SeqUtils.gc_fraction` (`gc=sum(count(x) for x in "CGScgs")`, `gc/length` ∈ [0,1], **empty ⇒ 0**, case-insensitive) + Wikipedia "GC skew" (GC/AT skew formulas, cumulative-skew origin/terminus — skew detail delegated to the sibling concept) + IUPAC codes (canonical A/C/G/T/U + degenerate/N; the 20 standard amino-acid letters). Hand-derived oracles `ATGC`→GC 0.5/skew 0, `GGGC`→GC 1.0/GC-skew 0.5/AT-skew 0 (a+t=0), `AAUUGGCC`→A/T/G/C/U 2/0/2/2/2, GC 0.5/AT-skew 1.0, `MKVLWA`→6 residues each 1. Corner cases: empty/null ⇒ all-zero composition; no G or C ⇒ zero-denominator skew 0.0; mixed case ⇒ case-insensitive; non-canonical letters (N, degenerate, X) tracked separately via `CountN`/`CountOther`, never folded into A/T/G/C/U. One ASSUMPTION (source-flagged): degenerate IUPAC codes excluded from GC/AT totals — Biopython's `gc_fraction` counts S toward GC and W toward the denominator, whereas the repo counts only A/T/G/C/U; the two agree **exactly** over the {A,T,G,C,U} alphabet (this unit's scope), differing only on degenerate symbols. NOTE: the Change History records SEQ-COMPOSITION-001 is a **duplicate/consolidated Registry entry** for the two composition methods already delivered under SEQ-STATS-001 (TestSpec §7) — a documentation/registry artifact over an existing implementation, not a new algorithm. Wrote source page [[seq-composition-001-evidence]] + concept [[base-composition]]; updated hub [[algorithm-validation-evidence]] (frontmatter sources +SEQ-COMPOSITION-001 + body evidence-link list; source_commit left at its prior value) and wiki/index.md (Sources + Concepts entries); added a reciprocal sibling bullet to [[nucleotide-composition-skew]] (navigation-only, frontmatter untouched — not derived from this source). Graph: one typed edge on the new concept, `relates_to concept:test-unit-registry` (the standard per-unit registration edge); no source-backed concept-to-concept edge (the source defines composition; it does not assert a typed relation to skew/CpG/complexity — those are wiki-navigational siblings). Contradictions: none — Biopython, Wikipedia, and IUPAC agree on formulas and alphabet; the degenerate-code handling is a documented scope simplification, not a conflict. Follow-ups: **SEQ-STATS-001 is not yet ingested** — when it is, it is the original home of these two composition methods and should enrich [[base-composition]] rather than create a new page; the GC-skew member (`gc-skew`) flagged on [[nucleotide-composition-skew]] is also still un-ingested.
   graph: +2 nodes, +1 typed edge (relates_to test-unit-registry)
- 2026-07-10 — ingest docs/Evidence/SEQ-DINUC-001-Evidence.md (Dinucleotide Analysis — dinucleotide frequency `f_XY = count/(N−1)`, the Karlin genomic-signature odds ratio `ρ_XY = f_XY/(f_X·f_Y)`, and codon frequencies; method `CalculateDinucleotideRatios`, Analysis area). CONTEXT check: read [[base-composition]] (SEQ-COMPOSITION-001, single-base layer), [[cpg-island-detection]] (EPIGEN-CPG-001, CpG O/E), [[k-mer-statistics]] (2-mer diversity), [[codon-usage-comparison]] + [[seq-codon-freq-001-evidence]] (codon frequency), and searched wiki for dinucleotide/k-mer/Karlin/genomic-signature. GENUINELY DISTINCT → NEW concept [[dinucleotide-relative-abundance]]: the Karlin & Burge (1995) dinucleotide relative-abundance / genomic signature (odds ratio over all 16 base-steps, `ρ=1` no-bias independence baseline, over/under thresholds `ρ≥1.23`/`ρ≤0.78`) is not covered anywhere — CpG O/E in [[cpg-island-detection]] is literally its `CG`-specialized case (same odds-ratio shape, differing only by the `N` vs `N−1` dinucleotide-frequency normalization, `N/(N−1)→1` for long sequences, a documented modeling choice), and [[base-composition]] is the single-base layer it builds `f_X`/`f_Y` on. Fully sourced: Karlin PMC126251 (`ρ_XY = f_XY/(f_X·f_Y)`, normalized frequencies, `r=1.0` no bias) + Karlin & Burge 1995 via MBE 19(6):964 (thresholds `ρ≤0.78`/`ρ≥1.23`, documentation-only — method returns raw ratio) + Gardiner-Garden & Frommer 1987 (CpG O/E = same shape normalized by `N`) + Kazusa CUTG (codon frequency = count/total non-overlapping in-frame triplets, non-ACGT excluded). Exact-rational oracles `ATGCGCGT`: GC/CG ρ=64/21≈3.0476, AT ρ=32/7≈4.5714, TG/GT ρ=32/21≈1.5238; f_XY GC=CG=2/7, AT=TG=GT=1/7; codon `ATGATGAAA` f0→ATG 2/3,AAA 1/3, f1→TGA 1.0. Corner cases: div-by-zero guard (absent constituent base ⇒ expected 0 ⇒ ratio 0), case-insensitive, null/empty/len<2→empty ratios+freqs, len<3→empty codons. Two documented modeling choices (Karlin `(N−1)` vs Gardiner-Garden `(N)` normalization; `U` counted as a fifth base for RNA), no source contradictions. The source's codon-frequency output is the SAME metric as SEQ-CODON-FREQ-001 → cross-linked to [[codon-usage-comparison]]/[[seq-codon-freq-001-evidence]], NOT re-derived. Wrote source page [[seq-dinuc-001-evidence]] + concept [[dinucleotide-relative-abundance]]; enriched [[base-composition]] (dinucleotide-composition bullet now names the Karlin odds ratio + new page; frontmatter +SEQ-DINUC-001 source + source_commit→HEAD) and [[cpg-island-detection]] (O/E paragraph now frames CpG O/E as the `CG`-specialized case of the general Karlin ρ + the N vs N−1 normalization note; frontmatter +SEQ-DINUC-001 source + source_commit→HEAD + updated→2026-07-10); updated wiki/index.md (Sources + Concepts entries). Graph: two typed edges on the new concept — `relates_to concept:test-unit-registry` (standard per-unit registration) and `relates_to concept:cpg-island-detection` (source explicitly states CpG O/E is the same odds-ratio shape specialized to CG, differing by the N/(N−1) factor). Contradictions: none — Karlin, Karlin & Burge, Gardiner-Garden & Frommer, and Kazusa CUTG agree on the odds-ratio definition, the independence baseline, and the codon count/total convention. Follow-ups: none blocking — a future `dinucleotide-frequency`-only or `CalculateDinucleotideFrequencies` unit would enrich this same page.
   graph: +2 nodes, +2 typed edges (relates_to test-unit-registry; relates_to cpg-island-detection)
- 2026-07-10 — ingest docs/Evidence/SEQ-ENTROPY-PROFILE-001-Evidence.md (Shannon Entropy Profile — a sliding-window per-symbol Shannon entropy profile: per fixed-width window `H = −Σ pᵢ log₂ pᵢ` bits over the window's **single-base (k=1) mono-nucleotide** composition, emitted as a standalone entropy profile; `SequenceComplexity` area). CONTEXT check: read the complexity/entropy-family concepts ([[windowed-sequence-complexity-profile]] SEQ-COMPLEX-WINDOW-001, [[k-mer-statistics]] k-entropy, [[sequence-complexity-compression-lempel-ziv]] LZ76, [[dust-low-complexity-score]] DUST, [[base-composition]]). NOT genuinely distinct → this is the **SAME per-window Shannon measure and window geometry** already documented as the Shannon component of the `ComplexityPoint` profile on [[windowed-sequence-complexity-profile]], just emitted **alone** (entropy channel only, a `List<double>` rather than a list of `ComplexityPoint`s, no paired linguistic complexity). Following the established second-entry-point precedent (k-entropy hosted on [[k-mer-statistics]] as both `AnalyzeKmers.Entropy` and `CalculateKmerEntropy`), ENRICHED [[windowed-sequence-complexity-profile]] with a "Second entry point" section rather than creating a new concept. It is the **character-level (k=1) counterpart** of the k-mer/block Shannon k-entropy (SEQ-COMPLEX-KMER-001, per-window mono-nucleotide n=4 vs whole-sequence L−k+1 k-mer distribution). Fully sourced: Shannon 1948 *Bell Syst Tech J* 27(3):379 (primary — `H=−Σp log p`, `0·log0=0`) + Wikipedia "Entropy (information theory)" (base-2 ⇒ bits/shannons, uniform max `log_b n`) + IntechOpen ch. 75997 (DNA 4-letter alphabet, max 2 bits = log₂4, sliding counter of width W). Hand-derived per-window oracles `AAAA`→0.0 / `AATT`→1.0 / `ATGC`→2.0 (=log₂4 uniform) / `AAAT`→0.8112781244591328 (3:1 skew) / `AATG`·`GCAA`→1.5 / `AAATTC`→1.4591479170272448; sliding profiles `AAATGC` w4 s1→[0.8112781244591328,1.5,2.0] / `AAATGCAA` w4 s2→[0.8112781244591328,2.0,1.5]; corner cases `windowSize>length`→empty profile, `==length`→single value, case-insensitive, bounds `0≤H≤log₂k≤2`. One non-value-affecting ASSUMPTION: mono-symbol (k=1) alphabet is the implementation's modelling choice (sources define H over any distribution; the k-mer/block generalization is the separate SEQ-COMPLEX-KMER-001 k-entropy). Wrote source page [[seq-entropy-profile-001-evidence]]; enriched concept [[windowed-sequence-complexity-profile]] (frontmatter +SEQ-ENTROPY-PROFILE-001 source + source_commit→HEAD + new graph edge; intro pointer + new "Second entry point" body section); updated hub [[algorithm-validation-evidence]] (frontmatter sources +SEQ-ENTROPY-PROFILE-001 + source_commit→HEAD + body evidence-link list) and wiki/index.md (Sources entry + note on the windowed-profile Concepts line). Graph: one typed edge on the enriched concept, `relates_to concept:test-unit-registry` (standard per-unit registration for the second entry point); no source-backed concept-to-concept edge (the source defines the entropy-profile method, it asserts no typed relation to siblings — those links are wiki-navigational). Contradictions: none — Shannon 1948, the Wikipedia information-theory definition, and the IntechOpen DNA application are mutually consistent; oracles follow directly. Follow-ups: Shannon character-entropy over a whole sequence (non-windowed) remains the last un-ingested named entropy member; the standalone linguistic-complexity unit SEQ-COMPLEX-001 still not separately ingested.
   graph: +1 node, +1 typed edge (relates_to test-unit-registry)
- 2026-07-10 — ingest docs/Evidence/SEQ-GC-ANALYSIS-001-Evidence.md (Comprehensive GC Analysis — the composite `GcAnalysisResult` bundling whole-sequence OverallGcContent `(G+C)/(A+T+G+C)×100` + OverallGcSkew `(G−C)/(G+C)` + OverallAtSkew `(A−T)/(A+T)` with a sliding-window GC%/GC-skew profile AND the population variance `σ²=Σ(xᵢ−μ)²/N` of each windowed series). CONTEXT check: read [[base-composition]] (SEQ-COMPOSITION-001, GC content), [[nucleotide-composition-skew]] (SEQ-ATSKEW-001, GC/AT skew — which flagged the `gc-skew` member as not-yet-ingested), [[windowed-sequence-complexity-profile]] (entropy/LC profile), [[centromere-analysis]] (GC-variability heuristic), and searched wiki for gc-content/gc-skew/windowed/isochore. The three overall scalars are ALREADY covered (GC content → base-composition; GC/AT skew → nucleotide-composition-skew) so they were NOT re-derived — cross-linked instead. GENUINELY DISTINCT layer = the **windowed GC-content/GC-skew profile** + the **compositional (population) variance** summaries (GcContentVariance/GcSkewVariance), the isochore-style GC-heterogeneity signal named in the ingest guidance → created ONE dedicated concept [[windowed-gc-profile-and-variance]], cross-linked as the composition sibling of the entropy/complexity [[windowed-sequence-complexity-profile]] and the explicit-profile version of the [[centromere-analysis]] GC-variability heuristic. Formulas fully sourced: Wikipedia GC-content (Madigan & Martinko 2003 *Brock* `GC%=(G+C)/(A+T+G+C)×100%`) + Wikipedia GC skew (Lobry 1996 *Mol Biol Evol* 13(5):660 + Grigoriev 1998 *NAR* 26(10):2286 cumulative-skew diagram; spectrum −1…+1, sign switch = replication origin/terminus) + Biopython `Bio.SeqUtils` (`GC_skew` "multiple windows along the sequence", zero-division ⇒ 0, ambiguous bases ignored; `gc_fraction` `[0,1]` = percentage ÷100) + Cuemath population variance (`σ²=Σ(xᵢ−μ)²/n`, worked `{12,13,12,14,19}`→μ=14/Σ=34/**6.8** anchor). Oracles `GGGCCAT`→GC% 71.42857…/GC-skew 0.2/AT-skew 0.0; windowed `GGCC` w2s2→windows GG(skew+1,GC%100)/CC(skew−1,GC%100)→`GcSkewVariance`=((1)²+(−1)²)/2=**1.0** / `GcContentVariance`=**0.0**; corner cases sequence-shorter-than-window ⇒ empty windowed lists ⇒ window-variances 0 while overall scalars still computed, no-G/C ⇒ skew 0/GC% 0, pure-G window ⇒ +1 / pure-C ⇒ −1, null DnaSequence ⇒ ArgumentNullException. Two documented ASSUMPTIONs, both labelling/estimator (neither correctness-affecting at the formula level): GC content reported as a **percentage ×100** (matching the repo `GcAnalysisResult`/`CalculateGcContent`) not Biopython's `[0,1]` fraction (differ only by 100), and window "variability" = **population variance ÷N** not sample variance ÷(N−1) (windows ARE the population, per Cuemath). Wrote source page [[seq-gc-analysis-001-evidence]] + concept [[windowed-gc-profile-and-variance]]; enriched [[base-composition]] (GC-variability bullet now names the dedicated windowed-GC concept + the ×100 vs `[0,1]` units note; frontmatter +SEQ-GC-ANALYSIS-001 source + source_commit→HEAD) and [[nucleotide-composition-skew]] (GC-skew-member sentence now records SEQ-GC-ANALYSIS-001 computes OverallGcSkew + windowed profile + variance; frontmatter +SEQ-GC-ANALYSIS-001 + source_commit→HEAD); updated hub [[algorithm-validation-evidence]] (frontmatter sources +SEQ-GC-ANALYSIS-001 + source_commit→HEAD + body evidence-link list + concept enumeration entry) and wiki/index.md (Sources + Concepts entries). Graph: three typed edges on the new concept — `relates_to concept:test-unit-registry` (standard per-unit registration) + `relates_to concept:base-composition` (GcAnalysisResult's OverallGcContent IS the GC-content quantity, windowed) + `relates_to concept:nucleotide-composition-skew` (its OverallGcSkew/OverallAtSkew ARE the skew family, windowed) — both concept edges source-backed by the single `GcAnalysisResult` computing all three together. Contradictions: none — GC-content (Brock), GC skew (Lobry/Grigoriev/Biopython), and the Cuemath population-variance anchor agree on formula, range, and estimator; the ×100-vs-fraction and ÷N-vs-÷(N−1) items are documented convention choices, not conflicts. Follow-ups: the standalone `gc-skew` unit flagged on [[nucleotide-composition-skew]] is still separately un-ingested (this composite realizes GC skew but is a distinct registry entry); a future dedicated windowed-GC or isochore unit would enrich [[windowed-gc-profile-and-variance]] rather than add a page.
   graph: +2 nodes, +3 typed edges (relates_to test-unit-registry; relates_to base-composition; relates_to nucleotide-composition-skew)
- 2026-07-10 — ingest docs/Evidence/SEQ-GC-PROFILE-001-Evidence.md (GC Content Profile — sliding-window GC% = `(G+C)/(A+T+G+C)×100` per fully-contained window; GC% only, no skew profile, no variance). CONTEXT: this is the STANDALONE windowed-GC%-content channel that the composite `GcAnalysisResult` (SEQ-GC-ANALYSIS-001) wraps — same measure + same window geometry, so per the ingest guidance ENRICHED the existing concept [[windowed-gc-profile-and-variance]] rather than creating a new one (mirrors the SEQ-ENTROPY-PROFILE-001 → windowed-sequence-complexity-profile precedent). Added a "Standalone entry point" section + intro pointer + graph edges. Sources: Wikipedia "GC-content" (`GC%=(G+C)/(A+T+G+C)×100%`, four-base denominator, undefined at A+T+G+C=0) + Biopython `Bio.SeqUtils.gc_fraction` ([0,1]×100=percentage; default `ambiguous="remove"` drops N from length; doctests `ACTGN` remove→0.50 / ignore→0.40 / RNA `GGAUCUUCGGAUCU`→0.50). Hand-derived window oracles (×100) `GGGG`→100 / `AAAA`→0 / `ATGC`→50 / `GGGA`→75 / `GCAT`→50 / `GGAN` (N excluded)→66.66666666666666; window count `⌊(n−w)/step⌋+1`, offsets 0/step/2·step, RNA U=non-GC, bounds `[0,100]`, `windowSize>length`/null/empty→empty profile, case-insensitive. One open ASSUMPTION (all-N window ⇒ A+T+G+C=0 div-by-zero → GC% 0, mirroring sibling `GcSkewCalculator`; not dictated by sources) — N-exclusion denominator is source-backed (`remove` mode), not assumed; no source contradictions. Wrote source page [[seq-gc-profile-001-evidence]]; enriched concept [[windowed-gc-profile-and-variance]] (frontmatter +SEQ-GC-PROFILE-001 source, source_commit→599fc94, +3 graph edges); updated hub [[algorithm-validation-evidence]] (frontmatter sources +SEQ-GC-PROFILE-001) and wiki/index.md.
   graph: +1 node, +3 typed edges (seq-gc-profile-001-evidence relates_to test-unit-registry; relates_to base-composition; + shared windowed-gc-profile-and-variance node)
- 2026-07-10 — ingest docs/Evidence/SEQ-HYDRO-001-Evidence.md (Hydrophobicity Analysis — the whole-sequence **Kyte-Doolittle GRAVY index** `GRAVY=Σ(kd value)/length` + the **unweighted sliding-window hydropathy profile** of exactly `N−W+1` window means). CONTEXT check: searched wiki for hydrophob/hydropath/kyte/doolittle/gravy/transmembrane and read the two closest concepts — [[transmembrane-helix-prediction]] (PROTMOTIF-TM-001, the ProteinMotif segment-calling unit) and [[intrinsic-disorder-prediction-top-idp]] (its `CalculateHydropathy` utility) — plus [[base-composition]] and the SEQ-\* composition siblings. GRAVY (a whole-sequence hydropathy **scalar**) and the raw hydropathy **profile** are NOT represented anywhere: the transmembrane concept only carries the W=19/threshold-1.6 **segment-calling** application of the same profile and explicitly frames itself as a distinct ProteinMotif unit. GENUINELY DISTINCT (different output = scalar+raw profile vs membrane segments; different family = SEQ-\* sequence-statistics vs ProteinMotif) → created ONE dedicated concept [[hydrophobicity-gravy-and-profile]] as the protein-property member of the SEQ-\* family (amino-acid analogue of the nucleotide [[base-composition]]), cross-linked as the profile that the segment-calling [[transmembrane-helix-prediction]] thresholds and that shares the kd scale with the disorder anchor's `CalculateHydropathy`. Formulas fully sourced: Biopython `ProtParamData.py` kd scale (20 residues I 4.5…R −4.5, attributed to Kyte & Doolittle 1982) + `ProtParam.py` (`gravy()`=total/length, `protein_scale` loop `range(N−W+1)`, default `edge=1.0`=unweighted mean) + Expasy ProtParam (GRAVY = Σ hydropathy ÷ residues, positive⇒hydrophobic) + GCAT/Davidson (window 9 surface / window 19 transmembrane "peaks with scores greater than 1.6") + alakazam (CRAN) `gravy` (KD default, 1982 primary). Hand-derived oracles `A`→1.8 / `AG`→0.7 / `FLIV`→3.825 / `RKDE`→−3.85; window-3 `FLIV`→[3.7, 4.1666666667] / `AG` (W>N)→empty; case-insensitive, empty/null → GRAVY 0/empty profile. One documented DEVIATION (unknown residues **skipped** — GRAVY divides by recognized count, profile treats as 0 — whereas Biopython `gravy()` raises `KeyError`; Kyte-Doolittle 1982 + Expasy define values only for the 20 standard residues and are silent on ambiguity codes/gaps, so neither rule is source-mandated — an API-shape/robustness choice, every in-alphabet value stays exactly source-conformant, algorithm doc §5.4). Wrote source page [[seq-hydro-001-evidence]] + concept [[hydrophobicity-gravy-and-profile]]; enriched [[transmembrane-helix-prediction]] (reciprocal cross-link naming this unit as the segment-calling layer over the SEQ-HYDRO-001 profile + frontmatter +SEQ-HYDRO-001 source + source_commit→HEAD); updated hub [[algorithm-validation-evidence]] (frontmatter sources +SEQ-HYDRO-001 + source_commit→HEAD + body source-link list [also added the previously-missing seq-gc-profile-001-evidence] + concept enumeration entry) and wiki/index.md (Sources + Concepts). Graph: two typed edges on the new concept — `relates_to concept:test-unit-registry` (standard per-unit registration) + `relates_to concept:transmembrane-helix-prediction` (same kd scale + same `edge=1.0` window mean; SEQ-HYDRO's GCAT/Davidson source supplies the TM unit's W=19/1.6 defaults, which threshold this profile). Contradictions: none — Biopython, Expasy, GCAT/Davidson and alakazam agree on the scale, the GRAVY formula and the unweighted-window-mean rule; the unknown-residue item is a documented convention choice, not a conflict. Follow-ups: SEQ-STATS-001 (which SEQ-COMPOSITION-001 consolidates) and the standalone `gc-skew` unit remain separately un-ingested; a future protein molecular-weight / pI QC unit (bio-qc composition family) would be a sibling concept.
   graph: +2 nodes, +2 typed edges (relates_to test-unit-registry; relates_to transmembrane-helix-prediction)
- 2026-07-10 — ingest docs/Evidence/SEQ-MW-001-Evidence.md (Molecular Weight Calculation — the average-isotopic **molecular mass in Daltons** of a **protein, DNA or RNA** sequence via one shared polymerization formula `weight = Σ (monomer mass) − (len−1)·water`). CONTEXT check: read the ingest-named siblings [[base-composition]] (SEQ-COMPOSITION-001, the nucleotide/residue count tally) and [[hydrophobicity-gravy-and-profile]] (SEQ-HYDRO-001, the protein GRAVY scalar), and searched wiki/ for molecular-weight / mass / Dalton / monophosphate / residue-mass / protein-property (only log/backlog mentions — no existing MW/mass concept). DECISION: MW is a genuinely distinct, well-defined biochemical property and a whole-sequence **scalar** spanning all three molecule types (protein+DNA+RNA) — NOT reducible to base-composition (a count) or hydrophobicity (protein-only average) — so created ONE dedicated concept [[molecular-weight]] as the **mass member of the SEQ-\* sequence-statistics family**, cross-linked as the mass-weighted view of the [[base-composition]] monomer tally and a whole-sequence-scalar sibling of the [[hydrophobicity-gravy-and-profile]] GRAVY scalar. Formula fully sourced: Biopython `Bio.SeqUtils.molecular_weight` (single-strand `weight = sum(weight_table[x]) − (len−1)·water`, one water removed per bond — peptide for protein, phosphodiester for nucleic acid; `water = 18.0153` avg / `18.010565` monoisotopic; single monomer ⇒ free-monomer mass, circular ⇒ `weight −= water`, double-strand adds the complement strand and is an **error for protein**; "only unambiguous letters", nucleotides assumed to carry a 5' phosphate) + Biopython `Bio.Data.IUPACData` average mass tables ("from PubChem": protein *free-amino-acid* A 89.0932/G 75.0666/W 204.2252/…; DNA *monophosphate* A 331.2218/C 307.1971/G 347.2212/T 322.2085; RNA A 347.2212/C 323.1965/G 363.2206/U 324.1813) + Expasy Compute pI/Mw (protein Mw = "addition of average isotopic masses of amino acids ... and ... one water molecule", Da — algebraically identical to `Σ free-amino-acid − (len−1)·water`) + Expasy ProtParam (delegates "as in Compute pI/Mw") + Expasy FindMod (residue-mass table Ala 71.0788… + H₂O 18.01524). Re-derived oracles (average tables, water 18.0153): `AGC`→protein **249.29** / DNA **949.61** / RNA **997.61** (Biopython docstring), single monomer (zero bonds ⇒ free-monomer mass) `G`→**75.0666** / `A` DNA→**331.2218** / `A` RNA→**347.2212**; empty/null→0; two-monomer subtracts exactly one water (bond-count invariant). One documented DEVIATION (unknown/ambiguous symbols **skipped** — contribute no mass and no bond, so the reported mass reflects only recognized monomers and no invented "average" mass is used — whereas Biopython *rejects* unknown letters; this is the SEQ-\* skip-unknown / non-throwing convention shared with [[hydrophobicity-gravy-and-profile]]'s GRAVY and [[base-composition]]'s CountOther routing, an API-shape/robustness choice, every in-alphabet value stays exactly source-conformant) + two API-shape ASSUMPTIONs (ToUpperInvariant case-fold; the skip resolution). Wrote source page [[seq-mw-001-evidence]] + concept [[molecular-weight]]; enriched [[base-composition]] (new "Molecular weight" sibling bullet = the mass-weighted view of the same tally; frontmatter +SEQ-MW-001 source + source_commit→HEAD) and [[hydrophobicity-gravy-and-profile]] (intro now names [[molecular-weight]] as the mass sibling summing the same 20 residues' Daltons; frontmatter +SEQ-MW-001 + source_commit→HEAD); updated hub [[algorithm-validation-evidence]] (frontmatter sources +SEQ-MW-001 + source_commit→HEAD + body source-link list + concept enumeration entry) and wiki/index.md (Sources + Concepts entries). Graph: two typed edges on the new concept — `relates_to concept:test-unit-registry` (standard per-unit registration) + `relates_to concept:base-composition` (the single-strand mass formula is the per-monomer mass sum over the same {A,T,G,C,U}+amino-acid alphabets base-composition counts — MW is the mass-weighted composition, sharing the case-fold + skip-unknown contract). Contradictions: none — Expasy and Biopython agree on the formula (protein `Σ residue + water` = `Σ free-amino-acid − (len−1)·water`) and on the average-mass tables to rounding; the skip-vs-reject item is a documented convention choice, not a conflict. Follow-ups: MW pairs with **pI / isoelectric point** in the Expasy Compute pI/Mw source but pI/charge is a **separate calculation not in this unit** — a future SEQ-pI/charge QC unit (bio-qc composition family) would be a sibling concept; the standalone `gc-skew` unit and SEQ-STATS-001 remain separately un-ingested.
   graph: +2 nodes, +2 typed edges (relates_to test-unit-registry; relates_to base-composition)
- 2026-07-10 — ingest docs/Evidence/SEQ-PI-001-Evidence.md (Isoelectric Point (pI) — the **pH at which a protein's net charge is zero**, found by **bisecting the Henderson–Hasselbalch net-charge function over [0,14]** to ±0.01; charge `Σ_basic +1/(1+10^(pH−pKa)) + Σ_acidic −1/(1+10^(pKa−pH))`, Moore 1985/Peptides `charge_pI.cpp`, basic {R,K,H,N-term} / acidic {D,E,C,Y,C-term}). Genuinely distinct SEQ-\* protein-property scalar (charge/pH, not mass or hydropathy) → NEW concept [[isoelectric-point]], the ExPASy Compute pI/Mw partner of [[molecular-weight]] (updated its Scope line to link the new unit). EMBOSS `Epk.dat` scale adopted (single-pKa-per-residue, not Bjellqvist); composition-only/order-independent ("no electrostatic interactions"). Sources EMBOSS `iep` + Peptides `charge_pI.cpp`/`charge()` + seqinr `computePI` + ExPASy Compute pI/Mw + Bjellqvist 1993. Oracles `A`/`AG`→6.10, `D`→3.75, `K`→9.70, `DDDD`→3.23, `KKKK`→11.27, `FLPVLAGLTPSIVPKLVCLLTKKC`→9.67 (net charge 3.037398/2.914112/0.7184524 @ pH 5/7/9), `ACDEFGHIKLMNPQRSTVWY`→7.36; pI∈[0,14], two assumptions (empty/null→7.0; scale=EMBOSS), no source contradictions.
   graph: +2 nodes, +2 typed edges (isoelectric-point relates_to test-unit-registry; relates_to molecular-weight)
- 2026-07-10 — ingest docs/Evidence/SEQ-REPLICATION-001-Evidence.md (Replication Origin Prediction — the classic ori-finder via the **cumulative GC-skew minimum**: integrate an *integer* running skew (G:+1, C:−1, A/T:0, `Skew_0=0`, `|Genome|+1` prefix values over `i∈[0,n]`) and read `PredictedOrigin` off the global **minimum**, `PredictedTerminus` off the **maximum**). CONTEXT check: read the two ingest-named siblings [[nucleotide-composition-skew]] (the scalar `(G−C)/(G+C)` skew family, which explicitly *anticipated* a cumulative-skew origin locator as a future application) and [[windowed-gc-profile-and-variance]] (which notes the cumulative version locates origin/terminus, Grigoriev 1998), plus grepped wiki/ for replication/origin/terminus/cumulative/gc-skew and read [[base-composition]]. DECISION: the *locating* algorithm is GENUINELY DISTINCT from those scalar/windowed skew statistics — it adds **cumulative integration + argmin/argmax extremum search** and uses the **integer** cumulative form of Rosalind's Minimum Skew Problem (BA1F), not the ratio-normalised scalar — so created ONE dedicated concept [[replication-origin-cumulative-skew]] as the skew family's *localisation* member, cross-linked from both anticipating siblings (each now names this concept as the realised locator). Sources: Rosalind BA1F Minimum Skew Problem (rank 3, canonical worked example — 100-nt sample min −4 @ positions 53,97, re-derived in-session exactly → `PredictedOrigin=53` first-minimizer tie-break) + Grigoriev 1998 *Nucleic Acids Res* 26(10):2286 (rank 1 — cumulative diagram min=origin/max=terminus, extrema ≈ half chromosome apart, leading strand G-rich) + Lobry 1996 (PMID 8676740 — skew sign flips at origin/terminus) + Wikipedia "GC skew" (max=terminal/min=origin). Diagram oracles `CCGGGG`→`0,−1,−2,−1,0,+1,+2` (min −2@2/max +2@6), `GGGCCC`→`0,+1,+2,+3,+2,+1,0`, `AATT`→flat. One documented ASSUMPTION — `IsSignificant` redefined threshold-free as `max>min` (the invented `amplitude>count×0.01` constant was **removed** as untraceable); flat diagram (no/balanced G/C) ⇒ origin=terminus=0/`IsSignificant=false`; no source contradictions. Wrote source page [[seq-replication-001-evidence]] + concept [[replication-origin-cumulative-skew]]; enriched [[nucleotide-composition-skew]] and [[windowed-gc-profile-and-variance]] (each: cross-link naming the new locator concept + frontmatter +SEQ-REPLICATION-001 source + source_commit→HEAD); updated hub [[algorithm-validation-evidence]] (frontmatter sources +SEQ-REPLICATION-001) and wiki/index.md (Sources + Concepts entries). Follow-ups: the standalone `gc-skew` unit and SEQ-STATS-001 remain separately un-ingested.
   graph: +1 node, +2 typed edges (replication-origin-cumulative-skew relates_to nucleotide-composition-skew; relates_to test-unit-registry)
- 2026-07-10 — ingest docs/Evidence/SEQ-RNACOMP-001-Evidence.md (RNA-specific per-base complement — `GetRnaComplementBase`, the IUPAC-complete RNA complement map). **CONTEXT CORRECTION:** the `RNACOMP` slug is RNA **complement**, NOT RNA composition — the ingest brief mis-described it as the U-alphabet composition tally and pointed at [[base-composition]]; that page was left untouched. The unit is RNA complement chemistry (A→U/U→A/G→C/C→G, T→A treating DNA T as U; reciprocal R↔Y/M↔K/D↔H/B↔V; self-complementary W/S/X/N; non-IUPAC pass-through; recognized bases uppercased — casing-only divergence from case-preserving Biopython). **NO new concept** — enriched [[rna-base-pairing]], which already owns the RNA base complement (`GetComplement`, RNA-PAIR-001): added a *SEQ-family full-IUPAC RNA complement* subsection (this SEQ-\* surface maps EVERY IUPAC code, vs the canonical+N/R/Y subset the RnaSecondaryStructure copy documents; RNA sibling of DNA `GetComplementBase`, SEQ-COMP-001 not yet ingested). Wrote source [[seq-rnacomp-001-evidence]] with the full-alphabet + ambiguity-code oracles. Updated hub [[algorithm-validation-evidence]] (frontmatter sources + source_commit→HEAD; body already lists rna-base-pairing). Sources: Biopython `IUPACData.py`/`Seq.py`/API docs + bioinformatics.org SMS + NC-IUB 1984 (Cornish-Bowden 1985); no source contradictions. index.md: +1 source (no new concept).
   graph: +1 node, +1 typed edge (relates_to test-unit-registry on rna-base-pairing, source seq-rnacomp-001-evidence)
- 2026-07-10 — ingest docs/Evidence/SEQ-SECSTRUCT-001-Evidence.md (Protein Secondary Structure Prediction — **Chou-Fasman conformational propensities** Pα/Pβ/Pt evaluated as a generic **sliding-window mean-propensity profile**). CONTEXT check: the "SEQ-SECSTRUCT" slug is **protein** secondary structure (Chou-Fasman), NOT RNA — searched wiki/concepts for chou/fasman/secondary-structure/helix/sheet/GOR/propensity and read the RNA-structure family (already covered by [[rna-secondary-structure-prediction]] etc., left untouched) and the protein siblings [[hydrophobicity-gravy-and-profile]] (SEQ-HYDRO-001) + [[transmembrane-helix-prediction]]; no existing protein-secondary-structure/Chou-Fasman concept. DECISION: created ONE dedicated concept [[protein-secondary-structure-chou-fasman]] as the conformational-propensity member of the SEQ-\* sequence-statistics family — same sliding-window-mean machinery as the hydropathy profile (N−W+1 unweighted window means, W>N→empty, unknown residues skip-and-excluded, case-insensitive), different per-residue value table (Pα/Pβ/Pt vs Kyte-Doolittle). The unit under test is the generic **profile, not** the classic 4-of-6-helix/3-of-5-sheet nucleation-extension state machine, so default W=7 is an API convenience not a CF constant. Sources: Wikipedia Chou–Fasman (rank 4) + Kelley lecture (rank 4) + CSB|SJU Jakubowski (rank 4, Lys Pα 1.16) + Przytycka NCBI lecture (rank 4, Lys 1.14) + ravihansa3000/ChouFasman reference impl (rank 3, integer table Lys 114) + BMC PMC1780123 (rank 1). Oracles single-residue window = that residue's tuple, `AE`→helix 1.465/sheet 0.60/turn 0.70, `AEV`→helix 1.330, `AXE` W3 averages only A+E. Three assumptions (profile-not-state-machine / W=7 API default; skip-unknown; the one contested value **Lys Pα resolved to 1.14** — two independent sources incl. reference impl vs single CSB|SJU 1.16), ~50–60% Q3 / 29-protein reliability caveat, no contradictions on the 19 uncontested residues. Wrote source [[seq-secstruct-001-evidence]] + concept [[protein-secondary-structure-chou-fasman]]; updated hub [[algorithm-validation-evidence]] (frontmatter sources +SEQ-SECSTRUCT-001; body source-link list + concept enumeration entry) and wiki/index.md (Sources + Concepts entries). Graph: two typed edges on the new concept — relates_to test-unit-registry (standard per-unit registration) + relates_to hydrophobicity-gravy-and-profile (shared SEQ-\* sliding-window mean-propensity mechanism). Follow-ups: a future segment-calling unit that thresholds this profile through the nucleation/extension state machine would be a distinct sibling; SEQ-STATS-001 remains un-ingested.
   graph: +2 nodes, +2 typed edges (relates_to test-unit-registry; relates_to hydrophobicity-gravy-and-profile)
- 2026-07-10 — ingest docs/Evidence/SEQ-STATS-001-Evidence.md (**the ORIGINAL sequence-composition-statistics umbrella** — nucleotide composition A/T/G/C/U/N/Other + Length, GC content `(G+C)/(A+T+G+C+U)`, GC skew `(G−C)/(G+C)`, AT skew `(A−T)/(A+T)`; coverage recs also exercise `SummarizeNucleotideSequence` aggregator + `CalculateAminoAcidComposition`). CONTEXT check per brief: read [[base-composition]], [[nucleotide-composition-skew]], [[molecular-weight]], [[dinucleotide-relative-abundance]], [[codon-usage-comparison]] and the sibling source [[seq-composition-001-evidence]] (whose Change History records that SEQ-COMPOSITION-001 is a **duplicate/consolidated Registry entry over the two composition methods already delivered under SEQ-STATS-001**, TestSpec §7). DECISION: **NO new concept** — every method this umbrella exercises is already synthesized on an existing page (nucleotide composition + GC content → [[base-composition]]; GC skew + AT skew → [[nucleotide-composition-skew]]; amino-acid composition → [[base-composition]]). Wrote source-summary [[seq-stats-001-evidence]] as the connective tissue tying the piecemeal-ingested composition family to its original home. ENRICHED: [[base-composition]] (added "Original home" paragraph + frontmatter +SEQ-STATS-001 source + source_commit→HEAD + a second relates_to test-unit-registry graph edge from seq-stats-001-evidence), [[nucleotide-composition-skew]] (noted both skews first delivered together under SEQ-STATS-001 before AT skew split into SEQ-ATSKEW-001 and both re-bundled by SEQ-GC-ANALYSIS-001; frontmatter +source + source_commit→HEAD), [[seq-composition-001-evidence]] (resolved its "(SEQ-STATS-001 is not yet ingested)" note → cross-link to the now-ingested umbrella). Updated hub [[algorithm-validation-evidence]] (frontmatter sources +SEQ-STATS-001 + source_commit→HEAD; body source-link enumeration) and wiki/index.md (+1 Sources entry). Sources: Biopython `gc_fraction`/`GC_skew` (rank 3; `gc/length`∈[0,1], empty⇒0, case-insensitive, `ZeroDivisionError`⇒skew 0.0) + Wikipedia "GC skew" (rank 4) + Lobry 1996 *MBE* 13(5):660 (rank 1, primary). Oracles `ATGC`/`GGGC`/`AAAT`/`GCCC`/`AAUUGGCC` (GC content + both skews incl. negative + zero-denominator cases); one ASSUMPTION (degenerate IUPAC codes excluded from GC/AT totals, agrees with Biopython over {A,T,G,C,U}), no source contradictions. Resolves the "SEQ-STATS-001 remains un-ingested" follow-up flagged by SEQ-COMPOSITION/SEQ-ATSKEW/SEQ-REPLICATION/SEQ-SECSTRUCT ingests.
   graph: +1 node (source seq-stats-001-evidence), +1 typed edge (base-composition relates_to test-unit-registry, source seq-stats-001-evidence)
- 2026-07-10 — ingest docs/Evidence/SEQ-SUMMARY-001-Evidence.md (**`SummarizeNucleotideSequence` — the per-sequence summary record**: bundles Length, nucleotide Composition (A/T/G/C/U/N), GcContent `(G+C)/(A+T+G+C+U)`, Shannon Entropy `−Σpᵢlog₂pᵢ` bits, linguistic Complexity, and MeltingTemperature (Wallace len<14 / GC-Marmur-Doty len≥14) into one object). CONTEXT check per brief: read [[base-composition]], [[nucleotide-composition-skew]], [[molecular-weight]], [[seq-stats-001-evidence]], [[seq-entropy-profile-001-evidence]], [[windowed-sequence-complexity-profile]] and [[primer-dimer-thermodynamics-tm]]; grepped wiki/ for SEQ-SUMMARY / SummarizeNucleotideSequence / SEQ-TM-001 / Wallace / Marmur / linguistic. DECISION: **NO new concept** — SEQ-SUMMARY-001 is a **pure field-by-field aggregation** whose contract is "each summary field equals the value its canonical per-metric method returns on the identical input", and every aggregated method is already synthesized on an existing page (composition/GC → [[base-composition]]; entropy → SEQ-ENTROPY-PROFILE-001/[[windowed-sequence-complexity-profile]]; linguistic complexity → [[windowed-sequence-complexity-profile]]; Tm → `CalculateMeltingTemperature` on [[primer-dimer-thermodynamics-tm]], the legacy Wallace/Marmur default; SEQ-TM-001 not yet separately ingested). Wrote source-summary [[seq-summary-001-evidence]]. ENRICHED (cross-reference correction): [[base-composition]] and [[seq-stats-001-evidence]] previously called `SummarizeNucleotideSequence` a "thin aggregation wrapper re-exposing the same GC content / counts" — this artifact shows it is broader (also carries Shannon entropy + linguistic complexity + Tm), so both now describe it as a full per-sequence statistics record linking to the new page (not a contradiction — the earlier pages named only the composition facet). [[base-composition]] frontmatter +SEQ-SUMMARY-001 source + source_commit→HEAD + a relates_to test-unit-registry graph edge (source seq-summary-001-evidence). Updated hub [[algorithm-validation-evidence]] (frontmatter sources +SEQ-SUMMARY-001 + source_commit→HEAD; body source-link enumeration) and wiki/index.md (+1 Sources entry). Sources: Biopython `gc_fraction` (rank 3) + `MeltingTemp` (`Tm_Wallace` `4(G+C)+2(A+T)` / `Tm_GC` Marmur-Doty `69.3+0.41·%GC−650/N`, rank 3) + Wikipedia "Entropy (information theory)" (rank 4, base-2⇒bits, max log₂n) + Wikipedia "Linguistic sequence complexity" (rank 4, Trifonov 1990, `C=∏Uᵢ`, 0<C<1). Oracles `ATGCATGC`→GcContent 0.5/Entropy 2.0/Tm 24.0/Complexity 0.83968253968…, 16-mer `ATGCATGCATGCATGC`→Tm 43.375 (GC branch, len≥14 `64.9+41·(8−16.4)/16`); empty/null⇒zero summary (Length/GcContent/Entropy/Complexity/Tm all 0), case-insensitive. One ASSUMPTION — the Tm formula-selection threshold `length<14` (= SEQ-TM-001 `ThermoConstants.WallaceMaxLength`) is non-correctness-affecting for the summary because its contract is equality with `CalculateMeltingTemperature` on the same flag. No source contradictions. Follow-up: the standalone Tm unit **SEQ-TM-001** (`CalculateMeltingTemperature` Wallace/Marmur-Doty length-dispatch) remains separately un-ingested — its own scalar Tm concept would be the natural home for the Wallace/Marmur formulas currently only sketched on [[primer-dimer-thermodynamics-tm]] as the legacy default.
   graph: +1 node (source seq-summary-001-evidence), +1 typed edge (base-composition relates_to test-unit-registry, source seq-summary-001-evidence)
- 2026-07-10 — ingest docs/Evidence/SEQ-THERMO-001-Evidence.md (**DNA Duplex Thermodynamics — the full nearest-neighbor ΔH°/ΔS°/ΔG°/Tm 4-tuple** of a Watson-Crick DNA duplex, a verbatim port of Biopython `Bio.SeqUtils.MeltingTemp.Tm_NN` over the **DNA_NN3 / Allawi & SantaLucia 1997** parameter set). CONTEXT check per brief: read [[primer-dimer-thermodynamics-tm]] (the PRIMER-family NN Tm engine, PRIMER-TM-001), [[base-composition]] + [[seq-summary-001-evidence]] (the SEQ family — `SummarizeNucleotideSequence` bundles the **legacy Wallace/Marmur** scalar Tm, SEQ-TM-001), and [[seq-stats-001-evidence]]; grepped wiki/ for tm/melting/thermodynamic/nearest-neighbor/SantaLucia/duplex/ΔG. DECISION: **NEW concept warranted** — this is a genuinely distinct **full-NN ΔH°/ΔS°/ΔG° engine** returning the whole thermodynamic tuple, and it is **neither** of the two Tm surfaces already in the wiki: (a) NOT the legacy Wallace/Marmur `%GC` rule of thumb (SEQ-TM-001, no thermodynamic model), and (b) NOT the 2004-unified PCR-primer engine [[primer-dimer-thermodynamics-tm]] — it uses the **1997 DNA_NN3** table with **per-terminus `init_A/T` (2.3,4.1) / `init_G/C` (0.1,−2.8)** counted off the two terminal bases, vs the primer engine's 2004 unified table with a fixed duplex-init (+0.2/−5.7) + terminal-A·T penalty (+2.2/+6.9). Created [[dna-duplex-nearest-neighbor-thermodynamics]] as the SEQ-\* full-tuple thermodynamics member. The bimolecular Tm `=(1000·ΔH°)/(ΔS°+R·ln k)−273.15` (R=1.987, `k=(dnac1−dnac2/2)·1e−9`; default 25 nM each ⇒ k=C_T/4 = the F=4 non-self-complementary factor) and the **salt correction `0.368·(len−1)·ln[Na⁺]`** are the SAME SantaLucia-1998 equations the primer engine uses (N/2 NN stacks = len−1). ΔG°₃₇ `= ΔH° − 310.15·ΔS°/1000`. Sources: Biopython `MeltingTemp.Tm_NN` (rank 3, port target; DNA_NN3 table retrieved verbatim, two-terminus init `ends=seq[0]+seq[-1]`, docstring `CGTTCCAAAGATGTGGGCATGAGCTTAC`→**60.32 °C** at 25 nM / 50 mM Na) + MELTING 5 User Guide (Dumousseau 2012, rank 3; Tm eq §4.2, F-factor 1/4 §4.3, default model "all97") + Wikipedia Nucleic acid thermodynamics (rank 4, cites SantaLucia 1998; per-terminus init + F corroboration). Oracles: `GCGC` (Na 0.05 M, C_T 250 nM) → ΔH −30.0 / ΔS −84.91 / ΔG −3.67 / **Tm −18.6 °C**; `ATCG` (one A/T + one G/C end) → ΔH −23.6 / ΔS −71.81; **length<2 ⇒ `(0,0,0,0)`** (the one API/edge-case convention), case-insensitive (`ToUpperInvariant`), non-self-complementary F=4 only. Zero correctness-affecting assumptions (every constant source-backed), **no source contradictions** (Biopython, MELTING 5, Wikipedia agree on the DNA_NN3 table + two-terminus init + Tm/salt equations). Wrote source [[seq-thermo-001-evidence]] + concept [[dna-duplex-nearest-neighbor-thermodynamics]]; ENRICHED [[primer-dimer-thermodynamics-tm]] (added a "third DNA-duplex Tm surface" cross-reference contrasting the DNA_NN3-1997 vs unified-2004 parameter vintage + init bookkeeping, noting the identical salt term; frontmatter +SEQ-THERMO-001 source + source_commit→HEAD); updated hub [[algorithm-validation-evidence]] (frontmatter sources +SEQ-THERMO-001 + source_commit→HEAD; body source-link enumeration) and wiki/index.md (Sources + Concepts entries). Graph: two typed edges on the new concept — `relates_to test-unit-registry` (standard per-unit registration) + `relates_to primer-dimer-thermodynamics-tm` (both are SantaLucia bimolecular NN Tm engines sharing the Tm+salt equations, differing only in parameter vintage + init bookkeeping). Follow-up: the standalone Wallace/Marmur scalar Tm unit **SEQ-TM-001** (`CalculateMeltingTemperature`) remains separately un-ingested — its own scalar-Tm concept would be the natural home for the Wallace/Marmur formulas currently only sketched as the legacy default on [[primer-dimer-thermodynamics-tm]].
   graph: +2 nodes, +2 typed edges (dna-duplex-nearest-neighbor-thermodynamics relates_to test-unit-registry; relates_to primer-dimer-thermodynamics-tm)
- 2026-07-10 — ingest docs/Evidence/SEQ-TM-001-Evidence.md (**`CalculateMeltingTemperature` — the standalone scalar melting temperature**, length-dispatched at `ThermoConstants.WallaceMaxLength = 14` between the **Wallace rule** `4(G+C)+2(A+T)` for short oligos <14 nt and the **Marmur-Doty GC formula** `64.9+41·(GC−16.4)/N` for len≥14). CONTEXT check per brief: read [[primer-dimer-thermodynamics-tm]], [[dna-duplex-nearest-neighbor-thermodynamics]] and [[seq-summary-001-evidence]]; grepped wiki/ for SEQ-TM-001/CalculateMeltingTemperature/melting-temperature/Wallace/Marmur; confirmed the backlog expects a `melting-temperature` concept slug for `docs/algorithms/{MolTools,Statistics}/Melting_Temperature.md`. DECISION: **NEW concept warranted** — this is the **third distinct Tm surface** and the **canonical home of the Wallace/Marmur formulas** (previously only sketched as the legacy default on [[primer-dimer-thermodynamics-tm]]), explicitly flagged as the natural home by the SEQ-SUMMARY-001 and SEQ-THERMO-001 ingests. Created [[melting-temperature]] as the legacy **%GC scalar** member of the SEQ-\* family — a rule of thumb carrying **no ΔH°/ΔS°/ΔG°**, no concentration term, no salt correction (distinct from both SantaLucia NN engines). **Consolidation note:** the evidence file is a duplicate/consolidated Registry entry over the two methods already delivered under SEQ-THERMO-001 (`CalculateMeltingTemperature` → this scalar concept + `CalculateThermodynamics` → [[dna-duplex-nearest-neighbor-thermodynamics]]); its NN/`Tm_NN`/SantaLucia-1998 datasets belong to that other concept, so the scalar %GC Tm is this ingest's genuinely new surface. Sources: Biopython `MeltingTemp` `Tm_Wallace` (`4(G+C)+2(A+T)`, "rule of thumb 14–20 nt", `ACGTTGCAATGCCGTA`→48.0) + `Tm_GC` (`A+B·%GC−C/N+salt−D·%mismatch`; valueset-1 Marmur-Doty `69.3+0.41·%GC−650/N`, valueset-2 QuikChange `81.5+…−675/N`, rank 3). Oracles Wallace `ACGTTGCAATGCCGTA`→48.0 / `ATGCATGC`→24.0, Marmur-Doty `GCGCGCGCGCATATATATAT`→51.78 / 16-mer `ATGCATGCATGCATGC`→43.375 (length-dispatch, cross-checked by SEQ-SUMMARY-001); empty/null→0, case-insensitive. One ASSUMPTION (14-nt dispatch boundary = `WallaceMaxLength`, non-correctness-affecting under an explicit `useWallaceRule` flag) + NN-facet default-C_T assumption (250 nM repo vs 50 nM Biopython, formula identical) + repo-vs-Biopython Marmur-Doty constant variant (64.9/−672.4 vs 69.3/−650, same `A+B·%GC−C/N` family); no source contradictions. Wrote source [[seq-tm-001-evidence]] + concept [[melting-temperature]]; ENRICHED [[primer-dimer-thermodynamics-tm]] (legacy-default sketch now links to the dedicated concept; frontmatter +SEQ-TM-001 + source_commit→HEAD), [[dna-duplex-nearest-neighbor-thermodynamics]] ("Relationship to the other Tm surfaces" legacy-Tm bullet now links [[melting-temperature]] + shows the actual repo Marmur-Doty formula; frontmatter +SEQ-TM-001 + source_commit→HEAD), [[seq-summary-001-evidence]] (`MeltingTemperature` row + "not yet ingested" note now point to [[melting-temperature]]; frontmatter +SEQ-TM-001 + source_commit→HEAD). Updated hub [[algorithm-validation-evidence]] (frontmatter sources +SEQ-TM-001 + source_commit→HEAD; body source-link enumeration + concept-enumeration entry) and wiki/index.md (Sources + Concepts entries). Graph: three typed edges on the new concept — relates_to test-unit-registry (standard per-unit registration) + alternative_to dna-duplex-nearest-neighbor-thermodynamics + alternative_to primer-dimer-thermodynamics-tm (the %GC rule of thumb vs the two SantaLucia NN engines). Resolves the SEQ-TM-001 un-ingested follow-up flagged by SEQ-SUMMARY-001 / SEQ-THERMO-001.
   graph: +2 nodes, +3 typed edges (melting-temperature relates_to test-unit-registry; alternative_to dna-duplex-nearest-neighbor-thermodynamics; alternative_to primer-dimer-thermodynamics-tm)
- 2026-07-10 — ingest docs/Evidence/SPLICE-ACCEPTOR-001-Evidence.md (Acceptor Site Detection — 3' splice site prediction). CONTEXT check per brief: searched wiki/concepts for splice/splicing/acceptor/donor/intron/exon/PWM/motif/position-weight — NO existing splicing concept; read [[known-motif-search]] (which names position-weight-matrix scanning as its non-exact counterpart) and the hub. DECISION: **NEW concept warranted** — splicing is a genuinely distinct topic and this is the first of the splicing family (acceptor/donor/branch point). Created [[splice-acceptor-site-prediction]] as the family anchor: the canonical **GU-AG rule** 3'ss with three cis-elements — **branch point** (`yUnAy` A ~18–50 nt upstream), **polypyrimidine tract** (C/U-rich, U2AF65), **AG + context** (`(Yn)NYAG|G`) — and **three scoring surfaces**: default PWM+PPT `FindAcceptorSites` (AcceptorPwm weights from Shapiro & Senapathy 1987, `(score/(count+1)+2)/4` [0,1] heuristic), opt-in branch-point `FindAcceptorBranchPoint` (18–40 nt window, Gao 2008 `yUnAy` conservation y@−3 0.790/U@−2 0.746/A@0 0.923/y@+1 0.751, found ≥ 0.8·3.210; Mercer 2015 corroboration), opt-in **MaxEntScan score3ss** `ScoreAcceptorMaxEnt` (Yeo & Burge 2004 23-nt max-entropy model, MIT-licensed maxentpy 82,560-record tables `Data/maxent_score3.txt`, canonical `...AGgga`→2.89 / 8.19 / −0.08 cross-checks), plus U12 **YCCAC** non-canonical scoring (Hall & Padgett 1994). Corner cases: < 20 nt → empty guard, no-AG → empty, cryptic intronic AG decoys, weak PPT skipped. Wrote source [[splice-acceptor-001-evidence]]. Cross-linked to the PWM/motif family: [[known-motif-search]] (reciprocal — added a one-line pointer from its PWM-branch mention), [[regulatory-element-detection]], [[consensus-sequence]]. Updated hub [[algorithm-validation-evidence]] (frontmatter sources +SPLICE-ACCEPTOR-001 + source_commit→HEAD; body source-link enumeration + concept-enumeration entry) and wiki/index.md (Sources + Concepts entries). Graph: two typed edges on the new concept — relates_to test-unit-registry (standard per-unit registration) + relates_to known-motif-search (medium; the AcceptorPwm scorer is the PWM branch of the degenerate/consensus motif family). Deviations/assumptions: all RESOLVED design decisions (PWM weights verified against Shapiro & Senapathy 1987; normalization = heuristic design choice; U12 YCCAC replaces old fixed 0.6; branch-point + MaxEntScan constants source-traceable). Licence flag (not buried): the bundled MaxEntScan table is the MIT-licensed maxentpy port, NOT the academic-terms Burge-lab Perl models (`Data/maxent_score3.LICENSE.md`). Contradictions: none — encyclopedic + Shapiro/Senapathy + Burge + Yeo/Burge + Gao/Mercer sources mutually consistent. Follow-ups: donor (5' splice site / GU) and branch-point-as-own-unit likely future SPLICE-* family members to enrich this anchor.
   graph: +2 nodes, +2 typed edges (relates_to test-unit-registry, relates_to known-motif-search)
- 2026-07-10 — ingest docs/Evidence/SPLICE-DONOR-001-Evidence.md (Donor (5') Splice Site Detection — the GU/GT dinucleotide at the exon-intron boundary, MAG|GURAGU consensus, MaxEntScan score5ss). CONTEXT check per brief: read the just-ingested sibling anchor [[splice-acceptor-site-prediction]] — confirmed it is scoped specifically to the **3' acceptor** (branch point + PPT + AG) and explicitly anticipates a donor member ("anchor for the splicing family (acceptor / donor / branch point)"). DECISION: **NEW parallel concept warranted** (not enrichment of the acceptor page) — created [[splice-donor-site-prediction]] as the donor/5'ss member. Structural asymmetry captured: the donor is a **single contiguous `MAG|GURAGU` motif** (positions −3..+6; 0-G/+1-U ~100%, −1-G ~80% for U1 snRNP) with **no branch point / PPT**, so only **two** scoring surfaces vs the acceptor's three — a default **IUPAC consensus match-fraction** scorer `FindDonorSites`/`ScoreDonorSite` (plainer than the acceptor's Shapiro/Senapathy PWM+PPT; GC donors self-penalize via the +1-U mismatch 8/9<9/9, so the old 0.7 GC penalty was removed) and opt-in **MaxEntScan score5ss** `ScoreDonorMaxEnt` (Yeo & Burge 2004 9-nt window). Key score5-vs-score3 contrast made explicit on both pages: score5 is a **single 4⁷=16,384-entry table with NO overlapping sub-windows** (GT at 0-based 3..4 removed, rest=`window[0:3]+window[5:9]`), vs score3's 82,560-record overlapping-sub-window factorisation. Canonical cross-checks `cagGTAAGT`→10.86 / `gagGTAAGT`→11.08 / `taaATAAGT`→−0.12. GC-AG (~0.5–1%) + U12 **AT-AC**/`ATATCC` non-canonical donors; guards <9 nt/no-GT/empty→empty. Wrote source [[splice-donor-001-evidence]] + concept [[splice-donor-site-prediction]]. Cross-linked the two splice pages reciprocally (acceptor anchor sentence + "Relation to other units" now name the donor; acceptor frontmatter +SPLICE-DONOR-001 source + source_commit→HEAD). Updated hub [[algorithm-validation-evidence]] (frontmatter sources +SPLICE-DONOR-001 + source_commit→HEAD; body source-link enumeration + concept-enumeration entry) and wiki/index.md (Sources + Concepts entries). Sources: Wikipedia RNA splicing / Spliceosome (rank 4, `MAG|GURAGU` / GU-AG rule / U1 GU binding) + Shapiro & Senapathy 1987 (rank 1, 5'ss PWM position frequencies) + Burge/Tuschl/Sharp 1999 (rank 1, GC-AG / U12 AT-AC / `ATATCC`) + Yeo & Burge 2004 (rank 1, MaxEntScan score5ss) with the MIT-licensed maxentpy 16,384-record table + Alberts MBoC. Licence flag (not buried): bundled `Data/maxent_score5.txt` is the MIT maxentpy port, NOT the academic-terms Burge-lab Perl models (`Data/maxent_score5.LICENSE.md`). All prior assumptions ELIMINATED/RESOLVED (PWM→IUPAC binary consensus weights, normalization→plain match fraction, GC 0.7 penalty removed). Contradictions: none — encyclopedic + Shapiro/Senapathy + Burge + Yeo/Burge mutually consistent. Follow-up: branch-point-as-its-own-unit remains the likely third splicing-family member to enrich the acceptor anchor.
   graph: +2 nodes, +2 typed edges (splice-donor-site-prediction relates_to test-unit-registry; relates_to splice-acceptor-site-prediction)
- 2026-07-10 — ingest docs/Evidence/SPLICE-PREDICT-001-Evidence.md (Gene Structure Prediction (Intron/Exon) — the composite/umbrella of the splicing family). CONTEXT check per brief: read both existing splicing concepts [[splice-donor-site-prediction]] (SPLICE-DONOR-001, GT/GU 5'ss) and [[splice-acceptor-site-prediction]] (SPLICE-ACCEPTOR-001, AG 3'ss + branch point + PPT) plus the hub. DECISION: this source is NOT a redundant re-statement of donor+acceptor — it is a genuinely distinct **composite** that adds intron/exon **pairing**, **exon typing/phase**, and the **spliced sequence** neither boundary unit produces → **NEW concept warranted**: [[gene-structure-prediction-intron-exon]] as the third, integrative splicing-family member. Pipeline: find donors (GT/GU) + acceptors (AG) → pair into introns under the **GT-AG rule** (Breathnach & Chambon 1981 `MAG|GURAGU`/`(Y)nNCAG|G`; Burge 1999 one-donor-one-acceptor), bounded by minIntronLength (default 60) / maxIntronLength, **greedy non-overlapping** selection → exons = the gaps, typed Initial/Internal/Terminal/Single (Gilbert 1978) with phase = `(Σ preceding exon lengths) mod 3` (Alberts 2002) → spliced sequence + overall score = mean of per-intron `(donor.Score+acceptor.Score)/2` in [0,1]. Datasets D1 two-exon 153 nt→1 intron/2 exons/70-nt spliced, D2 single-exon 50 nt→Single, D3 empty→0/0/score 0. Corner cases: 30-bp introns excluded by default min=60 (documented), no sites→single-exon, DNA T≡U. All 4 assumptions RESOLVED (greedy=design decision, phase=trivial, overall-score=definition, **default branch-point 0.3 magic constant removed** → fallback `(donor+acceptor)/2`). Wrote source [[splice-predict-001-evidence]] + concept [[gene-structure-prediction-intron-exon]]. Cross-linked: enriched both boundary concepts' "Relation to other units" (donor + acceptor now name the composite that consumes them; appended SPLICE-PREDICT-001 to their `sources:` lists per the established multi-source pattern). Updated hub [[algorithm-validation-evidence]] (frontmatter sources +SPLICE-PREDICT-001 + source_commit→HEAD; body enumeration entry) and wiki/index.md (Sources + Concepts entries). Sources: Wikipedia Intron/Exon/Gene-structure (rank 4) + Breathnach & Chambon 1981 (rank 1, GT-AG >99%, donor/acceptor consensus) + Shapiro & Senapathy 1987 (rank 1, PWM scoring) + Burge/Tuschl/Sharp 1999 (rank 1, U2/U12/GC-AG + pairing invariant). Contradictions: none — all sources mutually consistent on the GT-AG split-gene model. Follow-up: branch-point-as-its-own-unit still the likely next SPLICE-* member.
   graph: +2 nodes, +3 typed edges (gene-structure-prediction-intron-exon depends_on splice-donor-site-prediction; depends_on splice-acceptor-site-prediction; relates_to test-unit-registry)
- 2026-07-10 — ingest docs/Evidence/SV-BREAKPOINT-001-Evidence.md (Breakpoint Detection from Split (soft-clipped) Reads). CONTEXT check per brief: searched wiki/concepts + index for structural/variant/breakpoint/split-read/discordant/CNV/SV — found NO germline-SV concept; the closest existing pages are the ONCOLOGY read-evidence [[gene-fusion-detection-read-evidence]] (STAR-Fusion/Arriba transcript fusion calling) and [[chromothripsis-inference]] (ONCO-SV-001 copy-number pattern), plus the gene-order [[genome-rearrangement-breakpoint-distance]] — all distinct topics. DECISION: **NEW concept warranted** — this is a genuinely distinct germline structural-variant method (breakpoint localization from CIGAR soft-clip junctions + positional clustering) and the **first of the germline SV family** (StructuralVariantAnalyzer). Created [[breakpoint-detection-split-reads]] as the family anchor. The soft-clip signature: CIGAR **`S` consumes the read (SEQ) but not the reference** and sits only at a read end → leading `S` = **left / L-breakpoint**, trailing `S` = **right / R-breakpoint**, no-`S` read = no signal; the breakpoint is the **marginal point between clipped and matched sequence** (single-base resolution). Clustering rule: per-read junctions grouped under **positional tolerance** (ClipCrop "clustered within 5-base differences", 5 b), **minimum support** (SoftSearch ≥ x soft-clipped reads at a position, default 5 configurable to 2; support = # clipped reads in cluster) and **same clip-side + same chromosome**, with a **min clip length** filter (≤ 5 clipped bases dropped) and `RefineBreakpoint`. Repository `SplitRead` record stores the junction directly (`PrimaryPosition` = anchored SAM POS / `SupplementaryPosition` = junction coordinate / `ClipLength`). Oracles: 3 reads @ 5000 ± ≤ 5 b same chr support ≥ 2 → one breakpoint ~5000; 1 isolated read → none; 5000 vs 5100 (gap > tolerance) → two below-support groups → none. Sources: SAM/BAM Format Spec (samtools/hts-specs `SAMv1.tex`, rank 2 — CIGAR consume table, POS = leftmost reference-consuming op, SEQ = M/I/S/=/X so clips recoverable, S only at ends) + Tattini/D'Aurizio/Magi 2015 Front Bioeng Biotechnol review (rank 1 — anchored+imprecise split read, single-base SR resolution, cluster-of-split-reads) + ClipCrop Suzuki 2011 (rank 1 — junction=breakpoint, `31S69M`, L/R distinction, 5-base clustering) + SoftSearch Hart 2013 (rank 1 — ≥ x soft-clipped reads beginning at position y, orientation-consistent combining, min clip > 5). One ASSUMPTION (the sources fix the per-read junction + tolerance clustering but do NOT prescribe the cluster's reported coordinate statistic → this unit reports the rounded **mean** of member junctions, sub-tolerance only, mirroring the sibling `ClusterSplitReads`; does not change cluster membership or support count). Scope: operates on already-extracted `SplitRead` records, NOT raw BAM parsing; research-grade, not for clinical use. Wrote source [[sv-breakpoint-001-evidence]] + concept [[breakpoint-detection-split-reads]]. Cross-linked (prose [[mentions]], no unsupported typed edges): orthogonal to the oncology read-evidence rearrangement units [[gene-fusion-detection-read-evidence]] + [[chromothripsis-inference]], distinct from the gene-order [[genome-rearrangement-breakpoint-distance]] and chromosome-scale [[synteny-and-rearrangement-detection]]. Updated hub [[algorithm-validation-evidence]] (frontmatter sources +SV-BREAKPOINT-001 + source_commit→HEAD e0ce1587; body source-link enumeration + concept-enumeration entry) and wiki/index.md (Sources + Concepts entries). Graph: one typed edge on the new concept — relates_to test-unit-registry (standard per-unit registration; the source does not explicitly relate the method to another wiki concept's algorithm, so no further typed edges emitted per the source-supported-only rule). Contradictions: none — SAM spec + ClipCrop + SoftSearch + Tattini/Magi mutually consistent on the soft-clip-junction=breakpoint model, positional clustering within a small tolerance, and the min-support / min-clip-length filter. Follow-ups: likely future germline-SV siblings — discordant-pair clustering, SV candidate assembly + microhomology at the junction, CNV segmentation, SV genotyping/filtering/merging — would enrich this anchor.
   graph: +2 nodes, +1 typed edge (breakpoint-detection-split-reads relates_to test-unit-registry)
- 2026-07-10 — ingest docs/Evidence/SV-CNV-001-Evidence.md (Read-Depth Copy Number Variation Detection — windowed read depth → log2 ratio → integer copy number → deletion/duplication call). CONTEXT check per brief: read the SV anchor [[breakpoint-detection-split-reads]] (SV-BREAKPOINT-001) and searched wiki/concepts + index for cnv/copy-number/segmentation/depth/coverage — closest existing pages are the oncology [[copy-number-alteration-classification]] (ONCO-CNA-001, log2 → threshold-binned CNA states, which ALREADY names SV-CNV-001's DetectCNV/SegmentCopyNumber in prose), [[allele-specific-copy-number-ascat]] (BAF+purity), and the assembly [[coverage-depth-calculation]] (exact per-base depth). DECISION: read-depth CNV segmentation is a genuinely distinct method (aggregate depth signal, not split-read junctions) → NEW concept [[read-depth-cnv-segmentation]] as the read-depth member of the germline SV family, cross-linked to the split-read anchor. Pipeline: read depth ∝ CN (Yoon 2009) → windowed counting (non-overlapping fixed windows, one read per start, per-window mean RD) → optional GC correction r_i'=r_i·m/m_GC → log2 ratio log2(observed/reference depth), reference = overall median of window means (default self-reference) → integer CN = round(2·2^log2) (CNVkit "round", diploid ploidy 2) → del/dup call + segment merge (DetectCNV/SegmentCopyNumber). Anchors (verbatim CNVkit): log2(1/2)=−1→CN1 loss, 0→CN2 Neutral, log2(3/2)=+0.585→CN3 gain, +1→CN4 amp; CN clamped ≥ 0. Corner cases: zero-depth window → no-call (log2(0)=−∞ undefined, homozygous-deletion candidate, NOT −∞); NaN log2 → Neutral (CN 2). Two source-supported ASSUMPTIONs (reference baseline = overall median of non-zero window means; diploid ploidy 2). Sources: Yoon, Xuan, Makarov, Ye & Sebat 2009 (Genome Res 19(9):1586, rank 1, PMC2752127 — RD∝CN, windowed counting, GC formula) + CNVkit call.py + calling docs (rank 3, Talevich 2016 — n=ref_copies·2^log2, anchors, max(0,ncopies) clamp, nan→neutral). Wrote source [[sv-cnv-001-evidence]] + concept [[read-depth-cnv-segmentation]]. Cross-linked: enriched SV anchor [[breakpoint-detection-split-reads]] ("first germline-SV sibling") and oncology [[copy-number-alteration-classification]] (wikilinked the DetectCNV/SegmentCopyNumber prose to the new concept — bidirectional inbound link). Updated hub [[algorithm-validation-evidence]] (frontmatter sources +SV-CNV-001 + source_commit→HEAD 59811dac; body source-link + concept enumeration entries) and wiki/index.md (Sources + Concepts entries). Contradictions: none — Yoon 2009 (RD∝CN, windowing, GC) and CNVkit (log2→CN arithmetic, anchors, clamping) cover disjoint pipeline stages, mutually consistent. Follow-ups: remaining germline-SV siblings (discordant-pair clustering, breakpoint-junction assembly + microhomology, SV genotyping/filtering/merging) would enrich the anchor.
   graph: +2 nodes, +2 typed edges (read-depth-cnv-segmentation relates_to test-unit-registry; relates_to breakpoint-detection-split-reads)
- 2026-07-10 — ingest docs/Evidence/SV-DETECT-001-Evidence.md (Structural Variant Detection from Paired-End Mapping (PEM) signatures — discordant read pairs). CONTEXT check per brief: read the SV anchor [[breakpoint-detection-split-reads]] (SV-BREAKPOINT-001) and the read-depth sibling [[read-depth-cnv-segmentation]] (SV-CNV-001), plus the hub [[algorithm-validation-evidence]] and [[test-unit-registry]]. DECISION: discordant-read-pair (PEM) detection is a GENUINELY DISTINCT read-evidence signature from split reads (span + orientation geometry of a mate pair, not within-read clip junctions) → NEW dedicated concept [[discordant-pair-sv-detection]] as the third germline-SV read-evidence channel; cross-linked the split-read anchor (added the discordant-pair sibling to its relation paragraph). Method: concordant = FR proper pair (SAM FLAG 0x02, upstream + / downstream −; RF is everted/discordant, NOT proper); span-discordant iff insertSize < μ−c·σ OR > μ+c·σ (BreakDancer -c, default 3 s.d. of the empirical insert-size distribution); classify — span > μ+c·σ FR → Deletion, span < μ−c·σ → Insertion, FF/RR → Inversion, everted RF → tandem Duplication, mates on different chromosomes → Translocation; signature-then-cluster with min support ≥ -r read pairs (default 2). Datasets (mean=400/sd=50/c=3 ⇒ bounds [250,550]): FR span 5000 → Deletion, FR span 100 → Insertion, FF → Inversion, RF → Duplication, chr1≠chr2 → Translocation, FR span 400 → not discordant. One flagged ASSUMPTION — inter-chromosomal (CTX/translocation) precedence evaluated BEFORE orientation, since inversion (INV) is defined only for intra-chromosomal flipped pairs. Corner cases: an insertion larger than the insert size is invisible to PEM span and its sequence is not recovered (that is the split-read channel's job); below-support clusters dropped; cutoff bounds exact. Sources: Medvedev, Stanciu & Brudno 2009 (Nat Methods 6(11s):S13, rank 1 — signature catalogue: deletion larger span / insertion smaller span / inversion flipped orientation / linking = translocation) + BreakDancer README (Chen 2009, rank 3 — -c s.d. cutoff, -r min support 2, DEL/INS/INV/ITX/CTX codes) + BreakDancer protocol (Fan 2014, PMC3661775, rank 1 — six classes, 3–4 s.d. thresholds) + cureffi/BWA proper-pair (rank 4 — FR concordant, RF/FF/RR abnormal, FLAG 0x02) + DELLY 2012 / SVXplorer 2020 (rank 1 — RF ⇒ tandem duplication, FF/RR ⇒ inversion; "integrated paired-end and split-read analysis"). Wrote source [[sv-detect-001-evidence]] + concept [[discordant-pair-sv-detection]]; enriched SV anchor [[breakpoint-detection-split-reads]]; updated hub [[algorithm-validation-evidence]] (frontmatter sources +SV-DETECT-001 + source_commit→HEAD e525e311; body source-link + concept enumeration entries) and wiki/index.md (Sources + Concepts entries). Consumes already-mapped read-pair records, not raw BAM; research-grade, not for clinical use. Contradictions: none — Medvedev review, BreakDancer README/protocol, SAM/BWA proper-pair convention, and DELLY/SVXplorer mutually consistent on the discordant span/orientation signatures, the s.d. cutoff, FR-concordant/RF-everted, and signature-then-cluster with min support. Follow-ups: remaining germline-SV siblings — SV candidate assembly + microhomology at the junction, SV genotyping/filtering/merging — would further enrich the anchor.
   graph: +2 nodes, +2 typed edges (discordant-pair-sv-detection relates_to test-unit-registry; relates_to breakpoint-detection-split-reads)
- 2026-07-10 — ingest docs/Evidence/TRANS-CODON-001-Evidence.md (Codon Translation — codon → amino-acid via the NCBI genetic-code tables, `GeneticCode.Translate`, area Translation). CONTEXT check per brief: searched wiki/concepts + index for translation/codon/genetic-code/reading-frame/orf/transcription — the existing codon pages are all DOWNSTREAM of the table ([[relative-synonymous-codon-usage]], [[codon-adaptation-index]], [[effective-number-of-codons]], [[codon-optimization]], [[rare-codon-analysis]], [[codon-usage-comparison]]) and [[open-reading-frame-detection]] READS the standard table but does not model it; NO concept represented the genetic-code table itself. DECISION: **NEW concept warranted** — the foundational codon→AA lookup is genuinely distinct and is the table every codon-level operation reads. Created [[genetic-code-translation]] as the codon-family's base table. Content: `GeneticCode` (Seqeron.Genomics.Core) ships 4 of NCBI's 33 tables as static singletons — Standard(1)/VertebrateMitochondrial(2)/YeastMitochondrial(3)/BacterialPlastid(11), plus GetByTableNumber; codon→AA mappings + start/stop sets taken VERBATIM from the NCBI `AAs`/`Starts` strings (2024-09-23); `Translate`/`IsStartCodon`/`IsStopCodon`/`GetCodonsForAminoAcid`; DNA↔RNA T→U normalization, case-insensitive, stop→`'*'`, len≠3→ArgumentException; Met/AUG & Trp/UGG single-codon fixed points (the w≡1 fact behind CAI + optimization). Oracles: full 64-codon Standard table; table-2 diffs (AGA/AGG→`*`, AUA→M, UGA→W); table-3 diffs (CUN→T, AUA→M, UGA→W). Sources: NCBI Genetic Codes (transl_table 1/2/3/11, official spec) + Wikipedia Genetic/Start/Stop codon + Nirenberg-Matthaei 1961 / Crick 1968. **CONTRADICTION FLAGGED** (source-vs-implementation, API-contract layer not the code tables): the Evidence doc's Documented Corner Cases + Known Failure Modes promise `NNN`/"invalid nucleotide" → ArgumentException, but the actual `GeneticCode.Translate` returns `'X'` for any VALID IUPAC ambiguity codon (alphabet ACGURYMKSWBDHVN) and throws only for a non-IUPAC char (e.g. Z) — so NNN→'X', not an exception; recorded on both the concept and source pages for reconciliation. The mapping tables themselves match NCBI exactly (doc Deviations: None — accurate). Wrote source [[trans-codon-001-evidence]] + concept [[genetic-code-translation]]. Cross-linked (prose [[mentions]]): enriched [[open-reading-frame-detection]] (standard genetic code link) and [[codon-optimization]] (Met/Trp fixed-point link) — bidirectional inbound links to the new concept. Updated hub [[algorithm-validation-evidence]] (frontmatter sources +TRANS-CODON-001; body source-link enumeration; source_commit left per recent-ingest precedent) and wiki/index.md (Sources + Concepts entries). Graph: one typed edge on the new concept — relates_to test-unit-registry (per-unit registration); the ORF/optimization relationships are captured as auto-derived [[mentions]] edges rather than typed edges, since the TRANS-CODON source does not itself name those units (source-supported-only rule). Follow-ups: whole-sequence framed translation (`Translator`, six-frame) and the MCP `TranslateDna`/`TranslateRna` surface are adjacent units not yet ingested; the IUPAC-'X' vs promised-exception discrepancy is a candidate for a doc or test reconciliation.
   graph: +2 nodes, +1 typed edge (genetic-code-translation relates_to test-unit-registry)
- 2026-07-10 — ingest docs/Evidence/TRANS-DIFF-001-Evidence.md (Differential Expression — two-group RNA-seq DE: log2 fold change + Welch's unequal-variance t-test + Benjamini-Hochberg FDR; `TranscriptomeAnalyzer.CalculateFoldChange`/`FindDifferentiallyExpressed`, `Seqeron.Genomics.Annotation`). CONTEXT check per brief: searched wiki/concepts + index for differential/expression/rna-seq/transcriptome/fold-change/deseq/fdr/benjamini/multiple-testing/p-value — found NO transcriptome/DE concept and NO dedicated FDR/t-test statistics concept; the closest existing pages are the METAGENOMICS differential-abundance [[significant-taxa-detection]] (META-TAXA-001, non-parametric Mann–Whitney U, which explicitly notes it applies NO built-in BH/FDR — "the caller's responsibility"), the EPIGENETICS two-sample [[differentially-methylated-regions]] (Fisher's-exact + q-value gate), the ONCOLOGY single-sample [[expression-outlier-zscore-signature-score]] (z-score, not two-group), and the enrichment-downstream [[pathway-enrichment-ora]] — all distinct topics. DECISION: **NEW concept warranted** — RNA-seq two-group DE is a genuinely distinct method and the **first ingested unit of the Transcriptome/RNA-seq family** (siblings TPM/FPKM quantification, quantile normalization, PCA/clustering, alternative splicing not yet ingested). Created [[differential-expression]] as the family anchor. Content: **log2FC = log2((mean(treatment)+c)/(mean(control)+c))** pseudocount c=1, positive = up in treatment/condition 2 (DESeq2 Love 2014 + Science Park lesson); **Welch unequal-variance two-sample t-test** `t=(X̄₂−X̄₁)/√(s₁²/N₁+s₂²/N₂)`, unbiased (N−1) variances, Welch-Satterthwaite df, exact two-sided Student-t tail `p=I_{ν/(ν+t²)}(ν/2,½)` via regularized incomplete beta (Welch 1947 + Student-t CDF identity); **Benjamini-Hochberg FDR** across all genes reproducing R `p.adjust(method="BH")` `pmin(1,cummin(n/i·p[o]))[ro]` (BH 1995), monotone non-decreasing, adj p ∈ [raw p,1]; **two-criterion DE gate** `|log2FC| ≥ threshold (default 1.0) AND adjusted p < alpha (default 0.05)`, fail either ⇒ not significant. Oracles: UP `log2(41/11)=+1.8981204…` / DOWN −1.8981204… (exact negative) / FLAT 0; Welch {1,2,3} vs {7,8,9} → t=7.348469, ν=4, p=0.0018262607 (cross-checked vs SciPy `ttest_ind(equal_var=False)`); BH raw (0.001,0.4,0.5,0.9) → adj (0.004,0.6667,0.6667,0.9). Corner cases: <2 replicates → p=1 (variance undefined); se=0 → p=1 (equal means)/p=0 (unequal, t=±∞); zero mean → pseudocount. Three source-backed ASSUMPTIONs (pseudocount 1, <2-replicate p=1, se=0 conventions). Scope: a **simple two-group** estimator, NOT the full DESeq2 negative-binomial GLM (DESeq2 cited only for the log2FC definition/sign + BH-as-standard fact); research-grade, not for clinical use. Wrote source [[trans-diff-001-evidence]] + concept [[differential-expression]]. Cross-linked (prose [[mentions]] + one typed edge): `alternative_to` the non-parametric rank-sum [[significant-taxa-detection]] (parametric Welch-t on log2FC effect sizes + built-in BH vs rank-sum with no FDR); structurally akin to the Fisher's-exact two-sample [[differentially-methylated-regions]] two-criterion gate; distinct from the single-sample z-score [[expression-outlier-zscore-signature-score]]; DE gene lists feed [[pathway-enrichment-ora]]. Updated hub [[algorithm-validation-evidence]] (frontmatter sources +TRANS-DIFF-001 + source_commit→HEAD e00919fd; body source-link enumeration + concept-enumeration entry) and wiki/index.md (Sources + Concepts entries). Contradictions: none — DESeq2 (log2FC/BH), Welch 1947, the Student-t CDF regularized-incomplete-beta identity, and R `p.adjust` BH are mutually consistent. Follow-ups: the remaining Transcriptome/RNA-seq siblings (Expression_Quantification TPM/FPKM, quantile normalization/log2 transform, PCA/k-means clustering/co-expression, Alternative_Splicing PSI/deltaPSI — all present in `TranscriptomeAnalyzer`) would enrich this family anchor; a dedicated multiple-testing/BH or Welch-t statistics concept could be factored out if more DE-style units accrue.
   graph: +2 nodes, +2 typed edges (differential-expression relates_to test-unit-registry; alternative_to significant-taxa-detection)
- 2026-07-10 — ingest docs/Evidence/TRANS-EXPR-001-Evidence.md (Expression Quantification — RNA-seq TPM / FPKM / RPKM + quantile normalization; `TranscriptomeAnalyzer.CalculateTPM`/`CalculateFPKM`/`QuantileNormalize`, `Seqeron.Genomics.Annotation`; MCP `CalculateTpm`/`QuantileNormalize`). CONTEXT check per brief: read the family anchor [[differential-expression]] (TRANS-DIFF-001, just created) and searched wiki/concepts + index for tpm/fpkm/rpkm/expression/quantification/normalization/quantile — the only matches were oncology single-sample expression scoring [[expression-outlier-zscore-signature-score]] (z-score, not count-normalization) and the DE anchor itself; NO count-normalization concept exists. DECISION: **NEW concept warranted** — TPM/FPKM quantification + quantile normalization is a genuinely distinct method (within-/cross-sample count normalization) that sits UPSTREAM of the two-group DE test, a sibling not a duplicate. Created [[expression-quantification]] as the second Transcriptome/RNA-seq family unit. Content: **TPM** `(X_i/l_i)/Σ(X_j/l_j)·10⁶` — length-normalize into RPK THEN rescale so each sample sums to exactly 10⁶ (the sum-to-a-million invariant RPKM lacks; average TPM = 10⁶/#transcripts = const; Wagner/Kin/Lynch 2012, Zhao/Ye/Stanton 2020); **FPKM/RPKM** `X_i·10⁹/(l_i·N)` — per-kilobase + per-million-mapped-reads, RPKM (single-end) ≡ FPKM (paired-end) formula, `TPM=FPKM/Σ FPKM·10⁶` so FPKM does NOT sum to a constant (Pimentel 2014; Mortazavi 2008 for original RPKM); **quantile normalization** (Bolstad 2003) — sort each column, set each rank to the cross-column arithmetic rank mean, re-place at original positions, no reference sample; TIE rule: tied values get the mean of the rank means they would span. Oracles: TPM A(X=10,l=2000)/B(20,4000)/C(30,1000) → (125000,125000,750000)/Σ=10⁶; FPKM X=1000,l=2000,N=10⁶ → 500; quantile C1=(5,2,3,4)/C2=(4,1,4,2)/C3=(3,4,6,8), rank means 2.0/3.0/4.666…/5.666…, C2's two tied `4`s (rows A,C) → 5.166… (final 5.17). Corner cases: all-zero counts → Σ(X/l)=0 → TPM 0/0 undefined → emit 0 for all genes; non-positive length/N → FPKM 0, excluded from RPK; empty → empty (all three). Two source-backed ASSUMPTIONs (all-zero TPM → 0 degenerate convention; effective length = annotated length, `l̃_i=l_i` standard substitution). Scope: classic count-normalization layer, NOT an effective-length fragment model (kallisto/salmon) or TMM/median-of-ratios size factors (edgeR/DESeq2); TPM/RPKM are within-sample relative measures misused across samples/protocols (Zhao 2020); research-grade, not for clinical use. Wrote source [[trans-expr-001-evidence]] + concept [[expression-quantification]]. Cross-linked: reciprocal prose links with the anchor [[differential-expression]] (DE consumes normalized expression upstream; updated the anchor's closing line + downstream paragraph and its "future siblings" note) plus one typed `relates_to` edge to it; also `relates_to` [[test-unit-registry]]. Updated hub [[algorithm-validation-evidence]] (frontmatter sources +TRANS-EXPR-001 + source_commit→HEAD deb32560; body evidence-link enumeration) and wiki/index.md (Sources + Concepts entries). Contradictions: none — Wagner 2012, Zhao/Ye/Stanton 2020, Pimentel 2014 and Bolstad 2003 are mutually consistent. Follow-ups: remaining TranscriptomeAnalyzer siblings (PCA/k-means clustering/co-expression, alternative/differential splicing PSI/deltaPSI, RNA-seq library QC) not yet ingested.
   graph: +2 nodes, +2 typed edges (expression-quantification relates_to test-unit-registry; relates_to differential-expression)
- 2026-07-10 — ingest docs/Evidence/TRANS-PROT-001-Evidence.md (Whole-sequence protein translation — the `Translator` class: framed / six-frame translation + genetic-code-parameterized ORF finding, area Translation). CONTEXT check per brief: read [[genetic-code-translation]] (TRANS-CODON-001, the codon→AA table), [[open-reading-frame-detection]] (GENOMIC-ORF-001, ATG-only standard-code six-frame scanner), and searched wiki/concepts + index for translation/translator/six-frame/protein/reading-frame — the prior TRANS-CODON ingest had explicitly flagged "whole-sequence framed translation (`Translator`, six-frame)" as an adjacent not-yet-ingested unit that would ENRICH the codon-table concept. DECISION: **enrich, no new concept** — `Translator` is the sequence-level layer directly ABOVE the codon table it composes (via a `GeneticCode` parameter), so it belongs on [[genetic-code-translation]] rather than a standalone page; economical + strongly preferred by brief. Content added: framed translation (`frame` 0/1/2, else throws), six-frame dict keyed −3…+3 excluding 0 (3 fwd + 3 rev-comp, each 5'→3'), optional `toFirstStop`, DNA↔RNA T→U + case-insensitive + trailing partial codon untranslated, and `Translator.FindOrfs` (genetic-code-parameterized, min-length + both-strand) — deliberately NOT contract-equivalent to `GenomicAnalyzer.FindOpenReadingFrames` ([[open-reading-frame-detection]], ATG-only/standard-code) nor annotation-layer `GenomeAnnotator.FindOrfs`. Oracle: human insulin B chain (UniProt P01308 pos 25–54) 90-nt DNA → `FVNQHLCGSHLVEALYLVCGERGFFYTPKT` (30 aa); all 4 tables (1/2/3/11) verified codon-by-codon vs NCBI (2024-09-23). Sources: Wikipedia Translation-biology / Reading-frame / Open-reading-frame + NCBI Genetic Codes + UniProt P01308. **No contradictions** (Deviations: None). Wrote source [[trans-prot-001-evidence]]; enriched concept [[genetic-code-translation]] (frontmatter sources +TRANS-PROT-001 + source_commit bumped to HEAD; new "Whole-sequence framed translation" section; intro + Scope + Reference-sources updated; +1 typed graph edge relates_to test-unit-registry for the new unit). Updated hub [[algorithm-validation-evidence]] (frontmatter sources +TRANS-PROT-001; body source-link enumeration) and wiki/index.md (Sources entry + genetic-code-translation concept line). Follow-ups: the MCP `TranslateDna`/`TranslateRna` surface remains an adjacent not-yet-ingested wrapper unit.
   graph: +1 node, +1 typed edge (genetic-code-translation relates_to test-unit-registry for trans-prot-001-evidence)

## [2026-07-10] ingest | docs/Evidence/TRANS-SIXFRAME-001-Evidence.md → trans-sixframe-001-evidence (source)
   Same whole-sequence `Translator` six-frame surface as TRANS-PROT-001 — per the brief, ENRICHED
   [[genetic-code-translation]] rather than creating a redundant concept. Created source page
   trans-sixframe-001-evidence: reference-implementation angle (Biopython six_frame_translations =
   governing algorithm; EMBOSS transeq frame values 1/2/3/F/-1/-2/-3/R/6; EMBOSS getorf -find 1
   START→STOP; NCBI transl_table=1; Wikipedia Reading-frame). Distinctive new detail folded into the
   concept: the **reverse-frame numbering convention** — repo uses Biopython **independent-offset**
   (frame -k = reverse-complement translated 5'→3' at offset k-1, no correspondence to +1), the
   documented "alternative" to EMBOSS's phase-locked default; `ACTGG` frame -1 = `P` here vs `S` under
   EMBOSS = labelling convention, NOT a bug (flagged as such, no contradiction). Also clarified
   `Translator.FindOrfs` = START→STOP, start-residue included/stop excluded, minLength in **amino acids**
   (getorf -minsize is nucleotides). Oracles: 39-nt six-frame table; GGGATGAAACCCTAAGGG → MKP (start 3,
   end 14 incl.). Updated genetic-code-translation frontmatter (sources +TRANS-SIXFRAME-001,
   source_commit→HEAD 950ce49) + six-frame & ORF bullets + Reference-sources; wiki/index.md (Sources
   entry + concept line). Genetic-code-parameterized sibling of ATG-only [[open-reading-frame-detection]].
   graph: +1 node, +1 typed edge

## [2026-07-10] ingest | docs/Evidence/TRANS-SPLICE-001-Evidence.md → trans-splice-001-evidence (source) + alternative-splicing-psi (concept)
   TRANS-SPLICE-001 — RNA-seq alternative / differential splicing: event classification + Percent-
   Spliced-In (PSI, Ψ). Per brief, read the Transcriptome-family concepts ([[differential-expression]],
   [[expression-quantification]]) first and CONFIRMED this is RNA-seq read-quantification of splicing
   (PSI from inclusion/exclusion reads), NOT the genomic splice-site motif predictors
   [[splice-donor-site-prediction]] / [[splice-acceptor-site-prediction]] /
   [[gene-structure-prediction-intron-exon]]. DECISION: NEW concept [[alternative-splicing-psi]] —
   genuinely distinct from both the DE anchor (gene-level mean test) and the genomic splice family
   (sequence-motif scoring); anchored on the transcriptome family via [[differential-expression]], with
   a disambiguation cross-link (body + index) to the genomic splice family. `TranscriptomeAnalyzer.
   CalculatePSI` / `DetectAlternativeSplicing`: PSI Ψ = I/(I+S) (Wang 2008 / PMC3330053 / SUPPA2) +
   opt-in rMATS length-normalized ψ̂ = (I/lᵢ)/(I/lᵢ+S/lₛ) (Shen 2014); ΔPSI = splicing-level differential
   readout; five canonical AS classes (SE/RI/A5SS/A3SS/MXE, Wang 2008 = rMATS codes), event needs ≥2
   isoforms. Oracles: (80,20)→0.80; (80,20,200,100)→0.6666…. Corner cases: 0/0→NaN, S=0→1, I=0→0,
   0≤PSI≤1. Two source-backed ASSUMPTIONs (length normalization opt-in; forward strand). No source
   contradictions. Wrote source [[trans-splice-001-evidence]] + concept [[alternative-splicing-psi]];
   enriched anchor [[differential-expression]] (splicing-level counterpart cross-link); updated
   wiki/index.md (Sources + Concepts entries). Follow-ups: rMATS/SUPPA2 replicate ΔPSI significance +
   isoform-switching remain adjacent not-yet-ingested transcriptome units.
   graph: +2 nodes, +2 typed edges (alternative-splicing-psi relates_to test-unit-registry; relates_to differential-expression)

## [2026-07-10] ingest | VARIANT-ANNOT-001-Evidence.md → variant-annot-001-evidence (source) + variant-effect-annotation-vep (concept)
   VEP-style variant effect annotation: map an already-called variant to its functional consequence +
   Sequence-Ontology term/accession + IMPACT (HIGH/MODERATE/LOW/MODIFIER) via the Ensembl
   OverlapConsequence predicate system (Constants.pm consequence→IMPACT→rank table + VariationEffect.pm
   peptide predicates). IMPACT is stored on the term, not computed; most-severe = lowest rank
   (McLaren 2016). Coding engine translates ref/alt codons through the standard genetic code (NCBI table 1)
   and compares peptides: synonymous / missense / stop_gained / stop_lost; frameshift purely length-based
   |alt−ref| mod 3 ≠ 0; inframe ins/del = ×3 indel. Precedence stop_gained>missense, start_lost(rank 7)>
   coding-substitution. Oracles GAA→GTA missense · TTA→TTG synonymous · CAA→TAA stop_gained · TAA→CAA
   stop_lost · ATG→ATC start_lost · AC→A frameshift · A→ATTT inframe_insertion. Two assumptions (table 1
   only; single-codon SNV comparison). Wrote source [[variant-annot-001-evidence]] + NEW concept
   [[variant-effect-annotation-vep]]; enriched [[genetic-code-translation]] (added variant-annotation as a
   reader of the codon table); added to hub [[algorithm-validation-evidence]] (sources + body link + concept
   list); updated wiki/index.md (Sources + Concepts). No source contradictions; research-grade, not clinical.
   Follow-ups: variant calling (SNP/indel), pathogenicity/ACMG classification, and splice-region consequence
   scoring are adjacent variant-family units not yet ingested.
   graph: +2 nodes, +3 typed edges (variant-effect-annotation-vep relates_to test-unit-registry; depends_on genetic-code-translation; relates_to somatic-variant-calling-tumor-normal)

## [2026-07-10] ingest | VARIANT-CALL-001-Evidence.md → variant-call-001-evidence (source) + 1 concept
   Ingested germline variant calling (SNP/indel from a reference↔query global alignment + transition/
   transversion classification + Ti/Tv). Created source [[variant-call-001-evidence]] and NEW concept
   [[germline-variant-calling-snp-indel]] (the calling/detection member of the variant-analysis family:
   `CallVariantsFromAlignment` → `SequenceAligner.GlobalAlign`, SNP/Insertion/Deletion per aligned column,
   `"-"` gap sentinel + 0-based position in-memory, VCF padding/1-based POS only in serialized `ToVcfLines`;
   Ti/Tv with undefined `#Tv=0`→0). Cross-linked both siblings: enriched [[variant-effect-annotation-vep]]
   (named the germline caller as the upstream producer) and [[somatic-variant-calling-tumor-normal]] (added
   the germline reference↔query counterpart clause); both siblings' frontmatter gained the source path.
   Added to hub [[algorithm-validation-evidence]] (frontmatter source + roster body link); updated
   wiki/index.md (Sources + Concepts). Sources VCFv4.3 + Danecek 2011 + Tan 2015 + Collins & Jukes 1994;
   no source contradictions; research-grade, not clinical. Follow-ups: pathogenicity/ACMG classification and
   read-pileup (depth/genotype) calling remain adjacent variant-family units not yet ingested.
   graph: +2 nodes, +3 typed edges (germline-variant-calling-snp-indel relates_to test-unit-registry; relates_to variant-effect-annotation-vep; alternative_to somatic-variant-calling-tumor-normal)

## [2026-07-10] ingest | docs/Evidence/VARIANT-INDEL-001-Evidence.md → variant-indel-001-evidence (source)
   Indel detection (FindInsertions/FindDeletions — filters over the aligned-column caller). ENRICHED the
   existing concept [[germline-variant-calling-snp-indel]] rather than forking a new concept: it is the
   indel facet of the same VARIANT-CALL-001 caller. Added an "Indel detection" subsection (directional
   length invariant insertion⇒ALT>REF / deletion⇒REF>ALT, per-base multi-indel columns, minimal_representation
   CFTR/BRCA2 oracles) and expanded the normalization ASM with Tan 2015 Algorithm 1 (suffix-then-prefix
   trimming) + PharmCAT tandem-repeat left-shift. Concept frontmatter gained the source path + a
   relates_to test-unit-registry edge for the new unit. Updated wiki/index.md (Sources + Concepts).
   Sources VCFv4.3 + Tan 2015 (PMID 25701572) + minimal_representation (Minikel) + PharmCAT; no source
   contradictions; research-grade, not clinical. Follow-ups: indel left-normalization as a standalone
   method and read-pileup genotype calling remain adjacent units not yet ingested.
   graph: +1 node, +1 typed edge (germline-variant-calling-snp-indel relates_to test-unit-registry for variant-indel-001-evidence)

## [2026-07-10] ingest | docs/Evidence/VARIANT-SNP-001-Evidence.md → variant-snp-001-evidence (source)
   SNP detection (FindSnps alignment-based + FindSnpsDirect positional/Hamming-style). ENRICHED the
   existing concept [[germline-variant-calling-snp-indel]] rather than forking a new concept: it is the
   SNP facet of the same VARIANT-CALL-001 caller (mirrors the VARIANT-INDEL-001 enrichment). Added a
   "SNP detection" subsection (FindSnps = filter over the caller; FindSnpsDirect = Hamming-mismatch
   enumeration over equal-length sequences, SNP count = Hamming distance; equal-length/common-prefix
   precondition; REF==ALT is not a variant; case-insensitive; oracles ATGC→ATTC / AAAA→TGTA / VCFv4.3 §1.1
   G→A). Concept frontmatter gained the source path + a relates_to test-unit-registry edge for the new
   unit; lede + reference-sources paragraph updated to cite Acharya 2017 (PMC5410656, Hamming). Updated
   wiki/index.md (Sources + Concepts). Sources VCFv4.3 + Wikipedia/Futuyma Transversion + Acharya 2017 +
   Collins & Jukes 1994; no source contradictions; research-grade, not clinical. Follow-ups: read-pileup
   genotype calling and indel left-normalization remain adjacent units not yet ingested.
   graph: +1 node, +1 typed edge (germline-variant-calling-snp-indel relates_to test-unit-registry for variant-snp-001-evidence)

## [2026-07-10] ingest | docs/Validation/FINDINGS_REGISTER.md -> findings-register (source) + 1 concept
   GOVERNANCE. New source page for the validation-campaign disposition ledger (every note across
   86 per-unit reports triaged into FIXED-NOW/FEASIBLE/NOT-POSSIBLE/BY-DESIGN; green-washing detection;
   2026-06-12 snapshot SUPERSEDED by the 2026-06-24 re-validation reset). Created concept
   validation-findings-disposition (the A/B/C/D triage process + spec-not-impl green-washing remedy) with
   relates_to edges to validation-and-testing and build-quality-gate. ENRICHED validation-and-testing
   (Validation campaign paragraph now links findings-register + disposition; added the source path/commit).
   Did NOT force the algorithm-validation-evidence hub (this is a governance register, not per-algorithm
   Evidence). Updated wiki/index.md (Sources + Concepts). No contradictions flagged.
   graph: +2 nodes, +2 typed edges (validation-findings-disposition relates_to validation-and-testing, relates_to build-quality-gate)

## [2026-07-10] ingest | docs/Validation/LIMITATIONS.md -> limitations (source) + 1 concept
   GOVERNANCE. New source page for the validated operating-envelope document (what the library does
   NOT do; every row BY-DESIGN + ✅ CLEAN; three kinds irreducible/data-blocked/scope across ~13 units;
   research-vs-clinical disclaimer). Created concept operating-envelope-and-limitation-policy (the
   LimitationPolicy Strict/Moderate/Permissive runtime guard, minimum-access-mode table, and the
   limitation taxonomy) with relates_to edges to scientific-rigor and validation-and-testing. ENRICHED
   scientific-rigor (LimitationPolicy bullet now links the new concept + names the three modes),
   research-grade-limitations (LIMITATIONS.md sentence links the per-unit envelope), and
   validation-and-testing (operating-envelope-document phrase links [[limitations]] + the concept) —
   each with LIMITATIONS.md added to sources + source_commit bumped. Did NOT force the
   algorithm-validation-evidence hub (governance envelope, not per-algorithm Evidence). Updated
   wiki/index.md (Sources + Concepts). No contradictions flagged.
   graph: +2 nodes, +2 typed edges (operating-envelope-and-limitation-policy relates_to scientific-rigor, relates_to validation-and-testing)

## [2026-07-10] ingest | docs/Validation/VALIDATION_LEDGER.md -> validation-ledger (source), 3 pages enriched
   GOVERNANCE. New source page for the live per-unit validation status tracker — the ground-truth
   "where things stand" board (distinct from the test-unit-registry ID-scheme/spec and from the
   superseded findings-register snapshot). Captured the two-context protocol (fresh session/unit,
   implementer != validator, external primary sources, mutation checks), Stage A/B + State legend,
   and all three phases: Phase 1 (86 implemented CLEAN after the 2026-06-24 reset; 1 defect fixed
   PARSE-GENBANK-001; 19 re-reset 2026-06-25; +24 new campaign units) / Phase 2 (148 units, 13 genuine
   defects fixed) / Phase 3 (12 enhanced units, 1 latent defect PHYLO-NEWICK-001). Did NOT create a new
   concept (economical — enriched existing governance concepts instead) and did NOT force the
   algorithm-validation-evidence hub. ENRICHED test-unit-registry (added the live-status ledger vs
   registry-spec distinction; +ledger source/commit bump), validation-and-testing (campaign paragraph
   now links the ledger as ground truth; +ledger source/commit bump), validation-findings-disposition
   (ledger wikilinks in the live-status section), and the findings-register source page (back-link +
   Where-this-fits). Updated wiki/index.md (Sources). No contradictions flagged.
   graph: +1 node, +1 typed edge (validation-ledger supersedes findings-register)

## [2026-07-10] ingest | docs/Validation/VALIDATION_PROTOCOL.md → validation-protocol (source)
   Ingested the two-stage, one-session-per-unit validation METHODOLOGY doc that the ledger references.
   Created wiki/sources/validation-protocol.md: fresh-context-per-unit (implementer != validator), Stage A
   (validate description vs external primary sources — papers/textbooks/Biopython/EMBOSS/samtools/Rosalind,
   independent cross-check) BEFORE Stage B (validate code realises it); never fix code to a wrong spec; two
   completion end-states (✅ CLEAN / 🔧 LIMITED); verdict legend; report template; net10.0 green baseline
   4484/0 (2026-06-12); Phase-1 scope = 86 implemented units. Did NOT create a separate concept (economical —
   the source page is the methodology's canonical home) and did NOT force the algorithm-validation-evidence
   hub. ENRICHED validation-ledger (raw-path ref → [[validation-protocol]] wikilink + Where-this-fits bullet),
   validation-and-testing (campaign paragraph now names the protocol/Stage-A-before-B/end-states),
   validation-findings-disposition (green-washing remedy tied to the protocol's independent-external-source
   session), test-unit-registry ({UNIT-ID} is the handle a protocol session validates). Updated wiki/index.md
   (Sources). No contradictions flagged.
   graph: no typed edges (ontology has no relates_to/typed predicate targeting a source node, and the protocol
   neither supersedes nor contradicts anything; auto-derived sourced_from/mentions edges suffice).

## [2026-07-10] ingest | docs/Validation/reports/ALIGN-GLOBAL-001.md → align-global-001-report (source)
   First per-unit VALIDATION REPORT ingested (all prior sources were Evidence/governance docs). New source
   page for the Stage A/B validation write-up of ALIGN-GLOBAL-001 (Needleman–Wunsch): Stage A PASS (spec
   faithful to Wikipedia — border d·j/d·i, max-of-three recurrence, GCATGCG/GATTACA optimum 0, GapExtend-as-d),
   Stage B PASS (GlobalAlignCore+Traceback + cancellation overload, 13/13 cross-verification table, integer DP
   no overflow), State CLEAN, no defects; one documented non-defect (empty-DnaSequence overload returns
   score-0 empty alignment vs string overload's AlignmentResult.Empty). Kept it as a source summary distinct
   from the pre-impl [[align-global-001-evidence]]. ENRICHED concept [[global-alignment-needleman-wunsch]]
   (added report to sources + source_commit bump + a one-line CLEAN-verdict cross-link). Did NOT create a new
   concept (algorithm already represented) and did NOT force the algorithm-validation-evidence hub (that hub
   rosters Evidence artifacts, not reports). Updated wiki/index.md (Sources). No contradictions flagged.
   graph: no typed edges (source→source/concept report link is an auto-derived mentions/sourced_from edge; no
   new typed concept-to-concept predicate is warranted).

## [2026-07-10] ingest | docs/Validation/reports/ALIGN-LOCAL-001.md → align-local-001-report (source) + 1 concept
   New source page align-local-001-report (Stage A/B validation write-up, both PASS, State CLEAN, 7/7
   tests; hand-recomputed Wikipedia DP matrix, score 13). Created NEW concept local-alignment-smith-waterman
   — genuinely unrepresented (siblings global-alignment-needleman-wunsch and semi-global-alignment-fitting
   each had a concept; local existed only in their comparison tables; slug reserved in backlog). Tied the
   report to [[validation-ledger]]/[[validation-protocol]] (did NOT force the algorithm-validation-evidence
   hub — no Evidence artifact exists for this unit). Updated index.md (Sources + Concepts); moved the doc
   from backlog pending → covered. No contradictions flagged.
   graph: +2 nodes, +2 typed edges

## [2026-07-10] ingest | docs/Validation/reports/ALIGN-MULTI-001.md → align-multi-001-report (source)
   New source page align-multi-001-report (full re-validation write-up: Stage A/B both PASS, State CLEAN,
   96 MSA-family tests green, suite 18208/0; re-confirms star + progressive/UPGMA + consistency/T-Coffee,
   probe-verified GARFIELD relation 200→375, once-a-gap + zero-gap-DP + signal-add library). Kept distinct
   from the pre-impl [[align-multi-001-evidence]]. ENRICHED concept [[multiple-sequence-alignment]] (added
   report to sources + source_commit bump + a one-line re-validation CLEAN-verdict cross-link). Did NOT
   create a new concept (MSA/consensus/progressive/consistency already represented) and did NOT force the
   algorithm-validation-evidence hub (that hub rosters Evidence artifacts, not reports); tied the report to
   [[validation-ledger]]/[[validation-protocol]] from the source page. Updated wiki/index.md (Sources).
   One minor note (not a contradiction): the concept's numbered list of "three implementations" omits
   MultipleAlignProgressive (lists star/iterative/consistency), while the report treats progressive as a
   first-class named variant and iterative as the addendum sibling — left as a follow-up, not rewritten.
   graph: no typed edges (source→concept report link is an auto-derived mentions/sourced_from edge; no new
   typed concept-to-concept predicate is warranted).

## [2026-07-10] ingest | docs/Validation/reports/ALIGN-SEMI-001.md → align-semi-001-report (source) + 1 concept
   New source page align-semi-001-report (Stage A/B validation write-up: both PASS, State CLEAN, 17 canonical
   + property tests green; hand-recomputed fitting cases M1 4 / GAP 2 / MIX 1 / MAX 3, 11/11 cross-verification
   table, first-row-0 + no-zero-floor + argmax-last-row code evidence). Kept distinct from the pre-impl
   [[align-semi-001-evidence]] artifact. ENRICHED existing concept [[semi-global-alignment-fitting]] (added
   report to sources + source_commit bump + one-line re-validation CLEAN-verdict cross-link). Did NOT create a
   new concept (semi-global/fitting/overlap/glocal already represented) and did NOT force the
   algorithm-validation-evidence hub; tied the report to [[validation-ledger]]/[[validation-protocol]] from the
   source page. Updated wiki/index.md (Sources). No backlog move (validation reports are coverage-excluded, not
   a backlog slug). No contradictions flagged.
   graph: no typed edges (source→concept report link is an auto-derived mentions/sourced_from edge; no new
   typed concept-to-concept predicate is warranted).

## [2026-07-10] ingest | docs/Validation/reports/ALIGN-STATS-001.md → align-stats-001-report (source) + 1 concept
   New source page align-stats-001-report (Stage A/B two-stage validation write-up: both PASS-WITH-NOTES,
   State CLEAN, full suite 6536/0; validator independently re-derived the EMBOSS 65/90/9 → 43.6/60.4/6.0%
   numbers + the hand 9-col case, 7/7 cross-verification table, O(L) column-classify with `score.Mismatch>0`
   similarity rule as code evidence). Sole note on both stages is the rendering-only srspair `:`/`.` display
   simplification (no counted statistic affected). Kept distinct from the pre-impl [[align-stats-001-evidence]]
   artifact. ENRICHED existing concept [[alignment-statistics]] (added report to sources + source_commit bump +
   one-line PASS-WITH-NOTES/CLEAN verdict cross-link). Did NOT create a new concept (identity/similarity/gap/
   percent-identity/scoring-matrix already represented by [[alignment-statistics]]) and did NOT force the
   algorithm-validation-evidence hub; tied the report to [[validation-ledger]]/[[validation-protocol]] from the
   source page. Updated wiki/index.md (Sources). No backlog move (validation reports are coverage-excluded, not
   a backlog slug). No contradictions flagged.
   graph: no typed edges (source→concept report link is an auto-derived mentions/sourced_from edge; no new
   typed concept-to-concept predicate is warranted).

## [2026-07-10] ingest | docs/Validation/reports/ANNOT-CODING-001.md → annot-coding-001-report (source)
   New source page annot-coding-001-report (Stage A/B two-stage validation write-up: Stage A PASS-WITH-NOTES,
   Stage B PASS after in-session fix, State CLEAN, full suite 6561/0). Records a GENUINE defect found+fixed:
   the shipped CalculateCodingPotential had no both-zero branch and counted a both-in-both-tables-zero hexamer
   as a scored 0, contradicting canonical CPAT/lncScore FrameKmer.py which continues (not counted). Fix added
   the missing branch; the mirror-defect also lived in a code-echo test C1 (rewritten 0.6931->1.3863 sourced
   value) and three docs (algorithm doc 2.2 / Evidence / TestSpec). Discriminating cross-check: both-zero case
   scores ln4 = 1.3862943611198906, not the diluted value. Kept distinct from the pre-impl
   [[annot-coding-001-evidence]] artifact. ENRICHED existing concept [[coding-potential-hexamer-score]] (added
   report to sources + commit bump + one-line both-zero-defect cross-link on the exact branch). Did NOT create a
   new concept (coding-potential/CPAT/hexamer already represented) and did NOT force the
   algorithm-validation-evidence hub; tied the report to [[validation-ledger]]/[[validation-protocol]]. Updated
   wiki/index.md (Sources). No backlog move (validation reports are coverage-excluded, not a backlog slug). No
   contradictions flagged (the report itself corrects a prior doc/test transcription error, now consistent).
   graph: no typed edges (source->concept report link is an auto-derived mentions/sourced_from edge; no new
   typed concept-to-concept predicate is warranted).

## [2026-07-10] ingest | docs/Validation/reports/ANNOT-CODONUSAGE-001.md → annot-codonusage-001-report (source) + 1 concept
   New source page annot-codonusage-001-report (Stage A/B two-stage validation write-up for RSCU,
   GenomeAnnotator.GetCodonUsage at GenomeAnnotator.cs:922-992: Stage A PASS, Stage B PASS-WITH-NOTES, State
   CLEAN, full suite 6568/0). NO code defect — the formula n_i·x/Σx was confirmed VERBATIM against LIRMM/Rivals
   and cubar est_rscu, families/stops against NCBI table 1; the PASS-WITH-NOTES qualifier is solely three
   test-only coverage gaps closed in-session with zero code change (Trp single-codon=1.0, the 61-sense-codon
   set/stop exclusion, the empty-enumerable Array.Empty branch). Kept distinct from the pre-impl
   [[annot-codonusage-001-evidence]] artifact. ENRICHED existing concept [[relative-synonymous-codon-usage]]
   (added report to sources + commit bump 9ce49ba->987ea6c + one-line CLEAN-verdict cross-link). Did NOT create a
   new concept (RSCU already fully represented) and did NOT force the algorithm-validation-evidence hub; tied the
   report to [[validation-ledger]]/[[validation-protocol]]. Updated wiki/index.md (Sources). No backlog move
   (validation reports are coverage-excluded, not a backlog slug). No contradictions flagged.
   graph: no typed edges (source->concept report link is an auto-derived mentions/sourced_from edge; no new
   typed concept-to-concept predicate is warranted).

## [2026-07-10] ingest | docs/Validation/reports/ANNOT-GENE-001.md → annot-gene-001-report (source) + prokaryotic-gene-prediction-rbs (concept)
   New source page annot-gene-001-report (Stage A/B two-stage validation write-up for ORF-based prokaryotic
   gene prediction + Shine-Dalgarno RBS detection — GenomeAnnotator.PredictGenes / FindRibosomeBindingSites /
   FindRibosomeBindingSitesBothStrands / ScanStrandForShineDalgarno: Stage A PASS, Stage B PASS, State CLEAN,
   filtered suite GenomeAnnotator_Gene_Tests 39/0). NO code defect — an INDEPENDENT re-validation that
   re-fetched every consensus/spacing value (SD AGGAGG, anti-SD YACCUCCUUA, ~8 nt location, 5 nt optimal
   aligned spacing / Chen 1994) and HAND-re-derived every reverse-strand coordinate in Python without lifting
   repo expected values; the reverse-strand SD mapping forwardPos = len − hit.position − motifLen (scan the
   reverse complement, since the SD is an mRNA feature) was mutation-falsified (forwardPosition = hit.position
   → R1/R3/R4 fail), zero code change. UNLIKE the two prior report ingests, gene prediction + strand/spacing-
   aware SD RBS detection was genuinely UNREPRESENTED (open-reading-frame-detection is the Analysis-layer ORF
   sibling and explicitly excludes the annotation layer; regulatory-element-detection only lists a bare AGGAGG
   catalog string) — so CREATED a new concept [[prokaryotic-gene-prediction-rbs]] and added surgical back-links
   from [[open-reading-frame-detection]] (scope section) and [[regulatory-element-detection]] (SD row). Kept
   distinct from any annot-gene-001-evidence artifact; did NOT force the algorithm-validation-evidence hub;
   tied the report to [[validation-ledger]]/[[validation-protocol]]. Updated wiki/index.md (+1 source, +1
   concept). No backlog move — the Annotation/Gene_Prediction.md slug is resolved only when a concept lists that
   algorithm doc in sources:; this concept's source is the coverage-excluded validation report (per instruction
   #7), same as the ANNOT-CODING/CODONUSAGE report precedent. No contradictions flagged.
   graph: +2 nodes (source + concept), +2 typed edges (prokaryotic-gene-prediction-rbs relates_to
   test-unit-registry + alternative_to open-reading-frame-detection); body [[wikilinks]] mentions auto-derived.

## [2026-07-10] ingest | docs/Validation/reports/ANNOT-GFF-001.md → annot-gff-001-report (source) + bed-format-parsing (concept)
   New source page annot-gff-001-report (Stage A/B two-stage validation write-up for the ANNOTATION-LAYER GFF3
   I/O — GenomeAnnotator.ToGff3 / ParseGff3 / ComputeCdsPhases + Format/Encode/ParseGff3Attributes: Stage A
   PASS, Stage B PASS, State CLEAN; GFF3 fixture 46/46, full dotnet test Seqeron.sln 18783/0). NO code defect —
   a FRESH re-validation of the campaign export-fidelity fix (real source/score columns + per-transcript
   cumulative CDS phase on both strands), with SO GFF3 Spec v1.26 retrieved live and phases hand-recomputed
   against the canonical EDEN gene: plus-strand cds00003 = 0,1,1; cds00001 = 0,0,0,0; minus-strand (input
   order) = 2,2,0 — the load-bearing formula phase_i = (3 − Σ preceding lengths mod 3) mod 3, segments ordered
   5′→3′ (ascending start on +, DESCENDING on −). Kept DISTINCT from PARSE-GFF-001 / [[parse-gff-001-evidence]]
   (the FileIO GffParser, different code path + record type: GeneAnnotation 0-based half-open export vs
   GenomicFeature file-1-based parse). Did NOT create a new concept — the GFF3 9-column format is already
   anchored on [[bed-format-parsing]]; surgically updated that concept instead (added the report to sources: +
   a cross-link paragraph noting the second, annotation-layer GFF3 path and the CDS-phase formula). Did NOT
   force the algorithm-validation-evidence hub; tied the report to [[validation-ledger]] / [[validation-protocol]]
   / [[validation-and-testing]]. No annot-gff-001-evidence artifact exists. Updated wiki/index.md (+1 source).
   No backlog move — the Annotation/GFF3_IO.md (gff3-io) slug resolves only when a concept lists that ALGORITHM
   doc in sources:; this ingest covers the coverage-excluded validation report, same precedent as the prior
   ANNOT-* report ingests. No contradictions flagged.
   graph: +1 node (source annot-gff-001-report); no typed edges added (surgical sources:/cross-link only); body
   [[wikilinks]] mentions auto-derived.

## [2026-07-10] ingest | docs/Validation/reports/ANNOT-ORF-001.md → annot-orf-001-report (source) + open-reading-frame-detection (concept)
   New source page annot-orf-001-report (two-stage validation write-up for the ANNOTATION-LAYER ORF detection —
   GenomeAnnotator.FindOrfs / FindLongestOrfsPerFrame: Stage A PASS, Stage B PASS-WITH-NOTES, End state CLEAN;
   ORF filter 35/0, full project 18208/0). NO code defect — start set ATG+GTG/TTG, stops TAA/TAG/TGA, six
   frames, `minLength` in AMINO ACIDS, 0-based half-open stop-inclusive coords all confirmed vs Rosalind (live
   fetch) / Wikipedia / NCBI ORFfinder / Deonier 2005 / Claverie 1997; the authoritative Rosalind four-protein
   sample reproduced (incl. the reverse-strand MLLGSFRLIPKETLIQVAGSSPCNLS) + nested shared-stop ORF (start 24
   MGMTPRLGLESLLE nests start 30 MTPRLGLESLLE, shared end 69) verified in-code. Sole note = the non-canonical
   `requireStartCodon=false` run-off seeding path (outside standard-ORF scope). Kept DISTINCT from the ATG-only
   GenomicAnalyzer unit [[genomic-orf-001-evidence]] (different code path + start set + minLength unit) and from
   any annot-orf-001-evidence artifact (none exists). Did NOT create a new concept — [[open-reading-frame-detection]]
   already covers this; surgically updated it (added the report to sources: + a one-line verdict cross-link in the
   sibling-ORF-finder section). Did NOT force the algorithm-validation-evidence hub; tied the report to
   [[validation-ledger]] / [[validation-protocol]] / [[validation-and-testing]]. Updated wiki/index.md (+1 source).
   No contradictions flagged.
   graph: +1 node (source annot-orf-001-report); no typed edges added (surgical sources:/cross-link only); body
   [[wikilinks]] mentions auto-derived.

## [2026-07-10] ingest | docs/Validation/reports/ANNOT-PROM-001.md → annot-prom-001-report (source) + promoter-detection (NEW concept)
   New source page annot-prom-001-report (two-stage validation write-up for prokaryotic
   promoter motif detection — GenomeAnnotator.FindPromoterMotifs: Stage A PASS, Stage B PASS,
   End state CLEAN; promoter-motif filter 20/0 across 28 [TestCase]s, zero code change). No
   defect: −35 `TTGACA` / −10 `TATAAT` consensus, E. coli per-position probabilities (sums
   373 / 412), score = Σ matched-position p / Σ all-6 p ∈ [0,1] all confirmed vs Wikipedia
   Promoter-genetics / Pribnow-box (live fetch) + Pribnow 1975 / Harley & Reynolds 1987; all 8
   partial-variant scores hand-recomputed (full 6-mer + prefix-5/suffix-5/prefix-4, e.g.
   `TTGAC` 0.855 / `ATAAT` 0.813 / full boxes 1.000). Declared scope limits (17 bp spacing NOT
   enforced; exact-substring not PWM/HMM/mismatch-tolerant) are TestSpec-locked, not defects;
   `position` is a 0-based string index, not the biological TSS-relative negative coordinate.
   CREATED a new concept [[promoter-detection]] — genuinely unrepresented: this scored −10/−35
   detector is a DIFFERENT code path from the exact-hexamer catalog scan
   [[regulatory-element-detection]] (MOTIF-REGULATORY-001, GenomicAnalyzer.FindMotif), which
   the existing concept itself flagged as "un-ingested". Sourced the concept from both the
   report and docs/algorithms/Annotation/Promoter_Detection.md; surgically cross-linked
   regulatory-element-detection (replaced the "un-ingested" note with the sibling link). Did
   NOT force the algorithm-validation-evidence hub; tied the report to [[validation-ledger]] /
   [[validation-protocol]] / [[validation-and-testing]]. No annot-prom-001-evidence artifact
   exists. Backlog: moved `docs/algorithms/Annotation/Promoter_Detection.md` (promoter-detection)
   from pending Annotation (4→3) to covered-via-concept (75→76 covered, 170→169 pending), since
   the new concept lists that algorithm doc in sources:. Updated wiki/index.md (+1 source, +1
   concept). No contradictions flagged.
   graph: +2 nodes (source annot-prom-001-report, concept promoter-detection), +2 typed edges
   (promoter-detection --alternative_to--> regulatory-element-detection; --relates_to-->
   test-unit-registry); body [[wikilinks]] mentions auto-derived.

## [2026-07-10] ingest | docs/Validation/reports/ANNOT-REPEAT-001.md
   Created source-summary wiki/sources/annot-repeat-001-report.md (two-stage validation report
   for ANNOT-REPEAT-001 — repetitive element detection & classification; Stage A/B PASS, End
   state CLEAN, suite 6566/0, ledger row 58 / finding A17; one real defect found+fixed:
   ClassifyRepeat bidirectional-containment → one-directional element⊆query). Tied the report to
   [[validation-ledger]] / [[validation-protocol]] / [[validation-and-testing]] and the algorithm
   concept [[repetitive-element-detection]], NOT the algorithm-validation-evidence hub. Reused the
   existing concept (already covers ANNOT-REPEAT-001) — surgically added the report to its
   sources:, bumped source_commit, and cross-linked the defect fix in the deviation section; no
   new concept. Kept distinct from the pre-existing evidence page annot-repeat-001-evidence.
   Updated wiki/index.md (+1 source line). No contradictions flagged. No typed graph edges added
   (report is source-type; N/A).

## [2026-07-10] ingest | docs/Validation/reports/ASSEMBLY-CONSENSUS-001.md
   Created source-summary wiki/sources/assembly-consensus-001-report.md (two-stage validation
   report for ASSEMBLY-CONSENSUS-001 — consensus computation, SequenceAssembler.ComputeConsensus,
   the C of OLC; Stage A PASS-WITH-NOTES / Stage B PASS / End state CLEAN, full suite 6532/0, zero
   code or test change; validator re-ran Biopython 1.85 dumb_consensus reference, 10/10 datasets
   match). Tied the report to [[validation-ledger]] / [[validation-protocol]] and the algorithm
   concept [[consensus-sequence]], NOT the algorithm-validation-evidence hub. Reused the existing
   concept (already covers ASSEMBLY-CONSENSUS-001) — surgically added the report to its sources:,
   bumped source_commit, and noted the independent re-validation (10/10 Biopython match, two
   parameter-reachable default divergences carried, not defects). Kept distinct from the
   pre-existing evidence page assembly-consensus-001-evidence. Updated wiki/index.md (+1 source
   line). No contradictions flagged. No typed graph edges added (report is source-type; N/A).

## [2026-07-10] ingest | docs/Validation/reports/ASSEMBLY-CORRECT-001.md
   Created source-summary wiki/sources/assembly-correct-001-report.md (two-stage validation
   report for ASSEMBLY-CORRECT-001 — k-mer spectrum two-sided read error correction,
   SequenceAssembler.ErrorCorrectReads; Stage A PASS-WITH-NOTES / Stage B PASS / End state CLEAN,
   full suite 6535/0). Tied the report to [[validation-ledger]] / [[validation-protocol]] /
   [[validation-and-testing]] and the algorithm concept [[kmer-spectrum-error-correction]], NOT
   the algorithm-validation-evidence hub. Reused the existing concept (already covers
   ASSEMBLY-CORRECT-001) — surgically added the report to its sources:, bumped source_commit, and
   captured the two carried Stage-A notes (repo tests all covering k-mers vs Musket's
   leftmost+rightmost; `>=` vs `>` threshold) plus the fixed code-echoing test A22 (M4 rewritten
   to a genuine no-valid-correction case AAAAAAAA*3 + AACCAAAA, k=4; no algorithm defect). Kept
   distinct from the pre-existing evidence page assembly-correct-001-evidence. Updated
   wiki/index.md (+1 source line). No contradictions flagged. No typed graph edges added (report
   is source-type; N/A).

## [2026-07-10] ingest | docs/Validation/reports/ASSEMBLY-COVER-001.md
   Created source-summary wiki/sources/assembly-cover-001-report.md (two-stage validation report
   for ASSEMBLY-COVER-001 — coverage/depth, SequenceAssembler.CalculateCoverage→int[]; Stage A
   PASS-WITH-NOTES / Stage B PASS / End state CLEAN, full suite 6532→6533/0, zero code change).
   Tied to [[validation-ledger]] / [[validation-and-testing]] and the algorithm concept
   [[coverage-depth-calculation]] (assembly COVER anchor), NOT the algorithm-validation-evidence
   hub. Reused the existing concept — surgically added the report to its sources:, bumped
   source_commit, cross-linked the new report, and REFINED the boundary-clip rule: the report
   shows the min(pos+L,refLen) clip is dead/defensive (FindBestAlignment places reads only where
   they fit entirely, so the overhang partial-contribution case is unreachable — a near-contradiction
   with the concept's prior "contributes only its overlapping portion" wording, now reconciled).
   Fixed test gap noted (empty reference → []). Kept distinct from the pre-existing evidence page
   assembly-cover-001-evidence. Updated wiki/index.md (+1 source line). No backlog slug (none
   existed). No typed graph edges added (report is source-type; N/A).

## [2026-07-10] ingest | docs/Validation/reports/ASSEMBLY-DBG-001.md
   Created source-summary wiki/sources/assembly-dbg-001-report.md (two-stage validation report
   for ASSEMBLY-DBG-001 — de Bruijn graph assembly, SequenceAssembler.BuildDeBruijnGraph +
   AssembleDeBruijn / Hierholzer helpers; Stage A PASS / Stage B PASS-WITH-NOTES / End state
   CLEAN, full suite 6497/0, zero code change). Tied to [[validation-ledger]] /
   [[validation-and-testing]] and the algorithm concept [[de-bruijn-graph-assembly]] (assembly
   DBG anchor), NOT the algorithm-validation-evidence hub. Reused the existing concept —
   surgically added the report to its sources:, bumped source_commit to HEAD, and cross-linked
   the new report alongside the pre-existing evidence artifact. Only issue was a test-coverage
   gap (three untested-but-documented branches: disconnected-graph one-contig-per-component,
   MinContigLength filter, BuildDeBruijnGraph(null) guard), closed in-session with three
   exact-value tests; no algorithm defect, no contradictions. Kept distinct from the evidence
   page assembly-dbg-001-evidence. Updated wiki/index.md (+1 source line). No backlog slug
   (backlog tracks the algorithm doc, already covered; the report path is not a backlog entry).
   No typed graph edges added (report is source-type; N/A).

## [2026-07-10] ingest | docs/Validation/reports/ASSEMBLY-MERGE-001.md
   Created source-summary wiki/sources/assembly-merge-001-report.md (two-stage validation report
   for ASSEMBLY-MERGE-001 — contig merging, SequenceAssembler.MergeContigs(c1,c2,overlapLength);
   Stage A PASS / Stage B PASS / End state CLEAN, full suite 6529/0, zero code change). Tied to
   [[validation-ledger]] / [[validation-and-testing]] and the algorithm concept
   [[contig-merge-overlap-collapse]] (assembly MERGE anchor), NOT the algorithm-validation-evidence
   hub. Reused the existing concept — surgically added the report to its sources:, bumped
   source_commit to HEAD, and cross-linked the new report alongside the pre-existing evidence
   artifact assembly-merge-001-evidence. Merge = X+Y[l:], |merge|=|c1|+|c2|-l; single fallback
   covers l<=0 and l>min -> plain concat (suffixPrefixMatch guard); BAAB/AAABBBA trace to exact
   Langmead SCS/OLC printed strings, all 12 tests exact-value, HARD gate PASS. No defect, no
   contradictions. Updated wiki/index.md (+1 source line). No backlog slug (backlog tracks the
   algorithm doc, already covered; the report path is not a backlog entry). No typed graph edges
   added (report is source-type; N/A).

## [2026-07-10] ingest | docs/Validation/reports/ASSEMBLY-OLC-001.md
   Created source-summary wiki/sources/assembly-olc-001-report.md (two-stage validation report for
   ASSEMBLY-OLC-001 — Overlap-Layout-Consensus, SequenceAssembler.AssembleOLC + FindAllOverlaps
   (+cancellable) + FindOverlap; Stage A PASS / Stage B PASS-WITH-NOTES / End state CLEAN, full
   suite 6494/0, zero code change). Tied to [[validation-ledger]] / [[validation-and-testing]] and
   the algorithm concept [[overlap-layout-consensus-assembly]] (assembly OLC anchor), NOT the
   algorithm-validation-evidence hub. Reused the existing concept — surgically added the report to
   its sources:, bumped source_commit to HEAD, updated:, and cross-linked the new report alongside
   the pre-existing evidence artifact assembly-olc-001-evidence. Longest-suffix-prefix via
   descending scan (pos1=len1-L, pos2=0), case-insensitive identity fraction, no self-overlaps;
   greedy best-successor layout = superstring merge. Independent Python re-derived the 12-edge
   GTACGTACGAT set {4,5}, CTCTAGGCC len 6, chain->AAAAACCCCCGGGGGTTTTT, 7/8 identity gate. PASS-
   WITH-NOTES = two sourced intentional simplifications (greedy vs Hamiltonian-optimal; concat vs
   majority consensus). Only issue a test-coverage gap (MinContigLength discard), closed with
   exact-value test M5b. No defect, no contradictions. Updated wiki/index.md (+1 source line). No
   backlog slug (backlog tracks the algorithm doc, already covered; the report path is not a
   backlog entry). No typed graph edges added (report is source-type; N/A).

## [2026-07-10] ingest | docs/Validation/reports/ASSEMBLY-SCAFFOLD-001.md
   Created source-summary wiki/sources/assembly-scaffold-001-report.md (two-stage validation report
   for ASSEMBLY-SCAFFOLD-001 — scaffolding: SequenceAssembler.Scaffold(contigs, links,
   gapCharacter='N'); Stage A PASS / Stage B PASS / End state CLEAN, full suite 6529->6531/0, zero
   production-code change). Tied to [[validation-ledger]] / [[validation-and-testing]] and the
   algorithm concept [[scaffolding]] (assembly SCAFFOLD anchor), NOT the algorithm-validation-
   evidence hub. Reused the existing concept — surgically added the report path to its sources:,
   bumped source_commit to HEAD, updated:, and cross-linked the new report alongside the pre-
   existing evidence artifact assembly-scaffold-001-evidence. Concatenate-with-N-run construction
   (Jackman/ABySS 2.0 verbatim); gap rule gapLength = gapSize>0 ? gapSize : UnknownGapLength(100)
   (NCBI AGP unknown-size default, source-backed constant); `used` HashSet = one scaffold per
   contig, unreached -> length-1 ascending. Hand cross-check reproduced ACGTNNNTTGGNNCCAA (17),
   AAAA+100N+TTTT (108), gap-0 -> 100 N. HARD gate PASS (M3/M4 lock exactly 100 N, defeating a wrong
   Math.Max(1,gapSize)); two branch-coverage gaps closed with sourced tests (successor-already-
   placed -> AANNCCNNNGG; multi-forward-link first-declared tie-break -> ["AANNCC","GG"]). No defect,
   no contradictions; documented simplifications (link ranking / overlap resolution / orientation)
   out of scope. Updated wiki/index.md (+1 source line). No backlog slug (backlog tracks the
   algorithm doc, already covered; the report path is not a backlog entry). No typed graph edges
   added (report is source-type; N/A).

## [2026-07-10] ingest | docs/Validation/reports/ASSEMBLY-STATS-001.md
   Created source-summary wiki/sources/assembly-stats-001-report.md (two-stage validation report for
   ASSEMBLY-STATS-001 — assembly statistics: GenomeAssemblyAnalyzer.CalculateStatistics / CalculateNx
   (3-arg core + 2-arg delegate) / CalculateN50 / CalculateAuN / FindGaps + CalculateNxCurve wrapper;
   Stage A PASS / Stage B PASS / State CLEAN / test-quality PASS, full suite 6497/0, zero code or test
   change). Tied to [[validation-ledger]] / [[validation-and-testing]] and the algorithm concept
   [[assembly-statistics]] (assembly STATS anchor), NOT the algorithm-validation-evidence hub. Reused
   the existing concept — surgically added the report path to its sources:, bumped source_commit to
   HEAD, updated:, and cross-linked the new report alongside the pre-existing evidence artifact
   assembly-stats-001-evidence. Inclusive Nx boundary proven equal to QUAST s<=limit via
   cumulative*100 >= totalLength*threshold (long accum, no overflow); auN=Sigma l^2 / Sigma l
   (lh3/QUAST au_metric); FindGaps 0-based inclusive [Start,End] maximal N-runs + minGap filter. Hand
   cross-check reproduced Assembly A -> N50 70/L50 2, N90 30/L90 5, auN 57.586; Assembly B -> 50/3;
   {50,50} -> 50/1 (inclusive); minGap-5 gap filter. 23 tests, exact Is.EqualTo / Within(1e-10), no
   green-washing. No defect, no contradictions (Miller 2010 / Wikipedia / QUAST N50.py / Heng Li 2020
   agree). Non-defect notes: CalculateNxCurve wrapper + MedianLength upper-median (Assumption 2),
   empty->zeros vs QUAST None (Assumption 1). Updated wiki/index.md (+1 source line). No backlog slug
   (no assembly-stats entry in backlog). No typed graph edges added (report is source-type; N/A).

## [2026-07-10] ingest | docs/Validation/reports/ASSEMBLY-TRIM-001.md
   Created source page assembly-trim-001-report (validation report, distinct from the
   assembly-trim-001-evidence artifact). Verdict Stage A PASS / Stage B FAIL -> FIXED / State CLEAN;
   full suite 6535/0. One real defect fully fixed: original TrimEnd/TrimStart took the global-minimum
   partial sum with NO `s<0` early break and chained the 5' pass onto the 3'-surviving window (~19.5%
   of random reads diverged from cutadapt); fix adds the early break to both passes, runs each
   independently over the full read [0,n), and adds the start>=stop -> (0,0) drop rule. Python port vs
   cutadapt quality_trim_index = 0/900k mismatches. Tied to validation-ledger + validation-and-testing.
   Surgically corrected the concept quality-trimming-running-sum: its "5'-pass-on-surviving-window"
   wording and BWA-only early-break framing described the pre-fix (buggy) behaviour; updated to the
   validated both-passes-over-full-read + start>=stop-drop algorithm, added the report to sources,
   bumped source_commit. Updated wiki/index.md (+1 source line). No backlog slug (report path not a
   backlog row; Quality_Trimming algorithm doc already covered-via-concept). No new typed graph edges
   (report is source-type; N/A).

## [2026-07-10] ingest | docs/Validation/reports/CHROM-ALPHASAT-001.md → chrom-alphasat-001-report (source)
   Created wiki/sources/chrom-alphasat-001-report.md — validation report for CHROM-ALPHASAT-001
   (alpha-satellite monomer detection: ChromosomeAnalyzer.DetectAlphaSatellite / FindCenpBBoxes +
   171-bp / 17-bp-CENP-B-box constants). Stage A/B PASS, CLEAN, no code defect; one Stage-B test gap
   (non-ACGT excluded from AT denominator) closed. Tied to validation-ledger / validation-protocol /
   validation-and-testing / test-unit-registry; NOT forced onto algorithm-validation-evidence. Kept
   distinct from CHROM-CENT-001 (whole-centromere unit) — this is the narrow monomer-detection slice.
   Surgically updated concept centromere-analysis: added the report to sources, bumped source_commit,
   cross-linked the DetectAlphaSatellite bullet, +1 typed graph edge (chrom-alphasat-001-report
   relates_to test-unit-registry). Updated wiki/index.md (+1 source line). No backlog slug matched
   (alphasat/alpha-satellite absent from backlog).
   graph: +1 node, +1 typed edge

## 2026-07-10 — ingest docs/Validation/reports/CHROM-ANEU-001.md
   Created source page chrom-aneu-001-report (two-stage validation verdict for CHROM-ANEU-001,
   aneuploidy detection — `ChromosomeAnalyzer.DetectAneuploidy` / `IdentifyWholeChromosomeAneuploidy`,
   ChromosomeAnalyzer.cs:832–917): Stage A PASS / Stage B PASS / CLEAN, 31 passed 0 failed, zero code
   or test change. Distinct from the evidence artifact chrom-aneu-001-evidence (docs/Evidence). Tied to
   validation-ledger / validation-and-testing; concept anchor aneuploidy-detection.
   Surgically updated concept aneuploidy-detection: added the report to sources, bumped source_commit,
   cross-linked the report in the intro, +1 typed graph edge (chrom-aneu-001-report relates_to
   test-unit-registry). Updated wiki/index.md (+1 source line). No backlog slug matched (backlog row is
   the algorithm doc, already covered by the pre-existing concept).
   graph: +1 node, +1 typed edge

## 2026-07-10 — ingest docs/Validation/reports/CHROM-CENT-001.md
   Created source page chrom-cent-001-report (two-stage validation verdict for CHROM-CENT-001,
   centromere classification + α-satellite suprachromosomal-family assignment — new
   `ChromosomeAnalyzer.AssignSuprachromosomalFamily` / `LoadBundledAlphaSatelliteReference` + confirmation
   pass over Levan / DetectAlphaSatellite / DetectHigherOrderRepeat, ChromosomeAnalyzer.cs:1090/1187):
   re-validated 2026-06-26 after limitation-fix 887a9945 ADDED SF assignment (prior validation SUPERSEDED);
   Stage A PASS / Stage B PASS / CLEAN, 18860 passed 0 failed, zero code or test change. Bundled CC0 Dfam
   reference byte-verified (ALRb CENP-B box@126=B-type); SF rule = HOR period + A/B composition. Distinct
   from the evidence artifact chrom-cent-001-evidence (docs/Evidence) and from the narrow monomer-slice unit
   chrom-alphasat-001-report. Captured the LIMITED end-state: SF1-vs-SF2 unresolved + Sf1OrSf2Dimeric branch
   runtime-guarded to Permissive (LimitationPolicy) — CHROM-CENT-001 named in the operating-envelope doc.
   Tied to validation-ledger / validation-protocol / validation-and-testing / test-unit-registry; NOT forced
   onto algorithm-validation-evidence.
   Surgically updated concept centromere-analysis: added the report to sources, bumped source_commit to
   d0034a86, cross-linked the report + verdict in the intro, +1 typed graph edge (chrom-cent-001-report
   relates_to test-unit-registry). Updated wiki/index.md (+1 source line). No backlog slug matched (report
   is a generated per-run validation artifact, excluded from coverage).
   graph: +1 node, +1 typed edge

## [2026-07-10] ingest | docs/Validation/reports/CHROM-HOR-001.md → chrom-hor-001-report (source)
   Created wiki/sources/chrom-hor-001-report.md — validation report for CHROM-HOR-001 (higher-order
   repeat (HOR) detection: ChromosomeAnalyzer.DetectHigherOrderRepeat(sequence, monomerLength=171) →
   HorResult, ChromosomeAnalyzer.cs:751). Stage A/B PASS, CLEAN, no code defect; one Stage-B non-ACGT
   test gap closed (N tail dropped as partial monomer). HOR period = smallest k where ≥90% of k-spaced
   monomers ≥95% identical; unit k×171 bp, copy ⌊monomers/k⌋, HasHigherOrderStructure=period≥2 (k=1 1-mer
   is NOT a HOR); defining inter-HOR ≥ intra-HOR ordering confirmed vs McNulty&Sullivan 2018 / Rosandić
   2024 / Willard 1985 / Alkan 2007; independent k=4/m=7 hand cross-check (period 4, unit 684 bp, copy 7,
   inter 100% / intra 64.91%). NO new concept — HOR already synthesized in centromere-analysis; kept
   distinct from the monomer-slice unit chrom-alphasat-001-report (HOR out of scope there) and the
   whole-centromere unit chrom-cent-001-report. Tied to validation-ledger / validation-protocol /
   validation-and-testing / test-unit-registry; NOT forced onto algorithm-validation-evidence. Documented
   data-blocked boundary: suprachromosomal-family/HOR-family assignment not attempted (needs T2T-CHM13 HOR
   libraries). Surgically updated concept centromere-analysis: added the report to sources, bumped
   source_commit to 26fb94a8, cross-linked the DetectHigherOrderRepeat bullet, +1 typed graph edge
   (chrom-hor-001-report relates_to test-unit-registry). Updated wiki/index.md (+1 source line). No
   backlog slug matched (report is a generated per-run validation artifact, excluded from coverage).
   graph: +1 node, +1 typed edge

## [2026-07-10] ingest | docs/Validation/reports/CHROM-KARYO-001.md → chrom-karyo-001-report (source)
   Created wiki/sources/chrom-karyo-001-report.md — validation report for CHROM-KARYO-001 (karyotype
   analysis: ChromosomeAnalyzer.AnalyzeKaryotype(chromosomes, expectedPloidyLevel) + DetectPloidy(
   normalizedDepths, expectedDiploidDepth), ChromosomeAnalyzer.cs:136–241). Stage A/B both PASS, End
   state CLEAN, 36 exact-value tests 0 failed, zero code or test change. IMPORTANT scope clarification:
   the session prompt framed this around the Levan arm-ratio / centromeric-index classification
   (metacentric/submetacentric/subtelocentric/telocentric) — the report finds that is a DIFFERENT unit,
   CHROM-CENT-001 (CalculateArmRatio/ClassifyChromosomeByArmRatio/AnalyzeCentromere, centromere-analysis
   / chrom-cent-001-report) per ALGORITHMS_CHECKLIST_V2.md:73–74. CHROM-KARYO-001 is the
   karyotype/ploidy/aneuploidy unit. Absolute copy-count aneuploidy ladder (Nullisomy0…Pentasomy5) +
   ploidy/depth-ratio mapping confirmed vs Wikipedia Aneuploidy/Ploidy; 7-row DetectPloidy hand
   cross-check reproduced (ratio→round(ratio×2) clamped [1,8]; confidence uses the CLAMPED ploidy as
   reference → 0 conf out of range; true median). Boundary vs CHROM-ANEU-001: sibling
   IdentifyWholeChromosomeAneuploidy hardcodes a diploid baseline — correct for ANEU, not a KARYO
   divergence. Findings none. NO new concept — karyotype-analysis already synthesizes both algorithms;
   kept distinct from the evidence artifact chrom-karyo-001-evidence. Tied to validation-ledger /
   validation-protocol / validation-and-testing / test-unit-registry; NOT forced onto
   algorithm-validation-evidence. Surgically updated concept karyotype-analysis: added the report to
   sources, bumped source_commit to fcb5a4bc, cross-linked the report + Levan/CENT scope note in the
   intro, +1 typed graph edge (chrom-karyo-001-report relates_to test-unit-registry). Updated
   wiki/index.md (+1 source line). No backlog slug matched (report is a generated per-run validation
   artifact, excluded from coverage).
   graph: +1 node, +1 typed edge

2026-07-10 — Ingested docs/Validation/reports/CHROM-SYNT-001.md (validation report, synteny analysis
   — collinear blocks + rearrangement detection). Created wiki/sources/chrom-synt-001-report.md
   (Stage A PASS-WITH-NOTES / Stage B PASS / CLEAN, 19 synteny tests, zero code change; re-derived
   from fresh context, source unchanged since cb113ce). Enriched
   [[synteny-and-rearrangement-detection]]: added report to sources, bumped source_commit to
   7a7cdd29, cross-linked the report verdict alongside the evidence artifact in the intro. Updated
   wiki/index.md (+1 source line). Did not force the algorithm-validation-evidence hub; tied to
   validation-ledger/validation-protocol. No backlog slug matched (report path is coverage-excluded).
   graph: +0 nodes, +0 typed edges (report is source-type; concept already relates_to test-unit-registry)

## [2026-07-10] ingest | docs/Validation/reports/CHROM-TELO-001.md → chrom-telo-001-report (source)
   Created wiki/sources/chrom-telo-001-report.md — validation report for CHROM-TELO-001 (telomere
   analysis: ChromosomeAnalyzer.AnalyzeTelomeres(name, sequence, telomereRepeat="TTAGGG",
   searchLength=10000, minTelomereLength=500, criticalLength=3000) + EstimateTelomereLengthFromTSRatio(
   tsRatio, referenceRatio=1.0, referenceLength=7000) + constant HumanTelomereRepeat="TTAGGG",
   ChromosomeAnalyzer.cs:250–352). Validated 2026-06-24, Stage A PASS / Stage B PASS / End state ✅ CLEAN,
   build 0 warnings/errors, 33 Telomere tests pass, zero code change. TTAGGG (Moyzis 1988 / Meyne 1989,
   91 vertebrate species) + CCCTAA reverse-complement 5′ strand + Cawthon 2002 T/S proportionality
   confirmed vs first-sources; 3′=forward TTAGGG backward from terminus, 5′=CCCTAA forward from start,
   contiguous-tandem-till-<0.7-similarity counting (internal motif not counted), length=windows×repeatLen,
   purity=match/total, T/S=refLen×(tsRatio/refRatio). Independent hand cross-check reproduced 1200/purity-1.0,
   divergent TTAGGA→5/6, searchLen-600→600 truncation, T/S {1.5,0.5,2.0,1@ref2,0}→{10500,3500,14000,3500,0}.
   Findings none; in-spec note: 3′ scan phase-anchored to the terminus. NO new concept — telomere-analysis
   already synthesizes both algorithms; kept distinct from the evidence artifact chrom-telo-001-evidence
   (docs/Evidence). Tied to validation-ledger / validation-protocol / validation-and-testing /
   test-unit-registry; NOT forced onto algorithm-validation-evidence. Surgically updated concept
   telomere-analysis: added the report to sources, bumped source_commit to 9dfe8fee, cross-linked the
   report verdict alongside the evidence artifact in the intro, +1 typed graph edge (chrom-telo-001-report
   relates_to test-unit-registry). Updated wiki/index.md (+1 source line). No backlog slug matched (report
   is a generated per-run validation artifact, excluded from coverage; the algorithm doc is already
   covered by telomere-analysis).
   graph: +1 node, +1 typed edge

## [2026-07-10] ingest | docs/Validation/reports/CODON-CAI-001.md → codon-cai-001-report (source)
   Created wiki/sources/codon-cai-001-report.md — validation report for CODON-CAI-001 (Codon Adaptation
   Index — CAI, CodonOptimizer.CalculateCAI(codingSequence, table, excludeSingleCodonAminoAcids=false) +
   helper CalculateRelativeAdaptiveness + derived SingleCodonAminoAcids set, CodonOptimizer.cs:473–522,
   :131–144). Validated 2026-06-24 with 2026-06-25 re-validation, Stage A PASS / Stage B PASS / End state
   CLEAN; 34 CAI fixture tests + 18787 full Seqeron.Genomics.Tests pass, zero production-code change (4
   edge-case tests added for the 1e-6 zero-freq clamp + NaN no-data-AA skip). w=f/f_max geometric mean
   exp((1/L)·Σ ln w) confirmed vs Wikipedia + Sharp & Li 1987 (PMID 3547335) + Jansen 2003 (PMC2684136,
   verbatim single-codon-AA exclusion quote) + Kazusa; former D-A1 divergence resolved (opt-in
   excludeSingleCodonAminoAcids). Hand cross-checks reproduced to ≤1e-10 (CUAACU→0.17056, AGAAGG→0.07071,
   CUGCUA→0.28284; AUGUGG incl 1.0/excl 0.0, AUGCUACUA 0.18566/0.08; clamp CUACUG→0.001, CUA→1e-6,
   UUUCUG→1.0). NO new concept — codon-adaptation-index already synthesizes the algorithm; kept distinct
   from the evidence artifact codon-cai-001-evidence (docs/Evidence). Tied to validation-ledger /
   validation-and-testing / test-unit-registry; did NOT force the algorithm-validation-evidence hub.
   Surgically updated concept codon-adaptation-index: added the report to sources, bumped source_commit to
   01b6d4e5, cross-linked the report verdict alongside the evidence artifact in the intro, +1 typed graph
   edge (codon-cai-001-report relates_to test-unit-registry). Updated wiki/index.md (+1 source line). No
   backlog slug matched (report path is coverage-excluded per SCHEMA).
   graph: +1 node, +1 typed edge

## [2026-07-10] ingest | docs/Validation/reports/CODON-ENC-001.md → codon-enc-001-report (source)
   Created wiki/sources/codon-enc-001-report.md — validation report for CODON-ENC-001 (Effective
   Number of Codons — ENC/Nc, Wright 1990; CodonUsageAnalyzer.CalculateEnc(string) core +
   CalculateEnc(DnaSequence) delegate + private CalculateEncCore, CodonUsageAnalyzer.cs:274–360).
   Validated 2026-06-15, Stage A PASS-WITH-NOTES / Stage B PASS-WITH-NOTES / End state ✅ CLEAN; full
   dotnet test 6527 passed, 0 failed, zero production-code change. Formula (Eq. 1 F̂=(n·Σp²−1)/(n−1),
   Eq. 3 Nc=2+9/F̂₂+1/F̂₃+5/F̂₄+3/F̂₆, Eq. 4 within-class averaging, Eq. 5a Ile fallback, cap 61)
   confirmed verbatim vs Fuglsang 2004 (BBRC 317:957–964) + codonW (Peden thesis, "Nc not calculated"
   for empty class) + NCBI degeneracy partition (9 doublets/1 triplet/5 quartets/3 sextets/2 singlets).
   Independent Python reference reproduced to full double precision: M3=41.288461538461526,
   M5 Ile-absent Eq.5a=39.47394540942927, C1 2:1-bias=56.0, M1 one-codon-per-aa=20, M2 near-uniform
   cap=61. Fixed two code-echo tests (old M3=29.0/M5=40.4 asserted the unsourced full-count fallback);
   defect B1 (ClassContribution :357–360 returns raw codon count instead of declining, diverging from
   codonW; low severity, unreachable on real coding sequences) pinned by an explicitly LIBRARY-SPECIFIC
   M5b test (whole-class-absent→29.0). Note A2 = non-source lower clamp at 20. NO new concept —
   effective-number-of-codons already synthesizes the algorithm; kept distinct from the evidence
   artifact codon-enc-001-evidence (docs/Evidence). Tied to validation-ledger / validation-and-testing /
   test-unit-registry; did NOT force the algorithm-validation-evidence hub. Surgically updated concept
   effective-number-of-codons: added the report to sources, bumped source_commit to 816a85f7, updated
   date to 2026-07-10, cross-linked the report verdict + B1 divergence in the intro, +1 typed graph edge
   (codon-enc-001-report relates_to test-unit-registry). Updated wiki/index.md (+1 source line). No
   backlog slug matched (report path is coverage-excluded per SCHEMA).
   graph: +1 node, +1 typed edge

## [2026-07-10] ingest | docs/Validation/reports/CODON-OPT-001.md → codon-opt-001-report (source)
   New source-summary page for the two-stage validation report of CODON-OPT-001 (sequence/codon
   optimization, CodonOptimizer.OptimizeSequence). Stage A/B both PASS, State CLEAN, 59/0 tests, zero
   code change. Existing concept codon-optimization already synthesizes the algorithm; kept distinct
   from the evidence artifact codon-opt-001-evidence (docs/Evidence). Tied to validation-ledger /
   validation-and-testing / test-unit-registry; did NOT force the algorithm-validation-evidence hub.
   Surgically updated concept codon-optimization: added the report to sources, bumped source_commit to
   9dfa56d9, updated date to 2026-07-10, cross-linked the report PASS/CLEAN verdict in the intro. Updated
   wiki/index.md (+1 source line). No new typed graph edges (report is a source-summary; mentions
   auto-derived).
   graph: +1 node, +0 typed edges

## [2026-07-10] ingest | docs/Validation/reports/CODON-RARE-001.md → codon-rare-001-report (source)
   New source-summary page for the fresh 2026-06-25 two-stage re-validation of CODON-RARE-001 (rare
   codon detection — CodonOptimizer.FindRareCodons :663 + CalculateMinMaxProfile :720 +
   FindRareCodonClusters :825). Stage A PASS-WITH-NOTES / Stage B PASS / CLEAN, 45/45 tests, zero code
   change; %MinMax (Clarke&Clark 2008) and Sherlocc (Chartier 2012) formulas confirmed verbatim, all
   hand oracles reproduced live. Existing concept rare-codon-analysis already synthesizes the algorithm;
   kept distinct from the evidence artifact codon-rare-001-evidence (docs/Evidence). Tied to
   validation-and-testing (ledger) / test-unit-registry; did NOT force the algorithm-validation-evidence
   hub. Surgically updated concept rare-codon-analysis: added the report to sources, bumped source_commit
   to 8ce0af79, updated date to 2026-07-10, cross-linked the report verdict in the intro. Updated
   wiki/index.md (+1 source line). No backlog slug to move. No new typed graph edges (report is a
   source-summary; mentions auto-derived).
   graph: +1 node, +0 typed edges

## [2026-07-10] ingest | docs/Validation/reports/CODON-RSCU-001.md → codon-rscu-001-report (source)
   New source-summary page for the 2026-06-15 two-stage validation of CODON-RSCU-001 (Relative
   Synonymous Codon Usage + codon counting — CodonUsageAnalyzer.CalculateRscu :88 / CountCodons :37).
   Stage A PASS / Stage B PASS-WITH-NOTES / End state ✅ CLEAN, suite 6526/0, zero code change; formula
   n_i·x/Σx confirmed verbatim (LIRMM/GenomicSig/seqinr, Sharp Tuohy & Mosurski 1986), hand oracles
   reproduced live; PASS-WITH-NOTES = two documented test-coverage gaps closed (absent-family→0 guard,
   stop-codon 3-fold family), fixture 16→19, all test-only. Kept distinct from the same-measure sibling
   report annot-codonusage-001-report (GenomeAnnotator.GetCodonUsage) and from the evidence artifact
   codon-rscu-001-evidence (docs/Evidence). Existing concept relative-synonymous-codon-usage already
   synthesizes the measure; tied to validation-and-testing (ledger) / test-unit-registry; did NOT force
   the algorithm-validation-evidence hub. Surgically updated concept relative-synonymous-codon-usage:
   added the report to sources, bumped source_commit to e3c96b23, cross-linked both report verdicts in
   the intro. Updated wiki/index.md (+1 source line). No new typed graph edges (report is a
   source-summary; mentions auto-derived).
   graph: +1 node, +0 typed edges

## [2026-07-10] ingest | docs/Validation/reports/CODON-STATS-001.md → codon-stats-001-report (source)
   New source-summary page for the 2026-06-15 two-stage validation of CODON-STATS-001 (Codon Usage
   Statistics — the codon-family aggregation method CodonUsageAnalyzer.GetStatistics + CalculateCai,
   CodonUsageAnalyzer.cs:142/389). Stage A PASS-WITH-NOTES / Stage B PASS / End state ✅ CLEAN, suite
   6528/0, zero production-code change; CAI exp[(1/L)Σln w] + non-synonymous/stop exclusions, GC3s
   "synonymous 3rd-position GC excl. Met/Trp/stop" (Peden 1999 §1.8.2.1.3 verbatim), GC1/2/3 (EMBOSS
   cusp), RSCU n·x/Σx all confirmed against fetched sources (Wikipedia/seqinr/CodonW/Biopython/Kazusa);
   hand oracles reproduced live (√0.122=0.34928…, ∛(·)=0.011149…, S6=0.47706538…). PASS-WITH-NOTES = 3
   documented unit/edge choices (GC3s-as-percentage, skip-zero-w vs 0.01 floor, GC3s 6-fold subtlety);
   2 test-quality defects fixed in-session (bounds-only S6 strengthened to exact geometric mean; missing
   non-ACGT test added). Kept distinct from the evidence artifact codon-stats-001-evidence (docs/Evidence,
   the fuller aggregation description). Existing concept codon-adaptation-index synthesizes the co-canonical
   CalculateCai; tied to validation-and-testing (ledger) / test-unit-registry; did NOT force the
   algorithm-validation-evidence hub, and did NOT create a new codon-position-GC concept (positional GC /
   GC3s already documented on the evidence page — economical). Surgically updated concept
   codon-adaptation-index: added the report to sources, bumped source_commit to 518339cc, cross-linked the
   CODON-STATS-001 verdict in the intro. Updated wiki/index.md (+1 source line). No new typed graph edges
   (report is a source-summary; mentions auto-derived).
   graph: +1 node, +0 typed edges

## [2026-07-10] ingest | docs/Validation/reports/CODON-USAGE-001.md → codon-usage-001-report (source)
   New source-summary page for the 2026-06-24 two-stage validation of CODON-USAGE-001 (raw codon-usage
   table + TVD comparison — CodonOptimizer.CalculateCodonUsage :634 / CompareCodonUsage :657 /
   SplitIntoCodons :687, MolTools). Stage A PASS-WITH-NOTES / Stage B PASS / End state ✅ CLEAN, 22/22
   unit tests, full Seqeron.Genomics.Tests suite 18208/0, zero production-code change. Count(c)=raw
   Dictionary<codon,int> + Similarity=1−Σ|f₁−f₂|/2 (TVD=½·L¹, ∈[0,1]) confirmed vs EMBOSS cusp
   (Number/Frequency/Fraction), Kazusa row format, Wikipedia codon-usage-bias, TVD theory; hand oracles
   reproduced (ATGGCTGCTTAA→{AUG:1,GCU:2,UAA:1}; M9=0.75, M7=0.5, S6=2/3). Sole Stage-A note = scope
   framing (per-1000 frequency / per-family fraction / RSCU belong to CODON-RSCU-001 / CODON-STATS-001 /
   SEQ-CODON-FREQ-001), not a formula error; no defect, no code change. Kept distinct from the evidence
   artifact codon-usage-001-evidence (docs/Evidence). Existing concept codon-usage-comparison synthesizes
   the measure; tied to validation-ledger / validation-and-testing / test-unit-registry; did NOT force the
   algorithm-validation-evidence hub, and did NOT create a new concept (measure already represented —
   economical). Surgically updated concept codon-usage-comparison: added the report to sources, bumped
   source_commit to b0db43a8, cross-linked the CODON-USAGE-001 verdict in the intro. Updated wiki/index.md
   (+1 source line). No new typed graph edges (report is a source-summary; mentions auto-derived).
   graph: +1 node, +0 typed edges

## [2026-07-10] ingest | docs/Validation/reports/COMPGEN-ANI-001.md → compgen-ani-001-report (source)
   New source-summary page for the two-stage validation report of COMPGEN-ANI-001 (Average Nucleotide
   Identity — ANIb, ComparativeGenomics.CalculateANI / CalculateReciprocalAni). Independent re-validation
   of the 69c51fa0 limitations-campaign change (gapped fragment alignment via SequenceAligner.LocalAlign +
   reciprocal two-way ANI): Stage A/B both PASS, End state CLEAN, 20/20 ANI tests + 480/480 ~Comparative,
   zero code/test change; both prior PASS-WITH-NOTES resolved (minAlignableFraction now active; gapped +
   reciprocal replace ungapped/single-direction). Goris 2007 + pyani conventions confirmed verbatim; hand
   oracles G2 gapped 1.0>ungapped 0.875, R3 (1.0+1.0)/2=1.0. Kept distinct from the evidence artifact
   compgen-ani-001-evidence. Existing concept average-nucleotide-identity already represents the algorithm
   (economical — no new concept). Surgically updated that concept: added the report to sources, bumped
   source_commit to 205b259d, cross-linked the CLEAN verdict in the intro. Updated wiki/index.md (+1 source
   line). Tied to validation-ledger / validation-and-testing / test-unit-registry; did NOT force the
   algorithm-validation-evidence hub. No new typed graph edges (report is a source-summary; mentions
   auto-derived).
   graph: +1 node, +0 typed edges
## [2026-07-10] ingest | docs/Validation/reports/COMPGEN-CLUSTER-001.md → compgen-cluster-001-report (source)
   New source-summary page for the two-stage validation report of COMPGEN-CLUSTER-001 (Conserved Gene
   Clusters — common intervals of permutations, ComparativeGenomics.FindConservedClusters, cs:914–1021 +
   IsIntervalOf helper). Independent re-validation: Stage A PASS / Stage B PASS-WITH-NOTES / End state
   CLEAN — NO code defect and no code change; Stage B is with-notes solely because three weak test
   assertions were strengthened in-session (M3 Contains.Item→Is.EquivalentTo; M5 Does.Not.Contain→Is.Empty;
   S3 Does.Not.Contain/Contains→Is.EquivalentTo), all brute-forced, closing a green-wash gap. Scope
   clarified: the generic COG/OrthoMCL/MCL-clustering prompt notwithstanding, the sole method under test is
   the common-interval model (COG/OrthoMCL grouping → COMPGEN-ORTHO-001/RBH-001). Bui-Xuan/Habib/Paul 2013
   + Didier 2013 + Uno & Yagiura 2000 + Heber & Stoye 2001 confirmed verbatim; golden vector Id7 vs
   (7 2 1 3 6 4 5)→7 sets = paper Example 1, Didier {1,2,3,4} yes/{1,2} no, all brute-forced; full suite
   6605 passed/0 failed, build 0 errors. maxGap API-shape-only (strict gap-free; gene-teams not
   implemented). Kept distinct from evidence artifact compgen-cluster-001-evidence. Existing concept
   conserved-gene-clusters-common-intervals already represents the algorithm (economical — no new concept);
   surgically updated it: added the report to sources, bumped source_commit to 665dc336, cross-linked the
   CLEAN verdict in the intro. Updated wiki/index.md (+1 source line). Tied to validation-ledger /
   validation-and-testing / test-unit-registry; did NOT force the algorithm-validation-evidence hub. No new
   typed graph edges (report is a source-summary; mentions auto-derived).
   graph: +1 node, +0 typed edges
## [2026-07-10] ingest | docs/Validation/reports/COMPGEN-COMPARE-001.md → compgen-compare-001-report (source)
   New source-summary page for the two-stage validation report of COMPGEN-COMPARE-001 (comprehensive
   two-genome comparison — core/dispensable gene partition + overall syntenic-gene fraction,
   ComparativeGenomics.CompareGenomes, cs:765–810). CompareGenomes is an aggregator delegating to the
   already-validated sub-units COMPGEN-RBH-001 / SYNTENY-001 / REARR-001; it adds only three sourced
   pieces of unit logic (Tettelin 2005 pan-genome core/dispensable, Moreno-Hagelsieb 2008/Tatusov 1997
   RBH=shared gene, fraction-of-syntenic-genes metric + MCScanX Wang 2012). Independent re-validation:
   Stage A PASS / Stage B PASS / End state CLEAN — NO code defect and no code/spec/test change. Four
   invariants confirmed (Conserved=|orthologs|; core+specific_i=|genome_i|; OverallSynteny∈[0,1] via
   Math.Min(1.0,…); swap symmetry); 8 cases recomputed vs code (M1 1/1/1, M2 0/2/2, C1 2/0/0, M3
   5-collinear+1→5/1/1 & Synteny 5/6=0.8333 hand-traced, S1 3-collinear→Synteny 0, S2 symmetry, M4
   empty→all-0, Null→ArgumentException×2). Test-quality gate PASS: exact Is.EqualTo/Within(1e-10), prior
   permissive GreaterThan(-OrEqualTo) tests removed, full MUST/SHOULD/COULD + all-four-invariant coverage,
   honest green (full suite 6605 passed/0 failed, build 0 errors). One BY-DESIGN simplification inherited
   from COMPGEN-RBH-001 (alignment-free 5-mer Jaccard, id≥0.3/cov≥0.5, vs Tettelin 50%/50%); Stage-B
   nice-to-have notes only, not defects; no follow-ups. Existing concept genome-comparison-core-dispensable
   already represents the algorithm (economical — no new concept); surgically updated it: added the report
   to sources, bumped source_commit to 654fe336, cross-linked the CLEAN verdict in the intro. Updated
   wiki/index.md (+1 source line). Tied to validation-ledger / validation-and-testing / test-unit-registry;
   did NOT force the algorithm-validation-evidence hub. Kept distinct from evidence artifact
   compgen-compare-001-evidence. No new typed graph edges (report is a source-summary; mentions auto-derived).
   graph: +1 node, +0 typed edges

## [2026-07-10] ingest | docs/Validation/reports/COMPGEN-DOTPLOT-001.md → compgen-dotplot-001-report (source)
   New source-summary page for the two-stage validation report of COMPGEN-DOTPLOT-001 (dot plot —
   word-match / k-tuple dot matrix, ComparativeGenomics.GenerateDotPlot(seq1, seq2, wordSize=10,
   stepSize=1), cs:1169-1207). Validated 2026-06-16. Independent re-validation: Stage A PASS / Stage B
   PASS-WITH-NOTES / End state CLEAN — NO implementation defect and no product-code change. Match
   relation D={(i,j):A[i..i+w-1]=B[j..j+w-1]} case-insensitive (ToUpperInvariant), suffix-tree
   FindAllOccurrences yields ALL overlapping occurrences, eager wordSize/stepSize<=0→AOORE before
   iterator. Sources retrieved this session: Huttley TIB (k=1 rule, AGCGT/AT→{(0,0),(4,1)}, x=seq1/y=seq2),
   EMBOSS dottup (exact word match, default wordsize 10, noise/sensitivity trade-off), Wikipedia (Gibbs &
   McIntyre 1970, main diagonal). Ten cross-checks recomputed vs code all matched (incl. ACGTACGT self
   w=4→7-dot set, ACGT self w=1→exact main diagonal, disjoint→∅, default w=10→{(0,0)}). Two TEST-QUALITY
   fixes in-session (not code defects): M3 self-diagonal Is.SupersetOf→Is.EquivalentTo exact main diagonal
   (anti-green-wash), and new S3 locking the default wordSize=10 path. Honest green: full suite 6606
   passed/0 failed/0 skipped, build 0 errors. Axis-orientation / case-fold / non-positive-window-throws are
   documented decisions, not defects; no open follow-ups. Existing concept dot-plot-word-match already
   represents the algorithm (economical — no new concept); surgically updated it: added the report to
   sources, bumped source_commit to 37c54d6d, cross-linked the CLEAN verdict + the two in-session test
   fixes in the intro. Updated wiki/index.md (+1 source line). Tied to validation-ledger /
   validation-and-testing / test-unit-registry; did NOT force the algorithm-validation-evidence hub. Kept
   distinct from evidence artifact compgen-dotplot-001-evidence. No new typed graph edges (report is a
   source-summary; mentions auto-derived).
   graph: +1 node, +0 typed edges

## [2026-07-10] ingest | docs/Validation/reports/COMPGEN-ORTHO-001.md → compgen-ortho-001-report (source)
   New source-summary page for the two-stage validation report of COMPGEN-ORTHO-001 (ortholog detection
   by Reciprocal Best Hits + in-paralog identification, ComparativeGenomics.FindOrthologs / FindParalogs /
   FindReciprocalBestHits, cs:334–518). Validated 2026-06-15. Independent re-validation: Stage A
   PASS-WITH-NOTES / Stage B PASS-WITH-NOTES / End state CLEAN — NO code defect and no code change; the
   PASS-WITH-NOTES grades are honestly-documented alignment-free simplifications (A1 5-mer Jaccard replaces
   BLAST bit-score ranking; A2 FindParalogs within-genome mutual-best-hit proxy does not discriminate in-
   vs out-paralogs; A3 the 50% coverage gate is largely subsumed by the identity gate, a 200k-pair brute
   force found no separating input — a consequence of A1). FindOrthologs delegates to FindReciprocalBestHits
   (single source of truth, locked by S5); the historical non-reciprocity defect is already fixed (M2:
   a1↔b1 kept, b2 excluded). Fitch 1970 + Moreno-Hagelsieb 2008 (RBH def, ≥50% coverage, E≤1e-6) + Tatusov
   1997 + Remm 2001 confirmed; Python 5-mer Jaccard recomputation reproduced every asserted identity (1.0 /
   0.667 ranking / 0.0 rejection) and showed TtBlock vs GcBlock = 0.5 not 0.0; hand oracles M1–M6/S1
   reproduced exactly; full suite 6506 passed/0 failed, build 0 errors. Three in-session test-quality fixes
   (corrected the wrong "Jaccard 0.0" comment; added OrthologPair.Coverage assertions; direct S5/S6 tests
   for FindReciprocalBestHits). Existing concept ortholog-detection-reciprocal-best-hits already represents
   the algorithm (economical — no new concept); surgically updated it: added the report to sources, bumped
   source_commit to 0752d91e, cross-linked the CLEAN verdict + in-session fixes in the intro. Updated
   wiki/index.md (+1 source line). Tied to validation-ledger / validation-and-testing / test-unit-registry;
   did NOT force the algorithm-validation-evidence hub. Kept distinct from evidence artifact
   compgen-ortho-001-evidence. No new typed graph edges (report is a source-summary; mentions auto-derived).
   graph: +1 node, +0 typed edges

## [2026-07-10] ingest | docs/Validation/reports/COMPGEN-RBH-001.md → compgen-rbh-001-report (source)
   New source-summary page for the two-stage validation report of COMPGEN-RBH-001 (reciprocal best hits
   — the core between-genome ortholog primitive ComparativeGenomics.FindReciprocalBestHits, delegate
   FindOrthologs, private FindBestHit / CalculateSequenceSimilarity, cs:410–549; the RBH-only slice of
   COMPGEN-ORTHO-001, no in-paralog rule). Validated 2026-06-16. Independent re-validation: Stage A
   PASS-WITH-NOTES / Stage B PASS / End state CLEAN — NO code defect and no code change. The one Stage-A
   note is the documented alignment-free simplification (5-mer Jaccard replaces BLAST bit-score ranking,
   order-preserving on all datasets); Stage B is a clean PASS. Moreno-Hagelsieb 2008 (RBH def verbatim
   from the PubMed abstract; ≥50%/E≤1e-6 body quotes paywalled but the load-bearing definition
   independently confirmed) + Tatusov 1997 (COG mutually-consistent BeTs) + Best-Match-Graph literature
   confirm the symmetric requirement. Formula RBH(a,b) ⇔ bestHit(a→G2)=b ∧ bestHit(b→G1)=a with
   deterministic tie-break + both gates; independent Python 5-mer Jaccard reproduced self-match 1.0 /
   alignLen 14, superstring 0.667, and the new coverage case AAAAACCCCCGGGGG vs AAAAACCCCCTTTTT → shared
   6 / union 16 = 0.375, cov 6/11 = 0.5455. 13 unit tests; M2 excludes the non-reciprocal b2. Two
   test-coverage gaps closed in-session (test-surface only, impl already correct): M7 (coverage-gate
   rejection, exact 6/16 & 6/11 values; kept at default minCoverage 0.5, rejected at 0.6) and S4 (< k=5
   short-sequence similarity-0). Full suite 6605 passed / 0 failed / 1 skipped (unrelated MFE benchmark),
   build 0 warnings / 0 errors. Existing concept ortholog-detection-reciprocal-best-hits already
   represents the algorithm (RBH IS its core method — economical, no new concept); surgically updated it:
   added the report to sources, bumped source_commit to 00c5ea42, cross-linked the RBH-001 CLEAN verdict
   + the two in-session coverage-gap closures (M7/S4) in the intro. Updated wiki/index.md (+1 source
   line). Tied to validation-ledger / validation-and-testing / test-unit-registry; did NOT force the
   algorithm-validation-evidence hub. Kept distinct from evidence artifact compgen-rbh-001-evidence and
   sibling report compgen-ortho-001-report. No new typed graph edges (report is a source-summary;
   mentions auto-derived).
   graph: +1 node, +0 typed edges

## [2026-07-10] ingest | docs/Validation/reports/COMPGEN-REARR-001.md → compgen-rearr-001-report (source)
   New source-summary page for the two-stage validation report of COMPGEN-REARR-001 (genome
   rearrangement detection by breakpoints on a signed gene-order permutation,
   ComparativeGenomics.DetectRearrangements / ClassifyRearrangement, cs:581-723). Validated 2026-06-15.
   Independent re-validation: Stage A PASS-WITH-NOTES / Stage B PASS-WITH-NOTES / End state CLEAN — NO
   correctness defect in DetectRearrangements. The doc's "breakpoint iff y≠x+1" reduction PROVEN exact
   vs the full Hunter criterion (against beta=identity both clauses (x,y) and (-y,-x) collapse to
   y=x+1); Hunter alpha=(-2,-3,+1,+6,-5,-4)→b=6 reproduced incl. the (-5,-4) exclusion, plus independent
   hand cross-check (M2=0, M3=2, M4=2, S2=2, C1∈[0,n+1]) computed before running code. Iterator builds
   signedRank=sign×rank, relabels markers to 1..n, walks [0,…,n+1] emitting a breakpoint on curr≠prev+1;
   eager null checks split from the yield-iterator, <2 markers→yield break, dangling ortholog
   TryGetValue-guarded; ClassifyRearrangement re-parses TargetPosition and delegates to the emission-time
   ClassifyBoundary (stored Type == re-classification). Two test-coverage gaps found & fixed in-session
   (the only Stage-B defects): M9b (null genome2 throws ArgumentNullException) and M10
   (ClassifyRearrangement fallback returns stored Type on null/unparsable TargetPosition) → 16 REARR
   tests. Full unfiltered suite 6508 passed / 0 failed, build 0 errors, warning-free. Two documented
   notes (not defects): per-boundary Inversion/Transposition classifier is an intentionally simplified
   heuristic with no formal single-permutation basis (doc §5.3 / Evidence Assumption 3); d_BP=n−sim
   equality is telomere-convention-dependent while the unit reports/tests the extended-permutation count
   b(alpha). Sources Hunter Lecture 16 / Tannier PMC3887456 / Bafna–Pevzner 1998. Existing concept
   genome-rearrangement-breakpoint-distance already represents the algorithm (economical — no new
   concept); surgically updated it: added the report to sources, bumped source_commit to 4c3caf90,
   cross-linked the CLEAN verdict + the two coverage-gap fixes in the intro. Updated wiki/index.md (+1
   source line). Tied to validation-ledger / validation-and-testing / test-unit-registry; did NOT force
   the algorithm-validation-evidence hub. Kept distinct from evidence artifact compgen-rearr-001-evidence.
   No new typed graph edges (report is a source-summary; mentions auto-derived).
   graph: +1 node, +0 typed edges

## [2026-07-10] ingest | docs/Validation/reports/COMPGEN-REVERSAL-001.md → compgen-reversal-001-report (source)
   Two-stage validation report for COMPGEN-REVERSAL-001 (reversal/inversion distance = unsigned
   breakpoint lower bound ceil(b/2), CalculateReversalDistance, ComparativeGenomics.cs:840-880).
   Stage A PASS / Stage B PASS / End state CLEAN — no defect, no code or test change. Unsigned
   |Delta|!=1 = Hubotter 2020 Def 2.1; d>=b/2 (Corollary 2.1.1) => ceil(b/2); lower bound not exact HP
   distance (by design). Existing concept genome-rearrangement-breakpoint-distance already documents the
   algorithm (economical — no new concept); surgically enriched it: added the report to sources, bumped
   source_commit to e4a1444b, cross-linked the CLEAN verdict in the intro + the Unsigned-reversal
   section. Updated wiki/index.md (+1 source line). Tied to validation-ledger / validation-and-testing /
   test-unit-registry; did NOT force the algorithm-validation-evidence hub. Kept distinct from evidence
   artifact compgen-reversal-001-evidence and sibling compgen-rearr-001-report.
   No new typed graph edges (report is a source-summary; mentions auto-derived).
   graph: +1 node, +0 typed edges

## [2026-07-10] ingest | docs/Validation/reports/COMPGEN-SYNTENY-001.md → compgen-synteny-001-report (source)
   Two-stage validation report for COMPGEN-SYNTENY-001 (whole-genome syntenic-block detection, MCScanX
   collinearity DP scoring; ComparativeGenomics.FindSyntenicBlocks(g1, g2, orthologMap, minAnchors=5,
   maxGap=25) + VisualizeSynteny; ComparativeGenomics.cs:84-299). Stage A PASS-WITH-NOTES / Stage B
   PASS / End state CLEAN — no defect, no code change; tests strengthened +3. DP recurrence + constants
   (MatchScore 50, GapPenalty -1, NumberofGaps<25, report >250 i.e. >=5 pairs) confirmed verbatim vs
   Wang 2012 (Oxford Academic HTML, WebFetch) + wyp1125/MCScanX README; impl score=n*50-Sigma(|dpos2|-1)
   proven the single-monotone-chain closed form. Two Stage-A doc-only notes (MAX_GAPS paper=25 vs current
   tool default 20 = F-SYNTENY-001; "over 250" resolved to >=250 per the paper's own equivalence). +3
   tests close real gaps (S2b null genome2, S2c null orthologMap = full 3-arg contract; S4
   direction-switch flush branch/INV-3). Honest green 6504 passed/0 failed, warning-free. One BY-DESIGN
   note: greedy single-pass chaining == predecessor-DP for direction-consistent/gap-bounded/non-interleaved
   inputs (doc 5.2/5.3). Existing concept synteny-and-rearrangement-detection is the shared synteny anchor
   and already documents the MCScanX DP model (economical — no new concept); surgically enriched it: added
   the report to sources, bumped source_commit to 3d86b2b7, cross-linked the CLEAN verdict + the +3 tests
   in the MCScanX DP model section. Updated wiki/index.md (+1 source line). Tied to validation-ledger /
   validation-and-testing / test-unit-registry; did NOT force the algorithm-validation-evidence hub. Kept
   distinct from evidence artifact compgen-synteny-001-evidence and the chromosome-scale chrom-synt-001-report.
   No new typed graph edges (report is a source-summary; mentions auto-derived).
   graph: +1 node, +0 typed edges

## [2026-07-10] ingest | docs/Validation/reports/CRISPR-GUIDE-001.md → crispr-guide-001-report (source) + 1 concept
   Per-unit validation report for CRISPR-GUIDE-001 (MolTools — CRISPR gRNA design + on-target efficacy
   scoring). Stage A/B both PASS, State CLEAN, 54/54 CRISPR tests green, zero code/test change. CRISPR was
   genuinely unrepresented in the wiki (no prior crispr/gRNA/PAM/Doench concept), so created the new concept
   crispr-guide-rna-design as the CrisprDesigner/MolTools anchor — two layers: (1) composition heuristic
   DesignGuideRnas/EvaluateGuideRna, (2) learned on-target efficacy Doench-2014 Rule Set 1
   (CalculateOnTargetDoench2014, 30-mer 4+20+3+3, intercept 0.59763615, 70-entry table byte-identical to the
   re-downloaded CRISPOR reference) + Doench-2016 Rule Set 2 / Azimuth (CalculateOnTargetRuleSet2, sklearn-free
   GBRT, externally-derived oracles). Source-summary page records the re-grounding, three worked oracles
   (≤2e-8), NGG PAM guard, and by-design heuristic simplifications. Tied to validation-ledger /
   validation-and-testing / test-unit-registry; did NOT force the algorithm-validation-evidence hub. Covered
   backlog slug guide-rna-design (moved Guide_RNA_Design.md pending→covered; MolTools 17→16; status 76→77
   covered / 169→168 pending). Updated wiki/index.md (+1 source, +1 concept). Concept sources list both the
   report and docs/algorithms/MolTools/Guide_RNA_Design.md. No contradictions.
   graph: +2 nodes, +2 typed edges (crispr-guide-rna-design relates_to test-unit-registry, relates_to primer-dimer-thermodynamics-tm)

## [2026-07-10] ingest | docs/Validation/reports/CRISPR-OFF-001.md → crispr-off-001-report (source)
   Per-unit validation report for CRISPR-OFF-001 (MolTools — CRISPR off-target scoring: scan a genome for
   near-matches to a gRNA and score off-target risk). Two independent parts, both Stage A/B PASS, State
   CLEAN, zero production-code change: (1) MIT/Hsu-Zhang 2013 single-hit + aggregate specificity
   (CalculateMitHitScore / CalculateMitSpecificityScore; 20-element W byte-identical to CRISPOR hitScoreM;
   hitScore=Pi(1-W[i])*score2*score3*100, aggregate 100/(100+Sigma)*100; orientation index 0=PAM-distal,
   19=PAM-proximal seed; oracles perfect->100, mm@5->60.5, mm@13->14.9, mm@19->41.7, agg{60.5}->62.30529595
   reproduced from source; committed orientation-guard test), and (2) CFD Doench-2016 (CalculateCfdScore
   in [0,1]; 240-entry mismatch + 16-entry PAM matrices cross-checked 240/240 + 16/16 identical across
   CRISPOR + iGWOS; PAM GG=1.0/AG=0.259259/...; oracles perfect+GG->1.0, iGWOS doctests
   0.4635989007074176 / 0.5140384614450001). 71/71 off-target+CFD tests green; full suite 6812 passed/0
   failed, 0 warnings. Closed finding C7 (CFD was its last off-target residual). ENRICHED existing concept
   crispr-guide-rna-design (the CrisprDesigner/MolTools anchor) with a Layer-3 off-target section rather
   than creating a new concept (off-target is the complementary half of on-target guide selection on the
   same class): retitled to add off-target, added the report to sources, bumped source_commit to
   c763b50c, added a relates_to test-unit-registry typed edge for the off-target unit, updated intro +
   assumptions. New source-summary page crispr-off-001-report. Updated wiki/index.md (+1 source line,
   revised the concept line). Tied to validation-ledger / validation-and-testing / test-unit-registry; did
   NOT force the algorithm-validation-evidence hub. Kept distinct from the on-target sibling
   crispr-guide-001-report. No contradictions. Backlog slug off-target-analysis (Off_Target_Analysis.md)
   left pending — this report validates the MIT/CFD scoring layer, not that algorithm doc.
   graph: +1 node, +1 typed edge (crispr-guide-rna-design relates_to test-unit-registry for CRISPR-OFF-001)

## [2026-07-10] ingest | docs/Validation/reports/CRISPR-PAM-001.md → crispr-pam-001-report (source)
   Two-stage validation report for test unit CRISPR-PAM-001 (CRISPR PAM site detection — locate
   protospacer-adjacent motifs like SpCas9 NGG on both strands for gRNA design). Stage A PASS / Stage B
   PASS-WITH-NOTES / State CLEAN, 58/58 PAM tests green, zero production-code change. Validates
   CrisprDesigner.FindPamSites (string + DnaSequence overloads, identical results) + GetSystem: 7 systems
   with literature-confirmed PAM/orientation/guide-len (SpCas9 NGG/3'/20, SpCas9-NAG NAG, SaCas9
   NNGRRT(R in A,G)/21, Cas12a/AsCas12a TTTV(V in A,C,G)/5'/23, LbCas12a TTTV/24, CasX/Cas12e TTCN/20;
   Jinek 2012/Hsu 2013/Ran 2015/Zetsche 2015/Liu 2019). IUPAC via IupacHelper.MatchesIupac (NC-IUB 1984).
   Both-strand scan (forward CCN = revcomp of reverse NGG); 0-based coords, Position always forward-strand
   (reverse hit Position=len-i-pamLen, PamSequence revcomp); boundary targetStart>=0 && targetEnd<len.
   Hand oracle CCAACGTACGT...(len31)->0 fwd + 1 rev hit Position0/PamSequence "CCA"/target20 (M8); overlap
   AGGTGG->@20,@23; TTTV rejects TTTT, NNGRRT rejects C at R. Sole Stage-B note = reverse-strand
   PamSite.TargetStart indexes the revcomp string (different coord system than forward Position;
   XML-documented at :1035-1041, spec-invariant-consistent) — not a bug. ENRICHED existing concept
   crispr-guide-rna-design with a Layer 0 (PAM geometry) section rather than a new concept (PAM finding is
   the front end of guide design on the same CrisprDesigner class): added report to sources, bumped
   source_commit to 13507add, added a relates_to test-unit-registry typed edge for CRISPR-PAM-001, updated
   intro + Validation-status. New source-summary page crispr-pam-001-report. Updated wiki/index.md (+1
   source line, revised the concept line). Tied to validation-ledger / validation-and-testing /
   test-unit-registry; did NOT force the algorithm-validation-evidence hub. Kept distinct from the
   on-target/off-target siblings. No contradictions. Backlog slug pam-site-detection
   (PAM_Site_Detection.md algorithm doc) left pending — this report validates the unit, not that doc.
   graph: +1 node, +1 typed edge (crispr-guide-rna-design relates_to test-unit-registry for CRISPR-PAM-001)

## [2026-07-10] ingest | docs/Validation/reports/DISORDER-LC-001.md → disorder-lc-001-report (source) + 1 concept
   Two-stage validation report for test unit DISORDER-LC-001 (protein low-complexity region detection —
   SEG / Wootton & Federhen). This is SEG low-complexity detection despite the DisorderPredictor host class
   and the DISORDER-LC prefix — a distinct algorithm from TOP-IDP intrinsic-disorder-prediction. Stage A
   PASS-WITH-NOTES / Stage B PASS / End state CLEAN-FIXED, full suite 6612/0, zero production-code change.
   Validates DisorderPredictor.PredictLowComplexityRegions (:497) + private CalculateShannonEntropy (:280) +
   ClassifyLowComplexityType. H=-Sum p_i*log2(p_i) bits/residue + defaults W=12/K1=2.2/K2=2.5 confirmed
   verbatim vs GCG/Weizmann SEG help + ncbi-seg manpage + NCBI blast_seg.c (kSegWindow/kSegLocut/kSegHicut,
   s_Entropy); two-stage trigger(<=K1)/greedy-extend(whole segment <=K2)/merge/minLength scan. Detector
   re-implemented independently in Python — every expected value reproduced (M1 26Q->(0,25); M6
   20Q+60spacer+20A->(0,34),(67,99) as genuine H=2.4939<=2.5<2.6077 crossings with spacer H=3.585 gap
   preserved; S3 12Q->(0,11)). Test gate PASS: 3 coverage gaps closed with exact sourced values (M8 A/L-rich
   label branch A=L=50%, M9 custom triggerWindow W=4, M10 custom extensionThreshold K2=2.0), fixture 18->21.
   Stage-A notes (documented, non-blocking) = Shannon-entropy trigger stands in for WF eq-3 multinomial
   complexity + P0 optimization (disclosed algorithm-doc simplification, can differ from NCBI SEG on
   mixed-complexity edges, matches on homopolymer/biased inputs) + cosmetic X-rich/X/Y-rich label; neither
   moves segment boundaries. ENRICHED existing concept protein-low-complexity-seg (the SEG anchor already
   covering DISORDER-LC-001 + PROTMOTIF-LC-001) rather than a new concept: added report to sources, bumped
   source_commit to c9ed6cf3, added a relates_to test-unit-registry typed edge sourced from the report,
   updated intro to cite the report + validation-ledger. New source-summary page disorder-lc-001-report.
   Updated wiki/index.md (+1 source line). Tied to validation-ledger / validation-and-testing /
   test-unit-registry; did NOT force the algorithm-validation-evidence hub. Kept distinct from the evidence
   artifact disorder-lc-001-evidence. No contradictions. Backlog: ProteinPred SEG algorithm doc already
   covered via protein-low-complexity-seg (no pending slug to move).
   graph: +1 node, +1 typed edge (protein-low-complexity-seg relates_to test-unit-registry for DISORDER-LC-001, sourced from the report)

## [2026-07-10] ingest | docs/Validation/reports/DISORDER-MORF-001.md → disorder-morf-001-report (source)
   Two-stage validation report for test unit DISORDER-MORF-001 (MoRF / Molecular Recognition Feature
   prediction — the "dip within disorder" heuristic). Stage A PASS-WITH-NOTES / Stage B PASS / End state
   CLEAN, full unfiltered suite 6609/0, dotnet build 0 errors (4 pre-existing NUnit-analyzer warnings in
   unrelated files), zero code change. Validates DisorderPredictor.PredictMoRFs(sequence, minLength=10,
   maxLength=70) (DisorderPredictor.cs:615-671) over PredictDisorder window-21 normalized-TOP-IDP scores
   (constants MoRFOrderThreshold=0.5 :578 / MoRFMinLength=10 :584 / MoRFMaxLength=70 :589). MoRF = maximal
   ordered dip (per-residue disorder d<0.5) of length 10-70 flanked both sides by d>=0.5, score
   (0.5-mean d)/0.5 clamped [0,1]; every constant source-traced — 0.5 (Cheng/Oldfield PMC2570644), 10-70
   band (Mohan 2006 PMID16935303 / Wikipedia), TOP-IDP (Campen 2008 PMC2676888). Validator re-derived the
   smoothed profile independently in Python from source TOP-IDP raw values (not the repo): 25P+30L+25P dips
   [29,50] len22 mean0.362033 score0.275934; 25P+30I+25P dips [28,51] len24 score0.399608; both reproduce the
   locked test values, plus the empty cases (40L no-flank, 40P no-dip, 25P+16L+25P dip-8<10, 25P+95L+25P
   dip-87>70, 15L+30P terminal) and S1 two-dip (29,50)+(89,110). HARD test gate PASS (exact coords + scores
   Within(1e-6), monotonicity/bounds invariants alongside not instead of exact values, no green-washing).
   Stage-A notes: N1 = one bounded assumption (exact dip flank/run-length parameters live in Oldfield 2005's
   paywalled Methods -> documented qualitative criterion; load-bearing constants all source-traceable, not a
   correctness error); N2 = spec-prose M1 coordinate nit "20-34" -> 29-50/len22/score0.275934 corrected
   (doc-only, no code/test change); N3 = iota-subtype naming confirmed correct vs Mohan 2006 verbatim
   (Wikipedia less precise). ENRICHED the existing concept morf-prediction-dip-in-disorder (already the MoRF
   anchor) rather than a new concept — MoRF prediction is genuinely distinct but already represented: added
   report to sources, bumped source_commit to dc13c70f + updated 2026-07-10, updated intro to cite the report
   + validation-ledger. New source-summary page disorder-morf-001-report. Updated wiki/index.md (+1 source
   line). Tied to validation-ledger / validation-and-testing / test-unit-registry; did NOT force the
   algorithm-validation-evidence hub. Kept distinct from the evidence artifact disorder-morf-001-evidence.
   No contradictions. Backlog: the ProteinPred MoRF_Prediction algorithm-doc slug (morf-prediction) is a
   separate algorithm-doc reconciliation item, not this validation report (excluded from coverage) — no
   pending slug to move.
   graph: +1 node (disorder-morf-001-report source), +0 typed edges (report ties to validation-ledger /
   validation-and-testing / test-unit-registry via auto-derived mentions/sourced_from; concept already
   carries relates_to test-unit-registry + depends_on intrinsic-disorder-prediction-top-idp)

## [2026-07-10] ingest | docs/Validation/reports/DISORDER-PRED-001.md → disorder-pred-001-report (source)
   Two-stage validation report for test unit DISORDER-PRED-001 — the CORE per-residue intrinsic-disorder
   predictor DisorderPredictor.PredictDisorder (the shared PredictDisorder anchor the whole protein-disorder
   family reads from). Stage A PASS / Stage B PASS / End state CLEAN, full DisorderPredictor test family
   113/113 green (focused DisorderPredictor_DisorderPrediction_Tests 22/22), NO code changed this session, no
   divergences at either stage. Validates PredictDisorder (DisorderPredictor.cs:190) -> CalculatePerResidueScores
   (:227 centered truncated window) -> CalculateDisorderScore (:255 normalize (prop-TopIdpMin)/TopIdpRange +
   average), threshold score>=disorderThreshold (:242); constants TopIdpMin=-0.884(W) / TopIdpMax=0.987(P) /
   TopIdpRange=1.871 / TopIdpCutoff=0.542, default window 21; bundled CalculateHydropathy (Kyte-Doolittle 1982).
   Score S(aa)=(TOP-IDP(aa)-(-0.884))/1.871 averaged over window, disordered when >=0.542. All 20 Table 2
   values + cutoff 0.542 + window 21 + normalization + Dunker 8/8/4 classification (order {W,C,F,I,Y,V,L,N} /
   disorder {A,R,G,Q,S,P,E,K} / ambiguous {H,M,T,D}) re-confirmed LIVE vs Campen 2008 PMC2676888 (WebFetch of
   PMC full text; prediction eq I=-(<TOP-IDP>-0.542), positive=>ordered) and Wikipedia/Dunker 2001 PMID11381529.
   Hand-computed normalized scores reproduce spec M8b: W->0.0, I->0.398/1.871=0.21272 (ordered), E->1.620/1.871=
   0.86585 (disordered), P->1.0 (disordered). Edge semantics defined/sourced: empty->zeroed (INV-7), unknown
   residues skipped in mean (poly-X->0.0), termini truncated, case via ToUpperInvariant, scores in [0,1]. Tests
   assert exact externally-confirmed values (20 propensities M8, normalized W=0/I=0.2127/E=0.8660/P=1.0 M8b/S1,
   poly-I 0.0 / poly-E,P 1.0 M4-M6, hydropathy I=4.5/W=-0.9/E=-3.5 C4) with tight tolerances, not tautologies;
   division guarded count>0. Honest scoping preserved: composition heuristic AUC~0.65-0.72 (vs IUPred2A
   0.75-0.80, flDPnn 0.85-0.90), NOT an energy/ML predictor — no false IUPred-grade claim. Prior cosmetic
   ranking-string fix (S before K) already present in code+Evidence, no new divergence. Findings: NONE at both
   stages. SURGICALLY updated the existing concept intrinsic-disorder-prediction-top-idp (already the core
   PredictDisorder anchor) rather than a new concept — the algorithm is already fully represented: added the
   report to sources, bumped source_commit to 920bd895 + updated 2026-07-10, added a one-line CLEAN verdict
   cross-link to disorder-pred-001-report + validation-ledger in the intro. New source-summary page
   disorder-pred-001-report. Updated wiki/index.md (+1 source line). Tied to validation-ledger /
   validation-and-testing / test-unit-registry; did NOT force the algorithm-validation-evidence hub. Kept
   distinct from the evidence artifact disorder-pred-001-evidence and the sibling reports disorder-lc-001-report
   / disorder-morf-001-report. No contradictions.
   graph: +1 node (disorder-pred-001-report source), +0 typed edges (report ties to validation-ledger /
   validation-and-testing / test-unit-registry via auto-derived mentions/sourced_from; concept already carries
   relates_to test-unit-registry + relates_to protein-low-complexity-seg — no new typed edges added)

## [2026-07-10] ingest | docs/Validation/reports/DISORDER-PROPENSITY-001.md → disorder-propensity-001-report (source)
   Two-stage validation report for test unit DISORDER-PROPENSITY-001 — the O(1) per-residue propensity
   primitives beneath the sliding-window predictor: the raw TOP-IDP Table-2 lookup GetDisorderPropensity +
   Dunker order/disorder/ambiguous classification (IsDisorderPromoting / DisorderPromotingAminoAcids /
   OrderPromotingAminoAcids / AmbiguousAminoAcids), DisorderPredictor.cs:680-712 over the DisorderPropensity
   dict :85-107 and the three sets :111-121. Stage A PASS / Stage B PASS / End state CLEAN, full UNFILTERED
   suite 6609 passed / 0 failed / 0 skipped (14-test canonical fixture green), NO code/test/spec change this
   session, no divergences. Scope note: the session prompt's sliding-window smoothing / mean-disorder profile
   is NOT this unit — it lives in PredictDisorder under DISORDER-PRED-001; per Registry
   (ALGORITHMS_CHECKLIST_V2.md 3658-3676) + Method Index (5070-5073) this unit is scoped to the O(1) lookup +
   classification, validated as registered. GetDisorderPropensity returns the RAW un-normalized value
   (W-0.884...P+0.987), NOT the [0,1] normalized Si PredictDisorder uses. All 20 Table-2 values match exactly
   vs Campen 2008 PMC2676888 (WebFetch, Table 2 verbatim), min W=-0.884 / max P=0.987 confirmed; Dunker 8/8/4
   sets disorder {A,R,G,Q,S,P,E,K} / order {W,C,F,I,Y,V,L,N} / ambiguous {H,M,T,D} match, union 20 pairwise
   disjoint (INV-4). localCIDER (sequenceParameters.py:221) corroborates the Campen citation (not the raw
   numbers); PubMed 18991772 / USF DigitalCommons / Bentham confirm the rank string. Pitfall confirmed +
   correctly handled: rendered rank string places ...Q,K,S,E,P but Table-2 values give S=0.341<K=0.586 (by
   value ...Q,S,K,E,P) — presentation artifact, numeric values authoritative, both locked. Edge semantics:
   scale defined for 20 standard residues only, unknown->0.0 (GetValueOrDefault contract, not source value),
   case-fold via ToUpperInvariant — both impl contracts, honestly registered. HARD test-quality gate PASS:
   fixture hardcodes sourced literals (not code echoes) M1/M2 exact values + M6/M7/M8 exact set equivalence +
   counts; M3 disorder->true / M4 order->false / M5 ambiguous->false exercise both predicate branches; exact
   Is.EqualTo(...).Within(1e-10) + Is.EquivalentTo+Count, no Greater/AtLeast/Contains where exact known, no
   widened tolerance, no skipped tests; all 5 members + both branches + edge/contract cases (unknown X/Z/B/*
   ->0.0, lowercase, sorted property order) covered; honest green on full unfiltered suite. IsDisorderPromoting
   <=> membership in DisorderPromotingAminoAcids (M10) consistent. Findings: NONE at both stages. SURGICALLY
   updated the existing concept intrinsic-disorder-prediction-top-idp (already the shared PredictDisorder
   anchor + already hosts the DISORDER-PROPENSITY-001 primitives section) rather than a new concept — the
   primitives are already fully represented: added the report to sources, bumped source_commit to 540cb0df,
   added a one-line CLEAN verdict cross-link to disorder-propensity-001-report + validation-ledger in the
   primitives section. New source-summary page disorder-propensity-001-report. Updated wiki/index.md (+1 source
   line). Tied to validation-ledger / validation-and-testing / test-unit-registry; did NOT force the
   algorithm-validation-evidence hub. Kept distinct from the evidence artifact disorder-propensity-001-evidence
   and the sibling reports disorder-pred-001-report / disorder-lc-001-report / disorder-morf-001-report. No
   contradictions.
   graph: +1 node (disorder-propensity-001-report source), +0 typed edges (report ties to validation-ledger /
   validation-and-testing / test-unit-registry via auto-derived mentions/sourced_from; concept already carries
   its typed relates_to edges — no new typed edges added)

## 2026-07-10 — ingest docs/Validation/reports/DISORDER-REGION-001.md
   Ingested the per-unit VALIDATION REPORT for DISORDER-REGION-001 (disordered-region detection — the
   segment-calling layer collapsing the per-residue PredictDisorder TOP-IDP profile into contiguous IDR
   regions, plus opt-in MobiDB-lite v3 flavour typing). Re-validated 2026-06-25, supersedes 2026-06-24 pass
   (reset to pending after F4 ClassifyRegionFlavorMobiDbLite added). Verdict Stage A PASS-WITH-NOTES / Stage B
   PASS / End state CLEAN. Region-calling logic (IdentifyDisorderedRegions DisorderPredictor.cs:358, single-pass
   with explicit trailing-run branch -> no off-by-one; ClassifyDisorderedRegion :428 / CalculateConfidence
   :470) validated correct + sourced: consecutive-grouping score>=0.542, configurable minRegionLength default 5,
   0-based inclusive Start/End length=End-Start+1, Long IDR label at length>30 anchored to Ward 2004 DISOPRED2
   (>30 consecutive residues, 2.0% archaea/4.2% bacteria/33% eukaryotes). Independent Python recompute
   reproduced M2 P30->[0,29] / M6 W10+P20->[11,29] / S2 P20+W30->[0,18] / S5 W15+P20+W15->[16,33] / M14 two
   regions [(16,33),(51,68)] + homopolymer MeanScore/Confidence P30 1.0/1.0, E30 0.866/0.707, K30 0.786/0.532,
   S30 0.655/0.246 — all matching test assertions. Stage-B fidelity DEFECT found and FULLY FIXED this session:
   ClassifyRegionFlavorMobiDbLite computed f+ over {R,K} only, but MobiDB-lite v3 states.py translation table
   maps H->positive (f+=(R+K+H)/L) — HHHHHHHHAA (f+=0.8) returned WeaklyCharged instead of
   PositivePolyelectrolyte; the 16 prior flavour tests were all His-free so the gap was invisible. Fix: added H
   to positive count, corrected docstring (translation table cited), added F5b (HHHHHHHHAA->PPE) + F5c
   (HHHHDDDDAA->Polyampholyte, H balancing D); boundaries + default RegionType/Confidence untouched. Notes
   (PASS-WITH-NOTES) = disclosed first-principles labelling heuristics (enrichment 0.25 ~5x-random NOT Das&Pappu
   NCPR, priority Pro>Acidic>Basic>S/T>Long>Standard, Confidence=(mean-0.542)/(1-0.542)) affecting labels not
   boundaries, + no citable per-region confidence standard (MobiDB-lite/IUPred/PONDR report per-residue only;
   LimitationPolicy-guarded branch min access Permissive, throws under default Moderate); cosmetic Campen comment
   S/K-swap nit. Boundary fixture 24/24, flavour fixture 16/16 (incl. F5b/F5c), full unfiltered dotnet test
   Failed 0 (Genomics 18819). New source-summary page disorder-region-001-report. SURGICALLY enriched the
   existing anchor concept intrinsic-disorder-prediction-top-idp (already hosts the DISORDER-REGION-001
   region-detection section + already linked the evidence) rather than a new concept — added the report to
   sources, bumped source_commit to 8bd6a5e1, cross-linked disorder-region-001-report + the CLEAN/one-fix
   verdict + validation-ledger in the region-detection section. Updated wiki/index.md (+1 source line). Tied to
   validation-ledger / validation-and-testing / test-unit-registry; did NOT force the
   algorithm-validation-evidence hub. Kept distinct from the evidence artifact disorder-region-001-evidence and
   the sibling reports disorder-pred-001-report / disorder-propensity-001-report / disorder-lc-001-report /
   disorder-morf-001-report. No contradictions (the 0.25-vs-NCPR clarification is a disclosed heuristic
   distinction, already recorded on the evidence page, not a new source contradiction).
   graph: +1 node (disorder-region-001-report source), +0 typed edges (report ties to validation-ledger /
   validation-and-testing / test-unit-registry via auto-derived mentions/sourced_from; concept already carries
   its typed relates_to edges — no new typed edges added)
- 2026-07-10 — Ingested `docs/Validation/reports/EPIGEN-AGE-001.md` (per-unit two-stage validation report for EPIGEN-AGE-001, epigenetic age / DNAm clocks — Horvath 2013 multi-tissue 353-CpG + Horvath 2018 skin&blood 391-CpG + Levine 2018 PhenoAge 513-CpG; `EpigeneticsAnalyzer.CalculateEpigeneticAge`/`CalculateSkinBloodAge`/`CalculatePhenoAge`). Verdict Stage A PASS / Stage B PASS / End state CLEAN — no defect, no code change; 34/34 filtered tests pass; embedded 391 + 513 coefficient tables numerically identical to biolearn (`Horvath2.csv`/`PhenoAge.csv`, max abs diff 0.0), multi-tissue 353 byte-identical; anti.trafo (adult.age=20) byte-for-byte, strict `<` puts Y=0 on linear branch → 20; PhenoAge untransformed, locked by PA6 negative control (≠ anti.trafo(60.664)=1293.944). New source-summary page epigen-age-001-report. SURGICALLY updated the existing anchor concept epigenetic-age-horvath-clock (added the report to sources + bumped source_commit e90a7598→4416f109, cross-linked epigen-age-001-report + validation-ledger in the intro, and enriched the Scope section: report independently confirmed all three clocks vs biolearn). Updated wiki/index.md (+1 source line). Tied to validation-ledger / validation-and-testing / test-unit-registry; did NOT force the algorithm-validation-evidence hub. Kept distinct from the evidence artifact epigen-age-001-evidence. No contradictions. Backlog: EPIGEN-AGE-001's algorithm doc (`Epigenetic_Age_Estimation.md`) already covered-via-concept — no backlog move needed.
   graph: +1 node (epigen-age-001-report source), +0 typed edges (report ties to validation-ledger / validation-and-testing / test-unit-registry via auto-derived mentions/sourced_from; concept already carries its typed relates_to edge — no new typed edges added)

## [2026-07-10] lint | structural + staleness + coverage + graph pass (488→491 pages)
Ran the full lint suite. Staleness CLEAN, graph lint CLEAN (491 nodes, 0 issues). Findings triaged and fixed:
- **Broken wikilinks (4)** — false positives from escaped-pipe table cells (`[[slug\|Alias]]`) on rna-secondary-structure-prediction / rna-stem-loop-enumeration; all targets exist. Root cause was `WIKILINK_RE` in `scripts/wiki_lint.py` capturing the trailing `\` into the slug. Fixed the regex to exclude `\` from the target class and tolerate an optional leading backslash before the alias pipe (verified against all 4 flagged strings + plain-pipe / no-alias / mixed-row controls). Linter-only change; no wiki page edited for this.
- **Oversize (2)** — backlog.md (448) and algorithm-validation-evidence.md (466), both over the 400 soft cap and under the 800 hard cap; ACCEPTED as meta/hub pages (user decision), revisit near 800.
- **Coverage** — 174 uncovered → 168 (all `docs/algorithms/**`, backlog-tracked) after handling the 6 non-algorithm docs: ingested the 3 overviews (new source pages mcp-readme, skills-strategy, golden-skills-regression, cross-linked to each other + skill-layer / scientific-rigor / mcp-plan) and added 3 SCHEMA excludes for the ledgers (docs/mcp/MCP_STATUS.md, docs/mcp/traceability.md, docs/skills/golden/tasks.md). Zero non-algorithm docs now uncovered.
- Updated wiki/index.md (+3 source lines), wiki/concepts/skill-layer.md (wikilinked the strategy/golden/mcp-readme refs), and regenerated the compiled graph (491 nodes, 3734 edges).
   graph: +3 nodes (mcp-readme, skills-strategy, golden-skills-regression sources), +0 typed edges (new source pages carry mentions only; no `graph.relationships` blocks added).

## [2026-07-10] ingest | MCP tool catalog (batch, replaces per-tool ingestion)
Added `wiki/concepts/mcp-tool-catalog.md` (354 lines) — one durable reference mapping all **427 MCP tools across 11 servers** to the concept each thin wrapper delegates to (deterministic map from `tools/wiki-ingest/mcp_map.py`). Concept-grouped rows (`[[concept]] — tool_a, tool_b`) keep every tool→concept visible while staying compact; unmapped tools grouped by wrapped C# class per server. Coverage: **209 mapped / 218 unmapped**, **120 distinct concepts** referenced (all resolve; no broken wikilinks). No new concept page created — unmapped tools recorded as gaps (largest clusters: Parsers 39/41 whole-server gap, PrimerDesigner, GenomeAssemblyAnalyzer, StructuralVariantAnalyzer, RestrictionAnalyzer, RnaSecondaryStructure energy terms). Typed graph: catalog `relates_to` three-front-doors + skill-layer (both high, README-sourced). Updated wiki/index.md (+1 Concepts line).
- **Validation reports (12)** — ingested the 12 remaining `docs/Validation/reports/**` as full `[[<unit>-report]]` verdict pages (GENOMIC-REPEAT-001, META-RESIST-001, MIRNA-PAIR-001, MIRNA-TARGET-001, PANGEN-CORE-001, PARSE-GENBANK-001, PROTMOTIF-TM-001, RNA-PARTITION-001, SEQ-COMPLEX-COMPRESS-001, SV-CNV-001, SV-DETECT-001, TRANS-SPLICE-001), each wired into its existing concept (`sources:` + body wikilink) and listed under index `## Sources`; all end-state ✅ CLEAN. Complements the bulk `validation-verdicts` registry (these 12 were not in it). Also removed a stray tool-call tag block at the end of `discordant-pair-sv-detection.md`.
- **Oversize split (2)** — split the two soft-cap pages: `backlog.md` (448→149) moved its per-domain pending tables to the new [[backlog-pending]] (type index); `concepts/algorithm-validation-evidence.md` (468→376) moved its inline 190-item `[[..-evidence]]` enumeration to the new [[evidence-artifact-index]] (type synthesis, so the links still count as inbound — no evidence page orphaned). Both new pages indexed; wiki lint clean (0 oversize).

## [2026-07-11] ingest | docs/algorithms/Annotation/GFF3_IO.md → gff3-io (concept)
   New concept [[gff3-io]] synthesizing the ANNOT-GFF-001 **annotation-layer** GFF3 helper
   (`GenomeAnnotator.ParseGff3` / `ToGff3` / `EncodeGff3Value`, Implementation Status *Simplified*) —
   the second GFF3 path, kept distinct from the FileIO `GffParser` unit ([[parse-gff-001-evidence]]).
   Carries what is distinct: reduced `GenomicFeature` parse record (drops `seqid`/`source` → lossy
   round-trip), `GeneAnnotation`-only exporter, col-9 encoder, and the load-bearing **per-transcript
   cumulative CDS phase** `(3 − Σ preceding-lengths mod 3) mod 3` (5′→3′, both strands). Reused
   existing pages rather than duplicating the 9-column schema/dialect facts (on [[parse-gff-001-evidence]])
   and coordinate contrast (on [[bed-format-parsing]] / [[insdc-feature-location]]). Cross-linked from
   [[annot-gff-001-report]] (validation verdict), [[bed-format-parsing]], [[parse-gff-001-evidence]].
   Re-homed the `parse_gff3` / `to_gff3` MCP tools from [[bed-format-parsing]] onto [[gff3-io]] in
   [[mcp-tool-catalog]] (mapped 209→210, unmapped 218→217, concepts 120→121). Moved GFF3_IO from
   [[backlog-pending]] (Annotation 3→2) to [[backlog]] covered (77→78). Updated index. No contradictions.
   graph: +1 node, +1 typed edge

## [2026-07-13] ingest | docs/algorithms/Annotation/Gene_Prediction.md → prokaryotic-gene-prediction-rbs (concept, reconciliation)
   Reconciliation ingest: the primary spec was already fully synthesized by the existing concept
   [[prokaryotic-gene-prediction-rbs]] (created 2026-07-10, alongside the ANNOT-GENE-001 validation
   ingest), which lists `docs/algorithms/Annotation/Gene_Prediction.md` in `sources:` and is
   NOT stale (its `source_commit` == HEAD `ec9209f6`; the doc has no commits since). Verified the
   concept faithfully covers the doc's contract, INV-01…04, declared simplifications (no overlap/
   best-model resolution, promoter −10/−35 not integrated → [[promoter-detection]], forward-only
   legacy RBS, length-only score), SD model, complexity, and edge cases — no enrichment needed.
   Per the established project pattern, `docs/algorithms/**` primary specs are synthesized into
   concepts and cited in `sources:` rather than given a dedicated `wiki/sources/` page (no algorithm
   spec has one). Sole gap was the stale backlog (generated 2026-07-09, before the concept existed):
   moved Gene_Prediction.md from [[backlog-pending]] (Annotation 2→1) to [[backlog]] covered (78→79,
   pending 167→166). No new pages, no new graph edges, no contradictions.

## [2026-07-13] ingest | docs/algorithms/Annotation/ORF_Detection.md → [[open-reading-frame-detection]] (concept, enriched)
   Reconciliation ingest of the annotation-layer ORF primary spec (test unit ANNOT-ORF-001,
   `GenomeAnnotator.FindOrfs`). The existing concept [[open-reading-frame-detection]] already
   synthesized the sibling annotation finder (ATG/GTG/TTG starts, amino-acid `minLength`,
   `searchBothStrands`/`requireStartCodon` flags, `Translator.FindOrfs`, hard-coded standard code,
   nested-start-per-stop semantics) and cited its validation report [[annot-orf-001-report]] — but
   did NOT list the spec in `sources:` and lacked `FindLongestOrfsPerFrame`. Enriched the "sibling
   annotation ORF finder" section with `FindLongestOrfsPerFrame` (calls `FindOrfs(minLength:1)`,
   groups under signed frame keys 1..3 forward / -1..-3 reverse, longest-per-frame), default
   minLength 100, 0-based half-open `[Start,End)` span with terminal `*` retained in `ProteinSequence`,
   and pending-starts-per-frame; removed the now-false "not ingested here" note. Added
   `docs/algorithms/Annotation/ORF_Detection.md` to `sources:`, bumped `source_commit` → HEAD
   `4b5d552f`. Moved the row from [[backlog-pending]] (Annotation 1→0, section removed) to [[backlog]]
   covered (79→80, pending 166→165, 27→26 domains). Per the established pattern, no dedicated
   `wiki/sources/` page for the algorithm spec. No new concept, no new graph edges, no contradictions.

## [2026-07-13] ingest | docs/algorithms/CANONICAL_MAP.md → [[canonical-algorithm-map]] (source, meta/index)
   Ingested the META/INDEX doc (not a per-algorithm spec) as ONE concise source page
   `wiki/sources/canonical-algorithm-map.md` recording the project's algorithm-identity authority:
   alias→canonical test-unit IDs (SEQ-COMPOSITION-001→SEQ-STATS-001, SEQ-TM-001→SEQ-THERMO-001,
   GENOMIC-TANDEM-001→REP-TANDEM-001), folder-bucket normalization (MolTools / Population_Genetics /
   K-mer / RnaStructure — four merges complete; Motif_Analysis/Sequence_Comparison/Genomic_Analysis/
   Extended_Annotation/Extended_Assembly pending), one-canonical-doc-per-concept rules, and the
   retained legacy/baseline methods (UPGMA, JC/K2P, chi-sq HWE, Nussinov, OLC). Cross-linked to the
   coverage/registry infrastructure: [[backlog]] (identity-vs-coverage complement),
   [[algorithms-checklist-v2]], [[algorithm-validation-evidence]], [[validation-ledger]]. Listed under
   index Sources; updated the [[backlog]] Notes row (CANONICAL_MAP now points at the new page,
   README.md remains index-only). Per the meta-doc rule: NO per-algorithm concept pages created. No
   typed graph edges (source page). No contradictions.

## [2026-07-13] ingest | docs/algorithms/Chromosome_Analysis/Higher_Order_Repeat_Detection.md → enriched [[centromere-analysis]]
   Reused the existing chromosome centromere/satellite anchor rather than creating a new page: the
   HOR method (`DetectHigherOrderRepeat`) was already synthesized there (validated as CHROM-HOR-001).
   Added the primary spec to `sources:` (source_commit → HEAD b8c0053), enriched the HOR bullet with
   genuinely-distinct spec detail: the `HorResult` record struct, `monomerLength` default/validation,
   smallest-qualifying-period rule, O(M²·L²) cost with per-pair memoisation, and the not-implemented
   cascading/nested HOR-of-HORs decomposition (HORmon / alpha-CENTAURI). Resolved the last pending
   Chromosome_Analysis backlog row: moved to Covered via concept, updated counts (80→81 covered,
   165→164 pending, 26→25 domains) and removed the section from [[backlog-pending]]. No new page,
   no new typed edges (existing chrom-hor-001-report edge already links this slice). No contradictions.

## [2026-07-13] ingest | docs/algorithms/Codon/Codon_Usage_Statistics.md → new concept [[codon-usage-statistics]]
   Created a dedicated concept for the codon-usage family's aggregation/reporting view
   (CODON-STATS-001, `CodonUsageAnalyzer.GetStatistics` + `CalculateCai`). Rationale for a new page
   rather than enrichment: the component measures (RSCU/ENC/CAI) are each owned by existing concepts,
   but the **positional-GC block (GC1/GC2/GC3, GC3s, OverallGc)** and the one-pass bundle contract are
   genuinely distinct and unowned. Page synthesizes the bundle contract + positional GC — GC3s =
   %G/C at synonymous third positions with Met/Trp/stop excluded (Peden 1999 §1.8.2.1.3, the ATGGCA
   GC3=50% vs GC3s=0% contrast), GC3s-as-percentage (100× CodonW fraction) deviation, and the
   skip-zero-w CAI note (vs the [[codon-adaptation-index]] 1e-6 clamp; bundled E.coli/Kazusa tables
   never hit it). The CODON-STATS-001 Evidence and Validation-report artifacts were already ingested
   as source pages ([[codon-stats-001-evidence]] / [[codon-stats-001-report]]); this closes the
   remaining *algorithm-spec* coverage gap. Inbound links added from [[relative-synonymous-codon-usage]]
   and [[codon-usage-comparison]] (retargeted their "aggregation view" mentions from the evidence
   source page to the concept). Moved the last pending Codon backlog row to Covered via concept;
   updated counts (81→82 covered, 164→163 pending, 25→24 domains) and removed the Codon section from
   [[backlog-pending]]. Added [[codon-usage-statistics]] to index Concepts. No contradictions.
   graph: +1 node, +4 typed edges
- 2026-07-13 — ingest docs/algorithms/Complexity/DUST_Score.md (DUST Score algorithm SPEC — the primary per-algorithm spec for the triplet-frequency low-complexity masking score `S(x)=∑_t c(c−1)/2 / (L−2)`, unit SEQ-COMPLEX-DUST-001, status *Simplified*). CONTEXT check: the concept [[dust-low-complexity-score]] already exists (created 2026-07-10 from the SEQ-COMPLEX-DUST-001 Evidence artifact) and already synthesizes the formula, defaults (k=3, window 64, threshold 2.0/level 20), worked oracles, and corner cases. Per the ingest brief (prefer REUSE; do NOT create a redundant concept for a spec already covered), ENRICHED the existing concept rather than creating a new page. Added the doc's distinct **implementation surface** not previously captured: entry points `CalculateDustScore(DnaSequence,int)`/`(string,int)` + `MaskLowComplexity` in `SequenceComplexity.cs`; complexity O(L·wordSize) time / O(min(L,4^wordSize)) space; the "no suffix tree used" reuse-policy note; the **normalization-bug correction** (earlier code divided by `words−1`, over-scaling; fixed to divide by the word count for the `1/(L−2)` normalization); and the **SDUST symmetric perfect-interval masking rule is NOT implemented** (MaskLowComplexity is a fixed-window threshold scan, so masked boundaries may differ from dustmasker/sdust). Frontmatter: added docs/algorithms/Complexity/DUST_Score.md to sources (now spec + Evidence), source_commit→HEAD (be6d15b), updated→2026-07-13. Backlog reconciliation: moved the `docs/algorithms/Complexity/DUST_Score.md` row from [[backlog-pending]] (Complexity 4→3) to [[backlog]] Covered-via-concept table → [[dust-low-complexity-score]]; counts 82→83 covered, 163→162 pending. No new page ⇒ no index.md change (concept already listed). No graph change (no new nodes/edges — concept and its registry edge already exist). Sources agree with the already-recorded Evidence synthesis (Morgulis 2006, Li 2025 longdust, lh3/sdust); one new spec detail is the divisor-correction history, consistent with the `1/(L−2)` normalization already on the page — no contradictions. Follow-ups: none — the DUST algorithm is now covered by both its Evidence and its primary spec.
- 2026-07-13 — ingest docs/algorithms/Complexity/K-mer_Entropy.md (K-mer Entropy algorithm SPEC — the primary per-algorithm spec for `SequenceComplexity.CalculateKmerEntropy`, unit SEQ-COMPLEX-KMER-001, status *Production*; Shannon entropy `H=−Σ pᵢ log₂ pᵢ`, `pᵢ=nᵢ/(L−k+1)`, bits, over the overlapping k-mer distribution). CONTEXT check: the concept [[k-mer-statistics]] already synthesizes this exact unit in its "Second entry point" section (created from the SEQ-COMPLEX-KMER-001 Evidence artifact) — formula, contract, bounds, worked oracles (ACGT k=1/k=2, ATATAT, AAAA, AAACGT) all already present. Per the ingest brief (prefer REUSE; do NOT create a redundant concept for a spec already covered), ENRICHED the existing concept rather than creating a new `k-mer-entropy` page. Added the doc's distinct **implementation surface** not previously captured: entry points `CalculateKmerEntropy(DnaSequence,int)`/`(string,int)` → private `CalculateKmerEntropyCore(string,int)` in `SequenceComplexity.cs`; single linear scan with `Dictionary<string,int>`; complexity **O(N·k) time, O(D·k) space**; the "suffix tree evaluated and not used" reuse-policy note (mirrors DUST); the **generalization-of per-base Shannon entropy** comparison (`CalculateShannonEntropy` ≈ k=1 over 4 nt, max 2 bits, order-blind); the default **k=2** and unconstrained/no-IUPAC alphabet; and the **Not-implemented** items — normalised entropy `H/log₂N` and the entropy–rank ratio (Çakır arXiv:2511.05300), so cross-length values are not directly comparable and RC-equivalence/position are not modelled. Frontmatter: added docs/algorithms/Complexity/K-mer_Entropy.md to sources (now Evidence + K-mer_Statistics spec + this spec), source_commit→HEAD (1166783), updated→2026-07-13. Backlog reconciliation: moved the `docs/algorithms/Complexity/K-mer_Entropy.md` row from [[backlog-pending]] (Complexity 3→2) to [[backlog]] Covered-via-concept → [[k-mer-statistics]]; counts 83→84 covered, 162→161 pending. No new page ⇒ no index.md change (concept already listed). No graph change (no new nodes/edges — concept + registry edges already exist; the pending `shannon-entropy` reference kept as plain text, not a broken wikilink). Sources agree with the already-recorded Evidence synthesis (Li 2025 longdust, Entropy–Rank Ratio, Shannon 1948) — no contradictions. Follow-ups: none — k-mer entropy is now covered by both its Evidence and its primary spec.
- 2026-07-13 — ingest docs/algorithms/Complexity/Lempel_Ziv_Complexity.md (Lempel–Ziv Complexity algorithm SPEC — the primary per-algorithm spec for LZ76 complexity `c(S)`, its length-normalization `c/(n/log_b n)`, and `EstimateCompressionRatio`, unit SEQ-COMPLEX-COMPRESS-001, status *Production*). CONTEXT check: the concept [[sequence-complexity-compression-lempel-ziv]] already exists (created 2026-07-10 from the SEQ-COMPLEX-COMPRESS-001 Evidence + Validation-report artifacts) and already synthesizes the parsing rule, raw-count oracles, normalization worked example, homopolymer/b<2/trailing-component corner cases, and the family placement vs DUST/k-mer-entropy/SEG/repeat-detection. Per the ingest brief (prefer REUSE; do NOT create a redundant concept for a spec already covered), ENRICHED the existing concept rather than creating a new `lempel-ziv-complexity` page. Added the doc's distinct **implementation surface** not previously captured as a new "Implementation surface and spec cross-check" section: the three entry points `CalculateLempelZivComplexity`→int, `CalculateNormalizedLempelZivComplexity`→double, `EstimateCompressionRatio`→double thin delegate (registry-canonical name behind the `complexity_compression_ratio` MCP tool, INV-05) on `SequenceComplexity.cs` (Seqeron.Genomics.Analysis); the `HashSet<string>` single-pass factorization with the "suffix tree evaluated and not used" reuse note; complexity **O(n) amortized / O(n²) worst case, O(n) space**; and null/empty handling (null DnaSequence throws, null/empty string→0). Frontmatter: added docs/algorithms/Complexity/Lempel_Ziv_Complexity.md to sources (now spec + Evidence + Validation report), source_commit→HEAD (6cc4dea), updated→2026-07-13. Backlog reconciliation: moved the row from [[backlog-pending]] (Complexity 2→1) to [[backlog]] Covered-via-concept → [[sequence-complexity-compression-lempel-ziv]]; counts 84→85 covered, 161→160 pending (24 domains unchanged). No new page ⇒ no index.md change (concept already listed). No graph change (no new nodes/edges — concept + its registry edge already exist). CONTRADICTION flagged on the page: the spec (last reviewed 2026-06-14) still documents the **raw count** being returned for single-symbol input (b<2), but validation (2026-07-10) corrected the code to the reference clamp-to-2 normalized value (`"0"×16 → 1.25`); the read-only spec now lags shipped behavior on that one degenerate case, validated 1.25 is authoritative. Follow-ups: none — LZ76 complexity is now covered by its spec, Evidence, and Validation report.
- 2026-07-13 — ingest docs/algorithms/Complexity/Windowed_Complexity.md (Windowed Sequence Complexity algorithm SPEC — the primary per-algorithm spec for `SequenceComplexity.CalculateWindowedComplexity`, unit SEQ-COMPLEX-WINDOW-001, status *Production*; sliding-window driver emitting one `ComplexityPoint` per fully-contained window carrying per-window Shannon entropy + summation-form linguistic complexity). CONTEXT check: the concept [[windowed-sequence-complexity-profile]] already exists (created 2026-07-10 from the SEQ-COMPLEX-WINDOW-001 + SEQ-ENTROPY-PROFILE-001 Evidence artifacts) and already synthesizes the two per-window measures, window geometry / `ComplexityPoint` contract, worked oracles (ACGTACGT LC=23/29, AAAAAAAA LC=6/29, the L=24/w=8/s=8 → 3-window geometry), corner cases, and the second standalone Shannon-entropy-profile entry point. Per the ingest brief (prefer REUSE; do NOT create a redundant concept for a spec already covered), ENRICHED the existing concept rather than creating a new `windowed-complexity` page. Added the doc's distinct **implementation surface** not previously captured, as a new "Implementation surface" section: entry point `CalculateWindowedComplexity(DnaSequence, windowSize=64, stepSize=10)` on `SequenceComplexity.cs` (Seqeron.Genomics.Analysis) returning a lazily-evaluated `IEnumerable<ComplexityPoint>`; delegation to the shared `CalculateShannonEntropyCore`/`CalculateLinguisticComplexityCore` helpers (window values equal the standalone scalar metrics); complexity **O((L/s)·w²) time**, per-window space bounded by distinct-subword count; the "suffix tree evaluated and rejected" reuse note; and the **Not-implemented** linear-time suffix-tree profile of Troyanskaya et al. (2002) (direct per-window enumeration is authoritative for the ≤6 bounded word lengths). Frontmatter: added docs/algorithms/Complexity/Windowed_Complexity.md to sources (now spec + 2 Evidence), source_commit→HEAD (579178e), updated→2026-07-13. Backlog reconciliation: moved the row from [[backlog-pending]] (Complexity 1→0, section removed) to [[backlog]] Covered-via-concept → [[windowed-sequence-complexity-profile]]; counts 85→86 covered, 160→159 pending, 24→23 domains. No new page ⇒ no index.md change (concept already listed). No graph change (no new nodes/edges — concept + its two registry edges already exist). Spec agrees with the already-recorded Evidence synthesis (Shannon 1948, Troyanskaya 2002, Gabrielian & Bolshoy 1999) — no contradictions. Follow-ups: none — windowed complexity is now covered by both its Evidence and its primary spec; this closes the Complexity domain in the pending backlog.
- 2026-07-13 — ingest docs/algorithms/Extended_GC_Skew_Analysis/AT_Skew.md (AT Skew algorithm SPEC — the primary per-algorithm spec for `GcSkewCalculator.CalculateAtSkew`, unit SEQ-ATSKEW-001, status *Production*; the scalar strand-composition statistic `(A−T)/(A+T)`). CONTEXT check per brief (checked for existing GC-skew / AT-skew / strand-asymmetry / replication-origin concepts first): the concept [[nucleotide-composition-skew]] already exists (created 2026-07-10 from the SEQ-ATSKEW-001 / SEQ-STATS-001 / SEQ-GC-ANALYSIS-001 / SEQ-REPLICATION-001 Evidence artifacts) and already synthesizes the formula, `[−1,+1]` bounds, zero-denominator⇒0.0 convention, case-insensitive/ignore-non-A/T symbol handling, worked oracles, and the strand-asymmetry biology (Lobry 1996, Charneski 2011) with cross-links to [[replication-origin-cumulative-skew]] and [[windowed-gc-profile-and-variance]]. Per the brief (AT skew is the compositional sibling of GC skew — fold into the existing concept, do NOT create a redundant page), ENRICHED the existing concept rather than creating a new `at-skew` page. Added the doc's distinct **implementation surface** not previously captured, as a new "Implementation" section: both scalar skews live in `GcSkewCalculator` (Seqeron.Genomics.Analysis); AT skew has two overloads over a shared private core — `CalculateAtSkew(string)` (canonical; upper-cases, counts A/T via `string.Count`; null/empty⇒0) and `CalculateAtSkew(DnaSequence)` (forwards the normalized value object; null⇒ArgumentNullException); complexity **O(n) time / O(1) space**; the "suffix tree does not apply — a two-symbol count, not substring search" reuse-policy note; and the **global-scalar-only** scope (windowed/cumulative AT-skew profiles and AT-skew-based origin location deliberately not implemented → the same class's GC-skew-based `CalculateWindowedGcSkew`/`CalculateCumulativeGcSkew`/`PredictReplicationOrigin` for localization). Frontmatter: added docs/algorithms/Extended_GC_Skew_Analysis/AT_Skew.md to sources (now spec + 4 Evidence), source_commit→HEAD (4beb586), updated→2026-07-13. Backlog reconciliation: moved the row from [[backlog-pending]] (Extended_GC_Skew_Analysis 2→1) to [[backlog]] Covered-via-concept → [[nucleotide-composition-skew]]; counts 86→87 covered, 159→158 pending (23 domains unchanged — Comprehensive_GC_Analysis still pending in the domain). No new page ⇒ no index.md change (concept already listed). No graph change (no new nodes/edges — concept + its registry edge already exist). Spec agrees with the already-recorded Evidence synthesis (Lobry 1996, Charneski 2011, Wikipedia "GC skew", Biopython `GC_skew`) — no contradictions. Follow-ups: the sibling `docs/algorithms/Extended_GC_Skew_Analysis/Comprehensive_GC_Analysis.md` (expected slug `comprehensive-gc-analysis`, likely folds into [[windowed-gc-profile-and-variance]]) remains the one pending row in the domain.
- 2026-07-13 — ingest docs/algorithms/Extended_GC_Skew_Analysis/Comprehensive_GC_Analysis.md (Comprehensive GC Analysis algorithm SPEC — the primary per-algorithm spec for `GcSkewCalculator.AnalyzeGcContent`, unit SEQ-GC-ANALYSIS-001, status *Production*; the composite that bundles whole-sequence OverallGcContent `(G+C)/(A+T+G+C)×100` + OverallGcSkew `(G−C)/(G+C)` + OverallAtSkew `(A−T)/(A+T)` with a sliding-window GC%/GC-skew profile AND the population variance `σ²=Σ(xᵢ−μ)²/N` of each windowed series, returned as one `GcAnalysisResult` record). CONTEXT check per brief (checked existing GC-content/GC-skew/windowed concepts first): the concept [[windowed-gc-profile-and-variance]] already exists (created 2026-07-10 from the SEQ-GC-ANALYSIS-001 Evidence artifact) and already synthesizes the six outputs, formulas, window geometry (only fully-contained windows; shorter-than-window ⇒ empty lists), the ÷N population-variance choice, worked oracles (`GGGCCAT`→71.428…/0.2/0.0; `GGCC` w2s2→`GcSkewVariance` 1.0/`GcContentVariance` 0.0), and the two labelling/estimator assumptions. Per the brief (this is an aggregation/bundle folding into an existing GC concept — enrich, do NOT create a redundant page), ENRICHED the existing concept rather than creating a `comprehensive-gc-analysis` page; exact precedent = the 2026-07-13 AT_Skew spec ingest into [[nucleotide-composition-skew]]. Added the doc's distinct **implementation surface** not previously captured, as a new "Implementation" section: lives in `GcSkewCalculator` (Seqeron.Genomics.Analysis); one entry point two overloads over a shared core — `AnalyzeGcContent(DnaSequence, windowSize=1000, stepSize=100)` (canonical; validates non-null → `ArgumentNullException`; delegates to core) and `AnalyzeGcContent(string, …)` (parity; null/empty ⇒ zero result, empty lists, `SequenceLength=0`); the eight-field `GcAnalysisResult` record shape; window positions 0-based, `WindowStart`/`WindowEnd` inclusive, `Position`=midpoint `start+windowSize/2`, count `⌊(n−w)/step⌋+1` (INV-05); complexity **O(n + W·w) time / O(W) space** with each window recounted independently (no incremental sliding accumulator, matching the `CalculateWindowedGcSkew` cores); the suffix-tree-not-applicable reuse note (counting/aggregation only); and the out-of-scope note (no cumulative diagram / origin call — those belong to `PredictReplicationOrigin` SEQ-REPLICATION-001 / `CalculateCumulativeGcSkew`). Frontmatter: added docs/algorithms/Extended_GC_Skew_Analysis/Comprehensive_GC_Analysis.md to sources (now spec + 3 Evidence), source_commit→HEAD (7b4fd24), updated→2026-07-13. Backlog reconciliation: moved the row from [[backlog-pending]] to [[backlog]] Covered-via-concept → [[windowed-gc-profile-and-variance]]; removed the now-empty Extended_GC_Skew_Analysis pending section; counts 87→88 covered, 158→157 pending, 23→22 domains (Extended_GC_Skew_Analysis fully covered). No new page ⇒ no index.md change (concept already listed). No graph change (concept + its registry/base-composition/nucleotide-composition-skew edges already exist; spec adds no new typed edges). Hub [[algorithm-validation-evidence]] unchanged (already records SEQ-GC-ANALYSIS-001; spec docs are not added to the hub, per AT_Skew precedent). Spec agrees with the already-recorded Evidence synthesis (Lobry 1996, Grigoriev 1998, Brock/Wikipedia GC-content, Biopython `GC_skew`/`gc_fraction`, Cuemath population variance) — no contradictions. Follow-ups: none — this was the last pending row in the Extended_GC_Skew_Analysis domain (now fully covered); a future dedicated standalone `gc-skew` unit (still flagged on [[nucleotide-composition-skew]]) would enrich existing concepts rather than add a page.
- 2026-07-13 — ingest docs/algorithms/FileIO/BED_Parsing.md (BED Format Parsing algorithm SPEC — the primary per-algorithm spec for `BedParser`, unit PARSE-BED-001, status *Simplified*; tab-delimited 0-based half-open genomic-interval parser BED3→BED12 plus an interval toolkit). CONTEXT check per brief (there is very likely an existing `bed-format-parsing` concept — REUSE): the concept [[bed-format-parsing]] already exists (created 2026-07-10 from the PARSE-BED-001 Evidence artifact; also enriched by the ANNOT-GFF-001 report) and already synthesizes the coordinate system, BED3→BED12 column ladder, and all validation rules (chromStart≤chromEnd, score clamp 0–1000, strand ∈ {+,−,.}, Auto column-count lock, BED12 block constraints, header/comment skipping). Per the brief (prefer REUSE; create a new page ONLY if genuinely warranted — it was not), ENRICHED the existing concept rather than creating a redundant page. Added the doc's distinct **implementation surface** not previously captured, as two new sections: (1) "Parser surface & behavioral notes (`BedParser`)" — `Parse`/`ParseFile` O(n)/O(1), tab-first-then-whitespace(<3 fields) split, and the interval-toolkit table with the behavioral gotchas: `MergeOverlapping` O(r log r) merges *touching* intervals (`next.ChromStart <= current.ChromEnd`), `Intersect`/`Subtract` O(a·b) worst case, `ExpandIntervals` swaps upstream/downstream on negative strand, `CalculateCoverage` emits depth change-points not per-base rows, and the BED12 block helpers `ExpandBlocks`/`GetTotalBlockLength`/`GetIntrons` with exon expansion `exonStart = chromStart + blockStarts[i]`; (2) "`Auto` vs. explicit format modes" — `Auto` locks field count to the first data line (mixed-width files partially skipped), explicit `Bed3`/`Bed6`/`Bed12` do NOT currently force their nominal field count, display fields (thickStart/thickEnd/itemRgb) parsed but not semantically validated, and bigBed/full-UCSC-toolchain out of scope. Frontmatter: added docs/algorithms/FileIO/BED_Parsing.md to sources (now spec + Evidence + ANNOT-GFF-001 report), source_commit→HEAD (1edb4d1), updated→2026-07-13. Backlog reconciliation: moved the row from [[backlog-pending]] (FileIO 7→6) to [[backlog]] Covered-via-concept → [[bed-format-parsing]]; counts 88→89 covered, 157→156 pending (22 domains unchanged — six FileIO parsers still pending). No new page ⇒ no index.md change (concept already listed). No graph change (no new nodes/edges — concept + its test-unit-registry edge already exist; body [[wikilinks]] mentions auto-derived). Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub, per prior spec-ingest precedent). Spec agrees with the already-recorded Evidence synthesis (UCSC BED FAQ format1, BEDTools Quinlan & Hall 2010) — no contradictions. Follow-ups: the six sibling FileIO parser specs (EMBL/FASTA/FASTQ/GFF/GenBank/VCF) remain pending; each has an existing Evidence-derived concept to enrich by the same pattern.
- 2026-07-13 — ingest docs/algorithms/FileIO/EMBL_Parsing.md (EMBL Parsing algorithm SPEC — the primary per-algorithm spec for `EmblParser`, unit PARSE-EMBL-001, status *Simplified*; line-oriented EMBL/INSDC flat-file parser plus INSDC location parsing, feature/subsequence extraction, remote-aware assembly, and `ToGenBank` conversion). CONTEXT check per brief (there is very likely an existing Evidence-derived concept for EMBL/INSDC flat-file parsing — REUSE): the concept [[insdc-feature-location]] already exists (created 2026-07-10 from the PARSE-GENBANK-001 + PARSE-EMBL-001 Evidence + the GenBank Validation report) and already synthesizes the shared INSDC location-descriptor grammar (1-based-inclusive coords, `n..m`/`n^m`/`<`/`>` partials, `complement`/`join`/`order` + no-nested-join/order rule, remote refs), the operator-assembly semantics (`complement(join(a,b))==join(complement(b),complement(a))`), the offline-first caller-supplied remote resolver, and the worked oracles. Per the brief (prefer REUSE; enrich only genuinely-distinct implementation content; create a new page ONLY if warranted — it was not, the location grammar is deliberately shared with the GenBank cousin and the EMBL line-type/vocab detail already lives in the [[parse-embl-001-evidence]] source page), ENRICHED the existing concept rather than creating a redundant `embl-parsing` page; exact precedent = the 2026-07-13 BED_Parsing spec ingest into [[bed-format-parsing]]. Added the doc's distinct **EMBL-parser implementation surface** not previously captured, as a new "The EMBL parser surface (`EmblParser`)" section: the two-character line-code record shape (`ID`/`AC`/`SV`/`DE`/`KW`/`OS`/`OC`/`RN`/`RA`/`RT`/`RL`/`FT`/`SQ` + `//` terminator, INV-01), the `ID`-line grammar `accession; SV n; topology; mol_type; data_class; division; length BP` with `SV`→`AC` fallbacks, the `\n//` split + `ID`-prefixed-block filter + `GroupLinesByPrefix` concatenation, `AdditionalFields` preservation of non-consumed codes (`DT`/`DR`/`CC`/`OG`), `SQ` letters-only-uppercased extraction, the entry points `ParseLocation`/`GetFeatures`/`GetCDS`/`GetGenes`/`ExtractSequence`/`ResolveLocationSequence`/`ToGenBank` with their O(n)/O(features+references) complexity, the shared `SequenceFormatHelper.ParseLocationParts` location parser, and the *Simplified*-by-design scope (flatten qualifiers to a string dict with bare→"true", no round-trip serialization, skip malformed separators, no occurrence-count/`SQ`-composition validation). Frontmatter: added docs/algorithms/FileIO/EMBL_Parsing.md to sources (now 2 Evidence + GenBank report + this spec), source_commit→HEAD (29ba2bb), updated→2026-07-13. Backlog reconciliation: moved the row from [[backlog-pending]] (FileIO 6→5) to [[backlog]] Covered-via-concept → [[insdc-feature-location]]; counts 89→90 covered, 156→155 pending (22 domains unchanged — five FileIO parsers still pending). Updated the [[insdc-feature-location]] index.md line to note it now also synthesizes the EmblParser algorithm spec. No new page ⇒ concept already listed in index. No graph change (no new nodes/edges — concept + its registry edge already exist; body [[wikilinks]] mentions auto-derived). Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub, per prior spec-ingest precedent). Spec agrees with the already-recorded Evidence synthesis (EBI EMBL User Manual, INSDC Feature Table Definition) — no contradictions. Follow-ups: the four sibling FileIO parser specs (FASTA/FASTQ/GFF/VCF) plus GenBank_Parsing remain pending; each has an existing Evidence-derived concept to enrich by the same pattern.
- 2026-07-13 — ingest docs/algorithms/FileIO/FASTA_Parsing.md (FASTA Parsing algorithm SPEC — the primary per-algorithm spec for `FastaParser` (`Parse`/`ParseFile`/`ParseFileAsync` + `ToFasta`/`WriteFile`), unit PARSE-FASTA-001, status *Simplified*; the `>`-defline + sequence-line sequence-file parser). CONTEXT check per brief (there is very likely an existing Evidence-derived concept for FASTA parsing — REUSE): UNLIKE the BED/EMBL spec-ingest precedents, FASTA had **no** concept page — only the Evidence source page [[parse-fasta-001-evidence]], which had explicitly judged (from the format-facts artifact alone) that "no separate concept page is warranted." The primary spec is genuinely richer than that Evidence: it carries a full **Contract** (4 parameters, return table), **invariants INV-01/02/03**, a **complexity table**, concrete **entry points**, the **opt-in `SequenceAlphabet`/`FastaRecord` overloads** (default strict-DNA `FastaEntry`/`DnaSequence` vs opt-in `StrictDna`/`IupacNucleotide`+gap/`Rna`/`Protein`+`*` returning raw-string-preserving `FastaRecord`), and a deviations table — substantial implementation content not represented anywhere in the wiki, and the backlog explicitly expects a `fasta-parsing` slug. Per the brief ("Create a new concept page ONLY if genuinely warranted and not represented") and the FileIO per-format pattern (BED→[[bed-format-parsing]], EMBL→[[insdc-feature-location]]), CREATED the concept [[fasta-parsing]] (concepts/fasta-parsing.md). It synthesizes: the format (defline first-token=Id / rest=description, multi-line-or-single-line sequence lines with in-sequence whitespace ignored, multi-FASTA = concatenated records); the two-state line-oriented state machine (current header + `StringBuilder` buffer, emit-on-both-header-and-sequence rule at each `>`/EOF); INV-01/02/03; the `FastaParser`/`FastaEntry`/`DnaSequence` contract + entry-point/complexity tables (`Parse` O(n)/O(m), `ToFasta` O(n)/O(n), default width 80); the opt-in `SequenceAlphabet` overloads → `FastaRecord` (no `ToFasta`/`WriteFile` for the non-DNA record type); and the edge cases / intentional simplifications (null-or-empty→∅, header-without-sequence dropped, blank-line + in-sequence-whitespace skip, lowercase→upper, multi-space-defline single-leading-space, round-trip parse→write→parse). Inbound link: updated [[parse-fasta-001-evidence]] (source page) — softened its "no separate concept page is warranted" note to point the implementation surface at [[fasta-parsing]] (Evidence keeps the literature-traced format facts + character-set tables). Frontmatter on the new concept: sources = docs/algorithms/FileIO/FASTA_Parsing.md + docs/Evidence/PARSE-FASTA-001-Evidence.md; source_commit = HEAD (a84ee65); created/updated 2026-07-13. Backlog reconciliation: moved the row from [[backlog-pending]] (FileIO 5→4) to [[backlog]] Covered-via-concept → [[fasta-parsing]]; status line 90→91 covered, 155→154 pending (22 domains unchanged — four FileIO parsers still pending: FASTQ/GFF/GenBank/VCF). index.md: added the [[fasta-parsing]] Concepts line and retargeted the [[parse-fasta-001-evidence]] line's tail to reference the new concept. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub, per prior spec-ingest precedent). New page ⇒ new graph node concept:fasta-parsing + one typed edge relates_to concept:test-unit-registry (mirrors [[bed-format-parsing]]). Spec agrees with the already-recorded Evidence synthesis (Wikipedia FASTA format, NCBI BLAST/FASTA spec, Lipman & Pearson 1985, NC-IUB 1985) — no contradictions. Follow-ups: the four sibling FileIO parser specs (FASTQ/GFF/GenBank/VCF) remain pending; FASTQ/GFF/GenBank each have an existing Evidence-derived concept or sibling to enrich, VCF likewise.
  graph: +1 node, +1 typed edge
- 2026-07-13 — ingest docs/algorithms/FileIO/FASTQ_Parsing.md (FASTQ Parsing algorithm SPEC — the primary per-algorithm spec for `FastqParser`, unit PARSE-FASTQ-001, status *Simplified*; the 4-line-record FASTQ parser + Phred decode/encode, quality/length filtering, quality/adapter trimming, statistics, writing, and paired-end interleave/split helpers). CONTEXT check per brief (searched concepts/ + sources/ for existing fastq / quality-score / phred pages): found the encoding concept [[phred-quality-encoding]], the QC-stats concept [[fastq-quality-statistics]], the trimming concept [[quality-trimming-running-sum]], and the Evidence source page [[parse-fastq-001-evidence]] — but NO concept synthesizing the PARSER surface itself (state machine + `FastqParser` contract), the direct analog of [[fasta-parsing]]. Per the FileIO per-format pattern and the reserved `fastq-parsing` slug in [[backlog-pending]], created a NEW concept [[fastq-parsing]] (sources the algorithm spec) that owns the record state machine, tolerant-assembly rules, full contract surface, complexity, paired-end helpers, and simplifications, while explicitly delegating encoding math → [[phred-quality-encoding]] and QC stats → [[fastq-quality-statistics]] (consumed, not duplicated). Inbound links added from [[fasta-parsing]] (FASTQ sibling now points at the concept) and [[phred-quality-encoding]]. Moved the FileIO/FASTQ_Parsing row from [[backlog-pending]] to the covered table in [[backlog]] (FileIO 4→3; counts adjusted); added the index entry. No source contradictions. Follow-ups: FileIO pending still has GFF/GenBank/VCF parsing (reserved slugs gff-parsing/genbank-parsing/vcf-parsing).
  graph: +1 node, +1 typed edge
- 2026-07-13 — ingest docs/algorithms/FileIO/GFF_Parsing.md (GFF/GTF Parsing algorithm SPEC — the primary per-algorithm spec for the FileIO **`GffParser`**, unit PARSE-GFF-001, status *Simplified*; the fuller of the repo's two GFF3 code paths — GFF3 **and** GTF/GFF2 dialects, hierarchical gene models, `<8`-field rejection, filter/merge/extract/statistics utilities). CONTEXT check per brief (searched concepts/ + sources/ for an existing GFF-parsing concept): found the sibling annotation-layer concept [[gff3-io]] (ANNOT-GFF-001 — `GenomeAnnotator.ParseGff3`/`ToGff3`, a deliberately reduced lightweight helper) and the FileIO Evidence source page [[parse-gff-001-evidence]] (format-facts artifact for PARSE-GFF-001) — but **NO concept synthesizing the FileIO `GffParser` parser surface itself**, the direct analog of [[fasta-parsing]]/[[fastq-parsing]]. The Evidence page had explicitly judged "no dedicated concept page is warranted for this FileIO unit," but the primary spec is genuinely richer: full `GffParser` contract (Parse/ParseFile/reader overloads with the `GffFormat` enum), the line-oriented parse + attribute-dialect rules, the `Auto`-only-`##gff-version` gotcha, `GffRecord`/`GeneModel` shapes, `BuildGeneModels` hierarchy vocabulary + multi-`Parent` fan-out, complexity table, and the merge/extract/statistics/write helpers — substantial implementation content in the wiki nowhere, and the backlog explicitly reserved the `gff-parsing` slug. Per the brief (create a new `gff-parsing` concept following the FileIO per-format pattern, cross-linking — not duplicating — the sibling [[gff3-io]]) and the FASTA/FASTQ precedents, CREATED the concept [[gff-parsing]] (concepts/gff-parsing.md). It owns: the two-code-path disambiguation (FileIO GffParser vs annotation-layer [[gff3-io]]); the format recap (9-col, 1-based inclusive, dialect encodings deferred to Evidence); invariants INV-01..04; the line-oriented parse (skip blank/`#`/`##`, tab-split cols 1–8, `<8`-field reject, `.`→null score/phase, missing col-9→empty attr dict); the per-dialect attribute rules (GFF3 `;`+first-`=`+URL-unescape / GTF `;`+first-space+strip-quotes / GFF2 non-GTF branch); the **`Auto` conservative-detection gotcha** (only `##gff-version`, else GFF3, never infers GTF from syntax); the `GffParser`/`GffRecord`/`GeneModel` contract + entry-point + complexity tables (`Parse` O(n), `MergeOverlapping` O(n log n)); `BuildGeneModels` (`gene`→`mRNA`/`transcript`/`ncRNA`→`exon`/`CDS`/`*utr*` restricted vocabulary, comma-separated multi-`Parent` attaches child to each parent); and the edge cases / intentional simplifications (incl. the `<8`-vs-`<9`-field skip contrast with [[gff3-io]] and the lossy-round-trip writer). Deliberately does NOT re-derive the 9-column schema / dialect / RFC 3986 percent-escape tables / SO EDEN oracle — those stay on [[parse-gff-001-evidence]], cross-linked. Inbound links added: (1) [[parse-gff-001-evidence]] — replaced its "no dedicated concept page is warranted" note with a pointer to [[gff-parsing]] as the implementation-surface synthesis (Evidence keeps the literature-traced format facts); (2) [[gff3-io]] — its two FileIO-path references now cite [[gff-parsing]] alongside [[parse-gff-001-evidence]]. Frontmatter on the new concept: sources = docs/algorithms/FileIO/GFF_Parsing.md + docs/Evidence/PARSE-GFF-001-Evidence.md; source_commit = HEAD (c200712); created/updated 2026-07-13. Backlog reconciliation: moved the row from [[backlog-pending]] (FileIO 3→2, total 154→153) to [[backlog]] Covered-via-concept → [[gff-parsing]]; status line 92→93 covered, 153→152 pending (22 domains unchanged — GenBank/VCF FileIO parsers still pending). index.md: added the [[gff-parsing]] Concepts line and retargeted the [[parse-gff-001-evidence]] line's tail to reference the new concept. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub, per prior spec-ingest precedent). New page ⇒ new graph node concept:gff-parsing + one typed edge relates_to concept:test-unit-registry (mirrors [[fasta-parsing]]/[[fastq-parsing]]). Spec agrees with the already-recorded Evidence synthesis (Sequence Ontology GFF3 v1.26, UCSC GFF/GTF FAQ, Wikipedia GFF, WashU GTF2.2) — no contradictions. Follow-ups: the two sibling FileIO parser specs (GenBank/VCF) remain pending — GenBank folds into [[insdc-feature-location]]/[[parse-genbank-001-evidence]], VCF into [[parse-vcf-001-evidence]].
  graph: +1 node, +1 typed edge
- 2026-07-13 — ingest docs/algorithms/FileIO/GenBank_Parsing.md (GenBank Parsing algorithm SPEC — the primary per-algorithm spec for the FileIO **`GenBankParser`**, unit PARSE-GENBANK-001, status *Simplified*; the NCBI GenBank flat-file parser — keyword-section record model, INSDC feature locations, feature/subsequence extraction, and CDS translation). CONTEXT check per brief (GenBank shares the INSDC feature-location grammar with EMBL; there is very likely an existing concept — REUSE): the concept [[insdc-feature-location]] already exists (created 2026-07-10 from the PARSE-GENBANK-001 + PARSE-EMBL-001 Evidence + the GenBank Validation report; enriched 2026-07-13 with the EMBL parser surface from the EMBL spec) and already synthesizes the shared INSDC location-descriptor grammar (1-based-inclusive coords, `n..m`/`n^m`/`<`/`>` partials, `complement`/`join`/`order` + no-nested-join/order rule, remote refs), the operator-assembly semantics, the offline-first caller-supplied remote resolver, the worked oracles, AND the GenBank two-stage validation verdict paragraph (multi-line qualifier fix, Stage A/B PASS, CLEAN). The GenBank-specific format/testing facts (18 division codes, U49845 canonical record, defensive contracts) already live in the [[parse-genbank-001-evidence]] source page and the validation write-up in [[parse-genbank-001-report]]. Per the brief (REUSE the concept; enrich only genuinely-distinct GenBank-parser implementation content; do NOT create a new concept — it was not warranted; do NOT create a dedicated sources/ page for an algorithm spec) and the exact EMBL-spec-ingest precedent, ENRICHED the existing concept rather than creating a redundant `genbank-parsing` page. Added the doc's distinct **GenBankParser implementation surface** not previously captured, as a new "The GenBank parser surface (`GenBankParser`)" section (parallel to the EMBL one): the keyword-section-in-column-1 record shape (`LOCUS`/`DEFINITION`/`ACCESSION`/`VERSION`/`KEYWORDS`/`SOURCE`+`ORGANISM`/`REFERENCE`/`FEATURES`/`ORIGIN` + `//` terminator, INV-01), the `\n//` split + `LOCUS`-prefixed-block filter + header-keyword section parse with continuation-line merge + `AdditionalFields` preservation, the **fixed-column `LOCUS` line** (name 13–28 / length+bp/aa 30–40 / mol-type 45–47 / topology 56–63 / division 65–67 / date 69–79 DD-MMM-YYYY|DD-MMM-YY) contrasted with EMBL's single delimited `ID` line, the 18-code GenBank division set (PRI…ENV, overlapping-but-not-identical to EMBL's), `KEYWORDS`-`.`→empty / organism-from-`SOURCE` / taxonomy-from-indented-`ORGANISM`, the `ORIGIN` 60-base digit/whitespace-strip + uppercase normalization, the entry points `Parse`/`ParseFile`/`ParseLocation`/`GetFeatures`/`GetCDS`/`GetGenes`/`ExtractSequence`/`GetQualifier`/`TranslateCDS` with O(n) complexity and the `realStart=partStart-1`/`realEnd=partEnd` join+IUPAC-reverse-complement extraction via `FeatureLocationHelper`, the **GenBank-only `TranslateCDS`** decision rule (return existing `/translation` verbatim, else extract CDS nucleotides and translate with a built-in standard codon table, unknown codons → `X`; no alternative-genetic-code metadata — no EMBL analog), the shared `SequenceFormatHelper.ParseLocationParts` location parser, and the *Simplified*-by-design scope (flatten qualifiers to key/value strings, no exact round-trip, no `LOCUS`-length-vs-`ORIGIN` validation, uppercase-normalized case). Frontmatter: added docs/algorithms/FileIO/GenBank_Parsing.md to sources (now 2 Evidence + GenBank Validation report + EMBL spec + this GenBank spec), source_commit→HEAD (0a11c83), updated→2026-07-13. Backlog reconciliation: moved the `docs/algorithms/FileIO/GenBank_Parsing.md` row from [[backlog-pending]] (FileIO 2→1, total 153→152) to [[backlog]] Covered-via-concept → [[insdc-feature-location]]; status line 93→94 covered, 152→151 pending (22 domains unchanged — VCF the last pending FileIO parser). index.md: expanded the [[insdc-feature-location]] line to note it now synthesizes BOTH INSDC-dialect parser specs (GenBankParser + EmblParser). No new page ⇒ no new Concepts index line. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub, per prior spec-ingest precedent). No graph change (no new nodes/edges — concept:insdc-feature-location + its relates_to concept:test-unit-registry edge already exist; body [[wikilinks]] mentions auto-derived). Spec agrees with the already-recorded Evidence + Validation synthesis (NCBI GenBank Sample Record U49845, NCBI GenBank Overview, INSDC Feature Table Definition, Biopython Bio.GenBank.Scanner) — no contradictions. Follow-ups: `docs/algorithms/FileIO/VCF_Parsing.md` (slug `vcf-parsing`) is the one remaining pending FileIO parser — folds into [[parse-vcf-001-evidence]] by the same pattern.
- 2026-07-13 — ingest docs/algorithms/FileIO/VCF_Parsing.md (VCF Parsing algorithm SPEC — the primary per-algorithm spec for the FileIO **`VcfParser`**, unit PARSE-VCF-001, status *Simplified*; the Variant Call Format parser — `##`-header + 8-fixed-column records + optional FORMAT/genotype samples, variant classification, filtering, genotype/sample inspection, statistics + Ti/Tv, INFO helpers, and record writing — the LAST pending FileIO parser). CONTEXT check per brief (searched concepts/ + sources/ for an existing VCF-parsing concept): found the FileIO Evidence source page [[parse-vcf-001-evidence]] (format-facts artifact for PARSE-VCF-001) but **NO concept synthesizing the `VcfParser` parser surface itself** — the direct analog of [[fasta-parsing]]/[[fastq-parsing]]/[[gff-parsing]]. The Evidence page had explicitly judged "no dedicated concept page is warranted," but the primary spec is genuinely richer: the four-entry-point contract (`Parse`/`ParseWithHeader`/`ParseFile`/`ParseFileWithHeader` with distinct null/empty/missing behaviors — `ParseWithHeader(null)` throws, empty→default `VCFv4.3` header), `VcfRecord`/`VcfHeader` shapes, the line-oriented parse + sentinel normalization, the classification table (SNP/MNP/Ins/Del/Symbolic/Complex), the Ti/Tv formula, zygosity rules, the full filter/genotype/statistics/INFO helper surface + complexity table, and the intentional simplifications — substantial implementation content in the wiki nowhere, and the backlog reserved the `vcf-parsing` slug. Per the brief (create a new `vcf-parsing` concept following the FileIO per-format pattern) and the FASTA/FASTQ/GFF precedents, CREATED the concept [[vcf-parsing]] (concepts/vcf-parsing.md). It owns: the format recap (8 fixed columns CHROM/POS(1-based)/ID/REF/ALT/QUAL/FILTER/INFO + optional FORMAT/samples, sample parsing gated on `#CHROM` names); invariants INV-01..03; the four-entry-point table with null/empty/missing semantics; the line-oriented parse (tab-split, `<8`-field + non-integer-`POS` reject/skip-not-raise, header INFO/FORMAT/FILTER promotion + `OtherMetadata`) and sentinel normalization (`QUAL "."`→null, `FILTER "."`→empty array, valueless INFO flag→`"true"`); the `VcfParser`/`VcfRecord`/`VcfHeader` contract + entry-point + complexity tables; `ClassifyVariant` (REF/ALT-length rules with Symbolic override for `*`/`<…>`/breakend `[`/`]`); `CalculateTiTvRatio` (transitions AG/GA/CT/TC ÷ transversions, over ALL ALT alleles, WG≈2.0–2.1/exome≈3.0); zygosity (`GT` split on `/`|`|`, missing allele `.`→indeterminate, `IsHet` false when any allele missing); and the edge cases / intentional simplifications (`FilterPassing`=exactly `PASS` (`.` excluded), stringly-typed no-`Number`/`Type`-enforcement INFO/samples, lossy normalized writer, no BCF/BGZF). Deliberately does NOT re-derive the column/INFO/FORMAT field tables, symbolic/breakend/spanning-deletion grammar, the five 2026-03-11 spec-compliance corrections, or the 1000 Genomes/ClinVar oracle — those stay on [[parse-vcf-001-evidence]], cross-linked. Inbound link: replaced the Evidence page's "no dedicated concept page is warranted" note with a pointer to [[vcf-parsing]] as the implementation-surface synthesis (Evidence keeps the literature-traced format facts); bumped its `updated`→2026-07-13 (kept its `sources:` scoped 1:1 to the Evidence doc + `source_commit` unchanged, per the [[parse-gff-001-evidence]] precedent — the algorithm spec goes only in the concept's `sources:`). Frontmatter on the new concept: sources = docs/algorithms/FileIO/VCF_Parsing.md + docs/Evidence/PARSE-VCF-001-Evidence.md; source_commit = HEAD (a8fa2ab); created/updated 2026-07-13. Backlog reconciliation: moved the row from [[backlog-pending]] (FileIO 1→0 — section removed, domains 22→21, total 152→151) to [[backlog]] Covered-via-concept → [[vcf-parsing]]; status line 94→95 covered, 151→150 pending. index.md: added the [[vcf-parsing]] Concepts line (after [[gff-parsing]]). Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub, per prior spec-ingest precedent). New page ⇒ new graph node concept:vcf-parsing + one typed edge relates_to concept:test-unit-registry (mirrors [[fasta-parsing]]/[[fastq-parsing]]/[[gff-parsing]]). Spec agrees with the already-recorded Evidence synthesis (SAMtools HTS-specs VCFv4.3, Wikipedia VCF, Danecek et al. 2011) — no contradictions. Follow-ups: FileIO parser specs are now ALL covered (BED/EMBL/FASTA/FASTQ/GFF/GenBank/VCF); no FileIO rows remain in [[backlog-pending]].
  graph: +1 node, +1 typed edge
- 2026-07-13 — ingest docs/algorithms/K-mer/K-mer_Counting.md (K-mer Counting algorithm SPEC — the primary per-algorithm spec for the FOUNDATIONAL synchronous count primitive `KmerAnalyzer.CountKmers`, unit KMER-COUNT-001). CONTEXT check per brief (prefer REUSE; check [[k-mer-statistics]] + any k-mer-counting/frequency concept): the whole K-mer counting family already exists — [[asynchronous-kmer-counting]], [[both-strand-kmer-counting]], [[k-mer-statistics]], [[unique-and-mincount-kmers]], [[k-mer-generation]], [[k-mer-positions]], [[k-mer-euclidean-distance]] — BUT every one of them explicitly defers the shared `L−k+1` sync-count definition to a "future sync-count concept (not yet ingested)"; the canonical synchronous `CountKmers` primitive (KMER-COUNT-001) had NO owning concept. backlog-pending reserved the `k-mer-counting` slug. Per the brief (create only if genuinely distinct + unrepresented — it is: it is the primitive the entire family delegates to) CREATED the concept [[k-mer-counting]] (concepts/k-mer-counting.md). It owns: the shared sliding-window count model `Count(w)=Σ 1(S[i..i+k)=w)` over the `L−k+1` overlapping windows; the multi-surface entry-point table (canonical `string`, `DnaSequence` wrapper, span `CountKmersSpan`, cancellation-aware, async, both-strand — identical arithmetic, execution-only differences deferred to the sibling pages); invariants INV-01..04; the contract + alphabet-agnostic literal counting (non-ACGT/IUPAC `N` counted, so the 4^k bound is DNA-only); the **`k≤0` validation-order gotcha** (empty string short-circuits to empty dict BEFORE validating k → no throw; non-empty string k≤0 throws; span path validates k FIRST → span k≤0 always throws) — a genuinely-distinct behavior not on any sibling page; the O(n·k) EFFECTIVE complexity from per-window length-k string-key materialization + hashing (dictionary pressure), suffix-tree evaluated-not-used; worked oracles ATGG k=3→{ATG:1,TGG:1} / GTAGAGCTGT k=2·3·4→9·7/8·8/7·7 / AAAA k=2→{AA:3}; and the deviations (None for core count) / intentionally-simplified (no DNA-alphabet enforcement on raw-string/span, both-strand sums not canonical-collapses) / not-implemented (Jellyfish `-C`/Mash canonical collapsing, 2-bit/minimizer encodings). Deliberately does NOT re-derive the async cancellation contract, both-strand inversion-symmetry, or the AnalyzeKmers statistics — those stay on their sibling pages, cross-linked. Frontmatter: sources = docs/algorithms/K-mer/K-mer_Counting.md (a spec, NOT an Evidence file — no wiki/sources/ page created, per spec-ingest precedent); mcp_tools = count_kmers; source_commit = HEAD (a0600db); created/updated 2026-07-13. mcp-tool-catalog: moved `count_kmers` from KmerAnalyzer Unmapped → mapped [[k-mer-counting]] (Analysis 59→60 mapped, 32→31 unmapped). Backlog reconciliation: moved the row from [[backlog-pending]] (K-mer 3→2, total 151→150) to [[backlog]] Covered-via-concept → [[k-mer-counting]]; status line 95→96 covered, 150→149 pending (21 domains unchanged). index.md: added the [[k-mer-counting]] Concepts line (leading the K-mer counting cluster, before [[asynchronous-kmer-counting]]). Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub, per prior spec-ingest precedent). New page ⇒ new graph node concept:k-mer-counting + one typed edge relates_to concept:test-unit-registry (mirrors the sibling K-mer units). Spec agrees with the family's already-recorded synthesis (Wikipedia K-mer, Rosalind K-mer Composition, Compeau et al. 2011, Marçais & Kingsford 2011) — no contradictions. Follow-ups: two K-mer specs remain pending — K-mer_Frequency_Analysis.md (slug `k-mer-frequency-analysis`) and K-mer_Search.md (slug `k-mer-search`).
  graph: +1 node, +1 typed edge
- 2026-07-13 — ingest docs/algorithms/K-mer/K-mer_Frequency_Analysis.md (K-mer Frequency Analysis algorithm SPEC — the primary per-algorithm spec for the DISTRIBUTION-shaping trio `KmerAnalyzer.GetKmerFrequencies` + `GetKmerSpectrum` + `CalculateKmerEntropy`, unit KMER-FREQ-001). CONTEXT/decision per brief (prefer REUSE existing k-mer concepts; the backlog reserved slug `k-mer-frequency-analysis`; create a new concept ONLY if frequency analysis — relative frequencies / spectrum / distribution — is genuinely distinct from raw counting + statistics and unrepresented): checked [[k-mer-counting]] (owns the shared `L−k+1` count), [[k-mer-statistics]] (owns the scalar `KmerStatistics` summary + the Shannon k-entropy), and [[k-mer-euclidean-distance]] (consumes a frequency vector). Finding: the k-entropy is ALREADY covered (k-mer-statistics owns `H=−Σ p log₂ p` as `AnalyzeKmers.Entropy` and `SequenceComplexity.CalculateKmerEntropy`/SEQ-COMPLEX-KMER-001; this spec's `KmerAnalyzer.CalculateKmerEntropy` is a THIRD entry point to the identical formula), BUT `GetKmerFrequencies` (normalized frequency vector `f_i=c_i/Σc_j`) and `GetKmerSpectrum` (count-of-counts multiplicity histogram) are genuinely distinct, unrepresented OUTPUTS — no concept owns the frequency-map or the spectrum-producer. Per the brief, CREATED the concept [[k-mer-frequency-analysis]] (concepts/k-mer-frequency-analysis.md), cross-linking (NOT duplicating) the siblings. It owns: the three distribution views table (frequencies/spectrum/entropy with definitions), INV-01 Σf_i=1.0 / INV-02 spectrum Σ(count×mult)=L−k+1 / INV-03 0≤H≤log₂D; the contract (all delegate to `CountKmers` → case-insensitive, null/empty→empty/empty/0.0, k>L→empty/empty/0.0, k≤0 on non-empty→ArgumentOutOfRangeException, the shared `k≤0` validation-order deferred to [[k-mer-counting]]); the edge-case table (empty, k>L, single k-mer→{w:1.0}/{1:1}/0, homopolymer AAAA k=2→{AA:1.0}/{3:1}/0); O(n) time / O(u) space for all three on one CountKmers scan; a "the k-entropy is the same Shannon k-entropy owned by [[k-mer-statistics]]" section (third entry point, derivation NOT re-derived); and deviations/limitations — the ACCEPTED deviation that the spec described 4-decimal entropy rounding but the current source returns the raw unrounded double (`Math.Log2`, skips zero-freq terms), plus the observed-only / no-nᵏ-space-normalization / no-smoothing limitation. Deliberately does NOT re-derive the k-entropy formula, bounds, or oracles (those stay on [[k-mer-statistics]]). Relation section cross-links: built on [[k-mer-counting]] (shared count); `GetKmerFrequencies` is the exact frequency vector [[k-mer-euclidean-distance]] consumes for its L2 distance (that concept surfaces the `kmer_frequencies` MCP tool); `GetKmerSpectrum` histogram feeds error-detection / assembly QC → [[kmer-spectrum-error-correction]] (which surfaces the `kmer_spectrum` MCP tool); vs the scalar-summary [[k-mer-statistics]]. MCP catalog LEFT UNTOUCHED — `kmer_frequencies` and `kmer_spectrum` are already mapped to [[k-mer-euclidean-distance]] / [[kmer-spectrum-error-correction]] (not Unmapped, unlike the count_kmers precedent), and the brief did not request MCP remapping; the relationship is documented in prose instead, so no `mcp_tools:` frontmatter added to the new concept. Frontmatter: sources = docs/algorithms/K-mer/K-mer_Frequency_Analysis.md (a spec, NOT an Evidence file — no wiki/sources/ page created, per spec-ingest precedent); source_commit = HEAD (6b60958); created/updated 2026-07-13. Backlog reconciliation: moved the row from [[backlog-pending]] (K-mer 2→1, header 150→149 docs) to [[backlog]] Covered-via-concept → [[k-mer-frequency-analysis]]; status line 96→97 covered, pending 149→148, and the covered-note 151→150 (pre-existing count drift between the two files left otherwise as-is). index.md: added the [[k-mer-frequency-analysis]] Concepts line (before [[k-mer-counting]] in the K-mer cluster). Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub, per prior spec-ingest precedent). New page ⇒ new graph node concept:k-mer-frequency-analysis + one typed edge relates_to concept:test-unit-registry (mirrors the sibling K-mer units). Spec agrees with the family's already-recorded synthesis (Wikipedia K-mer/Entropy, Shannon 1948, Rosalind K-mer Composition, Teeling et al. 2004 TETRA, Chor et al. 2009) — no contradictions. Follow-ups: one K-mer spec remains pending — K-mer_Search.md (slug `k-mer-search`).
  graph: +1 node, +1 typed edge
- 2026-07-13 — ingest docs/algorithms/K-mer/K-mer_Search.md (K-mer Search algorithm SPEC — the SPECIAL-INTEREST-k-mer trio `KmerAnalyzer.FindMostFrequentKmers` + `FindUniqueKmers` + `FindClumps`, unit KMER-FIND-001; the LAST pending K-mer spec). CONTEXT/decision per brief (prefer REUSE existing k-mer concepts; the brief's parenthetical assumed k-mer-search ≈ "locating a query k-mer's positions" and expected a likely fold into [[k-mer-positions]]): on reading the actual spec, K-mer Search is NOT position-lookup — it is the "find k-mers of special interest" family (Rosalind BA1B frequent-words + BA1E (L,t) clumps), genuinely distinct from the where-index [[k-mer-positions]]. Checked all four suggested concepts ([[k-mer-counting]], [[k-mer-positions]], [[k-mer-statistics]], [[k-mer-frequency-analysis]]) plus [[unique-and-mincount-kmers]]. Finding: `FindUniqueKmers` (singletons, Count==1) is ALREADY owned by [[unique-and-mincount-kmers]] (KMER-UNIQUE-001), but `FindMostFrequentKmers` (frequent words, all ties at max count) and `FindClumps` (the (L,t) sliding-window clump finder) are genuinely distinct and UNREPRESENTED — no concept owns them. Per the brief's decision rule (create only if distinct + unrepresented — it is) and the reserved `k-mer-search` slug, CREATED the concept [[k-mer-search]] (concepts/k-mer-search.md). It owns: the three-operation core-model table (most-frequent {w:Count=max}, singletons {w:Count=1}, (L,t) clump); the clump sliding-window count-update algorithm (drop-left/add-right mutable dict + HashSet dedup); INV-01..04; and THE sharp edge = the TWO validation regimes (most-frequent/unique reuse `CountKmers` so k≤0 on non-empty THROWS ArgumentOutOfRangeException, whereas `FindClumps` NEVER throws — null/empty, k≤0, windowSize<k, windowSize>L, minOccurrences≤0 all → empty); complexity O(n) filters / O(n·(L−k+1)) worst clump; oracles ACGTTGCATGTCGCATGATGCATGAGAGCT k=4→{CATG,GCAT} each 3×, TGCA a (25,3) clump in the BA1E sample; deviations None; intentionally-simplified = FindClumps returns only qualifying patterns not their supporting windows (window/multiplicity traces NOT implemented). Deliberately does NOT re-derive the singleton total/distinct/unique terminology or the "unique≠KmerStatistics.UniqueKmers" gotcha — those stay on [[unique-and-mincount-kmers]], cross-linked. Frontmatter: sources = docs/algorithms/K-mer/K-mer_Search.md (a spec, NOT an Evidence file — no wiki/sources/ page created, per spec-ingest precedent); source_commit = HEAD (1a7f803); created/updated 2026-07-13. Inbound links (≥1 required, added as light body cross-refs WITHOUT altering their frontmatter/source_commit, since those pages are not re-synthesized from this spec): [[unique-and-mincount-kmers]] (most-frequent counterpart + note that its FindUniqueKmers is also documented by KMER-FIND-001) and [[k-mer-counting]] (added k-mer-search to its "filters over the shared count" list). Backlog reconciliation: moved the row from [[backlog-pending]] (removed the K-mer section entirely, header 149→148 docs / 21→20 domains) to [[backlog]] Covered-via-concept → [[k-mer-search]]; status line 97→98 covered, pending-tables reference 150→148 docs / 21→20 domains — THE K-MER DOMAIN IS NOW FULLY COVERED. index.md: added the [[k-mer-search]] Concepts line (after [[unique-and-mincount-kmers]] in the K-mer cluster). Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub, per prior spec-ingest precedent). New page ⇒ new graph node concept:k-mer-search + two typed edges (relates_to concept:test-unit-registry mirroring the sibling K-mer units, and relates_to concept:k-mer-counting for the filter-over-the-count relationship). Spec agrees with the family's already-recorded synthesis (Rosalind BA1B, Rosalind BA1E, Wikipedia K-mer) — no contradictions.
  graph: +1 node, +2 typed edges
- 2026-07-13 — ingest docs/algorithms/Metagenomics/PanGenome_Core_Accessory.md (Pan-Genome Construction algorithm SPEC — the primary per-algorithm spec for `PanGenomeAnalyzer.ConstructPanGenome`, unit PANGEN-CORE-001; the LAST pending Metagenomics spec). CONTEXT/decision per brief (REUSE the existing concept): the concept [[pan-genome-core-accessory-partition]] already exists and fully covers the spec's theory (occupancy core/accessory/unique with the fractional `occupancy/N ≥ coreFraction` Roary rule, Kislyuk fluidity φ, Heaps open/closed). No new page warranted. ENRICHED the concept surgically: added the spec (`docs/algorithms/Metagenomics/PanGenome_Core_Accessory.md`) as the FIRST entry in `sources:`, bumped source_commit 9957b47→24019ea and updated 2026-07-10→2026-07-13, and inserted a new "Implementation (`PanGenomeAnalyzer`)" section with the genuinely-distinct impl content NOT previously on the page: entry points (`ConstructPanGenome(genomes, identityThreshold=0.9, coreFraction=0.99)`, private `CalculateGenomeFluidity`/`DeterminePanGenomeType`/`EstimateHeapsDecayExponent`), the sharp `GetCoreGeneClusters` caveat (the standalone Registry `IdentifyCoreGenes` referent uses the FLOOR rule `occupancy ≥ floor(threshold·totalGenomes)`, NOT the fractional test — can admit borderline clusters `ConstructPanGenome` rejects at small N), the k=7 k-mer Jaccard clusterer + suffix-tree-evaluated-and-not-used note, single-dictionary-order Heaps α fit (zero-novelty floored to 1) → order-dependent borderline calls, no fluidity jackknife σ², and complexity (ConstructPanGenome O(g²·s)/O(g); fluidity O(N²·C); Heaps α O(N·C)) + INV-01..07 test coverage. Did NOT re-derive the fluidity formula, occupancy rule, or Heaps criterion (already on the page). NO wiki/sources/ page created (spec, not an Evidence/Validation report — per spec-ingest precedent). Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No new graph nodes/edges (concept reused). Backlog reconciliation: moved the row from [[backlog-pending]] (removed the Metagenomics section entirely, header 148→147 docs / 20→19 domains) to [[backlog]] Covered-via-concept → [[pan-genome-core-accessory-partition]]; status line 98→99 covered, 148→147 pending / 21→20 domains — THE METAGENOMICS DOMAIN IS NOW FULLY COVERED. index.md unchanged (concept already registered, line 485). Spec agrees with the concept's already-recorded synthesis (Tettelin 2005/2008, Kislyuk 2011, Page 2015 Roary, micropan) — no contradictions.
- 2026-07-13 — ingest docs/algorithms/MolTools/DNA_Dimer_Tm.md (DNA Self-/Hetero-Dimer Tm algorithm SPEC — the primary per-algorithm spec for `PrimerDesigner.FindMostStableDimer` + `CalculateDimerMeltingTemperature` / `CalculateSelfDimerMeltingTemperature`, unit PRIMER-TM-001 self-/hetero-dimer extension, status *Simplified*; a Primer3/`ntthal`-style thermodynamic alignment over the SantaLucia & Hicks 2004 unified NN model for the most-stable intermolecular duplex two oligos can form). CONTEXT check per brief (REUSE existing Tm / nearest-neighbor / primer-dimer / thermodynamics concepts): the concept [[primer-dimer-thermodynamics-tm]] already exists (PRIMER-TM-001, created 2026-07-10) and already comprehensively synthesizes this exact algorithm — the SantaLucia unified NN model, the bimolecular Tm Eq. 3 with the x=4/x=1 symmetry factor, the Eq. 5 Na⁺ salt correction, the full `CalculateDimerThermodynamicsNtthal` DP (mismatch/loop/bulge/overhang), and the SAME worked oracles this spec gives (GCGCGCGC self → 40.0906 °C; TGCATGCATG/CATGCATGCA → 25.6596 °C; GCGCATGCGC → 43.1572 °C). The two sibling NN-Tm surfaces [[dna-duplex-nearest-neighbor-thermodynamics]] (SEQ-THERMO-001) and [[melting-temperature]] (SEQ-TM-001 Wallace/Marmur-Doty) also exist and are already cross-linked from it. No new page warranted (dimer Tm is NOT distinct-and-unrepresented — it is the concept's core subject). Per the brief (ENRICH the closest existing concept; add doc to `sources:`; surgical str_replace only for genuinely-distinct dimer content) ENRICHED [[primer-dimer-thermodynamics-tm]] rather than creating a redundant `dna-dimer-tm` page. Added the spec (`docs/algorithms/MolTools/DNA_Dimer_Tm.md`) as the FIRST entry in `sources:`, bumped source_commit 52c02ee→d5b7080 and updated 2026-07-10→2026-07-13. Surgical addition: a new "Two-tier dimer API" paragraph in the dimer section capturing the spec's genuinely-distinct content not previously on the page — the public `FindMostStableDimer` `DimerResult` record (Strand1Start/Strand2Start 0-based 5′ indices, BasePairs, ΔH°/ΔS°/ΔG°37) keeps a SIMPLER gapless O(n·m) contiguous-Watson-Crick-run scorer (reads strand 2 in 3′→5′ order, scores each ≥2 bp run = init + Σ stacks + terminal-A·T penalty + Eq. 5 salt term, keeps the highest-Tm run) and therefore UNDERESTIMATES overhang/loop-stabilised dimers (ASM-02), whereas the Tm methods delegate to the full ntthal DP at machine-precision parity; the dimer defaults (50 mM Na⁺, 50 nM strand — Primer3/ntthal `dna_conc`, differing from the per-oligo Tm's 0.5 µM); the design note that the repository SUFFIX TREE was evaluated and rejected (score-based thermodynamic alignment over offsets, not exact-substring matching); and the dimer path being monovalent-salt only (no divalent Mg²⁺). Did NOT re-derive the NN model, Eq. 3/Eq. 5, the ntthal DP, or the shared oracles (already on the page). NO wiki/sources/ page created (spec, not an Evidence/Validation report — per spec-ingest precedent). Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No new graph nodes/edges (concept reused; no new typed relationships — the spec supports none beyond the already-recorded relates_to concept:test-unit-registry edges). Backlog reconciliation: moved the row from [[backlog-pending]] (MolTools 16→15, header 147→146 docs) to [[backlog]] Covered-via-concept → [[primer-dimer-thermodynamics-tm]]; status line 99→100 covered, 147→146 pending, pending-tables reference 148→147. index.md unchanged (concept already registered). Spec agrees with the concept's already-recorded synthesis (SantaLucia & Hicks 2004, Untergasser 2012 Primer3, primer3-py `thal.c`) — no contradictions. Follow-ups: MolTools has 15 specs still pending (DNA_Hairpin_Folding_Tm → folds into this same concept's hairpin section; NearestNeighbor_Salt_Corrected_Tm / LNA_Adjusted_Nearest_Neighbor_Tm → per-oligo Tm; Primer/Probe/Restriction specs → their own concepts).
- 2026-07-13 — ingest docs/algorithms/MolTools/DNA_Hairpin_Folding_Tm.md (DNA Hairpin Folding + secondary-structure hairpin Tm algorithm SPEC — the primary per-algorithm spec for `PrimerDesigner.FindMostStableHairpin` + `CalculateHairpinMeltingTemperature` + the `HairpinResult` record struct, unit PRIMER-TM-001 hairpin/secondary-structure Tm extension, status *Simplified*; the most-stable intramolecular stem+loop self-fold (MFE) over the SantaLucia & Hicks 2004 unified NN stem stacks + Table 4 hairpin-loop initiation, with the unimolecular concentration-independent Tm Eq. 11). CONTEXT/decision per brief (REUSE the existing concept's hairpin section — a prior partial assessment found it already thoroughly covers this spec): the concept [[primer-dimer-thermodynamics-tm]] already has an "Intramolecular hairpin self-folding" section that comprehensively synthesizes this exact spec — the reused `NnUnifiedParams` stem NN stacks, the Table 4 loop-initiation ΔG°37 by size (3→3.5 … 30→6.3) with loop ΔH°=0 / ΔS°=−ΔG°37·1000/310.15, the Jacobson-Stockmayer 2.44 large-loop extrapolation, the EXCLUDED bimolecular-init + terminal-A·T terms (unimolecular nucleation = loop init), the Eq. 11 concentration-independent Tm (Vallone & Benight 1999), the <3-nt steric floor, the tri/tetraloop bonus (now bundled via [[primer-tm-001-special-loop-evidence]]), and the SAME primary oracle GGGCTTTTGCCC → ΔH°=−25.8, ΔS°=−75.48486216346927, ΔG°37=−2.3883700000000054, Tm=68.6404 °C. No new page warranted (hairpin Tm is NOT distinct-and-unrepresented — it is an existing subsection's core subject; a new `dna-hairpin-folding-tm` page would duplicate it). Per the brief ENRICHED [[primer-dimer-thermodynamics-tm]] rather than creating a redundant page. Added the spec (`docs/algorithms/MolTools/DNA_Hairpin_Folding_Tm.md`) as the 2nd entry in `sources:` (after DNA_Dimer_Tm.md), bumped source_commit d5b7080→1f9b4d8 (updated already 2026-07-13). Surgical addition to the hairpin section: the genuinely-distinct IMPLEMENTATION content the section lacked — the legacy folder signature `FindMostStableHairpin(sequence, minStemLength=2, loopBonusDeltaG37=0)`, the single-stem-one-loop restriction (no bulges/internal/multibranch), the closing-pair-scan-then-maximal-inward-extension MFE mechanic, the loop<3 / stem<minStemLength (minStemLength must be ≥2 = ≥1 NN stack) rejection, the O(n³) worst-case complexity (O(n²) closing-pair scan × O(n) extension; primers ≤~40 nt), `CalculateHairpinMeltingTemperature` returning the Eq. 11 Tm or NaN, and the `HairpinResult` record fields (StemStart/StemEnd 0-based outermost pair, StemLength, LoopSize, ΔH°/ΔS°/ΔG°37). Did NOT re-derive the NN stem/loop model, Eq. 11, the tri/tetraloop bundling, or the shared oracle (already on the page). NO wiki/sources/ page created (spec, not an Evidence/Validation report — per spec-ingest precedent). Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No new graph nodes/edges (concept reused; the spec supports no new typed relationship beyond the already-recorded relates_to concept:test-unit-registry via the hairpin Evidence). Backlog reconciliation: moved the row from [[backlog-pending]] (MolTools 15→14, header 146→145 docs) to [[backlog]] Covered-via-concept → [[primer-dimer-thermodynamics-tm]]; status line 100→101 covered, 146→145 pending, pending-tables reference 147→146. index.md unchanged (concept already registered). Spec agrees with the concept's already-recorded synthesis (SantaLucia & Hicks 2004 Table 1/Table 4/Eqs 7–11, SantaLucia 1998 unified NN, Vallone & Benight 1999) — no contradictions. Follow-ups: MolTools has 14 specs still pending (DNA_Hairpin_Special_Loop_Bonus → the [[primer-tm-001-special-loop-evidence]] tri/tetraloop tables already summarized in this concept's hairpin section; NearestNeighbor_Salt_Corrected_Tm / LNA_Adjusted_Nearest_Neighbor_Tm → per-oligo Tm; Primer/Probe/Restriction specs → their own concepts).
- 2026-07-13 — ingest docs/algorithms/MolTools/DNA_Hairpin_Special_Loop_Bonus.md (DNA Hairpin Special Tri/Tetraloop Bonus algorithm SPEC — the primary per-algorithm spec for the bundled special-loop hairpin path, unit PRIMER-TM-001, status *Production*; full Primer3 ntthal intramolecular-hairpin DP that auto-applies the 16 triloop / 76 tetraloop stability-bonus tables, keyed on the full loop string incl. the closing pair). CONTEXT check: the concept [[primer-dimer-thermodynamics-tm]] already synthesizes this exact unit in its "Intramolecular hairpin self-folding" section (built from the PRIMER-TM-001-SPECIAL-LOOP Evidence, [[primer-tm-001-special-loop-evidence]]) — the bundled 16/76 tables, bsearch keying, ΔH cal/mol / ΔS cal/(K·mol) convention, and the GGGGCGAAAGCCCC/GGGCGAAGCCC oracles are all already present. Per the brief (REUSE; do NOT create a redundant concept for a spec already covered), ENRICHED the existing concept rather than creating a new dna-hairpin-special-loop-bonus page, and did NOT create a wiki/sources/ page (spec, not Evidence). Added the doc's distinct implementation surface not previously captured: the separate public entry point PrimerDesigner.CalculateHairpinThermodynamicsNtthal(sequence, sodiumMolar=0.05) → NtthalHairpin.Run returning HairpinThermodynamics? (ΔH/ΔS/ΔG37/Tm/BasePairs=N/2, null on no-structure/invalid); the full monomer ntthal DP (initMatrix2/fillMatrix2 + internal-loop/bulge modelling + calc_terminal_bp/END5 + tracebacku) at O(n³)/O(n²) vs the legacy single-stem folder; and the key distinction that the ntthal unimolecular Tm CARRIES the [Na⁺] salt correction — Tm = ΔH°/(ΔS° + (N/2−1)·saltCorrection) − 273.15, saltCorrection = 0.368·ln[Na⁺], monovalent-only (dv=dntp=0) — unlike the already-documented legacy Eq. 11 (no salt term). Frontmatter: added docs/algorithms/MolTools/DNA_Hairpin_Special_Loop_Bonus.md to sources, source_commit→HEAD (f8bd7ad), updated 2026-07-13. Backlog reconciliation: moved the row from [[backlog-pending]] (MolTools 14→13) to [[backlog]] Covered-via-concept → [[primer-dimer-thermodynamics-tm]]; counts 101→102 covered, 145→144 pending. No new page/index change, no typed edges added, no contradictions.
- 2026-07-13 — ingest docs/algorithms/MolTools/Hybridization_Probe_Design.md (Hybridization Probe Design algorithm SPEC — the primary per-algorithm spec for the GENERIC `ProbeDesigner.DesignProbes` designer + `DesignTilingProbes` + `CheckSpecificity` + the opt-in `EvaluateTaqManProbe`/`SelectTaqManStrand`, unit PROBE-DESIGN-001, status *Simplified*). CONTEXT/decision per brief (survey wiki/concepts for existing probe/primer/oligo concepts to REUSE; ENRICH the closest existing concept UNLESS probe design is genuinely distinct from primer design and unrepresented → then create a focused concept with ≥1 inbound link). Survey: [[taqman-probe-design-rules]] (PROBE-DESIGN-001, but scoped to the opt-in TaqMan 5'-nuclease hard rules — it treats the generic designer as "the unchanged default" and does NOT synthesize it), [[probe-offtarget-specificity-scan]] (PROBE-VALID-001, the standalone gapped-SW off-target scan), [[primer-dimer-thermodynamics-tm]] / [[primer3-weighted-penalty-objective]] (primer siblings), [[melting-temperature]]. Finding: the GENERIC hybridization-probe designer — candidate-window enumeration + the fixed additive-penalty ranking (baseline 1.0 minus GC/Tm/homopolymer/self-comp/secondary-structure/repeat/terminal-G/C penalties), application-specific defaults (Microarray/FISH/Northern/qPCR/Southern length·Tm·GC windows), tiling probes, prefix-sum GC, the suffix-tree specificity overload, and the Standard/Tiling/Antisense/LNA/MolecularBeacon probe types — is genuinely DISTINCT from primer design and essentially UNREPRESENTED (only mentioned in passing on the TaqMan page). The backlog even reserved the slug `hybridization-probe-design`. Per the brief's decision rule CREATED the concept [[hybridization-probe-design]] (concepts/hybridization-probe-design.md) rather than folding into the TaqMan page. It owns: the primer-vs-probe distinction; the 8-row additive-penalty table (score ranks heuristically, is NOT a hybridization probability); the two-regime Tm (Wallace <14 nt / salt-adjusted `81.5+16.6·log₁₀[Na⁺]+41·GC−600/N` via `CalculateSaltAdjustedTm`, not full NN); the application-defaults table; the algorithm + complexity (O(n×m) scan, prefix-sum GC O(1), CheckSpecificity O(m), DesignTilingProbes O(n)); the "specificity is a post-shortlist filter not a full rerank" sharp edge (maxProbes×5 shortlist; requireUnique=false + specificity 0 → a probe survives with final score 0); INV-01..03; edge cases; and the scope/simplifications (no thermodynamic binding model, no DB-alignment, MGB/LNA/dual-quencher chemistries not implemented — the LNA Tm-adjust is the separate PROBE-DESIGN-001 LNA variant). Deliberately did NOT re-derive the TaqMan hard rules or the off-target scan — those stay on their own pages, cross-linked. Inbound link (≥1 required): added a body [[hybridization-probe-design]] link to [[taqman-probe-design-rules]]'s "Scope" section (which already names the generic designer as the default) — a light cross-ref WITHOUT altering that page's frontmatter/source_commit (not re-synthesized from this spec). Frontmatter: sources = docs/algorithms/MolTools/Hybridization_Probe_Design.md (a spec, NOT an Evidence file — no wiki/sources/ page created, per spec-ingest precedent); source_commit = HEAD (acf2ff3); created/updated 2026-07-13. Backlog reconciliation: moved the row from [[backlog-pending]] (MolTools 13→12, header 144→143 docs) to [[backlog]] Covered-via-concept → [[hybridization-probe-design]]; status line 102→103 covered, 144→143 pending, pending-tables reference 144→143. index.md: added the [[hybridization-probe-design]] Concepts line (before [[taqman-probe-design-rules]] in the PROBE cluster). Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub, per prior spec-ingest precedent). New page ⇒ new graph node concept:hybridization-probe-design + three typed edges (relates_to concept:taqman-probe-design-rules — TaqMan is the opt-in extension of this generic designer per §2.2.1; relates_to concept:probe-offtarget-specificity-scan — the suffix-tree overload's CheckSpecificity + the sibling scan; relates_to concept:melting-temperature — the shared salt-adjusted/Wallace Tm helpers). Spec agrees with the family's already-recorded synthesis (Wikipedia Hybridization probe/FISH/DNA microarray/Nucleic acid thermodynamics, SantaLucia 1998, Breslauer 1986; TaqMan §2.2.1 vs the TaqMan Evidence) — no contradictions. Follow-ups: MolTools has 12 specs still pending (LNA/NN salt-corrected Tm → per-oligo Tm concepts; Off_Target_Analysis, PAM_Site_Detection → CRISPR concept; Primer_Design, Primer_Structure_Analysis, Primer3_Penalty_Objective → primer concepts; Probe_Validation → [[probe-offtarget-specificity-scan]]; the 3 Restriction specs).
  graph: +1 node, +3 typed edges
- 2026-07-13 — ingest docs/algorithms/MolTools/LNA_Adjusted_Nearest_Neighbor_Tm.md (LNA-Adjusted Nearest-Neighbour Tm algorithm SPEC — the primary per-algorithm spec for `PrimerDesigner.CalculateNearestNeighborThermodynamicsLna` + `CalculateMeltingTemperatureNNLna` and the qualitative `ProbeDesigner.EvaluateMgbProbeDesign` check, unit PROBE-DESIGN-001, status *Production*; an opt-in additive extension that adds the McTigue/Peterson/Kahn 2004 sequence-dependent LNA-DNA NN increments ΔΔH°/ΔΔS° to the SantaLucia DNA NN stack so the design Tm reflects the stabilization an internal LNA monomer confers). CONTEXT/decision per brief (survey wiki/concepts for existing Tm/NN/melting-temperature concepts — [[melting-temperature]], [[primer-dimer-thermodynamics-tm]], [[dna-duplex-nearest-neighbor-thermodynamics]], [[hybridization-probe-design]], [[taqman-probe-design-rules]] — and ENRICH the closest existing Tm concept UNLESS the LNA adjustment is genuinely distinct-and-wiki-worthy as its own page). Finding: the LNA adjustment is NOT a new algorithm — it reuses the SAME per-oligo NN engine (`CalculateNearestNeighborThermodynamics`) and the SAME bimolecular Eq. 3 Tm + salt corrections already synthesized in [[primer-dimer-thermodynamics-tm]]'s "Per-oligo design Tm and salt corrections" section, differing only by an ADDITIVE per-step increment for locked bases; the LNA methods live in the same `PrimerDesigner.cs`, and the prior sibling MolTools Tm specs (DNA_Dimer_Tm, DNA_Hairpin_Folding_Tm, DNA_Hairpin_Special_Loop_Bonus) all folded into this exact concept — the concept's own log entries even foreshadowed "LNA_Adjusted_Nearest_Neighbor_Tm → per-oligo Tm". A standalone `lna-adjusted-nearest-neighbor-tm` page would fragment the per-oligo Tm surface. Per the brief ENRICHED [[primer-dimer-thermodynamics-tm]] rather than creating a new page. Added the spec (`docs/algorithms/MolTools/LNA_Adjusted_Nearest_Neighbor_Tm.md`) as the 4th `sources:` entry, bumped source_commit f8bd7ad→59c4c07 (updated already 2026-07-13). Surgical addition: a new "## LNA-adjusted per-oligo Tm (base-modified probes)" section (placed right after the per-oligo section it extends) capturing the genuinely-distinct content — the LNA monomer chemistry (2′-O/4′-C methylene bridge, C3′-endo lock, "largest known increase in thermal stability of any modified DNA duplex", McTigue et al. 2004); the additive ΔH°/ΔS° increment model reusing the base NN engine + Eq. 3 (MELTING `enthalpy += lockedAcidValue`); the 32 verbatim McTigue-2004 LNA-DNA increments (16 5′-locked `X_L N` / 16 3′-locked `M X_L`, keyed by DNA step + locked position, from MELTING 5 `McTigue2004lockedmn.xml`, stored kcal/mol ÷1000); the entry points (`CalculateNearestNeighborThermodynamicsLna`, `CalculateMeltingTemperatureNNLna(…, saltMode=Owczarzy2004Monovalent)`) and the 0-based internal `lnaPositions` in (0,length−1) with set semantics; INV-01 (empty set = perfect-match NN) / INV-02 (internal LNA raises Tm) / INV-03 (terminal or out-of-range → null/NaN — McTigue never parameterised terminal LNA, ASM-01) / INV-04 (increment = verbatim McTigue value); null lnaPositions → ArgumentNullException, empty/<2 nt/non-ACGT → null/NaN; the accepted ~0.09 °C base-model offset vs MELTING `mct04` (base = SantaLucia unified, not McTigue's own reference set, ASM-02); the worked oracle CCATTGCTACC LNA@4 → base ΔH°=−80.8/ΔS°=−221.7 + TTL/AA(+2.326,+8.1) + TLG/AC(−1.540,−3.0) → Tm 63.528 °C (mct04 63.614), all-DNA 59.692 → +3.84 °C; LNA mismatch-discrimination + consecutive/terminal-LNA cooperative models (IDT 2012) out of scope; and the qualitative `ProbeDesigner.EvaluateMgbProbeDesign` 3′-MGB design-rule check (Kutyavin et al. 2000; 12–20 nt window, 3′ attachment, NO quantitative MGB ΔTm — empirical MGB-Eclipse residual; ArgumentNullException on null probe), cross-linked to [[hybridization-probe-design]] and [[taqman-probe-design-rules]]. Did NOT re-derive the base NN model, Eq. 3, or the salt corrections (already on the page). NO wiki/sources/ page created (spec, not an Evidence/Validation report — per spec-ingest precedent). Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No new graph nodes/edges (concept reused; the spec supports no new typed relationship beyond the already-recorded relates_to concept:test-unit-registry). Backlog reconciliation: moved the row from [[backlog-pending]] (MolTools 12→11, header 143→142 docs) to [[backlog]] Covered-via-concept → [[primer-dimer-thermodynamics-tm]]; status line 103→104 covered, 143→142 pending, domains aligned to 19 (matching backlog-pending). index.md unchanged (concept already registered, line 510). Spec agrees with the concept's already-recorded synthesis (SantaLucia 1998 unified NN, McTigue/Peterson/Kahn 2004, MELTING 5 / Dumousseau 2012, Kutyavin 2000) — no contradictions. Follow-ups: MolTools has 11 specs still pending (Melting_Temperature / NearestNeighbor_Salt_Corrected_Tm → the per-oligo Tm surfaces already synthesized here + [[melting-temperature]]; Off_Target_Analysis, PAM_Site_Detection → CRISPR concept; Primer_Design, Primer_Structure_Analysis, Primer3_Penalty_Objective → primer concepts; Probe_Validation → [[probe-offtarget-specificity-scan]]; the 3 Restriction specs).
- 2026-07-13 — ingest docs/algorithms/MolTools/Melting_Temperature.md (Melting Temperature Calculation algorithm SPEC — the primary per-algorithm spec for the MolTools primer-side scalar Tm `PrimerDesigner.CalculateMeltingTemperature` + `CalculateMeltingTemperatureWithSalt` over `ThermoConstants`, unit PRIMER-TM-001, status *Simplified*; the closed-form Wallace (length<14) / Marmur-Doty (length≥14) rule-of-thumb Tm plus an opt-in additive sodium correction). CONTEXT/decision per brief (REUSE the closest existing Tm concept — [[melting-temperature]] and [[primer-dimer-thermodynamics-tm]]). Finding: [[melting-temperature]] (SEQ-TM-001) is already the CANONICAL HOME of the Wallace + Marmur-Doty formulas — the exact core of this spec (identical `WallaceMaxLength=14`, A/T=2, G/C=4, Marmur-Doty 64.9/41.0/16.4 constants, same length-dispatch, same worked oracles GCGCGCGC→32°C and the ACGTACGT…→51.78°C Marmur-Doty case). The formulas themselves need no re-derivation. The backlog-pending row already reserved this doc for slug `melting-temperature`. Per the brief ENRICHED [[melting-temperature]] rather than creating a redundant `primer-melting-temperature` page or folding into [[primer-dimer-thermodynamics-tm]] (which is the NN-thermodynamics concept and explicitly points here as the Wallace/Marmur-Doty home). Added the spec as the 2nd `sources:` entry, bumped source_commit 52c02ee→b506d99, updated 2026-07-10→2026-07-13. Surgical addition: a new "## The MolTools primer-side twin (PRIMER-TM-001) and simple salt correction" section capturing the genuinely-distinct content NOT previously on the page (which had scoped the scalar Tm to SequenceStatistics/SEQ-TM-001 with "no salt correction"): the second primer-oriented implementation `PrimerDesigner.CalculateMeltingTemperature` sharing the same ThermoConstants; the Marmur-Doty `Math.Max(0,…)` ≥0 clamp (INV-04); only A/C/G/T contributing to the counted length (non-DNA ignored not rejected), null/empty/zero-count ⇒ 0; and — the key distinct feature — the OPT-IN simple additive sodium correction `PrimerDesigner.CalculateMeltingTemperatureWithSalt(primer, naConcentration=50 mM)` = `Tm_base + 16.6·log10([Na⁺]/1000)` (Owczarzy 2004 simple form, `ThermoConstants.SaltCoefficient=16.6`, rounded 1 dp), explicitly contrasted with SEQ-TM-001's no-salt scalar and the NN salt models (SantaLucia Eq.5 / Owczarzy 1/Tm) on [[primer-dimer-thermodynamics-tm]]; worked oracle 51.78°C@50 mM → +16.6·log10(0.05)=−21.6 → 30.2°C; the fixed 14-nt switch assumption (spec notes some literature uses ~17–20 bp). Did NOT re-derive the Wallace/Marmur-Doty formulas or their oracles (already on the page). NO wiki/sources/ page created (spec, not an Evidence/Validation report — per spec-ingest precedent). Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No new graph nodes/edges (concept reused; the spec supports no new typed relationship beyond the already-recorded alternative_to concept:primer-dimer-thermodynamics-tm / alternative_to concept:dna-duplex-nearest-neighbor-thermodynamics edges). Backlog reconciliation: moved the row from [[backlog-pending]] (MolTools 11→10) to [[backlog]] Covered-via-concept → [[melting-temperature]]; status line 104→105 covered, 142→141 pending. index.md unchanged (concept already registered). Spec agrees with the concept's already-recorded synthesis (Thein & Wallace 1986, Marmur & Doty 1962, Owczarzy 2004, SantaLucia 1998) — no contradictions; the primer twin's simple 16.6·log10 salt correction is a genuinely-new documented surface, not a conflict. Follow-ups: MolTools has 10 specs still pending (NearestNeighbor_Salt_Corrected_Tm → the NN Tm surfaces on [[primer-dimer-thermodynamics-tm]]; Off_Target_Analysis, PAM_Site_Detection → CRISPR; Primer_Design, Primer_Structure_Analysis, Primer3_Penalty_Objective → primer concepts; Probe_Validation → [[probe-offtarget-specificity-scan]]; the 3 Restriction specs).
- 2026-07-13 — ingest docs/algorithms/MolTools/NearestNeighbor_Salt_Corrected_Tm.md (Nearest-Neighbour Salt-Corrected Melting Temperature algorithm SPEC — the primary per-algorithm spec for the opt-in per-oligo design Tm `PrimerDesigner.CalculateMeltingTemperatureNN` + `CalculateNearestNeighborThermodynamics` and the mismatch/dangling-end `*Mismatch` path, unit PRIMER-TM-001, status *Production*; the SantaLucia 1998/2004 unified NN ΔH°/ΔS° → bimolecular Eq. 3 Tm with SantaLucia Eq. 5 entropy, Owczarzy 2004 monovalent + Owczarzy 2008 divalent Mg²⁺/dNTP salt corrections). CONTEXT/decision per brief (REUSE [[primer-dimer-thermodynamics-tm]] — the NN-thermodynamics home that already synthesizes this exact per-oligo NN engine and salt corrections; verify coverage, enrich only genuinely-distinct content, add the doc to sources). Finding: the concept's "Per-oligo design Tm and salt corrections" section already comprehensively covers this spec — the 10 NN stacks + duplex init (+0.2/−5.7) + terminal-A·T (+2.2/+6.9) + symmetry (−1.4), Eq. 3 Tm with x∈{1,4}, SantaLucia Eq. 5 (N=2·(L−1)), Owczarzy 2004 quadratic 1/Tm + Owczarzy 2008 divalent with R=√[Mg]/[Mon] regime + Ka=3×10⁴ dNTP chelation, the internal single-mismatch (DNA_IMM) + single dangling-end (DNA_DE) tables mirroring Biopython Tm_NN, default C_T=0.5 µM, the SEQ-THERMO-001 (1997 DNA_NN3) distinction, and the 35.8 °C / GCGCGC / ATGCATGC oracles. No new page warranted (per-oligo NN Tm is the existing section's core subject; a `nearestneighbor-salt-corrected-tm` page would duplicate it) — ENRICHED the concept per brief. Added the spec (`docs/algorithms/MolTools/NearestNeighbor_Salt_Corrected_Tm.md`) as the 1st `sources:` entry, bumped source_commit 59c4c07→11bded13 (updated already 2026-07-13). Surgical addition to the per-oligo section: the one genuinely-distinct implementation detail the section lacked — the spec §3.3 three-outcome PARAMETER contract, distinct from the sequence guard: `CalculateMeltingTemperatureNN`'s concentration parameters are domain-validated up front (Tm eval `R·ln(C_T/x)` and corrections `ln[Na⁺]`/`ln[Mg²⁺]` are undefined at non-positive args), so a non-positive strandConcentrationMolar (≤0/NaN), a non-positive sodiumMolar (≤0/NaN incl. zero salt whose ln(0)=−∞ would leak ≈−273.15 °C or a silent NaN), a negative magnesiumMolar, or a negative dntpMolar each throw ArgumentOutOfRangeException — every call resolves to exactly one of: a finite theory-correct Tm, a NaN sentinel (guarded non-computable sequence), or a documented ArgumentOutOfRangeException, never an undisciplined NaN/Inf leak. Did NOT re-derive the NN model, Eq. 3, the salt corrections, the mismatch/dangling-end tables, or the oracles (already on the page). NO wiki/sources/ page created (spec, not an Evidence/Validation report — per spec-ingest precedent). Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No new graph nodes/edges (concept reused; the spec supports no new typed relationship beyond the already-recorded relates_to concept:test-unit-registry via the PRIMER-TM-001-NN Evidence). Backlog reconciliation: moved the row from [[backlog-pending]] (MolTools 10→9) to [[backlog]] Covered-via-concept → [[primer-dimer-thermodynamics-tm]]; recounted actual pending rows = 139 (pre-existing headers 141/142/143 were mutually inconsistent overcounts) and corrected all three to 139 across 18 domains; status line 105→106 covered. index.md unchanged (concept already registered). Spec agrees with the concept's already-recorded synthesis (SantaLucia 1998 PNAS, SantaLucia & Hicks 2004, Owczarzy 2004/2008, Biopython MeltingTemp, Allawi/SantaLucia 1997/98, Peyret 1999, Bommarito 2000) — no contradictions. Follow-ups: MolTools has 9 specs still pending (Off_Target_Analysis, PAM_Site_Detection → CRISPR; Primer_Design, Primer_Structure_Analysis, Primer3_Penalty_Objective → primer concepts; Probe_Validation → [[probe-offtarget-specificity-scan]]; the 3 Restriction specs).
- 2026-07-13 — ingest docs/algorithms/MolTools/Off_Target_Analysis.md (Off-Target Analysis algorithm SPEC — the primary per-algorithm spec for the honest position-weighted heuristic `CrisprDesigner.FindOffTargets` + `CalculateSpecificityScore`, unit CRISPR-OFF-001, status *Simplified*; PAM-constrained near-match enumeration on both strands with a 5-per-seed / 2-per-non-seed positional mismatch penalty and a 0–100 specificity aggregate). CONTEXT/decision per brief (survey wiki/concepts for existing CRISPR/guide-RNA/off-target concepts to REUSE; enrich the closest UNLESS genuinely distinct-and-unrepresented). Finding: [[crispr-guide-rna-design]] already OWNS the CRISPR off-target surface as its "Layer 3 — off-target risk scoring (CRISPR-OFF-001)" and shares the exact test unit (CRISPR-OFF-001) and doc family, but it details only the learned MIT/Hsu-Zhang 2013 and CFD (Doench 2016) models — it referenced the honest heuristic `FindOffTargets`/`CalculateSpecificityScore`/`CalculateOffTargetScore` only in one passing "unchanged; the scored models were added on top" clause. Per the brief ENRICHED that concept rather than creating a redundant `off-target-analysis` page. Added `docs/algorithms/MolTools/Off_Target_Analysis.md` as the 5th `sources:` entry, bumped source_commit 13507add→36296b48, updated 2026-07-10→2026-07-13. Surgical addition inside Layer 3: a "Honest position-weighted heuristic" paragraph capturing the genuinely-distinct content — the `FindOffTargets(guide, genome, maxMismatches=3, systemType=SpCas9)` contract (both-strand PAM scan reusing Layer-0 geometry, yield on `0 < mismatches ≤ maxMismatches`; ArgumentNullException / ArgumentOutOfRangeException on maxMismatches∉[0,5] / ArgumentException on guide-length≠system); the hit record (Position/Sequence/Mismatches/MismatchPositions 0-based/IsForwardStrand/OffTargetScore); the **5-per-seed / 2-per-non-seed** positional penalty with the seed as the **12 bp PAM-proximal window** (last-12 Cas9/SaCas9, first-12 Cas12a) — flagged as a THIRD distinct notion of "seed" on this surface vs Layer-1's 10-nt design seed and the MIT/CFD 20-position vector; the documented subset systems table (SpCas9 NGG/20, SaCas9 NNGRRT/21, Cas12a TTTV/23); `CalculateSpecificityScore` using a **fixed mismatch cap of 4** (independent of the FindOffTargets default 3), summing penalties, returning `max(0, 100 − totalPenalty)`; INV-01 exact matches excluded / INV-02 mismatch bound / INV-03 clamp [0,100] (no off-targets ⇒ 100); and the scope limits (strict single-PAM matching misses alternate-PAM off-targets unless modeled as own system e.g. SpCas9-NAG; mismatches only, no bulges/gaps; no chromatin), with the spec's explicit redirect to CalculateCfdScore / CalculateMit* for published scoring. Did NOT re-derive the MIT/CFD models (already on the page). NO wiki/sources/ page created (spec, not an Evidence/Validation report — per spec-ingest precedent). Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No new graph nodes/edges (concept reused; the spec supports no new typed relationship beyond the already-recorded relates_to concept:test-unit-registry via the CRISPR-OFF-001 report). Backlog reconciliation: moved the row from [[backlog-pending]] (MolTools 9→8, header 139→138 docs) to [[backlog]] Covered-via-concept → [[crispr-guide-rna-design]]; status line 106→107 covered, 139→138 pending, pending-tables reference 139→138. index.md unchanged (concept already registered). Spec agrees with the concept's already-recorded synthesis (Hsu 2013, Fu 2013; the heuristic is explicitly a simplified honest ranking, not a conflict with the learned CFD/MIT models) — no contradictions. Follow-ups: MolTools has 8 specs still pending (PAM_Site_Detection → this CRISPR concept's Layer 0; Primer_Design, Primer_Structure_Analysis, Primer3_Penalty_Objective → primer concepts; Probe_Validation → [[probe-offtarget-specificity-scan]]; the 3 Restriction specs).
- 2026-07-13 — ingest docs/algorithms/MolTools/PAM_Site_Detection.md (PAM Site Detection algorithm SPEC — the primary per-algorithm spec for the CRISPR Layer-0 geometric front end `CrisprDesigner.FindPamSites` (DnaSequence + string overloads) + `GetSystem` + `IupacHelper.MatchesIupac`, unit CRISPR-PAM-001, status *Simplified*; IUPAC-aware both-strand PAM scan resolving each system's motif/orientation/guide-length and extracting the bounds-checked spacer). CONTEXT/decision per brief (the concept [[crispr-guide-rna-design]] already OWNS the CRISPR surface — PAM detection is its "Layer 0 — PAM site detection (CRISPR-PAM-001)" — REUSE it; verify coverage, enrich only genuinely-distinct implementation content, add the doc to sources; do NOT create a dedicated wiki/sources/ page for a spec). Finding: Layer 0 already comprehensively synthesizes this spec — the 7-system PAM/IUPAC/orientation/guide-len table, both-strand scan, IupacHelper.MatchesIupac (NC-IUB 1984), 0-based Position = forward-strand PAM start with reverse Position = len−i−pamLen + reverse-complemented PamSequence, the `targetStart ≥ 0 && targetEnd < len` boundary check, the reverse-strand TargetStart coordinate caveat (CrisprDesigner.cs:1035–1041), and a worked SpCas9 oracle. Per the brief ENRICHED that section rather than creating a redundant `pam-site-detection` page. Added `docs/algorithms/MolTools/PAM_Site_Detection.md` as the 6th `sources:` entry (after Off_Target_Analysis.md), bumped source_commit 36296b48→e2ce3072, updated already 2026-07-13. Surgical addition (one paragraph at the end of Layer 0): the genuinely-distinct content the section lacked — the explicit spacer-interval FORMULAS per orientation (PAM-after-target Cas9: targetStart = PamPos − guideLength, targetEnd = PamPos − 1; PAM-before-target Cas12a/CasX: targetStart = PamPos + pamLength, targetEnd = targetStart + guideLength − 1); the full `PamSite` output record (Position/PamSequence/TargetSequence/TargetStart/IsForwardStrand + resolved `System` metadata: name, PAM pattern, guide length, orientation, description); the O(n) time / O(k) yielded-sites space complexity; the overload contract distinction on INVALID input (DnaSequence throws ArgumentNullException on null vs the string overload returning empty for null/empty and upper-casing before scanning — for valid input both yield identical sites, refining the page's prior "identical results" phrasing); and the spec §5.3 scope note (sequence-pattern detector over the fixed 7-system table only — no cleavage-efficiency / chromatin / unlisted-Cas discovery). Did NOT re-derive the 7-system table, both-strand scan, IUPAC matching, coordinate conventions, or the oracle (already on the page). NO wiki/sources/ page created (spec, not an Evidence/Validation report — per spec-ingest precedent). Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No new graph nodes/edges (concept reused; the spec supports no new typed relationship beyond the already-recorded relates_to concept:test-unit-registry via the CRISPR-PAM-001 report). Backlog reconciliation: moved the row from [[backlog-pending]] (MolTools 8→7, header 138→137 docs) to [[backlog]] Covered-via-concept → [[crispr-guide-rna-design]]; status line 107→108 covered, 138→137 pending, pending-tables reference 138→137. index.md unchanged (concept already registered). Spec agrees with the concept's already-recorded synthesis (Wikipedia PAM/CRISPR, Jinek 2012, Zetsche 2015, Anders 2014) — no contradictions. Follow-ups: MolTools has 7 specs still pending (Primer_Design, Primer_Structure_Analysis, Primer3_Penalty_Objective → primer concepts; Probe_Validation → [[probe-offtarget-specificity-scan]]; the 3 Restriction specs).
- 2026-07-13 — ingest docs/algorithms/MolTools/Primer3_Penalty_Objective.md (Primer3 Weighted Penalty Objective (Per-Primer) algorithm SPEC — the primary per-algorithm spec for `PrimerDesigner.CalculatePrimer3Penalty` + `DefaultPrimer3Weights` / `DefaultPrimer3Optima`, unit PRIMER-TM-001, status *Production*; the de-facto field-standard `p_obj_fn` left/right-primer objective = a weight-gated one-sided-deviation weighted sum over Tm/size/GC%/self_any/self_end/num_ns that Primer3 minimises to select the best primer). CONTEXT/decision per brief (survey wiki/concepts for existing primer-design / Primer3 / primer-selection concepts to REUSE; ENRICH the closest UNLESS genuinely distinct-and-unrepresented). Finding: the concept [[primer3-weighted-penalty-objective]] ALREADY EXISTS (created 2026-07-10 from the PRIMER-TM-001 Evidence) and IS this exact algorithm — it already comprehensively synthesizes the p_obj_fn weighted-sum formula, the `|Tm−60|+|len−20|` default collapse (opt_size=20, opt_tm=60, temp/length weights 1, GC/self/num_ns weights 0, the `DEFAULT_OPT_GC_PERCENT=PR_UNDEFINED_INT_OPT` header-vs-manual-50 subtlety), sign-gated / weight-zero-short-circuit / non-negative / penalty=0-at-optimum structural properties, self_any/self_end as caller-supplied dpal inputs, and the worked oracles ((60,20)→0, (63,20)→3, (57,18)→5, (62.5,22)→4.5, plus non-default-weight cases). No new page warranted (the penalty objective is the concept's core subject; the reserved slug `primer3-penalty-objective` would duplicate it). Per the brief ENRICHED [[primer3-weighted-penalty-objective]]: added the SPEC (`docs/algorithms/MolTools/Primer3_Penalty_Objective.md`) as the FIRST `sources:` entry (before the Evidence file), bumped source_commit 92f89a5→540ba0d, updated 2026-07-10→2026-07-13. Surgical addition: a new "## Scope: which Primer3 terms this reproduces" section capturing the genuinely-distinct spec content the Evidence-derived page lacked — that the implementation reproduces only the SIX core left/right-primer terms as a deterministic O(1) exact weighted sum (not heuristic/probabilistic), with everything else deliberately OUT of scope: not-implemented per-primer terms (`*_TH` thermodynamic-alignment branch + `temp_cutoff`, `pos_penalty`, `end_stability`, `seq_quality`, `repeat_sim`, `template_mispriming`) and the entire pair-level `PRIMER_PAIR_*` objective (Tm-difference, product size, pair complementarity) — all default weight 0 so the DEFAULT selection still matches Primer3 exactly, the gap only appearing when a caller enables an unsupported weight; plus the legacy convenience `Score`/`CalculatePrimerScore` (used by `EvaluatePrimer`/`DesignPrimers`) kept UNCHANGED alongside the new validated `CalculatePrimer3Penalty`, and the Tm-term input coming from the SEQ-THERMO-001 routine ([[melting-temperature]] cross-link). Did NOT re-derive the formula, defaults, structural properties, or oracles (already on the page). NO wiki/sources/ page created (spec, not an Evidence/Validation report — per spec-ingest precedent). Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No new graph nodes/edges (concept reused; the spec supports no new typed relationship beyond the already-recorded relates_to concept:test-unit-registry / relates_to concept:primer-dimer-thermodynamics-tm edges). Backlog reconciliation: moved the row from [[backlog-pending]] (MolTools 7→6, header 137→136 docs) to [[backlog]] Covered-via-concept → [[primer3-weighted-penalty-objective]]; status line 108→109 covered, 137→136 pending. index.md unchanged (concept already registered, line 265). Spec agrees with the concept's already-recorded synthesis (Untergasser 2012, Koressaar 2007, Primer3 source `p_obj_fn`/`pr_set_default_global_args_2`, Primer3 manual §19) — no contradictions. Follow-ups: MolTools has 6 specs still pending (Primer_Design, Primer_Structure_Analysis → primer concepts; Probe_Validation → [[probe-offtarget-specificity-scan]]; the 3 Restriction specs).
- 2026-07-13 — ingest docs/algorithms/MolTools/Primer_Design.md (PCR primer-pair design algorithm SPEC, unit PRIMER-DESIGN-001, `PrimerDesigner.DesignPrimers`). CONTEXT/decision per brief: surveyed the three sibling primer concepts ([[primer3-weighted-penalty-objective]], [[primer-dimer-thermodynamics-tm]], [[hybridization-probe-design]]); the end-to-end primer-PAIR design PIPELINE is genuinely distinct and unrepresented — the penalty page explicitly defers pair ranking as "a separate (future) concern", the dimer page is thermodynamics only, and the probe page is a different reagent. Created a focused new concept [[primer-design]] (backlog already anticipated the `primer-design` slug). Content: the four-stage pipeline (search regions ≤200 bp up/downstream → candidate generation with reverse-complement-before-scoring → greedy independent per-side selection by the legacy additive `CalculatePrimerScore` heuristic → pair compatibility `|Tm_f−Tm_r|≤5 °C AND !HasPrimerDimer`), ProductSize formula (INV-03), the key distinction that DesignPrimers ranks with the legacy 100-based higher-is-better Score and NOT the validated lower-is-better [[primer3-weighted-penalty-objective]] `CalculatePrimer3Penalty`, DefaultParameters (18–25 bp/40–60 % GC/57–63 °C/OptimalTm 60/MaxHomopolymer 4/…), contract+invariants INV-01/02/03, O(n²) complexity, and the simplifications (greedy not global, fixed heuristic, no genome-wide specificity) + the accepted 57–63 vs 55–65 °C Tm deviation. Cross-links all three MolTools siblings + hub [[algorithm-validation-evidence]]. Inbound link added surgically to [[primer3-weighted-penalty-objective]] (its pair-level scope note now points to [[primer-design]] as the pipeline handling that deferred concern). NO wiki/sources/ page created (spec, not an Evidence/Validation report). Registered in index.md Concepts. Backlog: moved MolTools/Primer_Design.md to Covered-via-concept, removed from [[backlog-pending]] (MolTools 6→5), adjusted counts (pending 136→135). Added 2 typed relates_to graph edges (→primer3-weighted-penalty-objective, →primer-dimer-thermodynamics-tm), both explicitly source-supported.
   graph: +1 node, +2 typed edges
- 2026-07-13 — ingest docs/algorithms/MolTools/Primer_Structure_Analysis.md (Primer Structure Analysis algorithm SPEC — the primary per-algorithm spec for the primer secondary-structure QC screens on `PrimerDesigner`, unit PRIMER-STRUCT-001, status *Simplified*: five heuristics — `HasHairpinPotential` boolean stem-loop screen (length-dispatched nested-loop <100 bp / suffix-tree ≥100 bp), `HasPrimerDimer` boolean terminal 3′-complementarity, `Calculate3PrimeStability` terminal-5-mer NN ΔG (SantaLucia 1998 + Primer3 `PRIMER_MAX_END_STABILITY`), `FindLongestHomopolymer`, `FindLongestDinucleotideRepeat`). CONTEXT/decision per brief (survey primer concepts, REUSE closest — [[primer-design]], [[primer-dimer-thermodynamics-tm]], [[primer3-weighted-penalty-objective]] — enrich unless genuinely distinct + unrepresented). Finding: PRIMER-STRUCT-001 is a DISTINCT test unit — the boolean/scalar **screening** surface that deliberately "exposes discrete boolean or scalar quality signals rather than a full thermodynamic folding model", i.e. the low-fidelity ALTERNATIVE to the ntthal thermodynamic-Tm engine (PRIMER-TM-001, [[primer-dimer-thermodynamics-tm]]) and the per-candidate screens [[primer-design]]'s `EvaluatePrimer` consumes. No owning concept existed (primer-design only names the methods as inputs; the thermo concept explicitly is the full-fidelity path and contrasts itself against boolean screens); the 3′-stability step even uses a DISTINCT SantaLucia-1998 ΔG table, not the 2004 NnUnifiedParams. Per brief (create a focused concept when genuinely distinct + unrepresented) CREATED [[primer-structure-qc-screens]] (concepts/primer-structure-qc-screens.md), owning: the five-method table + thresholds (minStemLength=4/minLoopLength=3/minComplementarity=4), the length-dispatched hairpin strategy table + INV-01, the terminal-3′-window dimer rule (last min(8,len1,len2) bases vs revcomp), the SantaLucia-1998 3′-stability ΔG table + terminal-init (+0.98 G·C/+1.03 A·T) + INV-02 + GCGCG/TATAT extremes, the run heuristics INV-03/04, and the scope note (no full thermo folding; the thermo sibling is the separate PRIMER-TM-001 unit). Inbound links wired from [[primer-design]] (EvaluatePrimer's per-candidate screens) and [[primer-dimer-thermodynamics-tm]] (boolean-vs-thermo contrast in the dimer section). sources: the spec; source_commit 208c7e40. Backlog: moved the MolTools row to *Covered via concept*, MolTools pending 5→4, total 135→134, covered 109→110, note added; removed from backlog-pending. index.md: added the concept entry. Graph: +1 node, +2 typed edges (alternative_to concept:primer-dimer-thermodynamics-tm; relates_to concept:primer-design). No contradictions.
   graph: +1 nodes, +2 typed edges
- 2026-07-14 — ingest docs/algorithms/MolTools/Probe_Validation.md (Probe Validation algorithm SPEC, unit PROBE-VALID-001, `ProbeDesigner.ValidateProbe`). CONTEXT/decision per brief: REUSED the existing concept [[probe-offtarget-specificity-scan]] (the gapped Smith–Waterman opt-in scan is the improvement layered on this same unit) rather than create a new page. Enriched it with the genuinely-distinct DEFAULT `ValidateProbe(probe, references, maxMismatches=3, selfComplementarityThreshold=0.3)` surface it lacked: the specificity-score multiplicity penalty (h==0→0, h==1→1, h>1→1/h), pooled default OffTargetHits (on/off separation only in ScanOffTargetsGapped), self-complementarity fraction vs threshold, always-run hairpin secondary-structure flag, IsValid rule (no issues OR OffTargetHits≤1 && SelfComplementarity≤0.4), empty-probe structured-invalid vs null-throws behaviour, and the `CheckSuffix`/`CheckSpecificity(probe, ISuffixTree)` O(m) exact-hit helper. Clarified the default ValidateProbe is ungapped-Hamming and the SW scan is an opt-in supplement (behaviour unchanged). Added the spec to sources:, bumped source_commit→f54e8240. NO wiki/sources/ page created (spec, not an Evidence/Validation report). Backlog: moved to Covered-via-concept → [[probe-offtarget-specificity-scan]] (covered 110→111, pending 134→133); removed from backlog-pending (MolTools 4→3). No new graph nodes/edges (prose-only enrichment). No contradictions.
- 2026-07-14 — ingest docs/algorithms/MolTools/Restriction_Digest_Simulation.md (Restriction Digest Simulation — the primary per-algorithm SPEC for `RestrictionAnalyzer.Digest`/`GetDigestSummary`/`CreateMap`/`AreCompatible`/`FindCompatibleEnzymes`, RESTR family, status *Simplified*, Test Unit ID N/A). CONTEXT/decision: NEW focused concept [[restriction-digest-simulation]] rather than enriching the existing [[restriction-enzyme-filtering]] — digest simulation is a genuinely distinct operation (partitions a **target sequence** into fragments, builds a restriction map, tests end compatibility) vs. filtering's **enzyme-library metadata selection** scope. The filtering concept itself already flagged digest as a separate not-yet-ingested RESTR unit, so the two are siblings, not one page. Considered a shared `restriction-enzyme-analysis` home for all three MolTools restriction specs but rejected: filtering (library selection), digest (sequence→fragments), and site-detection (locate cuts) are three distinct algorithm surfaces, each warranting its own concept; the filtering anchor already exists as its own page. New page synthesizes: forward-strand-cut half-open partition `[0,c1),…,[ck,L)` → k+1 fragments (palindromic sites not double-counted), the fragment-sum invariant (Σ lengths = L, gel-checkable), DigestSummary (descending sizes), RestrictionMap (per-enzyme grouped positions, UniqueCutters/NonCutters, zero-enzymes→full-catalog scan), the AreCompatible blunt/overhang truth table (BamHI+BglII GATC / EcoRV+SmaI blunt / EcoRI+PstI not) + symmetry, no-cut→single-fragment / boundary-null-enzyme / no-zero-length edge cases, and the Simplified-scope limits (no gel migration / partial digest / methylation / circular DNA). Inbound link added from [[restriction-enzyme-filtering]] (its "Relation to the rest of MolTools" section now links the digest sibling). sources: the spec path; source_commit 0f959309. NO wiki/sources/ page created (spec, not an Evidence/Validation report). index.md: added concept entry after [[restriction-enzyme-filtering]]. Backlog: moved to Covered-via-concept → [[restriction-digest-simulation]] (covered 111→112, pending 133→132); removed from backlog-pending (MolTools 3→2). No typed graph edges added (prose-only concept; Test Unit N/A → no test-unit-registry edge; body wikilinks provide mentions edges). No contradictions with existing pages.
- 2026-07-14 — ingest docs/algorithms/MolTools/Restriction_Enzyme_Filtering.md (Restriction Enzyme Filtering algorithm SPEC — the primary per-algorithm spec for the four library-selection helpers `GetEnzymesByCutLength(length)` / `GetEnzymesByCutLength(min,max)` / `GetBluntCutters()` / `GetStickyCutters()` on `RestrictionAnalyzer`, unit RESTR-FILTER-001, status *Production*; pure metadata set-operations over the built-in enzyme table by recognition-length and blunt-vs-sticky end type — no sequence input). CONTEXT check per brief: the concept [[restriction-enzyme-filtering]] already exists (created 2026-07-10 from the RESTR-FILTER-001 Evidence artifact) and already synthesizes the two filter axes, the total blunt/sticky partition (INV-01), inclusive-range/boundary/empty-interval behavior, worked enzyme examples, and the SfiI interrupted-palindrome exclusion — but its `sources:` listed only the Evidence doc (so the spec row was still pending). Per the brief (REUSE the existing concept — do NOT create a new page or a `wiki/sources/` page for the spec), ENRICHED the existing concept rather than creating a redundant page. Added the doc's distinct **implementation surface** not previously captured, as a new "Implementation surface" section: the four total (never-throw, never-null) entry points on `RestrictionAnalyzer` (`Seqeron.Genomics.MolTools`, `RestrictionAnalyzer.cs`); the fixed static `Dictionary<string,RestrictionEnzyme>` library with per-record `RecognitionLength` and record-derived `IsBluntEnd` (`CutPositionForward == CutPositionReverse`), classifying end type from stored cut positions rather than re-deriving cleavage from sequence; lazy LINQ `Where` ⇒ deferred `IEnumerable<>` (insertion order; materialize with `ToList()`); complexity **O(e) time / O(1) extra space** single linear pass with the "suffix tree evaluated and inapplicable — no text to search" reuse note; and the **Not-implemented** overhang-direction (5'/3') / overhang-sequence filtering (use `AreCompatible` / `FindCompatibleEnzymes` on [[restriction-digest-simulation]]). Frontmatter: added docs/algorithms/MolTools/Restriction_Enzyme_Filtering.md to sources (now spec + Evidence), source_commit→HEAD (6a76515), updated→2026-07-14. Backlog reconciliation: moved the row from [[backlog-pending]] (MolTools 2→1) to [[backlog]] Covered-via-concept → [[restriction-enzyme-filtering]]; pending pointer count 132→131 (18 domains unchanged — Restriction_Site_Detection still pending in MolTools). No new page ⇒ no index.md change (concept already listed). No graph change (no new nodes/edges — concept + its test-unit-registry edge already exist; spec adds no new typed edges). Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub, per prior spec-ingest precedent). Spec agrees with the already-recorded Evidence synthesis (Wikipedia "Restriction enzyme" / "Sticky and blunt ends" / "List of restriction enzyme cutting sites", NEB/REBASE KpnI/EcoRI, PMC SfiI interrupted palindrome) — no contradictions. Follow-ups: the sibling `docs/algorithms/MolTools/Restriction_Site_Detection.md` (expected slug `restriction-site-detection`) remains the one pending row in the MolTools domain.

## [2026-07-14] ingest | docs/algorithms/MolTools/Restriction_Site_Detection.md → restriction-site-detection (concept)
   Created the RESTR family's **location** primitive concept — given a sequence + enzyme, find
   where the recognition sequence occurs and where it cuts. Distinct from its two siblings:
   [[restriction-enzyme-filtering]] (enzyme selection, no sequence) and
   [[restriction-digest-simulation]] (fragmentation, consumes the forward-strand cut positions
   detection produces). Both-strand IUPAC scan, cut-position formulas, palindrome double-report,
   overhang classification, RestrictionSite contract + invariants. Added inbound wikilinks from
   both sibling concepts (surgical). Moved the backlog row to *Covered via concept*; removed the
   MolTools section from [[backlog-pending]] — this CLOSES the MolTools domain (113 covered / 131
   pending across 17 domains). No typed graph edges (spec cites no cross-unit relationship beyond
   the sibling mentions, which auto-derive). No contradictions.

## [2026-07-14] ingest | docs/algorithms/Oncology/Allele_Specific_Copy_Number_Derivation.md → allele-specific-copy-number-ascat (concept, enriched)
   FIRST Oncology *algorithm-spec* reconciled (the domain's concept pages were all synthesized
   from Evidence docs, so every Oncology spec row was still pending). CONTEXT check: the concept
   [[allele-specific-copy-number-ascat]] (created 2026-07-09 from the ONCO-ASCAT-001 Evidence
   artifact) already synthesizes this exact unit — ASCAT nA/nB inversion + grid GoF, ASPCF
   penalised-LS segmentation, Battenberg two-state subclonal, and McGranahan/PICTograph/DeCiFering
   multiplicity+CCF. Per brief, ENRICHED it rather than create a redundant page. Added the primary
   spec's distinct surface not previously captured: (1) the reported minor-allele GoF vs the
   *selection* objective (adds the major-allele integer distance and breaks exact ties toward the
   lower ploidy ψ — the 2n-vs-4n parsimony convention); (2) the `OncologyAnalyzer` implementation
   surface — the dual segmentation path (greedy `SegmentAlleleSpecific` O(L) retained as an accepted
   deviation alongside global-optimum `SegmentAlleleSpecificAspcf` O(L²)/chrom), plus
   `FitPurityPloidy`/`DeriveMultiplicity`/`FitSubclonalCopyNumber`, and the *Simplified* status
   (fixed-grid fit, two-adjacent-state subclonal, no asmultipcf/WGD refit). Frontmatter: added the
   spec to sources (now spec + Evidence), source_commit→77e902ee, updated→2026-07-14. NO wiki/sources/
   page created (spec, not an Evidence/Validation report). No new page ⇒ no index.md change (concept
   already listed). Backlog: added Covered-via-concept row → [[allele-specific-copy-number-ascat]]
   (covered 113→114, pending 131→130); removed from [[backlog-pending]] (Oncology 37→36, total
   130→129). No new graph nodes/edges (prose-only enrichment; the concept's test-unit-registry edge
   already exists; body wikilinks provide mentions edges) ⇒ graph lint/extract skipped. Hub
   [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub, per prior
   spec-ingest precedent). No contradictions: the spec and the already-recorded Evidence synthesis
   agree stage-for-stage (Van Loo 2010 + ascat.runAscat.R, Nilsen 2012 + Ross 2021, Nik-Zainal 2012 /
   Battenberg, McGranahan 2016 / PICTograph / DeCiFering — each a disjoint stage).

## [2026-07-14] ingest | docs/algorithms/Oncology/Cancer_Cell_Fraction_Estimation.md → cancer-cell-fraction-clonal-clustering (concept, enriched)
   Second Oncology *algorithm-spec* reconciled. CONTEXT check: the concept
   [[cancer-cell-fraction-clonal-clustering]] (created 2026-07-09 from the ONCO-CCF-001 Evidence
   artifact) already synthesizes this exact unit — the McGranahan/Zheng/Tarabichi CCF closed form
   `CCF = f·(ρ·N_T + 2(1−ρ))/(ρ·m)`, the [0,1] reported cap with exposed RawCcf, and the
   deterministic Lloyd 1D k-means with quantile seeding + highest-centroid clonal rule. Per brief,
   ENRICHED it rather than create a redundant page. Added the primary spec's distinct implementation
   surface not previously captured: the two `OncologyAnalyzer` entry points with return shapes
   (`CcfEstimate.Ccf`/`RawCcf`; `CcfClustering.Centroids`/`Assignments` in input order/
   `ClonalClusterIndex`=k−1), their complexities (EstimateCcf O(1); ClusterCcfValues O(n·k·i) +
   O(n log n) sort, O(n+k) space), the exact exception types on each precondition, the quantile
   seed (j+0.5)/k and lower-index tie-break, and the suffix-tree-not-applicable note. Frontmatter:
   added the spec to sources (now spec + Evidence), source_commit→a78f8c60, updated→2026-07-14. NO
   wiki/sources/ page created (spec, not an Evidence/Validation report). No new page ⇒ no index.md
   change (concept already listed). Backlog: added Covered-via-concept row →
   [[cancer-cell-fraction-clonal-clustering]] (covered 114→115, pending 130→129); removed from
   [[backlog-pending]] (Oncology 36→35, total 129→128). No new graph nodes/edges (prose-only
   enrichment; the concept's existing typed edges are unchanged; body wikilinks provide mentions
   edges) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs
   are not added to the hub, per prior spec-ingest precedent). No contradictions: the spec and the
   already-recorded Evidence synthesis agree stage-for-stage.

## [2026-07-14] ingest | docs/algorithms/Oncology/Cancer_Variant_Annotation.md → [[cancer-variant-tier-classification-amp-asco-cap]] (concept enrich)
   Oncology reconciliation: the AMP/ASCO/CAP four-tier concept already existed (synthesized from
   ONCO-ANNOT-001-Evidence). Added the primary spec to its sources: and bumped source_commit to HEAD.
   Surgical enrich: new "Implementation surface (ONCO-ANNOT-001 spec)" section — CancerVariantAnnotationInput
   record fields, the three entry-point signatures with per-op complexity (O(1) classify / O(n) annotate /
   O(1) COSMIC lookup), the BenignPopulationMafThreshold=0.01 constant, ordinal COSMIC key equality +
   MatchCancerHotspots caller-supplied-set reuse, and the suffix-tree-not-applicable note. No wiki/sources/
   page (spec ≠ Evidence/Validation report). No new page ⇒ index.md unchanged (concept already listed).
   Backlog: added Covered-via-concept row (covered 115→116, pending 129→128); removed from
   [[backlog-pending]] (Oncology 35→34). No new graph nodes/edges (implementation prose only; the concept's
   existing typed edges are unchanged; body wikilinks provide mentions edges) ⇒ graph lint/extract skipped.
   Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No contradictions:
   the spec and the already-recorded Evidence synthesis agree on tiers, the 1% cutoff, and III↔IV discrimination.

## [2026-07-14] ingest | docs/algorithms/Oncology/Clinical_Actionability_Assessment.md → [[clinical-actionability-oncokb-levels]] (concept enrich)
   Oncology reconciliation: the OncoKB therapeutic-levels concept already existed (synthesized from
   ONCO-ACTION-001-Evidence). Added the primary spec to its sources: and bumped source_commit to HEAD.
   Surgical enrich: new "Implementation surface (ONCO-ACTION-001 spec)" section — the five OncologyAnalyzer
   entry points (AssessActionability O(n·k) / ClassifyActionabilityLevel O(k) / GetTherapyRecommendations
   O(k log k) / CompareLevels O(1) / IsStandardCare O(1)), the OncoKbLevel enum integer-order encoding of
   the combined R1>1>2>3A>3B>4>R2 order (CompareLevels = integer compare, no lookup table), the
   VariantActionabilityInput null-vs-empty associations contract, the four output fields
   (Highest{Sensitive,Resistance,Combined}Level + IsActionable), and the suffix-tree-not-applicable /
   Framework-status notes. No wiki/sources/ page (spec ≠ Evidence/Validation report). No new page ⇒
   index.md unchanged (concept already listed). Backlog: added Covered-via-concept row (covered 116→117,
   pending 128→127); removed from [[backlog-pending]] (Oncology 34→33, total 127→126). No new graph
   nodes/edges (implementation prose only; the concept's existing typed edges are unchanged; body wikilinks
   provide mentions edges) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged
   (spec docs are not added to the hub). No contradictions: the spec and the already-recorded Evidence
   synthesis agree on the seven levels, the three ordering axes, and the NotActionable empty-level semantics.

## [2026-07-14] ingest | docs/algorithms/Oncology/Clonal_Hematopoiesis_Filtering.md → [[clonal-hematopoiesis-cfdna-filtering]] (concept enrich)
   Oncology reconciliation: the CHIP-filtering concept already existed (ONCO-CHIP-001, synthesized from
   ONCO-CHIP-001-Evidence) and already covers the full primary spec — the three methods (IdentifyCHIPVariants,
   FilterCHIP with matched-WBC rule a + gene+VAF fallback rule b, CallVariantOrigin strict Bolton rule), the
   VAF≥0.02 threshold, the {DNMT3A,TET2,ASXL1,TP53,JAK2,SF3B1,SRSF2,PPM1D} panel, the φ=2.0/1.5 fold ratio,
   and matching worked oracles. Added the primary spec (docs/algorithms/Oncology/Clonal_Hematopoiesis_Filtering.md)
   as the FIRST entry in sources: (ahead of the Evidence doc), bumped source_commit 90f75a1→f0c2c13 and
   updated 2026-07-09→2026-07-14. Surgical enrich: new "Entry points and implementation" block capturing the
   spec's genuinely-distinct implementation content not previously on the page — all four static
   OncologyAnalyzer methods incl. the fourth IsCanonicalChipGene predicate, ChipVariant/WbcObservation record
   shapes, exact (chrom,1-based pos,ref,alt) locus key, the HashSet WBC-loci membership + Dictionary
   locus→best-observation data structures, O(n·g) / O(n+w) linear complexity, the suffix-tree-not-used note,
   and input-order preservation. No wiki/sources/ page (spec ≠ Evidence/Validation report). No new page ⇒
   index.md unchanged (concept already listed). Backlog: added Covered-via-concept row (covered 117→118,
   pending 127→126); removed from [[backlog-pending]] (Oncology 33→32, total 126→125). No new graph
   nodes/edges (implementation prose only; the concept's existing typed edges are unchanged; body wikilinks
   provide mentions edges) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged
   (spec docs are not added to the hub). No contradictions: the spec and the already-recorded Evidence
   synthesis agree on the three methods, all thresholds, the driver-gene panel, and the origin-call semantics.

## [2026-07-14] ingest | docs/algorithms/Oncology/Clonal_Subclonal_Classification.md → [[clonal-subclonal-classification-ccf-posterior]] (concept enrich)
   Oncology reconciliation: the CCF-posterior clonal/subclonal concept already existed (ONCO-CLONAL-001,
   synthesized from ONCO-CLONAL-001-Evidence) and already covers the full primary spec — the expected-VAF
   relation f(c)=αMc/(2(1−α)+αq) (Landau M=1 + Satas/DeCiFering multiplicity-general form), the uniform-prior
   Binomial posterior on the 100-point grid c∈[0.01,1], the verbatim classification rule clonal iff
   P(CCF>0.95)>0.5, the strict CCF>0.95 IdentifyClonalMutations threshold, the per-variant q-over-ploidy-scalar
   assumption, and matching worked oracles (A1/B2/C1/D/E). Added the primary spec
   (docs/algorithms/Oncology/Clonal_Subclonal_Classification.md) as the FIRST entry in sources: (ahead of the
   Evidence doc), bumped source_commit 7309394→9a7b5ef and updated 2026-07-09→2026-07-14. Surgical enrich: new
   "Implementation (ONCO-CLONAL-001 spec)" block capturing the spec's genuinely-distinct implementation content
   not previously on the page — OncologyAnalyzer.cs entry points ClassifyClonality(variants,purity) /
   IdentifyClonalMutations(ccfValues), the log-space Binomial with C(N,a) omitted (cancels under
   normalisation), the degenerate all-zero-posterior → flat-posterior fallback (stays subclonal), the
   suffix-tree-not-used note, and O(n·G)/O(n) + O(m) complexity. No wiki/sources/ page (spec ≠
   Evidence/Validation report). No new page ⇒ index.md unchanged (concept already listed). Backlog: added
   Covered-via-concept row (covered 118→119, pending 126→125); removed from [[backlog-pending]] (Oncology
   32→31, total 125→124). No new graph nodes/edges (implementation prose only; the concept's existing typed
   edges relates_to test-unit-registry / alternative_to cancer-cell-fraction-clonal-clustering / depends_on
   allele-specific-copy-number-ascat are unchanged; body wikilinks provide mentions edges) ⇒ graph
   lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub).
   No contradictions: the spec and the already-recorded Evidence synthesis agree on the model, both thresholds
   (0.95 CCF / 0.5 posterior), the grid, and the point-estimate variant.

## [2026-07-14] ingest | docs/algorithms/Oncology/Complex_Rearrangement_Classification.md → enriched [[chromothripsis-inference]]
   Oncology reconciliation (spec ≠ Evidence, so NO wiki/sources/ page). The concept synthesized earlier from
   ONCO-SV-001-Evidence already covered the spec's model (six hallmark criteria, oscillation counting +
   two-state hallmark, ≥10 first-pass screen, ≥6 SV floor, ≥7/4–6 confidence tiers, exponential-null CV>1
   clustering). Added the spec to sources: and bumped source_commit → ff6a38af. Surgical str_replace added a
   genuinely-distinct "Implementation (OncologyAnalyzer)" section: three entry points
   (CountCopyNumberStateOscillations / TestBreakpointClustering / ClassifyComplexRearrangement), complexity
   (O(n), O(m log m)), output-record fields, the six named decision constants, the k-transitions→k+1-segments
   rule, and the criteria-A/B-only gate (C–F / chromoplexy / BFB out of scope; suffix tree N/A). No new page ⇒
   index.md unchanged (concept already listed). Backlog: added Covered-via-concept row (covered 119→120,
   pending 125→124); removed from [[backlog-pending]] (Oncology 31→30, total 124→123). No new graph
   nodes/edges (implementation prose only; concept's existing typed edges relates_to test-unit-registry /
   copy-number-alteration-classification unchanged; body wikilinks provide mentions edges) ⇒ graph
   lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub).
   No contradictions: spec and the recorded Evidence synthesis agree on the model, thresholds, and worked datasets.
- 2026-07-14 — ingest docs/algorithms/Oncology/Copy_Number_Alteration_Classification.md (Copy-Number Alteration Classification SPEC — the primary per-algorithm spec for `OncologyAnalyzer.Log2RatioToCopyNumber`/`CallCopyNumber`/`ClassifyCopyNumber`/`ClassifyCopyNumbers`, unit ONCO-CNA-001, status Production; log2 copy ratio → absolute CN → five discrete CNA states via CNVkit hard thresholds). CONTEXT/decision per Oncology reconciliation brief: REUSED the existing concept [[copy-number-alteration-classification]] (created 2026-07-09 from the ONCO-CNA-001 Evidence artifact) rather than create a new page or a wiki/sources/ page (spec ≠ Evidence report). The concept already covered the model (n=2·2^log2, absolute_threshold binning, default cutoffs, GISTIC2 corroboration, NaN/boundary/purity corner cases, diploid-reference assumption, docs-vs-source threshold note); its sources: listed only the Evidence doc, keeping the spec row PENDING. Enriched with a genuinely-distinct "Implementation surface (ONCO-CNA-001)" section the Evidence page lacked: the four OncologyAnalyzer entry points (continuous / integer / full CopyNumberCall / order-and-length-preserving batch), the CopyNumberCall record fields, the validation contract (thresholds exactly-four strictly-ascending non-NaN → ArgumentException; ploidy>0 → ArgumentOutOfRangeException; null batch → ArgumentNullException; null thresholds → default fallback), and O(1)/O(m) complexity plus the single-region (no segmentation / no allele-specific) scope. Added the spec to sources: (now spec + Evidence), bumped source_commit→1555a132, updated→2026-07-14. Backlog: moved the row to [[backlog]] Covered-via-concept → [[copy-number-alteration-classification]] (covered 120→121, pending 124→123); removed from [[backlog-pending]] (Oncology 30→29, pending-total 123→122). No new page ⇒ index.md unchanged (concept already listed). No new graph nodes/edges (implementation prose only; concept's existing relates_to test-unit-registry edge and body wikilinks unchanged) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No contradictions: spec and the recorded Evidence synthesis agree on the model, thresholds, and state mapping.
## [2026-07-14] ingest | docs/algorithms/Oncology/Driver_Mutation_Detection.md → enriched [[driver-gene-classification-20-20-rule]]
   Oncology reconciliation (spec ≠ Evidence, so NO wiki/sources/ page). The concept synthesized earlier from
   ONCO-DRIVER-001-Evidence already covered the 20/20-rule model (Vogelstein 2013 / Tokheim 2020; f_OG > 0.20
   recurrent-missense → Oncogene, f_TSG > 0.20 truncating → TSG, strict `>`, ≥2 recurrence per Miller 2017,
   dual-pass dominant-fraction tie-break → Ambiguous, IDH1/dispersed-truncating oracles, companion ops
   MatchCancerHotspots / ScoreDriverPotential=max(f) / IdentifyDriverMutations⊆input). Added the spec to
   sources: (now spec + Evidence), bumped source_commit → a9f32c33, updated → 2026-07-14. Surgical str_replace
   added a genuinely-distinct "Implementation (per the algorithm spec)" section the Evidence-derived concept
   lacked: OncologyAnalyzer.cs location + four entry points, input-order preservation, 1-based per-position
   recurrence dictionary, ordinal (case-sensitive) gene match, O(1) (gene,position) hotspot set lookup with the
   suffix-tree-not-applicable note, ArgumentNullException on null inputs, and O(N) complexity (space O(P) /
   O(G+H)). No new page ⇒ index.md unchanged (concept already listed). Backlog: added Covered-via-concept row
   (covered 122→123, pending 122→121); removed from [[backlog-pending]] (Oncology 28→27, total 121→120). No new
   graph nodes/edges (implementation prose only; concept's existing relates_to test-unit-registry edge and body
   wikilinks unchanged) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs
   are not added to the hub). No contradictions: spec and the recorded Evidence synthesis agree on model,
   thresholds, and worked datasets.

- 2026-07-14 — ingest docs/algorithms/Oncology/CtDNA_Analysis.md (ctDNA Analysis SPEC — the primary per-algorithm spec for `OncologyAnalyzer.CtDnaDetectionProbability`/`ExpectedMutantMolecules`/`IsCtDnaDetected`/`CalculateTumorFraction`/`CalculateMeanVaf`/`HaploidGenomeEquivalents`, unit ONCO-CTDNA-001, status Production; Poisson limit-of-detection, tumour fraction = 2·VAF, mean VAF, mass→genome-equivalents). CONTEXT/decision per Oncology reconciliation brief: REUSED the existing concept [[ctdna-detection-and-tumor-fraction]] (created 2026-07-09 from the ONCO-CTDNA-001 Evidence artifact) rather than create a new page or a wiki/sources/ page (spec ≠ Evidence report). The concept already covered the whole model (Poisson p=1−e^(−ndk) + λ<3 low-burden regime, the τ≥0.95 ∧ λ≥1 detectability decision, TF=2·VAF diploid identity, mean VAF, 3.3pg⇒303 GE/ng conversion, worked oracles, detection range/background floor); its sources: listed only the Evidence doc, keeping the spec row PENDING. Enriched with genuinely-distinct implementation content the Evidence page lacked: a "Public surface (OncologyAnalyzer)" note (the six named entry points incl. ExpectedMutantMolecules returning λ directly, O(1) scalar / O(n) reporter-aggregate complexity, shared private CalculateVaf helper, direct 1−e^(−λ) with no Math.Expm1) and a "Not implemented" scope boundary (fragmentomics AnalyzeFragmentSizeDistribution absent — no BAM infra; CHIP-background filtering delegated to [[clonal-hematopoiesis-cfdna-filtering]] ONCO-CHIP-001, multi-variant MRD calling to [[tumor-informed-mrd-detection]] ONCO-MRD-001). Added the spec to sources: (now spec + Evidence), bumped source_commit→0ae8dfa5, updated→2026-07-14. Backlog: moved the row to [[backlog]] Covered-via-concept → [[ctdna-detection-and-tumor-fraction]] (covered 121→122, pending 123→122); removed from [[backlog-pending]] (Oncology 29→28, pending-total 122→121). No new page ⇒ index.md unchanged (concept already listed). No new graph nodes/edges (implementation prose only; concept's existing relates_to test-unit-registry / clonal-hematopoiesis-cfdna-filtering edges and body wikilinks unchanged) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No contradictions: spec and the recorded Evidence synthesis agree on the model, thresholds, constants, and worked datasets.

- 2026-07-14 — ingest docs/algorithms/Oncology/Focal_Amplification_Detection.md (Focal Amplification Detection SPEC — primary per-algorithm spec for `OncologyAnalyzer.DetectFocalAmplifications`/`IdentifyAmplifiedOncogenes`/`IsFocalAmplification`, unit ONCO-CNA-002, status Simplified; GISTIC2 length-based focal/broad split at broad_len_cutoff=0.98 + amplitude gate t_amp=0.1, arm→oncogene panel mapping). CONTEXT/decision per Oncology reconciliation brief: REUSED the existing concept [[focal-amplification-detection]] (created 2026-07-09 from the ONCO-CNA-002 Evidence artifact) rather than create a new page or a wiki/sources/ page (spec != Evidence report). The concept already covered the full model (two-part predicate, strict <0.98 boundary rule, six-gene oncogene panel with 12q→MDM2+CDK4, worked dataset, corner cases, both assumptions). Enriched with genuinely-distinct implementation content the Evidence page lacked: an "Implementation surface (ONCO-CNA-002)" section (the three OncologyAnalyzer entry points + thresholds default, CopyNumberArmSegment/FocalAmplificationThresholds types, End−Start length, Ordinal-ignore-case arm matching, ArgumentNull/Argument validation, O(n)/O(n+g) complexity with fixed panel g=6, suffix tree N/A, segmentation upstream in StructuralVariantAnalyzer.SegmentCopyNumber SV-CNV-001/ONCO-CNA-001, INV-01..04). Added the spec to sources: (now spec + Evidence), bumped source_commit→e8b2df0e, updated→2026-07-14. Backlog: moved the row to [[backlog]] Covered-via-concept → [[focal-amplification-detection]] (covered 123→124, pending 121→120); removed from [[backlog-pending]] (Oncology 27→26). No new page ⇒ index.md unchanged (concept already listed). No new graph nodes/edges (implementation prose only; concept's existing relates_to edges and body wikilinks unchanged) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No contradictions: spec and the recorded Evidence synthesis agree on model, thresholds, panel, and worked datasets.

- 2026-07-14 — ingest docs/algorithms/Oncology/Fusion_Breakpoint_Analysis.md (Fusion Breakpoint Analysis SPEC — primary per-algorithm spec for `OncologyAnalyzer.AnalyzeBreakpoint`/`PredictFusionProtein`, unit ONCO-FUSION-003, status Framework; Arriba site-vocabulary gating + AGFusion chimeric-CDS concat/translate/first-stop-truncation, codon-phase in-frame rule (fivePrimeCodingBases − threePrimeStartPhase) mod 3 == 0). CONTEXT/decision per Oncology reconciliation brief: REUSED the existing concept [[fusion-breakpoint-frame-and-protein-prediction]] (created 2026-07-10 from the ONCO-FUSION-003 Evidence artifact) rather than create a new page or a wiki/sources/ page (spec != Evidence report). The concept already covered the full model (four-state BreakpointFrameStatus, Arriba two-way vs AGFusion three-way framing, breakpoint-site gating, PredictFusionProtein steps, out-of-frame whole-codon trim, worked oracles). Enriched with genuinely-distinct implementation content the Evidence page lacked: an "Implementation surface" section (the two OncologyAnalyzer entry points + BreakpointAnalysis/FusionProteinPrediction record fields, offsets sourced from FusionBreakpoint, IsInFrame reuse from ONCO-FUSION-001, shared GeneticCode.Standard.Translate NCBI table 1, HasPrematureStop, O(1)/O(n) complexity, ArgumentNull/ArgumentOutOfRange validation, uppercase normalization, Framework status ⇒ caller-supplied CDS, slice/concat build so no suffix tree). Added the spec to sources: (now spec + Evidence), bumped source_commit→5465dd6b, updated→2026-07-14. Backlog: moved the row to [[backlog]] Covered-via-concept → [[fusion-breakpoint-frame-and-protein-prediction]] (covered 124→125, pending 120→119); removed from [[backlog-pending]] (Oncology 26→25, pending-total 119→118). No new page ⇒ index.md unchanged (concept already listed). No new graph nodes/edges (implementation prose only; concept's existing relates_to gene-fusion-detection-read-evidence / test-unit-registry edges and body wikilinks unchanged) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No hard contradiction; NOTE a benign framing nuance: the spec's §5.3 "intentionally simplified" describes AGFusion's "in-frame (with mutation)" middle class as mapped to InFrame (native-frame mid-codon junction, (b−p) mod 3 == 0), while the Evidence-derived concept describes the frameshifted-3′ complement case mapped to OutOfFrame ((b−p) mod 3 != 0) — these are two distinct junction geometries both handled correctly by the single codon-phase rule, not a rule disagreement.

- 2026-07-14 — lint | fixed 16 graph-provenance findings (5 typed edges re-sourced to existing source pages; 11 unsupported low-value typed edges removed while body wikilinks retain mentions), sharded the 547-line flat index into 6 bounded indexes, updated [[SCHEMA]], and rebuilt the graph (523 nodes / 4162 edges). Structural lint and graph lint clean.

- 2026-07-14 — lint follow-up | restored the 11 previously removed typed relationships at user request and rebuilt the graph (523 nodes / 4173 edges). The 5 provenance re-sources remain fixed; graph lint intentionally reports the restored 11 relationships because their `source:` values name concept pages rather than source pages.

- 2026-07-14 — ingest docs/algorithms/Oncology/HLA_Nomenclature_And_Allele_Specific_LOH.md (HLA nomenclature parsing + allele-specific HLA LOH SPEC — primary per-algorithm spec for `OncologyAnalyzer.ParseHlaAllele`/`TryParseHlaAllele`/`DetectHlaLoh`, unit ONCO-HLA-001, status Production; WHO IPD-IMGT/HLA nomenclature grammar validator + LOHHLA (McGranahan 2017) two-threshold LOH rule: allele CN < 0.5 AND allelic-imbalance paired-t p < 0.01, both strict). CONTEXT/decision per Oncology reconciliation brief: REUSED the existing concept [[hla-nomenclature-and-allele-specific-loh]] (created 2026-07-10 from the ONCO-HLA-001/ONCO-IMMUNE-001 Evidence artifacts) rather than create a new page or a wiki/sources/ page (spec != Evidence report). The concept already covered the full model (nomenclature grammar/fields/suffix set, validity rules, LOHHLA thresholds, worked oracles, over-calling guard, HomozygousLoss tie-break assumption, relation to allele-specific-copy-number-ascat). Enriched with genuinely-distinct implementation content the Evidence page lacked: an "API surface (implementation)" section (three OncologyAnalyzer entry points; HlaAllele/HlaLohResult/HlaAlleleCopyNumber records; LostAllele enum {None,Allele1,Allele2,Both}; ArgumentNull/Argument/Format exception contract + non-throwing Try wrapper; case-insensitive HLA-/suffix, upper-cased gene, leading zeros preserved; O(n) parse / O(1) LOH; suffix-only-when-trailing-letter; no suffix tree). Added the spec to sources: (now spec + 2 Evidence), bumped source_commit→305cb139, updated→2026-07-14. Backlog: moved the row to [[backlog]] Covered-via-concept → [[hla-nomenclature-and-allele-specific-loh]] (covered 126→127, pending 118→117); removed from [[backlog-pending]] (Oncology 24→23, pending-total 118→117). No new page ⇒ index.md unchanged (concept already listed). No new graph nodes/edges (implementation prose only; concept's existing relates_to test-unit-registry / allele-specific-copy-number-ascat edges and body wikilinks unchanged) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No contradictions — WHO nomenclature standard and LOHHLA paper/reference-R cover disjoint parts of the unit.

- 2026-07-14 — ingest docs/algorithms/Oncology/HRD_Score.md (HRD composite genomic-scar score SPEC — primary per-algorithm spec for `OncologyAnalyzer.CalculateHRDScore`/`ClassifyHRDStatus`/`DetectHRD`(3 overloads)/`CalculateHrdTaiScore`/`CalculateHrdLstScore`, unit ONCO-HRD-001, status Simplified; HRD = LOH + TAI + LST unweighted sum, HRD-high cutoff ≥42 inclusive (Telli 2016/Stewart 2022); LOH via DetectLOH/scarHRD calc.hrd, TAI via calc.ai_new even-ploidy path, LST via calc.lst; embedded UCSC cytoBand acen centromere table GRCh38/GRCh37). CONTEXT/decision per Oncology reconciliation brief: REUSED the existing concept [[homologous-recombination-deficiency-score]] (created 2026-07-10 from the ONCO-HRD-001 Evidence artifact) rather than create a new page or a wiki/sources/ page (spec != Evidence report). The concept already covered the full model (sum + 42 cutoff, LOH/TAI/LST definitions, TAI/LST derivation steps, centromere table, worked oracles, even-ploidy TAI assumption, relation to allele-specific-copy-number-ascat). Enriched with genuinely-distinct implementation content the Evidence page lacked: an "Implementation surface (OncologyAnalyzer)" section (the three-tier input model — count-based CalculateHRDScore/ClassifyHRDStatus/DetectHRD(HrdComponents); the caller-supplied-TAI/LST overload DetectHRD(segments,tai,lst) with LOH still derived; the all-derived DetectHRD(segments,genome) + standalone CalculateHrdTaiScore/CalculateHrdLstScore; INV-07 components identity; public constant HrdHighScoreThreshold=42; O(1)/O(n log n)/O(n²)-LST complexity; ArgumentOutOfRange/ArgumentNull/Argument exception contract). Added the spec to sources: (now spec + Evidence), bumped source_commit→9ccd313a, updated→2026-07-14. Backlog: moved the row to [[backlog]] Covered-via-concept → [[homologous-recombination-deficiency-score]] (covered 127→128, pending 117→116); removed from [[backlog-pending]] (Oncology 23→22, pending-total 117→116). No new page ⇒ index.md unchanged (concept already listed). No new graph nodes/edges (implementation prose only; no new wikilinks or typed edges) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No contradictions — the primary papers, scarHRD reference R, and UCSC/NCBI coordinate databases cover disjoint parts and agree on the sum, cutoff, and component rules.

- 2026-07-14 — ingest docs/algorithms/Oncology/Homozygous_Deletion_Detection.md (Homozygous/Deep Deletion Detection SPEC — primary per-algorithm spec for `OncologyAnalyzer.DetectHomozygousDeletions`/`IsHomozygousDeletion`/`IdentifyDeletedTumorSuppressors`, unit ONCO-CNA-003, status Production; homozygous deletion = total integer CN 0 = cBioPortal "−2" DeepDeletion (CNVkit absolute_threshold, log2 ≤ −1.1 default cutoffs), arm→tumour-suppressor panel TP53/RB1/CDKN2A/PTEN/BRCA1/BRCA2 by NCBI Gene cytoband). CONTEXT/decision per Oncology reconciliation brief: REUSED the existing concept [[homozygous-deletion-detection]] (created 2026-07-09 from the ONCO-CNA-003 Evidence artifact) rather than create a new page or a wiki/sources/ page (spec != Evidence report). The concept already covered the full model (CN-0 predicate with three-source convergence table, shallow-vs-deep, boundary-inclusive −1.1 rule, six-gene tumour-suppressor panel with 13q→RB1+BRCA2, worked dataset, order-preserving filter, both assumptions, total-CN-vs-allele-specific limitation). Enriched with genuinely-distinct implementation content the Evidence page lacked: an "Entry points" note (the three OncologyAnalyzer methods incl. the IsHomozygousDeletion single-segment predicate, reuse of ONCO-CNA-001 CallCopyNumber + ONCO-CNA-002 CopyNumberArmSegment/ValidateArmSegment, thresholds/ploidy optional params with defaults and constraints, O(n)/O(1) complexity) plus two corner cases (NaN log2 = CNVkit no-call → neutral reference CN, never reported; CallCopyNumber validation → ArgumentOutOfRange for ploidy ≤ 0 / Argument for non-ascending thresholds). Added the spec to sources: (now spec + Evidence), bumped source_commit→0472f265, updated→2026-07-14. Backlog: moved the row to [[backlog]] Covered-via-concept → [[homozygous-deletion-detection]] (covered 128→129, pending 116→115); removed from [[backlog-pending]] (Oncology 22→21, pending-total 116→115). No new page ⇒ index.md unchanged (concept already listed). No new graph nodes/edges (implementation prose only; concept's existing relates_to copy-number-alteration-classification / focal-amplification-detection / test-unit-registry edges and body wikilinks unchanged) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No contradictions: spec and the recorded Evidence synthesis agree on the CN-0 definition, cBioPortal/Cheng/CNVkit source convergence, panel, and worked datasets.

- 2026-07-14 — ingest docs/algorithms/Oncology/Immune_Infiltration_Estimation.md (Immune Infiltration Estimation SPEC — primary per-algorithm spec for `ImmuneAnalyzer.EstimateInfiltration`/`EstimateTumorPurity`/`DeconvoluteImmuneCells`(NNLS)/`DeconvoluteImmuneCellsNuSvr`(CIBERSORT ν-SVR)/`LoadSignatureMatrix`/`LoadBundledAbisSignatureMatrix`, unit ONCO-IMMUNE-001, status Simplified; three methods — ESTIMATE ssGSEA infiltration/purity, NNLS/LLSR deconvolution, CIBERSORT ν-SVR deconvolution). CONTEXT/decision per Oncology reconciliation brief: REUSED the existing concept [[immune-infiltration-deconvolution]] (created 2026-07-10 from the ONCO-IMMUNE-001 Evidence artifact) rather than create a new page or a wiki/sources/ page (spec != Evidence report). The concept already covered the full model (linear mixture m=S·f, ν-SVR sweep {0.25,0.5,0.75}+lowest-RMSE+zero-clip+sum-1, z-score standardisation, NNLS baseline, LM22/ABIS/default signature matrices + licences, Yoshihara cosine purity + Affymetrix domain + negative→NaN, MCP-counter, worked oracles, corner cases). Enriched with genuinely-distinct implementation content the Evidence page lacked: ν-SVR solver internals (SMO-style pairwise coordinate ascent on β_i=α_i−α_i*, step δ=(g_p−g_q)/(K_pp+K_qq−2K_pq) clipped to box |β_i|≤C and ν-budget Σ|β_i|≤Cνℓ, w=Σβ_i·x_i, C=1=NuSvrCost libsvm default, ≤200·n SMO iters, INV-NUSVR-04 dual constraints, O(|ν|·(n²+n·t·m)), non-finite→ArgumentException); the GSVA-style ssGSEA integral specifics (descending rank, hit weight ∝ rank^τ with τ=0.25, miss step −1/nMiss, integral not KS max-deviation, empty hit-set→0, O(N log N), un-normalised→relative purity vs opt-in absolute EstimateTumorPurity); edge cases (empty profile purity ≈0.8225, no-overlap all-zero/BestNu=0 branch, NNLS Lawson-Hanson maxIterations=1000, malformed TSV→FormatException). Added the spec to sources: (now spec + Evidence), bumped source_commit→e5e2f908, updated→2026-07-14. Backlog: moved the row to [[backlog]] Covered-via-concept → [[immune-infiltration-deconvolution]] (covered 129→130, pending 115→114); removed from [[backlog-pending]] (Oncology 21→20, pending-total 115→114). No new page ⇒ index.md unchanged (concept already listed). No new graph nodes/edges (implementation prose only; concept's existing relates_to test-unit-registry / expression-outlier-zscore-signature-score edges and body wikilinks unchanged) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No contradictions — the ESTIMATE/CIBERSORT/ν-SVR/GSVA/ABIS/LM22 sources cover complementary methods and agree on the shared linear-mixture and signature-scoring framing.

- 2026-07-14 — ingest docs/algorithms/Oncology/Known_Fusion_Database_Lookup.md (Known Fusion Database Lookup SPEC — primary per-algorithm spec for `OncologyAnalyzer.GetFusionAnnotation`/`MatchKnownFusions`/`KnownFusionMatch`/`FusionDesignationSeparator`, unit ONCO-FUSION-002, status Framework; HGNC gene-fusion designation `gene5p::gene3p` (Bruford et al. 2021 Leukemia) + directional lookup against a caller-supplied known-fusion set). CONTEXT/decision per Oncology reconciliation brief: REUSED the existing concept [[gene-fusion-nomenclature-known-fusion-lookup]] (created 2026-07-10 from the ONCO-FUSION-002 Evidence artifact) rather than create a new page or a wiki/sources/ page (spec != Evidence report); picked the DATABASE-LOOKUP fusion concept, distinct from the read-evidence caller [[gene-fusion-detection-read-evidence]] and the breakpoint/frame unit [[fusion-breakpoint-frame-and-protein-prediction]]. The concept already covered the full model (`::` separator, 5′-first directional order, HGNC symbols, Framework nature/caller-supplied set/no bundled Mitelman-COSMIC-ChimerDB, BCR::ABL1 + EML4::ALK worked examples, case-insensitivity, direction-matters + hyphen-read-through corner cases, validation contract). Enriched with genuinely-distinct implementation content the Evidence page lacked: an "Implementation (spec: ONCO-FUSION-002)" section — the four OncologyAnalyzer entry points incl. KnownFusionMatch record fields (Designation/IsKnown/Annotation?) and the FusionDesignationSeparator constant; the two-step matching mechanism (supplied dict comparer first → O(L) hash for OrdinalIgnoreCase, else single O(k·L) linear case-insensitive fallback scan so callers are not case/order-trapped); the deliberate suffix-tree rejection (exact dictionary-key lookup over short symbols, not substring search → hash map is correct). Added the spec to sources: (now spec + Evidence), bumped source_commit→9bf78352, updated→2026-07-14. Backlog: moved the row to [[backlog]] Covered-via-concept → [[gene-fusion-nomenclature-known-fusion-lookup]] (covered 130→131, pending 114→113); removed from [[backlog-pending]] (Oncology 20→19, pending-total 114→113). No new page ⇒ index.md unchanged (concept already listed). No new graph nodes/edges (implementation prose only; concept's existing relates_to gene-fusion-detection-read-evidence / test-unit-registry edges and body wikilinks unchanged) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No contradictions — the HGNC recommendation and the spec agree exactly on the `::` separator, 5′-first directionality, and the BCR::ABL1 worked example.

- 2026-07-14 — ingest docs/algorithms/Oncology/Loss_Of_Heterozygosity.md (Loss of Heterozygosity / HRD-LOH SPEC — primary per-algorithm spec for `OncologyAnalyzer.DetectLOH`/`CalculateHrdLohScore`/`CalculateLOHFraction`, unit ONCO-LOH-001, status Production; LOH segment = minor CN 0 & major CN ≠ 0, strict > 15 Mb length filter, whole-chromosome (chrDel) exclusion; HRD-LOH score = count of qualifying regions per Abkevich 2012 + scarHRD calc.hrd.R + oncoscanR score_loh). CONTEXT/decision per Oncology reconciliation brief: REUSED the existing concept [[loss-of-heterozygosity-detection]] (created 2026-07-10 from the ONCO-LOH-001 Evidence artifact) rather than create a new page or a wiki/sources/ page (spec != Evidence report). The concept already covered the full model (LOH-segment criterion, strict 15 Mb cut-off, homozygous-deletion exclusion, whole-chromosome chrDel exclusion, major-CN-cap-to-1 + oncoscanR ≤1 bp adjacency merge, length-weighted per-chromosome LOH fraction as API choice, worked oracle table, ordering invariance, upstream-segmentation assumption, relation to HRD/HLA-LOH siblings). Enriched with genuinely-distinct implementation content the Evidence page lacked: the third entry point `CalculateHrdLohScore` (direct score) + `DetectLOH` returning Regions+Score, the OncologyAnalyzer.cs location, complexity (DetectLOH O(n log n)/O(n), CalculateLOHFraction O(n)); and a correctness fix to the corner-cases list — split the conflated "empty/null → score 0" into empty→score 0 vs null (or null chromosome)→ArgumentNullException, added End≤Start / negative CN → ArgumentException, absent-chromosome fraction 0.0, and ordinal case-sensitive chromosome matching. Added the spec to sources: (now spec + Evidence), bumped source_commit→d9c52e6d, updated→2026-07-14. Backlog: moved the row to [[backlog]] Covered-via-concept → [[loss-of-heterozygosity-detection]] (covered 131→132, pending 113→112); removed from [[backlog-pending]] (Oncology 19→18, pending-total 113→112). No new page ⇒ index.md unchanged (concept already listed). No new graph nodes/edges (implementation prose only; concept's existing relates_to test-unit-registry / allele-specific-copy-number-ascat / homologous-recombination-deficiency-score edges and body wikilinks unchanged) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No contradictions — Abkevich 2012, scarHRD calc.hrd.R, and oncoscanR score_loh cover disjoint parts and agree on the LOH criterion, the strict 15 Mb cut-off, and the whole-chromosome exclusion.

- 2026-07-14 — ingest docs/algorithms/Oncology/MHC_Peptide_Binding_Classification.md (MHC-Peptide Binding Classification SPEC — primary per-algorithm spec for `OncologyAnalyzer.ClassifyBindingAffinity`/`ClassifyBindingRank`/`IsValidPeptideLength`/`ClassifyMhcBinding`/`PredictBindingHalfLifeBimas`/`PredictIc50Smm`/`PredictAndClassifySmm`/`LoadScoringMatrix` + `MhcflurryAffinityPredictor`, unit ONCO-MHC-001, status Framework; three layers — IC50/%Rank binder-tier classification (Strong/Weak/NonBinder, strict `<` cutoffs), opt-in matrix prediction (BIMAS product rule / SMM `IC50=50000^(1−score)`), and a ported MHCflurry 2.0 pan-allele class-I network). CONTEXT/decision per Oncology reconciliation brief: REUSED the existing concept [[mhc-peptide-binding-prediction]] (created 2026-07-10 from the ONCO-MHC-001 Evidence artifact) rather than create a new page or a wiki/sources/ page (spec != Evidence report). The concept already covered the full model (IC50 50/500 nM tiers, %Rank class I 0.5/2 & class II 2/10, strict-`<` boundary semantics, class I 8–14 / class II 13–25 lengths, BIMAS product & SMM additive rules with worked anchors, MHCflurry BLOSUM62 945/777 encoding + tanh/sigmoid + geometric-mean ensemble + <0.03% oracle parity, caller-supplied-trained-data packaging boundary analogized to CIBERSORT LM22, neoantigen/HLA-LOH relations). Enriched with genuinely-distinct implementation content the Evidence page lacked: the `LoadScoringMatrix(IEnumerable<string>)` parse format (a `CONST=<value>` line + `RESIDUE=VALUE` rows, `ArgumentNullException` on null input, `FormatException` on a malformed token / non-numeric value / multi-character residue key), the `PmhcScoringMatrix` record + `PmhcScoringMethod` enum, and INV-07 matrix-predictor contract (empty matrix or length≠rows → ArgumentException, null peptide → ArgumentNullException). Added the spec to sources: (now spec + Evidence), bumped source_commit→f44fa40d, updated→2026-07-14. Backlog: moved the row to [[backlog]] Covered-via-concept → [[mhc-peptide-binding-prediction]] (covered 132→133, pending 112→111); removed from [[backlog-pending]] (Oncology 18→17, pending-total 112→111). No new page ⇒ index.md unchanged (concept already listed). No new graph nodes/edges (implementation prose only; concept's existing relates_to test-unit-registry / hla-nomenclature-and-allele-specific-loh / immune-infiltration-deconvolution / neoantigen-peptide-generation edges and body wikilinks unchanged) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No contradictions — the NetMHCpan %Rank tiers, Sette/IEDB IC50 tiers, BIMAS/Parker & SMM matrix rules, and the MHCflurry network cover complementary methods and agree on the shared affinity/percentile framing.

- 2026-07-14 — ingest docs/algorithms/Oncology/MRD_Detection.md (Minimal/Molecular Residual Disease Detection SPEC — primary per-algorithm spec for `OncologyAnalyzer.DetectMRD`/`TrackVariantsOverTime`/`IsVariantDetected`/`EstimateInvarSignal`/`IntegratedMutantAlleleFractionV2`/`EstimateInvarSignalWithSize`/`FragmentSizeProfile`(+`FromKernelDensity`)/`SuppressOutlierLoci`/`EstimateLocusBackground`/`PassesBothStrandsFilter`, unit ONCO-MRD-001, status Simplified; tumour-informed ctDNA MRD: Signatera ≥2-of-16 panel positivity rule + read-pooled IMAF + panel Poisson p=1−e^(−nfm) reusing ONCO-CTDNA-001, and the INVAR2 GLRT stack — AF-weighted per-locus mixture q=p·g+e(1−p), EM p̂, LR=logL(p̂)−logL(0), background-subtracted IMAFv2, fragment-size-weighted with-RL GLRT, Bonferroni outlier suppression, control-derived background + both-strands filter, opt-in Gaussian-KDE size profile). CONTEXT/decision per Oncology reconciliation brief: REUSED the existing concept [[tumor-informed-mrd-detection]] (created 2026-07-10 from the ONCO-MRD-001 Evidence artifact) rather than create a new page or a wiki/sources/ page (spec != Evidence report). The concept already covered the full model (Signatera ≥2/16 rule + detected-of-16 table + HR 7.2 clinical signal, longitudinal first-positive, panel Poisson reuse with worked λ, IMAF/IMAFv2, INVAR mixture/EM/GLRT with synthetic-recovery oracles, size weighting + flat-profile sanity, outlier suppression + control-background + both-strands oracles, KDE size profile with Silverman bandwidth, corner cases, per-variant-detected assumption + resolved KDE note). Enriched with a genuinely-distinct "Implementation surface (ONCO-MRD-001)" section the Evidence page lacked: OncologyAnalyzer.cs location + the eleven entry points incl. the ONCO-CTDNA-001 CtDnaDetectionProbability reuse; the TumorMarker/InvarLocus/InvarMolecule inputs and MrdResult{Status/DetectedVariantCount/TrackedVariantCount/IMAF/DetectionProbability} + InvarSignalResult{IMAFv2/EstimatedTumorFraction p̂/LikelihoodRatio/Detected/LocusCount} outputs; parameter defaults (τ=2, r_min=1, n=0, detectionThreshold=0 with constraints); the validation contract (null→ArgumentNullException, empty→ArgumentException, out-of-range params→ArgumentOutOfRangeException, TrackVariantsOverTime delegates per timepoint); numerical details (negative reads clamped to 0 in IMAF, zero background floored to 1/R_i for finite logs, informative loci AF>0 & R>0, Lanczos lchoose/gamma); complexity (DetectMRD O(m), TrackVariantsOverTime O(T·m), EstimateInvarSignal O(m·I) I=200, IMAFv2 O(m)); suffix tree N/A (positional marker matching); and the CHIP scope boundary → ONCO-CHIP-001 [[clonal-hematopoiesis-cfdna-filtering]] (FilterCHIP). Added the spec to sources: (now spec + Evidence), bumped source_commit→5fcdcf5a, updated→2026-07-14. Backlog: moved the row to [[backlog]] Covered-via-concept → [[tumor-informed-mrd-detection]] (covered 133→134, pending 111→110); removed from [[backlog-pending]] (Oncology 17→16, pending-total 111→110). No new page ⇒ index.md unchanged (concept already listed). No new graph nodes/edges (implementation prose only; concept's existing relates_to test-unit-registry / depends_on ctdna-detection-and-tumor-fraction edges and body wikilinks — incl. the newly-referenced clonal-hematopoiesis-cfdna-filtering already linked elsewhere — unchanged) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No contradictions — Reinert 2019/Signatera, Natera/Avanzini Poisson, Wan 2020, the INVAR2 reference R, and Silverman KDE cover disjoint stages and agree with the recorded Evidence synthesis on the ≥2 rule, Poisson panel model, and GLRT/EM/IMAFv2/size/outlier/background formulas.

- 2026-07-14 — ingest docs/algorithms/Oncology/Microsatellite_Instability_Detection.md (MSI Detection SPEC — primary per-algorithm spec for `OncologyAnalyzer.CalculateMSIScore`/`ClassifyMSIStatus`/`ClassifyBethesdaPanel`/`DetectMSI`, unit ONCO-MSI-001, status Simplified; scoring-and-classification layer only: MSI score = unstable/valid loci with MSIsensor2 binary MSI-H cutoff ≥20% inclusive, plus categorical NCI/Bethesda 5-marker rule 0→MSS/1→MSI-L/≥2→MSI-H; upstream per-locus chi-square tumour-vs-normal call out of scope). CONTEXT/decision per Oncology reconciliation brief: REUSED the existing concept [[microsatellite-instability-detection]] (created 2026-07-10 from the ONCO-MSI-001 Evidence artifact) rather than create a new page or a wiki/sources/ page (spec != Evidence report). The concept already covered the full model (continuous MSIsensor/MSIsensor2 fraction + 20% inclusive cutoff + 3.5% dataset-specific note, per-site chi-square/FDR 0.05 upstream, Bethesda 5-marker panel table + BAT-25/BAT-26/D2S123/D5S346/D17S250, two-inputs-two-classifiers modelling choice with no fabricated continuous MSI-L band, worked oracles, corner cases, MSS-vs-MSI-L ambiguity, TMB/HRD/interpretation-layer relations). Enriched with a genuinely-distinct "Implementation surface (ONCO-MSI-001 spec)" section the Evidence page lacked: the four OncologyAnalyzer entry points with signatures, the MsiResult record struct{UnstableLoci/TotalLoci/Score/Status} + MsiStatus enum (MSS/MSI_Low/MSI_High, with ClassifyMSIStatus never returning MSI_Low), the full validation contract (ArgumentOutOfRangeException on totalLoci≤0 / unstableLoci<0 / unstableLoci>totalLoci / non-finite or out-of-[0,1] score / Bethesda marker violations; ArgumentNullException on null and ArgumentOutOfRangeException on empty flags in DetectMSI; inclusive 20% cutoff), and complexity (Classify*/CalculateMSIScore O(1), DetectMSI O(n)/O(1), suffix tree N/A). Added the spec to sources: (now spec + Evidence), bumped source_commit→335bdb80, updated→2026-07-14. Backlog: moved the row to [[backlog]] Covered-via-concept → [[microsatellite-instability-detection]] (covered 134→135, pending 110→109); removed from [[backlog-pending]] (Oncology 16→15, pending-total 110→109). No new page ⇒ index.md unchanged (concept already listed). No new graph nodes/edges (implementation prose only; concept's existing relates_to test-unit-registry edge and body wikilinks unchanged) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No contradictions — Niu et al. 2014 (MSIsensor), the niu-lab msisensor2 README, and Boland et al. 1998 (NCI Bethesda) cover disjoint parts and agree with the recorded Evidence synthesis on the unstable-loci fraction, the 20% inclusive cutoff, and the 5-marker marker-count rule.

- 2026-07-14 — ingest docs/algorithms/Oncology/Mutational_Process_Classification.md (Mutational Process Classification SPEC — primary per-algorithm spec for `OncologyAnalyzer.ClassifyMutationalProcess`/`GetMutationalProcess`, unit ONCO-SIG-004, status Production; maps signature-fitting exposures → active mutagenic aetiologies via normalize (Wᵢ=eᵢ/Σe) → deconstructSigs 6% presence cutoff (strict `<`, 0.06 retained) → COSMIC SBS→aetiology map (SBS1/5→Aging, SBS2/13→APOBEC, SBS4→Tobacco, SBS7a–d→UV, SBS6/15/20/26→MMRd) → per-process sum → dominant argmax). CONTEXT/decision per Oncology reconciliation brief: REUSED the existing concept [[mutational-process-classification]] (created 2026-07-10 from the ONCO-SIG-004 Evidence artifact) rather than create a new page or a wiki/sources/ page (spec != Evidence report), picking the PROCESS/aetiology-classification sibling over the fitting/extraction, exposure-bootstrap and SBS96-catalog concepts. The concept already covered the full model (normalize→cutoff→map→aggregate pipeline pinned per-source, worked SBS2/13/1/4 oracle, sub-cutoff-mass-dropped/multiple-active-processes/unmapped-label/Σ=0/custom-cutoff corner cases, per-process-summation + per-signature-cutoff-then-group assumptions, ONCO-SIG-001/002/003 family relations). Enriched with a genuinely-distinct "Implementation contract (ONCO-SIG-004 API)" section the Evidence page lacked: the two OncologyAnalyzer entry points with signatures, ActiveProcesses (IReadOnlyList<ProcessActivity>, descending contribution then process enum) + DominantProcess (MutationalProcess enum, Unknown when none) outputs, case-insensitive label lookup, the validation contract (null exposures/label→ArgumentNullException, negative/NaN exposure→ArgumentException, cutoff NaN or outside [0,1)→ArgumentOutOfRangeException, empty/zero-total→empty set+Unknown), complexity O(k log k)/O(k) (log k = ordering ≤5 processes), suffix-tree N/A, and the not-implemented confidence-based presence → ONCO-SIG-003 [[signature-exposure-bootstrap-confidence-intervals]]. Added the spec to sources: (now spec + Evidence), bumped source_commit→7783b8d6, updated→2026-07-14. Backlog: moved the row to [[backlog]] Covered-via-concept → [[mutational-process-classification]] (covered 135→136, pending 109→108); removed from [[backlog-pending]] (Oncology 15→14, pending-total 109→108). No new page ⇒ index.md unchanged (concept already listed). No new graph nodes/edges (implementation prose only; concept's existing relates_to test-unit-registry / depends_on mutational-signature-fitting-and-extraction edges and body wikilinks unchanged) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No contradictions — deconstructSigs (Rosenthal 2016) + whichSignatures.R, COSMIC SBS aetiology strings, and Alexandrov 2020 cover disjoint parts and agree with the recorded Evidence synthesis on normalized relative contributions, the 6% strict-`<` cutoff, and the SBS→aetiology map.

- 2026-07-14 — ingest docs/algorithms/Oncology/Mutational_Signature_Exposure_Bootstrap.md (Mutational Signature Exposure Bootstrap CI SPEC — primary per-algorithm spec for `OncologyAnalyzer.BootstrapExposures` + `BootstrapResampling` enum + `ExposureConfidenceInterval` record, unit ONCO-SIG-003, status Framework; parametric bootstrap of NNLS signature exposures: resample the observed catalog R times → refit each by NNLS (reuses ONCO-SIG-002 FitSignatures) → per-signature type-7 percentile CI, with two schemes — fixed-N Multinomial (sigminer default) and per-channel Poisson(observedₖ) (Senkin 2021 MSA variant, N not fixed)). CONTEXT/decision per Oncology reconciliation brief: REUSED the existing concept [[signature-exposure-bootstrap-confidence-intervals]] (created 2026-07-10 from the ONCO-SIG-003 Evidence artifact) rather than create a new page or a wiki/sources/ page (spec != Evidence report). The concept already covered the full model (resample→refit→percentile pipeline pinned per-source, type-7 quantile with worked oracles, multinomial-vs-Poisson fixed/unfixed-N distinction + the single-channel deterministic-collapse-vs-positive-width discriminator, N=0/zero-count/R=1 corner cases, ordering/non-negativity/determinism invariants + full validation contract, family relations to ONCO-SIG-001/002 and clinical-interpretation units). Enriched with a genuinely-distinct "Implementation (entry point and samplers)" section the Evidence page lacked: the `BootstrapExposures(catalog, signatures, replicates=1000, confidence=0.95, seed=42, resampling=Multinomial)` signature + `ExposureConfidenceInterval{PointEstimate,Mean,Lower,Upper,Confidence}` fields + `BootstrapResampling` enum default preserving byte-for-byte behaviour; the two sampler realisations (multinomial via sequential conditional-binomial construction with Bernoulli-sum Binomial — not R sample+table; Poisson via Knuth multiplication-of-uniforms, λ=0→0); complexity O(R·(N+NNLS(n,k))) + O(R log R) percentile sort; and the not-implemented MSA Gaussian σ=10% noise model / Bayesian credible intervals / presence p-values (sigfit/signeR/sigminer report_bootstrap_p_value). Added the spec to sources: (now spec + Evidence), bumped source_commit→cabf04f1, updated→2026-07-14. Backlog: moved the row to [[backlog]] Covered-via-concept → [[signature-exposure-bootstrap-confidence-intervals]] (covered 136→137, pending 108→107); removed from [[backlog-pending]] (Oncology 14→13, pending-total 108→107). No new page ⇒ index.md unchanged (concept already listed). No new graph nodes/edges (implementation prose only; concept's existing relates_to test-unit-registry / depends_on mutational-signature-fitting-and-extraction edges and body wikilinks unchanged) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No contradictions — Senkin 2021 (MSA), sigminer sig_fit_bootstrap, Huang/Wojtowicz/Przytycka 2018, Efron 1979, Hyndman & Fan 1996, Lawson & Hanson 1974 and Knuth cover complementary parts and agree with the recorded Evidence synthesis on the resample→NNLS-refit→percentile pipeline, the two count-resampling schemes, and the type-7 percentile CI.

- 2026-07-14 — ingest docs/algorithms/Oncology/Mutational_Signature_Extraction_NMF.md (De-novo Mutational-Signature Extraction via NMF SPEC — primary per-algorithm spec for `OncologyAnalyzer.ExtractSignatures` + `NmfObjective` enum + `SelectRank` + `MatchToReferenceSignatures`, unit ONCO-SIG-002, status Production; V≈W·H non-negative factorisation of a channels×samples count matrix into signatures W and exposures H via Lee & Seung multiplicative updates, Frobenius (Theorem 1) or KL/Poisson (Theorem 2, SigProfiler choice), L1-column-normalised signatures, Brunet-2004 consensus/cophenetic + Rousseeuw silhouette rank selection, cosine reference matching). CONTEXT/decision per Oncology reconciliation brief: REUSED the existing concept [[mutational-signature-fitting-and-extraction]] (created 2026-07-10 from the ONCO-SIG-002 Evidence artifact — the shared unit page already covering BOTH the NNLS refit and this NMF extraction) rather than create a new page or a wiki/sources/ page (spec != Evidence report). The concept already covered the full extraction model (V≈WH blind-source-separation framing, Frobenius + KL Lee & Seung updates, L1 normalisation + scale absorption, non-convex/local-optimum + permutation-scale ambiguity + ε denominator guard sharp edges, consensus/cophenetic + silhouette rank selection, greedy cosine COSMIC matching). Enriched with a genuinely-distinct "API surface and implementation notes" section the Evidence page lacked: the three OncologyAnalyzer entry points with signatures + defaults (ExtractSignatures maxIterations=10_000/tolerance=1e-10/seed=42, Frobenius default & preserved 5-arg overload; SelectRank runs=20/stabilityThreshold=0.80/minStability=0.20, KL default; MatchToReferenceSignatures), the SignatureExtractionResult/RankSelectionResult/RankStability/SignatureMatch records, opaque-channel V[channel][sample] contract, NmfEpsilon=1e-12 denominator+init floor, the relative-improvement convergence test (prevObj−obj)/max(prevObj,1) with tolerance=0 for planted-recovery tests, O(m·k·n)/O(I·m·k·n) complexity; and tightened the rank-selection rule to the spec's explicit "largest k with avg stability ≥0.80 and min ≥0.20, else highest-average-stability; cophenetic + mean reconstruction error are per-rank diagnostics not the selector" plus the no-embedded-COSMIC caller-supplied-reference note. Added the spec to sources: (now spec + Evidence), bumped source_commit→d69205e3, updated→2026-07-14. Backlog: moved the row to [[backlog]] Covered-via-concept → [[mutational-signature-fitting-and-extraction]] (covered 137→138, pending 107→106); removed from [[backlog-pending]] (Oncology 13→12, pending-total 108→107). No new page ⇒ index.md unchanged (concept already listed). No new graph nodes/edges (implementation prose only; concept's existing relates_to test-unit-registry / depends_on sbs96-mutational-signature-catalog edges and body wikilinks unchanged) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No contradictions — Lee & Seung 2001, Alexandrov 2013/2020, SigProfilerExtractor (Islam 2022), Brunet 2004 and Rousseeuw 1987 cover complementary parts and agree with the recorded Evidence synthesis on the V≈WH factorisation, the two objectives, L1 normalisation, and the consensus/silhouette rank rule.

- 2026-07-14 — ingest docs/algorithms/Oncology/Mutational_Signature_Fitting.md (Mutational Signature Fitting / NNLS Refitting + Cosine Similarity SPEC — primary per-algorithm spec for `OncologyAnalyzer.FitSignatures` + `CosineSimilarity` + `ReconstructCatalog`, unit ONCO-SIG-002, status Framework; the supervised REFIT half of the shared signature-deconvolution unit: min_x ‖S·x − d‖₂², x≥0 solved by Lawson-Hanson active set, cosine reconstruction gate ≥0.95, exposure proportions x/Σx, caller-supplied reference signatures — no embedded COSMIC). CONTEXT/decision per Oncology reconciliation brief: REUSED the existing concept [[mutational-signature-fitting-and-extraction]] (created 2026-07-10 from the ONCO-SIG-002 Evidence artifact — the shared unit page already covering BOTH this NNLS refit and the NMF extraction) rather than create a new page or a wiki/sources/ page (spec != Evidence report). The concept already covered the full refit model (NNLS objective + active-set clamp-and-refit, gradient w=Aᵀ(y−Ax), cosine similarity CS∈[0,1] with zero-norm→0.0, reconstruction S·x + 0.95 gate + proportion form, the three worked NNLS oracles including the constraint-binding x₁=−1→clamp→[0,0.5] case, family relations to ONCO-SIG-001/003/004). Enriched the "API surface and implementation notes" section — which previously listed only the EXTRACTION-layer entry points — with a genuinely-distinct **refit** layer the concept lacked: the three OncologyAnalyzer entry points (`CosineSimilarity(a,b)` zero-norm→0.0; `FitSignatures(catalog,signatures)`→`SignatureFitResult{Exposures,NormalizedExposures,Reconstruction,ReconstructionCosineSimilarity}`; `ReconstructCatalog(signatures,exposures)`), the dense passive-set normal-equations inner solve s_P=((S_P)ᵀS_P)⁻¹(S_P)ᵀd by Gaussian elimination with partial pivoting (ε=1e-12), the singular-passive-set/collinear-signatures→component-stays-0-not-throw behaviour, and the O(k³+k²·n) per outer iteration over ≤O(k) outer iterations complexity. Added the spec to sources: (now Evidence + Fitting spec + Extraction_NMF spec), bumped source_commit→c559752e, updated 2026-07-14. Backlog: moved the row to [[backlog]] Covered-via-concept → [[mutational-signature-fitting-and-extraction]] (covered 138→139, pending 106→105); removed from [[backlog-pending]] (Oncology 12→11, pending-total 107→106). No new page ⇒ index.md unchanged (concept already listed). No new graph nodes/edges (implementation prose only; concept's existing relates_to test-unit-registry / depends_on sbs96-mutational-signature-catalog edges and body wikilinks unchanged) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No contradictions — MutationalPatterns (Blokzijl 2018), deconstructSigs (Rosenthal 2016), Lawson-Hanson 1974 and iMutSig (Pan & Wang 2020) cover complementary parts and agree with the recorded Evidence synthesis on the NNLS objective, the active-set clamp-and-refit, the cosine reconstruction gate, and the proportion normalisation.

- 2026-07-14 — ingest docs/algorithms/Oncology/Neoantigen_Peptide_Generation.md (Neoantigen Candidate Peptide Window Generation SPEC — primary per-algorithm spec for `OncologyAnalyzer.GenerateNeoantigenPeptides` + `NeoantigenPeptide` record struct + `MhcClassIMinPeptideLength`/`MhcClassIMaxPeptideLength` constants, unit ONCO-NEO-001, status Framework; deterministic windowing step: for a somatic missense substitution at 1-based position p in protein P, enumerate every length-k window (class I default k=8..11, MHC-supported 8..14) of the mutant protein P′ that spans p — start s∈[max(1,p−k+1),min(p,L−k+1)] — each paired with the matched wild-type k-mer at identical coordinates (the agretope) plus the 0-based mutation offset; no MHC binding/IC50 scoring, no frameshift/indel/fusion neopeptides). CONTEXT/decision per Oncology reconciliation brief: REUSED the existing concept [[neoantigen-peptide-generation]] (created 2026-07-10 from the ONCO-NEO-001 Evidence artifact) rather than create a new page or a wiki/sources/ page (spec != Evidence report). The concept already covered the full model (mutation-spanning window rule with 0-based start bounds + interior count=k, mutant/WT agretope pairing + differential agretopic index framing, class I 8–11 default vs NetMHCpan 8–14 + 9-mer dominance + ProGeo-neo 21-mer ±10-flank equivalence, worked Y5C interior + M1V terminal-truncation + synonymous-no-peptide + short-protein corner cases, missense-only + binding-out-of-scope scope, relations to [[mhc-peptide-binding-prediction]] affinity gate and [[hla-nomenclature-and-allele-specific-loh]] presentation platform). Enriched with a genuinely-distinct "Implementation surface" section the Evidence page lacked: the `GenerateNeoantigenPeptides(wildTypeProtein, mutantResidue, mutationPosition, minLength=8, maxLength=11)`→IReadOnlyList<NeoantigenPeptide> entry point ordered length-then-start; the NeoantigenPeptide record-struct fields (Length/StartPosition 1-based/MutantPeptide/WildTypePeptide/MutationOffset 0-based) + the 8/11 class I constants; the full validation/exception contract (null→ArgumentNullException; empty protein, mutantResidue==WT, minLength<1, maxLength<minLength→ArgumentException; mutationPosition∉[1,L]→ArgumentOutOfRangeException; k>L silently skipped; opaque one-letter strings, no alphabet validation, case preserved); the six invariants INV-01..06 covered by OncologyAnalyzer_GenerateNeoantigenPeptides_Tests.cs; and the suffix-tree-N/A deviation (bounded arithmetic range, not multi-query exact match) + OncologyAnalyzer-over-NeoantigenPredictor naming. Added the spec to sources: (now spec + Evidence), bumped source_commit→ce89ed9d, updated→2026-07-14. Backlog: moved the row to [[backlog]] Covered-via-concept → [[neoantigen-peptide-generation]] (covered 139→140, pending 105→104); removed from [[backlog-pending]] (Oncology 11→10, pending-total 106→105). No new page ⇒ index.md unchanged (concept already listed). No new graph nodes/edges (implementation prose only; concept's existing relates_to test-unit-registry / mhc-peptide-binding-prediction / hla-nomenclature-and-allele-specific-loh edges and body wikilinks unchanged) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No contradictions — pVACtools (Hundal 2020), ProGeo-neo (Li 2020), NetMHCpan-4.0/4.1 (Jurtz 2017/DTU) and Wells 2020 (TESLA) cover complementary parts and agree with the recorded Evidence synthesis on the mutation-spanning-window rule, the matched-WT agretope pairing, and the class I length band.

- 2026-07-14 — ingest docs/algorithms/Oncology/SBS96_Trinucleotide_Context_Catalog.md (SBS-96 Trinucleotide Context Catalog SPEC — primary per-algorithm spec for `OncologyAnalyzer.ClassifySbsContext` + `EnumerateSbs96Channels` + `Build96ContextCatalog`, unit ONCO-SIG-001, status Production; the foundational SBS classification step: fold each single-base substitution onto the pyrimidine strand and tally into the 96 = 6×4×4 canonical COSMIC channels labelled `5'[REF>ALT]3'`). CONTEXT/decision per Oncology reconciliation brief: REUSED the existing concept [[sbs96-mutational-signature-catalog]] (created 2026-07-10 from the ONCO-SIG-001 Evidence artifact) rather than create a new page or a wiki/sources/ page (spec != Evidence report). The concept already covered the full model (6×4×4=96 channel definition, pyrimidine-strand reverse-complement fold rule with worked TGA/G>T→T[C>A]A oracle and a 7-row oracle table, partition invariants Σcounts=classifiable-variants, purine-fold / non-SBS / non-ACGT / ref==alt corner cases, family relations to the fitting/HRD/MSI/interpretation units). Enriched with a genuinely-distinct "Implementation surface (ONCO-SIG-001 spec)" section the Evidence page lacked: the three OncologyAnalyzer entry points with signatures + return types (ClassifySbsContext O(1); EnumerateSbs96Channels O(96) substitution-major deterministic order; Build96ContextCatalog O(n) with all 96 keys always present incl. zero counts for fixed 96-dim vector shape, Ordinal-keyed dict), the case-insensitive upper-casing + validation/exception contract (non-ACGT or ref==alt → ArgumentException; null → ArgumentNullException), and the constant-time-classification / suffix-tree-N/A note. Added the spec to sources: (now spec + Evidence), bumped source_commit→a0d6cd87, updated→2026-07-14. Backlog: moved the row to [[backlog]] Covered-via-concept → [[sbs96-mutational-signature-catalog]] (covered 140→141, pending 104→103); removed from [[backlog-pending]] (Oncology 10→9). No new page ⇒ index.md unchanged (concept already listed). No new graph nodes/edges (implementation prose only; concept's existing relates_to test-unit-registry edge and body wikilinks unchanged) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No contradictions — Alexandrov 2013, COSMIC SBS96, SigProfilerMatrixGenerator (Bergstrom 2019) cover complementary parts and agree with the recorded Evidence synthesis on the 6×4×4 channel space and the pyrimidine reverse-complement fold.

- 2026-07-15 — ingest docs/algorithms/Oncology/Sequencing_Artifact_Detection.md (Sequencing Artifact Detection SPEC — primary per-algorithm spec for `OncologyAnalyzer.ClassifyArtifact` + `CalculateGivScore` + `CalculateStrandBias` + `FilterArtifacts` + `DetectOxoGArtifacts`, unit ONCO-ARTIFACT-001, status Production; deterministic rule-based classifier: FFPE cytosine-deamination C>T/G>A vs OxoG 8-oxoG oxidation G>T/C>A disjoint substitution classes, GIV = read1Alt/read2Alt with DamagedGivThreshold 1.5 / UndamagedGivScore 1.0, and Phred-scaled two-sided Fisher strand bias FS = -10·log10(max(p, 1e-320)) on the [refFwd,refRev,altFwd,altRev] table). CONTEXT/decision per Oncology reconciliation brief: REUSED the existing concept [[sequencing-artifact-detection]] (created 2026-07-09 from the ONCO-ARTIFACT-001 Evidence artifact) rather than create a new page or a wiki/sources/ page (spec != Evidence report). The concept already covered the three-signal model (disjoint substitution classes, GIV ratio + 1.5/1.0 thresholds with Damage-estimator ≈2 acoustic-shearing note, GATK FisherStrand FS + MIN_PVALUE floor, worked oracles, FilterArtifacts output ⊆ input, GIV zero-denominator + zero-margin-table edge cases, BAM-parser-absent API-shape note, research-grade/not-for-clinical scope). Enriched with genuinely-distinct implementation content the Evidence page lacked: a new "Flag decision and entry points" section stating the exact flag rule (FFPE always flagged as a conservative pre-filter; OxoG flagged iff GIV > 1.5; other → never), the five OncologyAnalyzer entry points with signatures (ClassifyArtifact→ArtifactCall{Type,GivScore,StrandBiasPhred,IsArtifact}, CalculateGivScore, CalculateStrandBias, FilterArtifacts drops-flagged-keeps-order, DetectOxoGArtifacts returns flagged OxoG), the Lanczos log-gamma summed-hypergeometric two-sided Fisher (O(N) column sum, no suffix tree) + the three named constants; upgraded the [20,0,0,20] oracle from "FS large (>0)" to the exact p=1.4508889×10⁻¹¹ → FS=108.384; and added the GATK Mutect2 LearnReadOrientationModel/FilterMutectCalls (F1R2 OrientationBiasFilter) probabilistic-posterior pointer + F1R2/F2R1-tensor-vs-per-strand-count input contrast to Scope and limitations. Added the spec to sources: (now Evidence + spec), bumped source_commit→da317d1a, updated→2026-07-15. Backlog: moved the row to [[backlog]] Covered-via-concept → [[sequencing-artifact-detection]] (covered 141→142, pending 103→102); removed from [[backlog-pending]] (Oncology 9→8, pending-total 105→104). No new page ⇒ index.md unchanged (concept already listed). No new graph nodes/edges (implementation prose only; concept's existing relates_to test-unit-registry edge and body wikilinks unchanged) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No contradictions — Chen 2017 / Damage-estimator (oxidation GIV), Do & Dobrovic 2015 (FFPE deamination), Nature Methods 2017 (GIV thresholds) and GATK FisherStrand (strand bias) each cover a disjoint signal and agree with the recorded Evidence synthesis.
- 2026-07-15 — ingest docs/algorithms/Oncology/Somatic_Mutation_Calling.md (Somatic Mutation Calling SPEC — primary per-algorithm spec for `OncologyAnalyzer.CallSomaticMutations`/`Classify`/`FilterGermlineVariants`/`CalculateSomaticScore`, unit ONCO-SOMATIC-001, status Simplified; deterministic rule-based tumor-vs-matched-normal VAF classifier — Somatic/Germline/NotDetected via f_t≥τ_t (0.05) and f_n≤τ_n (0.01), somatic score max(0,f_t−f_n)). CONTEXT/decision per Oncology reconciliation brief: REUSED the existing concept [[somatic-variant-calling-tumor-normal]] (created 2026-07-10 from the ONCO-SOMATIC-001 Evidence artifact) rather than create a new page or a wiki/sources/ page (spec != Evidence report). The concept already covered the full model (Strelka S={f_t≠f_n} ref/ref state, Strelka2 continuous VAF, three thresholds with Yan-2021 5% WES LoD + 1% normal ceiling, Mutect2 germline filter + tumor-only ℓ_n=1 mode, max(0,f_t−f_n) documented-simplification score, seven worked oracles incl. boundary conventions, FilterGermlineVariants subset contract, LOH/low-purity/sub-5% corner cases, two flagged assumptions, research-grade/not-for-clinical scope). Enriched section 5 with genuinely-distinct implementation content the Evidence page lacked: the four `OncologyAnalyzer` entry points with the one-pass **O(n)** order-preserving `CallSomaticMutations` complexity, pure in-memory `VariantObservation` classification (no VCF parse, no suffix-tree), and the INV-05 uncovered-site `totalReads=0 ⇒ VAF=0` invariant plus the null/threshold-range/alt>total contract failures. Added the spec to sources: (now spec + Evidence×2), bumped source_commit→77ba0282, updated→2026-07-15. Backlog: moved the row to [[backlog]] Covered-via-concept → [[somatic-variant-calling-tumor-normal]] (covered 142→143, pending 102→101); removed from [[backlog-pending]] (Oncology 8→7, pending-total 104→103). No new page ⇒ index.md unchanged (concept already listed). No new graph nodes/edges (implementation prose only; concept's existing relates_to edges and body wikilinks unchanged) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No contradictions — Strelka/Strelka2 (Saunders 2012 / Kim 2018), Mutect2 (Benjamin 2019) and Yan 2021 each cover a disjoint facet and agree with the recorded Evidence synthesis; the spec's rule-based realization matches the concept's decision logic verbatim.
- 2026-07-15 — ingest docs/algorithms/Oncology/Tumor_Gene_Expression_Outlier.md (Tumor Gene Expression Outlier + Signature Score SPEC — primary per-algorithm spec for `OncologyAnalyzer.CalculateExpressionZScore`/`IdentifyOutlierGenes`/`CalculateSignatureScore`, unit ONCO-EXPR-001, status Framework; per-gene outlier z = (r−μ)/σ with sample SD divisor n−1 per cBioPortal NormalizeExpressionLevels.java, strict ±2 outlier rule z>+t Over / z<−t Under, and Lee et al. 2008 combined-z signature activity a = (Σz)/√k). CONTEXT/decision per Oncology reconciliation brief: REUSED the existing concept [[expression-outlier-zscore-signature-score]] (created 2026-07-10 from the ONCO-EXPR-001 Evidence artifact) rather than create a new page or a wiki/sources/ page (spec != Evidence report). The concept already covered the full model (z-score with (n−1) SD settled by the reference impl, diploid-vs-all base-population conventions, strict ±2 threshold, √k combined-z with k=1→z₁, zero-SD-throws behavioural deviation, n≥2 / empty-signature / scale-assumption edge cases, worked oracles {2,2,4,6,6}→μ=4,σ=2 with x=10→3.0 / x=8→2.0-boundary / x=−1→−2.5 and signature {3,1,−1,1}→2.0, research-grade scope, relation to the ssGSEA immune/stromal layer of [[immune-infiltration-deconvolution]]). Enriched with a genuinely-distinct "Implementation surface" section the Evidence page lacked: the three static `OncologyAnalyzer` entry points with signatures + defaults (CalculateExpressionZScore O(n); IdentifyOutlierGenes(sample, referenceCohorts, threshold=2.0) O(g·n) returning ExpressionOutlier{Gene,ZScore,Direction} in sample-dictionary iteration order; CalculateSignatureScore O(k)); the full validation/exception contract (null→ArgumentNullException; reference n<2 or σ=0, missing-cohort gene, empty signature→ArgumentException; non-positive threshold→ArgumentOutOfRangeException); case-sensitive dict keys + suffix-tree-N/A note; and the Framework no-bundled-cohort/signature/gene-set designation. Added the spec to sources: (now spec + 2 Evidence), bumped source_commit→c95a2e6e, updated→2026-07-15. Backlog: moved the row to [[backlog]] Covered-via-concept → [[expression-outlier-zscore-signature-score]] (covered 143→144, pending 101→100); removed from [[backlog-pending]] (Oncology 7→6, pending-total 103→102). No new page ⇒ index.md unchanged (concept already listed). No new graph nodes/edges (implementation prose only; concept's existing relates_to test-unit-registry edge and body wikilinks unchanged) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No contradictions — cBioPortal z-score spec + NormalizeExpressionLevels.java (n−1 SD, σ=0 fatal), cBioPortal FAQ (strict ±2), Lee et al. 2008 and GSVA/Hänzelmann 2013 (combined-z √k) cover complementary parts and agree with the recorded Evidence synthesis on the z-score, the strict outlier rule, and the √k signature denominator.
- 2026-07-15 — ingest docs/algorithms/Oncology/Tumor_Heterogeneity_Analysis.md (Tumor Heterogeneity Analysis SPEC — primary per-algorithm spec for `OncologyAnalyzer.CalculateITH`/`InferSubclones`/`AnalyzeHeterogeneity`, unit ONCO-HETERO-001, status Production; scalar ITH-summary layer — MATH = 100·1.4826·median(|f−median|)/median(f) per Mroz&Rocco 2013 / Mroz 2015 / maftools mathScore.R, Shannon H=−Σpᵢln pᵢ natural-log clonal diversity per Liu 2017/Shannon 1948, subclone count = occupied CCF clusters, subclonal fraction = #(CCF<0.95)/n with the Landau 2013 threshold). CONTEXT/decision per Oncology reconciliation brief: REUSED the existing concept [[intratumor-heterogeneity-metrics]] (created 2026-07-10 from the ONCO-HETERO-001 Evidence artifact) rather than create a new page or a wiki/sources/ page (spec != Evidence report). The concept already covered the full model (both MATH forms with the 1.4826 MAD-consistency factor + byte-for-byte maftools mapping, odd/even worked oracles 49.42 / 59.304, natural-log Shannon with {1.0}→0 / {0.5,0.5}→0.6931 / four-equal→ln4 oracles and the alpha-diversity cross-link, richness + 0.95 subclonal cut, zero-median/all-identical/single-clone corner cases, the two source-consistent assumptions cluster-size pᵢ + R even-count median, downstream relation to ONCO-CCF-001 clustering and ONCO-CLONAL-001 threshold). Enriched with a genuinely-distinct "Implementation surface (Seqeron.Genomics.Oncology)" section the Evidence page lacked: the three static OncologyAnalyzer entry points with signatures + complexity (CalculateITH O(n log n)/O(n) two median sorts with input-cloning non-mutation guarantee; InferSubclones O(n)/O(k) hash-set of labels; AnalyzeHeterogeneity O(n log n) returning HeterogeneityResult{MathScore,ShannonDiversity,SubcloneCount,SubclonalFraction} reusing ClusterCcfValues + ClonalCcfThreshold), the full validation/exception contract (null→ArgumentNullException; empty/non-finite/out-of-[0,1]/mismatched-length/zero-median→ArgumentException; clusterCount∉[1,count]→ArgumentOutOfRangeException), the suffix-tree-N/A note, and the not-implemented probabilistic PyClone/SciClone posterior-clustering boundary (deterministic k-means only). Added the spec to sources: (now Evidence + spec), bumped source_commit→e10eb245, updated→2026-07-15. Backlog: moved the row to [[backlog]] Covered-via-concept → [[intratumor-heterogeneity-metrics]] (covered 144→145); removed from [[backlog-pending]] (Oncology 6→5); also reconciled the pre-existing pending-total header drift to the true row count 98 in both [[backlog]] and [[backlog-pending]]. No new page ⇒ index.md unchanged (concept already listed in indexes/concepts.md). No new graph nodes/edges (implementation prose only; concept's existing depends_on ONCO-CCF-001 / relates_to ONCO-CLONAL-001 + test-unit-registry edges and body wikilinks unchanged) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No contradictions — Mroz&Rocco 2013, Mroz 2015, maftools mathScore.R, Liu 2017, Shannon 1948 and the reused Landau 2013 threshold cover complementary facets and agree with the recorded Evidence synthesis; the spec's formulas and 49.42/59.304 oracles match the concept verbatim.
- 2026-07-15 — ingest docs/algorithms/Oncology/Tumor_Mutational_Burden.md (Tumor Mutational Burden SPEC — primary per-algorithm spec for `OncologyAnalyzer.CalculateTMB(int,double)` / `CalculateTMB(IEnumerable<SomaticCall>,double)` / `ClassifyTMB(double)` + `TmbStatus` enum, unit ONCO-TMB-001, status Production; deterministic exact ratio TMB = mutationCount / targetRegionMb (mut/Mb) with FDA inclusive TMB-High ⇔ TMB ≥ 10 cutoff, 1.1 Mb FoundationOne 315-gene example denominator, filtering delegated upstream to ONCO-SOMATIC-001). CONTEXT/decision per Oncology reconciliation brief: REUSED the existing concept [[tumor-mutational-burden]] (created 2026-07-10 from the ONCO-TMB-001 Evidence artifact) rather than create a new page or a wiki/sources/ page (spec != Evidence report). The concept already covered the full model (mut/Mb formula with assay-footprint denominator incl. 1.1 Mb F1 / ~0.8 Mb F1CDx / ~30-40 Mb WES, synonymous-vs-nonsynonymous counting + pre-count driver/germline filtering per Chalmers 2017, FDA ≥10 inclusive TMB-High cutoff per Marcus 2021 cross-checked to Merino 2020 harmonization, the two-tier single-cutoff ASSUMPTION rejecting the unsourced 6/20 three-tier scheme, six worked oracles incl. the 11/1.1=10.0 boundary and 99/10=9.9 not-high, region=0-throws / non-negative-inputs / <0.5 Mb small-denominator-instability / monotonicity corner cases, relation to MSI + somatic caller + neoantigen/clinical-interpretation layers, research-grade/not-for-clinical scope). Enriched with a genuinely-distinct "Implementation surface (Seqeron.Genomics.Oncology)" section the Evidence page lacked: the three static OncologyAnalyzer entry points with signatures + complexity (CalculateTMB(int,double) O(1) division; the CalculateTMB(IEnumerable<SomaticCall>,double) overload counting only Somatic-status calls in one O(n) pass then delegating; ClassifyTMB O(1) returning the two-value TmbStatus enum), the full validation/exception contract (mutationCount<0 or targetRegionMb NaN/∞/≤0 → ArgumentOutOfRangeException with region=0 undefined-not-∞; null collection → ArgumentNullException then delegate; negative/non-finite tmb → ArgumentOutOfRangeException), and the suffix-tree-N/A + no-duplicated-filtering note. Added the spec to sources: (now spec + Evidence), bumped source_commit→eaf0bb76, updated→2026-07-15. Backlog: moved the row to [[backlog]] Covered-via-concept → [[tumor-mutational-burden]] (covered 145→146, pending 98→97); removed from [[backlog-pending]] (Oncology 5→4, pending-total 98→97). No new page ⇒ index.md unchanged (concept already listed in indexes/concepts.md). No new graph nodes/edges (implementation prose only; concept's existing relates_to test-unit-registry / somatic-variant-calling / microsatellite-instability edges and body wikilinks unchanged) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No contradictions — Chalmers 2017 (mut/Mb + synonymous counting + <0.5 Mb instability), Marcus 2021 (FDA ≥10 F1CDx cutoff) and Merino 2020 (harmonized unit/basis) cover complementary facets and agree with the recorded Evidence synthesis; the spec's ratio + inclusive ≥10 rule and 11/1.1=10.0 oracle match the concept verbatim.
- 2026-07-15 — ingest docs/algorithms/Oncology/Tumor_Phylogeny_Reconstruction.md (Tumor Phylogeny Reconstruction SPEC — primary per-algorithm spec for `OncologyAnalyzer.ReconstructPhylogeny(IReadOnlyList<CcfCluster>, double)` / `IdentifyTrunkMutations` / `IdentifyBranchMutations`, unit ONCO-PHYLO-001, status Simplified; deterministic constraint-satisfaction / perfect-phylogeny clonal-tree builder over CCF clusters — LICHeE Eq.2 lineage precedence u.CCFᵢ≥v.CCFᵢ−ε + presence pattern u=0⇒v=0 and Eq.5 per-node/per-sample sum rule Σchildren≤parent+ε, deepest-valid-ancestor / ascending-id deterministic tie-break, default ε=0). CONTEXT/decision per Oncology reconciliation brief: REUSED the existing concept [[tumor-phylogeny-clonal-tree-reconstruction]] (created 2026-07-10 from the ONCO-PHYLO-001 Evidence artifact) rather than create a new page or a wiki/sources/ page (spec != Evidence report). The concept already covered the full model (both ordering constraints in CCF form with LICHeE/PICTograph provenance, presence-pattern nesting rule, sum-rule branching-vs-nesting mechanics, trunk-vs-branch structural clonal/subclonal distinction, deepest-valid-ancestor tie-break + ε=0 as two source-consistent assumptions, three hand-derived worked datasets incl. the sum-rule chain-not-siblings case, the two ancestor≥descendant / per-node-sum invariants as test oracles + empty/single boundary cases, upstream ONCO-CCF-001 / sibling ONCO-CLONAL-001 / downstream ONCO-HETERO-001 relations, distance-based NJ/UPGMA out-of-scope contrast, research-grade/not-for-clinical scope). Enriched with a genuinely-distinct "Implementation surface (ONCO-PHYLO-001 spec)" section the Evidence page lacked: the three OncologyAnalyzer entry points with signatures + return shape (ReconstructPhylogeny→ClonalPhylogeny{RootId synthetic-normal-below-min-id CCF=1, Clusters input-order, Edges one ClonalEdge/non-root, SampleCount}; IdentifyTrunk/BranchMutations partitioning the clusters), the descending-total-CCF processing order + per-sample budget-debit construction step, the full validation/exception contract (same-nonzero-length CCF in [0,1] + unique ids, tolerance≥0 not-NaN; null→ArgumentNullException, ragged/NaN/out-of-range/dup-id→ArgumentException, negative/NaN tolerance→ArgumentOutOfRangeException, empty→root-only tree), and the O(n²·k) time / O(n·k) space complexity + suffix-tree-N/A (numeric constraint satisfaction, not substring search) note. Added the spec to sources: (now Evidence + spec), bumped source_commit→abca521a, updated→2026-07-15. Backlog: moved the row to [[backlog]] Covered-via-concept → [[tumor-phylogeny-clonal-tree-reconstruction]] (covered 146→147, pending 97→96); removed from [[backlog-pending]] (Oncology 4→3). No new page ⇒ index.md unchanged (concept already listed in indexes/concepts.md). No new graph nodes/edges (implementation prose only; concept's existing depends_on ONCO-CCF-001 / relates_to ONCO-CLONAL-001 + test-unit-registry edges and body wikilinks unchanged) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No contradictions — Popic 2015 (LICHeE, Genome Biology) and Zheng 2022 (PICTograph, Bioinformatics) state the same two constraints in VAF/CCF form and agree with the recorded Evidence synthesis; the spec's edge rule, sum rule, and worked walk-through match the concept verbatim.
- 2026-07-15 — ingest docs/algorithms/Oncology/Tumor_Ploidy_Estimation.md (Tumor Ploidy Estimation + Whole-Genome-Doubling SPEC — primary per-algorithm spec for `OncologyAnalyzer.EstimatePloidy` / `DetectWholeGenomeDoubling(segments, ReferenceGenome=GRCh38)` / `DetectWholeGenomeDoublingFromSuppliedLength` / `GetAutosomeLengths` / `GetAutosomalGenomeLength`, unit ONCO-PLOIDY-001, status Production; two post-segmentation genome-state summaries — average ploidy ψ = Σ(CN·L)/ΣL length-weighted mean total CN per Patchwork/Van Loo 2010 n-scale, and the facets-suite `is_genome_doubled` WGD rule frac_elevated_mcn = Σlen[mcn≥2 ∧ chr∈1..22]/autosomal_genome > 0.5 strict, mcn=tcn−lcn, reference-genome denominator 2,875,001,522 bp GRCh38 / 2,881,033,286 bp GRCh37). CONTEXT/decision per Oncology reconciliation brief: REUSED the existing concept [[tumor-ploidy-estimation-and-whole-genome-doubling]] (created 2026-07-10 from the ONCO-PLOIDY-001 Evidence artifact) rather than create a new page or a wiki/sources/ page (spec != Evidence report). The concept already covered the full model (Patchwork length-weighted ψ with the 3.0 worked oracle, n-scale 2n=diploid + >2.7n aneuploidy, the three load-bearing WGD details major-not-total CN / reference-autosome denominator / strict >0.5 with GRCh38 half=1,437,500,761 bp boundary, the legacy supplied-length overload, empty/invalid/boundary/all-1:1/chrXY/2:0-LOH corner cases, UCSC+Ensembl provenance, research-grade/not-for-clinical scope, and the distinct-from-ASCAT-joint-fit framing vs [[allele-specific-copy-number-ascat]] + purity counterpart [[tumor-purity-from-mutation-vaf]]). Enriched with a genuinely-distinct "Implementation surface (ONCO-PLOIDY-001 spec)" section the Evidence page lacked: the four OncologyAnalyzer entry points with signatures + defaults, single-pass O(n)/O(1) complexity + suffix-tree-N/A note, the shared AlleleSpecificSegment record + ValidateSegment helper reused from ONCO-LOH-001 / ONCO-HRD-001, the full validation/exception contract (null→ArgumentNullException; End≤Start or negative CN→ArgumentException; EstimatePloidy empty→ArgumentException vs DetectWholeGenomeDoubling empty→false; undefined ReferenceGenome→ArgumentOutOfRangeException), the INV-03 bounding property min CN≤ψ≤max CN, and the accepted registry deviation (scalar DetectWholeGenomeDoubling(ploidy) stub vs canonical per-segment method). Added the spec to sources: (now spec + Evidence), bumped source_commit→b6006db0, updated→2026-07-15. Backlog: moved the row to [[backlog]] Covered-via-concept → [[tumor-ploidy-estimation-and-whole-genome-doubling]] (covered 147→148, pending 96→95); removed from [[backlog-pending]] (Oncology 3→2, pending-total 96→95). No new page ⇒ index.md unchanged (concept already listed in indexes/concepts.md). No new graph nodes/edges (implementation prose only; concept's existing relates_to test-unit-registry / allele-specific-copy-number-ascat edges and body wikilinks unchanged) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No contradictions — Patchwork (length-weighted ψ), Van Loo 2010 (n-scale + >2.7n aneuploidy), Bielski 2018 + facets-suite is_genome_doubled (major-CN≥2 / >50% WGD rule) and UCSC/Ensembl (reference autosome lengths) each cover a disjoint facet and agree with the recorded Evidence synthesis; the spec's ψ=3.0 oracle and half=1,437,500,761 bp WGD boundary match the concept verbatim.
- 2026-07-15 — ingest docs/algorithms/Oncology/Tumor_Purity_Estimation.md (Tumor Purity Estimation SPEC — primary per-algorithm spec for `OncologyAnalyzer.EstimatePurityFromVaf(double)` / `EstimatePurityFromVAF(IEnumerable<VariantObservation>)` / `EstimatePurity(IEnumerable<PurityVariant>)`, unit ONCO-PURITY-001, status Production; closed-form CNAqc expected-VAF inversion — diploid-het special case ρ=2·VAF and allele-specific general inversion π=2v/[m+v(2−n_tot)] over the mixture denominator 2(1−π)+π·n_tot, median cross-variant aggregation, normal fixed at 2 copies per CNAqc/FACETS/ABSOLUTE). CONTEXT/decision per Oncology reconciliation brief: REUSED the existing concept [[tumor-purity-from-mutation-vaf]] (created 2026-07-10 from the ONCO-PURITY-001 Evidence artifact) rather than create a new page or a wiki/sources/ page (spec != Evidence report). The concept already covered the full model (CNAqc expected-VAF generative formula v=m·π·c/[2(1−π)+π·n_tot] with FACETS m*=mΦ+(1−Φ) and ABSOLUTE corroboration, the π=2v/(m·c+2v−v·n_tot) inversion + ρ=2·VAF diploid special case, the CNAqc worked oracles 60%→30% VAF and the 2:1 33%/66% two-peak multiplicity ambiguity, median aggregation, subclonal/low-purity/no-informative-variant corner cases, distinct-from-ASCAT-joint-fit framing vs [[allele-specific-copy-number-ascat]] and ploidy counterpart [[tumor-ploidy-estimation-and-whole-genome-doubling]], research-grade/not-for-clinical scope). Enriched with a genuinely-distinct "Implementation surface (ONCO-PURITY-001 spec)" section the Evidence page lacked: the three OncologyAnalyzer overloads with signatures — incl. the single-scalar EstimatePurityFromVaf(double) the concept had not named — the PurityVariant(vaf, Multiplicity, TumorTotalCopyNumber) record, the O(n log n)/O(n) collection complexity (median sort, O(1) per variant) vs O(1) scalar, the ONCO-VAF-001 CalculateVAF read-count delegation, and the full validation/exception contract (null→ArgumentNullException; empty→ArgumentException; VAF∉[0,1], diploid VAF>0.5, m<1, n_tot<1, or ρ∉[0,1]/non-positive denominator→ArgumentOutOfRangeException). Added the spec to sources: (now Evidence + spec), bumped source_commit→c031ee74, updated→2026-07-15. Backlog: moved the row to [[backlog]] Covered-via-concept → [[tumor-purity-from-mutation-vaf]] (covered 148→149, pending 95→94); removed from [[backlog-pending]] (Oncology 2→1, pending-total 95→94). No new page ⇒ index.md unchanged (concept already listed in indexes/concepts.md). No new graph nodes/edges (implementation prose only; concept's existing relates_to test-unit-registry / alternative_to allele-specific-copy-number-ascat edges and body wikilinks unchanged) ⇒ graph lint/extract skipped. Hub [[algorithm-validation-evidence]] unchanged (spec docs are not added to the hub). No contradictions — CNAqc (expected-VAF formula + worked peaks), FACETS (mixing-model denominator), ABSOLUTE (inverse purity/copy-number correction) cover complementary facets and agree with the recorded Evidence synthesis; the spec's ρ=2·VAF, allele-specific inversion, and 60%→30% / 2:1 π=1 oracles match the concept verbatim.
