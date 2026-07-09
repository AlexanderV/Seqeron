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
