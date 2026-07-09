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
