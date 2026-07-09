# Wiki Index

The catalog of all pages in this wiki. Each entry: a wikilink to the page and a one-line summary. The LLM reads this first when answering queries to identify candidate pages.

Keep summaries tight — one line each. The index is engineered to be cheap to read; a fat index defeats its purpose.

When this file exceeds ~300 lines or the wiki passes ~150 pages, shard into `wiki/indexes/<type>.md` and replace this file with a directory of shards. See the `scaling-playbook.md` reference in the `llm-wiki` skill for the migration procedure.

---

## Sources

- [[readme]] — the project's front-page `README.md`: what Seqeron is, its three entry points, and headline facts.
- [[algorithms-checklist-v2]] — the test-unit validation registry: 364 units (255 done / 109 proposed), DoD, evidence, coverage.
- [[advanced-testing-checklist]] — 2026-03-19 effectiveness/gap analysis rating the ten testing techniques by applicability, coverage, effort, and P0–P3 priority.
- [[mcp-methods-audit]] — a 2026-01-23 census of the C# API surface: 277 public static methods across 54 classes, with counting rules and the tool-count caveat.
- [[mcp-checklist]] — the SUPERSEDED (2026-07-01) MCP build tracker: a 12-server/241-tool plan never built; DoD gates, error-code catalog, two-tests-per-tool. Live status now in `docs/mcp/MCP_STATUS.md`.
- [[mcp-plan]] — the SUPERSEDED (2026-07-01) MCP design doc (sibling of the checklist): v3→v4 rationale, full 241-tool inventory across 12 servers, and the still-valid standards (naming, JSON-Schema 2020-12, error ranges, SemVer, 2-tests-per-tool, doc contract).
- [[mcp-prompt]] — the CURRENT one-tool-per-session MCP-completion subagent prompt: names the real 11 servers (server→project→file), the 3-part tool DoD (binding+tests+docs), and reverses the old no-business-asserts test policy. Ground truth: `docs/mcp/MCP_STATUS.md`.
- [[sonar-gate-plan]] — the Sonar-ratchet tracker: 66/66 SonarAnalyzer rules moved from report-only to blocking (or silenced-with-justification); Groups A–E, notable behaviour changes, and review-not-blind-fix on the S1244/S125 giants.
- [[align-global-001-evidence]] — evidence artifact for test unit ALIGN-GLOBAL-001 (Needleman–Wunsch): Wikipedia sources, the `GCATGCG`/`GATTACA` worked example, and the GapExtend-as-*d* / GapOpen-unused note.
- [[align-multi-001-evidence]] — evidence artifact for test unit ALIGN-MULTI-001 (MSA): main doc + MUSCLE and T-Coffee addenda; star/iterative/consistency aligners, SP-score and consensus design choices.
- [[align-semi-001-evidence]] — evidence artifact for test unit ALIGN-SEMI-001 (semi-global / fitting): Wikipedia + Rosalind SIMS/SMGB + Brudno 2003 sources, fitting-variant init/traceback, INV-1..5.
- [[align-stats-001-evidence]] — evidence artifact for test unit ALIGN-STATS-001 (identity/similarity/gaps + srspair formatting): EMBOSS needle/AlignFormats + BLAST + pseqsid sources, denominator-includes-gaps rule, three datasets.
- [[annot-coding-001-evidence]] — evidence artifact for test unit ANNOT-CODING-001 (CPAT hexamer coding-potential score): CPAT paper + FrameKmer.py + EMBOSS tcode sources, mean-ln-ratio formula, ±1 pseudo-scores, worked example 0.3466.
- [[annot-codonusage-001-evidence]] — evidence artifact for test unit ANNOT-CODONUSAGE-001 (RSCU): LIRMM + PMC2528880 + Sharp & Li 1986 + CodonU sources, n_i·x/Σx formula, Leu 3/1.5/1.5/0/0/0 example, pooled-counts + sense-codons-only, no CAI pseudocount.

## Concepts

- [[three-front-doors]] — one algorithm engine exposed through skills, the C# API, and MCP, with identical results.
- [[skill-layer]] — the Agent-Skill routing + discipline layer that keeps 427 tool schemas out of the model's context.
- [[layered-architecture]] — the strict up-only dependency layering (Levels 0–4) enforced by architecture tests.
- [[scientific-rigor]] — runtime honesty: `LimitationPolicy`, tool-only computation, provenance on every result.
- [[validation-and-testing]] — 22k+ tests across ten methodologies plus the per-unit validation campaign.
- [[build-quality-gate]] — the build-time static-analysis gate: SonarAnalyzer ratcheted to blocking under `TreatWarningsAsErrors`, fix-or-silence, review-not-blind-fix.
- [[test-unit-registry]] — the area-prefixed Test Unit ID scheme and per-unit record behind the validation effort.
- [[definition-of-done]] — the six-criterion acceptance bar (TestSpec, tests, ≥80% coverage, edge cases, CI, evidence) each unit must clear.
- [[algorithm-validation-evidence]] — the templated per-unit `docs/Evidence/` artifact pattern (sources, worked-example dataset, deviations) behind the "Evidence documented" DoD criterion.
- [[global-alignment-needleman-wunsch]] — the canonical DP global-alignment algorithm (`GlobalAlign`): linear-gap recurrence, O(nm), traceback, affine-as-extension.
- [[multiple-sequence-alignment]] — aligning 3+ sequences (NP-complete): Seqeron's star (`MultipleAlign`), iterative-refinement (MUSCLE), and consistency (T-Coffee) aligners; sum-of-pairs + consensus.
- [[semi-global-alignment-fitting]] — ends-free "glocal" alignment: fitting/query-in-reference variant (Rosalind SIMS); NW recurrence with free reference end gaps, traceback from max of last row.
- [[alignment-statistics]] — post-alignment metrics: percent identity/similarity/gaps (EMBOSS/BLAST convention, denominator includes gaps, positive-score ⇒ similar) and srspair three-line markup.
- [[coding-potential-hexamer-score]] — CPAT hexamer usage-bias coding-potential score: mean of ln(coding/noncoding) over in-frame hexamers (wordSize 6/step 3), sign convention, ±1 pseudo-scores; Fickett TESTCODE as not-implemented alternative.
- [[relative-synonymous-codon-usage]] — RSCU (Sharp & Li 1986) codon-usage-bias measure: n_i·x/Σx normalization, 1.0=no bias, bounded [0,n_i], Σ-over-family=n_i invariant; pooled counts, sense codons only; base anchor of the codon-usage family (vs CAI's 0.5 pseudocount).

## Gotchas

- [[research-grade-limitations]] — beta, not for clinical use; simplified-subset implementations; internal-only validation.

## Synthesis

(populated as query answers are filed back)
