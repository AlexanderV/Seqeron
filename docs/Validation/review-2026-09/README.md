# Full Algorithm Code Review — Campaign 2026-09

**Branch:** `claude/stoic-maxwell-0olr2z` (all sessions push here; never to another branch).
**Scope:** every implemented algorithm (all ☑ units in `ALGORITHMS_CHECKLIST_V2.md`, the
extra units listed in `docs/Validation/VALIDATION_LEDGER.md` → *New units*, and any public
algorithmic method in a batch's owned files that no unit covers).
**Method:** `docs/Validation/VALIDATION_PROTOCOL.md` (Stage A description → Stage B
implementation), plus the three campaign rules below.

**Baseline (2026-09-28, orchestrator):** `dotnet build Seqeron.sln` 0 errors; `Seqeron.Genomics.Tests` 20267 passed / 0 failed.

## Campaign rules (non-negotiable)

1. **Sourced, confirmed fixes only.** Every change must be justified by an authoritative
   source (primary paper, standard/spec, reference textbook) **and** confirmed numerically
   against a reference implementation or dataset (Biopython, ViennaRNA, primer3-py,
   scikit-bio, scipy, the tool's original published code / worked example, a published table).
   Tests lock the *sourced* values; tests are never bent to match code.
2. **No simplifications unless a full implementation is impossible or critically complex.**
   If a method is a heuristic/approximation of a published algorithm and there is *any*
   realistic way to implement the real algorithm (published model parameters, tables, or a
   reference implementation are obtainable), implement the real algorithm. Only when it is
   genuinely impossible (proprietary/trained model not obtainable, needs a multi-GB database,
   etc.) keep it — and then it must be declared honestly (XML doc + `docs/Validation/LIMITATIONS.md`
   entry proposal in the batch report) with the exact reason.
3. **No duplication.** If an algorithm depends on another algorithm, it must *call* the
   canonical implementation, not re-implement it (e.g. GC content, reverse complement,
   Hamming/edit distance, Tm, translation, Shannon entropy, k-mer counting, alignment, stats
   helpers like mean/variance/p-values in `Seqeron.Genomics.Infrastructure/StatisticsHelper.cs`).
   Also check the MCP tool wrappers (`src/Seqeron/Mcp/**`) that expose your classes: a wrapper
   must delegate, not re-implement.

## Environment for a fresh container

```bash
apt-get update -qq; DEBIAN_FRONTEND=noninteractive apt-get install -y -qq dotnet-sdk-10.0
pip install -q biopython primer3-py ViennaRNA scikit-bio scipy numpy   # reference implementations
export DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
dotnet build Seqeron.sln -c Debug            # must stay 0 errors
dotnet test tests/Seqeron/Seqeron.Genomics.Tests/Seqeron.Genomics.Tests.csproj -c Debug --no-build --filter FullyQualifiedName~<Class>
```

Network: `en.wikipedia.org`, `ncbi.nlm.nih.gov`, `europepmc`, `rosalind.info` are blocked for
curl/WebFetch. Use `WebSearch` (snippets), PyPI packages, `raw.githubusercontent.com`
(reference source code of original tools), and publisher sites that are reachable.
Record for each source what was actually opened.

## Ownership & concurrency

Batches run concurrently in separate sessions and push to the same branch.

- Edit **only the files your batch owns** (table below) + their tests/TestSpecs/Evidence/algorithm docs.
- Shared files (`Seqeron.Genomics.Infrastructure/*`, `Seqeron.Genomics.Core/LimitationPolicy.cs`,
  `OncologyAnalyzer.cs`, MCP servers): **additive changes only** (new members), no signature
  changes of existing public members. If a cross-batch change is needed (e.g. a duplicate whose
  canonical lives in another batch's file and must change), do **not** do it — record it in your
  report under *Cross-batch dedup requests*; the dedup phase handles it.
- Removing a duplicate by making **your own** file call an existing canonical method elsewhere is allowed.
- **Do not edit** `VALIDATION_LEDGER.md` / `FINDINGS_REGISTER.md` / `ALGORITHMS_CHECKLIST_V2.md`
  (the orchestrator consolidates). Write your results only to `docs/Validation/review-2026-09/<BATCH>.md`.
- Commit per unit (or per fix), then `git fetch origin claude/stoic-maxwell-0olr2z && git rebase origin/claude/stoic-maxwell-0olr2z && git push -u origin claude/stoic-maxwell-0olr2z`; on rejection repeat the fetch/rebase/push (network errors: retry with backoff 2/4/8/16 s).
- **Two test tiers** (measured 2026-09-28: `Fuzzing` 50 %, `Properties` 29 %, `Metamorphic` 8 % of test CPU time):
  - **Fast tier — before every per-unit push:** `dotnet build Seqeron.sln` 0 errors, then
    `dotnet test tests/Seqeron/Seqeron.Genomics.Tests/Seqeron.Genomics.Tests.csproj -c Debug --no-build --filter "FullyQualifiedName!~.Fuzzing.&FullyQualifiedName!~.Properties.&FullyQualifiedName!~.Metamorphic."`
    must be green (all Unit/Combinatorial/Algebraic/Snapshot/Mutation/Differential/Architecture tests, not only your class).
    Never disable, skip or `[Ignore]` a test — the slow tiers are deferred by the filter only.
  - **Heavy tier — once, when the whole batch is done (mandatory, before the final push):**
    1. Find every test in `Fuzzing/`, `Properties/`, `Metamorphic/` that exercises your owned classes/methods.
    2. **Update** those whose expectations the batch's fixes legitimately changed (with the sourced justification in the report) — never weaken an invariant to make it pass.
    3. **Add** new fuzz/property/metamorphic tests for new behaviour/invariants introduced by the batch (e.g. new algorithm paths, new edge-case contracts).
    4. Run the **full** `Seqeron.Genomics.Tests` (no filter) + the MCP test project of your server; all green.
    5. Record in the report a *Heavy-tier* section: tests updated / added / run counts.

## Session management inside a batch

Batch sessions must not hold every unit in one context. The batch lead:
1. Reads this README and the batch row, lists the units + owned files.
2. For **each unit**, spawns **one fresh subagent** (Agent tool) that does Stage A + Stage B,
   implements the fixes + tests, builds, runs the unit tests, and returns a short structured
   summary (verdicts, findings, fixes, files changed, sources opened, cross-batch requests).
   Run unit subagents **sequentially** (one build directory).
3. The lead commits + pushes after each unit, appends the summary to the batch report, and keeps
   only the summaries in its own context.
4. After all units **and after the completeness-audit loop below**: one extra fresh subagent sweeps the owned files for **duplication**
   (internal and vs. other canonical implementations) and for public algorithmic methods not
   covered by any unit; fixes what is in-ownership; records the rest.
5. Final full test run, final push, report finished.

## Definition of Done — no doable leftovers (mandatory)

A batch is **not finished** while anything doable remains. "Deferred", "follow-up", "optional",
"out of scope for this pass", "left for later", "could be implemented", "TODO" are **not** valid
end states. The only accepted leftover is one that is *impossible or critically complex*, and it
must carry concrete proof (what exactly is missing: a proprietary model/dataset, an unavailable
upstream binary, a redesign spanning files owned by another batch — with the exact blocker).
"Takes time", "large", "many tests to update", "beyond this unit's scope" are **not** proof.

Before writing the final report the lead runs a **completeness-audit loop**:
1. Spawn a fresh **auditor** subagent. It reads the batch report, `git diff <batch-start>..HEAD`
   of the owned files, and greps the owned code + report for: `TODO`, `FIXME`, `deferred`,
   `follow-up`, `optional`, `not implemented`, `NotImplemented`, `simplif`, `approximat`,
   `out of scope`, `LIMITED`, `left for`, `future`, `could `, `should `, `not yet`; it also greps
   every other `B*.md` report for cross-batch requests addressed to this batch, and checks every
   unit's Stage A/B findings have a matching fix. It returns a list of leftovers, each classified
   **DOABLE** or **BLOCKED (with proof)**. The auditor is adversarial: when in doubt → DOABLE.
2. The lead implements every DOABLE item (fresh subagent per item/unit, same rules: sourced,
   reference-checked, tests, fast tier, commit + push).
3. Repeat 1–2 until the auditor returns **zero DOABLE items**. Then run the heavy tier.
   **Order is mandatory: units → audit/finish loop → duplication sweep → final audit → heavy tier (once).**
   Leftovers are finished *before* the duplication sweep, so the sweep also covers the code the
   finish loop added (no second sweep). The sweep's own leftovers ("duplication kept because…")
   go through one final audit round (fix every DOABLE item, fast tier) before the heavy tier.
   The heavy tier (Fuzzing/Properties/Metamorphic: update + add + full run) runs **exactly once, at
   the very end**, after the last DOABLE item is done — never between audit rounds, never per
   leftover item. During the loop (and in finisher sessions) use only the fast tier per change.
   If the heavy tier itself exposes a defect, fix it (fast tier while fixing) and re-run only the
   affected heavy-tier classes, then one final full heavy run.
4. The final report and the final message end with a `## Leftovers` section: either
   `none`, or each BLOCKED item with its proof. Cross-batch requests *to other batches* are
   listed there as well (they are not leftovers of this batch only if the target file is owned
   by another batch).

**"LIMITATIONS proposals" / "Simplifications kept" / "Limitations proposed (not bugs)" are leftovers too.**
A documented limitation is not a finish line: if the fuller behaviour exists in a source or a
reference tool, implement it (as the default when the reference does so, otherwise as an extra
option/overload/parameter, keeping backward compatibility). Examples that are DOABLE, not
limitations: emitting the trailing partial window like Biopython; returning *all* extrema like
Rosalind BA1F; a thermodynamic (Turner) score instead of a pair count ("slow" → optimise or make
it opt-in); a mismatch/gap-tolerant variant like EMBOSS einverted ("a different algorithm" → add it
as a variant); a missing genetic-code parameter on an MCP tool ("signature change" → add an optional
parameter); an unsourced convenience metric → find the source or replace/rename it per a sourced
definition. Only items with a real BLOCKED proof may remain in the Limitations section, and each
must say why it is impossible or critically complex.

Never end the session with a question to the user or a list of "possible next steps": no human
is watching. If something is doable, do it.

## Batch report format (`docs/Validation/review-2026-09/<BATCH>.md`)

```markdown
# <BATCH> — <areas>
| Unit | Stage A | Stage B | State | Fix summary |
|---|---|---|---|---|
## Confirmed fixes
### F<n> — <unit>: <title>
- Defect (with minimal repro), source (what was opened), reference cross-check (numbers), fix, tests added, commit.
## Duplication removed
## Simplifications kept (with reason) / LIMITATIONS proposals
## Cross-batch dedup requests
## Test results (counts before/after)
```

## Batches

| Batch | Units | Owned files |
|---|---|---|
| B01 | SEQ-GC-001 SEQ-COMP-001 SEQ-REVCOMP-001 SEQ-VALID-001 SEQ-COMPLEX-001 SEQ-ENTROPY-001 SEQ-GCSKEW-001 SEQ-ATSKEW-001 SEQ-GC-ANALYSIS-001 SEQ-REPLICATION-001 SEQ-RNACOMP-001 | Core/SequenceExtensions.cs, Core/DnaSequence.cs, Core/RnaSequence.cs, Core/ISequence.cs, Core/IupacHelper.cs, Analysis/GcSkewCalculator.cs |
| B02 | TRANS-CODON-001 TRANS-PROT-001 TRANS-SIXFRAME-001 CODON-CAI-001 CODON-ENC-001 CODON-OPT-001 CODON-RARE-001 CODON-RSCU-001 CODON-STATS-001 CODON-USAGE-001 | Core/Translator.cs, Core/GeneticCode.cs, Core/ProteinSequence.cs, MolTools/CodonOptimizer.cs, MolTools/CodonUsageAnalyzer.cs |
| B03 | SEQ-STATS-001 SEQ-COMPOSITION-001 SEQ-MW-001 SEQ-PI-001 SEQ-HYDRO-001 SEQ-THERMO-001 SEQ-TM-001 SEQ-DINUC-001 SEQ-SECSTRUCT-001 SEQ-CODON-FREQ-001 SEQ-ENTROPY-PROFILE-001 SEQ-GC-PROFILE-001 SEQ-SUMMARY-001 | Analysis/SequenceStatistics.cs |
| B04 | SEQ-COMPLEX-COMPRESS-001 SEQ-COMPLEX-DUST-001 SEQ-COMPLEX-KMER-001 SEQ-COMPLEX-WINDOW-001 REP-STR-001 REP-TANDEM-001 REP-INV-001 REP-DIRECT-001 REP-PALIN-001 REP-APPROX-001 | Analysis/SequenceComplexity.cs, Analysis/RepeatFinder.cs |
| B05 | PAT-EXACT-001 PAT-APPROX-001 PAT-APPROX-002 PAT-APPROX-003 PAT-IUPAC-001 PAT-PWM-001 MOTIF-CONS-001 MOTIF-DISCOVER-001 MOTIF-SHARED-001 MOTIF-REGULATORY-001 MOTIF-GENERATE-001 | Analysis/MotifFinder.cs, Alignment/ApproximateMatcher.cs, src/SuffixTree/** |
| B06 | KMER-COUNT-001 KMER-FREQ-001 KMER-FIND-001 KMER-ASYNC-001 KMER-BOTH-001 KMER-DIST-001 KMER-GENERATE-001 KMER-POSITIONS-001 KMER-STATS-001 KMER-UNIQUE-001 | Analysis/KmerAnalyzer.cs |
| B07 | PRIMER-TM-001 PRIMER-DESIGN-001 PRIMER-STRUCT-001 PRIMER-NNTM-001 PRIMER-HAIRPIN-001 PRIMER-DIMER-001 PROBE-DESIGN-001 PROBE-VALID-001 PROBE-LNATM-001 PROBE-EVALUE-001 | MolTools/PrimerDesigner.cs, MolTools/ProbeDesigner.cs, MolTools/NtthalDimer.cs, MolTools/NtthalHairpin.cs, Infrastructure/ThermoConstants.cs |
| B08 | CRISPR-PAM-001 CRISPR-GUIDE-001 CRISPR-OFF-001 RESTR-FIND-001 RESTR-DIGEST-001 RESTR-FILTER-001 | MolTools/CrisprDesigner.cs, MolTools/AzimuthRuleSet2.cs, MolTools/RestrictionAnalyzer.cs |
| B09 | ALIGN-GLOBAL-001 ALIGN-LOCAL-001 ALIGN-SEMI-001 ALIGN-MULTI-001 ALIGN-STATS-001 GENOMIC-COMMON-001 GENOMIC-MOTIFS-001 GENOMIC-ORF-001 GENOMIC-REPEAT-001 GENOMIC-SIMILARITY-001 GENOMIC-TANDEM-001 | Alignment/SequenceAligner.cs, Alignment/AnchorBasedAligner.cs, Infrastructure/AlignmentTypes.cs, Analysis/GenomicAnalyzer.cs |
| B10 | ASSEMBLY-OLC-001 ASSEMBLY-DBG-001 ASSEMBLY-STATS-001 ASSEMBLY-MERGE-001 ASSEMBLY-SCAFFOLD-001 ASSEMBLY-COVER-001 ASSEMBLY-CONSENSUS-001 ASSEMBLY-TRIM-001 ASSEMBLY-CORRECT-001 | Alignment/SequenceAssembler.cs, Chromosome/GenomeAssemblyAnalyzer.cs |
| B11 | ANNOT-ORF-001 ANNOT-GENE-001 ANNOT-PROM-001 ANNOT-GFF-001 ANNOT-CODING-001 ANNOT-REPEAT-001 ANNOT-CODONUSAGE-001 SPLICE-DONOR-001 SPLICE-ACCEPTOR-001 SPLICE-PREDICT-001 SPLICE-MAXENT3-001 SPLICE-MAXENT5-001 | Annotation/GenomeAnnotator.cs, Annotation/SpliceSitePredictor.cs |
| B12 | RNA-STRUCT-001 RNA-STEMLOOP-001 RNA-ENERGY-001 RNA-PAIR-001 RNA-HAIRPIN-001 RNA-MFE-001 RNA-PSEUDOKNOT-001 RNA-DOTBRACKET-001 RNA-INVERT-001 RNA-PARTITION-001 RNA-ACCESS-001 RNA-PKPREDICT-001 RNA-PKRECURSIVE-001 | Analysis/RnaSecondaryStructure.cs, Core/Turner2004Parameters.cs |
| B13 | MIRNA-SEED-001 MIRNA-TARGET-001 MIRNA-PRECURSOR-001 MIRNA-PAIR-001 MIRNA-CONTEXT-001 MIRNA-PCT-001 MIRNA-CLASSIFY-001 MIRNA-CLEAVAGE-001 | Annotation/MiRnaAnalyzer.cs |
| B14 | PROTMOTIF-FIND-001 PROTMOTIF-PROSITE-001 PROTMOTIF-DOMAIN-001 PROTMOTIF-PATTERN-001 PROTMOTIF-SP-001 PROTMOTIF-TM-001 PROTMOTIF-CC-001 PROTMOTIF-LC-001 PROTMOTIF-COMMON-001 PROTMOTIF-HMM-001 | Analysis/ProteinMotifFinder.cs, Analysis/Plan7ProfileHmm.cs |
| B15 | DISORDER-PRED-001 DISORDER-REGION-001 DISORDER-MORF-001 DISORDER-PROPENSITY-001 DISORDER-LC-001 EPIGEN-CPG-001 EPIGEN-METHYL-001 EPIGEN-DMR-001 EPIGEN-BISULF-001 EPIGEN-CHROM-001 EPIGEN-AGE-001 | Analysis/DisorderPredictor.cs, Annotation/EpigeneticsAnalyzer*.cs |
| B16 | PHYLO-DIST-001 PHYLO-TREE-001 PHYLO-NEWICK-001 PHYLO-COMP-001 PHYLO-BOOT-001 PHYLO-STATS-001 POP-FREQ-001 POP-DIV-001 POP-HW-001 POP-FST-001 POP-LD-001 POP-SELECT-001 POP-ANCESTRY-001 POP-ROH-001 | Phylogenetics/PhylogeneticAnalyzer.cs, Population/PopulationGeneticsAnalyzer.cs |
| B17 | META-CLASS-001 META-PROF-001 META-ALPHA-001 META-BETA-001 META-BIN-001 META-FUNC-001 META-RESIST-001 META-PATHWAY-001 META-TAXA-001 META-CHECKM-001 META-TETRA-001 | Metagenomics/MetagenomicsAnalyzer.cs, Metagenomics/TaxonomyTree.cs |
| B18 | PANGEN-CORE-001 PANGEN-CLUSTER-001 PANGEN-HEAP-001 PANGEN-MARKER-001 COMPGEN-SYNTENY-001 COMPGEN-ORTHO-001 COMPGEN-REARR-001 COMPGEN-RBH-001 COMPGEN-COMPARE-001 COMPGEN-REVERSAL-001 COMPGEN-CLUSTER-001 COMPGEN-ANI-001 COMPGEN-DOTPLOT-001 | Metagenomics/PanGenomeAnalyzer.cs, Analysis/ComparativeGenomics.cs |
| B19 | CHROM-TELO-001 CHROM-CENT-001 CHROM-KARYO-001 CHROM-ANEU-001 CHROM-SYNT-001 CHROM-ALPHASAT-001 CHROM-HOR-001 TRANS-EXPR-001 TRANS-DIFF-001 TRANS-SPLICE-001 | Chromosome/ChromosomeAnalyzer.cs, Annotation/TranscriptomeAnalyzer.cs |
| B20 | PARSE-FASTA-001 PARSE-FASTQ-001 PARSE-BED-001 PARSE-VCF-001 PARSE-GFF-001 PARSE-GENBANK-001 PARSE-EMBL-001 QUALITY-PHRED-001 QUALITY-STATS-001 | Seqeron.Genomics.IO/** |
| B21 | VARIANT-CALL-001 VARIANT-SNP-001 VARIANT-INDEL-001 VARIANT-ANNOT-001 SV-DETECT-001 SV-BREAKPOINT-001 SV-CNV-001 | Annotation/VariantCaller.cs, Annotation/VariantAnnotator.cs, Annotation/StructuralVariantAnalyzer.cs |
| B22 | ONCO-SOMATIC-001 ONCO-VAF-001 ONCO-DRIVER-001 ONCO-ARTIFACT-001 ONCO-ANNOT-001 ONCO-TMB-001 ONCO-MSI-001 ONCO-HRD-001 ONCO-LOH-001 | Oncology/OncologyAnalyzer.SomaticCalling.cs, .DriversArtifactsAnnotation.cs, .TmbMsi.cs, .HrdLoh.cs |
| B23 | ONCO-SIG-001 ONCO-SIG-002 ONCO-SIG-003 ONCO-SIG-004 ONCO-FUSION-001 ONCO-FUSION-002 ONCO-FUSION-003 | Oncology/OncologyAnalyzer.Signatures.cs, .Fusions.cs |
| B24 | ONCO-CNA-001 ONCO-CNA-002 ONCO-CNA-003 ONCO-PURITY-001 ONCO-PLOIDY-001 ONCO-ASCAT-001 ONCO-CLONAL-001 ONCO-CCF-001 ONCO-PHYLO-001 ONCO-HETERO-001 | Oncology/OncologyAnalyzer.CopyNumberPloidy.cs, .Clonality.cs, .PhylogenyHeterogeneity.cs |
| B25 | ONCO-CTDNA-001 ONCO-MRD-001 ONCO-CHIP-001 ONCO-HLA-001 ONCO-ACTION-001 ONCO-SV-001 ONCO-EXPR-001 | Oncology/OncologyAnalyzer.CtdnaMrdChip.cs, .ClinicalMisc.cs |
| B26 | ONCO-NEO-001 ONCO-MHC-001 ONCO-IMMUNE-001 MHC-NN-001 MHC-MATRIX-001 IMMUNE-NUSVR-001 | Oncology/OncologyAnalyzer.Neoantigen.cs, Oncology/MhcflurryAffinityPredictor.cs, Oncology/ImmuneAnalyzer.cs |

Paths are relative to `src/Seqeron/Algorithms/Seqeron.Genomics.<Project>/`.

## Phase 2 (after all batches) — cross-batch deduplication

A dedicated session collects every *Cross-batch dedup request* from the batch reports, runs a
whole-repo duplication sweep, and removes each confirmed duplicate by routing callers to the
canonical implementation (behaviour-preserving, tests green). Then the orchestrator consolidates
the ledger, findings register and wiki.
