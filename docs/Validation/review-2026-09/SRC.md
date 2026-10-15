# SRC — sources & citations (2026-10-14/15)

Not a batch: no algorithm units. Scope: QUEUE "Consolidation task — re-anchor citations to reachable sources" and
"Network access opened … two BLOCKED items re-opened" (RO1, RO2). Network check at start: `web.expasy.org`,
`www.kazusa.or.jp`, `en.wikipedia.org`, `www.addgene.org` answer 200 (payload checked); Wikipedia's `api.php` throttles
(HTTP 429), so pages were fetched as raw wikitext (`index.php?action=raw`, 200).

## A. RO2 — ProtScale exponential weight-variation model (SRC for B03) — DONE

- Fetched `https://web.expasy.org/protscale/protscale-doc.html`: it defines the model only in words ("makes the weights
  decrease exponentially from the central position to the window edge") and its W 7 / edge 10 % example
  (10 %, 12 %, 39 %, 100 %) is **stale** — the live tool prints and uses 0.10, 0.13, 0.29, 1.00.
- Equation derived from the live tool (`cgi-bin/protscale/protscale.pl`; printed weight tables + verbose 3-dp scores,
  34 runs, W 3–21, edge 0–99.5 %): **w(0) = 1, w(d) = (e − 0.01) + (1 − e)·(0.01/(1 − e))^(d/R)**, R = (W − 1)/2,
  max deviation 6.9e-5 (= rounding). Five full KD profiles agree within 0.000496.
- Implemented additively: `HydrophobicityProfileOptions.WeightModel` (`HydropathyWeightModel.Linear` default |
  `Exponential`, Standard convention only), MCP `hydrophobicity_profile` `weightModel`. 14 NUnit cases + 1 MCP test, all
  locked to ProtScale output. Evidence: `docs/Evidence/SEQ-HYDRO-001-Evidence.md` "SRC RO2". Commit 4bb88d80.
  B03.md leftover struck + Post-completion entry; QUEUE R6 done.

## C. RO1 — Kazusa 316407 per-thousand column (SRC for B02) — DONE

FIN-B02 had landed (B02.md ends with `## Leftovers`), so this was done in place. `EColiK12.PerThousand` = Kazusa
species=316407 "/1000" column verbatim ("Escherichia coli W3110 [gbbct]: 4332 CDS's (1372057 codons)"; Σ number
1 372 057, Σ /1000 = 1000.0; the preset's 64 fractions equal Kazusa's fraction column). Per F39's sourced design
(%MinMax reference takes per-thousand usage) the default E. coli %MinMax now equals CHARMING `calculateMinMax`
(commit 2c75d1b): `CUGAGA` w2 36.470588… → 58.10113519091846 (CHARMING 58.1); `AGAAGAAGA` w3 −86.3636… →
−88.91170431211498; the `CodonFrequencies` basis keeps the old values. Sherlocc now accepts `EColiK12` by default.
Tests: Codon tier 854/854, Codon heavy 304/304, MCP MolTools 411/411. Evidence `CODON-RARE-001-Evidence.md`
"SRC RO1". Commits a9e8c94d (+ hash records 0606d2b5, c49b251c). B02.md leftover (b) struck; QUEUE R7 done.

## B. Citation audit (Wikipedia / Addgene)

### Census (script, not by reading files)

| | count |
|---|---|
| tracked `*.md` + `*.cs` files outside `wiki/` citing them | 421 |
| URL occurrences | 793 |
| distinct URLs | 205 (200 Wikipedia, 5 Addgene) |
| distinct URLs answering HTTP 200 on 2026-10-14 | 205 / 205 (3 via redirect) |

(The orchestrator's 428 / 790 / 206 differs by the counting regex: the first pass of mine truncated URLs at `,` / `'` /
`)` and swallowed one `</summary>`; corrected. 5 occurrences moved while B08 edited its files during the audit.)

Classification of the 788 occurrences of the first census pass (verdicts carried over to the corrected URLs, see
the audit file): **NUMBER 173, DEFINITION 79, BIBLIOGRAPHY 536.** Every NUMBER occurrence and every distinct
DEFINITION URL was checked against the fetched page by six verification passes (claim context vs page text, numbers
grepped, primaries looked up); full per-occurrence table with page quotes and primaries:
`docs/Evidence/SRC-Citation-Audit.md`.

### Verdicts (by occurrence)

| class | SUPPORTED | NO_NUMBER → definition OK | SUPPORTED_DEF | NOT_STATED / NOT_SUPPORTED_DEF | LINK-OK |
|---|---|---|---|---|---|
| NUMBER (173) | 97 | 55 | — | 21 | — |
| DEFINITION (79) | 12 | — | 61 | 6 | — |
| BIBLIOGRAPHY (536) | — | — | — | — | 536 |

So **27 occurrences (3.4 %) cited a page that does not state what was attributed to it**; plus 13 "verbatim"/detail
mismatches inside otherwise supported blocks (fixed or marked). No number that drives code or a test was wrong — in
every defect the value itself was right and rests on a primary already cited (or now cited) in the same document; the
defect was the attribution.

### Re-anchoring of numbers that drive code or tests

93 SUPPORTED occurrences back a number that drives code or a test expectation. For each the verification pass recorded
the primary source and where it states the number (paper / reference-tool source) — column "primary" in
`SRC-Citation-Audit.md`; Wikipedia stays as a secondary citation, now marked "re-fetched 2026-10-14" in 98 Evidence
sections. Explicitly re-anchored in this session:
- **PRIMER-DESIGN-001** (SRC for B07): Wikipedia "Primer (molecular biology)" states only 18–24 nt. Addgene states
  verbatim "Length of 18-24 bases; 40-60% G/C content; … Tm of 50-60°C; Primer pairs should have a Tm within 5°C".
  Primary defaults now in the TestSpec: Primer3 `libprimer3.c` `pr_set_default_global_args_1` (primer3-py 2.0.3 sdist):
  size 18/20/27, Tm 57/60/63, GC 20/80, poly-X 5, GC clamp 0, max end GC 5, pair ΔTm 100; EMBOSS 6.6.0 `eprimer3.acd`:
  minsize 18, optsize 20, maxsize 27, mintm 57, opttm 60, maxtm 63, mingc 20, ogcpercent 50, maxgc 80, maxdifftm 100,
  gcclamp 0, maxpolyx 5. The library's 40–60 % GC and ≤ 5 °C defaults rest on Addgene (a lab guideline, now the only
  source cited for them); MaxLength 25 is a library convention (between Primer3 27 and Addgene 24), stated as such.
- NMF (`OncologyAnalyzer.Signatures.cs`, SRC for B24): 10 000 iterations = SigProfilerExtractor `sigpro.py` l.341–343
  `min_nmf_iterations` default (max 1 000 000); ε floor practice = scikit-learn `_nmf.py` `EPSILON` (float32 eps);
  1e-12 stated as an implementation choice.
- Tajima's D ≈ 0.273 tests (POP-DIV): recomputed with Tajima 1989 eq. 38 (0.27345); Wikipedia stops at d = 0.08.

### Confirmed fixes (attribution defects; the numbers were right)

- F1 PRIMER-DESIGN-001: 40–60 % GC, Tm 50–60 °C, pair ΔTm ≤ 5 °C attributed to Wikipedia → Addgene + Primer3/EMBOSS anchors.
- F2 ANNOT-REPEAT-001: STR 1–6 bp and "two or more copies" (drive `maxUnitLength = 6` + a MUST test) → Tóth 2000, Benson 1999.
- F3 CHROM-CENT-001: 171 bp α-satellite (drives `AlphaSatelliteMonomerLength`) → PMC6121732, not the Centromere page.
- F4 SPLICE-DONOR-001: MAG|GURAGU −3..+6 (drives the donor window) → Shapiro & Senapathy 1987 / Yeo & Burge 2004.
- F5 SPLICE-ACCEPTOR-001: BPS 18–50 nt → Mercer 2015 / Gao 2008 (code uses 18–40); AT-AC → Burge 1999; U2AF65/35 marked.
- F6 ONCO-IMMUNE-001: LM22 547 genes × 22 cell types (drives ImmuneAnalyzer) → Newman 2015 / Chen 2018; ssGSEA → Yoshihara 2013.
- F7 SEQ-THERMO-001: self-complementary x = 1 → SantaLucia 1998 / Biopython `Tm_NN`; page formula quoted as on page.
- F8 POP-DIV tests (Genomics + MCP): D ≈ 0.273 → Tajima 1989 eq. 38; MCP comment's "a₁ ≈ 0.096" (that is θ_W) fixed.
- F9 OncologyAnalyzer.Signatures: NMF quote not on page; 10 000 / ε re-sourced (above).
- F10 CHROM-KARYO-001: "Nullisomy: 0 copies" (code maps 0 → "Nullisomy") not on the Aneuploidy page → re-sourced to Wikipedia "Nullisomic" ("lack of both the normal chromosomal pairs … (2n-2)").
- F11 PARSE-FASTQ-001: '@'/'I' heuristic not on page and stale → points to the revised FastQC `PhredEncoding.java` rule.
- F12 PARSE-GFF-001: "first 8 fields" → page says 7 (col 8 differs GFF2/GTF vs GFF3).
- F13 PARSE-GENBANK-001: "34+ trillion bp (Oct 2024)" → 53.9 trillion bases (release 271.0, Apr 2026), as the page now says.
- F14 META-PROF-001: sum-to-1 marked definitional (not on page). F15 META-ALPHA-001: Pielou J → Pielou 1966.
- F16 PROTMOTIF-DOMAIN-001: zinc finger "~50 aa" (wrong; PF00096 ~23–28) removed; h-region 7–15 → von Heijne 1985 (page: 5–16).
- F17 COMPGEN-COMPARE-001: 50.2 % (Acropora) marked ScienceDirect-only, not re-verifiable.
- F18 VARIANT-SNP-001: "four transversions per base" → two per base (four pair types); list marked derived.
- F19 SEQ-COMPLEX-COMPRESS-001: "verbatim" LZ76 pseudocode had two inverted conditions → corrected to the page.
- F20 ALIGN-GLOBAL-001, PROBE-VALID-001, SPLICE-PREDICT-001, MIRNA-PAIR-001, MOTIF-REGULATORY-001 (TATAWAW), META-BIN-001,
  TRANS-DIFF-001, PHYLO-DIST-001, META-TAXA test comment: quotes not on the page corrected or marked paraphrase/derived.

Commits: 0a56b38a, e19f148f, bd16257d, 61214e0d, 5a92af35, 6af1e9b9, 6aaf270e, 552e8b04 (fixer log in the audit
file's verdict notes). No number or test assertion changed in B.

### What stays "inherited, not re-verified"

Nothing in the audited set: all 205 URLs were re-fetched and every NUMBER/DEFINITION claim checked. Out-of-scope
occurrences (B02, B08 files) got verdicts but no edits — see Leftovers.

## Tests

- A: targeted `SequenceStatistics|Hydrophob` 795/795; MCP Analysis 433/433.
- C: Codon tier 854/854; Codon Fuzzing/Properties/Metamorphic 304/304; MCP MolTools 411/411.
- B: touched classes 110/110 + MCP Population 3/3.
- PrimerDesigner XML-doc edit first added a line and shifted the `docs/mcp` line anchors (5 `DocRef_LineAnchor` failures
  in MCP MolTools); fixed by keeping the edit on one line (MCP MolTools 413/413).
- **Final full fast tier** (after all commits): `dotnet build Seqeron.sln` 0 errors; `Seqeron.Genomics.Tests` fast filter
  **18415/18415**; MCP Analysis 433/433, MCP MolTools 413/413, MCP Population 56/56. Heavy tier: C ran the Codon
  Fuzzing/Properties/Metamorphic classes (304/304); A adds an option whose default path is unchanged; B is comment/doc only.

## Leftovers

DOABLE: none in SRC's scope.

Out of scope for this session (verified, listed — not edited):
- **B02** (finished; doc-only, registered as QUEUE R8): `docs/Evidence/CODON-OPT-001-Evidence.md:21` "0 < CAI ≤ 1"
  attributed to Wikipedia "Codon Adaptation Index", which states no range (the bound follows from w_i = f_i/max f ≤ 1;
  primary Sharp & Li 1987); `CODON-STATS-001-Evidence.md:25` log form exp[(1/L)Σ ln w_i] not on the page (identity);
  `CODON-USAGE-001-Evidence.md:17` RSCU / cosine-similarity bullets not on the cited page.
- **B08** (running; owns its citations): `docs/Evidence/RESTR-FILTER-001-Evidence.md:41` per-enzyme 4/6/8-bp site
  lengths cited to Wikipedia "Restriction site", which is an index page ("usually short (3 to 8 bp)", no enzyme table);
  `tests/.../RestrictionAnalyzer_FindSites_Tests.cs:~871` cases labelled "Wikipedia" include BglII and ApaI, which are not
  in the page's Examples table (the other 19 enzymes and cut positions match). Both values are correct per Biopython
  `Bio.Restriction` (B08's oracle); only the label is wrong.
- `wiki/**` not touched (orchestrator re-ingests).
