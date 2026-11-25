# B11 A9 hand-off: U12 log-odds typing (owner run required)

Status: BLOCKED on obtaining the published numeric U12 matrices; nothing was reclassified (README Definition of Done).
Date of attempt: 2026-11-25. Audit item: B11.md "Audit worklist (round 1)" A9.

## Current behaviour (to be replaced)
`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs`
- `U12DonorConsensus` (:256) and `ScoreU12DonorSite` (:1056): fraction of consensus positions matched (/RTATCCTY).
- `U12BranchSiteConsensus = "UUCCUUAAC"` (:342), `ScoreU12AcceptorSite` (:1078): best `ConsensusMatchFraction` of the branch site 8-21 nt upstream of the 3' ss.
- `DefaultU12MinSiteScore = 0.75` (:1105): unpublished cut-off on a match fraction, used by `PredictIntrons`, `PredictGeneStructure` (:1354 path) and `DetectAlternativeSplicing` (:1499).

## Why blocked (every source tried)
Published numeric tables needed: U12 5' ss and branch-site frequency matrices (Sharp & Burge 1997 Cell 91:875; Burge, Padgett & Sharp 1998 Mol Cell 2:773; Burge, Tuschl & Sharp 1999 RNA World 2nd ed.; Levine & Durbin 2001 NAR 29:4006 - supplement frequency tables; Alioto 2007 NAR 35:D110; Sheth 2006 NAR 34:3955; Moyer 2020 NAR 48:7066).
- Levine & Durbin 2001: `pmc.ncbi.nlm.nih.gov/articles/PMC60238/` opened (brief's PMC55920 is a different, unrelated paper). Text states the matrices are in Supplementary Material (`nar_29_19_4006__gke532table1.xls`, `...table2.xls`; authors' site www.sanger.ac.uk/Users/rd/U12/); the xls links return 404 and Europe PMC reports the article "not open access" for supplementaryFiles. Published decision rule IS available in the text: log-odds in bits, U12 donor > 9 bits and branch > 6 bits (AT-AC and GT-AG scans), branch-to-3'ss distance 8-21 nt, matrices trained on 48 non-redundant U12 introns with genomic-background pseudocounts (reused by Levine's thesis, `ftp.sanger.ac.uk/pub/resources/theses/levine/chapter4.pdf`, opened with pdftotext: no matrix tables).
- Sharp & Burge 1997: `cell.com/cell/fulltext/S0092-8674(00)80479-1` HTTP 403; not in PMC / Europe PMC.
- Burge 1998 / Burge 1999: not open access (not reachable).
- Alioto 2007 U12DB `PMC1635337` full text via Europe PMC REST: no tables; `genome.crg.cat/datasets/u12/` gives only the consensus RTATCCTTT / TTCCTTRAY and an FTP SQL dump (`ftp://genome.imim.es/pub/software/u12/u12db_v1_0.sql.gz`, GPL, not reachable here and not a printed paper table).
- Sheth 2006 `PMC1557818` (full text + supplement zip, figures only), Moyer 2020 `PMC7367187` (+ supplement PDF: 107 lines, no matrices), Olthof 2019 `PMC6717393` (supplement xlsx = intron lists), Taxonomy of introns 2024 `PMC11347168` (supplement: PWMs shown as logos only), Turunen 2013 review `PMC3584512` (no matrix table), Madan 2015 `PMC4349895`: no printed numeric U12 matrices. The SpliceRack / intronIC matrices are GPL data files and must not be bundled.
- `rest.ensembl.org`, Europe PMC REST full text (curl) work; ncbi/pmc supplementary xls 404.

## Planned implementation (when matrices are in hand)
Owner supplies the numbers printed in the papers (Sharp & Burge 1997 Tables / Burge 1999 tables, or Levine & Durbin 2001 supplementary Table 1/2 frequency tables) as a plain CSV (position, A, C, G, U counts or frequencies) in `docs/Evidence/SPLICE-PREDICT-001/u12_matrices.csv` with citation.
1. `SpliceSitePredictor.U12LogOddsDonor(string window)` / `U12LogOddsBranch(string window)`: log2(f_U12(b,pos) / f_bg(b)) summed over positions, pseudocounts from background as in Levine & Durbin (genomic base frequencies); background = uniform 0.25 unless the paper states the genome frequencies. Matrices as `private static readonly double[,]` with a citation comment.
2. Typing rule: `PredictIntrons` types an intron U12 iff donor log-odds > 9 bits and branch log-odds > 6 bits (Levine & Durbin 2001), or, if the transcribed matrix differs, the threshold published with that matrix; log-odds against the existing U2 model (S&S / MaxEnt) as the optional alternative, `U12 vs U2 log-odds > 0` (Sharp & Burge 1997 classification). New overload parameters `u12MinDonorBits`, `u12MinBranchBits` with these defaults; the match-fraction overload `PredictIntrons(..., u12MinSiteScore)` kept as `[Obsolete]`-free legacy only if the existing tests still need it (otherwise removed, update callers at :1122, :1354, :1499).
3. `FindDonorSites` / `FindAcceptorSites` U12 branches (:503, :638) report log-odds converted to the existing 0-1 `Score` through a documented monotone map, or expose bits in a new `Bits` property (additive record field).
4. Update TestSpec SPLICE-PREDICT-001 / SPLICE-DONOR-001 / SPLICE-ACCEPTOR-001, Evidence, `docs/algorithms/Splicing/*`, MCP `predict_introns`/`find_donor_sites` docs.

## Oracle
Independent Python (`scratchpad/oracle_u12.py`, not committed): load the same CSV, compute log-odds with numpy, score random sequences with planted /RTATCCTY donors and UUCCUUAAC branch sites, compare to the C# values through a small console driver (pattern: `drv` project referencing `Seqeron.Genomics.Annotation.csproj`, as used for A10). Pass criterion: |bits difference| < 1e-9 on >= 1,000 random windows, identical U12 calls on >= 300 random sequences.
Published examples: fetch U12 introns from Ensembl REST (`rest.ensembl.org/lookup/symbol/homo_sapiens/<SYMBOL>?expand=1`, `/sequence/id/<ENSG>?type=genomic`) for genes named in Levine & Durbin / U12DB (e.g. SCN4A, SCN5A, CACNA1, MAPK1/ERK2 intron U12 AT-AC, TMEM...); report the type call per known U12 intron and the fraction of U2 introns (HBB, TP53, CFTR canonical introns) wrongly typed U12 (expected 0).

## Tests to add
Exact bits for a consensus window and a one-mismatch window (locked to the CSV), U12 call on a planted AT-AC intron, no U12 call for U2 introns of HBB/TP53, threshold parameter boundary (= threshold is not U12), legacy overload behaviour.

## Verification
```
export DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
dotnet build Seqeron.sln -c Debug 2>&1 | grep -E " error |Error\(s\)" | sort -u | tail    # only CS8978
dotnet test tests/Seqeron/Seqeron.Genomics.Tests -c Debug --no-build --filter "FullyQualifiedName~Splic|FullyQualifiedName~SpliceSite|FullyQualifiedName~Transcriptome"
dotnet test tests/Seqeron/Seqeron.Mcp.Annotation.Tests -c Debug --no-build
```
