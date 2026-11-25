# Licences and assemblies (B09 WP79)

Status: permanent record. Every third-party program or library whose source code was read for a port anywhere in the repository (all review batches, not only B09), its licence, what the licence demands, and where the derived code lives.
Maintained together with `THIRD-PARTY-NOTICES.md` (repository root), the `LICENSE` file of every copyleft assembly and the architecture guard `tests/Seqeron/Seqeron.Genomics.Tests/Architecture/LicenceSeparationTests.cs`.
I am not a lawyer: this is the technical separation the owner asked for ("variant 1 - separate GPL assemblies"), not legal advice.

## 1. The decision and the layout

* The ports made by reading the source code of a program are derivative works of that program. The repository's main assemblies are MIT.
* Ports of GPL / LGPL programs live in **separate assemblies named by their subject** (never "Gpl" / "Copyleft"), one assembly per subject area, each with its own `LICENSE` file, an SPDX header (`// SPDX-License-Identifier: ...`, origin note) in every `.cs` file, `<PackageLicenseExpression>` in the csproj, its own test project(s) (moved tests and oracle data; the root namespace of the test projects stays `Seqeron.Genomics.Tests` so embedded resource names are unchanged) and, for the MCP tools, an MCP server project of its own with the same tool names as before.
* **The MIT assemblies do not depend on the copyleft assemblies** (no project reference, no `InternalsVisibleTo`, no assembly reference, no reflection string, no MCP tool reference). The copyleft assemblies depend on the MIT ones; the shared building blocks they need (cancellation scope, input-contract helpers, `PairScoring`, `BuildConsensus`, ...) were made `public` in the MIT assemblies.
* Permissive programs (MAFFT, MUSCLE, SPOA, minimap2 / ksw2, Mash, Biopython, parasail, NCBI toolkit, maxentpy, CNVkit, MHCflurry, MobiDB-lite, LICHeE, ...) stay in the MIT assemblies; their notices and licence texts are in `THIRD-PARTY-NOTICES.md`.
* Programs with a custom / non-commercial / unknown licence are NOT decided here: section 5 "NEEDS OWNER DECISION".
* Guard: `LicenceSeparationTests` (no MIT project / assembly / source file references a copyleft assembly; every `.cs` file of a copyleft project and of its test project carries the SPDX header of its licence; every copyleft project has a LICENSE file).

## 2. Copyleft assemblies

| Assembly | SPDX | Source files | Lines | Test project(s) / MCP server | Holds |
|---|---|---|---|---|---|
| `Seqeron.Genomics.ClustalW` | LGPL-3.0-or-later | 23 | 11,685 | `tests/Seqeron/Seqeron.Genomics.ClustalW.Tests` | ClustalW 2.1 (`ClustalWAligner`: pairwise / quick / profile / multiple stages, trees, command line, output formats, matrices) and the guide-tree progressive / iterative MSA family (`MultipleAlignProgressive`, `MultipleAlignIterative`, `ProgressiveMsaOptions`, guide trees, distances) |
| `Seqeron.Genomics.Kalign` | GPL-3.0-or-later | 4 | 3,115 | `tests/Seqeron/Seqeron.Genomics.Kalign.Tests` | Kalign 3.3.5 (`KalignAligner`, `KalignProgram`, MSF / Clustal / FASTA I/O) |
| `Seqeron.Genomics.PartialOrderAlignment` | GPL-2.0-or-later | 4 | 2,424 | `tests/Seqeron/Seqeron.Genomics.PartialOrderAlignment.Tests` | POA 2.0 (`PoaAligner`, `PoaProgram`, score matrices, PO / PIR / CLUSTAL output) |
| `Seqeron.Genomics.TCoffee` | GPL-2.0-or-later | 22 | 11,530 | `tests/Seqeron/Seqeron.Genomics.TCoffee.Tests` | T-Coffee 13.41 (`TCoffeeAligner`: libraries, extension, proba / pair methods, structure glue, tree, evaluate, program options) incl. the consistency aligner `MultipleAlignConsistency` (`ConsistencyMsaOptions`) and ProbCons 1.12 (`ProbConsAligner`, `ProbConsProgram`) |

| MCP server | SPDX | Source files | Tool docs | Tests |
|---|---|---|---|---|
| `Seqeron.Mcp.MultipleAlignment` | GPL-3.0-or-later | 12 | 25 (docs/mcp/tools/multiplealignment) | `tests/Seqeron/Seqeron.Mcp.MultipleAlignment.Tests` |

The MIT assemblies are listed in the root `LICENSE` file context: every other project under `src/` and `tests/` is MIT.

## 3. Inventory of the ported programs

### 3.1 Copyleft (separate assemblies)

| Program | Version | Licence (as the upstream states it) | Class | Where | Licence text taken from |
|---|---|---|---|---|---|
| EMBOSS (applications: needle, water, stretcher, matcher, supermatcher, needleall, merger, megamerger, seqmatchall, wordmatch, diffseq, est2genome, dotmatcher, showalign, infoalign, distmat, plotcon, consambig, fuzznuc, fuzzpro, fuzztran, patmatmotifs, patmatdb, prosextract, etandem, equicktandem, getorf, tcode, checktrans, plotorf, transeq, seqret and the other sequence / feature programs) | 6.6.0 | GNU General Public License, version 2 or (at your option) any later version (file headers of the application sources, e.g. emboss/needle.c, dotmatcher.c) | strong copyleft (GPL-2.0-or-later) | planned: Seqeron.Genomics.Emboss.* (programs, by area); NOT YET SEPARATED - currently Alignment / Analysis assemblies (section 6) | EMBOSS 6.6.0 source (header of each file) and /usr/share/doc/emboss/copyright (Debian emboss 6.6.0+dfsg-12ubuntu2) |
| EMBOSS libraries AJAX and NUCLEUS (sequence / feature / report readers and writers, ACD, USA, databases, alignment output formats, PCRE, zlib deflate) | 6.6.0 | GNU Lesser General Public License, version 2.1 or (at your option) any later version (headers of ajax/*.c, nucleus/*.c); the embedded PCRE (University of Cambridge licence) and the zlib deflate (zlib licence) keep their own terms | weak copyleft (LGPL-2.1-or-later) | planned: Seqeron.Genomics.Emboss.SequenceInput (libraries); NOT YET SEPARATED - currently Infrastructure / Alignment / Analysis (section 6) | EMBOSS 6.6.0 source headers; Debian copyright file (License: emboss = AJAX / NUCLEUS LGPL, applications GPL; PCRE_license; zlib) |
| T-Coffee | 13.41.0 (Debian 13.45.0.4846264+really13.41.0.28bdc39+dfsg-1 sources) | GNU General Public License, version 2 or (at your option) any later version (file headers, Debian copyright: GPL-2+) | strong copyleft (GPL-2.0-or-later) | Seqeron.Genomics.TCoffee | T-Coffee source file headers; /usr/share/doc/t-coffee/copyright |
| ProbCons (hosted in the T-Coffee assembly because the port runs on the pair-HMM tables embedded in T-Coffee) | 1.12 | Public domain (ProbCons README: 'PROBCONS has been made freely available as PUBLIC DOMAIN software') | public domain | Seqeron.Genomics.TCoffee (ProbConsAligner*, ProbConsProgram) | /usr/share/doc/probcons/copyright |
| Kalign | 3.3.5 | GNU General Public License, version 3 or (at your option) any later version | strong copyleft (GPL-3.0-or-later) | Seqeron.Genomics.Kalign | Kalign 3.3.5 source (COPYING) and /usr/share/doc/kalign/copyright |
| ClustalW | 2.1 (Debian clustalw 2.1+lgpl-7) | GNU Lesser General Public License, version 3 or (at your option) any later version | weak copyleft (LGPL-3.0-or-later) | Seqeron.Genomics.ClustalW (also hosts the guide-tree progressive / iterative MSA family built on ClustalW's distances and trees) | /usr/share/doc/clustalw/copyright (License: LGPL-3+) |
| POA (Partial Order Alignment) | 2.0 | GNU General Public License, version 2 or (at your option) any later version | strong copyleft (GPL-2.0-or-later) | Seqeron.Genomics.PartialOrderAlignment | POA 2.0 source (COPYING) and /usr/share/doc/poa/copyright |


### 3.2 Permissive and public domain (kept in the MIT assemblies, NOTICE)

| Program | Version | Licence | Class | Where | Licence text taken from |
|---|---|---|---|---|---|
| MAFFT | 7.505 | MAFFT licence (BSD-style: retain the notice, reproduce it in binary distributions, no endorsement) | permissive | NOTICE (THIRD-PARTY-NOTICES.md) | /usr/share/doc/mafft/copyright |
| MUSCLE | 3.8.1551 | Public domain dedication ('MUSCLE is public domain software') | public domain | NOTICE | /usr/share/doc/muscle3/copyright |
| SPOA | 4.0.8 | MIT | permissive | NOTICE | upstream LICENSE (Robert Vaser) |
| minimap2 / ksw2 | 2.26 | MIT | permissive | NOTICE | upstream LICENSE.txt (Heng Li, Dana-Farber, Broad Institute) |
| Mash (core) | 2.3 | BSD 3-Clause (Battelle National Biodefense Institute); kseq.h MIT; MurmurHash3 public domain; robin_hood and Cap'n Proto MIT | permissive | NOTICE | /usr/share/doc/mash/copyright |
| Biopython | 1.88 | Biopython License Agreement (dual with BSD 3-Clause for some files) | permissive | NOTICE | Biopython LICENSE.rst (installed package) |
| parasail | 1.3.4 | Battelle BSD-style licence (PyPI: BSD) | permissive (custom BSD-style wording) | NOTICE | upstream LICENSE.md |
| NCBI C++ Toolkit / BLAST+ (dustmasker, ORFfinder, BLAST statistics) | 2.12.0 | Public domain (United States Government work) | public domain | NOTICE | /usr/share/doc/ncbi-blast+/copyright |
| HMMER 3 / Easel (profile HMM text format facts, E-value definitions) | 3.4 | BSD 3-Clause (HMMER), BSD 2-Clause (Easel) | permissive | none required (file-format facts and published formulas only) | upstream LICENSE |
| maxentpy (MaxEntScan tables) | master | MIT | permissive | NOTICE | upstream LICENSE |
| CNVkit (cnvlib/call.py arithmetic) | master | Apache License 2.0 | permissive | NOTICE (+ Apache-2.0 text in LICENSES.md appendix) | upstream LICENSE |
| MHCflurry 2.0 (network architecture, pseudosequences, weights) | 2.1.5 | Apache License 2.0 | permissive | NOTICE (+ Apache-2.0 text) | upstream LICENSE |
| MobiDB-lite (thresholds) | v3 | Apache License 2.0 | permissive | NOTICE (+ Apache-2.0 text) | upstream LICENSE |
| LICHeE (PHYNetwork constraint network) | master | MIT (Copyright (c) 2014 Victoria Popic) | permissive | NOTICE | upstream LICENSE |
| Pfam profile HMMs (data) | - | CC0 1.0 | public domain | NOTICE | Pfam licence statement |
| ABIS-Seq signature matrix (data) | - | CC BY 4.0 | permissive (attribution) | NOTICE | Monaco et al. 2019 |


### 3.3 Programs of the other batches (read for a port, a verbatim table or a behaviour check)

Evidence = the batch report / Evidence document of the unit (`docs/Evidence/<UNIT>-Evidence.md`, `docs/Validation/review-2026-09/B*.md`) and the header of the C# file. "port" = the C# code follows the source of the program (identifiers, structure, constants). "consulted" = the source was read to confirm behaviour or a formula; the batch states that no code was copied (or the file is an independent implementation of a published algorithm). Licence texts were fetched from the upstream repositories (`raw.githubusercontent.com`, neutral User-Agent `Seqeron-review/1.0`) or taken from the package metadata.

| Program | Upstream licence (as stated) | What the repository did (evidence) | Class | Action |
|---|---|---|---|---|
| primer3 / primer3-py `thal.c` (ntthal) | GNU GPL version 2 (LICENSE of primer3-org/primer3 and libnano/primer3-py) | `MolTools/NtthalDimer.cs`, `NtthalHairpin.cs`: "faithful port" / "bit-faithful port of primer3-py 2.3.1 thal.c" (PRIMER-TM-001-SPECIAL-LOOP-Evidence) | strong copyleft | NOT YET SEPARATED (assembly planned, see section 6) |
| DECIPHER 3.9.4 `ConsensusSequence` | GPL-3 (DESCRIPTION of the Bioconductor package) | `Analysis/MotifFinder.DecipherConsensus.cs`: "a line-by-line port of DECIPHER 3.9.4" | strong copyleft | NOT YET SEPARATED (assembly planned, see section 6) |
| RSAT `oligo-analysis` v1.169, `dyad-analysis` v1.78, `RSA.disco.lib` (rsa-tools/rsat-code 10043f2) | GNU Affero General Public License v3 (LICENSE.txt of rsa-tools/rsat-code) | `Analysis/MotifFinder.OligoAnalysis*.cs`, `MotifFinder.DyadAnalysis.cs`, `MotifFinder.SharedMotifsDegenerate.cs`, `MotifFinder.OligoSequenceType.cs`: "Ported from rsa-tools/rsat-code" | strong copyleft (network clause) | NOT YET SEPARATED (assembly planned, see section 6) |
| Ckmeans.1d.dp (R package) | LGPL (>= 3) (DESCRIPTION) | `Oncology/OncologyAnalyzer.Clonality.cs`: "a line-by-line port" of the optimal 1-D k-means | weak copyleft | NOT YET SEPARATED (assembly planned, see section 6) |
| ASCAT / Battenberg (VanLoo-lab) | GNU GPL v3 (LICENSE of VanLoo-lab/ascat; Battenberg is GPL-3 as well) | `Oncology/OncologyAnalyzer.CopyNumberPloidy.cs`: ASCAT `ascat.aspcf`, `ascat.runAscat`, `determine_copynumber`, Battenberg clonal-state arithmetic ported function by function (ONCO-ASCAT-001-Evidence, F58 ...) | strong copyleft | NOT YET SEPARATED (assembly planned, see section 6) |
| micropan (R package) | GPL-2 (DESCRIPTION) | `Metagenomics/PanGenomeAnalyzer.cs`: lower / upper bounds of the fitted Heaps' law parameters "copied verbatim" (numeric constants) | strong copyleft (two constants) | two numeric constants: listed in section 5 |
| SkewIT | GNU GPL v3 (LICENSE) | `Analysis/GcSkewCalculator.cs`: the decision steps of `skewi.py` are described "verbatim"; the genus threshold table is NOT embedded (caller-supplied, B01 F31) | strong copyleft (algorithm steps only) | algorithm steps only (no table): listed in section 5 |
| Tandem Repeats Finder 4.10 | GNU Affero GPL v3 | `Analysis/RepeatFinder.cs`: source read "for behaviour only; nothing transliterated"; the simulated tables are NOT shipped (REP-APPROX-001-Evidence, B04 F43 / F56) | consulted | none (no code copied) |
| MISA | GNU GPL (misa.pl) | `RepeatFinder.cs` doc: the "reject false type motifs" rule is the behaviour of MISA; independent implementation (REP-STR-001-Evidence) | consulted | none |
| miRDeep2 | GNU GPL | MIRNA-PRECURSOR-001-Evidence: "no GPL miRDeep2 code or weights were ported, copied, or linked" | consulted | none |
| CRISPOR / CHOPCHOP / cas-offinder | various (CRISPOR: academic licence on the web server code; no LICENSE file at crisporWebsite) | `MolTools/CrisprDesigner.cs`: published scores re-implemented; constants cited, no code ported (REP / CRISPR evidence) | consulted | none; NEEDS OWNER DECISION only if a later audit shows copied code (listed in section 5) |
| Azimuth Rule Set 2 (Microsoft Research) | BSD 3-Clause ("Copyright (c) 2015, Microsoft Research") | `MolTools`: Rule Set 2 model reproduced; oracle vectors from `scripts/azimuth/extract_azimuth_model.py` | permissive | NOTICE entry (Azimuth) |
| TargetScan 7 scripts (`targetscan_70_BL_PCT.pl`) | no LICENSE file in the nsoranzo/targetscan mirror | `Annotation/MiRnaAnalyzer.cs`: PCT arithmetic "retrieved verbatim this session" (comment) | unknown | NEEDS OWNER DECISION (section 5) |
| HMMER 3 / Easel | BSD 3-Clause / BSD 2-Clause | `Analysis/Plan7ProfileHmm.cs`: file-format facts "verbatim from the HMMER User's Guide" | permissive | none beyond NOTICE (facts) |
| NCBI C++ Toolkit (dustmasker `symdust.cpp`, ORFfinder `COrf`) | public domain (US Government work) | `Analysis/SequenceComplexity.Dustmasker.cs`, `Core/OrfScanner.cs` rules | public domain | NOTICE |
| MobiDB-lite v3 | Apache-2.0 | `Analysis/DisorderPredictor.cs`: thresholds verbatim | permissive | NOTICE (+ Apache-2.0 text) |
| CNVkit `cnvlib/call.py` | Apache-2.0 | `Chromosome|Oncology/CopyNumberMath.cs`: "Ported from CNVkit" | permissive | NOTICE (+ Apache-2.0 text) |
| MHCflurry 2.0 | Apache-2.0 | `Oncology/MhcflurryAffinityPredictor.cs`: "faithful C# port", weights + pseudosequences embedded | permissive | NOTICE (+ Apache-2.0 text) |
| LICHeE | MIT | `Oncology/OncologyAnalyzer.PhylogenyHeterogeneity.cs`: "Port of LICHeE's PHYNetwork" | permissive | NOTICE |
| maxentpy (MaxEntScan tables) | MIT | embedded `Data/maxent_score3.txt`, score5 factorisation (SPLICE-ACCEPTOR-001-Evidence) | permissive | NOTICE |
| AGFusion | MIT (per the project) | `Oncology/OncologyAnalyzer.Fusions.cs`: model table "verbatim from AGFusion model.py" | permissive | NOTICE (see section 5, LICENSE file not fetched: HTTP 404) |
| BBMap / BBDuk | BBTools licence (custom) | `Analysis/SequenceComplexity.cs`: entropy masking port | custom | NEEDS OWNER DECISION (section 5) |
| Pfam, ABIS-Seq, REBASE (via Biopython), Cosmic / ClinVar-style lookup tables | CC0, CC BY 4.0, Biopython licence (REBASE allenz columns: Copyright Richard J. Roberts - redistribution question, see B08 report) | embedded data resources | data | NOTICE; REBASE allenz columns: NEEDS OWNER DECISION |
| scipy, numpy, R (stats), Jellyfish, CD-HIT, bwa, seqtk, sdust, GATK, Ensembl VEP, maftools, sigminer, deconstructSigs, CheckM, QUAST, InterProScan, ViennaRNA, pknotsRG, ChouFasman, seqinr, ape, micropan (rest), epigeneticclock, methylKit, INVAR, facets, msisensor, LOHHLA, oncokb, STAR-Fusion, arriba, breakdancer, cBioPortal | various | cited as oracles or as the published source of a formula; the batch reports record no ported code | consulted / oracle | none |



## 4. Black-box oracles (run, output recorded; NOT derivative)

These programs were executed (installed Debian / Ubuntu packages, or user-run Docker packages) and their output compared with the ports; no source code of them was ported by the entries below. The oracle generators are the scripts under `scripts/alignment/generate_*_oracle.py` and `docs/Validation/review-2026-09/handoff/`; the recorded outputs are the TSV files under the `TestData` directories of the test projects.

| Program | Version / package |
|---|---|
| EMBOSS 6.6.0 binaries (all programs of the table above, run and their output recorded) | emboss 6.6.0+dfsg-12ubuntu2 |
| t-coffee 13.41 | t-coffee 13.45.0.4846264+really13.41.0.28bdc39+dfsg-1 |
| kalign 3.3.5 | kalign 1:3.3.5-1 |
| clustalw 2.1 (also a patched build in the user-run Docker package handoff/A1-10d) | clustalw 2.1+lgpl-7 |
| muscle 3.8.1551 | muscle3 3.8.1551-3 |
| mafft 7.505 | mafft 7.505-1 |
| probcons 1.12 | probcons 1.12-14 |
| poa 2.0, spoa 4.0.8, abpoa 1.5.1 | poa 2.0+20060928-9; spoa 4.0.8-1; abpoa 1.5.1-1build1 |
| minimap2 2.26 | minimap2 2.26+dfsg-1build1 |
| MUMmer 3.23 | mummer 3.23+dfsg-8 |
| Mash 2.3 | mash 2.3+dfsg-6 |
| fasta3 / fasta36 36.3.8i | fasta3 36.3.8i.14-Nov-2020-1 |
| NCBI BLAST+ 2.12.0, ORFfinder | ncbi-blast+ 2.12.0+ds-4build2 |
| LAST 1542 (Docker package handoff/A2-19c-last; MAFFT --lastpair) | last-align 1542-1 |
| TM-align 20190822, MUSTANG 3.2.4, SAP 1.1.3, DaliLite v3 / v5 (T-Coffee structure methods; run in user Docker packages, output parsed by the T-Coffee glue port) | tm-align 20190822+dfsg-2build1; mustang 3.2.4-1; SAP 1.1.3 (github.com/mathbio-nimr-mrc-ac-uk/SAP); DaliLite (ekhidna2.biocenter.helsinki.fi) |
| exonerate 2.4.0, hhsuite 3.3.0, DIALIGN / DIALIGN-TX, PRANK, Clustal Omega, FSA, AMAP, ProDA (installed; formats / comparison only, nothing ported) | see /usr/share/doc |
| Biopython 1.88, parasail 1.3.4, scipy, numpy (Python reference implementations) | pip |


## 5. NEEDS OWNER DECISION (custom / non-commercial / unknown licence; nothing was moved or removed for these)

The code stays where it is until the owner decides (relicense request, replace by an independent implementation, move to a separate assembly, or remove). Evidence for each item is the licence text quoted.

| Program | Version | Licence text / evidence | Class | Files |
|---|---|---|---|---|
| FASTA 3.4 (University of Virginia; the alignment half and scoring ported for the MAFFT `--fastapair` / `--fastswpair` / `--fastaparttree` options and the `fasta34` oracle) | 3.4 (fasta34t*; ships with MAFFT's `fasta34` dependency) | Custom: 'Copyright 1988 ... 1995 William R. Pearson and the University of Virginia. All rights reserved. The FASTA program and documentation may not be sold or incorporated into a commercial product, in whole or in part, without written consent of William R. Pearson and the University of Virginia.' (COPYRIGHT file of the FASTA 3.4 source); smith_waterman_sse2.* additionally 'copyright (c) 2006 by Michael Farrar ... may not be sold or incorporated into a commercial product ... without written consent'. | custom / commercial-use restriction (NOT compatible with MIT or GPL as written) | src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.MafftFasta34*.cs, SequenceAligner.MafftDndfast7.cs (calls), MafftProgram*.cs paths that use them; tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAligner_MafftFasta34*_Tests.cs; docs F314 / F323 / F325 / F346-F349. The later FASTA 36 releases (fasta3 36.3.8i in Debian) are Apache-2.0, but the code ported here is the FASTA 3.4 source. |
| MUMmer (nucmer, promer, mummer, show-coords, show-snps, show-aligns, show-diff, show-tiling, delta-filter, dnadiff, mummerplot, mapview, nucmer2xfig, repeat-match, exact-tandems, combineMUMs, gaps, delta2blocks, mummer-annotate, delta2maf, run-mummer1/3) | 3.23 | Artistic License (the Artistic License 1.0 as shipped with Debian: 'Copyright (c) 2003 by Stefan Kurtz and The Institute for Genomic Research' - Debian copyright file License: Artistic) | custom (Artistic 1.0: modified versions must be marked and the standard version made available; neither permissive nor GPL) | src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/Mummer*.cs, Nucmer*.cs, Promer*.cs, LibstdcxxSort.cs, AnchorBasedAligner.cs (chaining reference only); MCP tools AnalysisTools.Mummer*.cs; docs F150, F199-F213 |
| Mash Open Bloom Filter (bloom_filter.hpp, Arash Partow) used by `mash dist -b` / screen port | Mash 2.3 | Common Public License Version 1.0 (Debian copyright file of mash, Files: src/mash/bloom_filter.hpp) | weak copyleft, not GPL / LGPL (CPL-1.0: modifications of the program must be released under the CPL) | src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MashBloomFilter.cs (F197) |
| BBTools / BBMap (BBDuk.maskLowEntropy, EntropyTracker; entropy masking port) | 40.02 | BBTools licence: BSD 3-Clause text of the University of California / Lawrence Berkeley National Laboratory plus the 'Enhancements' paragraph (a perpetual licence for improvements that are made available) | custom (BSD-style with an additional grant clause) | src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs (entropy masking) |
| parasail (score API reproduced) | 1.3.4 | Battelle's BSD-style licence text (custom wording, listed by PyPI as BSD) | permissive-style custom wording (kept in the MIT assemblies and in NOTICE; listed for completeness) | src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ParasailAligner.cs |
| PCRE 3.x embedded in EMBOSS (ajax/pcre) | EMBOSS 6.6.0 | PCRE licence (University of Cambridge): free use with an acknowledgement sentence in the documentation, altered versions marked | custom (BSD-style with acknowledgement); lives in the EMBOSS motif assembly under that assembly's copyleft licence | src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/EmbossPcre.cs |
| libstdc++ std::sort order (LibstdcxxSort.cs: reproduction of the introsort of GCC's libstdc++ to match the order of MUMmer / minimap2 outputs) | GCC 13 | GNU General Public License v3 with the GCC Runtime Library Exception (libstdc++) | strong copyleft with runtime exception; the algorithm is the published introsort - the file reproduces the order, not the text | src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/LibstdcxxSort.cs |
| glibc random() / rand() (GlibcRandom.cs: TYPE_3 additive feedback generator, to reproduce `rand()` driven outputs) | glibc 2.39 | GNU Lesser General Public License v2.1 or later | weak copyleft; the generator is the published algorithm (r[i] = r[i-3] + r[i-31]) - the file reproduces the algorithm, not the text | src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/GlibcRandom.cs |


## 6. Status of the separation (what is done, what is registered)

Done in this work package (commit "B09 WP79 ..."): Kalign, POA, ClustalW (+ the guide-tree progressive / iterative MSA family), T-Coffee (+ the consistency aligner and ProbCons) are separate assemblies with LICENSE files, SPDX headers, tests and an MCP server (`Seqeron.Mcp.MultipleAlignment`, 25 tools; same tool names); the MIT assemblies were cleaned of every reference to them and `LicenceSeparationTests` enforces it.
Shared MIT building blocks made `public` for the copyleft assemblies: `AlignmentCancellation`, `SequenceAligner.RequireProgramInput` / `RequireProgramCall` / `RequireNoNullElements` / `RequireCells(..., owner)`, `BuildConsensus`, `PairScoring` and its helper types, `ParallelFor`, `LinearAlignCore` and the other helpers the compiler reported (they carry XML documentation).
Behaviour changes to MIT APIs forced by the separation (documented in B09.md F363): `MsaOptions` keeps only the generic members (`GapModel`, `ExactMaxStateCells`, `ExactCarrilloLipman`, `ExactLowerBoundHint`, `MaxDegreeOfParallelism`, `CancellationToken`); the guide-tree members moved to `ProgressiveMsaOptions` (ClustalW assembly) and `ConsistencyMsaOptions` (T-Coffee assembly); `MultipleAlignExact` uses the star alignment as its Carrillo-Lipman lower bound unless `ExactLowerBoundHint` is given (the optimum does not depend on it).

Registered, not yet separated (each is a straight application of the same procedure; the lead schedules them): see B09.md "Leftovers" of WP79 - EMBOSS (6.6.0, the largest port: assemblies by area SequenceInput / PairwiseAlignment / Motifs / Orf / Repeats / ReportFormats), and the ports of other batches that the inventory above classifies as copyleft (primer3 ntthal, DECIPHER ConsensusSequence, RSAT oligo / dyad analysis, Ckmeans.1d.dp, ASCAT / Battenberg). Until they are moved they are the only copyleft-derived code that sits in an MIT assembly; `LICENSES.md` is the list the guard test will be extended with.

## Appendix A. Apache License 2.0 (text)

Applies to CNVkit, MHCflurry and MobiDB-lite (see `THIRD-PARTY-NOTICES.md`); the copyright lines of the individual projects are in the notices.

```

                                 Apache License
                           Version 2.0, January 2004
                        http://www.apache.org/licenses/

   TERMS AND CONDITIONS FOR USE, REPRODUCTION, AND DISTRIBUTION

   1. Definitions.

      "License" shall mean the terms and conditions for use, reproduction,
      and distribution as defined by Sections 1 through 9 of this document.

      "Licensor" shall mean the copyright owner or entity authorized by
      the copyright owner that is granting the License.

      "Legal Entity" shall mean the union of the acting entity and all
      other entities that control, are controlled by, or are under common
      control with that entity. For the purposes of this definition,
      "control" means (i) the power, direct or indirect, to cause the
      direction or management of such entity, whether by contract or
      otherwise, or (ii) ownership of fifty percent (50%) or more of the
      outstanding shares, or (iii) beneficial ownership of such entity.

      "You" (or "Your") shall mean an individual or Legal Entity
      exercising permissions granted by this License.

      "Source" form shall mean the preferred form for making modifications,
      including but not limited to software source code, documentation
      source, and configuration files.

      "Object" form shall mean any form resulting from mechanical
      transformation or translation of a Source form, including but
      not limited to compiled object code, generated documentation,
      and conversions to other media types.

      "Work" shall mean the work of authorship, whether in Source or
      Object form, made available under the License, as indicated by a
      copyright notice that is included in or attached to the work
      (an example is provided in the Appendix below).

      "Derivative Works" shall mean any work, whether in Source or Object
      form, that is based on (or derived from) the Work and for which the
      editorial revisions, annotations, elaborations, or other modifications
      represent, as a whole, an original work of authorship. For the purposes
      of this License, Derivative Works shall not include works that remain
      separable from, or merely link (or bind by name) to the interfaces of,
      the Work and Derivative Works thereof.

      "Contribution" shall mean any work of authorship, including
      the original version of the Work and any modifications or additions
      to that Work or Derivative Works thereof, that is intentionally
      submitted to Licensor for inclusion in the Work by the copyright owner
      or by an individual or Legal Entity authorized to submit on behalf of
      the copyright owner. For the purposes of this definition, "submitted"
      means any form of electronic, verbal, or written communication sent
      to the Licensor or its representatives, including but not limited to
      communication on electronic mailing lists, source code control systems,
      and issue tracking systems that are managed by, or on behalf of, the
      Licensor for the purpose of discussing and improving the Work, but
      excluding communication that is conspicuously marked or otherwise
      designated in writing by the copyright owner as "Not a Contribution."

      "Contributor" shall mean Licensor and any individual or Legal Entity
      on behalf of whom a Contribution has been received by Licensor and
      subsequently incorporated within the Work.

   2. Grant of Copyright License. Subject to the terms and conditions of
      this License, each Contributor hereby grants to You a perpetual,
      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
      copyright license to reproduce, prepare Derivative Works of,
      publicly display, publicly perform, sublicense, and distribute the
      Work and such Derivative Works in Source or Object form.

   3. Grant of Patent License. Subject to the terms and conditions of
      this License, each Contributor hereby grants to You a perpetual,
      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
      (except as stated in this section) patent license to make, have made,
      use, offer to sell, sell, import, and otherwise transfer the Work,
      where such license applies only to those patent claims licensable
      by such Contributor that are necessarily infringed by their
      Contribution(s) alone or by combination of their Contribution(s)
      with the Work to which such Contribution(s) was submitted. If You
      institute patent litigation against any entity (including a
      cross-claim or counterclaim in a lawsuit) alleging that the Work
      or a Contribution incorporated within the Work constitutes direct
      or contributory patent infringement, then any patent licenses
      granted to You under this License for that Work shall terminate
      as of the date such litigation is filed.

   4. Redistribution. You may reproduce and distribute copies of the
      Work or Derivative Works thereof in any medium, with or without
      modifications, and in Source or Object form, provided that You
      meet the following conditions:

      (a) You must give any other recipients of the Work or
          Derivative Works a copy of this License; and

      (b) You must cause any modified files to carry prominent notices
          stating that You changed the files; and

      (c) You must retain, in the Source form of any Derivative Works
          that You distribute, all copyright, patent, trademark, and
          attribution notices from the Source form of the Work,
          excluding those notices that do not pertain to any part of
          the Derivative Works; and

      (d) If the Work includes a "NOTICE" text file as part of its
          distribution, then any Derivative Works that You distribute must
          include a readable copy of the attribution notices contained
          within such NOTICE file, excluding those notices that do not
          pertain to any part of the Derivative Works, in at least one
          of the following places: within a NOTICE text file distributed
          as part of the Derivative Works; within the Source form or
          documentation, if provided along with the Derivative Works; or,
          within a display generated by the Derivative Works, if and
          wherever such third-party notices normally appear. The contents
          of the NOTICE file are for informational purposes only and
          do not modify the License. You may add Your own attribution
          notices within Derivative Works that You distribute, alongside
          or as an addendum to the NOTICE text from the Work, provided
          that such additional attribution notices cannot be construed
          as modifying the License.

      You may add Your own copyright statement to Your modifications and
      may provide additional or different license terms and conditions
      for use, reproduction, or distribution of Your modifications, or
      for any such Derivative Works as a whole, provided Your use,
      reproduction, and distribution of the Work otherwise complies with
      the conditions stated in this License.

   5. Submission of Contributions. Unless You explicitly state otherwise,
      any Contribution intentionally submitted for inclusion in the Work
      by You to the Licensor shall be under the terms and conditions of
      this License, without any additional terms or conditions.
      Notwithstanding the above, nothing herein shall supersede or modify
      the terms of any separate license agreement you may have executed
      with Licensor regarding such Contributions.

   6. Trademarks. This License does not grant permission to use the trade
      names, trademarks, service marks, or product names of the Licensor,
      except as required for reasonable and customary use in describing the
      origin of the Work and reproducing the content of the NOTICE file.

   7. Disclaimer of Warranty. Unless required by applicable law or
      agreed to in writing, Licensor provides the Work (and each
      Contributor provides its Contributions) on an "AS IS" BASIS,
      WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
      implied, including, without limitation, any warranties or conditions
      of TITLE, NON-INFRINGEMENT, MERCHANTABILITY, or FITNESS FOR A
      PARTICULAR PURPOSE. You are solely responsible for determining the
      appropriateness of using or redistributing the Work and assume any
      risks associated with Your exercise of permissions under this License.

   8. Limitation of Liability. In no event and under no legal theory,
      whether in tort (including negligence), contract, or otherwise,
      unless required by applicable law (such as deliberate and grossly
      negligent acts) or agreed to in writing, shall any Contributor be
      liable to You for damages, including any direct, indirect, special,
      incidental, or consequential damages of any character arising as a
      result of this License or out of the use or inability to use the
      Work (including but not limited to damages for loss of goodwill,
      work stoppage, computer failure or malfunction, or any and all
      other commercial damages or losses), even if such Contributor
      has been advised of the possibility of such damages.

   9. Accepting Warranty or Additional Liability. While redistributing
      the Work or Derivative Works thereof, You may choose to offer,
      and charge a fee for, acceptance of support, warranty, indemnity,
      or other liability obligations and/or rights consistent with this
      License. However, in accepting such obligations, You may act only
      on Your own behalf and on Your sole responsibility, not on behalf
      of any other Contributor, and only if You agree to indemnify,
      defend, and hold each Contributor harmless for any liability
      incurred by, or claims asserted against, such Contributor by reason
      of your accepting any such warranty or additional liability.

   END OF TERMS AND CONDITIONS

   APPENDIX: How to apply the Apache License to your work.

      To apply the Apache License to your work, attach the following
      boilerplate notice, with the fields enclosed by brackets "[]"
      replaced with your own identifying information. (Don't include
      the brackets!)  The text should be enclosed in the appropriate
      comment syntax for the file format. We also recommend that a
      file or class name and description of purpose be included on the
      same "printed page" as the copyright notice for easier
      identification within third-party archives.

   Copyright [yyyy] [name of copyright owner]

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

```
