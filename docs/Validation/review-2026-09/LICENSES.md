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
| `Seqeron.Genomics.AlignmentConsensus` | GPL-3.0-only | 1 | 437 | `tests/Seqeron/Seqeron.Genomics.AlignmentConsensus.Tests` | (see the F-entry of the batch report) |
| `Seqeron.Genomics.ClustalW` | LGPL-3.0-or-later | 23 | 11,685 | `tests/Seqeron/Seqeron.Genomics.ClustalW.Tests` | ClustalW 2.1 (`ClustalWAligner`: pairwise / quick / profile / multiple stages, trees, command line, output formats, matrices) and the guide-tree progressive / iterative MSA family (`MultipleAlignProgressive`, `MultipleAlignIterative`, `ProgressiveMsaOptions`, guide trees, distances) |
| `Seqeron.Genomics.Kalign` | GPL-3.0-or-later | 4 | 3,115 | `tests/Seqeron/Seqeron.Genomics.Kalign.Tests` | Kalign 3.3.5 (`KalignAligner`, `KalignProgram`, MSF / Clustal / FASTA I/O) |
| `Seqeron.Genomics.PartialOrderAlignment` | GPL-2.0-or-later | 4 | 2,424 | `tests/Seqeron/Seqeron.Genomics.PartialOrderAlignment.Tests` | POA 2.0 (`PoaAligner`, `PoaProgram`, score matrices, PO / PIR / CLUSTAL output) |
| `Seqeron.Genomics.TCoffee` | GPL-2.0-or-later | 22 | 11,530 | `tests/Seqeron/Seqeron.Genomics.TCoffee.Tests` | T-Coffee 13.41 (`TCoffeeAligner`: libraries, extension, proba / pair methods, structure glue, tree, evaluate, program options) incl. the consistency aligner `MultipleAlignConsistency` (`ConsistencyMsaOptions`) and ProbCons 1.12 (`ProbConsAligner`, `ProbConsProgram`) |
| `Seqeron.Genomics.Emboss.SequenceFormats` | GPL-2.0-or-later (library files LGPL-2.1-or-later) | 107 | 50,551 | `tests/Seqeron/Seqeron.Genomics.Emboss.SequenceFormats.Tests` | EMBOSS 6.6.0 AJAX / NUCLEUS sequence, feature, report and alignment readers / writers (the lowest EMBOSS layer: ACD, USA, databases, DRCAT / EDAM data, deflate), seqret / union and the other sequence programs, transeq conventions |
| `Seqeron.Genomics.Emboss.FeatureTools` | GPL-2.0-or-later (library files LGPL-2.1-or-later) | 8 | 2,615 | `tests/Seqeron/Seqeron.Genomics.Emboss.FeatureTools.Tests` | EMBOSS report / feature-output formats and the feature programs (showfeat, extractfeat, twofeat, ...) |
| `Seqeron.Genomics.Emboss.Repeats` | GPL-2.0-or-later | 6 | 1,484 | `tests/Seqeron/Seqeron.Genomics.Emboss.Repeats.Tests` | EMBOSS wordmatch, seqmatchall, etandem, equicktandem |
| `Seqeron.Genomics.Emboss.PairwiseAlignment` | GPL-2.0-or-later | 23 | 7,252 | `tests/Seqeron/Seqeron.Genomics.Emboss.PairwiseAlignment.Tests` | EMBOSS needle, water, stretcher, matcher, supermatcher, needleall, merger, megamerger, est2genome, diffseq, dotmatcher, showalign, distmat, plotcon, wordfinder |
| `Seqeron.Genomics.Emboss.Motifs` | GPL-2.0-or-later (library files LGPL-2.1-or-later) | 11 | 4,123 | `tests/Seqeron/Seqeron.Genomics.Emboss.Motifs.Tests` | EMBOSS fuzznuc, fuzzpro, fuzztran, patmatmotifs, patmatdb, prosextract, marscan, consambig, the embedded PCRE |
| `Seqeron.Genomics.Emboss.OpenReadingFrames` | GPL-2.0-or-later | 6 | 1,555 | `tests/Seqeron/Seqeron.Genomics.Emboss.OpenReadingFrames.Tests` | EMBOSS getorf, plotorf, tcode, checktrans (the EMBOSS ORF conventions) |
| `Seqeron.Genomics.OptimalClustering` | LGPL-3.0-or-later | 1 | 638 | `tests/Seqeron/Seqeron.Genomics.OptimalClustering.Tests` | Ckmeans.1d.dp 4.3.6 (optimal one-dimensional k-means of cancer cell fractions, BIC choice of k); the MIT heterogeneity analysis has its own optimal 1-D k-means |
| `Seqeron.Genomics.AlleleSpecificCopyNumber` | GPL-3.0-only | 7 | 6,619 | `tests/Seqeron/Seqeron.Genomics.AlleleSpecificCopyNumber.Tests` | ASCAT 3.2.0 (ASPCF / asmultipcf segmentation, purity / ploidy fit, multiplicity), Battenberg 3.0.0 (phased-BAF segmentation, subclonal states, merging, callSubclones), the R Mersenne-Twister generator |
| `Seqeron.Genomics.PurityPeakQc` | GPL-3.0-or-later | 3 | 1,857 | `tests/Seqeron/Seqeron.Genomics.PurityPeakQc.Tests` | CNAqc 1.1.5 peak-based purity QC (analyze_peaks, subclonal expectations), BMix 1.0.0 binomial mixtures, R kmeans (Hartigan-Wong / MacQueen) |
| `Seqeron.Genomics.OverrepresentedWords` | AGPL-3.0-only | 6 | 2,988 | `tests/Seqeron/Seqeron.Genomics.OverrepresentedWords.Tests` | RSAT oligo-analysis v1.169, dyad-analysis v1.78, RSA.disco.lib (background models, word / dyad significance, shared-motif statistics) |
| `Seqeron.Genomics.OligoThermodynamics` | GPL-2.0-or-later | 13 | 8,390 | `tests/Seqeron/Seqeron.Genomics.OligoThermodynamics.Tests` | Primer3 / primer3-py 2.3.1: ntthal (thal.c hairpins, dimers, end stability), dpal.c alignment-mode complementarity, masker.c template masking, libprimer3 primer / probe picking (`ThermodynamicPrimerDesigner`, `ThermodynamicProbeDesigner`, `PrimerParameters`, `PrimerPairOptions`, `ProbeParameters`, mispriming libraries) |
| `Seqeron.Genomics.Last` | GPL-3.0-or-later | 8 | 3,215 | `tests/Seqeron/Seqeron.Genomics.Last.Tests` | LAST 1256 (`LastPrograms`: `lastdb` and `lastal` for DNA - YASS subset seed, adaptive seeds, gapless and gapped x-drop extension, the ambiguity scores of the matrix, the Karlin-Altschul parameters of ungapped alignments, TAB and MAF output; the programs that `mafft --lastpair` starts; no E-values: the ALP library is not ported) |

| MCP server | SPDX | Source files | Tool docs | Tests |
|---|---|---|---|---|
| `Seqeron.Mcp.AlignmentConsensus` | GPL-3.0-only | 3 | 1 (docs/mcp/tools/alignmentconsensus) | `tests/Seqeron/Seqeron.Mcp.AlignmentConsensus.Tests` |
| `Seqeron.Mcp.Last` | GPL-3.0-or-later | 3 | 1: mafft_command_last (docs/mcp/tools/last) | `tests/Seqeron/Seqeron.Mcp.Last.Tests` |
| `Seqeron.Mcp.MultipleAlignment` | GPL-3.0-or-later | 12 | 25 (docs/mcp/tools/multiplealignment) | `tests/Seqeron/Seqeron.Mcp.MultipleAlignment.Tests` |
| `Seqeron.Mcp.EmbossSequenceFormats` | GPL-2.0-or-later | 9 | 42 (docs/mcp/tools/embosssequenceformats) | `tests/Seqeron/Seqeron.Mcp.EmbossSequenceFormats.Tests` |
| `Seqeron.Mcp.EmbossFeatureTools` | GPL-2.0-or-later | 4 | 9 (docs/mcp/tools/embossfeaturetools) | `tests/Seqeron/Seqeron.Mcp.EmbossFeatureTools.Tests` |
| `Seqeron.Mcp.EmbossRepeats` | GPL-2.0-or-later | 4 | 3 (docs/mcp/tools/embossrepeats) | `tests/Seqeron/Seqeron.Mcp.EmbossRepeats.Tests` |
| `Seqeron.Mcp.EmbossPairwiseAlignment` | GPL-2.0-or-later | 13 | 16 (docs/mcp/tools/embosspairwisealignment) | `tests/Seqeron/Seqeron.Mcp.EmbossPairwiseAlignment.Tests` |
| `Seqeron.Mcp.EmbossMotifs` | GPL-2.0-or-later | 9 | 12 (docs/mcp/tools/embossmotifs) | `tests/Seqeron/Seqeron.Mcp.EmbossMotifs.Tests` |
| `Seqeron.Mcp.EmbossOpenReadingFrames` | GPL-2.0-or-later | 5 | 4 (docs/mcp/tools/embossopenreadingframes) | `tests/Seqeron/Seqeron.Mcp.EmbossOpenReadingFrames.Tests` |
| `Seqeron.Mcp.OverrepresentedWords` | AGPL-3.0-only | 4 | 3: oligo_analysis, dyad_analysis, shared_motifs_significance (docs/mcp/tools/overrepresentedwords) | `tests/Seqeron/Seqeron.Mcp.OverrepresentedWords.Tests` |
| `Seqeron.Mcp.OligoThermodynamics` | GPL-2.0-or-later | 4 | 11: design_primers, evaluate_primer, generate_primer_candidates, primer_dimer, design_probes, design_tiling_probes, design_antisense_probes, design_molecular_beacon, design_probes_primer3, validate_probe, analyze_oligo (docs/mcp/tools/oligothermodynamics) | `tests/Seqeron/Seqeron.Mcp.OligoThermodynamics.Tests` |

The MIT assemblies are listed in the root `LICENSE` file context: every other project under `src/` and `tests/` is MIT.

## 3. Inventory of the ported programs

### 3.1 Copyleft (separate assemblies)

| Program | Version | Licence (as the upstream states it) | Class | Where | Licence text taken from |
|---|---|---|---|---|---|
| EMBOSS (applications: needle, water, stretcher, matcher, supermatcher, needleall, merger, megamerger, seqmatchall, wordmatch, diffseq, est2genome, dotmatcher, showalign, infoalign, distmat, plotcon, consambig, fuzznuc, fuzzpro, fuzztran, patmatmotifs, patmatdb, prosextract, etandem, equicktandem, getorf, tcode, checktrans, plotorf, transeq, seqret and the other sequence / feature programs) | 6.6.0 | GNU General Public License, version 2 or (at your option) any later version (file headers of the application sources, e.g. emboss/needle.c, dotmatcher.c) | strong copyleft (GPL-2.0-or-later) | SEPARATED (WP80, F364): `Seqeron.Genomics.Emboss.SequenceFormats` / `.FeatureTools` / `.Repeats` / `.PairwiseAlignment` / `.Motifs` / `.OpenReadingFrames` (by subject area, GPL-2.0-or-later) and the MCP servers `Seqeron.Mcp.Emboss*` (six, 86 tools) | EMBOSS 6.6.0 source (header of each file) and /usr/share/doc/emboss/copyright (Debian emboss 6.6.0+dfsg-12ubuntu2) |
| EMBOSS libraries AJAX and NUCLEUS (sequence / feature / report readers and writers, ACD, USA, databases, alignment output formats, PCRE, zlib deflate) | 6.6.0 | GNU Lesser General Public License, version 2.1 or (at your option) any later version (headers of ajax/*.c, nucleus/*.c); the embedded PCRE (University of Cambridge licence) and the zlib deflate (zlib licence) keep their own terms | weak copyleft (LGPL-2.1-or-later) | SEPARATED (WP80, F364): the library files (per-file SPDX LGPL-2.1-or-later) live in the same Emboss assemblies, chiefly `Seqeron.Genomics.Emboss.SequenceFormats` and `.FeatureTools`; the LGPL version is read from the SPDX / file headers of the EMBOSS 6.6.0 sources (2.1 or later) | EMBOSS 6.6.0 source headers; Debian copyright file (License: emboss = AJAX / NUCLEUS LGPL, applications GPL; PCRE_license; zlib) |
| T-Coffee | 13.41.0 (Debian 13.45.0.4846264+really13.41.0.28bdc39+dfsg-1 sources) | GNU General Public License, version 2 or (at your option) any later version (file headers, Debian copyright: GPL-2+) | strong copyleft (GPL-2.0-or-later) | Seqeron.Genomics.TCoffee | T-Coffee source file headers; /usr/share/doc/t-coffee/copyright |
| ProbCons (hosted in the T-Coffee assembly because the port runs on the pair-HMM tables embedded in T-Coffee) | 1.12 | Public domain (ProbCons README: 'PROBCONS has been made freely available as PUBLIC DOMAIN software') | public domain | Seqeron.Genomics.TCoffee (ProbConsAligner*, ProbConsProgram) | /usr/share/doc/probcons/copyright |
| Kalign | 3.3.5 | GNU General Public License, version 3 or (at your option) any later version | strong copyleft (GPL-3.0-or-later) | Seqeron.Genomics.Kalign | Kalign 3.3.5 source (COPYING) and /usr/share/doc/kalign/copyright |
| LAST (lastdb, lastal) | 1256 (gitlab.com/mcfrith/last; the version of the oracle) | GNU General Public License, version 3 or (at your option) any later version (README: "LAST is distributed under the GNU General Public License, either version 3 of the License, or (at your option) any later version"; COPYING.txt = GPL-3; newer files carry SPDX-License-Identifier: GPL-3.0-or-later; Debian copyright of last-align: GPL-3+) | strong copyleft (GPL-3.0-or-later) | Seqeron.Genomics.Last (WP85, F376). The MAFFT side of `--lastpair` (`pairlocalalign -R`, `readlastresx`, `block2reg`, `putlocalhom_last`) is MAFFT code (BSD-style MAFFT licence) and stays in the MIT assembly, which only knows the interface `ILastPrograms` | LAST 1256 source (README.rst, COPYING.txt, file headers) and /usr/share/doc/last-align/copyright |
| ClustalW | 2.1 (Debian clustalw 2.1+lgpl-7) | GNU Lesser General Public License, version 3 or (at your option) any later version | weak copyleft (LGPL-3.0-or-later) | Seqeron.Genomics.ClustalW (also hosts the guide-tree progressive / iterative MSA family built on ClustalW's distances and trees) | /usr/share/doc/clustalw/copyright (License: LGPL-3+) |
| POA (Partial Order Alignment) | 2.0 | GNU General Public License, version 2 or (at your option) any later version | strong copyleft (GPL-2.0-or-later) | Seqeron.Genomics.PartialOrderAlignment | POA 2.0 source (COPYING) and /usr/share/doc/poa/copyright |


### 3.2 Permissive and public domain (kept in the MIT assemblies, NOTICE)

| Program | Version | Licence | Class | Where | Licence text taken from |
|---|---|---|---|---|---|
| MAFFT | 7.505 | MAFFT licence (BSD-style: retain the notice, reproduce it in binary distributions, no endorsement) | permissive | NOTICE (THIRD-PARTY-NOTICES.md) | /usr/share/doc/mafft/copyright |
| CONTRAfold (inference: inside / outside, posterior base-pairing probabilities, default parameters) | 2.02 | BSD 3-clause ("Copyright (c) 2006, Chuong Do"; manual.tex of contrafold_v2_02.tar.gz, section "License (BSD)"; contra.stanford.edu) | permissive | `Seqeron.Genomics.Alignment` (`ContrafoldProgram.cs`; WP86, F377), NOTICE | manual.tex of contrafold_v2_02.tar.gz (sha256 9e265726...) |
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
| primer3 / primer3-py `thal.c` (ntthal) | GNU GPL version 2 (LICENSE of primer3-org/primer3 and libnano/primer3-py) | `MolTools/NtthalDimer.cs`, `NtthalHairpin.cs`: "faithful port" / "bit-faithful port of primer3-py 2.3.1 thal.c" (PRIMER-TM-001-SPECIAL-LOOP-Evidence) | strong copyleft | SEPARATED (WP80, F368): `Seqeron.Genomics.OligoThermodynamics` (GPL-2.0-or-later; source headers of thal.c / dpal.c: "either version 2 of the License, or (at your option) any later version") with the MCP server `Seqeron.Mcp.OligoThermodynamics` (11 tools); it also holds the dpal.c / masker.c / libprimer3 ports found in `PrimerDesigner` / `ProbeDesigner` (`PrimerDesigner.Masking.cs`, `.MisprimingLibrary.cs`, `.BoundAndPosition.cs`, `.SequenceQuality.cs`, `.TemplateMispriming.cs`, the Primer3 picker and the structure screens). The MIT `PrimerDesigner` keeps the nearest-neighbour Tm / GC / dimer / hairpin screens written from the published models |
| DECIPHER 3.9.4 `ConsensusSequence` | GPL-3 (DESCRIPTION of the Bioconductor package) | `Analysis/MotifFinder.DecipherConsensus.cs`: "a line-by-line port of DECIPHER 3.9.4" | strong copyleft | SEPARATED: `Seqeron.Genomics.AlignmentConsensus` (GPL-3.0-only), MCP `Seqeron.Mcp.AlignmentConsensus` |
| RSAT `oligo-analysis` v1.169, `dyad-analysis` v1.78, `RSA.disco.lib` (rsa-tools/rsat-code 10043f2) | GNU Affero General Public License v3 (LICENSE.txt of rsa-tools/rsat-code) | `Analysis/MotifFinder.OligoAnalysis*.cs`, `MotifFinder.DyadAnalysis.cs`, `MotifFinder.SharedMotifsDegenerate.cs`, `MotifFinder.OligoSequenceType.cs`: "Ported from rsa-tools/rsat-code" | strong copyleft (AGPL-3.0-only, network clause) | SEPARATED (WP80, F367): `Seqeron.Genomics.OverrepresentedWords` (AGPL-3.0-only) and the MCP server `Seqeron.Mcp.OverrepresentedWords` (oligo_analysis, dyad_analysis, shared_motifs_significance); LICENSE.txt of rsa-tools/rsat-code = AGPL v3 (the perl scripts carry no per-file licence header) |
| Ckmeans.1d.dp (R package) | LGPL (>= 3) (DESCRIPTION) | `Oncology/OncologyAnalyzer.Clonality.cs`: "a line-by-line port" of the optimal 1-D k-means | weak copyleft | SEPARATED (WP80, F365): `Seqeron.Genomics.OptimalClustering` (LGPL-3.0-or-later); the MIT `AnalyzeHeterogeneity` uses an independently written optimal 1-D k-means (Bellman / Fisher dynamic program) |
| ASCAT / Battenberg (VanLoo-lab) | GNU GPL v3 (LICENSE of VanLoo-lab/ascat; Battenberg is GPL-3 as well) | `Oncology/OncologyAnalyzer.CopyNumberPloidy.cs`: ASCAT `ascat.aspcf`, `ascat.runAscat`, `determine_copynumber`, Battenberg clonal-state arithmetic ported function by function (ONCO-ASCAT-001-Evidence, F58 ...) | strong copyleft | SEPARATED (WP80, F366): `Seqeron.Genomics.AlleleSpecificCopyNumber` (GPL-3.0-only: DESCRIPTION of ASCAT 3.2.0 and Battenberg 3.0.0 "License: GPL-3"). CNAqc 1.1.5 ("License: GPL (>= 3) + file LICENSE"), BMix 1.0.0 ("GPL (>= 3)") and the R `kmeans` port found in the same file: `Seqeron.Genomics.PurityPeakQc` (GPL-3.0-or-later); the R Mersenne-Twister generator is in the ASCAT assembly |
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
| LAST (Docker package handoff/A2-19c-last; MAFFT --lastpair; run 2 builds 1256 / 1250 / 1282 / 1411 / 1542 from gitlab.com/mcfrith/last and selects the first fast one, named in the returned out/versions.txt; the audit machine has 1542) | last-align 1542-1 |
| CONTRAfold 2.02 (built in the user Docker package handoff/A2-19c-rna, `-O3 -DNDEBUG`) and DAFS 0.0.3 (built from the public release; run as a black box: `mafft --contrafold` / `--dafs` runs and probes recorded; the CONTRAfold inference is ported, see 3.2) | contrafold_v2_02.tar.gz (contra.stanford.edu); satoken/dafs v0.0.3 (github.com) |
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



### 5.1 FASTA 3.4: where MAFFT's `fasta34` comes from (facts established in WP80)

* The `mafft` script (Debian 7.505-1, `/usr/bin/mafft` lines 1630-1642) does `FASTA_4_MAFFT=$(which fasta34)` and stops with "This option requires the fasta34 program (FASTA version x.xx or higher) installed in your PATH" when it is missing. **MAFFT does not contain or build `fasta34`**: `dpkg -L mafft` lists only `mafft`, its helper binaries under `/usr/lib/mafft/lib/mafft` (disttbfast, dndfast7, splittbfast, ...) and two ruby scripts; Debian `mafft` Depends: libc6, perl (Recommends blast2, libwww-perl, lynx, ruby); its copyright file covers MAFFT's own BSD-style licence only.
* The Debian package `fasta3` (36.3.8i.14-Nov-2020-1, installed in this environment; `dpkg -S /usr/bin/fasta36` = fasta3) installs `fasta36` and **no `fasta34`**; its copyright file says "License: Apache-2.0" (Pearson / University of Virginia, 1996-2014) for the FASTA 36 sources. So the Apache-2.0 sources are NOT what `fasta34` is.
* The `fasta34` of the oracle runs (WP75, F346-F349; ALIGN-MULTI-001-Evidence "WP75", handoff `A2-19c-fastapair`) was built by the user from **`fasta-34.26.5.shar.Z` = the FASTA 3.4 sources** (binary sha256 `7491909b...`), i.e. the sources whose `COPYRIGHT` file restricts commercial use (quoted in the table above). The archive was only parsed, never executed, by the review.
* The source files that the port's headers and the Evidence cite (all of FASTA 3.4, release 34.26.5): `comp_lib.c`, `dropnfa.c` (`init_work`, `do_fasta`, `savemax`, `spam`, `sconn`, `dmatch` / `FLOCAL_SIGN` ...), `mshowbest.c`, `scaleswn.c` (`proc_hist_a`, `find_za`, `zs_to_E`, `scale_scores`), `compacc.c` (`sortbeste`), `alt_parms.h`, `apam.c`, `upam.h` (BLOSUM50), `uascii.h`, `initfa.c`, `doinit.c`, `mmgetaa.c` (header of `SequenceAligner.MafftFasta34.cs`; Evidence "Sources" list). `smith_waterman_sse2.*` (Farrar, with its own "may not be sold" wording) is listed in the table row as a further restricted file of that tree.
* Consequence: the ported code is derived from FASTA **3.4** sources with the restricted COPYRIGHT, not from the Apache-2.0 FASTA 36 sources; the later FASTA 36 release (Apache-2.0) contains successors of the same files but the port and its oracle (`fasta34` 34.26.5 output, `-z3 -m10 -b500000 ...`) have not been compared with them.

### 5.2 Options per item (nothing below was done; each option names its consequence)

| Item | keep + notice | separate assembly | remove | rewrite |
|---|---|---|---|---|
| FASTA 3.4 (`SequenceAligner.MafftFasta34*.cs`, `MafftProgram*` paths using `--fastapair` / `--fastswpair` / `--fastaparttree`) | keep: the MIT package then carries code under a no-commercial-use licence - the notice states it, but "MIT" is no longer true for the whole package | `Seqeron.Genomics.FastaPearson`-style assembly under the FASTA licence (not GPL): the MIT MAFFT code would lose `dndfast7`-based options (or take a hook as for ASCAT) | drop the three MAFFT options and their oracle tests (about 4 F entries of tests) | re-derive from the Apache-2.0 FASTA 36 sources and re-verify against the `fasta34` oracle: probable but unproven equality (scaling, `-z 3`, `savemax`) |
| MUMmer 3.23 (Artistic 1.0: modified versions marked, standard version available) | keep: add the marking sentence and the Artistic text to NOTICE | separate assembly under Artistic 1.0 (neither MIT nor GPL) - the MIT `AnchorBasedAligner` stays | drop `Mummer*.cs`, `Nucmer*.cs`, `Promer*.cs` and the MCP tools | rewrite from the published papers (MUM / anchor chaining), results would differ from the oracle |
| Mash Bloom filter (`MashBloomFilter.cs`, CPL-1.0) | keep + notice: CPL requires that changes to the filter source are released under the CPL | separate assembly under CPL-1.0 (the `mash dist -b` / screen option hooks to it) | remove `-b` / the Bloom-filter screen | rewrite the Bloom filter (a textbook structure); the false-positive pattern of the oracle would need re-checking |
| BBTools (BBDuk entropy masking) | keep + notice: BSD-3 text + the "Enhancements" paragraph (improvements must be made available to the licensor) | separate assembly under the BBTools licence | remove `maskLowEntropy` | rewrite from the published entropy definition |
| TargetScan 7 scripts (`MiRnaAnalyzer.cs`: PCT arithmetic taken verbatim; no LICENSE file in the mirror) | keep + notice and ask the authors for written terms (the code stays MIT with an unknown-licence snippet) | separate assembly: needs a licence to name first | remove the PCT term of the context score | rewrite the PCT arithmetic from Friedman et al. 2009 / Agarwal et al. 2015 and keep a differential test against the verbatim values
| REBASE `allenz` columns (organism, methylation; Copyright Richard J. Roberts, redistribution question of B08) | keep + notice and ask REBASE for permission | data in a separate package | ship only the enzymes table derived from Biopython (sites / cuts) and drop the organism / methylation columns | rebuild organism / methylation from the primary literature
| parasail / PCRE wording | keep + notice (BSD-style / PCRE acknowledgement sentence already in NOTICE) | n/a (PCRE sits in the EMBOSS Motifs assembly) | n/a | n/a |
| `LibstdcxxSort.cs` (GPL-3 + runtime exception) / `GlibcRandom.cs` (LGPL-2.1+) | keep: the files reproduce published algorithms (introsort; r[i] = r[i-3] + r[i-31]) and not the library text | move to the Alignment copyleft siblings if the owner reads them as derived | replace by `Array.Sort` / `System.Random` (the MUMmer / minimap2 tie order and the `rand()` driven oracles then differ) | already independent: re-implement from Musser 1997 / Knuth if wording is a concern |
| GISTIC2 ziggurat deconstruction (`OncologyAnalyzer.CopyNumberPloidy.cs`, ported from the MATLAB sources `perform_deconstruction.m`, `atomic_zigg_deconstruction.m`, `normalize_by_arm_length.m`; found in WP80, licence NOT verified offline) | keep + notice once the Broad licence of GISTIC2 is read | separate assembly named by meaning (copy-number deconstruction) | remove the `ziggurat` features (F50-F52) | rewrite from Mermel et al. 2011 and the published ziggurat description |
| R-reproducing numerics in MIT `StatisticsHelper` (`BandwidthNrd0`, `GaussianKernelDensity` = R `density`, `SampleQuantileType7`, `RMean`, `RSampleVariance`, `ExtendedPrecision*`, `BrentRoot`, `NormalQuantile`, `Digamma`, `RegularizedIncompleteBeta`: written from the documented formulas, checked against R; `PeakPick` follows peakPick 0.11, MIT) | keep: formulas of R `stats` are not copied code (R is GPL-2 / GPL-3) | n/a | n/a | n/a |

## 6. Status of the separation (what is done, what is registered)

WP79 (commit "B09 WP79 ..."): Kalign, POA, ClustalW (+ the guide-tree progressive / iterative MSA family), T-Coffee (+ the consistency aligner and ProbCons), DECIPHER ConsensusSequence are separate assemblies with LICENSE files, SPDX headers, tests and MCP servers.
WP80 (this work package), all done and green:

* F364 EMBOSS 6.6.0: six assemblies by subject area (`Seqeron.Genomics.Emboss.SequenceFormats`, `.FeatureTools`, `.Repeats`, `.PairwiseAlignment`, `.Motifs`, `.OpenReadingFrames`), six test projects and six MCP servers (86 tools, same names). "Emboss" is the identity of the derived work (precedent: Kalign, ClustalW, T-Coffee); the area suffix carries the meaning; the sequence / feature / report / alignment libraries are in the lowest assembly because programs of every area need the writers. The MIT seams (translator conventions, TopN tie order, column marks, infoalign / distmat / alignment statistics, string pattern overloads, `AlignmentTools.ExternalFormatWriter` hook for the EMBOSS alignment formats) were cut.
* F365 Ckmeans.1d.dp: `Seqeron.Genomics.OptimalClustering` (LGPL-3.0-or-later) + tests; the MIT heterogeneity analysis has its own optimal 1-D k-means.
* F366 ASCAT 3.2.0 / Battenberg 3.0.0 / CNAqc 1.1.5 / BMix 1.0.0 / R kmeans: `Seqeron.Genomics.AlleleSpecificCopyNumber` (GPL-3.0-only) and `Seqeron.Genomics.PurityPeakQc` (GPL-3.0-or-later) + tests (no MCP tools exist for them; the oncology API is C#-only).
* F367 RSAT oligo-analysis / dyad-analysis: `Seqeron.Genomics.OverrepresentedWords` (AGPL-3.0-only) + `Seqeron.Mcp.OverrepresentedWords` (3 tools) + tests; the MIT exact PWM p-value takes an explicit `PwmBackgroundChain` (own `FromWordFrequencies` builder, so `pwm_score_pvalue` stays on the MIT server).
* F368 Primer3 (ntthal / dpal / masker / libprimer3): `Seqeron.Genomics.OligoThermodynamics` (GPL-2.0-or-later) + `Seqeron.Mcp.OligoThermodynamics` (11 tools: design_primers, evaluate_primer, generate_primer_candidates, primer_dimer, design_probes, design_tiling_probes, design_antisense_probes, design_molecular_beacon, design_probes_primer3, validate_probe, analyze_oligo) + tests. The MIT `PrimerDesigner` / `ProbeDesigner` keep the nearest-neighbour Tm, GC, simple dimer / hairpin and probe screens.
* F376 (WP85) LAST 1256: `Seqeron.Genomics.Last` (GPL-3.0-or-later) + `Seqeron.Genomics.Last.Tests` (the 724 recorded `lastal` calls and the 204 `mafft --lastpair` runs of the user's Docker package). The MIT `Seqeron.Genomics.Alignment` has the MAFFT side and the interface `ILastPrograms` and does not reference the GPL assembly: a caller gives an installation of `mafft` an implementation (`MafftInstallation.Last`); without one the script stops with its "Install LAST" text, so the MIT MCP tool `mafft_command` is unchanged.
* `LicenceSeparationTests` covers every new assembly automatically (auto-discovery by `PackageLicenseExpression`), and has new tests: the MIT MCP servers expose no copyleft tool name, and the expected copyleft assemblies exist by name.

Shared MIT building blocks made `public` for the copyleft assemblies carry XML documentation; a few helpers got argument guards because the reflection fuzz tests enumerate public static methods.

Still registered (owner decisions, section 5): FASTA 3.4, MUMmer, Mash Bloom filter, BBTools, TargetScan, REBASE allenz, parasail / PCRE wording, `LibstdcxxSort` / `GlibcRandom`, GISTIC2 (new, licence not verified), R-reproducing numerics. Nothing of section 5 was moved or removed.


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
