# Review provenance & re-review queue — campaign review-2026-09

> **Generated** by `tools/validation/review_provenance.py` — do not hand-edit.
>
> Every unit below was reviewed by an autonomous agent session (Stage A + Stage B per
> `docs/Validation/VALIDATION_PROTOCOL.md`), not by a human reviewer. This table exists so the
> whole campaign can be **re-run selectively** once the review tooling is reworked: take the
> units whose `Re-review` box is unchecked (or all of them, per batch) and replay the protocol.
> The reviewing agent and its session are recorded *per commit* on the campaign branch:
> `git log --format='%h %s %(trailers:key=Co-Authored-By,valueonly) %(trailers:key=Claude-Session,valueonly)'`.

Campaign start: `baa72cf06~1`. Commits scanned: 224. Units: 258.

| Unit | Batch | Stage A | Stage B | State | Review commits | Re-review |
|---|---|---|---|---|---|---|
| SEQ-GC-001 | B01 | 🟡 PASS-WITH-NOTES | ✅ PASS | ✅ CLEAN | `fb7e448b4` | ☐ |
| SEQ-COMP-001 | B01 | ✅ PASS | 🟡 PASS-WITH-NOTES | ✅ CLEAN | `abf556998` | ☐ |
| SEQ-REVCOMP-001 | B01 | 🟡 PASS-WITH-NOTES | ✅ PASS | ✅ CLEAN | `029231323` | ☐ |
| SEQ-VALID-001 | B01 | 🟡 PASS-WITH-NOTES | 🟡 PASS-WITH-NOTES | ✅ CLEAN | `3ce8dd469` | ☐ |
| SEQ-COMPLEX-001 | B01 | 🟡 PASS-WITH-NOTES | ✅ PASS | ✅ CLEAN | `e8c076e48` `6c4e882df` | ☐ |
| SEQ-ENTROPY-001 | B01 | 🟡 PASS-WITH-NOTES | ❌ FAIL → fixed | ✅ CLEAN | `cc4787160` | ☐ |
| SEQ-GCSKEW-001 | B01 | 🟡 PASS-WITH-NOTES | ❌ FAIL → fixed | ✅ CLEAN | `2b2088c86` | ☐ |
| SEQ-ATSKEW-001 | B01 | ✅ PASS | ✅ PASS | ✅ CLEAN | `7959bf735` | ☐ |
| SEQ-GC-ANALYSIS-001 | B01 | 🟡 PASS-WITH-NOTES | ❌ FAIL → fixed | ✅ CLEAN | `ef4011485` | ☐ |
| SEQ-REPLICATION-001 | B01 | 🟡 PASS-WITH-NOTES | 🟡 PASS-WITH-NOTES | ✅ CLEAN | `05b94361e` | ☐ |
| SEQ-RNACOMP-001 | B01 | ✅ PASS | ✅ PASS | ✅ CLEAN | `0d6901259` | ☐ |
| TRANS-CODON-001 | B02 | PASS-WITH-NOTES | FAIL→fixed | FIXED | `53ef1d76b` | ☐ |
| TRANS-PROT-001 | B02 | PASS-WITH-NOTES | FAIL→fixed | FIXED | `5052a4c2a` | ☐ |
| TRANS-SIXFRAME-001 | B02 | PASS-WITH-NOTES | FAIL→fixed | FIXED | `32158911a` | ☐ |
| CODON-CAI-001 | B02 | PASS-WITH-NOTES | FAIL→fixed | FIXED | `be84d5132` | ☐ |
| CODON-ENC-001 | B02 | PASS-WITH-NOTES | FAIL→fixed | FIXED | `c24dc9ea8` | ☐ |
| CODON-OPT-001 | B02 | PASS-WITH-NOTES (docs corrected) | FAIL→fixed | FIXED | `e02301bd7` | ☐ |
| CODON-RARE-001 | B02 | PASS-WITH-NOTES | FAIL→fixed | FIXED (1 open item) | `26e3c3479` | ☐ |
| CODON-RSCU-001 | B02 | PASS-WITH-NOTES | FAIL→fixed | FIXED | `1f3508045` | ☐ |
| CODON-STATS-001 | B02 | PASS-WITH-NOTES | FAIL→fixed | FIXED | `141521133` | ☐ |
| CODON-USAGE-001 | B02 | PASS-WITH-NOTES | FAIL→fixed | FIXED | `dd5ee2294` | ☐ |
| SEQ-STATS-001 | B03 | — | — | not reviewed | `dd9067805` | ☐ |
| SEQ-COMPOSITION-001 | B03 | — | — | not reviewed | `dd9067805` | ☐ |
| SEQ-MW-001 | B03 | PASS-with-notes | FAIL→FIXED | FIXED | `ceb6a35d3` | ☐ |
| SEQ-PI-001 | B03 | FAIL→FIXED | FAIL→FIXED | FIXED | `92f55a0fb` | ☐ |
| SEQ-HYDRO-001 | B03 | PASS-with-notes | FAIL→FIXED | FIXED | `8220b7c0b` | ☐ |
| SEQ-THERMO-001 | B03 | — | — | not reviewed | `654629fcf` | ☐ |
| SEQ-TM-001 | B03 | — | — | not reviewed | `654629fcf` | ☐ |
| SEQ-DINUC-001 | B03 | PASS-with-notes | PASS (+F13) | FIXED (dedup) | `d663ba4b1` | ☐ |
| SEQ-SECSTRUCT-001 | B03 | FAIL→FIXED | PASS (profile) + new method | FIXED | `4e587193b` | ☐ |
| SEQ-CODON-FREQ-001 | B03 | PASS-with-notes | FAIL→FIXED | FIXED | `f92e88697` | ☐ |
| SEQ-ENTROPY-PROFILE-001 | B03 | PASS-with-notes | FAIL→FIXED | FIXED (dedup) | `0acd18e87` | ☐ |
| SEQ-GC-PROFILE-001 | B03 | PASS-with-notes | FAIL→FIXED | FIXED (dedup) | `0acd18e87` | ☐ |
| SEQ-SUMMARY-001 | B03 | PASS-with-notes | FAIL→FIXED | FIXED | `92f0dac1e` | ☐ |
| SEQ-COMPLEX-COMPRESS-001 | B04 | — | — | not reviewed | `a699874c7` | ☐ |
| SEQ-COMPLEX-DUST-001 | B04 | FAIL → corrected | FAIL → fixed | FIXED | `c1e8e97d9` `6f78c824b` `569961111` | ☐ |
| SEQ-COMPLEX-KMER-001 | B04 | PASS-WITH-NOTES (citation corrected) | PASS | CLEAN (dedup + docs) | `d3a45f7a0` `d05364ddd` `5713c7d47` `bb9b23b28` | ☐ |
| SEQ-COMPLEX-WINDOW-001 | B04 | — | — | not reviewed | `5713c7d47` `1e70d7c6e` | ☐ |
| REP-STR-001 | B04 | PASS-WITH-NOTES (reporting convention sourced: MISA, pytrf/Krait, maximal repetitions) | FAIL → fixed | FIXED | `d3a45f7a0` `f1abd8660` `4176b592a` | ☐ |
| REP-TANDEM-001 | B04 | — | — | not reviewed | `d3a45f7a0` `bb4742e4b` `a05770adb` | ☐ |
| REP-INV-001 | B04 | PASS-WITH-NOTES (reporting rule sourced: EMBOSS palindrome; einverted scoring not promised) | FAIL → fixed | FIXED | `105bc08ba` `850e996bc` | ☐ |
| REP-DIRECT-001 | B04 | FAIL → corrected (reporting convention sourced: MUMmer `repeat-match` maximal pairs, Gusfield 1997 §7.12, REPuter) | FAIL → fixed | FIXED | `d05364ddd` `f1abd8660` `aa32cf499` `8ce6a4d86` | ☐ |
| REP-PALIN-001 | B04 | PASS-WITH-NOTES (REVP "every reverse palindrome, length 4–12, any order" re-sourced; ACGT-only alphabet added; TestSpec TATA/ATAT label swap fixed) | FAIL → fixed | FIXED | `c3cc345dc` | ☐ |
| REP-APPROX-001 | B04 | FAIL → corrected (TRF 4.10.0 compiled as oracle: % stats between ADJACENT copies, copy number = aligned consensus columns, local wraparound DP, N never matches, overlapping periods reported; sum-of-heads table re-derived exactly) | FAIL → fixed | FIXED | `6a2140638` `210c1389c` `92e85db36` `5b6f73986` `ffe2267a3` `77b758499` | ☐ |
| PAT-EXACT-001 | B05 | PASS-WITH-NOTES | PASS-WITH-NOTES | FIXED | `802099bd5` | ☐ |
| PAT-APPROX-001 | B05 | PASS | PASS-WITH-NOTES | CLEAN | `f37dded8e` | ☐ |
| PAT-APPROX-002 | B05 | PASS | PASS-WITH-NOTES | FIXED | `72efb3066` | ☐ |
| PAT-APPROX-003 | B05 | PASS | FAIL→FIXED | FIXED | `fcfde7ead` | ☐ |
| PAT-IUPAC-001 | B05 | PASS | PASS-WITH-NOTES | CLEAN | `43c39ba69` | ☐ |
| PAT-PWM-001 | B05 | PASS-WITH-NOTES | PASS-WITH-NOTES | FIXED | `4021a5ce4` | ☐ |
| MOTIF-CONS-001 | B05 | PASS-WITH-NOTES | PASS-WITH-NOTES | FIXED | `ca5910b18` | ☐ |
| MOTIF-DISCOVER-001 | B05 | PASS-WITH-NOTES | FAIL→FIXED | FIXED | `2a9d8e666` | ☐ |
| MOTIF-SHARED-001 | B05 | PASS-WITH-NOTES | PASS-WITH-NOTES | FIXED | `82ad132b3` | ☐ |
| MOTIF-REGULATORY-001 | B05 | PASS-WITH-NOTES | FAIL→FIXED | FIXED | `1950c3f20` | ☐ |
| MOTIF-GENERATE-001 | B05 | PASS-WITH-NOTES | FAIL→FIXED | FIXED | `a0b526961` | ☐ |
| KMER-COUNT-001 | B06 | PASS-with-notes | PASS-with-notes | ☑ reviewed; F13 (audit round 1, WP4); F26 (audit round 2, WP8) | `f9be81380` `e0d4e473b` | ☐ |
| KMER-FREQ-001 | B06 | PASS-with-notes | PASS | ☑ reviewed; F12 (audit round 1, WP3); F26 (audit round 2, WP8) | `244520ea3` | ☐ |
| KMER-FIND-001 | B06 | PASS-with-notes | PASS-with-notes | ☑ reviewed; F1 (efficiency) fixed; F26 (audit round 2, WP8) | `f06151b71` | ☐ |
| KMER-ASYNC-001 | B06 | PASS-with-notes | PASS-with-notes | ☑ reviewed; F2, F3 fixed; F14 (audit round 1, WP4) | `3488d2ea8` | ☐ |
| KMER-BOTH-001 | B06 | PASS-with-notes | PASS | ☑ reviewed, no defect; F25 (audit round 2, WP8) | `6228d165f` | ☐ |
| KMER-DIST-001 | B06 | PASS-with-notes | PASS | ☑ reviewed; F11 (WP2), F15 (WP4), F19–F22 (WP6), F23–F24 (WP7), F27 (WP8), F28–F31 (WP9), F32–F36 (WP10), F37–F38 (WP11) | `1d9d04daf` `5cf05a992` `f7b4d57df` `7d8e4688f` `da9c8f535` `c74372f1b` | ☐ |
| KMER-GENERATE-001 | B06 | PASS-with-notes | PASS-with-notes | ☑ reviewed; F4 (robustness) fixed | `5ef210c6d` | ☐ |
| KMER-POSITIONS-001 | B06 | PASS-with-notes | PASS-with-notes | ☑ reviewed; F5 (efficiency) fixed | `a502f56a6` | ☐ |
| KMER-STATS-001 | B06 | PASS-with-notes | PASS-with-notes | ☑ reviewed; F6, F7, F8 fixed | `e6498131f` | ☐ |
| KMER-UNIQUE-001 | B06 | PASS-with-notes | PASS-with-notes | ☑ reviewed; F9 fixed | `9cbf226fd` | ☐ |
| PRIMER-TM-001 | B07 | ⚠ (salt formula mis-described, Owczarzy mis-cited) → fixed | ⚠ → fixed | FIXED | `a9f618373` | ☐ |
| PRIMER-DESIGN-001 | B07 | ⚠ (Primer3 Tm window applied to Marmur–Doty Tm; greedy pair selection declared although real algorithm implementable; `target_end` documented inclusive) → fixed | ❌ → fixed | FIXED | `4147d446a` | ☐ |
| PRIMER-STRUCT-001 | B07 | ⚠ (dimer model described as "terminal complementarity" but implemented as parallel complementarity; 3′-stability "<5 nt → 0" and N handling not Primer3; heuristic hairpin/dimer screens declared as the structure screen although Primer3's thermodynamic one is implementable) → fixed | ❌ → fixed | FIXED | `0ce1d25a4` | ☐ |
| PRIMER-NNTM-001 | B07 | ⚠ (NN4 = SantaLucia & Hicks 2004 described as "identical to SantaLucia 1998"/"SantaLucia (1998) unified"; mismatch path described as mirroring Biopython `Tm_NN` although it had no terminal-mismatch table; TestSpec "CLEAN" + wrong fixture path; basic Tm documented as ignoring U) → fixed | ❌ → fixed | FIXED | `270918321` | ☐ |
| PRIMER-HAIRPIN-001 | B07 | ⚠ (ntthal docs: "Not implemented: divalent salt", `saltCorrection = 0.368·ln([Na⁺])`, validated only vs primer3-py 2.3.0 at dv = 0; legacy folder doc pointed to UNAFold for loop terms that the ntthal path implements; MCP `evaluate_primer` described a "hairpin potential" although it reports Primer3 ntthal Tm values) → fixed | ❌ → fixed | FIXED | `dbd90d50e` | ☐ |
| PRIMER-DIMER-001 | B07 | ⚠ (TestSpec claimed ΔG ≤ 0 or no structure; docs said divalent not modelled / ANY only / "≥ 2 contiguous bp" / C_T 0.5 µM; validated only vs primer3-py 2.3.0 at dv = 0) → fixed | ❌ → fixed | FIXED | `4fbf8d386` | ☐ |
| PROBE-DESIGN-001 | B07 | ⚠ (probe Tm documented as "Wallace / salt-adjusted" helper with no stated conditions although presets cite NN-based tools; qPCR preset Tm 68–72 / GC 40–60 vs ABI 68–70 / 30–80; self-structure only sequence heuristics although Primer3 ntthal screen implementable; ε260 "nearest-neighbor (simplified)" was a mononucleotide sum; validated "CLEAN" with Tm 49.35 on a non-NN scale) → fixed | ❌ → fixed | FIXED | `19bdf7300` | ☐ |
| PROBE-VALID-001 | B07 | ⚠ (Kane et al. 2000 cited only for the > 75 % identity threshold — its > 15-nt contiguous-identity criterion missing, identity threshold applied as ≥; `ValidateProbe` self-structure still the heuristic fold-back fraction / fixed-loop stem although DesignProbes uses Primer3's ntthal screen; 1/N "specificity" presented as if sourced; `IsValid` lenient rule "≤ 1 hit ∧ fraction ≤ 0.4 ⇒ valid despite issues" unsourced; `maxMismatches = 3` justified by CRISPR guide data; "duplex-Tm off-target modelling not implemented" although implementable; status "Simplified", TestSpec "Not Started") → fixed | ❌ → fixed | FIXED | `2c6b7b9d3` | ☐ |
| PROBE-LNATM-001 | B07 | ⚠ (McTigue increments documented as applied to the SantaLucia & Hicks 2004 set and the 0.087 °C gap to MELTING `mct04` excused as a "base-set choice"; Owczarzy et al. 2011 single/consecutive/mismatch LNA model declared unobtainable/absent although its parameters and reference implementation are in MELTING 5; terminal/self-complementary handling and TestSpec "CLEAN" vs a non-reproduced worked example) → fixed | ❌ → fixed | FIXED | `576eee030` | ☐ |
| PROBE-EVALUE-001 | B07 | ⚠ (TestSpec/report cited blast_stat's 2/−3 non-affine megablast row λ 0.55 / K 0.21 as the "published ungapped" value; K declared to have "no citable closed form" although BLAST+ computes it (`BlastKarlinLHtoK`); no gapped/edge-effect model although the doc frames hits as BLAST-grade) → fixed | ❌ → fixed | FIXED | `f9795ab41` | ☐ |
| CRISPR-PAM-001 | B08 | PASS-with-notes | FAIL→FIXED | FIXED | — | ☐ |
| CRISPR-GUIDE-001 | B08 | — | — | not reviewed | — | ☐ |
| CRISPR-OFF-001 | B08 | — | — | not reviewed | — | ☐ |
| RESTR-FIND-001 | B08 | — | — | not reviewed | — | ☐ |
| RESTR-DIGEST-001 | B08 | — | — | not reviewed | — | ☐ |
| RESTR-FILTER-001 | B08 | — | — | not reviewed | — | ☐ |
| ALIGN-GLOBAL-001 | B09 | PASS-WITH-NOTES | PASS-WITH-NOTES | FIXED | `5e6c6f1ea` | ☐ |
| ALIGN-LOCAL-001 | B09 | PASS-WITH-NOTES | PASS-WITH-NOTES | FIXED | `4299f18a0` | ☐ |
| ALIGN-SEMI-001 | B09 | PASS-WITH-NOTES | PASS | FIXED | `1cd574217` | ☐ |
| ALIGN-MULTI-001 | B09 | — | — | not reviewed | — | ☐ |
| ALIGN-STATS-001 | B09 | — | — | not reviewed | — | ☐ |
| GENOMIC-COMMON-001 | B09 | — | — | not reviewed | — | ☐ |
| GENOMIC-MOTIFS-001 | B09 | — | — | not reviewed | — | ☐ |
| GENOMIC-ORF-001 | B09 | — | — | not reviewed | — | ☐ |
| GENOMIC-REPEAT-001 | B09 | — | — | not reviewed | — | ☐ |
| GENOMIC-SIMILARITY-001 | B09 | — | — | not reviewed | — | ☐ |
| GENOMIC-TANDEM-001 | B09 | — | — | not reviewed | — | ☐ |
| ASSEMBLY-OLC-001 | B10 | PASS-with-notes | FAIL-fixed | ✅ fixed | `1d719ff92` | ☐ |
| ASSEMBLY-DBG-001 | B10 | ⏳ | ⏳ | pending | — | ☐ |
| ASSEMBLY-STATS-001 | B10 | ⏳ | ⏳ | pending | — | ☐ |
| ASSEMBLY-MERGE-001 | B10 | ⏳ | ⏳ | pending | — | ☐ |
| ASSEMBLY-SCAFFOLD-001 | B10 | ⏳ | ⏳ | pending | — | ☐ |
| ASSEMBLY-COVER-001 | B10 | ⏳ | ⏳ | pending | — | ☐ |
| ASSEMBLY-CONSENSUS-001 | B10 | ⏳ | ⏳ | pending | — | ☐ |
| ASSEMBLY-TRIM-001 | B10 | ⏳ | ⏳ | pending | — | ☐ |
| ASSEMBLY-CORRECT-001 | B10 | ⏳ | ⏳ | pending | — | ☐ |
| ANNOT-ORF-001 | B11 | 🟡 PASS-WITH-NOTES | ❌→fixed | ✅ CLEAN | `8690148fd` | ☐ |
| ANNOT-GENE-001 | B11 | 🟡 PASS-WITH-NOTES | ❌→fixed | ✅ CLEAN (RBS score heuristic declared) | `b7712f246` | ☐ |
| ANNOT-PROM-001 | B11 | — | — | not reviewed | — | ☐ |
| ANNOT-GFF-001 | B11 | — | — | not reviewed | — | ☐ |
| ANNOT-CODING-001 | B11 | — | — | not reviewed | — | ☐ |
| ANNOT-REPEAT-001 | B11 | — | — | not reviewed | — | ☐ |
| ANNOT-CODONUSAGE-001 | B11 | — | — | not reviewed | — | ☐ |
| SPLICE-DONOR-001 | B11 | — | — | not reviewed | — | ☐ |
| SPLICE-ACCEPTOR-001 | B11 | — | — | not reviewed | — | ☐ |
| SPLICE-PREDICT-001 | B11 | — | — | not reviewed | — | ☐ |
| SPLICE-MAXENT3-001 | B11 | — | — | not reviewed | — | ☐ |
| SPLICE-MAXENT5-001 | B11 | — | — | not reviewed | — | ☐ |
| RNA-STRUCT-001 | B12 | — | — | not reviewed | — | ☐ |
| RNA-STEMLOOP-001 | B12 | — | — | not reviewed | — | ☐ |
| RNA-ENERGY-001 | B12 | — | — | not reviewed | — | ☐ |
| RNA-PAIR-001 | B12 | PASS-WITH-NOTES | FAIL (fixed) | FIXED | `440b2e505` | ☐ |
| RNA-HAIRPIN-001 | B12 | — | — | not reviewed | — | ☐ |
| RNA-MFE-001 | B12 | — | — | not reviewed | — | ☐ |
| RNA-PSEUDOKNOT-001 | B12 | — | — | not reviewed | — | ☐ |
| RNA-DOTBRACKET-001 | B12 | PASS-WITH-NOTES | FAIL (fixed) | FIXED | `0ed8d4b4c` | ☐ |
| RNA-INVERT-001 | B12 | — | — | not reviewed | — | ☐ |
| RNA-PARTITION-001 | B12 | — | — | not reviewed | — | ☐ |
| RNA-ACCESS-001 | B12 | — | — | not reviewed | — | ☐ |
| RNA-PKPREDICT-001 | B12 | — | — | not reviewed | — | ☐ |
| RNA-PKRECURSIVE-001 | B12 | — | — | not reviewed | — | ☐ |
| MIRNA-SEED-001 | B13 | PASS-WITH-NOTES | FAIL→FIXED | FIXED | `579f163a8` | ☐ |
| MIRNA-TARGET-001 | B13 | PASS-WITH-NOTES | FAIL→FIXED | FIXED | `937d2ae66` | ☐ |
| MIRNA-PRECURSOR-001 | B13 | — | — | not reviewed | — | ☐ |
| MIRNA-PAIR-001 | B13 | — | — | not reviewed | — | ☐ |
| MIRNA-CONTEXT-001 | B13 | — | — | not reviewed | — | ☐ |
| MIRNA-PCT-001 | B13 | — | — | not reviewed | — | ☐ |
| MIRNA-CLASSIFY-001 | B13 | — | — | not reviewed | — | ☐ |
| MIRNA-CLEAVAGE-001 | B13 | — | — | not reviewed | — | ☐ |
| PROTMOTIF-FIND-001 | B14 | — | — | not reviewed | `566c54697` | ☐ |
| PROTMOTIF-PROSITE-001 | B14 | — | — | not reviewed | `de4ddcdb1` | ☐ |
| PROTMOTIF-DOMAIN-001 | B14 | — | — | not reviewed | — | ☐ |
| PROTMOTIF-PATTERN-001 | B14 | — | — | not reviewed | — | ☐ |
| PROTMOTIF-SP-001 | B14 | — | — | not reviewed | — | ☐ |
| PROTMOTIF-TM-001 | B14 | — | — | not reviewed | — | ☐ |
| PROTMOTIF-CC-001 | B14 | — | — | not reviewed | — | ☐ |
| PROTMOTIF-LC-001 | B14 | — | — | not reviewed | — | ☐ |
| PROTMOTIF-COMMON-001 | B14 | — | — | not reviewed | — | ☐ |
| PROTMOTIF-HMM-001 | B14 | — | — | not reviewed | — | ☐ |
| DISORDER-PRED-001 | B15 | 🟡 PASS-WITH-NOTES | 🟡 PASS-WITH-NOTES | ✅ CLEAN | `06ce61d02` | ☐ |
| DISORDER-REGION-001 | B15 | 🟡 PASS-WITH-NOTES | 🟡 PASS-WITH-NOTES | ✅ CLEAN | `663a912f0` | ☐ |
| DISORDER-MORF-001 | B15 | — | — | not reviewed | — | ☐ |
| DISORDER-PROPENSITY-001 | B15 | — | — | not reviewed | — | ☐ |
| DISORDER-LC-001 | B15 | — | — | not reviewed | — | ☐ |
| EPIGEN-CPG-001 | B15 | — | — | not reviewed | — | ☐ |
| EPIGEN-METHYL-001 | B15 | — | — | not reviewed | — | ☐ |
| EPIGEN-DMR-001 | B15 | — | — | not reviewed | — | ☐ |
| EPIGEN-BISULF-001 | B15 | — | — | not reviewed | — | ☐ |
| EPIGEN-CHROM-001 | B15 | — | — | not reviewed | — | ☐ |
| EPIGEN-AGE-001 | B15 | — | — | not reviewed | — | ☐ |
| PHYLO-DIST-001 | B16 | PASS-with-notes (L=0 case wrongly specified as 0) | FIXED | FIXED | `b45d0b765` | ☐ |
| PHYLO-TREE-001 | B16 | — | — | not reviewed | — | ☐ |
| PHYLO-NEWICK-001 | B16 | — | — | not reviewed | — | ☐ |
| PHYLO-COMP-001 | B16 | — | — | not reviewed | — | ☐ |
| PHYLO-BOOT-001 | B16 | — | — | not reviewed | — | ☐ |
| PHYLO-STATS-001 | B16 | — | — | not reviewed | — | ☐ |
| POP-FREQ-001 | B16 | — | — | not reviewed | — | ☐ |
| POP-DIV-001 | B16 | — | — | not reviewed | — | ☐ |
| POP-HW-001 | B16 | — | — | not reviewed | — | ☐ |
| POP-FST-001 | B16 | — | — | not reviewed | — | ☐ |
| POP-LD-001 | B16 | — | — | not reviewed | — | ☐ |
| POP-SELECT-001 | B16 | — | — | not reviewed | — | ☐ |
| POP-ANCESTRY-001 | B16 | — | — | not reviewed | — | ☐ |
| POP-ROH-001 | B16 | — | — | not reviewed | — | ☐ |
| META-CLASS-001 | B17 | PASS | PASS (after F1) | FIXED | `19c25e642` | ☐ |
| META-PROF-001 | B17 | PASS-WITH-NOTES (Evidence doc invariants corrected) | PASS | CLEAN | `13632d173` | ☐ |
| META-ALPHA-001 | B17 | PASS (Chao 1984→1987 attribution for bias-corrected form corrected) | FAIL → FIXED | FIXED | `d589b2913` | ☐ |
| META-BETA-001 | B17 | — | — | not reviewed | — | ☐ |
| META-BIN-001 | B17 | — | — | not reviewed | — | ☐ |
| META-FUNC-001 | B17 | — | — | not reviewed | — | ☐ |
| META-RESIST-001 | B17 | — | — | not reviewed | — | ☐ |
| META-PATHWAY-001 | B17 | — | — | not reviewed | — | ☐ |
| META-TAXA-001 | B17 | — | — | not reviewed | — | ☐ |
| META-CHECKM-001 | B17 | — | — | not reviewed | — | ☐ |
| META-TETRA-001 | B17 | — | — | not reviewed | — | ☐ |
| PANGEN-CORE-001 | B18 | PASS-WITH-NOTES | FAIL→fixed | FIXED | `dbafec315` | ☐ |
| PANGEN-CLUSTER-001 | B18 | — | — | not reviewed | — | ☐ |
| PANGEN-HEAP-001 | B18 | — | — | not reviewed | — | ☐ |
| PANGEN-MARKER-001 | B18 | — | — | not reviewed | — | ☐ |
| COMPGEN-SYNTENY-001 | B18 | — | — | not reviewed | — | ☐ |
| COMPGEN-ORTHO-001 | B18 | — | — | not reviewed | — | ☐ |
| COMPGEN-REARR-001 | B18 | — | — | not reviewed | — | ☐ |
| COMPGEN-RBH-001 | B18 | — | — | not reviewed | — | ☐ |
| COMPGEN-COMPARE-001 | B18 | — | — | not reviewed | — | ☐ |
| COMPGEN-REVERSAL-001 | B18 | — | — | not reviewed | — | ☐ |
| COMPGEN-CLUSTER-001 | B18 | — | — | not reviewed | — | ☐ |
| COMPGEN-ANI-001 | B18 | — | — | not reviewed | — | ☐ |
| COMPGEN-DOTPLOT-001 | B18 | — | — | not reviewed | — | ☐ |
| CHROM-TELO-001 | B19 | PASS-WITH-NOTES | FAIL → fixed | FIXED | `e6279d018` | ☐ |
| CHROM-CENT-001 | B19 | PASS-WITH-NOTES | FAIL → fixed | FIXED (localisation heuristic kept as declared LIMITATION) | `d5ec859ea` | ☐ |
| CHROM-KARYO-001 | B19 | — | — | not reviewed | — | ☐ |
| CHROM-ANEU-001 | B19 | — | — | not reviewed | — | ☐ |
| CHROM-SYNT-001 | B19 | — | — | not reviewed | — | ☐ |
| CHROM-ALPHASAT-001 | B19 | — | — | not reviewed | — | ☐ |
| CHROM-HOR-001 | B19 | — | — | not reviewed | — | ☐ |
| TRANS-EXPR-001 | B19 | — | — | not reviewed | — | ☐ |
| TRANS-DIFF-001 | B19 | — | — | not reviewed | — | ☐ |
| TRANS-SPLICE-001 | B19 | — | — | not reviewed | — | ☐ |
| PARSE-FASTA-001 | B20 | PASS-with-fixes | PASS-with-fixes | 🔧 fixed | `c57912072` | ☐ |
| PARSE-FASTQ-001 | B20 | 🟡 PASS-with-fixes | ❌→🔧 fixed | ✅ CLEAN (🔧 7 fixes) | `99d8b890f` | ☐ |
| PARSE-BED-001 | B20 | — | — | not reviewed | — | ☐ |
| PARSE-VCF-001 | B20 | — | — | not reviewed | — | ☐ |
| PARSE-GFF-001 | B20 | — | — | not reviewed | — | ☐ |
| PARSE-GENBANK-001 | B20 | — | — | not reviewed | — | ☐ |
| PARSE-EMBL-001 | B20 | — | — | not reviewed | — | ☐ |
| QUALITY-PHRED-001 | B20 | — | — | not reviewed | — | ☐ |
| QUALITY-STATS-001 | B20 | — | — | not reviewed | — | ☐ |
| VARIANT-CALL-001 | B21 | PASS-WITH-NOTES (doc wrongly said indels not left-aligned) | FAIL→FIXED | FIXED+LIMITATION | `d31c1d7d5` | ☐ |
| VARIANT-SNP-001 | B21 | PASS-WITH-NOTES (doc claimed case-insensitive; unequal-length truncation treated as assumption) | FAIL→FIXED | FIXED | `a8a684882` | ☐ |
| VARIANT-INDEL-001 | B21 | — | — | not reviewed | — | ☐ |
| VARIANT-ANNOT-001 | B21 | — | — | not reviewed | — | ☐ |
| SV-DETECT-001 | B21 | — | — | not reviewed | — | ☐ |
| SV-BREAKPOINT-001 | B21 | — | — | not reviewed | — | ☐ |
| SV-CNV-001 | B21 | — | — | not reviewed | — | ☐ |
| ONCO-SOMATIC-001 | B22 | PASS-with-fixes | PASS | ☑ fixed | `cb637f86b` | ☐ |
| ONCO-VAF-001 | B22 | — | — | not reviewed | — | ☐ |
| ONCO-DRIVER-001 | B22 | — | — | not reviewed | — | ☐ |
| ONCO-ARTIFACT-001 | B22 | — | — | not reviewed | — | ☐ |
| ONCO-ANNOT-001 | B22 | — | — | not reviewed | — | ☐ |
| ONCO-TMB-001 | B22 | — | — | not reviewed | — | ☐ |
| ONCO-MSI-001 | B22 | — | — | not reviewed | — | ☐ |
| ONCO-HRD-001 | B22 | — | — | not reviewed | — | ☐ |
| ONCO-LOH-001 | B22 | — | — | not reviewed | — | ☐ |
| ONCO-SIG-001 | B23 | ⚠️ (channel-order doc gap) | ✅ (192/192 vs SPMG port) | FIXED (docs + dedup; no numeric change) | `42cf5cf05` | ☐ |
| ONCO-SIG-002 | B23 | — | — | not reviewed | `63d7fac4f` | ☐ |
| ONCO-SIG-003 | B23 | — | — | not reviewed | — | ☐ |
| ONCO-SIG-004 | B23 | — | — | not reviewed | — | ☐ |
| ONCO-FUSION-001 | B23 | — | — | not reviewed | — | ☐ |
| ONCO-FUSION-002 | B23 | — | — | not reviewed | — | ☐ |
| ONCO-FUSION-003 | B23 | — | — | not reviewed | — | ☐ |
| ONCO-CNA-001 | B24 | PASS-WITH-NOTES | FAIL → fixed | FIXED | `94fd19159` | ☐ |
| ONCO-CNA-002 | B24 | PASS-WITH-NOTES | PASS-WITH-NOTES (after F4) | LIMITED | `bc86794e8` | ☐ |
| ONCO-CNA-003 | B24 | PASS | FAIL → fixed | FIXED | `c35c48565` | ☐ |
| ONCO-PURITY-001 | B24 | PASS | FAIL (in B22 file) | LIMITED | `f4ebddf88` | ☐ |
| ONCO-PLOIDY-001 | B24 | FAIL → corrected | FAIL → fixed | FIXED | `b175277f8` | ☐ |
| ONCO-ASCAT-001 | B24 | FAIL → corrected | FAIL → fixed | FIXED | `9d5713566` | ☐ |
| ONCO-CLONAL-001 | B24 | PASS-WITH-NOTES | FAIL → fixed | FIXED | `5a94942d2` | ☐ |
| ONCO-CCF-001 | B24 | FAIL → corrected | FAIL → fixed | FIXED | `7e40c60d5` | ☐ |
| ONCO-PHYLO-001 | B24 | PASS-WITH-NOTES | FAIL → fixed | FIXED | `3bb86f375` | ☐ |
| ONCO-HETERO-001 | B24 | PASS-WITH-NOTES | FAIL → fixed | FIXED | `3e5c55ff7` | ☐ |
| ONCO-CTDNA-001 | B25 | 🟡 (3.3 pg/GE field convention vs 3.205 pg primary; kept, documented) | ❌→✅ | FIXED | `f889d5915` | ☐ |
| ONCO-MRD-001 | B25 | — | — | not reviewed | — | ☐ |
| ONCO-CHIP-001 | B25 | — | — | not reviewed | — | ☐ |
| ONCO-HLA-001 | B25 | — | — | not reviewed | — | ☐ |
| ONCO-ACTION-001 | B25 | — | — | not reviewed | — | ☐ |
| ONCO-SV-001 | B25 | — | — | not reviewed | — | ☐ |
| ONCO-EXPR-001 | B25 | — | — | not reviewed | — | ☐ |
| ONCO-NEO-001 | B26 | PASS-WITH-NOTES (default length 8–14 vs docs 8–11 aligned; stop-gain unspecified) | FAIL → FIXED | FIXED | `6bbad85da` `93fc76621` | ☐ |
| ONCO-MHC-001 | B26 | — | — | not reviewed | — | ☐ |
| ONCO-IMMUNE-001 | B26 | — | — | not reviewed | — | ☐ |
| MHC-NN-001 | B26 | — | — | not reviewed | — | ☐ |
| MHC-MATRIX-001 | B26 | — | — | not reviewed | — | ☐ |
| IMMUNE-NUSVR-001 | B26 | — | — | not reviewed | — | ☐ |

## Batches with units not yet recorded as reviewed

- B03: 4 unit(s) without a verdict row in `B03.md`
- B04: 3 unit(s) without a verdict row in `B04.md`
- B08: 5 unit(s) without a verdict row in `B08.md`
- B09: 8 unit(s) without a verdict row in `B09.md`
- B11: 10 unit(s) without a verdict row in `B11.md`
- B12: 11 unit(s) without a verdict row in `B12.md`
- B13: 6 unit(s) without a verdict row in `B13.md`
- B14: 10 unit(s) without a verdict row in `B14.md`
- B15: 9 unit(s) without a verdict row in `B15.md`
- B16: 13 unit(s) without a verdict row in `B16.md`
- B17: 8 unit(s) without a verdict row in `B17.md`
- B18: 12 unit(s) without a verdict row in `B18.md`
- B19: 8 unit(s) without a verdict row in `B19.md`
- B20: 7 unit(s) without a verdict row in `B20.md`
- B21: 5 unit(s) without a verdict row in `B21.md`
- B22: 8 unit(s) without a verdict row in `B22.md`
- B23: 6 unit(s) without a verdict row in `B23.md`
- B25: 6 unit(s) without a verdict row in `B25.md`
- B26: 5 unit(s) without a verdict row in `B26.md`
