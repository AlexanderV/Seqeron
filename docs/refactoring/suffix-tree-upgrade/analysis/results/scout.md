# Completeness scouting (from the K11/K12 agent) - NOT scored
Date 2026-11-02. [O] = I opened the URL this session; [S] = URL seen only in a search result (not opened; verify before use).

## 1. Competitors a sceptical reviewer would expect (missing from K1..K12)
1. Movi 2 (move-structure pangenome index, 2025): default index 56% smaller than Movi 1, PML queries ~35x faster than SPUMONI; HPRC R2 (~460 haplotypes). Biggest r-index rival. [O] https://pmc.ncbi.nlm.nih.gov/articles/PMC12632922/ (Movi 1: iScience 2024 [S] https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11696632/)
2. Move-r (SEA 2024) and Move-rb (arXiv 2609.30089, Sep 2026, bidirectional r-index with approximate matching, said to beat br-index, b-move, columba): the move-structure answer to C3/C12. [O] https://arxiv.org/abs/2609.30089 (ID format unusual, fetched via WebFetch only); Move-r [S] https://drops.dagstuhl.de/storage/00lipics/lipics-vol301-sea2024/LIPIcs.SEA.2024.1/LIPIcs.SEA.2024.1.pdf
3. Columba + b-move (bidirectional FM / bidirectional move, optimal search schemes, up to 13 errors, 1M 151bp reads on hg38 in 24 s at 0 errors; RLC variant 5.6x less RAM on 6,115-strain pangenome; AGPL). SeqAn3 (K2) is its only peer in C3. [O] https://pmc.ncbi.nlm.nih.gov/articles/PMC12724072/
4. Suffixient arrays (Cenzato et al., 2024/25): tiny SA subset that locates one occurrence and all MEMs by binary search, I/O-model efficient, 3.5 ns/char, claimed 10-100x faster than r-index. Possible successor to the r-index for MEMs. [O] https://arxiv.org/abs/2407.18753
5. grlBWT (and block_RLBWT): fastest BWT builder in ropebwt3's benchmark (human100 in 8.3 h but ~2 TB scratch disk). Must sit in the C10/C11/C14 comparison. [O, link in README] https://github.com/ddiazdom/grlBWT
6. sr-index (subsampled r-index, CPM 2021) and r-index-f (SEA 2022): ropebwt3's author calls sr-index "probably the better solution" for 14.6T-symbol data; r-index-f in the same paper's untested list. [O] https://arxiv.org/html/2409.00613v2
7. copMEM2 (Bioinformatics 2023): multithreaded sampling MEM finder; human-vs-mouse MEMs >=50 in 59 s, 10.4 GB, 12 threads; compared with E-MEM, bfMEM, copMEM. Modern replacement for the sparseMEM/essaMEM/E-MEM notes under K4; paper does not mention MUMmer4. [O] https://pmc.ncbi.nlm.nih.gov/articles/PMC10209524/
8. gsufsort (AMB 2020): linear-time SA/LCP/BWT/document array for string collections (gSACA-K); direct rival to K5 for generalized SA. [S] https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7507297/
9. GBWT + GCSA2 (vg, Bioinformatics 2020): haplotype-aware pangenome index; GBWT of 1000 Genomes and 108,070 TOPMed haplotypes. Reference for "pangenome-scale" C18. [S] https://pmc.ncbi.nlm.nih.gov/articles/PMC7223266
10. k-mer/BWT-based indexes: SBWT + Themisto (colored k-mer index), Fulgor, Metagraph. ropebwt3's paper benchmarks the last two against FM-indexes (metagraph human100 build 16.9 h, 251 GB), so they do compete on exact lookup. [S] https://github.com/algbio/Themisto ; benchmark [O] https://pmc.ncbi.nlm.nih.gov/articles/PMC11646566/
11. GPU: only research code found, nothing mature at genome scale. Mainz multi-GPU SA thesis (2023; the 3.44 s human-genome figure is unverified), BarraCUDA FM-index aligner, haystackfm (Rust/WebGPU FM-index). Treat C15's GPU clause as a gap for every tool. [S] https://openscience.ub.uni-mainz.de/handle/20.500.12030/9652 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC3278344 ; https://docs.rs/crate/haystackfm/0.2.0
12. Compressed suffix trees 2021-26: no new library surfaced; sdsl-lite CST (K1) and the r-index suffix-tree theory (JACM 2020, no released tree API) remain the reference. [O] https://arxiv.org/abs/1809.02792

## 2. Capabilities that distinguish tools but are not in C1..C24
1. Bidirectional extension (extend left AND right, FMD / b-move): prerequisite for SMEMs and search schemes; C3/C8 only hint at it. [O] https://raw.githubusercontent.com/lh3/ropebwt3/v3.10/fm-index.h (rb3_fmd_extend)
2. Search schemes (optimal partitioning for k errors) as a distinct mode from backtracking: Columba/SeqAn3. [O] https://pmc.ncbi.nlm.nih.gov/articles/PMC12724072/
3. Streaming/online matching statistics with bounded latency and RAM (patterns never buffered): PHONI vs MONI differ here, C5 does not capture it. [O] https://arxiv.org/pdf/2011.05610
4. Locate with sampled-SA tradeoff (sample rate s vs time/space; ropebwt3 expected O(s/m) for m similar genomes; r-index sampling at run boundaries). [O] https://raw.githubusercontent.com/lh3/ropebwt3/v3.10/NEWS.md (release 3.8)
5. Top-k / document listing / document counting (C6 only partly covers it; SPUMONI has a document array, Gagie et al. 2017 give repetitive-collection indexes). [S] https://users.dcc.uchile.cl/~gnavarro/abstracts/ir17.html
6. Coordinate translation: global offset to (contig, offset, strand) plus names/lengths sidecar; essential for genomics, absent from the C list. [O] https://raw.githubusercontent.com/lh3/ropebwt3/v3.10/fm-index.h (rb3_pos_t, rb3_ssa_t)
7. Text-side ambiguity/masking policy (N, IUPAC, lowercase in the INDEXED text); C2 covers only the pattern. MONI strips non-ACGTN, ropebwt3 stores ACGTN. [O] https://pmc.ncbi.nlm.nih.gov/articles/PMC8892979
8. Index merge / union of two indexes (not the same as C17 append): ropebwt3 'merge', BWT merge in the paper. [O] https://raw.githubusercontent.com/lh3/ropebwt3/v3.10/main.c
9. Canonical/deterministic serialization: ropebwt3 README says the same BWT need not give the same FMR file, so a content hash (C20 integrity) is not reproducible. [O] https://raw.githubusercontent.com/lh3/ropebwt3/v3.10/README.md
10. Random access to SA / ISA / LCP / Phi (LCE queries) on the compressed index: JACM r-index theory; relevant to C8/C9. [O] https://arxiv.org/abs/1809.02792
11. Affine-gap alignment directly on the index (BWA-SW over a DAWG) versus k-error enumeration; C3 conflates the two. [O] https://raw.githubusercontent.com/lh3/ropebwt3/v3.10/README.md
12. Classification workloads (pseudo-matching lengths, host depletion, adaptive sampling) as a first-class query; SPUMONI. [O] https://pmc.ncbi.nlm.nih.gov/articles/PMC8237286/
13. Cache/latency behaviour per query character (move structure: O(1) steps, up to 30x faster backward search) as a sub-measure of C12. [S] https://www.biorxiv.org/content/10.1101/2023.11.04.565615.full.pdf

## 3. Caveats found while scoring K11/K12 (relevant to the whole matrix)
- Independent head-to-head data exist: ropebwt3 paper Tables 3-4 (build + query for ropebwt3, grlBWT, pfp-thresholds, MONI, Movi, metagraph, fulgor). Use it to calibrate C10-C12 for all BWT-class tools. [O] https://pmc.ncbi.nlm.nih.gov/articles/PMC11646566/
- Licence conflict: MONI repo LICENSE is GPL-3.0 but the bioconda recipe says MIT. [O] https://bioconda.github.io/recipes/moni/README.html
- Big-BWT upstream is GitLab (gitlab.com/manzai/Big-BWT, README file is readme.md); MONI/SPUMONI point to a GitHub copy (alshai/Big-BWT). Big-BWT switched to libsais in Aug 2026 commits. [O] https://gitlab.com/manzai/Big-BWT
