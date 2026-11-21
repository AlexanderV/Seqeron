# minimap2 2.26: base-level alignment, paired reads and the program options

**Units:** ALIGN-MULTI-001 (audit A3-9c2, A5-9b, A5-9c; B09 F187 - F190).
**Code:** `Minimap2Chaining.Align.cs` (`align.c`), `Minimap2Chaining.Pe.cs` (`pe.c`, `mm_seg_gen`), `Minimap2Chaining.Pipeline.cs` (`mm_map_frag`), `Minimap2Chaining.Output.cs` (`format.c`), `Minimap2Chaining.Index.cs` (`--junc-bed`, `--alt`), `Minimap2Program.cs` (`main.c`, `options.c`, `bseq.c` / `kseq.h`,
the pipeline of `map.c`), `Ksw2Aligner.Splice.cs` (`ksw_exts2_sse`), `Ksw2Aligner.Local.cs` (`ksw_ll_i16`). MCP: `minimap2_map`, `ksw2_extend` (variants `exts2`, `ll`).
**Reference:** minimap2 v2.26-r1175 (Li 2018, Bioinformatics 34:3094; Li 2021, Bioinformatics 37:4572); sources read at raw.githubusercontent.com/lh3/minimap2/v2.26: `align.c`, `hit.c`, `pe.c`, `map.c`, `format.c`, `main.c`, `options.c`, `bseq.c`, `bseq.h`, `index.c`, `lchain.c`, `seed.c`, `splitidx.c`, `ksw2.h`, `ksw2_ll_sse.c`, `ksw2_exts2_sse.c`, `ksw2_extz2_sse.c`,
`ksw2_extd2_sse.c`, `ketopt.h`, `mmpriv.h`, `minimap.h`, `kseq.h` (the in-memory reader is a port of `kseq_read` / `ks_getuntil2`). The installed program (Debian 2.26-r1175, one SSE4.1 code path, no CPU dispatch: no `cpuid`, no AVX, `pblendvb` / `pmaxsb` present) is the oracle.

## 1. What is ported

`Minimap2Program.Run(args, files)` is `minimap2 [options] ref.fa query.fa [query2.fa]` on in-memory files: it returns what the program writes to the standard output (PAF or SAM), to the standard error and to the `-o` file, and its exit status. Underneath, `Minimap2Chaining.Map` (one read) and
`Minimap2Chaining.MapFragment` (one or more reads of a fragment) are `mm_map_frag`; `Minimap2Chaining.FormatFragment` is the output step of the pipeline (`mm_write_paf3` / `mm_write_sam3`).

### 1.1 Base-level alignment (`-c`, `-a`, `--cs`, `--MD`; `align.c`)

`MinimapMapOptions.Cigar` switches `mm_align_skeleton` on. Per chain (`mm_align1`):
1. the anchors are squeezed (`mm_squeeze_a`), the first / last anchors that look like noise are dropped (`mm_fix_bad_ends`; in splice mode `mm_fix_bad_ends_splice` with a `ksw_ll` extension score of the end seed; `--no-end-flt` switches it off, short reads take the longest exact stretch `mm_max_stretch`), seeds in long indel regions are ignored (`mm_filter_bad_seeds`,
   `mm_filter_bad_seeds_alt`: `MM_SEED_IGNORE`, `MM_SEED_LONG_JOIN`);
2. the extension regions on both sides are computed from the neighbouring seeds of other chains and the `-g` / `-A` / `-O` / `-E` limits (`rs0 qs0 re0 qe0`);
3. the left extension (reversed sequences, `KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT | KSW_EZ_REV_CIGAR`), the gap filling between the retained anchors (first pass `KSW_EZ_APPROX_MAX`, then, when `mm_test_zdrop` sees a score drop, a second pass with the exact Z-drop; code 2 = the drop may hide an inversion, tested by
   `ksw_ll_i16`) and the right extension call `mm_align_pair`: `ksw_extz2_sse` (equal gap costs), `ksw_extd2_sse` (dual affine `-O q,q2 -E e,e2`) or, with `--splice`, `ksw_exts2_sse` (donor / acceptor scores, `-C`, `--junc-bonus`, `--junc-bed`, `-u`, `-J`); the short-read mode aligns the max stretch without gaps;
4. a Z-drop that does not extend beyond `min_cnt` anchors truncates the chain, otherwise `mm_split_reg` splits it (`zd:i:`) and `mm_align1_inv` tries to rescue an inversion between the two parts (`tp:A:I` / `i`);
5. `mm_update_extra` left-aligns the indels (`mm_fix_cigar`: shifts, `5I6D7I` -> `5I` + `6D`, leading I / D removed), counts matches, block length, ambiguous bases (`nn`) and the best-scoring segment (`ms`), optionally rewrites `M` as `=` / `X` (`--eqx`).

After the loop `mm_filter_regs` (`--max-clip-ratio`, `-s`), `mm_update_dp_max` (queries of at least 500 bases with similar best alignments), `mm_hit_sort` by `ms`, `mm_set_parent` with the DP maxima (`dp_max2`, `n_sub`), `mm_select_sub` and `mm_set_sam_pri`; `mm_set_mapq` uses `ms` / `ms2` and the identity. Output tags: `NM ms AS nn ts tp cm s1 s2 de zd rl cg cs MD`; SAM:
flags 0x1 0x2 0x4 0x8 0x10 0x20 0x40 0x80 0x100 0x800, clipping `S` / `H` (`-Y`, `--secondary-seq`), `SA`, `RG`, `FI`, mate fields and TLEN, `CG:B:I` for more than 65 535 operations (`-L`), header `@HD @SQ @RG @PG` with the permuted command line.

### 1.2 Kernels

`Ksw2Aligner.ExtendS2` (`ksw_exts2_sse`, new), `LocalScore` (`ksw_ll_qinit` + `ksw_ll_i16`, new; the zero-scored padding lanes of the striped layout can hold the maximum, so the end coordinate can lie beyond the query, which `mm_align1_inv` then reads as bytes before the region), and the F125 ports `ExtendZ2` / `ExtendD2`
are run lane by lane on the byte layout of the C code. `ExtendD2` / `ExtendD` accept a second extension penalty 0 (`-E 2,0`; the earlier check refused it). The Alignment project is built with `Optimize = true` (the lane emulation is about four times slower unoptimised).

### 1.3 Paired reads (`map.c`, `hit.c`, `pe.c`)

Two query files (or an interleaved file in fragment mode: consecutive reads with the same name after removing `/N`) are one fragment of `n_segs` reads. `collect_minimizers` sketches every read with its segment id and offsets the positions; `comput_sc` forbids a step inside a read that skips more than the
reference gap (`n_seg > 1`), pays the cheaper penalty (or a bonus when the reads overlap) between reads; the re-chaining with `max_occ` runs when the best chain does not cover all reads; `mm_select_sub_multi` keeps secondary chains by the distance on the reference (`pri1 = 0.2`, `pri2 = 0.7`); `mm_seg_gen` splits the fragment chains into
the chains of the reads; every read gets `mm_set_parent`, the alignment and `mm_set_mapq`; with `-c` (also `-a`) `mm_pair` pairs the hits (`pe_bonus` 33, `dp_thres`, orientation key `rid << 32 | rs << 1 | (s ^ rev)`), lifts the paired hits to primary, raises the mapping quality (`6.02 * dd / a - 4.343 * ln(n_sub)`), sets `proper_frag`
(SAM 0x2) and `mm_set_pe_thru`. `pe_ori` (`-x sr`: 1 = forward-reverse; the default 0 = forward-forward) reverse-complements the first / second read before mapping and flips the hits back. 2.26 has no `--pe-ori` option.

### 1.4 The other options (`main.c`)

Parsed as `ketopt` does (permuted arguments, `--opt=value` / `--opt value`, unique prefixes of long options, `--cs` / `--rmq` with `=` only), `-x` applied first and cumulatively over the current options (`-x splice -x asm5` keeps the splice flags and fails the `--rmq` check), `mm_check_opt` messages and exit
status, `-v` verbosity of warnings and errors, `-V`, `-h`, a missing query file or option argument, `-R` read groups, `-o FILE`, `--junc-bed`, `--alt` / `--alt-drop`, `--qstrand`, `--no-pairing`, `--secondary=`, `--paf-no-hit`, `--sam-hit-only`, `-y`, `-Q`, `--max-qlen`, `--seed`, `--no-hash-name`, `--heap-sort=`, `--dual=`, `--frag=`, `--end-seed-pen`, `--min-dp-len`, `--cap-sw-mem`,
`--min-occ-floor`, `--score-N`, `--max-chain-iter`, `--chain-gap-scale`, `--chain-skip-scale`, `-u`, `-J`, `-C`, `--junc-bonus`, `--splice-flank=`, `--end-bonus`, `-z`, `-O`, `-E`, `-s` (default 80: the program sets `min_dp_max` once from the defaults 40 x 2, it does not follow `-m` / `-A`). The reader is a port of `kseq_read`
(CRLF blank lines become a one-character sequence, a `@` or `+` at the start of a sequence line ends it, a malformed FASTQ record stops the file with a warning).

**`-I` and `--split-prefix`.** The reference is read in mini-batches of `min(50 Mbases, -I)` sequences and a part ends when the bases read exceed `-I`; every part is indexed and mapped on its own (`mid_occ` is computed from the first part and reused), the output of a part follows the previous one and a SAM header has no `@SQ` lines (warning at `-v 2`).
`--split-prefix` keeps the hits of all parts and merges them (`merge_hits`: the target ids shifted, `mm_update_dp_max`, sort, `mm_set_parent`, `mm_select_sub`, `mm_set_mapq`, `mm_pair` with the first part's fragment gap); the SAM `@SQ` lines are then written after `@PG`. The merge reads the query files again with `n_fp` = the number of query
files, so two files without the fragment mode pair the wrong reads with the stored hits: the port reproduces that.

**No effect on the standard output** (measured: oracle cases carry them): `-t`, `-K` (a mini-batch of the query reader is never split inside a fragment below 1 Mbase), `-d`, `--bucket-bits`, `--no-kalloc`, `--cap-kalloc`, `--print-qname`, `--print-seeds`, `--print-chains`, `--print-aln-seq` (the debug log goes to the standard error), `-2`, `--idx-no-seq` (an error with `-c`).

## 2. Verification

| Oracle | What | Result |
|---|---|---|
| `minimap2_align_oracle.jsonl.gz` | whole stdout of 2 700 command lines: 6 presets x 150, random alignment options, SAM, short reads, splice, repeat families, Z-drop / inversion rescue | 2 700 / 2 700 equal |
| `minimap2_pair_oracle.jsonl.gz` | 2 100 command lines of the paired mode (two files, interleaved, FR / FF, PAF with / without `-c`, SAM, options) | 2 100 / 2 100 equal |
| `minimap2_options_oracle.jsonl.gz` | 2 403 command lines: index parts / split-prefix, chaining and seeding options, output flags, `--alt`, `--qstrand`, all-vs-all, unusual input, a 610 kb query with `CG:B:I`, option errors (stdout + stderr at `-v 2` + exit status + `-o` files) | 2 403 / 2 403 equal |
| `ksw2_c_oracle.tsv` | the C `ksw_ll_i16`, `ksw_exts2_sse`, `ksw_extz2_sse`, `ksw_extd2_sse` (v2.26, `gcc -msse4.1`) on 4 x 500 random inputs, every field of `ksw_extz_t` and the CIGAR | 2 000 / 2 000 equal |

`valgrind -q --track-origins=yes` on two cases of every group of the three program oracles: no invalid access, no use of uninitialised values. The program has one SSE4.1 code path (no run-time dispatch), so its output does not depend on the CPU.

## 3. Differences and limits

- The `[M::...]` progress log of the program is not reproduced; `-d` writes no index file; a k above 28, a window of 256 or more and similar inputs that make the C program abort raise `ArgumentOutOfRangeException`.
- The warnings of the reader that the program prints per mini-batch are printed once per part.
- Memory: the traceback of the ksw2 SIMD kernels is limited by `Ksw2Aligner.MaxTracebackBytes` (1 GiB; the C code has no limit).
