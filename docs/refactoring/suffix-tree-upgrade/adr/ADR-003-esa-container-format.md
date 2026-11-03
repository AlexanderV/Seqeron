# ADR-003: ESA persistent container (a separate format kind next to v6)

Status: **proposed** (work package S1-0, 2026-11-03). Contract: [ADR-001](ADR-001-index-contract.md). Construction and LCP encodings: [ADR-002](ADR-002-sa-construction.md). Paths are repo-relative, `file:line` at commit `b0922125`. **UNVERIFIED** = not checked against code or a measurement.

## 1. Decisions and constraints

* The ESA container is a **separate kind** with its own magic; the v6 suffix-tree format and its loader stay untouched and readable (owner decision; PROCESS.md rule 2).
* Managed only: files are read through `MemoryMappedFile` (BCL) and plain `FileStream`/`RandomAccess`; **no P/Invoke and no `madvise`/prefetch hints** (the v6 provider uses `LibraryImport`/`posix_madvise` for advisory prefetch, `MappedFileStorageProvider.Prefetch.cs:23-76`; the ESA reader does not copy that).
* Positions are 64-bit in the contract ([ADR-001](ADR-001-index-contract.md) section 8), so the container carries a position width of 4 or 8 bytes instead of fixing `int` (v6 fixes `int32` text length and node count, `PersistentConstants.cs:57,63`).
* Immutable, read-only, zero-copy mapped access, thread-safe concurrent readers. Same file content for the same input and options (deterministic writer: no timestamps, no random ids) so that S5 can assert byte-identical serial and parallel builds.

## 2. Identification

| Item | Value |
|---|---|
| Kind name | **SESA** ("Seqeron enhanced suffix array"), format version 1.0 |
| Magic (8 bytes at offset 0) | `89 53 45 53 41 0D 0A 1A` = `0x89 'S' 'E' 'S' 'A' CR LF SUB` (PNG-style: high bit set, line-ending bytes catch text-mode corruption). As a little-endian `int64`: `0x1A0A0D4153455389` |
| Recommended file extension | `.sesa` (v6 documentation uses `.suffix`, `SuffixTree.Persistent/README.md`; no extension is enforced anywhere) |
| v6 magic for comparison | `PersistentConstants.MAGIC_NUMBER = 0x5452454558494646` (`PersistentConstants.cs:33`), read with `ReadInt64` and compared at `PersistentStorageHeaderValidation.cs:92,117` |
| Documentation defect found | the v6 constant is commented and documented as ASCII `"SUFFIXTR"` (`PersistentConstants.cs:33`, README table) but decodes to `"TREEXIFF"` (big-endian) / bytes `46 46 49 58 45 45 52 54` ("FFIXEERT", little-endian write; the write path was not inspected, UNVERIFIED). The two magics differ in every byte, so dispatch is unaffected. The text is an S0 documentation item, not touched here |

## 3. File layout

```
[ header page: 4096 bytes = 128-byte header + section table (up to 62 entries) + zero padding ]
[ sections, each starting at a 64-byte aligned offset; sections >= 1 MiB at a 4096-byte aligned offset ]
```
All integers little-endian. Elements are naturally aligned. No pointers inside the file, only indices, so the file is position independent. Section order written by v1: TEXT, ALPHABET, SA, LCP, LCP_OVERFLOW, CHILD, BUCKET, NEXTSYM, ISA, RMQ, BUILD_INFO (gaps are zero padding).

### Header (128 bytes)

| Offset | Size | Field | Meaning |
|---|---|---|---|
| 0 | 8 | `MAGIC` | section 2 (written **last**, section 6) |
| 8 | 2 | `VERSION_MAJOR` | 1 |
| 10 | 2 | `VERSION_MINOR` | 0 |
| 12 | 4 | `ENDIAN_PROBE` | `0x01020304` stored little-endian (bytes `04 03 02 01`) |
| 16 | 8 | `TEXT_LEN` | `n` (uint64) |
| 24 | 8 | `NODE_COUNT` | logical `NodeCount` of the terminated tree (ADR-001 section 4); equals the Ukkonen tree's value |
| 32 | 8 | `FILE_SIZE` | expected file length, validated on open (as v6 `SIZE`, `PersistentConstants.cs:73`) |
| 40 | 8 | `SECTION_TABLE_OFFSET` | 128 |
| 48 | 4 | `SECTION_COUNT` | <= 62 |
| 52 | 4 | `SECTION_ENTRY_SIZE` | 64 |
| 56 | 8 | `LRS_DEPTH` | longest repeated substring length (0 if none) |
| 64 | 8 | `LRS_RANK` | first rank `r` with `LCP[r] == LRS_DEPTH`, `-1` if none (O(1) `LongestRepeatedSubstring`, like v6 `DEEPEST_NODE`/`LRS_DEPTH`, `PersistentConstants.cs:87,93`) |
| 72 | 4 | `POS_WIDTH` | 4 or 8: byte width of SA, ISA, CHILD and BUCKET entries; 4 requires `n + 1 <= 2^32 - 1` |
| 76 | 4 | `SIGMA` | number of distinct code units in the text |
| 80 | 32 | `TEXT_SHA256` | SHA-256 of the text as UTF-16LE code units (same payload as `SuffixTreeSerializer`, `README.md` "Serialization Format v2"), zeros if not requested |
| 112 | 8 | `BUILD_FLAGS` | bit 0: SA is sentinel-augmented (`SA[0] = n`); always 1 in v1 |
| 120 | 4 | reserved | 0 |
| 124 | 4 | `HEADER_CRC32C` | CRC-32C over bytes `[0,124)` and the whole section table (CRC field treated as zero) |

### Section entry (64 bytes, at `128 + 64*i`)

| Offset | Size | Field | Meaning |
|---|---|---|---|
| 0 | 4 | `KIND` | section 4 |
| 4 | 4 | `FLAGS` | bit 0 `REQUIRED` (a reader that does not know the kind must refuse the file), bit 1 `HAS_CRC` |
| 8 | 8 | `OFFSET` | absolute file offset, aligned as above |
| 16 | 8 | `LENGTH` | payload bytes (padding excluded) |
| 24 | 8 | `COUNT` | number of elements |
| 32 | 4 | `ELEM_WIDTH` | bytes per element (0 for non-array payloads) |
| 36 | 4 | `CODEC` | kind specific |
| 40 | 8 | `PARAM0` | kind specific |
| 48 | 8 | `PARAM1` | kind specific |
| 56 | 4 | `CRC32C` | CRC-32C of the payload when `HAS_CRC`, else 0 |
| 60 | 4 | reserved | 0 |

## 4. Section kinds (v1)

| Kind | Name | Required | Content |
|---|---|---|---|
| 1 | `TEXT` | yes | `COUNT = n`. Codec 1: one byte per code unit, identity (code units < 256: ASCII, Latin-1, DNA, protein). Codec 2: dense order-preserving rank bytes, needs `ALPHABET` (any text with <= 256 distinct code units). Codec 3: UTF-16LE `ushort`. Codec 4: reserved (2-bit packed DNA, S2) |
| 2 | `ALPHABET` | if codec 2 or `BUCKET` | sorted distinct code units, `uint16` each, `COUNT = SIGMA` |
| 3 | `SA` | yes | `COUNT = n + 1`, width `POS_WIDTH`, `SA[0] = n`, codec 1 plain |
| 4 | `LCP` | yes | `COUNT = n + 1`. Codec 1: `U8` (escape value 255), codec 2: `U16` (escape 65535), codec 3: `U32`, codec 4: `U64`. `PARAM0` = number of overflow entries. `LCP[0] = LCP[1] = 0` |
| 5 | `LCP_OVERFLOW` | iff `PARAM0 > 0` | sorted pairs `(rank, value)`, each field `POS_WIDTH` bytes, for the entries equal to the escape value |
| 6 | `CHILD` | yes | child table, width `POS_WIDTH`, `COUNT = n + 1`. Codec 1: Abouelhoda-style single-array packing of `up/down/nextlIndex`; the exact bit layout is fixed in WP 1b.2 (**UNVERIFIED**; prototype P3 compares it with a 3-array reference). Further layouts get new codec ids without a format bump |
| 7 | `BUCKET` | no | first-`k`-symbol table: `PARAM0 = k`, `PARAM1 = number of dense symbol codes`, `COUNT = sigma^k + 1` first ranks, width `POS_WIDTH`. Default `k` = largest with `sigma^k * POS_WIDTH <= n/2` bytes (DNA, 10 Mbp: `k = 10`, 4 MB) |
| 8 | `ISA` | no (navigation pack) | `COUNT = n + 1`, width `POS_WIDTH`; `ISA[SA[i]] = i` for text positions, `ISA[n] = 0` |
| 9 | `RMQ` | no (navigation pack) | range-minimum / previous- and next-smaller helper over `LCP`: `PARAM0` = block size, codec 1 = block minima + sparse table, in-block scan (about 0.3-0.5 B/char, UNVERIFIED) |
| 10 | `NEXTSYM` | no (fast child lookup) | `NEXTSYM[i] = text[SA[i] + LCP[i]]` (dense code), 1 byte (2 for codec 3), `COUNT = n + 1` |
| 11 | `BUILD_INFO` | no | UTF-8 `key=value` lines: tool version, git commit of the writer, SA-IS variant, options. No timestamps |
| 12-4095 | reserved | - | S3 document map, S4 reversed-text link, S7 BWT / wavelet / samples, S2 packed text |
| >= 2^28 | private | - | experiments; always written with `REQUIRED` clear |

Notes:
* The navigation pack (`ISA`, `RMQ`) is optional because the suffix-link decision of [ADR-001](ADR-001-index-contract.md) section 5 makes it needed only by the matching-statistics family; when absent the engine builds it lazily in memory (never writes to a read-only mapping).
* Everything dense is a plain array so an engine can wrap the mapping as `ReadOnlySpan<T>` for `count <= int.MaxValue` or as a pointer-plus-`long` view (`SaView<T>`) beyond.

## 5. Size budget (index bytes per input character, text included)

| Section | Random DNA 10 Mbp (this gate) | Worst case |
|---|---|---|
| `TEXT` | 1.0 | 2.0 (UTF-16 codec 3) |
| `SA` (4-byte) | 4.0 | 4.0 (8.0 with `POS_WIDTH` 8) |
| `LCP` | 1.0 (U8, escapes ~0) | 4.0 (`a^n`, U32) |
| `CHILD` | 4.0 | 4.0 |
| `BUCKET` | <= 0.4 | <= 0.5 by construction |
| **Core** | **about 10.4** | about 14.5 |
| `ISA` | 4.0 | 4.0 |
| `RMQ` | about 0.4 | about 0.4 |
| **Full (Core + navigation pack)** | **about 14.8** | about 18.9 |
| + `NEXTSYM` | +1.0 = about 15.8 (still <= 16) | +1.0 |

The 16 B/char gate is defined on random DNA (PLAN.md:83: expected about 10 for the core). The table shows Core about 10.4 and Full about 14.8-15.8 on the gate input (estimates from array sizes, UNVERIFIED until the harness measures real files); adversarial inputs exceed 16 and are reported per class, not hidden. Compared with v6: 56.2 B/char for the same input (`BENCHMARK.md:23`).

## 6. Integrity, validation and trust

| Level | Cost | Checks |
|---|---|---|
| `Header` (default `Open`) | O(sections), cold open of the order of milliseconds like v6 (4.84 ms, `BENCHMARK.md:177`) | magic, version major, endian probe, header CRC, `FILE_SIZE == actual length`, every section inside the file, aligned, non-overlapping, `COUNT * ELEM_WIDTH == LENGTH`, required sections present and consistent with `TEXT_LEN`, `POS_WIDTH`, `SIGMA`, codec/alphabet agreement |
| `Sections` | O(file size) (CRC-32C by hardware intrinsics where available, slice-by-8 managed fallback; both give the same value) | + every `HAS_CRC` payload |
| `Structure` | O(n) | + `SA` is a permutation of `0..n` with `SA[0] = n` (bitmap); `LCP` recomputed with the Phi pass of ADR-002 section 4 and compared; adjacent-rank order check using `LCP` (`text[SA[i-1]+LCP[i]] < text[SA[i]+LCP[i]]` or the shorter suffix first); child table round trip; `ISA` inverse of `SA`; `RMQ` spot checks |

Access modes: `Checked` (default) clamps/validates every index read from `SA`, `ISA`, `CHILD`, `LCP_OVERFLOW` against its section bounds before use (a corrupt body then throws `InvalidDataException`, never reads outside the mapping); `Trusted` skips the checks and is opt-in for files the process wrote itself or verified at `Structure` level. This closes the weakness recorded in the K0 gaps for v6 ("query hot loops read raw pointers with Debug.Assert-only bounds ... Load validates header only"). Throughput gates in [PROGRESS-S1.md](../PROGRESS-S1.md) are measured in the default mode and `Trusted` is reported next to it (the overhead of `Checked` is UNVERIFIED).

`SuffixTreeSerializer.CalculateLogicalHash` accepts any `ISuffixTree` (`SuffixTreeSerializer.cs:30-47`); it must give the same value for an ESA and for the Ukkonen tree of the same text when requirement R-HASH of [ADR-001](ADR-001-index-contract.md) section 4 holds. `TEXT_SHA256` is the hash of the text alone and is a different value.

Writer atomicity: sections are streamed to `<path>.tmp` with the first 8 bytes zero; after all sections, the header (with CRCs) is written, the file is flushed, the magic is written **last**, then the file is renamed over `<path>`. A crash leaves either no file or a file without magic, which the dispatcher rejects as "unknown index kind". The writer builds sections directly from the build buffers ([ADR-002](ADR-002-sa-construction.md) section 3); no second in-memory copy.

## 7. Versioning and compatibility

* `VERSION_MAJOR` changes only for incompatible changes to the header or to the meaning of existing kinds; a reader refuses a larger major (`NotSupportedException` with both versions).
* `VERSION_MINOR` changes for additive features: new optional section kinds, new codecs of optional sections, new header bits in reserved space. A reader of the same major opens any minor and ignores unknown **non-required** sections; unknown **required** sections or unknown codecs of required sections are refused.
* The header page is 4096 bytes so up to 62 sections fit without moving anything.
* Endianness: little-endian only. A big-endian host throws `NotSupportedException` at open (zero-copy mapping cannot byte-swap). The v6 code makes the same assumption implicitly through raw pointer reads (`PersistentSuffixTreeNavigator.cs:139-140`, UNVERIFIED for big-endian hosts).

## 8. Load dispatch and relation to v6

New entry point, additive, outside the v6 classes: `SuffixIndexFactory.Open(string path, EsaOpenOptions? options = null)` (final name in WP 1b.7).
1. Read the first 8 bytes with `RandomAccess.Read` (no mapping yet).
2. If they equal the v6 magic bytes (the `int64` `0x5452454558494646` stored by the v6 writer) -> `PersistentSuffixTreeFactory.LoadPersistent(path)` (`PersistentSuffixTreeFactory.cs:130-146`, unchanged). Returns `PersistentSuffixTree` as `ISuffixTree`.
3. If they equal the SESA magic -> map the file, validate at the requested level, return `EsaSuffixTree`.
4. Otherwise -> `InvalidDataException("unknown suffix index kind")`.

Relation to v6, all of it tested (WP 1b.7):
* `PersistentSuffixTreeFactory.Load(esaFile)` keeps throwing its existing `InvalidOperationException("... Magic number mismatch.")` (`PersistentStorageHeaderValidation.cs:117-118`); the test pins this so nobody "fixes" the v6 loader by accident. `SuffixIndexFactory.Open(v6File)` returns the v6 tree. Existing v6 tests are unchanged and green.
* v6 stays the format of the online Ukkonen builder (and of S6 online append); SESA is the format of every static engine. No binary conversion exists or is planned: to convert, build an ESA from the v6 tree's text (`ISuffixTreeSearch.Text`) in seconds instead of the 14.9 s Ukkonen build (`BENCHMARK.md:22`).
* The logical export format v2 (`SuffixTreeSerializer`, text + hash) remains the portable exchange format for both; `Import` rebuilds with Ukkonen.
* v6 and SESA can be opened side by side in one process and queried through the same `ISuffixTree`, which is what the differential fuzz harness does.

## 9. Mapping and lifetime

`MemoryMappedFile.CreateFromFile(path, FileMode.Open, null, 0, MemoryMappedFileAccess.Read)`; one view of the whole file acquired with `SafeMemoryMappedViewHandle.AcquirePointer` and released in `Dispose`; sections are exposed as typed views over that pointer. All accessors check `ObjectDisposedException` at the public boundary only (not per element). Files larger than the address space are an S2 topic (out-of-core is untested in v6 as well, `BENCHMARK.md:188`). OS read-ahead is left to the platform: cold behaviour (first 100,000 queries after a cold load: 0.488 s vs 0.346 s warm for v6, `BENCHMARK.md:179`) is measured in WP 1b.9 and a managed pre-touch of sections (sequential read of the mapping) is the only warm-up allowed.

## 10. Risks and open points

* The `CHILD` bit packing is not designed here (P3). The container only fixes the section and the codec id.
* Whether `Checked` access costs more than a few percent is UNVERIFIED; if it does, the default stays `Checked` and the gate is reported in both modes.
* `LRS_RANK` / `NODE_COUNT` depend on the child-table pass; the writer computes them while building the table.
* Decision for the owner: extension and kind name (`.sesa`, "SESA"); 4096-byte header page versus a 256-byte header (the page costs 4 KiB, irrelevant above 1 MB).
