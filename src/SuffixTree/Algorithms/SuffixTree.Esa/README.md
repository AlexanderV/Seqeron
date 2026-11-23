# SuffixTree.Esa

An enhanced-suffix-array (ESA) backend for the suffix-tree API of this repository, for .NET 10. The index is a suffix array, its LCP array and a child table over
the text; it is node for node the terminated suffix tree of the Ukkonen implementation (same node counts, same `Traverse` event stream, same logical hash,
checked by the tests), so every `ISuffixTree` consumer works unchanged. It adds a 64-bit navigation and suffix-array surface (`ISuffixIndex`,
`ISuffixTreeNavigation`, `ISuffixArrayAccess`), and a persistent container, **SESA**, that is built by an external-memory pipeline with bounded RAM and read in
place through a memory mapping.

**Status.** Stage S1 of the suffix-tree upgrade (`docs/refactoring/suffix-tree-upgrade/`), measured and reviewed, not yet released as a package: no `PackageId`, no tag.
**There is no public entry that builds a container yet** (section "Building a container today"). Everything below is x86-64 Linux, .NET 10 unless it says otherwise.

**When to reach for it.** The in-memory `SuffixTree` (Ukkonen) takes about 169 bytes per character and 14.6 s for 10 Mbp of DNA; `SuffixTree.Persistent` (v6) takes 14.9 s and 56 bytes per
character on disk. The ESA takes 10 bytes per character (plus the 2-byte-per-character string it shares with the caller) and about one second per 10 Mbp in memory, answers counts and locates
within 1.15x of the fastest tools measured (count 647k against 747k per second, locate 990k against 1.14M patterns per second), and as a container is opened in milliseconds and read from disk with no managed heap. Use the Ukkonen tree when the text must grow online, or the text is a few
kilobytes and nothing else matters. Use `SuffixTree.Persistent` when your program must build a file through the public API today. Use the ESA for everything that is built once and queried a lot,
for texts that do not fit the heap (build the container offline), and for code that wants the navigation API.

| | `SuffixTree` (Ukkonen, in memory) | `SuffixTree.Persistent` (v6, mmap) | `SuffixTree.Esa` in memory | `SuffixTree.Esa` container (SESA, mmap) |
|---|---|---|---|---|
| build 10 Mbp random DNA, 1 thread | 14.6 s | 14.9 s | **1.138 s** [1.118-1.174] (public `CreateEsa` 0.841 s [0.832-0.907]) | 2.33 s [2.04-2.39] at 10^7, 4 threads; 17.9 s [15.9-19.4] at 10^8, 4 threads |
| size per character | 169 B (+2 B text) | 56.2 B (file) | **10.00 B** (+2 B for the caller's string it shares) | 9.25 B (core) / 13.25 B (with the Psi section) on disk |
| peak memory of the build | 1,761 MB | 1,179 MB (46 MB anonymous) | 153.8 MB (15.4 B/char) | RssAnon 253 MB [247-256] at 10^8 (286 MB on genome-like text; gate 320 MB), whatever n |
| count of 20-mers, 1 thread | 589,810 /s | 388,371 /s | 647,176 /s [611,792-824,474] | 1.00-1.03 of the in-memory tree on DNA, 4 threads (0.93 on a 20-letter text) |
| locate of 1,000 8-mers | 49,653 patterns/s | 39,943 patterns/s | 990,228 patterns/s [876,565-1,627,933] (`FindAllOccurrences`); 2,173,239 /s with the lazy `Locate` | as in memory |
| alphabet | UTF-16 | UTF-16 (ASCII stored in 1 B) | UTF-16 | **bytes only (0-255)** |
| longest text | 1,073,741,823 | 1,073,741,823 | 1,073,741,823 (a `string` cannot be longer) | 2,147,483,646 (a 32-bit surface; 64 bits is stage S2) |

Sources: the Ukkonen and v6 columns are `docs/refactoring/suffix-tree-upgrade/BENCHMARK.md` section 1 (2026-11-03); the ESA columns are the gate tables of ADR-004 section 5a "1b.9" (A1, A2, A3, A4, A6, C),
rows `analysis/bench/results_1b9_rows.tsv` (lines 9, 12, 22, 44, 202, 205, 206) and `results_1b9_chain.txt`; five fresh processes per row, median [min-max], 4 vCPU Intel Xeon under KVM, 15 GB. The machine
has two timing states (1.2-1.4x between groups of fresh processes), so read a ratio of less than 1.4 between two rows of different days as "the same". 10 Mbp is `dna_10m` of `analysis/bench`.

---

## Features

| Feature | Details |
|---|---|
| **`ISuffixTree` unchanged** | `Contains`, `FindAllOccurrences`, `CountOccurrences`, `LongestRepeatedSubstring`, distinct substrings, longest common substring, exact-match anchors, maximal exact / unique matches, maximal repeated pairs, `Traverse`, `PrintTree`; identical answers to the Ukkonen tree (a differential fuzz of 111,227 cases in 10 minutes, 0 failures). |
| **Long navigation surface** | `IndexLocus` (a suffix-array interval and a depth), `TryExtendRight` / `TryExtendLeft`, lazy paged `Locate`, `TryParent`, `TryChild`, `GetChildren`, `TrySuffixLink`, `EdgeSymbol`, `SuffixAt`, `RankOf`, `LcpAt`, `Lcp`, `Lca`. `Support(IndexOperation)` says which operations are native and which build something first. |
| **Container (SESA)** | One file, memory-mapped read-only, read in place; CRC-32C per section; three validation levels; `Checked` / `Trusted` access; deterministic bytes. |
| **Concurrent readers** | Any number of threads over one tree; `Dispose` while calls are in flight is safe (section "Threads and disposal"). |
| **Managed code only** | `SuffixTree.Esa` has no native binaries and no P/Invoke of its own. (The v6 loader that `SuffixIndexFactory.Open` can dispatch to, in the referenced `SuffixTree.Persistent`, has advisory `posix_madvise` / `PrefetchVirtualMemory` imports that a SESA file never reaches.) `IsAotCompatible` and `IsTrimmable` are set; a console program that builds an in-memory tree, opens a container at the `Structure` level, walks to a parent and calls the dispatcher published with `PublishAot` (linux-x64) without a trim or AOT warning and ran correctly (speed not measured). |

---

## Quick start

All samples below were compiled and run against this repository (see the footer), in a console project with implicit usings and `using SuffixTree; using SuffixTree.Esa; using SuffixTree.Indexing; using System.Text;`.

### An in-memory index of a string

```csharp
using EsaSuffixTree tree = EsaSuffixTreeFactory.CreateEsa("abracadabra");   // the tree shares the string, it does not copy it

tree.Contains("bra");                                   // true
tree.CountOccurrences("a");                             // 5
int[] where = tree.FindAllOccurrences("abra").Order().ToArray();   // [0, 7]: the library returns suffix-array order, sort if you need ascending
tree.LongestRepeatedSubstring();                        // "abra"
```

`EsaSuffixTreeFactory.Create(string)` returns the same tree typed as `ISuffixTree`, which is not `IDisposable`: use `CreateEsa` (or cast) to write `using`.

### Walk the tree: extend a match, count it, page through its positions

```csharp
IndexLocus locus = tree.Root;
foreach (char c in "abr")
{
    // never pass the same variable as the locus and as the result (see "Limits and guarantees", item 10)
    if (!tree.TryExtendRight(locus, c, out IndexLocus next)) { /* the string does not occur */ break; }
    locus = next;
}
long count = tree.Count(locus);                          // 2 occurrences of "abr"
Span<long> page = stackalloc long[8];
int got = tree.Locate(locus, skip: 0, page);             // up to 8 positions; call again with skip += got for the next page (suffix-array order)

Span<IndexLocus> kids = stackalloc IndexLocus[8];
int kidCount = tree.GetChildren(tree.Root, kids);        // the terminator child ($, key -1) first, then ascending by first symbol
int first = tree.EdgeSymbol(kids[1], 0);                 // 'a'
tree.TryExtendRight(tree.Root, "abra".AsSpan(), out IndexLocus abra);
tree.IsExplicitNode(abra);                               // true: "abra" is followed by different symbols
tree.TrySuffixLink(abra, out IndexLocus link);           // the node of "bra"; the first call builds the navigation pack (about 0.1 s per 10 Mbp)
```

### Matching operations

```csharp
foreach (var m in tree.FindMaximalExactMatches("cadabrax", minLength: 4))
    Console.WriteLine($"text@{m.PositionInText} query@{m.PositionInQuery} length {m.Length}");   // "cadabra" and "abra"
string common = tree.LongestCommonSubstring("xcadabrax");      // "cadabra"
long distinct = tree.CountDistinctSubstrings();                 // 54
```

These need suffix links, so the first call builds the navigation pack (the inverse suffix array and an LCP pyramid, 4.13 B/char, once, thread-safe).

### Open a container

```csharp
using EsaSuffixTree file = EsaSuffixTreeFactory.Open("genome.sesa");        // Header level, Checked access: the defaults
Console.WriteLine($"{file.TextLength} symbols, {file.CountOccurrences("GATTACA")} x GATTACA");

using var strict = EsaSuffixTreeFactory.Open("genome.sesa",
    new EsaOpenOptions { Validation = EsaValidationLevel.Structure, VerifyTextSha256 = true });   // re-derive everything once

using var fast = EsaSuffixTreeFactory.Open("genome.sesa",
    new EsaOpenOptions { Validation = EsaValidationLevel.Sections, Access = EsaAccessMode.Trusted });

ISuffixTree any = SuffixIndexFactory.Open("genome.sesa");                  // a SESA container or a v6 persistent file, told apart by the magic
((IDisposable)any).Dispose();
```

`SuffixIndexFactory.Open` returns `ISuffixTree`, which is not `IDisposable`: cast to dispose it. A v6 file ignores the options.

**The levels** (`EsaValidationLevel`, each includes the cheaper ones):

| level | what it checks | cost | what it does not catch |
|---|---|---|---|
| `Header` (default) | magic, version, endian probe, header CRC-32C, file size, every section inside the file, aligned, not overlapping, element counts and widths, the alphabet, the layout of the packed text | O(sections): 13-15 ms for the first `Open` of a fresh process (JIT), 7.4-8.6 ms with ReadyToRun, 13 us warm; 20 ms with a cold page cache | any damage to a section body |
| `Sections` | also the CRC-32C of every section payload and that every byte between the sections is zero | one sequential pass: 33 ms for the 132 MB container of 10 Mbp, about 0.6 s at 10^8 | a body that was damaged and then re-checksummed |
| `Structure` | also the arrays against each other and the text: the suffix array is a permutation and sorted, the LCP array and the child table equal the recomputed ones, the statistics, the Psi array, the symbol counts (and with `VerifyTextSha256` the text hash) | 1.5 s at 10 Mbp, 22-27 s at 10^8; it grows as n^2 / chunk (one pass over the suffix array per 4 Mi positions) | an adversarial container that is consistent in every check |

**The modes** (`EsaAccessMode`): `Checked` (default) validates every index that is read from the file before it is used and bounds every walk: a damaged body is an `InvalidDataException`, never a read outside the mapping, a hang or a
crash. `Trusted` leaves those checks out of the hot loops (0-6 % faster by the gate measurement): **a corrupt file may crash the process, hang it or return garbage**. Use `Trusted` for a file this process wrote itself or after a
`Structure` open. See "Limits and guarantees".

### Texts that are not bytes

A container holds **one byte per symbol** (code units 0-255, at most 256 distinct, packed in ceil(log2 sigma) bits). A UTF-8 file is therefore a text of its bytes, and a pattern has to be the same bytes:

```csharp
using EsaSuffixTree utf8 = EsaSuffixTreeFactory.Open("cafes.sesa");
string needle = Encoding.Latin1.GetString(Encoding.UTF8.GetBytes("café"));   // the two bytes C3 A9 as two chars
utf8.CountOccurrences(needle);       // counts the UTF-8 "café"
utf8.CountOccurrences("café");       // 0: U+00E9 as one char is the single byte E9, which does not occur
utf8.CountOccurrences("cafć");  // 0: a code unit above 255 never occurs in a container
```

The in-memory tree takes any `string` (UTF-16 code units, lone surrogates and U+0000 included; a surrogate pair is two symbols) and any alphabet up to 65,536 distinct code units.

### Threads and disposal

```csharp
using EsaSuffixTree shared = EsaSuffixTreeFactory.Open("genome.sesa");
Parallel.For(0, 8, _ => { for (int i = 0; i < 1000; i++) shared.CountOccurrences("GATTACA"); });   // one tree, any number of readers

shared.Dispose();
shared.Contains("A");   // ObjectDisposedException (before any argument check), also from a thread that was already inside another call
```

Every read is safe for concurrent callers (the navigation pack is built once under a lock; the longest repeated substring is cached with a compare-exchange). `Dispose` is idempotent, does not wait, and
**a call that was in flight when `Dispose` ran completes with its answer**; the mapping is released by the last call to leave. A call that starts after `Dispose` throws `ObjectDisposedException`. The
results that read the mapping later (`Text`, `ISuffixIndex.Text`, the sequence of `EnumerateSuffixes`) count themselves as a call at every access. The protocol and its proof are in ADR-004 section 5a "1b.9g".

---

## Building a container today

`SesaBuilder` is `internal`: **a program that references this library cannot build a container.** The ways that exist are the benchmark application and the hand-off package, both in this repository:

```
dotnet run -c Release --project apps/SuffixTree.Esa.Bench -- sesa build --input text.bin --work /big/disk/work --keep true
#   -> /big/disk/work/out.sesa   (the input is a file of bytes, one symbol per byte; out.sesa, out.sesa.tmp and out.sesa.tmp.d in the work directory are deleted first)
dotnet run -c Release --project apps/SuffixTree.Esa.Bench -- sesa resume build --input text.bin --work /big/disk/work --resume true --keep true   # the same with Resume on
analysis/s1b-handoff/s1b.sh        # the whole chain (generate, build, verify, Structure level, queries) for 10^6 .. 3x10^9 symbols; S1B-RUNBOOK.md next to it
```

What a build needs: four worker threads by default; a free disk of 25-26 bytes per character during the build (the container reserves 13.25-14.25, the scratch up to 12); a page cache of about the size of the packed text plus the
difference-cover rank file (about 0.4-0.5 bytes per character) or the build slows by 1.7x at 0.7 of that floor and 5x at 0.5 (ADR-004 "1b.9d"); a peak anonymous memory of at most 320 MB at every n. A failed build removes what it made (an existing
target stays as it was) unless `Resume` is on, which keeps `state.json` and the finished stages for the next call.

---

## The container format in one paragraph

A SESA file is a header page of 4,096 bytes (a 128-byte header with the magic `89 'S' 'E' 'S' 'A' CR LF SUB` written last, the text length, node count, longest repeated substring, the text's SHA-256 and a CRC-32C of
the header page and the section table, then a table of up to 62 entries of 64 bytes) followed by sections, each at a 64-byte aligned offset (4,096 from 1 MiB), all little endian with no pointers: the packed TEXT (b bits per
symbol, big-endian bit stream, 64 bytes of zero padding), the ALPHABET, the symbol COUNTS, the suffix array SA (4 bytes per rank, `SA[0] = n`), the LCP array (1, 2 or 4 bytes per entry with an overflow list), the child table
CHILD (4 bytes per rank, layout A, in a tagged or a packed codec), the optional PSI array and BUILD_INFO, each with its own CRC-32C. The same input and options give the same bytes. The full definition, the validation rules and the
amendments are in [ADR-003](../../../../docs/refactoring/suffix-tree-upgrade/adr/ADR-003-esa-container-format.md).

---

## Performance

Measured on one machine (4 vCPU Intel Xeon under KVM, 15 GB, Ubuntu 24.04, .NET 10.0.403, Release, defaults, tiered compilation and dynamic PGO on). Five fresh processes per row, median [min-max]. Warm unless stated.

| what | condition | figure | row |
|---|---|---|---|
| in-memory build, SA + LCP + child table | 10 Mbp random DNA, 1 thread | 1.138 s [1.118-1.174]; public `CreateEsa` 0.841 s [0.832-0.907]; 36 cold builds of the day 0.927 s; same-session libsais (SA only) 0.344 s | `results_1b9_rows.tsv:9,44`, `results_1b9_libsais.tsv` |
| in-memory build at 50 Mbp | random DNA | 5.23 s [4.52-6.82] | `results_1b9_rows.tsv:173` |
| size in memory | `dna_10m` | 10.00 B/char (text 1 + SA 4 + LCP 1 + child 4); `rep_10m` 11.00; `a^n` or `abab...` 13.00 (the LCP array needs 4 bytes per entry; measured by the reviewer at 5 Mbp) | `:22`, ADR-004 1b.9 A2 |
| container build, production chain | 10^8 symbols, 4 threads, fsync on | 17.9 s [15.9-19.4]; RssAnon 253 MB [247-256]; 13.25 B/char; 10 of 10 chains passed their 11 checks | `results_1b9_chain.txt` |
| count of 20-mers | 10 Mbp DNA, 1 thread / 4 threads, one shared tree | 647,176 /s [611,792-824,474] / 3,098,222 /s (3.96x of 783,122 on the same instrument) | `:202`, `results_1b9_abapp.txt` |
| locate of 8-mers | 1,000 patterns, 1 thread / 4 threads | 990,228 patterns/s (`FindAllOccurrences`), 2,173,239 (lazy `Locate`) / 4,841,674 | `:205,206` |
| MEM / MUM / anchors of a 100 kbp query | 10 Mbp text | 0.109 / 0.083 / 0.083 s first call, 0.075 / 0.051 / 0.056 s steady | ADR-004 1b.9 A8 |
| file-backed against in-memory tree | 4 threads, `Checked` | DNA 1.02-1.03 (count), 1.00 (locate); 20-letter text 0.934 (count), 0.967 (locate) | ADR-004 1b.9 table "file against memory" |
| `Trusted` against `Checked` | file-backed tree | 0-6 % faster | ADR-004 1b.9 A7 |
| first `Open` of a fresh process | 132 MB container | 15.2 ms [13.1-15.4] (JIT); 8.5 ms [7.4-8.6] with ReadyToRun; 12.9 us for the next `Open` in the same process | `results_1b9_open_r2r.txt` |
| first queries of a process | `Checked`, mapped | the first 10,000-20,000 run as tier-0 code, 10-40 ms in all | ADR-004 1b.9g |

Not measured: arm64, Windows, macOS, a text above 10^8 symbols in memory, the cold-cache latency of a query on a mapped file, the speed of a NativeAOT build, the same on a second machine.

---

## Limits and guarantees

What the code promises, what it does not, and where each statement is pinned. "ADR" is `docs/refactoring/suffix-tree-upgrade/adr/ADR-004-external-memory-construction.md` section 5a unless named.

| # | limit or guarantee | why | pinned by |
|---|---|---|---|
| 1 | **`Trusted` access can crash the process, hang it or return garbage on a corrupt file, by contract.** In 1,500 random damages of a container, 47 killed the process and 4 hung it in `Trusted` mode; the same 1,500 were refused or answered without a fault in `Checked` mode. | the hot loops have no bounds checks | `EsaAccessMode.Trusted` documentation; ADR "1b.8" (the corruption fuzz reads damaged files in `Checked` mode only) |
| 2 | **`Checked` never reads outside a section**: a damaged body is an `InvalidDataException`, at every validation level. It can still return a wrong answer when the damage is in range (a swapped suffix-array entry): only `Structure` catches that. | every value read from the file that is used as an index is validated first | `SesaReaderCorruptionTests`, the corruption lane of the harness (105,044 cases), ADR "1b.7f" |
| 3 | **What each level detects** is in the table above. `Structure` is a check for accidental damage: a container that is consistent in every check cannot be told from a valid one. | the checks compare the arrays with each other | `SesaVerifier` documentation, ADR "1b.7f" threat model |
| 4 | **A file that is truncated, replaced in place or deleted and rewritten while it is mapped can fault the process** (a bus error, `AccessViolationException`), in either mode. Replace a container by renaming a new file over it, never by writing into it. | a memory mapping reads the file system's pages; a bound check cannot see the file shrink | ADR "1b.7f" limit (d) |
| 5 | **Longest text: 2,147,483,646 symbols** through `EsaSuffixTree` (a 32-bit surface; `NotSupportedException` above it naming the long API). The `Structure` level takes at most 2,147,483,645. `NodeCount` is an `int` and throws `NotSupportedException` when a tree has more nodes. The builder itself accepts up to 4,294,967,292 symbols, but no reader opens a container above 2,147,483,646. `ISuffixTree` positions are `int`; the 64-bit stage is S2. | positions and ranks are 32-bit in the arrays | `EsaSuffixTree.MaxTextLength`, `SesaHeaderReader`, ADR "1b.9f" |
| 6 | **Containers hold byte texts only** (code units below 256, at most 256 distinct). A container whose ALPHABET has a letter above 255 is refused as damaged; there is no encoding option in the build commands, so a UTF-8 or UTF-16 text is indexed as its bytes. A query with a code unit above 255 finds nothing in a container (count 0, `Contains` false, `TryExtendRight` false). The in-memory tree and the old v6 format keep UTF-16. | the packed text stores b <= 8 bits per symbol | `SesaHeaderReader` ("the packed text holds code units below 256 only"), ADR-003 amendment of 2026-11-18 |
| 7 | **A locus is checked in O(1) and not exactly.** A tuple that cannot be a locus of this text (out of range, depth not matching its interval) is an `ArgumentException`; a locus of another text of the same size, or a hand-made one that has the shape of a locus, is accepted and gives a wrong answer. Pass loci only to the tree that produced them. | the contract leaves interpretation to the engine and an exact check costs a table walk | `EsaSuffixIndex.CheckLocus` and `Resolve`; demonstrated in `analysis/review/1b10_review.md` |
| 8 | **`Dispose` of an in-memory tree releases nothing at once.** It raises the flag (every later call throws `ObjectDisposedException`); the arrays and, once built, the pinned navigation pack are freed when the tree is unreachable and its finalizers have run. Only a tree over a container frees its mapping at `Dispose` (when the last call has left). | the arrays are managed memory that a call in flight reads without a lock | `EsaSuffixTree.Dispose` documentation; `FacadeLifetimeTests` |
| 9 | **No `fallocate`, no reserved disk.** On Linux the build reserves the container with `fallocate` (through `FileStreamOptions.PreallocationSize`), so a full disk at the start is an `IOException`, and it checks the free space (the container plus 12 bytes per character of scratch) before it starts. Where the platform or the file system cannot reserve (the BCL ignores the request on some file systems; macOS and Windows take the explicit path, not run here), the file is sparse and a disk that fills during the build is a fault of the process, not an exception. The growth of the LCP section is reserved by writing zeros on every platform. | the BCL reserves blocks only for a file it creates | ADR "1b.9d", `FullDisk1b9dTests`, `Reservation1b9dTests` |
| 10 | **Do not pass one variable as both the locus and the result** of `TryExtendRight` or `TryExtendLeft` (`tree.TryExtendRight(l, c, out l)`): the methods clear `result` before they read `locus`, so the call returns a wrong answer (`false` for a string that occurs) or throws `ArgumentException`. Write the result to a new variable. The other navigation members are not affected. | `in` and `out` parameters alias the same variable | `analysis/review/1b10_review.md` (R01) |
| 11 | **`Resume` is off by default.** Without it an interrupted build leaves `<target>.tmp` and `<target>.tmp.d`, and the next build starts from the beginning after deleting them. With it the build writes `state.json`, `<target>.lock` and keeps the finished stages (a resume saves 70-94 % of a restart after a kill from the end of the LCP stage on); a failed resumable build keeps 13 bytes per character on the disk. **Two builds of one target without `Resume` are not supported**: the second deletes the first one's temporary file, and the first one's rename can publish the second one's half-written file. There is no directory `fsync`: a power loss can cost the newest stage or the final rename, never a wrong container (argued, not tested on a device). | the BCL cannot open a directory for `fsync`; the lock is only taken by resumable builds | ADR "1b.9e", "1b.9e.1", debt register G(1), G(4) |
| 12 | **The builder is internal** (section "Building a container today"). | S1 delivered the engine; the public entry is S8 | `SesaBuilder` (internal); `analysis/review/1b10_review.md` |
| 13 | **What ran where.** Measured and tested on x86-64 Linux only. The arm64 CRC-32C path was compiled and checked against the table path in software, never run on arm64; Windows and macOS were not run; the guard-page tests need Linux `mmap` / `mprotect` (they are skipped elsewhere); a big-endian host is refused (`NotSupportedException`); a 32-bit process was not tried. The hand-off package for the owner's Mac (arm64, 3x10^9 symbols) is delivered, its run is pending. | no other machine was available | ADR "1b.9f"; `analysis/s1b-handoff/S1B-RUNBOOK.md` |
| 14 | **The first `Open` of a fresh process is JIT-bound**: 13-15 ms against the 10 ms gate of the plan (7.4-8.6 ms with ReadyToRun of `SuffixTree.Esa.dll`, 13 us for the next `Open`). The gate is recorded as not met on the plain JIT. | the reader code is compiled at its first call | `analysis/bench/results_1b9_open_r2r.txt`, ADR "1b.9" A9 |
| 15 | **The gate at four threads on a 20-letter text is not met as written**: the file-backed tree counts at 0.934 of the in-memory tree (locate 0.967). The facade is not the cause (it is 0.99-1.015 of the index behind it): the 5-bit packed symbol is decoded at every step of the descent. DNA is at 1.02-1.03. | a text of 20 letters needs 5 bits per symbol in the file | ADR "1b.9" open point 1, "1b.9g.1" |
| 16 | **The reader protocol has a cost and platform assumptions.** Calls publish themselves with plain stores (no atomic, no fence); `Dispose` runs one process-wide memory barrier, which interrupts every core that runs a thread of the process (0.1 us when idle, 7-11 us with threads spinning, milliseconds when a core is descheduled), so opening and disposing thousands of containers per second is not free. Its safety rests on `Interlocked.MemoryBarrierProcessWide` being a real barrier (Linux `membarrier`, kernel 4.14 or later, or the runtime's own fallback). The x64 half of the proof was run; the weak-memory (arm64) half is an argument. `Dispose` does not wait: a reader that never returns holds the mapping. | the asymmetric protocol moves the cost from every call to `Dispose` | ADR "1b.9g" (lemmas 1-3, "Limits of the protocol") |
| 17 | **Cold-cache latency is not measured.** Every throughput figure is warm; the first queries of a cold mapped file pay page faults. | not asked of the gate | ADR "1b.9" limits |
| 18 | **Memory is 10 B/char in memory plus the caller's string; the navigation pack adds 4.13 B/char** the first time a suffix-link operation runs (`Parent`, `SuffixLink`, `ExtendLeft`, `RankOf`, `Lcp`, `Lca`, anchors, MEM, MUM, longest common substring). A repetitive text costs more (`rep_10m` 11.00, `a^n` 13.00 B/char in memory, 17.1 with the pack). `Support(IndexOperation)` reports `Degraded` until the pack exists. | the pack is built on first use | `EsaSuffixIndex` documentation, ADR-001 section 6 |
| 19 | **Result order.** `FindAllOccurrences` and `Locate` return positions in suffix-array order, not ascending; the representative of `LongestRepeatedSubstring` on a tie is the lexicographically smallest, and the text position reported for a longest common substring or an anchor is the first occurrence the navigator offers (another valid one than the Ukkonen tree's, for example). | the suffix array is the order | `ISuffixTreeSearch` remarks; differential tests compare these members as sets or by property |
| 20 | **The `Structure` open can write a scratch file**: when the lcp-interval stack of the structure pass gets deeper than 4 Mi entries (a text like `a^n` of more than 4 million symbols) it spills to `esa-verify-*.stack` in the system temporary directory, and removes it at the end. | the stack must not grow with the depth | `SesaVerifierTuning`, `SesaVerifier` |
| 21 | **Caller errors on a container-backed `Text` are reported as damage.** `file.Text[file.Text.Length]`, `Text.Substring` and `Text.Slice` with an out-of-range argument throw `InvalidDataException("... the file is damaged")` where the in-memory tree throws `IndexOutOfRangeException` / `ArgumentOutOfRangeException`. | the packed-text view uses one check for both | `analysis/review/1b10_review.md` (R06) |

---

Samples compiled and run on commit 04620cd5 (the head of branch `claude/stoic-maxwell-0olr2z` when the samples were run; `src/` and `apps/` are those of `12495c3e`), .NET 10.0.403, x86-64 Linux: the six `csharp` blocks of this file were extracted and run verbatim, in order, against containers built with the Bench command.
Tests: `tests/SuffixTree/SuffixTree.Esa.Tests` (default set 4,945 passed + 37 skipped in 14 min 31 s, Release; categories Exhaustive, Fuzz, Vectors and the checked build `-p:EsaBounds=true` are run by hand, see `docs/refactoring/suffix-tree-upgrade/PROGRESS-S1.md`).
Design and evidence: `docs/refactoring/suffix-tree-upgrade/` (`PLAN.md`, `BENCHMARK.md`, `adr/ADR-001` to `ADR-004`).
