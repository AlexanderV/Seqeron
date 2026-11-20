# Mash text commands — `bounds`, `dist`, `triangle`

Implementation: `Seqeron.Genomics.Analysis/MashCommands.cs` (`MashCommands.Bounds`, `Dist`, `Triangle`) on top of the MinHash sketching of `KmerAnalyzer` (`CreateMinHashSketch`, `CompareMinHashSketches`: Ondov et al. 2016, Genome Biol 17:132; Mash 2.3 `Sketch.cpp`, `hash.cpp`, `MinHashHeap.cpp`, `CommandDistance.cpp`). MCP tools `mash_bounds`, `mash_dist`, `mash_triangle`.
Audit item A5-10 (B09 F174). Sources read: marbl/Mash v2.3 `CommandBounds.cpp`, `CommandDistance.cpp`, `CommandTriangle.cpp` (raw.githubusercontent.com, User-Agent `Seqeron-review/1.0`).

## What is ported

| Command | Options | Output |
|---|---|---|
| `mash bounds` | `-k` (1-32), `-p` (0-1) | the two tables (Mash distance, screen distance) of the error of the lower Jaccard bound for sketch sizes 100 ... 10^6 and distances 0.05 ... 0.4: the smallest `x` with `BinomialCdf(x, s, m2j) > (1-p)/2` (the program scans `x` linearly with `gsl_cdf_binomial_P`; the port bisects the monotone CDF with `StatisticsHelper.BinomialCdf`), `j2m = -ln(2je/(1+je))/k` (screen: `1 - je^(1/k)`), printed as `j2m - dist` |
| `mash dist` | `-k -s -S -n -i -t -v -d` | one line per (query, reference) pair that passes `-v` (p-value) and `-d` (distance): `reference TAB query TAB distance TAB p-value TAB shared/total`; `-t`: the `#query` header and a row of distances per query (blank when a pair fails) |
| `mash triangle` | `-k -s -S -n -i -E -v -d` | Phylip lower-triangular matrix (`TAB count`, then `name TAB d ...`) or, with `-E` (implied by `-v` / `-d`), the edge list |

Numbers are printed as C++ `cout << double` does (`%g`, six significant digits, `inf` / `nan`). With `-i` a sequence shorter than `k` is not sketched (no row for it).

## Verification

`scripts/alignment/generate_mash_commands_oracle.py` runs the installed `mash` 2.3 on 220 random cases per command (families of related genomes 0.5-25 % divergent, unrelated ones, 1-3 records per file, `k` 4-32, sketch sizes 5-2000, seeds, `-n`, `-i`, `-t`, `-E`, `-v`, `-d`; for `dist` the references are sketched first with `mash sketch -o`), the whole stdout is compared as text (`TestData/Alignment/mash_commands_oracle.jsonl.gz`, `MashCommands_Oracle_Tests`): 660 / 660 equal.

## Not ported (registered, B09.md A5-10)

`mash screen` (containment of query sketches in a read mixture), the reads mode of `sketch` / `dist` (`-r -m -b -c -g`: noise filter, Bloom filter, coverage, genome size), the alphabet options (`-a`, `-z`, `-Z`; protein sketches) and `mash taxscreen` / `paste` / `info` (sketch file I/O).
