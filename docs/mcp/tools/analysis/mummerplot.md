# mummerplot

MUMmer mummerplot: the gnuplot script and the .fplot / .rplot / .hplot / .filter files of alignments, byte for byte.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `mummerplot` |
| **Method ID** | `MummerPlot.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

MUMmer 3.23 mummerplot (perl script 3.5): the gnuplot plot of the alignments of a nucmer / promer .delta, a .cluster (nucmer --nodelta), a mummer -c output (.out; use -c -b -L) or a show-tiling file (.tiling; coverage plot): the data files <prefix>.fplot (forward hits sorted by percent similarity, 'sR sQ sim' per end of an alignment, blocks separated by two empty lines), <prefix>.rplot (reverse hits), <prefix>.hplot (highlighted hits: breakpoints with breakLength = -b, SNP positions with snp = -S through show-snps) and the script <prefix>.gp (terminal x11 / postscript / png with --x11 --postscript --png or terminal = -t, size small / medium / large, ranges, tics of several sequences, grid, labels, the line styles, the plot or splot command; color / coverage / title / xRange / yRange / layout options) - and, for delta input with filter (-f) or layout (-l, --fat; needs rFile / qFile), the filtered delta file <prefix>.filter (delta-filter -r -q). For delta and cluster input the sequence files named in the first line of the file are read (files: name -> FASTA text) to get the lengths and offsets of the sequences, rFile / qFile (-R / -Q) name files with an ordered set of sequences (a multi-FASTA file, or lines 'id length [+|-]'), idR / idQ (-r / -q) plot one sequence. The script then runs gnuplot: gnuplotCommand is the command line it runs (the rendering - .ps / .png / the x11 window - is not part of the text output); stdout holds the line of gnuplot --version (gnuplotVersion; gnuplotStatus 0 = found, -1 = not found, other = outdated), messages the progress lines and warnings (multiple terminals / sizes, -f -l -S with non-delta input, ports of non-x11 plots, overlaid sequences, 'Forking mouse listener' and the connection warnings of the listener of the x11 plot with rPort / qPort). A die of the script (invalid terminal / size / range, empty or unrecognised match file, no alignment data, an unreadable -R / -Q file, a line that cannot be parsed) is reported on messages with a non-zero exitCode (the status the program returns is the stale errno of its last system call: the tool returns 2). When equal lengths or offsets make perl's hash order decide the layout of several sequences the program's own output differs from run to run; the port uses the order of first appearance. Byte for byte equal to the installed mummerplot on 330 + runs (every input type, every option) and 1 000 more random runs.

## Core Documentation Reference

- Source: [MummerPlot.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/MummerPlot.cs)
- Algorithm doc: [Mummer_Package_Programs.md](../../../algorithms/Alignment/Mummer_Package_Programs.md)
- Oracle: [generate_mummerplot_oracle.py](../../../../scripts/alignment/generate_mummerplot_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `match` | string | Yes | The text of the match file (.delta, .cluster, mummer -c .out or .tiling). |
| `files` | object | No | The files the script reads besides the match file, by name: the sequence files named in the first line of a delta / cluster file and the -R / -Q files (FASTA text or lists of 'id length [+\|-]'). |
| `matchName` | string | No | The file name of the match file as given to the script (printed in the messages; default out.delta). (default: out.delta) |
| `breakLength` | integer | No | -b: highlight alignments with breakpoints further than this many nucleotides from the nearest sequence end. |
| `color` | boolean | No | --color / --nocolor: colour the plot lines by percent similarity (true) or turn all colour off (false). |
| `coverage` | boolean | No | -c / --nocoverage: draw a reference coverage plot (default for tiling input). |
| `filter` | boolean | No | -f: show only the best alignments of a delta file (delta-filter -r -q), delta input only. |
| `layout` | boolean | No | -l: lay a delta multiplot out in an intelligible fashion (needs rFile / qFile), delta input only. |
| `prefix` | string | No | -p: prefix of the output files (default out). |
| `reverseVideo` | boolean | No | -rv: reverse video for x11 plots. (default: false) |
| `idR` | string | No | -r: plot this reference sequence id. |
| `idQ` | string | No | -q: plot this query sequence id. |
| `rFile` | string | No | -R: file (a key of files) with an ordered set of reference sequences. |
| `qFile` | string | No | -Q: file (a key of files) with an ordered set of query sequences. |
| `rPort` | integer | No | -rport: port to send the reference id and position to on a mouse double click (x11 plots). |
| `qPort` | integer | No | -qport: port for the query. |
| `size` | string | No | -s: 'small' (default), 'medium' or 'large'. |
| `small` | boolean | No | --small. (default: false) |
| `medium` | boolean | No | --medium. (default: false) |
| `large` | boolean | No | --large. (default: false) |
| `snp` | boolean | No | -S: highlight the SNPs (delta input; reads the sequence files). |
| `terminal` | string | No | -t: terminal 'x11' (default), 'postscript' or 'png'. |
| `x11` | boolean | No | --x11. (default: false) |
| `postscript` | boolean | No | --postscript. (default: false) |
| `png` | boolean | No | --png. (default: false) |
| `title` | string | No | --title: gnuplot plot title. |
| `xRange` | string | No | -x: x range '[min:max]'. |
| `yRange` | string | No | -y: y range '[min:max]'. |
| `fat` | boolean | No | --fat: lay out by the fattest alignment only (implies layout). (default: false) |
| `gnuplotStatus` | integer | No | The status of 'gnuplot --version': 0 found (default), -1 not found, other = outdated. |
| `gnuplotVersion` | string | No | What 'gnuplot --version' prints (default 'gnuplot 6.0 patchlevel 0'). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `files` | object | The files of the script by name (.fplot, .rplot, .hplot, .gp, .filter). |
| `stdout` | string | Standard output (the version line of gnuplot). |
| `messages` | string | Standard error of the script. |
| `exitCode` | integer | The exit status (0, or non-zero after a die). |
| `gnuplotCommand` | string | The gnuplot command line the script runs. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Match file text cannot be null or empty |
| 1001 | ERROR: Invalid terminal type, <name> |
| 1001 | ERROR: Could not read <file>, Unrecognized file type |

## Examples

### Example 1: A PNG plot of one alignment

**User Prompt:**
> Write the gnuplot script that plots this delta file as a PNG.

**Expected Tool Call:**
```json
{"tool": "mummerplot", "arguments": {"match": "ref.fa qry.fa\nNUCMER\n>ref qry 260 240\n11 250 1 240 3 3 0\n131\n-51\n0\n", "matchName": "x.delta", "files": {"ref.fa": ">ref demo\nggcccagtccagatcctcggaagtcccaactaacgaataagtagatccttctaaatagta\ngtatacgaatgtaatcaacatattaccgcacccaactttgattgggtcataaacaaacat\ncatgtttctggaggcggtttgtgagaccgacaaacgcgatagagtaatcttataaatcag\ngcaaatactcggtcctttctacgggtcgcggccgcgtggtcataggcagttttatgtttc\ntttgaggttaatgggaggca\n", "qry.fa": ">qry demo\nagatcctcggaagtcccaactaacgaataagtagatccttctaaatagtagtatacgaat\nataatcaacatattaccgcacccaactttgattgggtcataaacaaacatcatgtttctg\ngaggcggttttgagaccgacaaacgcgatagagtaatcttataaatcaggcaaatactcg\nggtcctttctacgggtcgcggccgcgtggtcataggcagttttatgtttctttgaggtta\n"}, "prefix": "demo", "png": true}}
```

**Response:**
```json
{"files": {"demo.fplot": "#-- forward hits sorted by %sim\n0 0 0\n0 0 0\n\n\n11 1 98.7551867219917\n250 240 98.7551867219917\n\n\n", "demo.rplot": "#-- reverse hits sorted by %sim\n0 0 0\n0 0 0\n\n\n", "demo.gp": "set terminal png tiny size 800,800\nset output \"demo.png\"\nset size 1,1\nset grid\nunset key\nset border 15\nset tics scale 0\nset xlabel \"ref\"\nset ylabel \"qry\"\nset format \"%.0f\"\nset mouse format \"%.0f\"\nset mouse mouseformat \"[%.0f, %.0f]\"\nset xrange [1:260]\nset yrange [1:240]\nset style line 1  lt 1 lw 3 pt 6 ps 1\nset style line 2  lt 3 lw 3 pt 6 ps 1\nset style line 3  lt 2 lw 3 pt 6 ps 1\nplot \\\n \"demo.fplot\" title \"FWD\" w lp ls 1, \\\n \"demo.rplot\" title \"REV\" w lp ls 2\n"}, "stdout": "gnuplot 6.0 patchlevel 0\n", "messages": "Reading delta file x.delta\nWriting plot files demo.fplot, demo.rplot\nWriting gnuplot script demo.gp\nRendering plot demo.png\n", "exitCode": 0, "gnuplotCommand": "gnuplot demo.gp"}
```

mummerplot -p demo --png x.delta: one forward alignment at 98.755 percent similarity; the axes are the reference and query lengths.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
