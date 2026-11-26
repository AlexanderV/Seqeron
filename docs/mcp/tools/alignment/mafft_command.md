# mafft_command

MAFFT 7.505 as a command line (the script mafft: options parsed as the script does, the pipeline of ported programs, FASTA / CLUSTAL / PHYLIP output, tree and distance files, the program's standard error); equal to the installed program on the committed command lines.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `mafft_command` |
| **Method ID** | `MafftProgram.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The script mafft of MAFFT 7.505 (Katoh & Standley 2013; mafft [options] in.fa > out) as a function of its command line and the text of the files it names (MafftProgram.Run): the options are parsed exactly as the script does (every option and the checks of its argument; unknown options and missing files print the script's usage text and return 1; --version, --help), the strategy is derived as the script does (the default FFT-NS-2, --retree n, --maxiterate n, --localpair / --globalpair / --genafpair = L-INS / G-INS / E-INS, --oldgenafpair, --weighti x (0 = the pair library without constraints), --auto (the rule of the script), --parttree / --dpparttree / --fastaparttree with --partsize / --groupsize, --fastapair / --fastswpair (F-INS-1 / H-INS-1: dndfast7 runs the ported FASTA 3.4 'fasta34 [-A] -z3 -m10 [-n] -Q -b500000 -E500000 -d500000' for every sequence against the following ones; its scores are the distances (hat2), its local alignments the pair library (hat3); a nucleotide set ends as the program does - 'count (fasta34 -z 3) = n' when a sequence has no alignment, the segmentation fault of dndfast7 when the alignment of the other strand is cut away), --treein file / --pileup / --randomchain / --large / --memsavetree / --memsavetreex / --stepadd / --youngestlinkage (the guide trees of loadtree, createchain and compacttree; with a pair library the program nodepair, --lhlimit n), --memsave / --nomemsave (MSalignmm, Falign_udpari_long; --auto for 20 000 sequences or more), --skipiterate x (dvtditr -E: fixed subalignments), --add file / --addprofile file (an alignment plus new sequences: disttbfast or tbfast -K -I, with --keeplength, --mapout, --compactmapout and --smoothing, with every pair library and guide tree above; with --maxiterate dvtditr -K), --core / --coreext / --corethr / --corewin (setcore: needs the file hat2 of the iterative methods or --distout), --aamatrix file (loadaamtx), --merge table (subalignments that keep their columns: the constrained UPGMA tree, the checks of the groups; disttbfast, tbfast and dvtditr), --text / --textmatrix file (an alphabet of 256 characters: every byte except blank, gap and line end is a letter; the entries of the file replace the scores of the matrix; CLUSTAL output marks identical columns only), --seed file / --seedtable file (aligned seed sequences: multi2hat3s writes the seeds and the blocks of every pair of them, disttbfast -E 1, tbfast -J and dvtditr -l align the input to them; several --seed options, --anysymbol), --addfragments / --addfull / --addlong / --addtotop new.fa (the program addsingle: the new sequences join the alignment one by one, with --6merpair / --10merpair / --6mermultipair / --10mermultipair / --multipair (the default) / --hybridpair / --localpair; --treein, --keeplength, --mapout, --treeout), --hybridpair without a fragment option (Hybrid-INS-1: pairlocalalign -Y writes the blocks of the old sequences against the added ones - nothing without --add -, disttbfast -E 1 -T -y the 6-mer distances of all sequences (hat2), tbfast loads hat3 and hat2; every run of tbfast with --scoreout ends with a segmentation fault of the program (sumofpairsscore after freeconstants) and --hybridpair --thread n without --add one of pairlocalalign -Y: both are refused), --codonpos file / --codonscore file (the codon part of Falign), --adjustdirection / --adjustdirectionaccurately (makedirectionlist and setdirection: sequences of the other strand get _R_ in front of their titles and are reverse complemented; with --merge the subalignments must agree), --anchors file (external anchors: readexternalanchors, pickpairanch, Falign_givenanchors; --anchors with an anchor behind the last residue of its sequence is refused, the program reads outside a table there), --sp (the dynamic program D__align / D__align_variousdist of Dalignmm.c: the gap penalty of a cell is weighted by the probability that the gap is new, from the table of the gap lengths of both groups; the program prints its progress and a check of the score to the standard error and stops with exit status 1 when the check fails; the programs splittbfast and addsingle do not know -@ and stop with 'illegal option @'; --sp with --add, which reads the matrix before its first row, and --maxiterate with --scoreout, which reads a freed block, end with a segmentation fault of the program and are refused), --shiftpenalty x (with --globalpair: the dynamic programs may jump to the best earlier cell at that price), --unalignlevel x / --allowshift (with --globalpair: the scoring matrix of a step is shifted by the offset of the distance, in the iterative refinement one matrix per class of distance), --nofft / --fft, --6merpair / --10merpair, --averagelinkage / --minimumlinkage / --mixedlinkage x, --unweight, --nuc / --amino, --op / --ep / --exp / --opdist, --bl / --jtt / --kimura / --tm / --kappa / --fmodel / --nwildcard / --nzero, --similaritylevel, --gop / --gep / --gexp / --lop / --lep / --lexp / --LOP / --LEXP, --leavegappyregion, --termgappenalty, --minimumweight, --randomseed, --noscore, --bunkatsunashi, --outnum, --scoreout, --thread / --threadtb / --threadit, --anysymbol / --preservecase, --quiet, --progress file, --out file) and run on the ported programs (replaceu, disttbfast, tbfast + pairlocalalign, splittbfast, dndpre, dvtditr, restoreu, f2cl) in the script's order. The output formats are the program's: FASTA of 60 residues per line (--inputorder, --reorder), --clustalout (CLUSTAL with the strategy name and the marks of f2cl), --phylipout, --namelength n; --treeout / --nodeout / --distout write the guide tree and distance files next to the input (name.tree, name.hat2) and are returned in files, as is --out. The arguments are those of the program, the LAST one is the input file; every file the arguments name is given as text in files (file name -> content); the file '-' is standard input (stdin). Options whose binaries are not ported yet (see the registered list in docs/Validation/review-2026-09/B09.md) are refused with a NotSupported error that names the option; options for programs that are not installed in the reference installation (--lastpair / --lastmultipair need lastdb / lastal and are registered with a Docker package in B09.md; --larapair, --foldalignpair, --scarnapair, --dafspair, --mccaskill and, without the parameter contrafold, --contrafold) give the script's own message 'Install ...' and return 1, as the reference does; with contrafold = true, --contrafold with a pair library runs (B09.md F377). Without arguments the program starts an interactive wizard (refused). The standard error equals the installed Debian mafft 7.505-1 byte for byte (the progress text with several threads, which differs from run to run in the program, is written as of one thread); standard output and files equal it on every case of scripts/alignment/generate_mafft_command_oracle.py and generate_mafft_compacttree_oracle.py (B09.md F315 ff.). A guide tree file with NaN or infinite branch lengths gives NaN weights: the dynamic programs then read cells of the trace-back matrix that an earlier alignment left behind, and the port keeps that matrix as the program does (B09.md F334). Inputs on which the program itself reads or writes outside its buffers, and --distout with the trees that name a random temporary directory or a process id in the program's message, are refused (B09.md F326 ff.); so is a protein run in which fasta34 takes a query for DNA (more than 85 percent of ACGTUN in more than 20 residues): the program then reads a score into an unassigned variable (B09.md F349).

## Core Documentation Reference

- Source: [MafftProgram.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/MafftProgram.cs)
- Algorithm: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md)
- Oracle: [generate_mafft_command_oracle.py](../../../../scripts/alignment/generate_mafft_command_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `arguments` | array<string> | Yes | The command-line arguments of mafft, one per element, the input file LAST (for example ["--localpair", "--maxiterate", "1000", "--clustalout", "in.fa"]). To split a command string like a shell, use the C# API MafftProgram.SplitCommandLine. |
| `files` | object | No | The text of the files the arguments name: file name -> content (the input FASTA, --add / --merge / --seed / --aamatrix / --treein files ...). |
| `stdin` | string | No | Text of standard input (read when the input file is '-'). |
| `contrafold` | boolean | No | Run as the Debian mafft 7.505-1 with CONTRAfold 2.02 installed beside contrafoldwrap (MafftInstallation.DebianWithContrafold): --contrafold with --globalpair / --localpair / --genafpair then runs contrafoldwrap, the ported CONTRAfold (single-precision inside / outside algorithm, posteriors >= 0.01) and the RNA step of tbfast (--weightr / --weightm / --rop); false (the default) is the plain Debian installation, which answers --contrafold with the script's 'Install CONTRAfold' text. (default: false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `exitCode` | integer | 0 = success, 1 = the script or a program stops with an error (usage text, unknown option, missing file, unknown letter ...). |
| `stdout` | string | The alignment in the output format (empty with --out). |
| `stderr` | string | The messages of the script and of the programs it runs (empty with --quiet), or the usage text. |
| `files` | object | The files the run wrote: --out, name.tree (--treeout / --nodeout), name.hat2 (--distout), name.notused (--sparsepickup), the --progress file. |

## Errors

| Code | Message |
|------|---------|
| 1001 | arguments is null or empty |
| 1002 | an option whose program is not ported yet (NotSupported, the message names the option) |
| 1003 | the file '-' is read from stdin |
| 1004 | a refused command line is an exit status 1 with the script's own message in stderr |

## Examples

### Example 1: Three proteins, CLUSTAL output

**User Prompt:**
> Align these three sequences with MAFFT and give the CLUSTAL format.

**Expected Tool Call:**
```json
{"tool": "mafft_command", "arguments": {"arguments": ["--quiet", "--clustalout", "ex.fa"], "files": {"ex.fa": ">a first\nMKVLAAGIVGLLLA\n>b\nMKVLAGIVGLLA\n>c\nMKVIAAGLVGLLLA\n"}}}
```

**Response:**
```json
{"exitCode": 0, "stdout": "CLUSTAL format alignment by MAFFT FFT-NS-2 (v7.505)\n\n\na               MKVLAAGIVGLLLA\nb               MKVL-AGIVGLLA-\nc               MKVIAAGLVGLLLA\n                ***: **:****  \n", "stderr": "", "files": {}}
```

The output equals `mafft --quiet --clustalout ex.fa` of the installed Debian mafft 7.505-1.

## See Also

- [pairwise_aligner](pairwise_aligner.md)
- [alignment_operation](alignment_operation.md)
