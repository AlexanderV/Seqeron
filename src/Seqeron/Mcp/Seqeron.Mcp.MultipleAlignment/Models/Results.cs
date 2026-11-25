// SPDX-License-Identifier: GPL-3.0-or-later
// This file is part of Seqeron.Mcp.MultipleAlignment, licensed under the terms stated in LICENSE (this directory).
namespace Seqeron.Mcp.MultipleAlignment.Tools;

/// <summary>Kalign 3.3.5 alignment: the rows of the sequences that were aligned, which input sequences those are and the alphabet the program detected.</summary>
public record KalignResultDto(
    string[] AlignedSequences,
    string Consensus,
    int TotalScore,
    double ExactTotalScore,
    int[] KeptIndices,
    string Alphabet,
    int[][] GuideTree);

/// <summary>Kalign 3.3.5 command: exit code, standard output, the files the program wrote (name to text) and what it printed to the error stream.</summary>
public record KalignCommandResultDto(
    int ExitCode,
    string Stdout,
    Dictionary<string, string> Files,
    string[] Messages);

/// <summary>POA 2.0 alignment: the rows (sequences in the program's output order, then the consensus sequences of -hb), the input index of every row (-1 for a consensus), and the PO, PIR and CLUSTAL texts.</summary>
public record PoaResultDto(
    string[] Names,
    string[] AlignedSequences,
    int[] Order,
    string Po,
    string Pir,
    string Clustal,
    string Stdout);
