// SPDX-License-Identifier: GPL-3.0-or-later
// This file is part of Seqeron.Mcp.Last, licensed under the terms stated in LICENSE (this directory).
using System.ComponentModel;
using ModelContextProtocol.Server;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Last;

namespace Seqeron.Mcp.Last.Tools;

/// <summary>The outcome of a run of the script mafft: exit status, standard output, standard error and the files the run wrote.</summary>
public sealed record MafftCommandLastResult(int ExitCode, string Stdout, string Stderr, Dictionary<string, string> Files);

public partial class LastTools
{
    /// <summary>The installation of the oracle (Debian mafft 7.505-1) with LAST 1256 in /opt/last/bin.</summary>
    private static readonly MafftInstallation WithLast = MafftInstallation.Debian with { Last = new LastPrograms() };

    [McpServerTool(Name = "mafft_command_last", Title = "MAFFT 7.505 with LAST 1256 - command line (--lastpair, --lastmultipair)", ReadOnly = true)]
    [Description("The script mafft of MAFFT 7.505 as a function of its command line and the text of the files it names (MafftProgram.Run), with the programs lastdb and lastal of LAST 1256 (ported, GPL) installed: --lastpair aligns every pair of nucleotide sequences with lastal (Last-INS-1). pairlocalalign -R of MAFFT 7.505 overflows its heap when three or more sequences are related (valgrind: invalid writes behind the rows of its table of local homologies; the run then crashes): such runs are refused with an error, the runs with two sequences and the runs without hits are answered. --lastmultipair stops with the messages of the script (it needs --addfragments, and dndpre does not know -Q).")]
    public static MafftCommandLastResult MafftCommandLast(
        [Description("The command-line arguments of mafft, one per element, the input file LAST (for example [\"--lastpair\", \"--clustalout\", \"in.fa\"]).")] string[] arguments,
        [Description("The text of the files the arguments name: file name -> content (the input FASTA, --addfragments / --add / --merge / --seed files ...).")] Dictionary<string, string>? files = null,
        [Description("Text of standard input (read when the input file is '-').")] string? stdin = null,
        CancellationToken cancellationToken = default)
    {
        if (arguments == null || arguments.Length == 0)
            throw new ArgumentException("arguments must contain the command-line arguments of mafft (the input file last).", nameof(arguments));
        var result = MafftProgram.Run(arguments, WithLast, files, stdin, cancellationToken);
        return new MafftCommandLastResult(result.ExitCode, result.Stdout, result.Stderr, new Dictionary<string, string>(result.Files));
    }
}
