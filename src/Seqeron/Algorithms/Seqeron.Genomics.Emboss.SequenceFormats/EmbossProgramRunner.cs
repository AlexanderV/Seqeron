// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Analysis;

/// <summary>What an EMBOSS program run leaves: the standard output, the error stream, the exit status and the files it wrote (B09 WP61, F280).</summary>
/// <param name="Stdout">What the program wrote to its standard output (<c>stdout</c> as an output file name).</param>
/// <param name="Stderr">What the library printed on the error stream.</param>
/// <param name="ExitCode">0, or 1 when the program died.</param>
/// <param name="Files">The files the program wrote, by name.</param>
public sealed record EmbossRunResult(string Stdout, string Stderr, int ExitCode, IReadOnlyDictionary<string, string> Files)
{
    public static string Join(IEnumerable<string> messages)
    {
        var sb = new StringBuilder();
        foreach (string m in messages)
            sb.Append(m).Append('\n');
        return sb.ToString();
    }
}

/// <summary>
/// The EMBOSS 6.6.0 programs of this port run from their command line (B09 WP61, F280): the arguments are resolved by the ACD library (<see cref="EmbossAcdArgs"/>), the values are read in the order of the ACD file
/// (the first sequence is read before the later values are checked), a bad value stops the program with the library's messages, and the result is what the program leaves (<see cref="EmbossRunResult"/>).
/// </summary>
public static class EmbossProgramRunner
{
    private static readonly Dictionary<string, Func<EmbossAcdArgs, IReadOnlyDictionary<string, string>, EmbossRunResult>> Programs = new(StringComparer.Ordinal)
    {
        ["revseq"] = EmbossSequencePrograms.Revseq,
        ["transeq"] = EmbossSequencePrograms.Transeq,
        ["extractseq"] = EmbossSequencePrograms.Extractseq,
        ["cutseq"] = EmbossSequencePrograms.Cutseq,
        ["degapseq"] = EmbossSequencePrograms.Degapseq,
        ["trimseq"] = EmbossSequencePrograms.Trimseq,
        ["skipseq"] = EmbossSequencePrograms.Skipseq,
        ["nthseq"] = EmbossSequencePrograms.Nthseq,
        ["notseq"] = EmbossSequencePrograms.Notseq,
        ["pasteseq"] = EmbossSequencePrograms.Pasteseq,
        ["trimest"] = EmbossSequencePrograms.Trimest,
        ["newseq"] = EmbossSequencePrograms.Newseq,
        ["yank"] = EmbossSequencePrograms.Yank,
        ["seqcount"] = EmbossSequencePrograms.Seqcount,
        ["splitter"] = EmbossSequencePrograms.Splitter,
        ["maskseq"] = EmbossSequencePrograms.Maskseq,
        ["maskambignuc"] = EmbossSequencePrograms.Maskambignuc,
        ["maskambigprot"] = EmbossSequencePrograms.Maskambigprot,
        ["trimspace"] = EmbossSequencePrograms.Trimspace,
        ["nospace"] = EmbossSequencePrograms.Nospace,
        ["notab"] = EmbossSequencePrograms.Notab,
        ["noreturn"] = EmbossSequencePrograms.Noreturn,
        ["shuffleseq"] = EmbossSequencePrograms.Shuffleseq,
        ["msbar"] = EmbossSequencePrograms.Msbar,
        ["makenucseq"] = EmbossSequencePrograms.Makenucseq,
        ["makeprotseq"] = EmbossSequencePrograms.Makeprotseq,
        ["nthseqset"] = EmbossSequencePrograms.Nthseqset,
        ["listor"] = EmbossSequencePrograms.Listor,
        ["infoseq"] = EmbossSequencePrograms.Infoseq,
        ["wordcount"] = EmbossSequencePrograms.Wordcount,
        ["seqxref"] = EmbossSequencePrograms.Seqxref,
        ["showorf"] = EmbossSequencePrograms.Showorf,
        ["prettyseq"] = EmbossSequencePrograms.Prettyseq,
        ["drfinddata"] = EmbossSequencePrograms.Drfinddata,
        ["dbtell"] = EmbossSequencePrograms.Dbtell,
        ["showdb"] = EmbossSequencePrograms.Showdb,
    };

    /// <summary>The programs that can be run.</summary>
    public static IReadOnlyCollection<string> ProgramNames => Programs.Keys;

    /// <summary>Runs <paramref name="program"/> with the arguments of its command line (<c>-sequence in.dat -outseq stdout -minsize 30</c>); <paramref name="files"/> are the files the program may read, by name.</summary>
    /// <exception cref="ArgumentException">An unknown program, or a qualifier the library would refuse.</exception>
    /// <exception cref="InvalidOperationException">The input is one on which the reference program never stops or dies by a fault.</exception>
    public static EmbossRunResult Run(string program, IReadOnlyList<string> args, IReadOnlyDictionary<string, string> files) => Run(program, args, files, new EmbossRunEnvironment());

    /// <summary>
    /// Runs <paramref name="program"/> as <see cref="Run(string, IReadOnlyList{string}, IReadOnlyDictionary{string, string})"/> with the clock of <paramref name="environment"/>: the date the writers print and the seed of the
    /// random programs (<c>msbar</c>, <c>shuffleseq</c>, <c>makenucseq</c>, <c>makeprotseq</c> take the microseconds of the clock as the seed).
    /// </summary>
    public static EmbossRunResult Run(string program, IReadOnlyList<string> args, IReadOnlyDictionary<string, string> files, EmbossRunEnvironment environment) =>
        Run(Programs, program, args, files, environment);

    /// <summary>Runs <paramref name="program"/> of the table <paramref name="programs"/> (the programs of one assembly: the other assemblies of the port keep the programs of their own area in tables of their own) as <see cref="Run(string, IReadOnlyList{string}, IReadOnlyDictionary{string, string}, EmbossRunEnvironment)"/>.</summary>
    public static EmbossRunResult Run(IReadOnlyDictionary<string, Func<EmbossAcdArgs, IReadOnlyDictionary<string, string>, EmbossRunResult>> programs, string program, IReadOnlyList<string> args, IReadOnlyDictionary<string, string> files, EmbossRunEnvironment environment)
    {
        ArgumentNullException.ThrowIfNull(programs);
        ArgumentNullException.ThrowIfNull(program);
        ArgumentNullException.ThrowIfNull(environment);
        ArgumentNullException.ThrowIfNull(args);
        ArgumentNullException.ThrowIfNull(files);
        if (!programs.TryGetValue(program, out var run))
            throw new ArgumentException($"The EMBOSS program '{program}' is not ported.", nameof(program));
        using var scope = EmbossRunEnvironment.Enter(environment);
        return run(EmbossAcdArgs.Parse(program, args), files);
    }
}
