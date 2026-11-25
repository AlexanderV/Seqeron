// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.Motifs, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.Emboss.Motifs;

/// <summary>The EMBOSS 6.6.0 programs of this assembly run from their command line (see <see cref="EmbossProgramRunner"/>).</summary>
public static class EmbossMotifProgramRunner
{
    private static readonly Dictionary<string, Func<EmbossAcdArgs, IReadOnlyDictionary<string, string>, EmbossRunResult>> Programs = new(StringComparer.Ordinal)
    {
        ["marscan"] = EmbossMotifPrograms.Marscan,
        ["patmatdb"] = EmbossMotifPrograms.Patmatdb,
        ["patmatmotifs"] = EmbossMotifPrograms.Patmatmotifs,
        ["prosextract"] = EmbossMotifPrograms.Prosextract,
    };

    /// <summary>The programs of this assembly that can be run.</summary>
    public static IReadOnlyCollection<string> ProgramNames => Programs.Keys;

    /// <summary>Runs <paramref name="program"/> with the arguments of its command line; <paramref name="files"/> are the files the program may read, by name.</summary>
    /// <exception cref="ArgumentException">An unknown program, or a qualifier the library would refuse.</exception>
    public static EmbossRunResult Run(string program, IReadOnlyList<string> args, IReadOnlyDictionary<string, string> files) => Run(program, args, files, new EmbossRunEnvironment());

    /// <summary>Runs <paramref name="program"/> with the clock of <paramref name="environment"/>.</summary>
    public static EmbossRunResult Run(string program, IReadOnlyList<string> args, IReadOnlyDictionary<string, string> files, EmbossRunEnvironment environment) =>
        EmbossProgramRunner.Run(Programs, program, args, files, environment);
}
