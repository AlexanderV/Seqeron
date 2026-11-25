// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Emboss.Motifs;
using Seqeron.Genomics.Emboss.OpenReadingFrames;
using Seqeron.Genomics.Emboss.PairwiseAlignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>Runs a program of the EMBOSS port by name over the runners of every assembly (the programs of the port are spread over the assemblies by area).</summary>
internal static class EmbossAnyProgramRunner
{
    public static EmbossRunResult Run(string program, IReadOnlyList<string> args, IReadOnlyDictionary<string, string> files, EmbossRunEnvironment environment)
    {
        if (EmbossMotifProgramRunner.ProgramNames.Contains(program))
            return EmbossMotifProgramRunner.Run(program, args, files, environment);
        if (EmbossOrfProgramRunner.ProgramNames.Contains(program))
            return EmbossOrfProgramRunner.Run(program, args, files, environment);
        if (EmbossPairwiseProgramRunner.ProgramNames.Contains(program))
            return EmbossPairwiseProgramRunner.Run(program, args, files, environment);
        return EmbossProgramRunner.Run(program, args, files, environment);
    }
}
