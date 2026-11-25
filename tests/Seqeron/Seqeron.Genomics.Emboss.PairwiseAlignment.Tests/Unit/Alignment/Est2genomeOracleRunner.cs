// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.PairwiseAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Emboss.PairwiseAlignment;
using System.Globalization;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>Maps the command line of a recorded <c>est2genome</c> run (TestData/Alignment/est2genome_oracle.jsonl.gz, scripts/alignment/generate_est2genome_oracle.py) onto <see cref="EmbossEst2genomeRequest"/>.</summary>
public static class Est2genomeOracleRunner
{
    private static readonly CultureInfo Ci = CultureInfo.InvariantCulture;

    /// <summary>The output of the port, or null when the program is expected to fail (crash / usage error).</summary>
    public static EmbossEst2genomeOutput? Run(IReadOnlyList<string> a, string estFasta, string genomeFasta)
    {
        var req = new EmbossEst2genomeRequest(estFasta, genomeFasta);
        var qualifiers = new List<string>();
        int i = 0;
        while (i < a.Count)
        {
            string opt = a[i++];
            switch (opt)
            {
                case "-match": req = req with { Match = int.Parse(a[i++], Ci) }; break;
                case "-mismatch": req = req with { Mismatch = int.Parse(a[i++], Ci) }; break;
                case "-gappenalty": req = req with { GapPenalty = int.Parse(a[i++], Ci) }; break;
                case "-intronpenalty": req = req with { IntronPenalty = int.Parse(a[i++], Ci) }; break;
                case "-splicepenalty": req = req with { SplicePenalty = int.Parse(a[i++], Ci) }; break;
                case "-minscore": req = req with { MinScore = int.Parse(a[i++], Ci) }; break;
                case "-reverse": req = req with { Reverse = true }; break;
                case "-nousesplice": req = req with { UseSplice = false }; break;
                case "-mode": req = req with { Mode = a[i++] }; break;
                case "-nobest": req = req with { Best = false }; break;
                case "-space": req = req with { Space = float.Parse(a[i++], Ci) }; break;
                case "-shuffle": req = req with { Shuffle = int.Parse(a[i++], Ci) }; break;
                case "-seed": req = req with { Seed = int.Parse(a[i++], Ci) }; break;
                case "-align": req = req with { Align = true }; break;
                case "-width": req = req with { Width = int.Parse(a[i++], Ci) }; break;
                case "-sbegin1" or "-send1" or "-sbegin2" or "-send2": qualifiers.Add(opt + " " + a[i++]); break;
                case "-sreverse1" or "-slower1" or "-supper1" or "-sreverse2" or "-slower2" or "-supper2": qualifiers.Add(opt); break;
                default: throw new ArgumentException(opt);
            }
        }

        if (qualifiers.Count > 0)
            req = req with { SequenceQualifiers = string.Join(' ', qualifiers) };
        try
        {
            return EmbossEst2genomePrograms.Run(req);
        }
        catch (Exception e) when (e is ArgumentException or InvalidOperationException)
        {
            return null;
        }
    }
}
