// SPDX-License-Identifier: LGPL-3.0-or-later
// Derived from ClustalW 2.1 (http://www.clustal.org/), Copyright (C) Des Higgins, Julie Thompson, Toby Gibson; LGPL-3.0-or-later.
// Part of Seqeron.Genomics.ClustalW, licensed under the GNU Lesser General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.ClustalW;
using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-14d, B09 F297): ClustalW 2.1 <c>-CLUSTERING=UPGMA</c> together with <c>-ITERATION=alignment | tree</c> (<c>-NUMITER</c>), <c>-PROFILE1 -PROFILE2 -PROFILE</c> and
/// <c>-SEQUENCES</c> - the code paths of <c>TreeInterface::getWeightsForProfileAlignUPGMA</c> and of the iteration with UPGMA tree weights - against the program's own rows.
/// Oracle (TestData/Alignment/clustalw_upgma_iteration_oracle.tsv, scripts/alignment/generate_clustalw_upgma_iteration_oracle.py): the unmodified Debian clustalw 2.1+lgpl-7 source built at -O0 by the user
/// (docker package handoff/A2-14-upgma-iteration; the installed binary and the -O2 build segmentation-fault on every UPGMA run) on 800 sets (400 DNA, 400 protein: 250 align-alignment, 200 align-tree,
/// 175 profile, 175 sequences; -NOWEIGHTS -MAXDIV -GAPOPEN -GAPEXT -MATRIX -NUMITER -QUICKTREE); the profiles of the profile commands are alignments written by the program itself (NJ).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_ClustalWUpgmaIteration_Tests
{
    private sealed record Row(string Id, string Mode, string Type, string[] Profile1, string[] Profile2, string Options, string[] Rows);

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.clustalw_upgma_iteration_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var list = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            list.Add(new Row(f[0], f[1], f[2], f[3].Split(';'), f[4].Length == 0 ? [] : f[4].Split(';'), f[5], f[6].Split(';')));
        }

        return list;
    });

    private static double D(string s) => double.Parse(s, CultureInfo.InvariantCulture);

    private static ClustalWProfileMatrix ProfileMatrix(string name) => name switch
    {
        "GONNET" => ClustalWProfileMatrix.Gonnet,
        "BLOSUM" => ClustalWProfileMatrix.Blosum,
        "PAM" => ClustalWProfileMatrix.Pam,
        "ID" => ClustalWProfileMatrix.Identity,
        _ => throw new ArgumentException(name),
    };

    private static ClustalWIteration IterationOf(Row r) => r.Mode == "align-tree" ? ClustalWIteration.Tree : r.Mode == "align-alignment" ? ClustalWIteration.Alignment : ClustalWIteration.None;

    private static ClustalWOptions DnaOptions(Row r)
    {
        var o = new ClustalWOptions { Clustering = ClustalWClustering.Upgma, Iteration = IterationOf(r) };
        foreach (var kv in r.Options.Length == 0 ? [] : r.Options.Split(','))
        {
            var p = kv.Split('=');
            switch (p[0])
            {
                case "n1": break;
                case "numiter": o = o with { NumIterations = int.Parse(p[1]) }; break;
                case "gapopen": o = o with { GapOpen = D(p[1]) }; break;
                case "gapext": o = o with { GapExtend = D(p[1]) }; break;
                case "maxdiv": o = o with { DivergenceCutoff = int.Parse(p[1]) }; break;
                case "noweights": o = o with { NoWeights = true }; break;
                case "quicktree": o = o with { PairwiseMode = ClustalWPairwiseMode.Quick }; break;
                default: throw new InvalidOperationException($"unknown oracle option {p[0]}");
            }
        }

        return o;
    }

    private static ClustalWProteinOptions ProteinOptions(Row r)
    {
        var o = new ClustalWProteinOptions { Clustering = ClustalWClustering.Upgma, Iteration = IterationOf(r) };
        foreach (var kv in r.Options.Length == 0 ? [] : r.Options.Split(','))
        {
            var p = kv.Split('=');
            switch (p[0])
            {
                case "n1": break;
                case "numiter": o = o with { NumIterations = int.Parse(p[1]) }; break;
                case "gapopen": o = o with { GapOpen = D(p[1]) }; break;
                case "gapext": o = o with { GapExtend = D(p[1]) }; break;
                case "maxdiv": o = o with { DivergenceCutoff = int.Parse(p[1]) }; break;
                case "noweights": o = o with { NoWeights = true }; break;
                case "matrix": o = o with { Matrix = ProfileMatrix(p[1]) }; break;
                case "quicktree": o = o with { PairwiseMode = ClustalWPairwiseMode.Quick }; break;
                default: throw new InvalidOperationException($"unknown oracle option {p[0]}");
            }
        }

        return o;
    }

    private static string[] Run(Row r)
    {
        if (r.Mode.StartsWith("align", StringComparison.Ordinal))
            return (r.Type == "dna"
                ? ClustalWAligner.MultipleAlignClustalW(r.Profile1, DnaOptions(r))
                : ClustalWAligner.MultipleAlignClustalW(r.Profile1, ProteinOptions(r))).AlignedSequences.ToArray();
        if (r.Type == "dna")
        {
            var o = DnaOptions(r);
            return (r.Mode == "profile"
                ? ClustalWAligner.AlignClustalWProfiles(r.Profile1, r.Profile2, o)
                : ClustalWAligner.AlignClustalWSequencesToProfile(r.Profile1, r.Profile2, o)).AlignedSequences.ToArray();
        }
        else
        {
            var o = ProteinOptions(r);
            return (r.Mode == "profile"
                ? ClustalWAligner.AlignClustalWProfiles(r.Profile1, r.Profile2, o)
                : ClustalWAligner.AlignClustalWSequencesToProfile(r.Profile1, r.Profile2, o)).AlignedSequences.ToArray();
        }
    }

    [Test]
    public void Oracle_Has800Sets_OfEveryModeAndBothTypes()
    {
        var rows = Oracle.Value;
        Assert.Multiple(() =>
        {
            Assert.That(rows, Has.Count.EqualTo(800));
            Assert.That(rows.Count(r => r.Type == "dna"), Is.EqualTo(400));
            Assert.That(rows.Count(r => r.Mode == "align-alignment"), Is.EqualTo(250));
            Assert.That(rows.Count(r => r.Mode == "align-tree"), Is.EqualTo(200));
            Assert.That(rows.Count(r => r.Mode == "profile"), Is.EqualTo(175));
            Assert.That(rows.Count(r => r.Mode == "sequences"), Is.EqualTo(175));
        });
    }

    [TestCase("align-alignment")]
    [TestCase("align-tree")]
    [TestCase("profile")]
    [TestCase("sequences")]
    public void UpgmaWithIterationAndProfileCommands_Rows_EqualClustalW21(string mode)
    {
        int bad = 0, n = 0;
        var first = new List<string>();
        foreach (var r in Oracle.Value.Where(x => x.Mode == mode))
        {
            n++;
            string[] got = Run(r);
            if (!got.SequenceEqual(r.Rows))
            {
                bad++;
                if (first.Count < 3)
                    first.Add($"set {r.Id} ({r.Mode} {r.Type}, {r.Options}): expected {string.Join('/', r.Rows)}, got {string.Join('/', got)}");
            }
        }

        Assert.That(n, Is.GreaterThan(150));
        Assert.That(bad, Is.Zero, string.Join("\n", first));
    }

    /// <summary>The oracle is not vacuous: with the neighbour-joining guide tree instead of UPGMA the same options give other rows on a good share of the sets (so the UPGMA weights decide).</summary>
    [Test]
    public void UpgmaRows_DifferFromTheNeighbourJoiningRows_OnManySets()
    {
        int differs = 0, n = 0;
        foreach (var r in Oracle.Value.Where(x => x.Mode is "align-alignment" or "align-tree"))
        {
            n++;
            string[] nj = r.Type == "dna"
                ? ClustalWAligner.MultipleAlignClustalW(r.Profile1, DnaOptions(r) with { Clustering = ClustalWClustering.NeighborJoining }).AlignedSequences.ToArray()
                : ClustalWAligner.MultipleAlignClustalW(r.Profile1, ProteinOptions(r) with { Clustering = ClustalWClustering.NeighborJoining }).AlignedSequences.ToArray();
            if (!nj.SequenceEqual(r.Rows)) differs++;
        }

        Assert.That(n, Is.EqualTo(450));
        Assert.That(differs, Is.GreaterThan(20), "sets where UPGMA and NJ give different rows");
    }
}
