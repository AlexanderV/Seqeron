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
/// ALIGN-MULTI-001 (audit A2-14, B09 F101): the ClustalW 2.1 <c>-CLUSTERING=UPGMA -ALIGN</c> multiple alignment (<c>SequenceAligner.ClustalWUpgma.cs</c>, <c>ClustalWOptions.Clustering</c>,
/// <c>ClustalWProteinOptions.Clustering</c>) against the program's own rows.
/// Oracle (TestData/Alignment/clustalw_upgma_oracle.tsv, scripts/alignment/generate_clustalw_upgma_oracle.py): the unmodified Debian clustalw 2.1+lgpl-7 source built at -O0 by the user (docker
/// package handoff/A2-14-clustalw-upgma; the installed binary and the -O2 build of the same source segmentation-fault on every UPGMA run: undefined behaviour of the program) on 600 sets (300 DNA,
/// 300 protein, 2-9 sequences; -NOWEIGHTS -MAXDIV -GAPOPEN -GAPEXT -MATRIX -PWMATRIX -QUICKTREE).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_ClustalWUpgmaOracle_Tests
{
    private sealed record Row(string Id, string Type, string[] Seqs, string Options, string[] Rows, string Dnd);

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.clustalw_upgma_oracle.tsv";
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
            list.Add(new Row(f[0], f[1], f[2].Split(';'), f[3], f[4].Split(';'), f[5]));
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

    private static ClustalWPairwiseMatrix PairwiseMatrix(string name) => name switch
    {
        "GONNET" => ClustalWPairwiseMatrix.Gonnet,
        "BLOSUM" => ClustalWPairwiseMatrix.Blosum,
        "PAM" => ClustalWPairwiseMatrix.Pam,
        "ID" => ClustalWPairwiseMatrix.Identity,
        _ => throw new ArgumentException(name),
    };

    private static ClustalWOptions DnaOptions(Row r)
    {
        var o = new ClustalWOptions { Clustering = ClustalWClustering.Upgma };
        foreach (var kv in r.Options.Length == 0 ? [] : r.Options.Split(','))
        {
            var p = kv.Split('=');
            switch (p[0])
            {
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
        var o = new ClustalWProteinOptions { Clustering = ClustalWClustering.Upgma };
        foreach (var kv in r.Options.Length == 0 ? [] : r.Options.Split(','))
        {
            var p = kv.Split('=');
            switch (p[0])
            {
                case "gapopen": o = o with { GapOpen = D(p[1]) }; break;
                case "gapext": o = o with { GapExtend = D(p[1]) }; break;
                case "maxdiv": o = o with { DivergenceCutoff = int.Parse(p[1]) }; break;
                case "noweights": o = o with { NoWeights = true }; break;
                case "matrix": o = o with { Matrix = ProfileMatrix(p[1]) }; break;
                case "pwmatrix": o = o with { Pairwise = o.Pairwise with { Matrix = PairwiseMatrix(p[1]) } }; break;
                case "quicktree": o = o with { PairwiseMode = ClustalWPairwiseMode.Quick }; break;
                default: throw new InvalidOperationException($"unknown oracle option {p[0]}");
            }
        }

        return o;
    }

    private static string[] Align(Row r) => r.Type == "dna"
        ? ClustalWAligner.MultipleAlignClustalW(r.Seqs, DnaOptions(r)).AlignedSequences.ToArray()
        : ClustalWAligner.MultipleAlignClustalW(r.Seqs, ProteinOptions(r)).AlignedSequences.ToArray();

    [Test]
    public void Oracle_Has600Sets_OfBothTypes()
    {
        var rows = Oracle.Value;
        Assert.Multiple(() =>
        {
            Assert.That(rows, Has.Count.EqualTo(600));
            Assert.That(rows.Count(r => r.Type == "dna"), Is.EqualTo(300));
            Assert.That(rows.Count(r => r.Type == "protein"), Is.EqualTo(300));
        });
    }

    [Test]
    public void MultipleAlignClustalW_Upgma_Rows_EqualClustalW21_OnAllOracleSets()
    {
        int bad = 0;
        var first = new List<string>();
        foreach (var r in Oracle.Value)
        {
            string[] got = Align(r);
            if (!got.SequenceEqual(r.Rows))
            {
                bad++;
                if (first.Count < 3)
                    first.Add($"set {r.Id} ({r.Type}, {r.Options}): expected {string.Join('/', r.Rows)}, got {string.Join('/', got)}");
            }
        }

        Assert.That(bad, Is.Zero, string.Join("\n", first));
    }

    [Test]
    public void GuideTreeText_EqualsClustalW21_NewtreeFile_OnAllOracleSets()
    {
        // the .dnd the program writes (RootedTreeOutput::printPhylipTree; sequences s000, s001, ...): topology, join order and the join heights (distance / 2) of the UPGMA
        int bad = 0;
        var first = new List<string>();
        foreach (var r in Oracle.Value)
        {
            var names = Enumerable.Range(0, r.Seqs.Length).Select(i => $"s{i:000}").ToArray();
            string got = r.Type == "dna"
                ? ClustalWAligner.ClustalWUpgmaGuideTreeText(r.Seqs, DnaOptions(r), names)
                : ClustalWAligner.ClustalWUpgmaGuideTreeText(r.Seqs, ProteinOptions(r), names);
            if (got.Replace('\n', '|') != r.Dnd)
            {
                bad++;
                if (first.Count < 3)
                    first.Add($"set {r.Id}: expected {r.Dnd}, got {got.Replace('\n', '|')}");
            }
        }

        Assert.That(bad, Is.Zero, string.Join("\n", first));
    }
}
