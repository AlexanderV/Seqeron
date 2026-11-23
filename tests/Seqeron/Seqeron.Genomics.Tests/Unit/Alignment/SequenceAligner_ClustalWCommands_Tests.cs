using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-14b, B09 F104): ClustalW 2.1 guide-tree files and the row order. <see cref="SequenceAligner.RunClustalW(IReadOnlyList{ClustalWSequence}, ClustalWOptions?, ClustalWRunOptions?)"/>,
/// <see cref="SequenceAligner.RunClustalWProfile(IReadOnlyList{ClustalWSequence}, IReadOnlyList{ClustalWSequence}, ClustalWOptions?, ClustalWRunOptions?)"/> and
/// <see cref="SequenceAligner.RunClustalWSequences(IReadOnlyList{ClustalWSequence}, IReadOnlyList{ClustalWSequence}, ClustalWOptions?, ClustalWRunOptions?)"/> against the installed <c>clustalw</c> 2.1:
/// the guide-tree file of <c>-ALIGN -NEWTREE</c> / <c>-NEWTREE</c> (<c>printPhylipTree</c>), the alignment along a user tree (<c>-USETREE</c>, <c>Tree::readTree</c>), the profile commands with
/// <c>-USETREE1/2</c> and the written trees of both profiles, and the order of the rows (<c>-OUTORDER</c>). Oracle: TestData/Alignment/clustalw_commands_oracle.tsv
/// (scripts/alignment/generate_clustalw_commands_oracle.py), 3 300 sets.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_ClustalWCommands_Tests
{
    private sealed record Row(string Id, string Command, string Type, string[] Names, string[] Seqs, int P1, string Options, string OutOrder, string Tree1, string Tree2,
        string[] OutNames, string[] OutRows, string TreeOut1, string TreeOut2, string Error);

    private static string Unescape(string s)
    {
        var sb = new System.Text.StringBuilder();
        bool escaped = false;
        foreach (char c in s)
        {
            if (escaped)
            {
                sb.Append(c switch { 'n' => '\n', 't' => '\t', _ => c });
                escaped = false;
            }
            else if (c == '\\')
            {
                escaped = true;
            }
            else
            {
                sb.Append(c);
            }
        }

        return sb.ToString();
    }

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.clustalw_commands_oracle.tsv";
        using var stream = TestDataResource.Open(name)
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
            list.Add(new Row(f[0], f[1], f[2], f[3].Split(';'), f[4].Split(';'), int.Parse(f[5], CultureInfo.InvariantCulture), f[6], f[7], Unescape(f[8]), Unescape(f[9]),
                f[10].Length == 0 ? [] : f[10].Split(';'), f[11].Length == 0 ? [] : f[11].Split(';'), Unescape(f[12]), Unescape(f[13]), f[14]));
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

    private static ClustalWRunResult Run(Row r)
    {
        bool treesOnly = r.Command.EndsWith("-trees", StringComparison.Ordinal);
        string command = treesOnly ? r.Command[..^6] : r.Command;
        var all = r.Names.Zip(r.Seqs, (n, s) => new ClustalWSequence(n, s)).ToList();
        var run = new ClustalWRunOptions { OutputOrder = r.OutOrder == "input" ? ClustalWOutputOrder.Input : ClustalWOutputOrder.Aligned };
        run = command switch
        {
            "newtree" => run with { GuideTreeOnly = true },
            "usetree" => run with { UseTree = r.Tree1 },
            "profile" => run with { UseTree1 = r.Tree1.Length > 0 ? r.Tree1 : null, UseTree2 = r.Tree2.Length > 0 ? r.Tree2 : null, GuideTreeOnly = treesOnly },
            "sequences" => run with { UseTree = r.Tree1.Length > 0 ? r.Tree1 : null, GuideTreeOnly = treesOnly },
            _ => run,
        };

        var opts = r.Options.Length == 0 ? [] : r.Options.Split(',');
        if (r.Type == "dna")
        {
            var o = new ClustalWOptions();
            foreach (var kv in opts)
            {
                var p = kv.Split('=');
                switch (p[0])
                {
                    case "gapopen": o = o with { GapOpen = D(p[1]) }; break;
                    case "gapext": o = o with { GapExtend = D(p[1]) }; break;
                    case "maxdiv": o = o with { DivergenceCutoff = int.Parse(p[1], CultureInfo.InvariantCulture) }; break;
                    case "noweights": o = o with { NoWeights = true }; break;
                    case "negative": o = o with { NegativeMatrix = true }; break;
                    case "quick": o = o with { PairwiseMode = ClustalWPairwiseMode.Quick }; break;
                    case "transweight": o = o with { TransitionWeight = float.Parse(p[1], CultureInfo.InvariantCulture) }; break;
                    case "pwgapopen": o = o with { PairwiseGapOpen = D(p[1]) }; break;
                    case "pwgapext": o = o with { PairwiseGapExtend = D(p[1]) }; break;
                    default: throw new InvalidOperationException($"unknown oracle option {p[0]}");
                }
            }

            return command switch
            {
                "align" or "newtree" or "usetree" => SequenceAligner.RunClustalW(all, o, run),
                "profile" => SequenceAligner.RunClustalWProfile(all.Take(r.P1).ToList(), all.Skip(r.P1).ToList(), o, run),
                _ => SequenceAligner.RunClustalWSequences(all.Take(r.P1).ToList(), all.Skip(r.P1).ToList(), o, run),
            };
        }
        else
        {
            var o = new ClustalWProteinOptions();
            foreach (var kv in opts)
            {
                var p = kv.Split('=');
                switch (p[0])
                {
                    case "gapopen": o = o with { GapOpen = D(p[1]) }; break;
                    case "gapext": o = o with { GapExtend = D(p[1]) }; break;
                    case "maxdiv": o = o with { DivergenceCutoff = int.Parse(p[1], CultureInfo.InvariantCulture) }; break;
                    case "noweights": o = o with { NoWeights = true }; break;
                    case "negative": o = o with { NegativeMatrix = true }; break;
                    case "quick": o = o with { PairwiseMode = ClustalWPairwiseMode.Quick }; break;
                    case "matrix": o = o with { Matrix = ProfileMatrix(p[1]) }; break;
                    case "pwmatrix": o = o with { Pairwise = o.Pairwise with { Matrix = PairwiseMatrix(p[1]) } }; break;
                    case "pwgapopen": o = o with { Pairwise = o.Pairwise with { GapOpen = D(p[1]) } }; break;
                    case "pwgapext": o = o with { Pairwise = o.Pairwise with { GapExtend = D(p[1]) } }; break;
                    case "nopgap": o = o with { NoPreferredGapPenalties = true }; break;
                    case "nohgap": o = o with { NoHydrophilicGapPenalties = true }; break;
                    case "gapdist": o = o with { GapSeparationDistance = int.Parse(p[1], CultureInfo.InvariantCulture) }; break;
                    default: throw new InvalidOperationException($"unknown oracle option {p[0]}");
                }
            }

            return command switch
            {
                "align" or "newtree" or "usetree" => SequenceAligner.RunClustalW(all, o, run),
                "profile" => SequenceAligner.RunClustalWProfile(all.Take(r.P1).ToList(), all.Skip(r.P1).ToList(), o, run),
                _ => SequenceAligner.RunClustalWSequences(all.Take(r.P1).ToList(), all.Skip(r.P1).ToList(), o, run),
            };
        }
    }

    /// <summary>Compares one oracle set; returns null when equal, else the first difference.</summary>
    private static string? Compare(Row r)
    {
        bool treesOnly = r.Command.EndsWith("-trees", StringComparison.Ordinal);
        string command = treesOnly ? r.Command[..^6] : r.Command;
        ClustalWRunResult res;
        try
        {
            res = Run(r);
        }
        catch (ClustalWTreeFileException e)
        {
            // Undefined: the tree puts the program into undefined behaviour (it hangs, corrupts its heap or frees a node twice); the oracle shows a crash for most of them and
            // an ordinary outcome for a few, so this library refuses all of them. Hang / Crash rows: any refusal is accepted; the others: the very refusal.
            if (e.Problem == ClustalWTreeFileProblem.Undefined)
            {
                return null;
            }

            if (r.Error is "Hang" or "Crash")
                return null;
            return r.Error == e.Problem.ToString() ? null : $"expected {(r.Error.Length == 0 ? "a result" : r.Error)}, got {e.Problem}";
        }

        if (r.Error.Length > 0)
        {
            // the program's crashes on trees with single-child nodes depend on the heap; a result is not claimed to differ
            return r.Error is "Hang" or "Crash" ? null : $"expected {r.Error}, got a result";
        }

        if (command == "newtree")
            return res.GuideTree == r.TreeOut1 ? null : $"tree {res.GuideTree} vs {r.TreeOut1}";
        if (command == "profile")
        {
            int p2 = r.Names.Length - r.P1;
            if (r.Tree1.Length == 0 && r.P1 >= 2 && res.GuideTree1 != r.TreeOut1)
                return $"tree1 {res.GuideTree1} vs {r.TreeOut1}";
            if (r.Tree2.Length == 0 && p2 >= 2 && res.GuideTree2 != r.TreeOut2)
                return $"tree2 {res.GuideTree2} vs {r.TreeOut2}";
        }
        else if ((command == "align" || (command == "sequences" && r.Tree1.Length == 0)) && res.GuideTree != r.TreeOut1)
        {
            return $"tree {res.GuideTree} vs {r.TreeOut1}";
        }

        if (treesOnly)
            return res.Alignment is null ? null : "an alignment was produced for a trees-only run";
        if (res.Alignment is null)
            return "no alignment";
        if (!res.Alignment.GetOutputNames().SequenceEqual(r.OutNames))
            return $"order {string.Join(';', res.Alignment.GetOutputNames())} vs {string.Join(';', r.OutNames)}";
        return res.Alignment.GetOutputRows().SequenceEqual(r.OutRows) ? null : "rows differ";
    }

    private static void AssertAll(Func<Row, bool> filter, int minimum)
    {
        var rows = Oracle.Value.Where(filter).ToList();
        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(minimum));
        var failures = new List<string>();
        foreach (var r in rows)
        {
            string? diff = Compare(r);
            if (diff != null)
                failures.Add($"set {r.Id} ({r.Command}, {r.Type}, {r.Options}): {diff}");
        }

        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(3)) + $"\n{failures.Count} / {rows.Count} sets differ");
    }

    /// <summary>
    /// The oracle rows on which the PROGRAM is in undefined behaviour, by set id and by what the port does (every other row is compared exactly by the five <c>AssertAll</c> tests). The program hangs, corrupts its
    /// heap or frees a node twice on a user tree that has single-child nodes or names the same sequence twice (<c>-USETREE</c>, ClustalW 2.1 <c>src/tree/</c>), and the oracle shows a crash for most of them
    /// and an ordinary outcome for a few; the port refuses every such tree (<see cref="ClustalWTreeFileProblem.Undefined"/>) and accepts any refusal on a Hang / Crash row. Format: <c>id: what the port does [what the oracle recorded]</c>;
    /// the list is the exact set (generated by running every row, see the test), no count tolerance.
    /// </summary>
    private static readonly string[] UndefinedBehaviourRows =
    [
        "1429/usetree: refused NotCompatible [Crash]",
        "1488/usetree: refused Undefined [Hang]",
        "1512/usetree: refused Undefined [Crash]",
        "1541/usetree: refused Undefined [Hang]",
        "1567/usetree: refused Undefined [Hang]",
        "1704/usetree: refused NotCompatible [Crash]",
        "1707/usetree: refused Undefined [NameNotFound]",
        "1722/usetree: refused Undefined [Hang]",
        "1760/usetree: refused Undefined [Crash]",
        "1779/usetree: refused NotCompatible [Crash]",
        "1848/usetree: result [Crash]",
        "1893/usetree: refused Undefined [Hang]",
        "1950/usetree: result [Crash]",
        "2153/profile: refused NotCompatible [Crash]",
        "2192/profile: result [Crash]",
        "2365/profile: refused Undefined [NameNotFound]",
        "2436/profile: refused Undefined [Hang]",
        "2477/profile: refused Undefined [Hang]",
        "2499/profile: refused Undefined [Hang]",
        "2863/sequences: refused NotCompatible [Crash]",
        "2951/sequences: refused NotCompatible [Crash]",
        "2968/sequences: refused Undefined [result]",
        "3158/sequences: refused Undefined [Hang]",
        "3160/sequences: refused NotCompatible [Crash]",
    ];

    [Test]
    public void UndefinedBehaviourRows_AreExactlyTheEnumeratedOnes()
    {
        var observed = new List<string>();
        foreach (var r in Oracle.Value)
        {
            string oracle = r.Error.Length == 0 ? "result" : r.Error;
            try
            {
                Run(r);
                if (r.Error is "Hang" or "Crash")
                    observed.Add($"{r.Id}/{r.Command}: result [{oracle}]");
            }
            catch (ClustalWTreeFileException e)
            {
                if (e.Problem == ClustalWTreeFileProblem.Undefined)
                    observed.Add($"{r.Id}/{r.Command}: refused Undefined [{oracle}]");
                else if (r.Error is "Hang" or "Crash")
                    observed.Add($"{r.Id}/{r.Command}: refused {e.Problem} [{oracle}]");
            }
        }

        Assert.That(observed, Is.EqualTo(UndefinedBehaviourRows), "observed:\n" + string.Join("\n", observed.Select(x => "\"" + x + "\",")));
    }

    [Test]
    public void Oracle_Has300SetsOfEachFeatureAndType()
    {
        var rows = Oracle.Value;
        Assert.Multiple(() =>
        {
            foreach (string type in new[] { "dna", "protein" })
            {
                Assert.That(rows.Count(r => r.Type == type && r.Command == "align"), Is.GreaterThanOrEqualTo(300), $"align {type}");
                Assert.That(rows.Count(r => r.Type == type && r.Command == "newtree"), Is.GreaterThanOrEqualTo(300), $"newtree {type}");
                Assert.That(rows.Count(r => r.Type == type && r.Command == "usetree"), Is.GreaterThanOrEqualTo(300), $"usetree {type}");
                Assert.That(rows.Count(r => r.Type == type && r.Command.StartsWith("profile", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(300), $"profile {type}");
                Assert.That(rows.Count(r => r.Type == type && r.Command.StartsWith("sequences", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(300), $"sequences {type}");
            }

            Assert.That(rows.Count(r => r.Error.Length > 0 && r.Command == "usetree"), Is.GreaterThanOrEqualTo(100), "refusals of user trees");
            foreach (string kind in new[] { "WrongFormat", "NotCompatible", "NameNotFound", "UnrootedWithoutDistances" })
                Assert.That(rows.Any(r => r.Error == kind), Is.True, kind);
            Assert.That(rows.Count(r => r.OutOrder == "input"), Is.GreaterThan(300));
            Assert.That(rows.Count(r => r.OutOrder == "aligned"), Is.GreaterThan(300));
        });
    }

    [Test]
    public void AlignWithNewTree_RowsOrderAndGuideTree_EqualClustalW21_OnAllOracleSets() => AssertAll(r => r.Command == "align", 660);

    [Test]
    public void NewTreeOnly_GuideTree_EqualClustalW21_OnAllOracleSets() => AssertAll(r => r.Command == "newtree", 660);

    [Test]
    public void UseTree_RowsOrderAndRefusals_EqualClustalW21_OnAllOracleSets() => AssertAll(r => r.Command == "usetree", 660);

    [Test]
    public void ProfileCommand_RowsAndProfileTrees_EqualClustalW21_OnAllOracleSets() => AssertAll(r => r.Command.StartsWith("profile", StringComparison.Ordinal), 660);

    [Test]
    public void SequencesCommand_RowsAndGuideTree_EqualClustalW21_OnAllOracleSets() => AssertAll(r => r.Command.StartsWith("sequences", StringComparison.Ordinal), 660);

    [Test]
    public void FiveDnaSequences_GuideTreeFileAndAlignedOrder_AreTheProgramsOutput()
    {
        // clustalw -INFILE=a.fa -TYPE=DNA -ALIGN -OUTPUT=CLUSTAL -NEWTREE=t.dnd (installed 2.1+lgpl-7): t.dnd and the row order of the alignment file
        var seqs = new[]
        {
            new ClustalWSequence("s1", "ACGTACGTTAGCTAGCTAGGATCGATCGAT"),
            new ClustalWSequence("s2", "ACGTACGATAGCTAGCTAGGATCGTTCGAT"),
            new ClustalWSequence("s3", "ACGAACGTTAGCTTGCTAGGATCGATCAAT"),
            new ClustalWSequence("s4", "TCGTACGTTAGCTAGCAAGGATCGATCGAT"),
            new ClustalWSequence("s5", "ACGTACGTTAGGTAGCTAGGTTCGATCGAA"),
        };
        var r = SequenceAligner.RunClustalW(seqs);
        Assert.Multiple(() =>
        {
            Assert.That(r.GuideTree, Is.EqualTo("(\n(\n(\ns1:0.00000,\ns2:0.06667)\n:0.00000,\ns3:0.10000)\n:0.00000,\ns4:0.06667,\ns5:0.10000);\n"));
            Assert.That(r.Alignment!.GetOutputNames(), Is.EqualTo(new[] { "s4", "s5", "s3", "s1", "s2" }));
            Assert.That(r.Alignment.GetOutputRows()[0], Is.EqualTo("TCGTACGTTAGCTAGCAAGGATCGATCGAT"));
            Assert.That(r.Alignment!.Names, Is.EqualTo(new[] { "s1", "s2", "s3", "s4", "s5" }), "rows keep the input order");
        });
        var input = SequenceAligner.RunClustalW(seqs, null, new ClustalWRunOptions { OutputOrder = ClustalWOutputOrder.Input });
        Assert.That(input.Alignment!.GetOutputNames(), Is.EqualTo(new[] { "s1", "s2", "s3", "s4", "s5" }));
        Assert.That(input.Alignment.Rows, Is.EqualTo(r.Alignment!.Rows));
    }

    [Test]
    public void TwoSequences_GuideTree_IsOneLineWithHalfTheDistance()
    {
        var r = SequenceAligner.RunClustalW([new ClustalWSequence("a", "ACGTACGTAC"), new ClustalWSequence("b", "ACGTACGAAC")]);
        Assert.That(r.GuideTree, Is.EqualTo("(a:0.05,b:0.05);\n"));
    }

    [Test]
    public void UserTree_NamesAreMatchedCaseInsensitively_AndRootedTreeNeedsNoDistances()
    {
        var seqs = new[] { new ClustalWSequence("Alpha", "ACGTACGTTAGCTAGCTAGG"), new ClustalWSequence("Beta", "ACGTACGATAGCTAGCTAGG"), new ClustalWSequence("Gamma", "ACGAACGTTAGCTTGCTAGG") };
        var withDistances = SequenceAligner.RunClustalW(seqs, null, new ClustalWRunOptions { UseTree = "(alpha:0.1,(BETA:0.05,gamma:0.05):0.02);" });
        var without = SequenceAligner.RunClustalW(seqs, null, new ClustalWRunOptions { UseTree = "(alpha,(beta,gamma));" });
        Assert.Multiple(() =>
        {
            Assert.That(withDistances.Alignment!.Rows, Has.Count.EqualTo(3));
            Assert.That(without.Alignment!.Rows.Select(x => x.Replace("-", "")), Is.EqualTo(seqs.Select(s => s.Residues)));
            Assert.That(withDistances.GuideTree, Is.Null, "a user tree is not rewritten");
        });
    }

    [TestCase("", ClustalWTreeFileProblem.WrongFormat)]
    [TestCase("a:0.1,b:0.1,c:0.1;", ClustalWTreeFileProblem.WrongFormat)]
    [TestCase("(a:0.1,b:0.1);", ClustalWTreeFileProblem.NotCompatible)]
    [TestCase("(a:0.1,(b:0.1,c:0.1):0.1,d:0.1);", ClustalWTreeFileProblem.NotCompatible)]
    [TestCase("(a:0.1,(b:0.1,x:0.1):0.1);", ClustalWTreeFileProblem.NameNotFound)]
    [TestCase("(a,b,c);", ClustalWTreeFileProblem.UnrootedWithoutDistances)]
    [TestCase("(a:0.1,b:0.1,c", ClustalWTreeFileProblem.Undefined)]
    public void UserTree_TheProgramsRefusals_AreTypedExceptions(string tree, ClustalWTreeFileProblem expected)
    {
        var seqs = new[] { new ClustalWSequence("a", "ACGTACGTTAGCTAGCTAGG"), new ClustalWSequence("b", "ACGTACGATAGCTAGCTAGG"), new ClustalWSequence("c", "ACGAACGTTAGCTTGCTAGG") };
        var ex = Assert.Throws<ClustalWTreeFileException>(() => SequenceAligner.RunClustalW(seqs, null, new ClustalWRunOptions { UseTree = tree }));
        Assert.That(ex!.Problem, Is.EqualTo(expected));
    }

    [Test]
    public void InvalidSequences_AreRefused()
    {
        Assert.Throws<ArgumentException>(() => SequenceAligner.RunClustalW([new ClustalWSequence("a", "ACGT")]));
        Assert.Throws<ArgumentException>(() => SequenceAligner.RunClustalW([new ClustalWSequence("a", "ACGT"), new ClustalWSequence("a", "ACGA")]));
        Assert.Throws<ArgumentException>(() => SequenceAligner.RunClustalW([new ClustalWSequence("a b", "ACGT"), new ClustalWSequence("c", "ACGA")]));
        Assert.Throws<ArgumentException>(() => SequenceAligner.RunClustalW([new ClustalWSequence("a", "ACGT"), new ClustalWSequence("", "ACGA")]));
        Assert.Throws<ArgumentException>(() => SequenceAligner.RunClustalW([new ClustalWSequence("a", "ACGT"), new ClustalWSequence("b", "")]));
        Assert.Throws<ArgumentException>(() => SequenceAligner.RunClustalWProfile([], [new ClustalWSequence("a", "ACGT")]));
    }

    [Test]
    public void ClustalWFixed_RoundsTheExactBinaryValueHalfToEven_LikePrintf()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceAligner.ClustalWFixed(0.015625, 5), Is.EqualTo("0.01562"), "exact tie, even digit stays");
            Assert.That(SequenceAligner.ClustalWFixed(0.046875, 5), Is.EqualTo("0.04688"), "exact tie, odd digit rounds up");
            Assert.That(SequenceAligner.ClustalWFixed(0.06667, 5), Is.EqualTo("0.06667"));
            Assert.That(SequenceAligner.ClustalWFixed(-0.000001, 5), Is.EqualTo("-0.00000"));
            Assert.That(SequenceAligner.ClustalWFixed(9.999996, 5), Is.EqualTo("10.00000"));
            Assert.That(SequenceAligner.ClustalWFixed(0.0, 5), Is.EqualTo("0.00000"));
        });
    }
}
