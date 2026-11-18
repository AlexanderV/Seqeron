using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A1-10b, B09 F56): <see cref="SequenceAligner.MultipleAlignClustalW(IEnumerable{string}, ClustalWOptions?, ScoringMatrix?)"/>,
/// the step-by-step port of ClustalW 2.1's multiple-alignment stage for DNA (neighbour-joining guide tree, tree-derived
/// sequence weights, delayed divergent sequences, position-specific gap penalties, Myers-Miller profile alignment).
///
/// Oracle (TestData/Alignment/clustalw_msa_full_oracle.tsv, scripts/alignment/generate_clustalw_msa_full_oracle.py): the
/// installed <c>clustalw</c> 2.1 <c>-TYPE=DNA -ALIGN -OUTPUT=FASTA -OUTORDER=INPUT</c> with the option variants -GAPOPEN, -GAPEXT,
/// -TRANSWEIGHT, -MAXDIV, -NOWEIGHTS, -NEGATIVE, -QUICKTREE (-KTUPLE, -WINDOW, -TOPDIAGS, -PAIRGAP), -PWGAPOPEN, -PWGAPEXT on
/// 1 000 DNA sets (2-25 sequences, up to 600 nt; mutated families, low-complexity, duplicates, tiny and IUB-ambiguity sets).
/// The binary computes its own pairwise stage, guide tree, weights and profile alignments: every output row must be identical.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_ClustalWMultiple_Tests
{
    private sealed record Row(string Id, string[] Seqs, string[] Rows, string Options);

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        string name = "Seqeron.Genomics.Tests.TestData.Alignment.clustalw_msa_full_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var list = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var f = line.Split('\t');
            list.Add(new Row(f[0], f[1].Split(';'), f[2].Split(';'), f.Length > 3 ? f[3] : ""));
        }

        return list;
    });

    /// <summary>The clustalw command line semantics of the oracle's option string (-KTUPLE sets the pair gap to ktuple + 4, -PAIRGAP overrides).</summary>
    private static ClustalWOptions Parse(string s)
    {
        var o = new ClustalWOptions();
        if (s.Length == 0)
            return o;
        foreach (var kv in s.Split(','))
        {
            var p = kv.Split('=');
            switch (p[0])
            {
                case "gapopen": o = o with { GapOpen = double.Parse(p[1], System.Globalization.CultureInfo.InvariantCulture) }; break;
                case "gapext": o = o with { GapExtend = double.Parse(p[1], System.Globalization.CultureInfo.InvariantCulture) }; break;
                case "pwgapopen": o = o with { PairwiseGapOpen = double.Parse(p[1], System.Globalization.CultureInfo.InvariantCulture) }; break;
                case "pwgapext": o = o with { PairwiseGapExtend = double.Parse(p[1], System.Globalization.CultureInfo.InvariantCulture) }; break;
                case "transweight": o = o with { TransitionWeight = float.Parse(p[1], System.Globalization.CultureInfo.InvariantCulture) }; break;
                case "maxdiv": o = o with { DivergenceCutoff = int.Parse(p[1]) }; break;
                case "noweights": o = o with { NoWeights = true }; break;
                case "negative": o = o with { NegativeMatrix = true }; break;
                case "quick": o = o with { PairwiseMode = ClustalWPairwiseMode.Quick }; break;
                case "ktuple": o = o with { Quick = o.Quick with { KTuple = int.Parse(p[1]), PairGap = int.Parse(p[1]) + 4 } }; break;
                case "window": o = o with { Quick = o.Quick with { Window = int.Parse(p[1]) } }; break;
                case "topdiags": o = o with { Quick = o.Quick with { TopDiagonals = int.Parse(p[1]) } }; break;
                case "pairgap": o = o with { Quick = o.Quick with { PairGap = int.Parse(p[1]) } }; break;
                default: throw new InvalidOperationException($"unknown oracle option {p[0]}");
            }
        }

        return o;
    }

    private static string[] Align(string[] seqs, ClustalWOptions? o = null) => SequenceAligner.MultipleAlignClustalW(seqs, o).AlignedSequences.ToArray();

    [Test]
    public void MultipleAlignClustalW_Rows_EqualClustalW21_OnAllOracleSets()
    {
        var rows = Oracle.Value;
        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(1000));
        int bad = 0;
        string? first = null;
        foreach (var r in rows)
        {
            string[] got = Align(r.Seqs, Parse(r.Options));
            if (!got.SequenceEqual(r.Rows) && ++bad == 1)
                first = $"set {r.Id} ({r.Options}): expected {string.Join('/', r.Rows)}, got {string.Join('/', got)}";
        }

        Assert.That(bad, Is.Zero, first);
    }

    [Test]
    public void MultipleAlignClustalW_OracleCoversEveryOptionAndTheSpecialPaths()
    {
        var rows = Oracle.Value;
        string Has(string key) => key;
        Assert.Multiple(() =>
        {
            foreach (string key in new[] { "gapopen", "gapext", "transweight", "maxdiv", "noweights", "negative", "quick", "ktuple", "window", "topdiags", "pairgap", "pwgapopen", "pwgapext" })
                Assert.That(rows.Count(r => r.Options.Split(',').Any(o => o.Split('=')[0] == Has(key))), Is.GreaterThan(20), key);
            Assert.That(rows.Count(r => r.Options.Length == 0), Is.GreaterThan(200), "default runs");
            Assert.That(rows.Count(r => r.Seqs.Length == 2), Is.GreaterThan(30), "two sequences (rooted two-leaf tree)");
            Assert.That(rows.Count(r => r.Seqs.Length >= 14), Is.GreaterThan(30), "large sets");
            Assert.That(rows.Count(r => r.Seqs.Any(s => s.Any(c => "ACGT".IndexOf(c) < 0))), Is.GreaterThan(60), "IUB ambiguity codes");
            Assert.That(rows.Count(r => r.Seqs.Distinct().Count() < r.Seqs.Length), Is.GreaterThan(20), "duplicate sequences (tied distances)");
        });
    }

    /// <summary>
    /// The oracle discriminates every mechanism: for each switch there are many sets whose rows change when the mechanism is turned
    /// off, so equality with clustalw is not an accident of easy inputs.
    /// </summary>
    [Test]
    public void MultipleAlignClustalW_EachMechanism_ChangesRowsOnManyOracleSets()
    {
        var sets = Oracle.Value.Where(r => r.Options.Length == 0 && r.Seqs.Length >= 3).Take(150).ToList();
        int noWeights = 0, delayed = 0, trans = 0, negative = 0, gap = 0, quick = 0;
        foreach (var r in sets)
        {
            var baseRows = Align(r.Seqs);
            if (!Align(r.Seqs, new ClustalWOptions { NoWeights = true }).SequenceEqual(baseRows)) noWeights++;
            if (!Align(r.Seqs, new ClustalWOptions { DivergenceCutoff = 0 }).SequenceEqual(baseRows)) delayed++;
            if (!Align(r.Seqs, new ClustalWOptions { TransitionWeight = 1f }).SequenceEqual(baseRows)) trans++;
            if (!Align(r.Seqs, new ClustalWOptions { NegativeMatrix = true }).SequenceEqual(baseRows)) negative++;
            if (!Align(r.Seqs, new ClustalWOptions { GapOpen = 5, GapExtend = 1 }).SequenceEqual(baseRows)) gap++;
            if (!Align(r.Seqs, new ClustalWOptions { PairwiseMode = ClustalWPairwiseMode.Quick }).SequenceEqual(baseRows)) quick++;
        }

        Assert.Multiple(() =>
        {
            Assert.That(noWeights, Is.GreaterThan(5), "sequence weights");
            Assert.That(delayed, Is.GreaterThan(5), "delayed divergent sequences");
            Assert.That(trans, Is.GreaterThan(3), "transition weight");
            Assert.That(negative, Is.GreaterThan(5), "negative matrix");
            Assert.That(gap, Is.GreaterThan(20), "gap penalties");
            Assert.That(quick, Is.GreaterThan(3), "fast pairwise stage");
        });
    }

    /// <summary>Hand run through the binary: clustalw -TYPE=DNA -ALIGN -OUTPUT=FASTA on these three sequences.</summary>
    [Test]
    public void MultipleAlignClustalW_ThreeSequenceExample_EqualsHandRun()
    {
        string[] seqs = ["ACGTACGTTTGACCAGTAGGCATCGAT", "ACGTACGATTGACCAGTAGCCATCGAT", "GGGGGGAAAATTTTCCCC"];
        var expected = Oracle.Value.First(r => r.Seqs.Length == 3 && r.Options.Length == 0);
        Assert.That(Align(expected.Seqs), Is.EqualTo(expected.Rows));
        var rows = Align(seqs);
        Assert.Multiple(() =>
        {
            Assert.That(rows.Select(r => r.Replace("-", "")), Is.EqualTo(seqs));
            Assert.That(rows.Select(r => r.Length).Distinct().Count(), Is.EqualTo(1));
        });
    }

    [Test]
    public void MultipleAlignClustalW_Invariants_RowsSpellInputs_EqualLength_InputOrder()
    {
        foreach (var r in Oracle.Value.Take(300))
        {
            var rows = Align(r.Seqs, Parse(r.Options));
            Assert.That(rows.Length, Is.EqualTo(r.Seqs.Length));
            Assert.That(rows.Select(x => x.Length).Distinct().Count(), Is.EqualTo(1));
            for (int i = 0; i < rows.Length; i++)
                Assert.That(rows[i].Replace("-", ""), Is.EqualTo(r.Seqs[i].ToUpperInvariant()), $"set {r.Id} row {i}");
        }
    }

    [Test]
    public void MultipleAlignClustalW_Contracts()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceAligner.MultipleAlignClustalW(Array.Empty<string>()).AlignedSequences, Is.Empty);
            var one = SequenceAligner.MultipleAlignClustalW(new[] { "ACGT" });
            Assert.That(one.AlignedSequences, Is.EqualTo(new[] { "ACGT" }));
            Assert.That(one.TotalScore, Is.Zero);
            Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignClustalW((IEnumerable<string>)null!));
            Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignClustalW((IEnumerable<DnaSequence>)null!));
            Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignClustalW(new[] { "ACGT", "ACQT" }));
            Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignClustalW(new[] { "ACGT", "" }));
            Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignClustalW(new[] { "ACGT", "ACGA" }, new ClustalWOptions { GapOpen = -1 }));
            Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignClustalW(new[] { "ACGT", "ACGA" }, new ClustalWOptions { PairwiseGapExtend = double.NaN }));
            Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignClustalW(new[] { "ACGT", "ACGA" }, new ClustalWOptions { TransitionWeight = -0.5f }));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.MultipleAlignClustalW(new[] { "ACGT", "ACGA" }, new ClustalWOptions { PairwiseMode = (ClustalWPairwiseMode)9 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.MultipleAlignClustalW(new[] { "ACGT", "ACGA" }, new ClustalWOptions { PairwiseMode = ClustalWPairwiseMode.Quick, Quick = new ClustalWQuickOptions { KTuple = 7 } }));
        });
    }

    [Test]
    public void MultipleAlignClustalW_DnaSequenceOverload_EqualsStringOverload_AndLowerCase()
    {
        string[] seqs = ["ACGTACGTTAGCTAGCTAGGATCGATCG", "ACGTACGATAGCTAGCTGGATCGATCA", "TTGACCGTAGGCATCGATCGGATAGCAG", "ACGTACGTTAGCTAGGATCGATCGAAA"];
        var viaString = SequenceAligner.MultipleAlignClustalW(seqs);
        var viaDna = SequenceAligner.MultipleAlignClustalW(seqs.Select(s => new DnaSequence(s)));
        var lower = SequenceAligner.MultipleAlignClustalW(seqs.Select(s => s.ToLowerInvariant()));
        Assert.Multiple(() =>
        {
            Assert.That(viaDna.AlignedSequences, Is.EqualTo(viaString.AlignedSequences));
            Assert.That(lower.AlignedSequences, Is.EqualTo(viaString.AlignedSequences));
            Assert.That(viaString.TotalScore, Is.EqualTo(SequenceAligner.SumOfPairsScore(viaString.AlignedSequences, SequenceAligner.SimpleDna, MsaGapModel.Affine)));
            Assert.That(viaString.Consensus.Length, Is.EqualTo(viaString.AlignedSequences[0].Length));
        });
    }

    [Test]
    public void MultipleAlignClustalW_DefaultOptions_AreClustalWDnaDefaults()
    {
        var o = new ClustalWOptions();
        Assert.Multiple(() =>
        {
            Assert.That(o.GapOpen, Is.EqualTo(15.0));
            Assert.That(o.GapExtend, Is.EqualTo(6.66));
            Assert.That(o.PairwiseGapOpen, Is.EqualTo(15.0));
            Assert.That(o.PairwiseGapExtend, Is.EqualTo(6.66));
            Assert.That(o.TransitionWeight, Is.EqualTo(0.5f));
            Assert.That(o.DivergenceCutoff, Is.EqualTo(30));
            Assert.That(o.NoWeights, Is.False);
            Assert.That(o.NegativeMatrix, Is.False);
            Assert.That(o.PairwiseMode, Is.EqualTo(ClustalWPairwiseMode.Full));
        });
        Assert.That(Align(["ACGTACGT", "ACGTTCGT", "ACGAACGT"]), Is.EqualTo(Align(["ACGTACGT", "ACGTTCGT", "ACGAACGT"], o)));
    }

    /// <summary>
    /// ClustalW stores -GAPEXT as a float: 100 * 6.66f = 665.99998 and the profile extension penalty is (int) 665, not 666
    /// (the pairwise stage multiplies in single precision and gets 666). The oracle fixes the rows; this documents the consequence:
    /// 6.66 and 6.65 are NOT the same parameter for ClustalW, but 6.66 and 6.6599998 are.
    /// </summary>
    [Test]
    public void MultipleAlignClustalW_GapExtend_IsStoredAsSinglePrecision()
    {
        string[] seqs = ["CCAACCAACCAAGGTT", "ACACTCCGCCCGCACTCAGAACTTT", "TAATTGG", "CCAACCGGTTAAGGCC"];
        Assert.That(Align(seqs, new ClustalWOptions { GapExtend = 6.66 }), Is.EqualTo(Align(seqs, new ClustalWOptions { GapExtend = 6.6599998 })));
    }
}
