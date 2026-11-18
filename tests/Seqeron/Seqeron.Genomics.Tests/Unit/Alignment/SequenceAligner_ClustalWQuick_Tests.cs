using System.Globalization;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A1-15b, B09 F55): ClustalW's FAST pairwise stage (<c>-QUICKTREE</c>, k-tuple diagonal scoring,
/// <see cref="SequenceAligner.ClustalWQuickPairwiseAlign(string, string, ClustalWQuickOptions?, int?)"/>) and
/// <see cref="MsaDistanceMode.ClustalWQuick"/>.
///
/// Oracle (TestData/Alignment/clustalw_quick_oracle.tsv, scripts/alignment/generate_clustalw_quick_oracle.py): the installed
/// <c>clustalw</c> 2.1 <c>-TYPE=DNA -QUICKTREE -ALIGN</c> "Sequences (i:j) Aligned. Score: X" (%lg of the double percent
/// score) / "Not Aligned" on 8 000+ pairs of 600+ sets with -KTUPLE / -WINDOW / -TOPDIAGS / -PAIRGAP / -SCORE variants, IUB-N
/// sequences, repeat-rich sets that reach ClustalW's "(Partial alignment)" bound.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_ClustalWQuick_Tests
{
    private static IEnumerable<string[]> ReadTsv(string resource)
    {
        string name = "Seqeron.Genomics.Tests.TestData.Alignment." + resource;
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        reader.ReadLine(); // header
        string? line;
        while ((line = reader.ReadLine()) is not null)
            if (line.Length > 0)
                yield return line.Split('\t');
    }

    private static readonly Lazy<List<string[]>> Rows = new(() => ReadTsv("clustalw_quick_oracle.tsv").ToList());

    private static ClustalWQuickOptions Opts(string[] f) => new()
    {
        KTuple = int.Parse(f[1]),
        Window = int.Parse(f[2]),
        TopDiagonals = int.Parse(f[3]),
        PairGap = int.Parse(f[4]),
        Percent = f[5] == "1",
    };

    [Test]
    public void ClustalWQuick_PrintedScore_MatchesClustalW21Oracle()
    {
        var rows = Rows.Value;
        Assert.That(rows.Count, Is.GreaterThan(8000));
        int bad = 0;
        string? first = null;
        foreach (var f in rows)
        {
            var r = SequenceAligner.ClustalWQuickPairwiseAlign(f[7], f[8], Opts(f), int.Parse(f[6]));
            bool ok;
            if (f[9] == "NA")
                ok = !r.Aligned;
            else
            {
                double expected = double.Parse(f[9], CultureInfo.InvariantCulture);
                // %lg prints 6 significant digits
                ok = r.Aligned && Math.Abs(r.Score - expected) <= 5.1e-6 * Math.Max(1.0, Math.Abs(expected));
            }

            if (!ok && ++bad == 1)
                first = $"{string.Join('/', f.Take(7))} {f[7]} / {f[8]}: expected {f[9]}, got {r.Score} (aligned {r.Aligned}, fragment {r.FragmentScore})";
        }

        Assert.That(bad, Is.Zero, first);
    }

    /// <summary>The oracle rows that ran with ClustalW's DNA defaults (2 / 4 / 4 / 5, percent) must be reproduced by a parameterless call.</summary>
    [Test]
    public void ClustalWQuick_DefaultOptions_MatchClustalWDefaultRuns()
    {
        var rows = Rows.Value.Where(f => f[1] == "2" && f[2] == "4" && f[3] == "4" && f[4] == "5" && f[5] == "1").ToList();
        Assert.That(rows.Count, Is.GreaterThan(1500));
        int bad = 0;
        foreach (var f in rows)
        {
            var r = SequenceAligner.ClustalWQuickPairwiseAlign(f[7], f[8], null, int.Parse(f[6]));
            bool ok = f[9] == "NA" ? !r.Aligned : r.Aligned && Math.Abs(r.Score - double.Parse(f[9], CultureInfo.InvariantCulture)) <= 5.1e-6 * Math.Max(1.0, r.Score);
            if (!ok) bad++;
        }

        Assert.That(bad, Is.Zero);
    }

    [Test]
    public void ClustalWQuick_OracleCoversTheParameterSpaceAndTheEdgeExits()
    {
        var rows = Rows.Value;
        Assert.Multiple(() =>
        {
            Assert.That(rows.Select(f => f[0]).Distinct().Order(), Is.EqualTo(new[] { "acgt", "iupac", "long", "rep" }));
            Assert.That(rows.Select(f => f[1]).Distinct().Count(), Is.EqualTo(4), "k-tuple 1..4");
            Assert.That(rows.Count(f => f[5] == "0"), Is.GreaterThan(500), "-SCORE=ABSOLUTE");
            Assert.That(rows.Count(f => f[9] == "NA"), Is.GreaterThan(500), "Not Aligned");
            Assert.That(rows.Count(f => f[4] != (int.Parse(f[1]) + 4).ToString()), Is.GreaterThan(500), "explicit -PAIRGAP");
            int partial = rows.Count(f => SequenceAligner.ClustalWQuickPairwiseAlign(f[7], f[8], Opts(f), int.Parse(f[6])).PartialAlignment);
            Assert.That(partial, Is.GreaterThan(0), "the oracle must reach ClustalW's '(Partial alignment)' exit");
        });
    }

    /// <summary>Hand run through the binary: `clustalw -QUICKTREE -ALIGN` prints 92.5926 / 27.7778 / 27.7778.</summary>
    [Test]
    public void ClustalWQuick_HandRunExample_PrintsSameScores()
    {
        string a = "ACGTACGTTTGACCAGTAGGCATCGAT", b = "ACGTACGATTGACCAGTAGCCATCGAT", c = "GGGGGGAAAATTTTCCCC";
        int longest = 27;
        Assert.Multiple(() =>
        {
            Assert.That(SequenceAligner.ClustalWQuickPairwiseAlign(a, b, null, longest).Score, Is.EqualTo(92.5926).Within(5e-5));
            Assert.That(SequenceAligner.ClustalWQuickPairwiseAlign(a, c, null, longest).Score, Is.EqualTo(27.7778).Within(5e-5));
            Assert.That(SequenceAligner.ClustalWQuickPairwiseAlign(b, c, null, longest).Score, Is.EqualTo(27.7778).Within(5e-5));
        });
    }

    [Test]
    public void ClustalWQuick_ResultIsSelfConsistent_DistanceFromScore()
    {
        foreach (var f in Rows.Value.Take(500))
        {
            var r = SequenceAligner.ClustalWQuickPairwiseAlign(f[7], f[8], Opts(f), int.Parse(f[6]));
            double calc = r.FragmentScore == 0 ? 0.0
                : f[5] == "1" ? (double)r.FragmentScore / Math.Min(f[7].Length, f[8].Length) * 100.0 : r.FragmentScore;
            Assert.That(r.Score, Is.EqualTo(calc).Within(1e-12));
            Assert.That(r.Distance, Is.EqualTo((100.0 - calc) / 100.0).Within(1e-12));
        }
    }

    [Test]
    public void ClustalWQuick_Defaults_AreClustalWDnaDefaults()
    {
        var o = new ClustalWQuickOptions();
        Assert.Multiple(() =>
        {
            Assert.That(o.KTuple, Is.EqualTo(2));
            Assert.That(o.Window, Is.EqualTo(4));
            Assert.That(o.TopDiagonals, Is.EqualTo(4));
            Assert.That(o.PairGap, Is.EqualTo(5));
            Assert.That(o.Percent, Is.True);
        });

        // clustalw -KTUPLE=2 alone sets the pair gap to 6 (CommandLineParser), the default without -KTUPLE is 5: the printed
        // score of this pair is 27.5862 with pair gap 5 and 25.8621 with pair gap 6 (binary).
        const string a = "GACTCCTTGTCCACTTTGTTGAGAGGCAAGCAAAAAGCATTATGTACCCCACTCCCCT", b = "GGCCCCGTACACCTTTGGGCGTCGTGGCTAGAGGGTAAGTTTTAAGGAAACCATCCGCT";
        Assert.That(SequenceAligner.ClustalWQuickPairwiseAlign(a, b).Score, Is.EqualTo(27.5862).Within(5e-5));
        Assert.That(SequenceAligner.ClustalWQuickPairwiseAlign(a, b, o with { PairGap = 6 }).Score, Is.EqualTo(25.8621).Within(5e-5));
    }

    [Test]
    public void ClustalWQuick_Contracts()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => SequenceAligner.ClustalWQuickPairwiseAlign((string)null!, "A"));
            Assert.Throws<ArgumentException>(() => SequenceAligner.ClustalWQuickPairwiseAlign("ACGQ", "ACGT"));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.ClustalWQuickPairwiseAlign("ACGT", "ACGT", new ClustalWQuickOptions { KTuple = 0 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.ClustalWQuickPairwiseAlign("ACGT", "ACGT", new ClustalWQuickOptions { KTuple = 5 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.ClustalWQuickPairwiseAlign("ACGT", "ACGT", new ClustalWQuickOptions { Window = 0 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.ClustalWQuickPairwiseAlign("ACGT", "ACGT", new ClustalWQuickOptions { TopDiagonals = 0 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.ClustalWQuickPairwiseAlign("ACGT", "ACGT", new ClustalWQuickOptions { PairGap = 0 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.ClustalWQuickPairwiseAlign("ACGTACGT", "ACGT", null, 3));
            var empty = SequenceAligner.ClustalWQuickPairwiseAlign("", "ACGT");
            Assert.That((empty.Distance, empty.Aligned, empty.FragmentScore), Is.EqualTo((1.0, false, 0)));
            var shortSeq = SequenceAligner.ClustalWQuickPairwiseAlign("A", "A");
            Assert.That(shortSeq.Aligned, Is.False, "shorter than the k-tuple: no fragment");
        });
    }

    #region MSA with ClustalWQuick distances

    private static IEnumerable<DnaSequence> Dna(params string[] s) => s.Select(x => new DnaSequence(x));

    [Test]
    public void MultipleAlign_QuickDistance_ProducesValidAlignments_AndDiffersFromFullOnSomeSets()
    {
        var rnd = new Random(9);
        int differs = 0;
        for (int t = 0; t < 60; t++)
        {
            string baseSeq = new(Enumerable.Range(0, rnd.Next(20, 70)).Select(_ => "ACGT"[rnd.Next(4)]).ToArray());
            var seqs = Enumerable.Range(0, rnd.Next(3, 7)).Select(_ =>
            {
                var chars = baseSeq.ToList();
                for (int m = 0; m < rnd.Next(1, 12); m++)
                {
                    int p = rnd.Next(chars.Count);
                    if (rnd.Next(3) == 0 && chars.Count > 10) chars.RemoveAt(p); else chars[p] = "ACGT"[rnd.Next(4)];
                }

                return new string(chars.ToArray());
            }).ToArray();

            var quick = SequenceAligner.MultipleAlignProgressive(Dna(seqs), null, new MsaOptions { Distance = MsaDistanceMode.ClustalWQuick, GuideTree = MsaGuideTree.NeighborJoining });
            var full = SequenceAligner.MultipleAlignProgressive(Dna(seqs), null, new MsaOptions { Distance = MsaDistanceMode.ClustalWAffine, GuideTree = MsaGuideTree.NeighborJoining });
            for (int i = 0; i < seqs.Length; i++)
                Assert.That(quick.AlignedSequences[i].Replace("-", ""), Is.EqualTo(seqs[i]));
            Assert.That(quick.AlignedSequences.Select(r => r.Length).Distinct().Count(), Is.EqualTo(1));
            if (!quick.AlignedSequences.SequenceEqual(full.AlignedSequences)) differs++;
        }

        Assert.That(differs, Is.GreaterThan(0), "the fast distances must be able to change the guide tree");
    }

    [Test]
    public void MultipleAlign_QuickDistance_UsesQuickPairwiseParameters()
    {
        var seqs = new[]
        {
            "ACGTACGTTTGACCAGTAGGCATCGATGGA", "ACGTACGATTGACCAGTAGCCATCGATGGA", "TTGACCAGTAGGCATCGATCCCAAGT",
            "ACGTAGGATTGACCGTAGCCATCGATGGA", "GGGTTTACCAGTAGGCATCGATGGAAACCT",
        };
        string[] Run(ClustalWQuickOptions q) => SequenceAligner.MultipleAlignProgressive(Dna(seqs), null,
            new MsaOptions { Distance = MsaDistanceMode.ClustalWQuick, QuickPairwise = q }).AlignedSequences.ToArray();

        // Default options reproduce the explicit ClustalW defaults.
        Assert.That(Run(new ClustalWQuickOptions { KTuple = 2, Window = 4, TopDiagonals = 4, PairGap = 5 }), Is.EqualTo(Run(new ClustalWQuickOptions())));
        // Option validation reaches the guide-distance stage.
        Assert.Throws<ArgumentOutOfRangeException>(() => Run(new ClustalWQuickOptions { KTuple = 9 }));
    }

    [Test]
    public void MultipleAlign_DefaultDistance_IsUnchanged()
    {
        var seqs = new[] { "ACGTACGTAC", "ACGTTCGTAC", "ACGTACGAAC" };
        Assert.That(new MsaOptions().Distance, Is.EqualTo(MsaDistanceMode.Linear));
        Assert.That(SequenceAligner.MultipleAlignProgressive(Dna(seqs)).AlignedSequences,
            Is.EqualTo(SequenceAligner.MultipleAlignProgressive(Dna(seqs), null, new MsaOptions()).AlignedSequences));
    }

    #endregion
}
