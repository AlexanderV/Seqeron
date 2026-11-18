namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A1-15, B09 F51): ClustalW's own pairwise stage as the guide-tree distance of the progressive,
/// iterative and consistency aligners (<see cref="MsaDistanceMode.ClustalWAffine"/>) and the public
/// <see cref="SequenceAligner.ClustalWPairwiseAlign(string, string, double, double)"/>.
///
/// Oracle (TestData/Alignment/clustalw_pairwise_oracle.tsv, clustalw_msa_oracle.tsv; scripts/alignment/generate_clustalw_pairwise_oracle.py):
/// the installed <c>clustalw</c> 2.1 (<c>-TYPE=DNA -ALIGN</c> prints "Sequences (i:j) Aligned. Score: N" with
/// N = (int)(100 * identities / MIN(len1, len2)) of FullPairwiseAlign's alignment) on 3 600 random / low-complexity / mutated
/// DNA and IUB-ambiguity pairs; the identity count and score columns come from a line-by-line Python port of
/// FullPairwiseAlign.cpp that was asserted equal to the binary on every row. The MSA rows come from an independent Python
/// progressive / iterative implementation driven by that distance. ClustalW's pairwise stage is an affine local
/// (Smith-Waterman-Gotoh) alignment, end = first row-major maximum, Myers-Miller midpoint tie rule.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_ClustalWDistance_Tests
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

    private static readonly Lazy<List<string[]>> Pairs = new(() => ReadTsv("clustalw_pairwise_oracle.tsv").ToList());
    private static readonly Lazy<List<string[]>> MsaRows = new(() => ReadTsv("clustalw_msa_oracle.tsv").ToList());

    #region Pairwise stage vs clustalw 2.1

    [Test]
    public void ClustalWPairwise_PrintedScoreIdentitiesAndScore_MatchClustalW21Oracle()
    {
        var rows = Pairs.Value;
        Assert.That(rows.Count, Is.GreaterThan(3000));
        int bad = 0;
        string? first = null;
        foreach (var f in rows)
        {
            var r = SequenceAligner.ClustalWPairwiseAlign(f[1], f[2]);
            bool ok = r.ClustalWScore == int.Parse(f[3]) && r.Identities == int.Parse(f[4]) && r.Score == int.Parse(f[5]);
            if (!ok && ++bad == 1) first = $"{f[1]} / {f[2]}: expected {f[3]}/{f[4]}/{f[5]}, got {r.ClustalWScore}/{r.Identities}/{r.Score}";
        }

        Assert.That(bad, Is.Zero, first);
    }

    [Test]
    public void ClustalWPairwise_OracleCoversAcgtIupacAndLongPairs()
    {
        var kinds = Pairs.Value.GroupBy(f => f[0]).ToDictionary(g => g.Key, g => g.Count());
        Assert.Multiple(() =>
        {
            Assert.That(kinds["acgt"], Is.GreaterThan(2000));
            Assert.That(kinds["iupac"], Is.GreaterThan(500));
            Assert.That(kinds["long"], Is.GreaterThan(10));
        });
    }

    /// <summary>Three pairs run through the binary by hand: clustalw -ALIGN prints Score 91, 40 and 26.</summary>
    [TestCase("ACGTACGTTAGCTAGCTAGGATCGATCGATTTAGC", "ACGTACGATAGCTAGCTGGATCGATCAATTTAGCAA", 91, 32)]
    [TestCase("ACGTACGTTAGCTAGCTAGGATCGATCGATTTAGC", "TTGACCGTAGGCATCGATCGGATAGCAGCT", 40, 12)]
    [TestCase("ACGTACGATAGCTAGCTGGATCGATCAATTTAGCAA", "TTGACCGTAGGCATCGATCGGATAGCAGCT", 26, 8)]
    public void ClustalWPairwise_HandRunExamples_PrintedScore(string a, string b, int printed, int identities)
    {
        var r = SequenceAligner.ClustalWPairwiseAlign(a, b);
        Assert.Multiple(() =>
        {
            Assert.That(r.ClustalWScore, Is.EqualTo(printed));
            Assert.That(r.Identities, Is.EqualTo(identities));
            Assert.That(r.PercentIdentity, Is.EqualTo(100.0 * identities / Math.Min(a.Length, b.Length)).Within(1e-12));
            Assert.That(r.Distance, Is.EqualTo(1.0 - (double)identities / Math.Min(a.Length, b.Length)).Within(1e-12));
        });
    }

    [Test]
    public void ClustalWPairwise_Alignment_IsLocalRegionConsistentWithCounts()
    {
        foreach (var f in Pairs.Value.Where(f => f[0] == "acgt").Take(600))
        {
            var r = SequenceAligner.ClustalWPairwiseAlign(f[1], f[2]);
            if (r.Score == 0)
            {
                Assert.That(r.AlignedSequence1, Is.Empty);
                Assert.That((r.Start1, r.End1, r.Start2, r.End2), Is.EqualTo((-1, -1, -1, -1)));
                continue;
            }

            Assert.Multiple(() =>
            {
                Assert.That(r.AlignedSequence1.Length, Is.EqualTo(r.AlignedSequence2.Length));
                Assert.That(r.AlignedSequence1.Replace("-", ""), Is.EqualTo(f[1][r.Start1..(r.End1 + 1)]));
                Assert.That(r.AlignedSequence2.Replace("-", ""), Is.EqualTo(f[2][r.Start2..(r.End2 + 1)]));
                int same = 0;
                for (int i = 0; i < r.AlignedSequence1.Length; i++)
                    if (r.AlignedSequence1[i] != '-' && r.AlignedSequence1[i] == r.AlignedSequence2[i]) same++;
                Assert.That(same, Is.EqualTo(r.Identities));
            });
        }
    }

    /// <summary>The path score recomputed from the rows with ClustalW's matrix (+1000 match, -900 mismatch) and gaps 3000 + 666 k equals the DP score.</summary>
    [Test]
    public void ClustalWPairwise_AlignmentRescoresToScore()
    {
        foreach (var f in Pairs.Value.Where(f => f[0] == "acgt").Take(800))
        {
            var r = SequenceAligner.ClustalWPairwiseAlign(f[1], f[2]);
            if (r.Score == 0) continue;
            long score = 0;
            char open = '\0';
            for (int i = 0; i < r.AlignedSequence1.Length; i++)
            {
                char x = r.AlignedSequence1[i], y = r.AlignedSequence2[i];
                if (x == '-' || y == '-')
                {
                    char kind = x == '-' ? 'a' : 'b';
                    score -= (kind == open ? 0 : 3000) + 666;
                    open = kind;
                }
                else
                {
                    score += x == y ? 1000 : -900;
                    open = '\0';
                }
            }

            Assert.That(score, Is.EqualTo(r.Score), $"{f[1]} / {f[2]}");
        }
    }

    [Test]
    public void ClustalWPairwise_EmptySequence_DistanceOne()
    {
        foreach (var (a, b) in new[] { ("", "ACGT"), ("ACGT", ""), ("", "") })
        {
            var r = SequenceAligner.ClustalWPairwiseAlign(a, b);
            Assert.Multiple(() =>
            {
                Assert.That(r.Distance, Is.EqualTo(1.0));
                Assert.That(r.PercentIdentity, Is.Zero);
                Assert.That(r.Identities, Is.Zero);
                Assert.That(r.ClustalWScore, Is.Zero);
            });
        }
    }

    [Test]
    public void ClustalWPairwise_IdenticalAndLowerCase_AreFullIdentity()
    {
        var r = SequenceAligner.ClustalWPairwiseAlign("acgtacgtac", "ACGTACGTAC");
        Assert.Multiple(() =>
        {
            Assert.That(r.Identities, Is.EqualTo(10));
            Assert.That(r.ClustalWScore, Is.EqualTo(100));
            Assert.That(r.Score, Is.EqualTo(10_000));
        });
    }

    [Test]
    public void ClustalWPairwise_DnaSequenceOverload_EqualsStringOverload()
    {
        var a = new DnaSequence("ACGTACGTTAGCTAGC");
        var b = new DnaSequence("ACGTTCGATAGCTGC");
        Assert.That(SequenceAligner.ClustalWPairwiseAlign(a, b), Is.EqualTo(SequenceAligner.ClustalWPairwiseAlign(a.Sequence, b.Sequence)));
    }

    [Test]
    public void ClustalWPairwise_InvalidInput_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.ClustalWPairwiseAlign((string)null!, "ACGT"));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.ClustalWPairwiseAlign("ACGT", (string)null!));
        Assert.Throws<ArgumentException>(() => SequenceAligner.ClustalWPairwiseAlign("ACGT", "AC1T"));
        Assert.Throws<ArgumentException>(() => SequenceAligner.ClustalWPairwiseAlign("ACGT", "ACGT", -1, 6.66));
        Assert.Throws<ArgumentException>(() => SequenceAligner.ClustalWPairwiseAlign("ACGT", "ACGT", 15, double.NaN));
    }

    /// <summary>Gap penalties are parameters: a free gap (0/0) lets the aligner bridge the insertion that the defaults refuse.</summary>
    [Test]
    public void ClustalWPairwise_GapPenaltiesAreHonoured()
    {
        const string a = "ACGTACGTAC" + "TTTTTTTTTT" + "GGCATGCATG";
        const string b = "ACGTACGTAC" + "GGCATGCATG";
        var strict = SequenceAligner.ClustalWPairwiseAlign(a, b);
        var free = SequenceAligner.ClustalWPairwiseAlign(a, b, 0, 0);
        Assert.Multiple(() =>
        {
            Assert.That(free.Identities, Is.EqualTo(20));
            Assert.That(free.Score, Is.EqualTo(20_000));
            Assert.That(strict.Score, Is.LessThan(free.Score));
        });
    }

    #endregion

    #region Guide-tree distance mode

    private static ScoringMatrix Scoring(string[] f) => new(int.Parse(f[0]), int.Parse(f[1]), 0, int.Parse(f[2]));

    [Test]
    public void ProgressiveAndIterative_ClustalWDistance_MatchIndependentPython()
    {
        var options = new MsaOptions { Distance = MsaDistanceMode.ClustalWAffine };
        int badProg = 0, badIter = 0, differFromLinear = 0;
        foreach (var f in MsaRows.Value)
        {
            var seqs = f[3].Split(';').Select(s => new DnaSequence(s)).ToList();
            var scoring = Scoring(f);
            var prog = SequenceAligner.MultipleAlignProgressive(seqs, scoring, options);
            var iter = SequenceAligner.MultipleAlignIterative(seqs, scoring, options);
            if (!prog.AlignedSequences.SequenceEqual(f[5].Split(';'))) badProg++;
            if (!iter.AlignedSequences.SequenceEqual(f[6].Split(';'))) badIter++;
            if (!prog.AlignedSequences.SequenceEqual(SequenceAligner.MultipleAlignProgressive(seqs, scoring).AlignedSequences)) differFromLinear++;
        }

        Assert.Multiple(() =>
        {
            Assert.That(badProg, Is.Zero, "progressive rows vs independent Python");
            Assert.That(badIter, Is.Zero, "iterative rows vs independent Python");
            Assert.That(differFromLinear, Is.GreaterThan(5), "the oracle must discriminate the two distance modes");
        });
    }

    [Test]
    public void LinearDistanceMode_IsTheDefault_AndUnchanged()
    {
        foreach (var f in MsaRows.Value.Take(40))
        {
            var seqs = f[3].Split(';').Select(s => new DnaSequence(s)).ToList();
            var scoring = Scoring(f);
            var def = new MsaOptions();
            Assert.Multiple(() =>
            {
                Assert.That(def.Distance, Is.EqualTo(MsaDistanceMode.Linear));
                Assert.That(SequenceAligner.MultipleAlignProgressive(seqs, scoring, def).AlignedSequences,
                    Is.EqualTo(SequenceAligner.MultipleAlignProgressive(seqs, scoring).AlignedSequences));
                Assert.That(SequenceAligner.MultipleAlignIterative(seqs, scoring, def).AlignedSequences,
                    Is.EqualTo(SequenceAligner.MultipleAlignIterative(seqs, scoring).AlignedSequences));
                Assert.That(SequenceAligner.MultipleAlignConsistency(seqs, scoring, null).AlignedSequences,
                    Is.EqualTo(SequenceAligner.MultipleAlignConsistency(seqs, scoring).AlignedSequences));
            });
        }
    }

    /// <summary>The distance mode only chooses the guide tree: every result still spells its inputs, rows are equal length.</summary>
    [Test]
    public void ConsistencyWithClustalWDistance_RowsSpellInputs()
    {
        var options = new MsaOptions { Distance = MsaDistanceMode.ClustalWAffine };
        foreach (var f in MsaRows.Value.Take(60))
        {
            var raw = f[3].Split(';');
            var seqs = raw.Select(s => new DnaSequence(s)).ToList();
            var r = SequenceAligner.MultipleAlignConsistency(seqs, Scoring(f), options);
            Assert.That(r.AlignedSequences.Select(x => x.Replace("-", "")), Is.EqualTo(raw));
            Assert.That(r.AlignedSequences.Select(x => x.Length).Distinct().Count(), Is.EqualTo(1));
        }
    }

    [Test]
    public void InvalidDistanceMode_Throws()
    {
        var seqs = new[] { new DnaSequence("ACGT"), new DnaSequence("ACGA"), new DnaSequence("ACGG") };
        var bad = new MsaOptions { Distance = (MsaDistanceMode)99 };
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.MultipleAlignProgressive(seqs, null, bad));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.MultipleAlignIterative(seqs, null, bad));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.MultipleAlignConsistency(seqs, null, bad));
    }

    [Test]
    public void TwoSequences_ClustalWDistance_SameAsLinearResult()
    {
        var seqs = new[] { new DnaSequence("ACGTACGTAC"), new DnaSequence("ACGTTCGTAC") };
        var options = new MsaOptions { Distance = MsaDistanceMode.ClustalWAffine };
        Assert.That(SequenceAligner.MultipleAlignProgressive(seqs, null, options).AlignedSequences,
            Is.EqualTo(SequenceAligner.MultipleAlignProgressive(seqs).AlignedSequences));
    }

    #endregion
}
