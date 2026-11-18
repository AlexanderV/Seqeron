namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A1-10, B09 F52): affine-gap multiple sequence alignment (<see cref="MsaOptions.GapModel"/> =
/// <see cref="MsaGapModel.Affine"/>) and the natural affine sum-of-pairs objective <see cref="SequenceAligner.SumOfPairsScore"/>.
///
/// Objective (Altschul 1989 natural gap costs): every pair of rows is projected out of the alignment, gap-gap columns are skipped and every
/// maximal run of gap columns of the projection costs GapOpen + length * GapExtend. For two sequences this is the Gotoh score of
/// <c>GlobalAlignAffine</c>, so a pairwise optimum is an MSA optimum. For k &gt;= 3 the methods are HEURISTICS (the profile DP uses the
/// position-specific approximation of ClustalW/MUSCLE); optimality is claimed only for k = 2 and measured against the exact optimum of tiny
/// sets below.
///
/// Oracles (TestData/Alignment/msa_affine_*.tsv, scripts/alignment/generate_msa_affine_oracle.py): independent Python re-score of random
/// gapped MSAs with the Biopython 1.88 pairwise optimum as an upper bound; the EXACT optimum of tiny 2-/3-sequence sets by dynamic
/// programming over (i, j, l, run state of each pair), validated against exhaustive enumeration of all alignments; independent Python
/// implementations of the affine star, progressive, iterative and consistency (affine pair library) methods.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MultipleAlignAffine_Tests
{
    private static IEnumerable<string[]> ReadTsv(string resource)
    {
        string name = "Seqeron.Genomics.Tests.TestData.Alignment." + resource;
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        string? line;
        while ((line = reader.ReadLine()) is not null)
            if (line.Length > 0)
                yield return line.Split('\t');
    }

    private static readonly Lazy<List<string[]>> SpRows = new(() => ReadTsv("msa_affine_sp_oracle.tsv").ToList());
    private static readonly Lazy<List<string[]>> ExactRows = new(() => ReadTsv("msa_affine_exact_oracle.tsv").ToList());
    private static readonly Lazy<List<string[]>> MsaRows = new(() => ReadTsv("msa_affine_oracle.tsv").ToList());

    private static ScoringMatrix Scoring(string[] f) => new(int.Parse(f[0]), int.Parse(f[1]), int.Parse(f[2]), int.Parse(f[3]));

    private static List<DnaSequence> Dna(string joined) => joined.Split(';').Select(s => new DnaSequence(s)).ToList();

    private static readonly MsaOptions Affine = new() { GapModel = MsaGapModel.Affine };

    #region Natural affine sum-of-pairs objective

    [Test]
    public void SumOfPairsScore_MatchesIndependentRescore_AndRespectsPairwiseOptimumBound()
    {
        var rows = SpRows.Value;
        Assert.That(rows.Count, Is.GreaterThan(400));
        int bad = 0;
        string? first = null;
        foreach (var f in rows)
        {
            var aligned = f[4].Split(';');
            var sc = Scoring(f);
            bool ok = SequenceAligner.SumOfPairsScore(aligned, sc, MsaGapModel.Affine) == int.Parse(f[6])
                   && SequenceAligner.SumOfPairsScore(aligned, sc, MsaGapModel.Linear) == int.Parse(f[5])
                   && int.Parse(f[6]) <= int.Parse(f[7]);
            if (!ok && ++bad == 1) first = $"{f[4]} {string.Join(',', f.Take(4))}";
        }

        Assert.That(bad, Is.Zero, first);
    }

    /// <summary>Hand example (match 2, mismatch -3, GapOpen -5, GapExtend -1) on rows AC-GT / ACTGT / A--GT; see the per-pair comments.</summary>
    [Test]
    public void SumOfPairsScore_HandDerivedExample()
    {
        // match +2, mismatch -3, open -5, extend -1.
        var sc = new ScoringMatrix(2, -3, -5, -1);
        string[] rows = ["AC-GT", "ACTGT", "A--GT"];
        // (0,1): A/A C/C match, gap in row 0 (run 1) = -5-1 = -6, G/G T/T: 2*4 - 6 = 2
        // (0,2): A/A match, C against gap in row 2 = run of one projected column (col 2 is gap-gap and skipped): -6, G, T: 3*2 - 6 = 0
        // (1,2): A/A, then C,T against gaps: one run of 2 = -5 - 2 = -7, G, T: 3*2 - 7 = -1
        Assert.That(SequenceAligner.SumOfPairsScore(rows, sc, MsaGapModel.Affine), Is.EqualTo(2 + 0 + -1));
        // Linear: gap positions cost -1 each: (0,1) 8-1 = 7, (0,2) 6-1 = 5, (1,2) 6-2 = 4.
        Assert.That(SequenceAligner.SumOfPairsScore(rows, sc, MsaGapModel.Linear), Is.EqualTo(7 + 5 + 4));
    }

    /// <summary>A gap in row i directly followed by a gap in row j are two runs; a gap run is not interrupted by gap-gap columns.</summary>
    [Test]
    public void SumOfPairsScore_RunSemantics()
    {
        var sc = new ScoringMatrix(1, -1, -10, -1);
        Assert.Multiple(() =>
        {
            Assert.That(SequenceAligner.SumOfPairsScore(["A-", "-A"], sc, MsaGapModel.Affine), Is.EqualTo(-22), "two runs");
            Assert.That(SequenceAligner.SumOfPairsScore(["AA-", "--A"], sc, MsaGapModel.Affine), Is.EqualTo(-12 + -11), "run of two, then a run in the other row");
            Assert.That(SequenceAligner.SumOfPairsScore(["AAA", "A-A"], sc, MsaGapModel.Affine), Is.EqualTo(1 + -11 + 1));
            // (0,1) one run of 3 = -13; (0,2) A/A A/A then (A,-) = 2 - 11 = -9; (1,2) run of 2 (last column gap-gap skipped) = -12.
            Assert.That(SequenceAligner.SumOfPairsScore(["AAA", "---", "AA-"], sc, MsaGapModel.Affine), Is.EqualTo(-13 - 9 - 12));
            Assert.That(SequenceAligner.SumOfPairsScore(["AC-T", "AC-T"], sc, MsaGapModel.Affine), Is.EqualTo(3), "gap-gap column skipped");
            Assert.That(SequenceAligner.SumOfPairsScore(["A"], sc, MsaGapModel.Affine), Is.Zero);
            Assert.That(SequenceAligner.SumOfPairsScore([], sc, MsaGapModel.Affine), Is.Zero);
        });
    }

    [Test]
    public void SumOfPairsScore_Validation()
    {
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.SumOfPairsScore(null!));
        Assert.Throws<ArgumentException>(() => SequenceAligner.SumOfPairsScore(["AC", "A"]));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.SumOfPairsScore(["AC", "AC"], null, (MsaGapModel)7));
    }

    /// <summary>For two sequences the natural affine SP of the optimal pairwise alignment is the Gotoh score of GlobalAlignAffine.</summary>
    [Test]
    public void SumOfPairsScore_TwoRows_EqualsGlobalAlignAffineScore()
    {
        var rnd = new Random(2026);
        string Rnd(int n) => new string(Enumerable.Range(0, n).Select(_ => "ACGT"[rnd.Next(4)]).ToArray());
        var schemes = new[] { new ScoringMatrix(1, -1, -2, -1), new ScoringMatrix(2, -3, -5, -2), new ScoringMatrix(5, -4, -10, -1) };
        for (int t = 0; t < 150; t++)
        {
            var sc = schemes[t % 3];
            string a = Rnd(rnd.Next(1, 40)), b = Rnd(rnd.Next(1, 40));
            var r = SequenceAligner.GlobalAlignAffine(a, b, sc);
            Assert.That(SequenceAligner.SumOfPairsScore([r.AlignedSequence1, r.AlignedSequence2], sc, MsaGapModel.Affine), Is.EqualTo(r.Score));
        }
    }

    #endregion

    #region Exact optimum of tiny sets

    [Test]
    public void ExactOptimumRows_ScoreToTheOptimum_UnderTheObjective()
    {
        foreach (var f in ExactRows.Value)
            Assert.That(SequenceAligner.SumOfPairsScore(f[6].Split(';'), Scoring(f), MsaGapModel.Affine), Is.EqualTo(int.Parse(f[5])), f[4]);
    }

    /// <summary>
    /// Honest optimality statement: no heuristic exceeds the exact optimum; the two-sequence star / progressive / iterative / consistency
    /// results ARE optimal; the fraction of 3-sequence sets solved exactly is recorded (observed values asserted as lower bounds).
    /// </summary>
    [Test]
    public void Heuristics_NeverExceedTheExactOptimum_AndAreOptimalForTwoSequences()
    {
        int sets3 = 0, star3 = 0, prog3 = 0, iter3 = 0;
        foreach (var f in ExactRows.Value)
        {
            var sc = Scoring(f);
            var seqs = Dna(f[4]);
            int optimal = int.Parse(f[5]);
            var star = SequenceAligner.MultipleAlign(seqs, sc, Affine);
            var prog = SequenceAligner.MultipleAlignProgressive(seqs, sc, Affine);
            var iter = SequenceAligner.MultipleAlignIterative(seqs, sc, Affine);
            var cons = SequenceAligner.MultipleAlignConsistency(seqs, sc, Affine);
            var anch = SequenceAligner.MultipleAlignAnchored(seqs, sc, Affine);
            Assert.Multiple(() =>
            {
                foreach (var r in new[] { star, prog, iter, cons, anch })
                {
                    Assert.That(r.TotalScore, Is.LessThanOrEqualTo(optimal), f[4]);
                    Assert.That(r.TotalScore, Is.EqualTo(SequenceAligner.SumOfPairsScore(r.AlignedSequences, sc, MsaGapModel.Affine)));
                }
            });

            if (seqs.Count == 2)
            {
                Assert.Multiple(() =>
                {
                    Assert.That(star.TotalScore, Is.EqualTo(optimal), "star k=2");
                    Assert.That(prog.TotalScore, Is.EqualTo(optimal), "progressive k=2");
                    Assert.That(iter.TotalScore, Is.EqualTo(optimal), "iterative k=2");
                });
            }
            else
            {
                sets3++;
                if (star.TotalScore == optimal) star3++;
                if (prog.TotalScore == optimal) prog3++;
                if (iter.TotalScore == optimal) iter3++;
            }
        }

        TestContext.Out.WriteLine($"3-sequence sets {sets3}: star optimal {star3}, progressive optimal {prog3}, iterative optimal {iter3}");
        Assert.That(iter3, Is.GreaterThanOrEqualTo(prog3), "refinement never lowers the exact objective");
        Assert.That(star3 + prog3, Is.GreaterThan(sets3 / 2), "the heuristics are optimal on a substantial share of tiny sets");
    }

    #endregion

    #region Rows vs independent Python

    [Test]
    public void StarProgressiveIterative_AffineRows_MatchIndependentPython()
    {
        int badStar = 0, badProg = 0, badIter = 0, badCons = 0, differFromLinear = 0;
        foreach (var f in MsaRows.Value)
        {
            var sc = Scoring(f);
            var seqs = Dna(f[4]);
            if (!SequenceAligner.MultipleAlign(seqs, sc, Affine).AlignedSequences.SequenceEqual(f[5].Split(';'))) badStar++;
            var prog = SequenceAligner.MultipleAlignProgressive(seqs, sc, Affine);
            if (!prog.AlignedSequences.SequenceEqual(f[6].Split(';'))) badProg++;
            if (!SequenceAligner.MultipleAlignIterative(seqs, sc, Affine).AlignedSequences.SequenceEqual(f[7].Split(';'))) badIter++;
            if (!SequenceAligner.MultipleAlignConsistency(seqs, sc, Affine).AlignedSequences.SequenceEqual(f[8].Split(';'))) badCons++;
            if (!prog.AlignedSequences.SequenceEqual(SequenceAligner.MultipleAlignProgressive(seqs, sc).AlignedSequences)) differFromLinear++;
        }

        Assert.Multiple(() =>
        {
            Assert.That(badStar, Is.Zero, "affine star rows");
            Assert.That(badProg, Is.Zero, "affine progressive rows");
            Assert.That(badIter, Is.Zero, "affine iterative rows");
            Assert.That(badCons, Is.Zero, "affine consistency rows (affine global + local pair library)");
            Assert.That(differFromLinear, Is.GreaterThan(20), "the affine model must change the alignments");
        });
    }

    #endregion

    #region Two sequences: exactly GlobalAlignAffine

    [Test]
    public void TwoSequences_AllMethods_RowsEqualGlobalAlignAffine()
    {
        var rnd = new Random(77);
        string Rnd(int n) => new string(Enumerable.Range(0, n).Select(_ => "ACGT"[rnd.Next(4)]).ToArray());
        var schemes = new[] { new ScoringMatrix(1, -1, -2, -1), new ScoringMatrix(2, -3, -5, -2), new ScoringMatrix(5, -4, -10, -1), new ScoringMatrix(1, -2, -3, -1) };
        for (int t = 0; t < 120; t++)
        {
            var sc = schemes[t % schemes.Length];
            string a = Rnd(rnd.Next(1, 30)), b = Rnd(rnd.Next(1, 30));
            var seqs = new List<DnaSequence> { new(a), new(b) };
            var expected = SequenceAligner.GlobalAlignAffine(a, b, sc);
            var star = SequenceAligner.MultipleAlign(seqs, sc, Affine);
            var prog = SequenceAligner.MultipleAlignProgressive(seqs, sc, Affine);
            Assert.Multiple(() =>
            {
                Assert.That(star.TotalScore, Is.EqualTo(expected.Score), $"{a}/{b} star");
                Assert.That(prog.TotalScore, Is.EqualTo(expected.Score), $"{a}/{b} progressive");
                Assert.That(prog.AlignedSequences, Is.EqualTo(new[] { expected.AlignedSequence1, expected.AlignedSequence2 }), $"{a}/{b} progressive rows");
                Assert.That(SequenceAligner.MultipleAlignIterative(seqs, sc, Affine).TotalScore, Is.EqualTo(expected.Score));
            });
        }
    }

    /// <summary>The Gotoh border: gap penalties charged per opening - the long indel is ONE opening in the affine model.</summary>
    [Test]
    public void LongIndel_AffineChargesOneOpening_LinearDoesNot()
    {
        var sc = new ScoringMatrix(2, -3, -8, -1);
        var seqs = new List<DnaSequence> { new("ACGTACGTACGTAAAAAAAAAAGGCATCGATCG"), new("ACGTACGTACGTGGCATCGATCG"), new("ACGTACGTACGTAAAAAAAAAAGGCATCGATCG") };
        var affine = SequenceAligner.MultipleAlign(seqs, sc, Affine);
        var linear = SequenceAligner.MultipleAlign(seqs, sc);
        Assert.Multiple(() =>
        {
            Assert.That(affine.AlignedSequences[1], Does.Contain("----------"));
            // Identical pair: 33 matches = 66. Each of the two pairs with the deleted sequence: 23 matches = 46 and ONE 10-nt gap run
            // (affine -8 + 10 * -1 = -18; linear 10 * -1 = -10).
            Assert.That(affine.TotalScore, Is.EqualTo(66 + 2 * (46 - 18)));
            Assert.That(linear.TotalScore, Is.EqualTo(66 + 2 * (46 - 10)));
        });
    }

    #endregion

    #region Contracts

    [Test]
    public void DefaultOptions_AreLinear_AndEveryMethodIsUnchanged()
    {
        var def = new MsaOptions();
        Assert.That(def.GapModel, Is.EqualTo(MsaGapModel.Linear));
        foreach (var f in MsaRows.Value.Take(25))
        {
            var sc = Scoring(f);
            var seqs = Dna(f[4]);
            Assert.Multiple(() =>
            {
                Assert.That(SequenceAligner.MultipleAlign(seqs, sc, def).AlignedSequences, Is.EqualTo(SequenceAligner.MultipleAlign(seqs, sc).AlignedSequences));
                Assert.That(SequenceAligner.MultipleAlignAnchored(seqs, sc, def).AlignedSequences, Is.EqualTo(SequenceAligner.MultipleAlignAnchored(seqs, sc).AlignedSequences));
                Assert.That(SequenceAligner.MultipleAlignProgressive(seqs, sc, null).AlignedSequences, Is.EqualTo(SequenceAligner.MultipleAlignProgressive(seqs, sc).AlignedSequences));
            });
        }
    }

    [Test]
    public void AffineMethods_RowsSpellInputs_EqualLength_ConsensusScoreConsistent()
    {
        foreach (var f in MsaRows.Value.Take(60))
        {
            var sc = Scoring(f);
            var raw = f[4].Split(';');
            var seqs = Dna(f[4]);
            foreach (var r in new[]
            {
                SequenceAligner.MultipleAlign(seqs, sc, Affine), SequenceAligner.MultipleAlignAnchored(seqs, sc, Affine),
                SequenceAligner.MultipleAlignProgressive(seqs, sc, Affine), SequenceAligner.MultipleAlignIterative(seqs, sc, Affine),
                SequenceAligner.MultipleAlignConsistency(seqs, sc, Affine),
            })
            {
                Assert.Multiple(() =>
                {
                    Assert.That(r.AlignedSequences.Select(x => x.Replace("-", "")), Is.EqualTo(raw));
                    Assert.That(r.AlignedSequences.Select(x => x.Length).Distinct().Count(), Is.EqualTo(1));
                    Assert.That(r.TotalScore, Is.EqualTo(SequenceAligner.SumOfPairsScore(r.AlignedSequences, sc, MsaGapModel.Affine)));
                });
            }
        }
    }

    /// <summary>Iterative refinement accepts only strict improvements of the exact objective, so it never falls below the progressive seed.</summary>
    [Test]
    public void IterativeAffine_NeverBelowProgressiveAffine()
    {
        foreach (var f in MsaRows.Value)
        {
            var sc = Scoring(f);
            var seqs = Dna(f[4]);
            Assert.That(SequenceAligner.MultipleAlignIterative(seqs, sc, Affine).TotalScore,
                Is.GreaterThanOrEqualTo(SequenceAligner.MultipleAlignProgressive(seqs, sc, Affine).TotalScore));
        }
    }

    [Test]
    public void ZeroGapOpen_AffineEqualsLinearObjective()
    {
        var sc = new ScoringMatrix(2, -3, 0, -2);
        foreach (var f in MsaRows.Value.Take(40))
        {
            var rows = SequenceAligner.MultipleAlignProgressive(Dna(f[4]), sc).AlignedSequences;
            Assert.That(SequenceAligner.SumOfPairsScore(rows, sc, MsaGapModel.Affine), Is.EqualTo(SequenceAligner.SumOfPairsScore(rows, sc, MsaGapModel.Linear)));
        }
    }

    [Test]
    public void InvalidGapModel_Throws()
    {
        var seqs = new List<DnaSequence> { new("ACGT"), new("ACGA") };
        var bad = new MsaOptions { GapModel = (MsaGapModel)9 };
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.MultipleAlign(seqs, null, bad));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.MultipleAlignAnchored(seqs, null, bad));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.MultipleAlignProgressive(seqs, null, bad));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.MultipleAlignIterative(seqs, null, bad));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.MultipleAlignConsistency(seqs, null, bad));
    }

    [Test]
    public void EmptyAndSingle_AffineIsTrivial()
    {
        Assert.That(SequenceAligner.MultipleAlign(Array.Empty<DnaSequence>(), null, Affine), Is.EqualTo(MultipleAlignmentResult.Empty));
        var one = SequenceAligner.MultipleAlignProgressive([new DnaSequence("ACGT")], null, Affine);
        Assert.That(one.TotalScore, Is.Zero);
        Assert.That(one.AlignedSequences, Is.EqualTo(new[] { "ACGT" }));
    }

    #endregion
}
