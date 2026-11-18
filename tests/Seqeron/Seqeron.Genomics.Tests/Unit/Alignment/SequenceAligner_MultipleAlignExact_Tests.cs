namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A1-11, B09 F60): EXACT sum-of-pairs multiple alignment <see cref="SequenceAligner.MultipleAlignExact"/>
/// (k-dimensional dynamic program, natural affine gap states, Carrillo–Lipman 1988 bound).
///
/// Oracle (TestData/Alignment/msa_exact_oracle.tsv, scripts/alignment/generate_msa_exact_oracle.py): the exact optimum of 804 random sets
/// (k = 2..5, linear and natural affine objective) by two INDEPENDENT Python programs (lattice DP with 2^k - 1 neighbours; DP over the
/// position vector and the run state of every pair), both asserted equal to the exhaustive enumeration of every alignment of the tiny sets and
/// (k = 2) to Biopython 1.88 PairwiseAligner. The C# program must equal that optimum with and without the Carrillo–Lipman bound, return rows
/// that score to it under <see cref="SequenceAligner.SumOfPairsScore"/>, and equal GlobalAlign / GlobalAlignAffine for k = 2. The heuristic aligners
/// never exceed it; how often each reaches it is recorded (Evidence file).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MultipleAlignExact_Tests
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

    private sealed record Row(MsaGapModel Model, ScoringMatrix Scoring, string[] Seqs, int Optimum, string[] OptimalRows)
    {
        public List<DnaSequence> Dna() => Seqs.Select(s => new DnaSequence(s)).ToList();
        public MsaOptions Options(bool prune = true) => new() { GapModel = Model, ExactCarrilloLipman = prune };
        public override string ToString() => $"{Model} {Scoring.Match}/{Scoring.Mismatch}/{Scoring.GapOpen}/{Scoring.GapExtend} {string.Join(';', Seqs)}";
    }

    private static readonly Lazy<List<Row>> Oracle = new(() => ReadTsv("msa_exact_oracle.tsv").Select(f => new Row(
        f[0] == "affine" ? MsaGapModel.Affine : MsaGapModel.Linear,
        new ScoringMatrix(int.Parse(f[1]), int.Parse(f[2]), int.Parse(f[3]), int.Parse(f[4])),
        f[5].Split(';'),
        int.Parse(f[6]),
        f[7].Split(';'))).ToList());

    // ------------------------------------------------------------------------------------------------ oracle

    [Test]
    public void Oracle_OptimalRows_ScoreToTheOptimum_UnderTheObjective()
    {
        Assert.That(Oracle.Value.Count, Is.GreaterThanOrEqualTo(800));
        foreach (var r in Oracle.Value)
            Assert.That(SequenceAligner.SumOfPairsScore(r.OptimalRows, r.Scoring, r.Model), Is.EqualTo(r.Optimum), r.ToString());
    }

    [Test]
    public void ExactOptimum_EqualsTheIndependentOracle_WithAndWithoutCarrilloLipman()
    {
        foreach (var r in Oracle.Value)
        {
            foreach (bool prune in new[] { true, false })
            {
                var res = SequenceAligner.MultipleAlignExactWithStatistics(r.Dna(), r.Scoring, r.Options(prune));
                var a = res.Alignment;
                Assert.Multiple(() =>
                {
                    Assert.That(a.TotalScore, Is.EqualTo(r.Optimum), $"{r} prune={prune}");
                    Assert.That(SequenceAligner.SumOfPairsScore(a.AlignedSequences, r.Scoring, r.Model), Is.EqualTo(r.Optimum), $"{r} rows score");
                    Assert.That(a.AlignedSequences.Select(x => x.Replace("-", "")), Is.EqualTo(r.Seqs), $"{r} rows are the input sequences");
                    Assert.That(a.AlignedSequences.Select(x => x.Length).Distinct().Count(), Is.EqualTo(1), r.ToString());
                    Assert.That(res.VisitedCells, Is.LessThanOrEqualTo(res.LatticeCells), r.ToString());
                    if (!prune || r.Seqs.Length == 2)
                        Assert.That(res.VisitedCells, Is.EqualTo(res.LatticeCells).Or.Not.EqualTo(0L), r.ToString());
                });
            }
        }
    }

    [Test]
    public void TwoSequences_AreGlobalAlign_AndGlobalAlignAffine_RowForRow()
    {
        int n = 0;
        foreach (var r in Oracle.Value.Where(x => x.Seqs.Length == 2 && x.Seqs.All(s => s.Length > 0)))
        {
            var seqs = r.Dna();
            var single = r.Model == MsaGapModel.Affine ? SequenceAligner.GlobalAlignAffine(seqs[0], seqs[1], r.Scoring) : SequenceAligner.GlobalAlign(seqs[0], seqs[1], r.Scoring);
            var exact = SequenceAligner.MultipleAlignExact(seqs, r.Scoring, r.Options());
            Assert.That(exact.AlignedSequences, Is.EqualTo(new[] { single.AlignedSequence1, single.AlignedSequence2 }), r.ToString());
            Assert.That(exact.TotalScore, Is.EqualTo(single.Score), r.ToString());
            n++;
        }

        Assert.That(n, Is.GreaterThanOrEqualTo(40));
    }

    /// <summary>Honest optimality statement: no heuristic exceeds the exact optimum; the share of sets each reaches it is recorded.</summary>
    [Test]
    public void Heuristics_NeverExceedTheExactOptimum_AndTheirOptimalShareIsRecorded()
    {
        foreach (var model in new[] { MsaGapModel.Linear, MsaGapModel.Affine })
        {
            var sets = Oracle.Value.Where(x => x.Model == model && x.Seqs.Length >= 3 && x.Seqs.Length <= 4 && x.Seqs.Sum(s => s.Length) <= 30).ToList();
            int star = 0, prog = 0, iter = 0, cons = 0, anch = 0, any = 0;
            foreach (var r in sets)
            {
                var seqs = r.Dna();
                var opt = new MsaOptions { GapModel = model };
                int[] scores =
                [
                    SequenceAligner.MultipleAlign(seqs, r.Scoring, opt).TotalScore,
                    SequenceAligner.MultipleAlignProgressive(seqs, r.Scoring, opt).TotalScore,
                    SequenceAligner.MultipleAlignIterative(seqs, r.Scoring, opt).TotalScore,
                    SequenceAligner.MultipleAlignConsistency(seqs, r.Scoring, opt).TotalScore,
                    SequenceAligner.MultipleAlignAnchored(seqs, r.Scoring, opt).TotalScore,
                ];
                Assert.That(scores, Is.All.LessThanOrEqualTo(r.Optimum), r.ToString());
                if (scores[0] == r.Optimum) star++;
                if (scores[1] == r.Optimum) prog++;
                if (scores[2] == r.Optimum) iter++;
                if (scores[3] == r.Optimum) cons++;
                if (scores[4] == r.Optimum) anch++;
                if (scores.Max() == r.Optimum) any++;
            }

            TestContext.Out.WriteLine($"{model}: {sets.Count} sets (k = 3..4, total length <= 30); optimal: star {star}, progressive {prog}, iterative {iter}, consistency {cons}, anchored {anch}, best-of-five {any}");
            Assert.That(sets.Count, Is.GreaterThan(100));
            Assert.That(iter, Is.GreaterThanOrEqualTo(prog), "refinement never lowers the objective");
            Assert.That(star + prog + iter, Is.GreaterThan(sets.Count / 2), "the heuristics are optimal on a substantial share of tiny sets");
        }
    }

    // ------------------------------------------------------------------------------------------------ brute force

    /// <summary>Exhaustive enumeration of every column sequence (written here, independent of the program under test).</summary>
    private static int BruteForce(string[] seqs, ScoringMatrix sc, MsaGapModel model)
    {
        int k = seqs.Length;
        var moves = Enumerable.Range(1, (1 << k) - 1).ToArray();
        int best = int.MinValue;
        var cols = new List<int>();
        void Rec(int[] pos)
        {
            if (Enumerable.Range(0, k).All(i => pos[i] == seqs[i].Length))
            {
                var rows = Enumerable.Range(0, k).Select(i =>
                {
                    var p = 0;
                    return new string(cols.Select(m => (m >> i & 1) != 0 ? seqs[i][p++] : '-').ToArray());
                }).ToArray();
                best = Math.Max(best, SequenceAligner.SumOfPairsScore(rows, sc, model));
                return;
            }

            foreach (int m in moves)
            {
                bool ok = true;
                for (int i = 0; i < k; i++)
                    if ((m >> i & 1) != 0 && pos[i] >= seqs[i].Length) ok = false;
                if (!ok)
                    continue;
                var np = (int[])pos.Clone();
                for (int i = 0; i < k; i++) np[i] += m >> i & 1;
                cols.Add(m);
                Rec(np);
                cols.RemoveAt(cols.Count - 1);
            }
        }

        Rec(new int[k]);
        return best;
    }

    [Test]
    public void ExactOptimum_EqualsExhaustiveEnumeration_OfTinySets_LinearAndAffine()
    {
        var rnd = new Random(11);
        var schemes = new[] { new ScoringMatrix(1, -1, -2, -1), new ScoringMatrix(2, -3, -5, -2), new ScoringMatrix(5, -4, -10, -1), new ScoringMatrix(1, -2, 0, -3), new ScoringMatrix(1, 0, -1, 0) };
        int checkedSets = 0;
        for (int t = 0; t < 60; t++)
        {
            int k = 3 + t % 2;
            int maxLen = k == 3 ? 4 : 3;
            var seqs = Enumerable.Range(0, k).Select(_ => new string(Enumerable.Range(0, rnd.Next(0, maxLen + 1)).Select(__ => "ACGT"[rnd.Next(rnd.Next(2) == 0 ? 2 : 4)]).ToArray())).ToArray();
            var sc = schemes[rnd.Next(schemes.Length)];
            foreach (var model in new[] { MsaGapModel.Linear, MsaGapModel.Affine })
            {
                int brute = BruteForce(seqs, sc, model);
                var dna = seqs.Select(s => new DnaSequence(s)).ToList();
                foreach (bool prune in new[] { true, false })
                    Assert.That(SequenceAligner.MultipleAlignExact(dna, sc, new MsaOptions { GapModel = model, ExactCarrilloLipman = prune }).TotalScore,
                        Is.EqualTo(brute), $"{model} prune={prune} {string.Join(';', seqs)} {sc}");
                checkedSets++;
            }
        }

        Assert.That(checkedSets, Is.EqualTo(120));
    }

    // ------------------------------------------------------------------------------------------------ pruning, guard, edge cases

    [Test]
    public void CarrilloLipman_PrunesSimilarSequences_AndKeepsTheOptimum()
    {
        string core = "ACGTTGCAAGCTTCGATCGGATCCATGCAGTCAAGTTGCA";   // 40 nt
        var seqs = new[] { core, core.Remove(10, 1), core.Insert(25, "G"), core.Remove(30, 2) }.Select(s => new DnaSequence(s)).ToList();
        var sc = new ScoringMatrix(2, -3, -5, -2);
        foreach (var model in new[] { MsaGapModel.Linear, MsaGapModel.Affine })
        {
            var pruned = SequenceAligner.MultipleAlignExactWithStatistics(seqs, sc, new MsaOptions { GapModel = model });
            Assert.That(pruned.LatticeCells, Is.GreaterThan(2_000_000));
            Assert.That(pruned.VisitedCells, Is.LessThan(pruned.LatticeCells / 100), $"{model}: the bound must prune a similar set to a thin corridor");
            Assert.That(pruned.LowerBound, Is.LessThanOrEqualTo(pruned.Alignment.TotalScore));
            Assert.That(SequenceAligner.SumOfPairsScore(pruned.Alignment.AlignedSequences, sc, model), Is.EqualTo(pruned.Alignment.TotalScore));
            foreach (var h in new[]
            {
                SequenceAligner.MultipleAlign(seqs, sc, new MsaOptions { GapModel = model }),
                SequenceAligner.MultipleAlignProgressive(seqs, sc, new MsaOptions { GapModel = model }),
            })
                Assert.That(h.TotalScore, Is.LessThanOrEqualTo(pruned.Alignment.TotalScore));
        }
    }

    [Test]
    public void CarrilloLipman_AgreesWithTheFullLattice_OnMediumSets()
    {
        // 3 x ~40 nt (~70 000 cells) and 4 x ~12 nt (~30 000 cells): full lattice vs pruned, both objectives
        var rnd = new Random(5);
        var sc = new ScoringMatrix(1, -1, -3, -1);
        for (int t = 0; t < 12; t++)
        {
            int k = 3 + t % 2, len = k == 3 ? 40 : 12;
            string baseSeq = new(Enumerable.Range(0, len).Select(_ => "ACGT"[rnd.Next(4)]).ToArray());
            var seqs = Enumerable.Range(0, k).Select(_ => new DnaSequence(Mutate(baseSeq, rnd))).ToList();
            foreach (var model in new[] { MsaGapModel.Linear, MsaGapModel.Affine })
            {
                var full = SequenceAligner.MultipleAlignExactWithStatistics(seqs, sc, new MsaOptions { GapModel = model, ExactCarrilloLipman = false });
                var pruned = SequenceAligner.MultipleAlignExactWithStatistics(seqs, sc, new MsaOptions { GapModel = model });
                Assert.That(pruned.Alignment.TotalScore, Is.EqualTo(full.Alignment.TotalScore), $"{model} {string.Join(';', seqs.Select(s => s.Sequence))}");
                Assert.That(full.VisitedCells, Is.EqualTo(full.LatticeCells));
                Assert.That(pruned.VisitedCells, Is.LessThanOrEqualTo(full.VisitedCells));
            }
        }
    }

    private static string Mutate(string s, Random rnd)
    {
        var sb = new System.Text.StringBuilder();
        foreach (char c in s)
        {
            double r = rnd.NextDouble();
            if (r < 0.1)
                sb.Append("ACGT"[rnd.Next(4)]);      // substitution
            else if (r >= 0.17)                      // 0.1 <= r < 0.17: deletion (nothing appended)
            {
                sb.Append(c);
                if (r < 0.24)
                    sb.Append("ACGT"[rnd.Next(4)]);  // insertion after the residue
            }
        }

        return sb.Length == 0 ? "A" : sb.ToString();
    }

    [Test]
    public void IdenticalSequences_ScoreIsPairCountTimesLengthTimesMatch()
    {
        string s = "ACGTACGGTCAGTCAGGTCAAGCTAGCTAGGATCGATCGGCTAGCATGCATCGA";
        var seqs = Enumerable.Range(0, 5).Select(_ => new DnaSequence(s)).ToList();
        var sc = new ScoringMatrix(3, -2, -4, -1);
        var res = SequenceAligner.MultipleAlignExactWithStatistics(seqs, sc, new MsaOptions { GapModel = MsaGapModel.Linear });
        Assert.That(res.Alignment.TotalScore, Is.EqualTo(10 * s.Length * 3));
        Assert.That(res.Alignment.AlignedSequences, Is.All.EqualTo(s));
        Assert.That(res.VisitedCells, Is.EqualTo(s.Length + 1), "a single diagonal corridor");
    }

    [Test]
    public void ComplexityGuard_ThrowsClearException_InsteadOfExhaustingMemory()
    {
        var seqs = new[] { "ACGTACGTACGTACGT", "TTGACCATTGACCATT", "GGCATTACGGCATTAC" }.Select(s => new DnaSequence(s)).ToList();
        var ex = Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignExact(seqs, null, new MsaOptions { ExactCarrilloLipman = false, ExactMaxStateCells = 1000 }));
        Assert.That(ex!.Message, Does.Contain("complexity limit"));
        Assert.That(ex.Message, Does.Contain("ExactMaxStateCells"));
        Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignExact(seqs, null, new MsaOptions { GapModel = MsaGapModel.Affine, ExactCarrilloLipman = false, ExactMaxStateCells = 5000 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.MultipleAlignExact(seqs, null, new MsaOptions { ExactMaxStateCells = 0 }));

        // the same set solves within the default budget
        Assert.DoesNotThrow(() => SequenceAligner.MultipleAlignExact(seqs));
    }

    [Test]
    public void TooManySequences_AndNullAndTrivialInputs()
    {
        var many = Enumerable.Range(0, 9).Select(_ => new DnaSequence("ACGT")).ToList();
        Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignExact(many));
        Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignExact(many.Take(7), null, new MsaOptions { GapModel = MsaGapModel.Affine }));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignExact(null!));
        Assert.That(SequenceAligner.MultipleAlignExact(Array.Empty<DnaSequence>()), Is.EqualTo(MultipleAlignmentResult.Empty));
        var one = SequenceAligner.MultipleAlignExact(new[] { new DnaSequence("ACGT") });
        Assert.That(one.AlignedSequences, Is.EqualTo(new[] { "ACGT" }));
        Assert.That(one.TotalScore, Is.Zero);
        // empty sequences are legal and align against gaps
        var withEmpty = SequenceAligner.MultipleAlignExact(new[] { new DnaSequence("ACG"), new DnaSequence("ACG"), new DnaSequence("") });
        Assert.That(withEmpty.AlignedSequences, Is.EqualTo(new[] { "ACG", "ACG", "---" }));
    }

    [Test]
    public void Result_IsDeterministic_AndExactMsaResultCarriesTheBound()
    {
        var seqs = new[] { "ACGTTGCA", "ACGTGCA", "ACTTGCAA", "ACGTTCA" }.Select(s => new DnaSequence(s)).ToList();
        var sc = new ScoringMatrix(2, -3, -5, -2);
        var a = SequenceAligner.MultipleAlignExactWithStatistics(seqs, sc, new MsaOptions { GapModel = MsaGapModel.Affine });
        var b = SequenceAligner.MultipleAlignExactWithStatistics(seqs, sc, new MsaOptions { GapModel = MsaGapModel.Affine });
        Assert.That(a.Alignment.AlignedSequences, Is.EqualTo(b.Alignment.AlignedSequences));
        Assert.That(a.LowerBound, Is.LessThanOrEqualTo(a.Alignment.TotalScore));
        Assert.That(a.StateValues, Is.EqualTo(a.VisitedCells * 75), "k = 4: 75 weak orders of the last-move times");
        var lin = SequenceAligner.MultipleAlignExactWithStatistics(seqs, sc);
        Assert.That(lin.StateValues, Is.EqualTo(lin.VisitedCells));
        Assert.That(a.Alignment.Consensus.Length, Is.EqualTo(a.Alignment.AlignedSequences[0].Length));
    }
}
