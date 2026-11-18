using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A1-14, B09 F61): protein / substitution-matrix multiple sequence alignment — every multiple aligner
/// (<c>MultipleAlign</c> center star, <c>MultipleAlignAnchored</c>, <c>MultipleAlignProgressive</c>, <c>MultipleAlignIterative</c>,
/// <c>MultipleAlignConsistency</c>, <c>MultipleAlignExact</c>) accepts <see cref="ProteinSequence"/> / string sequences with a
/// <see cref="MatrixScoring"/> (BLOSUM62, PAM250, EDNAFULL, ...; linear or affine, possibly fractional gap scores) and shares the dynamic programs,
/// merge, profile and objective code of the DNA overloads (no second implementation).
///
/// Oracles (TestData/Alignment/msa_matrix_*.tsv, scripts/alignment/generate_msa_matrix_oracle.py):
///  * <c>msa_matrix_sp_oracle.tsv</c>: independent Python re-score (exact rationals) of 420 random gapped protein / IUPAC-DNA MSAs, with the
///    Biopython 1.88 PairwiseAligner(substitution_matrix=...) optimum of every pair as an upper bound of the sum of pairs;
///  * <c>msa_matrix_exact_oracle.tsv</c>: the exact optimum of 308 small sets (k = 2: Biopython optimum; k = 3, 4: dynamic program over the
///    k-dimensional lattice, validated against exhaustive enumeration of every alignment for the tiniest sets);
///  * <c>msa_matrix_oracle.tsv</c>: independent Python implementations (written from the algorithm descriptions) of the star, anchored, progressive,
///    iterative and consistency aligners under linear / affine, integer / fractional gap scores: every row must be equal (238 sets, k = 2..6, protein
///    BLOSUM62 / PAM250 and IUPAC DNA EDNAFULL); k = 2 results are asserted equal to the Biopython optimum there.
/// Plus: the matrix path with a +1 / -1 substitution matrix reproduces the scalar <see cref="ScoringMatrix"/> aligners row for row.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MultipleAlignMatrix_Tests
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

    private static readonly Lazy<List<string[]>> SpRows = new(() => ReadTsv("msa_matrix_sp_oracle.tsv").ToList());
    private static readonly Lazy<List<string[]>> ExactRows = new(() => ReadTsv("msa_matrix_exact_oracle.tsv").ToList());
    private static readonly Lazy<List<string[]>> MsaRows = new(() => ReadTsv("msa_matrix_oracle.tsv").ToList());

    private static double D(string s) => double.Parse(s, CultureInfo.InvariantCulture);

    private static SubstitutionMatrix Matrix(string name) => name switch
    {
        "BLOSUM62" => SubstitutionMatrix.Blosum62,
        "PAM250" => SubstitutionMatrix.Pam250,
        "EDNAFULL" => SubstitutionMatrix.EdnaFull,
        _ => throw new ArgumentException(name),
    };

    private static MatrixScoring Scoring(string matrix, string open, string ext) => new(Matrix(matrix), D(open), D(ext));

    private static MsaOptions Options(string model) => new() { GapModel = model == "affine" ? MsaGapModel.Affine : MsaGapModel.Linear };

    private static MsaGapModel Model(string model) => model == "affine" ? MsaGapModel.Affine : MsaGapModel.Linear;

    private const double Eps = 1e-9;

    #region Sum-of-pairs objective

    [Test]
    public void MatrixSumOfPairsScore_MatchesIndependentRescore_AndRespectsPairwiseOptimumBound()
    {
        var rows = SpRows.Value;
        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(420));
        int bad = 0;
        string? first = null;
        foreach (var f in rows)
        {
            var aligned = f[3].Split(';');
            var sc = Scoring(f[0], f[1], f[2]);
            double linear = SequenceAligner.MatrixSumOfPairsScore(aligned, sc, MsaGapModel.Linear);
            double affine = SequenceAligner.MatrixSumOfPairsScore(aligned, sc, MsaGapModel.Affine);
            bool ok = Math.Abs(linear - D(f[4])) < Eps && Math.Abs(affine - D(f[5])) < Eps
                   && linear <= D(f[6]) + Eps && affine <= D(f[7]) + Eps;
            if (!ok && ++bad == 1) first = $"{f[3]} {f[0]} {f[1]} {f[2]}: {linear}/{affine} vs {f[4]}/{f[5]}";
        }

        Assert.That(bad, Is.Zero, first);
    }

    /// <summary>Hand example, BLOSUM62: M/M = 5, K/K = 5, V/V = 4, L/L = 4; the gap column of the pair costs the gap score.</summary>
    [Test]
    public void MatrixSumOfPairsScore_HandDerivedExample()
    {
        string[] rows = ["MK-V", "MKLV"];
        // Linear gap -4: 5 + 5 - 4 + 4 = 10. Affine open -10, extend -1: one run of length 1 = -11: 5 + 5 - 11 + 4 = 3.
        Assert.That(SequenceAligner.MatrixSumOfPairsScore(rows, new MatrixScoring(SubstitutionMatrix.Blosum62, -10, -4), MsaGapModel.Linear), Is.EqualTo(10));
        Assert.That(SequenceAligner.MatrixSumOfPairsScore(rows, new MatrixScoring(SubstitutionMatrix.Blosum62, -10, -1), MsaGapModel.Affine), Is.EqualTo(3));
        // EMBOSS needle numbers gapopen 10 / gapextend 0.5 (a gap of length 1 costs 10): -10 instead of -11 -> 5 + 5 - 10 + 4 = 4; fractional extension 0.5 for a run of 2.
        Assert.That(SequenceAligner.MatrixSumOfPairsScore(rows, MatrixScoring.Emboss(SubstitutionMatrix.Blosum62, 10, 0.5), MsaGapModel.Affine), Is.EqualTo(4));
        Assert.That(SequenceAligner.MatrixSumOfPairsScore(["MK--V", "MKLLV"], MatrixScoring.Emboss(SubstitutionMatrix.Blosum62, 10, 0.5), MsaGapModel.Affine),
            Is.EqualTo(5 + 5 - 10.5 + 4));
        // Lower-case rows are compared case-insensitively; gap-gap columns are neutral; one row has nothing to score.
        Assert.That(SequenceAligner.MatrixSumOfPairsScore(["mk-v", "mklv"], new MatrixScoring(SubstitutionMatrix.Blosum62, -10, -4)), Is.EqualTo(10));
        Assert.That(SequenceAligner.MatrixSumOfPairsScore(["MK-V", "MK-V"], new MatrixScoring(SubstitutionMatrix.Blosum62, -10, -4)), Is.EqualTo(14));
        Assert.That(SequenceAligner.MatrixSumOfPairsScore(["MKV"], new MatrixScoring(SubstitutionMatrix.Blosum62, -10, -4)), Is.Zero);
    }

    [Test]
    public void MatrixSumOfPairsScore_Validation()
    {
        var sc = new MatrixScoring(SubstitutionMatrix.Blosum62, -10, -1);
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => SequenceAligner.MatrixSumOfPairsScore(null!, sc));
            Assert.Throws<ArgumentNullException>(() => SequenceAligner.MatrixSumOfPairsScore(["MK"], null!));
            Assert.Throws<ArgumentException>(() => SequenceAligner.MatrixSumOfPairsScore(["MK", "MKV"], sc), "different lengths");
            Assert.Throws<ArgumentException>(() => SequenceAligner.MatrixSumOfPairsScore(["MK", "M1"], sc), "residue outside the matrix alphabet");
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.MatrixSumOfPairsScore(["MK", "MK"], sc, (MsaGapModel)9));
        });
    }

    #endregion

    #region Exact optimum

    [Test]
    public void MultipleAlignExact_EqualsTheOracleOptimum_AndItsRowsAchieveIt()
    {
        var rows = ExactRows.Value;
        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(300));
        Assert.That(rows.Select(f => f[4].Split(';').Length).Distinct().Order(), Is.EqualTo(new[] { 2, 3, 4 }));
        var failures = new List<string>();
        foreach (var f in rows)
        {
            var seqs = f[4].Split(';');
            var sc = Scoring(f[1], f[2], f[3]);
            var r = SequenceAligner.MultipleAlignExact(seqs, sc, Options(f[0]));
            double optimum = D(f[5]);
            double check = SequenceAligner.MatrixSumOfPairsScore(r.AlignedSequences, sc, Model(f[0]));
            bool ok = Math.Abs(r.ExactTotalScore - optimum) < Eps && Math.Abs(check - optimum) < Eps
                   && r.TotalScore == (int)Math.Round(optimum, MidpointRounding.AwayFromZero)
                   && r.AlignedSequences.Select(s => s.Replace("-", "")).SequenceEqual(seqs.Select(s => s.ToUpperInvariant()));
            if (!ok) failures.Add($"{f[0]} {f[1]} {f[2]} {f[3]} {f[4]}: got {r.ExactTotalScore}/{check}, expected {optimum}");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count}; first: {failures.FirstOrDefault()}");
    }

    [Test]
    public void MultipleAlignExact_Rows_AreTheOracleRows()
    {
        // The optimal ROWS are not unique in general; for these sets the deterministic result is checked only through the optimum above.
        // Here: k = 2 is exactly the pairwise global alignment under the same MatrixScoring (rows included).
        foreach (var f in ExactRows.Value.Where(f => f[4].Split(';').Length == 2))
        {
            var seqs = f[4].Split(';');
            var sc = Scoring(f[1], f[2], f[3]);
            var r = SequenceAligner.MultipleAlignExact(seqs, sc, Options(f[0]));
            var pair = f[0] == "affine" ? SequenceAligner.GlobalAlignAffine(seqs[0], seqs[1], sc) : SequenceAligner.GlobalAlign(seqs[0], seqs[1], sc);
            Assert.That(r.AlignedSequences, Is.EqualTo(new[] { pair.AlignedSequence1, pair.AlignedSequence2 }), $"{f[0]} {f[4]}");
            Assert.That(r.ExactTotalScore, Is.EqualTo(pair.ExactScore).Within(Eps));
        }
    }

    [Test]
    public void MultipleAlignExact_PruningDoesNotChangeTheOptimum_AndStatisticsAreScaled()
    {
        var sc = MatrixScoring.Emboss(SubstitutionMatrix.Blosum62, 10, 0.5);
        string[] seqs = ["MKVLAAGIVG", "MKVLSAGIVG", "MKVAAGLVG"];
        var pruned = SequenceAligner.MultipleAlignExactWithStatistics(seqs, sc, new MsaOptions { GapModel = MsaGapModel.Affine });
        var full = SequenceAligner.MultipleAlignExactWithStatistics(seqs, sc, new MsaOptions { GapModel = MsaGapModel.Affine, ExactCarrilloLipman = false });
        Assert.Multiple(() =>
        {
            Assert.That(pruned.Alignment.ExactTotalScore, Is.EqualTo(full.Alignment.ExactTotalScore).Within(Eps));
            Assert.That(pruned.ScoreScale, Is.EqualTo(2), "gap scores 9.5 / 0.5 -> halves");
            Assert.That(pruned.VisitedCells, Is.LessThanOrEqualTo(pruned.LatticeCells));
            Assert.That(pruned.LowerBound / (double)pruned.ScoreScale, Is.LessThanOrEqualTo(pruned.Alignment.ExactTotalScore + Eps));
        });
    }

    #endregion

    #region Heuristic aligners: independent Python rows

    private static readonly string[] Methods = ["star", "anchored", "progressive", "iterative", "consistency"];

    private static MultipleAlignmentResult Run(string method, string[] seqs, bool protein, MatrixScoring sc, MsaOptions opt)
    {
        if (protein)
        {
            var p = seqs.Select(s => new ProteinSequence(s)).ToList();
            return method switch
            {
                "star" => SequenceAligner.MultipleAlign(p, sc, opt),
                "anchored" => SequenceAligner.MultipleAlignAnchored(p, sc, opt),
                "progressive" => SequenceAligner.MultipleAlignProgressive(p, sc, opt),
                "iterative" => SequenceAligner.MultipleAlignIterative(p, sc, opt),
                _ => SequenceAligner.MultipleAlignConsistency(p, sc, opt),
            };
        }

        return method switch
        {
            "star" => SequenceAligner.MultipleAlign(seqs, sc, opt),
            "anchored" => SequenceAligner.MultipleAlignAnchored(seqs, sc, opt),
            "progressive" => SequenceAligner.MultipleAlignProgressive(seqs, sc, opt),
            "iterative" => SequenceAligner.MultipleAlignIterative(seqs, sc, opt),
            _ => SequenceAligner.MultipleAlignConsistency(seqs, sc, opt),
        };
    }

    [Test]
    public void Aligners_EveryRowEqualsTheIndependentPythonImplementation()
    {
        var rows = MsaRows.Value;
        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(238));
        Assert.That(rows.Select(f => f[0]).Distinct().Order(), Is.EqualTo(new[] { "affine", "linear" }));
        Assert.That(rows.Select(f => f[1]).Distinct().Order(), Is.EqualTo(new[] { "BLOSUM62", "EDNAFULL", "PAM250" }));
        var failures = new List<string>();
        int checkedRows = 0;
        foreach (var f in rows)
        {
            var seqs = f[4].Split(';');
            var sc = Scoring(f[1], f[2], f[3]);
            var opt = Options(f[0]);
            bool protein = f[1] != "EDNAFULL";
            for (int m = 0; m < Methods.Length; m++)
            {
                string expected = f[5 + m];
                if (expected.Length == 0)
                    continue; // the anchored oracle is skipped for sets with a repeated anchor text (arbitrary text occurrence in the suffix tree)
                var r = Run(Methods[m], seqs, protein, sc, opt);
                checkedRows++;
                if (string.Join(';', r.AlignedSequences) != expected)
                    failures.Add($"{Methods[m]} {f[0]} {f[1]} {f[2]} {f[3]} {f[4]}: got {string.Join(';', r.AlignedSequences)} expected {expected}");
                else if (Math.Abs(r.ExactTotalScore - SequenceAligner.MatrixSumOfPairsScore(r.AlignedSequences, sc, Model(f[0]))) > Eps)
                    failures.Add($"{Methods[m]} TotalScore {f[4]}");
            }
        }

        Assert.That(checkedRows, Is.GreaterThanOrEqualTo(1150));
        Assert.That(failures, Is.Empty, $"{failures.Count} of {checkedRows}; first: {failures.FirstOrDefault()}");
    }

    [Test]
    public void TwoSequences_EveryMethodIsAnOptimalPairwiseAlignment()
    {
        foreach (var f in MsaRows.Value.Where(f => f[4].Split(';').Length == 2))
        {
            var seqs = f[4].Split(';');
            var sc = Scoring(f[1], f[2], f[3]);
            var opt = Options(f[0]);
            var pair = f[0] == "affine" ? SequenceAligner.GlobalAlignAffine(seqs[0], seqs[1], sc) : SequenceAligner.GlobalAlign(seqs[0], seqs[1], sc);
            foreach (var method in new[] { "star", "progressive", "iterative" })
            {
                var r = Run(method, seqs, f[1] != "EDNAFULL", sc, opt);
                Assert.That(r.ExactTotalScore, Is.EqualTo(pair.ExactScore).Within(Eps), $"{method} {f[0]} {f[1]} {f[4]}");
            }

            Assert.That(SequenceAligner.MultipleAlign(seqs, sc, opt).AlignedSequences, Is.EqualTo(new[] { pair.AlignedSequence1, pair.AlignedSequence2 }), "star rows = the pairwise alignment");
        }
    }

    [Test]
    public void ConsistencyEngineTCoffeeLinked_AcceptsProteinAndProducesAValidAlignment()
    {
        var f = MsaRows.Value.First(r => r[1] == "BLOSUM62" && r[4].Split(';').Length >= 4);
        var seqs = f[4].Split(';');
        var sc = Scoring(f[1], f[2], f[3]);
        var linked = SequenceAligner.MultipleAlignConsistency(seqs, sc, new MsaOptions { GapModel = Model(f[0]), ConsistencyEngine = MsaConsistencyEngine.TCoffeeLinked });
        Assert.Multiple(() =>
        {
            Assert.That(linked.AlignedSequences.Select(s => s.Replace("-", "")), Is.EqualTo(seqs));
            Assert.That(linked.AlignedSequences.Select(s => s.Length).Distinct().Count(), Is.EqualTo(1));
            Assert.That(linked.ExactTotalScore, Is.EqualTo(SequenceAligner.MatrixSumOfPairsScore(linked.AlignedSequences, sc, Model(f[0]))).Within(Eps));
        });
    }

    #endregion

    #region Scalar / matrix equivalence (the matrix path is the same code)

    /// <summary>A +1 / -1 matrix over ACGT with integer gaps is the scalar <see cref="ScoringMatrix"/>: every aligner returns the same rows.</summary>
    [Test]
    public void MatrixPath_WithMatchMismatchMatrix_ReproducesTheScalarAligners()
    {
        const string alphabet = "ACGT";
        var table = new int[4, 4];
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                table[i, j] = i == j ? 2 : -3;
        var matrix = new SubstitutionMatrix("M2-3", alphabet, table);
        var scalar = new ScoringMatrix(2, -3, -5, -2);
        var matrixScoring = new MatrixScoring(matrix, -5, -2);

        var rnd = new Random(14);
        for (int t = 0; t < 40; t++)
        {
            int k = rnd.Next(3, 6);
            string baseSeq = new(Enumerable.Range(0, rnd.Next(10, 30)).Select(_ => alphabet[rnd.Next(4)]).ToArray());
            var seqs = Enumerable.Range(0, k).Select(_ => Mutate(rnd, baseSeq)).ToArray();
            var dna = seqs.Select(s => new DnaSequence(s)).ToList();
            foreach (var model in new[] { MsaGapModel.Linear, MsaGapModel.Affine })
            {
                var opt = new MsaOptions { GapModel = model };
                Assert.Multiple(() =>
                {
                    AssertSame(SequenceAligner.MultipleAlign(dna, scalar, opt), SequenceAligner.MultipleAlign(seqs, matrixScoring, opt), "star");
                    AssertSame(SequenceAligner.MultipleAlignAnchored(dna, scalar, opt), SequenceAligner.MultipleAlignAnchored(seqs, matrixScoring, opt), "anchored");
                    AssertSame(SequenceAligner.MultipleAlignProgressive(dna, scalar, opt), SequenceAligner.MultipleAlignProgressive(seqs, matrixScoring, opt), "progressive");
                    AssertSame(SequenceAligner.MultipleAlignIterative(dna, scalar, opt), SequenceAligner.MultipleAlignIterative(seqs, matrixScoring, opt), "iterative");
                    AssertSame(SequenceAligner.MultipleAlignConsistency(dna, scalar, opt), SequenceAligner.MultipleAlignConsistency(seqs, matrixScoring, opt), "consistency");
                    if (k <= 4)
                        AssertSame(SequenceAligner.MultipleAlignExact(dna.Take(3), scalar, opt), SequenceAligner.MultipleAlignExact(seqs.Take(3), matrixScoring, opt), "exact");
                });
            }
        }
    }

    private static string Mutate(Random rnd, string s)
    {
        var sb = new System.Text.StringBuilder();
        foreach (char c in s)
        {
            double r = rnd.NextDouble();
            if (r < 0.08) sb.Append("ACGT"[rnd.Next(4)]);
            else if (r < 0.13) sb.Append("");
            else if (r < 0.18) { sb.Append(c); sb.Append("ACGT"[rnd.Next(4)]); }
            else sb.Append(c);
        }

        return sb.Length == 0 ? "A" : sb.ToString();
    }

    private static void AssertSame(MultipleAlignmentResult scalar, MultipleAlignmentResult matrix, string what)
    {
        Assert.That(matrix.AlignedSequences, Is.EqualTo(scalar.AlignedSequences), what);
        Assert.That(matrix.TotalScore, Is.EqualTo(scalar.TotalScore), what);
        Assert.That(matrix.Consensus, Is.EqualTo(scalar.Consensus), what);
    }

    #endregion

    #region Contract

    private static readonly string[] Globins =
    [
        "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF",
        "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKV",
        "MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDIAAKYKELGYQG",
    ];

    [Test]
    public void Overloads_ProteinAndStringAndCaseAreEquivalent_AndScoreIsReal()
    {
        var sc = MatrixScoring.Emboss(SubstitutionMatrix.Blosum62, 10, 0.5);
        var opt = new MsaOptions { GapModel = MsaGapModel.Affine };
        var protein = SequenceAligner.MultipleAlignProgressive(Globins.Select(g => new ProteinSequence(g)).ToList(), sc, opt);
        var text = SequenceAligner.MultipleAlignProgressive(Globins, sc, opt);
        var lower = SequenceAligner.MultipleAlignProgressive(Globins.Select(g => g.ToLowerInvariant()).ToList(), sc, opt);
        Assert.Multiple(() =>
        {
            Assert.That(text.AlignedSequences, Is.EqualTo(protein.AlignedSequences));
            Assert.That(lower.AlignedSequences, Is.EqualTo(protein.AlignedSequences), "lower-case input, upper-case rows");
            Assert.That(text.AlignedSequences.Select(s => s.Replace("-", "")), Is.EqualTo(Globins));
            Assert.That(text.ExactTotalScore, Is.EqualTo(SequenceAligner.MatrixSumOfPairsScore(text.AlignedSequences, sc, MsaGapModel.Affine)).Within(Eps));
            Assert.That(text.TotalScore, Is.EqualTo((int)Math.Round(text.ExactTotalScore, MidpointRounding.AwayFromZero)));
            Assert.That(text.Consensus.Length, Is.EqualTo(text.AlignedSequences[0].Length));
            Assert.That(text.Consensus.Replace("-", ""), Is.Not.Empty);
            // The same sequences through the scalar path are impossible (letters outside ACGT); the matrix path accepts them.
            Assert.That(text.AlignedSequences.Select(s => s.Length).Distinct().Count(), Is.EqualTo(1));
        });
    }

    [Test]
    public void ConsensusIsTheColumnMajority_TiesPreferResiduesThenSmallerLetter()
    {
        var sc = new MatrixScoring(SubstitutionMatrix.Blosum62, -10, -1);
        var r = SequenceAligner.MultipleAlign(["MKVL", "MKVL", "MRVL", "MKVL"], sc);
        Assert.That(r.AlignedSequences, Is.EqualTo(new[] { "MKVL", "MKVL", "MRVL", "MKVL" }));
        Assert.That(r.Consensus, Is.EqualTo("MKVL"));
        // two sequences with a different residue in one column: tie -> the smaller letter (K < R)
        var tie = SequenceAligner.MultipleAlign(["MKV", "MRV"], sc);
        Assert.That(tie.Consensus, Is.EqualTo("MKV"));
    }

    [Test]
    public void TrivialInputs_AreHandledLikeTheDnaOverloads()
    {
        var sc = new MatrixScoring(SubstitutionMatrix.Blosum62, -10, -1);
        foreach (var call in new Func<IEnumerable<string>, MultipleAlignmentResult>[]
        {
            s => SequenceAligner.MultipleAlign(s, sc),
            s => SequenceAligner.MultipleAlignAnchored(s, sc),
            s => SequenceAligner.MultipleAlignProgressive(s, sc),
            s => SequenceAligner.MultipleAlignIterative(s, sc),
            s => SequenceAligner.MultipleAlignConsistency(s, sc),
            s => SequenceAligner.MultipleAlignExact(s, sc),
        })
        {
            Assert.That(call([]).AlignedSequences, Is.Empty);
            var one = call(["mkv"]);
            Assert.That(one.AlignedSequences, Is.EqualTo(new[] { "MKV" }));
            Assert.That(one.TotalScore, Is.Zero);
            Assert.That(one.Consensus, Is.EqualTo("MKV"));
        }
    }

    [Test]
    public void InvalidInput_IsRejectedWithClearExceptions()
    {
        var sc = new MatrixScoring(SubstitutionMatrix.Blosum62, -10, -1);
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlign((IEnumerable<string>)null!, sc));
            Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlign((IEnumerable<ProteinSequence>)null!, sc));
            Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlign(["MK", "MV"], null!));
            Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlign(["MK", null!], sc));
            Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlign(["MK", "M-V"], sc), "'-' is the gap, not a residue");
            Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignProgressive(["MK", "M1V"], sc), "digit outside the alphabet");
            Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignExact(["MK", "MV"], new MatrixScoring(SubstitutionMatrix.Blosum62, 1, -1)), "positive gap open");
            Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignConsistency(["MK", "MV"], new MatrixScoring(SubstitutionMatrix.Blosum62, -1, double.NaN)));
            Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignIterative(["MK", "MV"], new MatrixScoring(SubstitutionMatrix.Blosum62, -1, -0.0000001)), "finer than 1e-6");
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.MultipleAlignIterative(["MK", "MV"], sc, null, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.MultipleAlign(["MK", "MV"], sc, new MsaOptions { GapModel = (MsaGapModel)5 }));
            Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignExact(Enumerable.Repeat("MKV", 9), sc), "more than 8 sequences");
        });
    }

    [Test]
    public void ScoreIsExactForFractionalGaps_AndTotalScoreIsRoundedHalfAwayFromZero()
    {
        // Two identical residues and one extra: with gap extension 0.5 the star alignment of ["MKV", "MV"] puts K against a gap.
        var sc = new MatrixScoring(SubstitutionMatrix.Blosum62, 0, -0.5);
        var r = SequenceAligner.MultipleAlign(["MKV", "MV"], sc);
        // M/M 5 + (K vs gap) -0.5 + V/V 4 = 8.5
        Assert.That(r.AlignedSequences, Is.EqualTo(new[] { "MKV", "M-V" }));
        Assert.That(r.ExactTotalScore, Is.EqualTo(8.5));
        Assert.That(r.TotalScore, Is.EqualTo(9));
        // Cheap gaps beat the W/C mismatch (-2): M/M 5 + two gap columns of -0.5 = 4 (two gap columns), not 5 - 2 = 3.
        var cheap = SequenceAligner.MultipleAlign(["MW", "MC"], new MatrixScoring(SubstitutionMatrix.Blosum62, 0, -0.5));
        Assert.That(cheap.AlignedSequences.Select(r => r.Replace("-", "")), Is.EqualTo(new[] { "MW", "MC" }));
        Assert.That(cheap.AlignedSequences[0].Length, Is.EqualTo(3), "two gap columns (co-optimal arrangements exist)");
        Assert.That(cheap.ExactTotalScore, Is.EqualTo(4.0));
        var dear = SequenceAligner.MultipleAlign(["MW", "MC"], new MatrixScoring(SubstitutionMatrix.Blosum62, 0, -2.5));
        Assert.That(dear.AlignedSequences, Is.EqualTo(new[] { "MW", "MC" }));
        Assert.That(dear.ExactTotalScore, Is.EqualTo(3.0));
    }

    #endregion
}
