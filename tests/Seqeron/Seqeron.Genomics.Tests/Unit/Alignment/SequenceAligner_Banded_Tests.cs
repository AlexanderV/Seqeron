using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-GLOBAL-001 / ALIGN-LOCAL-001 / ALIGN-SEMI-001 (audit R14-12, B09 F305): banded alignment (Ukkonen 1985; Chao, Pearson &amp; Miller 1992) and X-drop on the core aligners, locked against the full dynamic programs on random pairs:
/// a band that covers the matrix gives the full result field for field (global, local and semi-global; linear and affine; scalar and substitution-matrix scoring); a narrow band never beats the full optimum and gives the full rows when
/// the full alignment lies inside it; <see cref="SequenceAligner.GlobalAlignBandedExact(string, string, ScoringMatrix?, bool, int)"/> is exact with a proof. Fast-tier property tests (seeded random pairs).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-GLOBAL-001")]
public class SequenceAligner_Banded_Tests
{
    private static string Dna(Random r, int n) => new([.. Enumerable.Range(0, n).Select(_ => "ACGT"[r.Next(4)])]);

    private static string Protein(Random r, int n) => new([.. Enumerable.Range(0, n).Select(_ => "ACDEFGHIKLMNPQRSTVWY"[r.Next(20)])]);

    private static string Mutate(Random r, string s, string alphabet, double rate)
    {
        var sb = new System.Text.StringBuilder();
        foreach (char c in s)
        {
            double x = r.NextDouble();
            if (x < rate / 3)
                sb.Append(alphabet[r.Next(alphabet.Length)]);
            else if (x >= rate / 3 && x < rate * 2 / 3)
            {
                // deletion: the letter is dropped
            }
            else if (x < rate)
                sb.Append(c).Append(alphabet[r.Next(alphabet.Length)]);
            else
                sb.Append(c);
        }

        return sb.Length == 0 ? "A" : sb.ToString();
    }

    private static IEnumerable<(string A, string B, ScoringMatrix S)> Pairs(int seed, int count, bool protein = false)
    {
        var r = new Random(seed);
        for (int i = 0; i < count; i++)
        {
            string a = protein ? Protein(r, r.Next(3, 40)) : Dna(r, r.Next(3, 45));
            string b = r.Next(4) == 0 ? (protein ? Protein(r, r.Next(3, 40)) : Dna(r, r.Next(3, 45))) : Mutate(r, a, protein ? "ACDEFGHIKLMNPQRSTVWY" : "ACGT", 0.2);
            var s = new ScoringMatrix(r.Next(1, 5), -r.Next(1, 4), -r.Next(0, 6), -r.Next(1, 4));
            yield return (a, b, s);
        }
    }

    private static BandedAlignOptions Whole(string a, string b, AlignmentType t, bool affine, ScoringMatrix? s, EndGaps ends = EndGaps.None) =>
        new() { Type = t, Affine = affine, Scoring = s, EndGaps = ends, DiagonalLow = -a.Length, DiagonalHigh = b.Length };

    [Test]
    public void BandCoveringTheMatrix_EqualsTheFullPrograms_GlobalLocalSemiGlobal_LinearAndAffine()
    {
        int checkedCases = 0;
        foreach (var (a, b, s) in Pairs(14012, 120))
        {
            Assert.That(SequenceAligner.BandedAlign(a, b, Whole(a, b, AlignmentType.Global, false, s)), Is.EqualTo(SequenceAligner.GlobalAlign(a, b, s)), $"global linear {a} / {b}");
            Assert.That(SequenceAligner.BandedAlign(a, b, Whole(a, b, AlignmentType.Global, true, s)), Is.EqualTo(SequenceAligner.GlobalAlignAffine(a, b, s)), $"global affine {a} / {b}");
            Assert.That(SequenceAligner.BandedAlign(a, b, Whole(a, b, AlignmentType.Local, false, s)), Is.EqualTo(SequenceAligner.LocalAlign(a, b, s)), $"local linear {a} / {b}");
            Assert.That(SequenceAligner.BandedAlign(a, b, Whole(a, b, AlignmentType.Local, true, s)), Is.EqualTo(SequenceAligner.LocalAlignAffine(a, b, s)), $"local affine {a} / {b}");
            foreach (var ends in new[] { EndGaps.Fitting, EndGaps.FittingSequence2, EndGaps.All, EndGaps.OverlapSuffix1Prefix2 })
            {
                Assert.That(SequenceAligner.BandedAlign(a, b, Whole(a, b, AlignmentType.SemiGlobal, false, s, ends)), Is.EqualTo(SequenceAligner.SemiGlobalAlign(a, b, s, ends)), $"semi linear {ends} {a} / {b}");
                Assert.That(SequenceAligner.BandedAlign(a, b, Whole(a, b, AlignmentType.SemiGlobal, true, s, ends)), Is.EqualTo(SequenceAligner.SemiGlobalAlignAffine(a, b, s, ends)), $"semi affine {ends} {a} / {b}");
            }

            checkedCases++;
        }

        Assert.That(checkedCases, Is.EqualTo(120));
    }

    [Test]
    public void BandCoveringTheMatrix_EqualsTheFullPrograms_WithASubstitutionMatrix()
    {
        var matrix = SubstitutionMatrix.Get("BLOSUM62");
        foreach (var (a, b, _) in Pairs(14013, 80, protein: true))
        {
            var ms = new MatrixScoring(matrix, -10, -0.5);
            Assert.That(SequenceAligner.BandedAlign(a, b, new BandedAlignOptions { Type = AlignmentType.Global, Affine = true, MatrixScoring = ms, DiagonalLow = -a.Length, DiagonalHigh = b.Length }), Is.EqualTo(SequenceAligner.GlobalAlignAffine(a, b, ms)), $"{a} / {b}");
            Assert.That(SequenceAligner.BandedAlign(a, b, new BandedAlignOptions { Type = AlignmentType.Local, Affine = true, MatrixScoring = ms, DiagonalLow = -a.Length, DiagonalHigh = b.Length }), Is.EqualTo(SequenceAligner.LocalAlignAffine(a, b, ms)), $"{a} / {b}");
            var linear = new MatrixScoring(matrix, 0, -1.5);
            Assert.That(SequenceAligner.BandedAlign(a, b, new BandedAlignOptions { Type = AlignmentType.Global, MatrixScoring = linear, DiagonalLow = -a.Length, DiagonalHigh = b.Length }), Is.EqualTo(SequenceAligner.GlobalAlign(a, b, linear)), $"{a} / {b}");
        }
    }

    // the diagonals j - i the path of an alignment visits
    private static (int Low, int High) PathBand(AlignmentResult r)
    {
        int i = r.StartPosition1, j = r.StartPosition2, low = j - i, high = j - i;
        for (int k = 0; k < r.AlignedSequence1.Length; k++)
        {
            if (r.AlignedSequence1[k] != '-')
                i++;
            if (r.AlignedSequence2[k] != '-')
                j++;
            low = Math.Min(low, j - i);
            high = Math.Max(high, j - i);
        }

        return (low, high);
    }

    [Test]
    public void NarrowBand_NeverBeatsTheFullOptimum_AndEqualsItWhenTheFullPathIsInside()
    {
        int inside = 0, outside = 0;
        foreach (var (a, b, s) in Pairs(14014, 250))
        {
            foreach (bool affine in new[] { false, true })
            {
                var full = affine ? SequenceAligner.GlobalAlignAffine(a, b, s) : SequenceAligner.GlobalAlign(a, b, s);
                int d = b.Length - a.Length;
                for (int w = 0; w <= 6; w += 3)
                {
                    var band = SequenceAligner.BandedAlign(a, b, GlobalBandOf(a, b, s, affine, w));
                    Assert.That(band.Score, Is.LessThanOrEqualTo(full.Score), $"{a} / {b} w={w}");
                    var (lo, hi) = PathBand(full);
                    if (lo >= Math.Min(0, d) - w && hi <= Math.Max(0, d) + w)
                    {
                        inside++;
                        Assert.That(band, Is.EqualTo(full), $"{a} / {b} w={w} affine={affine}");
                    }
                    else
                    {
                        outside++;
                    }
                }
            }
        }

        Assert.That(inside, Is.GreaterThan(300));
        Assert.That(outside, Is.GreaterThan(30));
    }

    private static BandedAlignOptions GlobalBandOf(string a, string b, ScoringMatrix s, bool affine, int w) =>
        new() { Type = AlignmentType.Global, Affine = affine, Scoring = s, DiagonalLow = Math.Min(0, b.Length - a.Length) - w, DiagonalHigh = Math.Max(0, b.Length - a.Length) + w };

    [Test]
    public void LocalBand_AroundTheSeedDiagonal_EqualsTheFullLocalAlignment_WhenInside()
    {
        var r = new Random(14015);
        int equal = 0;
        for (int t = 0; t < 120; t++)
        {
            string core = Dna(r, 30);
            string a = Dna(r, r.Next(5, 20)) + core + Dna(r, r.Next(5, 20));
            string b = Dna(r, r.Next(5, 20)) + Mutate(r, core, "ACGT", 0.1) + Dna(r, r.Next(5, 20));
            var s = new ScoringMatrix(2, -3, -4, -2);
            foreach (bool affine in new[] { false, true })
            {
                var full = affine ? SequenceAligner.LocalAlignAffine(a, b, s) : SequenceAligner.LocalAlign(a, b, s);
                if (full.AlignedSequence1.Length == 0)
                    continue;
                var (lo, hi) = PathBand(full);
                var band = affine ? SequenceAligner.LocalAlignAffineBanded(a, b, s, lo - 2, hi + 2) : SequenceAligner.LocalAlignBanded(a, b, s, lo - 2, hi + 2);
                Assert.That(band, Is.EqualTo(full), $"{a} / {b}");
                equal++;
            }
        }

        Assert.That(equal, Is.GreaterThan(150));
    }

    [Test]
    public void GlobalAlignBandedExact_EqualsTheFullProgram_WithAProof_AndFewCells()
    {
        long cells = 0, full = 0;
        foreach (var (a, b, s) in Pairs(14016, 200))
        {
            foreach (bool affine in new[] { false, true })
            {
                var exact = SequenceAligner.GlobalAlignBandedExact(a, b, s, affine, initialBandWidth: 2);
                var reference = affine ? SequenceAligner.GlobalAlignAffine(a, b, s) : SequenceAligner.GlobalAlign(a, b, s);
                Assert.That(exact.Alignment, Is.EqualTo(reference), $"{a} / {b} affine={affine}");
                cells += exact.Cells;
                full += (long)(a.Length + 1) * (b.Length + 1) * exact.Rounds;
            }
        }

        // similar sequences are proven in a narrow band: far fewer cells than rounds of the whole matrix
        Assert.That(cells, Is.LessThan(full));
        var longA = Dna(new Random(1), 4000);
        var longB = Mutate(new Random(2), longA, "ACGT", 0.02);
        var s2 = new ScoringMatrix(2, -3, -5, -2);
        var big = SequenceAligner.GlobalAlignBandedExact(longA, longB, s2, affine: true, initialBandWidth: 8);
        Assert.That(big.Proven, Is.True);
        Assert.That(big.Cells, Is.LessThan(4000L * 4000 / 4));
        Assert.That(big.Alignment.Score, Is.EqualTo(SequenceAligner.GlobalAlignAffine(longA, longB, s2).Score));
        var matrix = new MatrixScoring(SubstitutionMatrix.Get("BLOSUM62"), -10, -0.5);
        string pa = Protein(new Random(3), 300), pb = Mutate(new Random(4), pa, "ACDEFGHIKLMNPQRSTVWY", 0.1);
        Assert.That(SequenceAligner.GlobalAlignBandedExact(pa, pb, matrix, affine: true, initialBandWidth: 4).Alignment, Is.EqualTo(SequenceAligner.GlobalAlignAffine(pa, pb, matrix)));
    }

    [Test]
    public void XDrop_WithALargeX_EqualsTheBand_AndASmallXNeverBeatsIt()
    {
        foreach (var (a, b, s) in Pairs(14017, 150))
        {
            var options = Whole(a, b, AlignmentType.Local, true, s);
            var full = SequenceAligner.BandedAlign(a, b, options);
            Assert.That(SequenceAligner.BandedAlign(a, b, options with { XDrop = 100_000 }), Is.EqualTo(full));
            Assert.That(SequenceAligner.BandedAlign(a, b, options with { XDrop = 3 }).Score, Is.LessThanOrEqualTo(full.Score));
        }
    }

    [Test]
    public void Refusals()
    {
        Assert.Throws<ArgumentException>(() => SequenceAligner.BandedAlign("ACGT", "ACGTACGT", new BandedAlignOptions { DiagonalLow = -1, DiagonalHigh = 1 }));
        Assert.Throws<ArgumentException>(() => SequenceAligner.BandedAlign("ACGT", "ACGT", new BandedAlignOptions { DiagonalLow = 1, DiagonalHigh = 3 }));
        Assert.Throws<ArgumentException>(() => SequenceAligner.BandedAlign("ACGT", "ACGT", new BandedAlignOptions { DiagonalLow = 2, DiagonalHigh = 1 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.GlobalAlignBanded("ACGT", "ACGT", null, -1));
        Assert.Throws<ArgumentException>(() => SequenceAligner.BandedAlign(new string('A', 20000), new string('A', 20000), new BandedAlignOptions { DiagonalLow = -20000, DiagonalHigh = 20000 }));
        var empty = SequenceAligner.LocalAlignBanded("AAAA", "CCCC", null, -3, 3);
        Assert.That(empty, Is.EqualTo(SequenceAligner.LocalAlign("AAAA", "CCCC")));
    }
}
