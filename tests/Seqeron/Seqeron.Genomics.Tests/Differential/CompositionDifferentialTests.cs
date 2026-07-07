// 08_DIFFERENTIAL_TESTING rows 1-7 (Composition). Each test runs the production method against an
// INDEPENDENT oracle written from a different first principle (LINQ counting / a literal IUPAC lookup
// table / a regex / a hand-derived published value), never a copy of the implementation's code path.
// A divergence between the two implementations on identical input fails the test.

using System.Text.RegularExpressions;

namespace Seqeron.Genomics.Tests.Differential;

[TestFixture]
public class CompositionDifferentialTests
{
    private const double Tol = 1e-12;

    // ---- Row 1: SEQ-GC-001 — span switch vs LINQ Count-based oracle (exact GC%) ----

    // Independent oracle: GC% = 100 * |{G,C,g,c}| / |{A,T,G,C,U + lowercase}|, computed with LINQ.
    private static double GcPercentOracle(string seq)
    {
        int gc = seq.Count(c => "GgCc".Contains(c));
        int valid = seq.Count(c => "AaTtGgCcUu".Contains(c));
        return valid == 0 ? 0.0 : (double)gc / valid * 100.0;
    }

    [Test]
    [Category("SEQ-GC-001")]
    [TestCase("GCGC")]
    [TestCase("ATAT")]
    [TestCase("AcGtUu")]            // mixed case incl U
    [TestCase("ACGTNNNN---xyz")]    // ambiguous/other chars excluded from num AND denom
    [TestCase("")]                  // empty -> 0
    [TestCase("NNNN")]              // no valid nucleotides -> 0
    [TestCase("GGGGGGGCCCCCAAATTT")]
    public void GcContent_MatchesLinqCountOracle(string seq)
    {
        Assert.That(seq.AsSpan().CalculateGcContent(), Is.EqualTo(GcPercentOracle(seq)).Within(Tol));
    }

    // ---- Row 2: SEQ-COMP-001 — switch complement vs literal IUPAC lookup-table oracle ----

    // Independent oracle: a dictionary of the IUPAC DNA complement pairs (recognized base -> UPPERCASE
    // complement; unrecognized char passes through verbatim, preserving case).
    private static readonly Dictionary<char, char> DnaComplementTable = new()
    {
        ['A'] = 'T', ['T'] = 'A', ['G'] = 'C', ['C'] = 'G', ['U'] = 'A',
        ['R'] = 'Y', ['Y'] = 'R', ['S'] = 'S', ['W'] = 'W', ['K'] = 'M',
        ['M'] = 'K', ['B'] = 'V', ['D'] = 'H', ['H'] = 'D', ['V'] = 'B', ['N'] = 'N',
    };

    private static char ComplementOracle(char c)
    {
        char up = char.ToUpperInvariant(c);
        return DnaComplementTable.TryGetValue(up, out char comp) ? comp : c;
    }

    private static string ComplementOracle(string seq) => new string(seq.Select(ComplementOracle).ToArray());

    [Test]
    [Category("SEQ-COMP-001")]
    [TestCase("ACGT")]
    [TestCase("acgtu")]                 // case-insensitive, recognized -> uppercase
    [TestCase("RYSWKMBDHVN")]           // full IUPAC ambiguity alphabet
    [TestCase("ACGT-N.x")]              // gaps and non-IUPAC pass through verbatim
    public void Complement_MatchesLookupTableOracle(string seq)
    {
        Span<char> dest = new char[seq.Length];
        Assert.That(seq.AsSpan().TryGetComplement(dest), Is.True);
        Assert.That(new string(dest), Is.EqualTo(ComplementOracle(seq)));
    }

    // ---- Row 3: SEQ-REVCOMP-001 — optimized revcomp vs (complement then reverse) oracle ----

    private static string ReverseComplementOracle(string seq)
    {
        // Two independent steps: complement each base (lookup table), then reverse the whole string.
        var complemented = ComplementOracle(seq).ToCharArray();
        Array.Reverse(complemented);
        return new string(complemented);
    }

    [Test]
    [Category("SEQ-REVCOMP-001")]
    [TestCase("ACGT")]
    [TestCase("AACCGGTT")]
    [TestCase("acgtRYN")]
    [TestCase("ACGT-N.x")]
    [TestCase("A")]
    [TestCase("")]
    public void ReverseComplement_MatchesComplementThenReverseOracle(string seq)
    {
        Span<char> dest = new char[seq.Length];
        Assert.That(seq.AsSpan().TryGetReverseComplement(dest), Is.True);
        Assert.That(new string(dest), Is.EqualTo(ReverseComplementOracle(seq)));
    }

    // ---- Row 4: SEQ-VALID-001 — span-loop validation vs regex oracle (same bool) ----

    private static readonly Regex DnaRegex = new("^[ACGTacgt]*$", RegexOptions.Compiled);

    [Test]
    [Category("SEQ-VALID-001")]
    [TestCase("ACGT")]
    [TestCase("acgt")]
    [TestCase("ACGTN")]     // N is NOT valid DNA here (only ACGT) -> false
    [TestCase("ACGU")]      // U is RNA -> false for DNA
    [TestCase("ACG T")]     // space -> false
    [TestCase("")]          // empty -> true (vacuously)
    [TestCase("12345")]
    public void IsValidDna_MatchesRegexOracle(string seq)
    {
        Assert.That(seq.AsSpan().IsValidDna(), Is.EqualTo(DnaRegex.IsMatch(seq)));
    }

    // ---- Row 5: SEQ-COMPLEX-001 — linguistic complexity vs hand-derived REF values ----

    // The published linguistic-complexity definition (mean over k of |distinct k-mers| / min(4^k, n-k+1))
    // computed BY HAND for tiny sequences, independent of the production code path:
    //   "AAAA"  -> k1 1/4, k2 1/3, k3 1/2, k4 1/1  -> mean 0.520833...
    //   "ACGT"  -> all maximal at every k          -> mean 1.0
    //   "ATATAT"-> k1 2/4, k2 2/5, k3 2/4, k4 2/3, k5 2/2, k6 1/1 -> mean 0.677777...
    [Test]
    [Category("SEQ-COMPLEX-001")]
    [TestCase("AAAA", (0.25 + 1.0 / 3.0 + 0.5 + 1.0) / 4.0)]
    [TestCase("ACGT", 1.0)]
    [TestCase("ATATAT", (0.5 + 0.4 + 0.5 + 2.0 / 3.0 + 1.0 + 1.0) / 6.0)]
    public void LinguisticComplexity_MatchesHandDerivedReference(string seq, double expected)
    {
        Assert.That(SequenceStatistics.CalculateLinguisticComplexity(seq), Is.EqualTo(expected).Within(1e-12));
    }

    // ---- Row 6: SEQ-ENTROPY-001 — optimized Shannon vs naive histogram + closed-form REF ----

    // Independent histogram oracle: group letters (uppercased), -Σ p·log2 p.
    private static double ShannonOracle(string seq)
    {
        var letters = seq.ToUpperInvariant().Where(char.IsLetter).ToList();
        if (letters.Count == 0) return 0.0;
        return -letters.GroupBy(c => c)
                       .Select(g => (double)g.Count() / letters.Count)
                       .Sum(p => p * Math.Log2(p));
    }

    [Test]
    [Category("SEQ-ENTROPY-001")]
    [TestCase("ACGT")]      // 4 equiprobable symbols -> exactly 2.0 bits
    [TestCase("AAAA")]      // 1 symbol -> 0
    [TestCase("AACC")]      // 2 equiprobable -> 1.0 bit
    [TestCase("AA11")]      // digits filtered by IsLetter -> entropy of "AA" = 0
    [TestCase("")]          // empty -> 0
    [TestCase("AAACGTTTGC")]
    public void ShannonEntropy_MatchesHistogramOracle(string seq)
    {
        Assert.That(SequenceStatistics.CalculateShannonEntropy(seq), Is.EqualTo(ShannonOracle(seq)).Within(1e-12));
    }

    [Test]
    [Category("SEQ-ENTROPY-001")]
    public void ShannonEntropy_KnownClosedFormAnchors()
    {
        Assert.That(SequenceStatistics.CalculateShannonEntropy("ACGT"), Is.EqualTo(2.0).Within(Tol));
        Assert.That(SequenceStatistics.CalculateShannonEntropy("AACC"), Is.EqualTo(1.0).Within(Tol));
        Assert.That(SequenceStatistics.CalculateShannonEntropy("AAAA"), Is.EqualTo(0.0).Within(Tol));
    }

    // ---- Row 7: SEQ-GCSKEW-001 — windowed GC skew vs per-window independent count oracle ----

    // Independent scalar oracle: (G - C) / (G + C) counted with LINQ.
    private static double GcSkewOracle(string seq)
    {
        int g = seq.Count(c => c == 'G' || c == 'g');
        int c = seq.Count(ch => ch == 'C' || ch == 'c');
        return g + c == 0 ? 0.0 : (double)(g - c) / (g + c);
    }

    [Test]
    [Category("SEQ-GCSKEW-001")]
    [TestCase("GGGC")]
    [TestCase("GCCC")]
    [TestCase("CCGG")]
    [TestCase("AATT")]   // no G/C -> 0
    [TestCase("")]
    public void GcSkew_Scalar_MatchesCountOracle(string seq)
    {
        Assert.That(GcSkewCalculator.CalculateGcSkew(seq), Is.EqualTo(GcSkewOracle(seq)).Within(Tol));
    }

    [Test]
    [Category("SEQ-GCSKEW-001")]
    public void GcSkew_Windowed_MatchesPerWindowOracle()
    {
        const string seq = "GGGCCCGGGC"; // 10 nt
        const int window = 4, step = 2;
        var points = GcSkewCalculator.CalculateWindowedGcSkew(seq, window, step).ToList();

        // Independent enumeration of the windows and their skew.
        var expected = new List<(int start, int end, int pos, double skew)>();
        for (int i = 0; i + window <= seq.Length; i += step)
        {
            string w = seq.Substring(i, window);
            expected.Add((i, i + window - 1, i + window / 2, GcSkewOracle(w)));
        }

        Assert.That(points.Count, Is.EqualTo(expected.Count));
        for (int k = 0; k < points.Count; k++)
        {
            Assert.That(points[k].WindowStart, Is.EqualTo(expected[k].start), $"start[{k}]");
            Assert.That(points[k].WindowEnd, Is.EqualTo(expected[k].end), $"end[{k}]");
            Assert.That(points[k].Position, Is.EqualTo(expected[k].pos), $"pos[{k}]");
            Assert.That(points[k].GcSkew, Is.EqualTo(expected[k].skew).Within(Tol), $"skew[{k}]");
        }
    }
}
