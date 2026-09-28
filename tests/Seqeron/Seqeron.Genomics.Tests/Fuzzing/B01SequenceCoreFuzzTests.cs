namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Heavy-tier fuzz tests for the contracts introduced by review batch B01
/// (docs/Validation/review-2026-09/B01.md, F2–F13). Random UTF-16 garbage (incl. surrogate halves,
/// '\0', U+017F 'ſ') is fed to the B01 surfaces; every call must either return a theory-correct value or
/// throw the documented ArgumentException family — never hang, never leak a raw runtime exception.
///
/// Test Units: SEQ-VALID-001, SEQ-COMPLEX-001, SEQ-ENTROPY-001, SEQ-GCSKEW-001, SEQ-GC-ANALYSIS-001, B01-SWEEP.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
[Category("Composition")]
public class B01SequenceCoreFuzzTests
{
    private const int Seed = 20260928;

    private static string Garbage(Random rng, int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
        {
            chars[i] = rng.Next(6) switch
            {
                0 => (char)rng.Next(0, 0x10000),          // any UTF-16 unit, incl. surrogate halves
                1 => (char)rng.Next(0, 128),              // any ASCII
                2 => "ſKı\u0000￿"[rng.Next(5)],
                _ => "ACGTUNacgtunRYSWKMBDHV-."[rng.Next(24)],
            };
        }
        return new string(chars);
    }

    private static bool IsAsciiLetterIn(char c, string upperSet) =>
        upperSet.Contains(c) || (c is >= 'a' and <= 'z' && upperSet.Contains((char)(c - 32)));

    /// <summary>Validators never throw on garbage and agree with an independent per-char reference.</summary>
    [Test]
    [CancelAfter(60_000)]
    public void Validators_Garbage_MatchReference()
    {
        var rng = new Random(Seed);
        for (int trial = 0; trial < 3000; trial++)
        {
            string s = Garbage(rng, rng.Next(0, 40));
            var span = s.AsSpan();
            int refDna = Array.FindIndex(s.ToCharArray(), c => !IsAsciiLetterIn(c, "ACGT"));
            int refRna = Array.FindIndex(s.ToCharArray(), c => !IsAsciiLetterIn(c, "ACGU"));
            span.IndexOfInvalidDna().Should().Be(refDna);
            span.IndexOfInvalidRna().Should().Be(refRna);
            span.IsValidDna().Should().Be(refDna < 0);
            span.IsValidRna().Should().Be(refRna < 0);
            span.IsValidIupacDna().Should().Be(s.All(c => IsAsciiLetterIn(c, "ACGTRYSWKMBDHVN")));
            span.IsValidIupacRna().Should().Be(s.All(c => IsAsciiLetterIn(c, "ACGURYSWKMBDHVN")));

            var iupac = new IupacDnaSequence(s);
            iupac.IsValid().Should().Be(s.All(c => IsAsciiLetterIn(c, "ACGTURYSWKMBDHVN") || c is '-' or '.'));
            iupac.Length.Should().Be(s.Length, "ASCII-only upper-casing never changes length");
        }
    }

    /// <summary>
    /// LC string overload on arbitrary alphabets: the suffix-tree path (m &gt; 12) must equal direct distinct-substring
    /// counting over the same upper-cased text (terminator never collides with '\0' or U+FFFF); result finite, ≥ 0.
    /// </summary>
    [Test]
    [CancelAfter(120_000)]
    public void LinguisticComplexity_GarbageAlphabet_SuffixTreeEqualsDirectCount()
    {
        var rng = new Random(Seed + 1);
        for (int trial = 0; trial < 150; trial++)
        {
            string s = Garbage(rng, rng.Next(13, 70));
            int m = rng.Next(13, 90);
            string u = s.ToUpperInvariant();
            int n = u.Length, mm = Math.Min(m, n);
            long obs = 0, pos = 0;
            for (int i = 1; i <= mm; i++)
            {
                var set = new HashSet<string>(StringComparer.Ordinal);
                for (int j = 0; j + i <= n; j++) set.Add(u.Substring(j, i));
                obs += set.Count;
                pos += Math.Min(i < 31 ? 1L << (2 * i) : long.MaxValue, n - i + 1);
            }
            double lc = SequenceComplexity.CalculateLinguisticComplexity(s, m);
            double.IsFinite(lc).Should().BeTrue();
            lc.Should().Be((double)obs / pos, $"trial {trial}, m={m}");
        }
    }

    /// <summary>Shannon entropy on garbage: finite, in [0, 2], equal to the entropy of its A/C/G/T/U projection.</summary>
    [Test]
    public void ShannonEntropy_Garbage_BoundedAndProjectionInvariant()
    {
        var rng = new Random(Seed + 2);
        for (int trial = 0; trial < 2000; trial++)
        {
            string s = Garbage(rng, rng.Next(0, 60));
            double h = SequenceComplexity.CalculateShannonEntropy(s);
            h.Should().BeInRange(0, 2.0 + 1e-12);
            string projected = new(s.Where(c => "ACGTUacgtu".Contains(c)).ToArray());
            h.Should().Be(SequenceComplexity.CalculateShannonEntropy(projected));
        }
    }

    /// <summary>
    /// GC skew / analysis on garbage with random (possibly invalid) window/step: sizes &lt; 1 throw
    /// ArgumentOutOfRangeException eagerly; otherwise results are finite and bounded and enumeration terminates.
    /// </summary>
    [Test]
    [CancelAfter(60_000)]
    public void GcSkewAndAnalysis_GarbageAndRandomSizes_GuardedAndBounded()
    {
        var rng = new Random(Seed + 3);
        for (int trial = 0; trial < 1500; trial++)
        {
            string? s = rng.Next(20) == 0 ? null : Garbage(rng, rng.Next(0, 80));
            int w = rng.Next(-3, 12), step = rng.Next(-3, 12);

            Func<List<GcSkewPoint>> windowed = () => GcSkewCalculator.CalculateWindowedGcSkew(s!, w, step).ToList();
            Func<List<CumulativeGcSkewPoint>> cumulative = () => GcSkewCalculator.CalculateCumulativeGcSkew(s!, w).ToList();
            Func<GcAnalysisResult> analysis = () => GcSkewCalculator.AnalyzeGcContent(s!, w, step);

            if (w < 1 || step < 1)
            {
                windowed.Should().Throw<ArgumentOutOfRangeException>();
                analysis.Should().Throw<ArgumentOutOfRangeException>();
            }
            else
            {
                windowed().Should().OnlyContain(p => p.GcSkew >= -1 && p.GcSkew <= 1);
                var r = analysis();
                double.IsFinite(r.OverallGcContent).Should().BeTrue();
                r.OverallGcContent.Should().BeInRange(0, 100);
                r.WindowedGcContent.Should().OnlyContain(p => p.GcContent >= 0 && p.GcContent <= 100);
                r.GcContentVariance.Should().BeGreaterThanOrEqualTo(0);
            }

            if (w < 1) cumulative.Should().Throw<ArgumentOutOfRangeException>();
            else cumulative().Should().OnlyContain(p => double.IsFinite(p.CumulativeGcSkew));
        }
    }

    /// <summary>
    /// QualitySequence on random quality strings: accepted iff length matches and every char ∈ [offset, 126];
    /// otherwise ArgumentException (never clamping/padding). Accepted strings round-trip; Phred on random p.
    /// </summary>
    [Test]
    public void QualitySequence_RandomQualityStrings_AcceptedIffValid()
    {
        var rng = new Random(Seed + 4);
        for (int trial = 0; trial < 3000; trial++)
        {
            int n = rng.Next(0, 20);
            int offset = rng.Next(3) == 0 ? 64 : 33;
            int qLen = rng.Next(4) == 0 ? rng.Next(0, 22) : n;
            string q = new(Enumerable.Range(0, qLen).Select(_ => rng.Next(3) == 0
                ? (char)rng.Next(0, 300)
                : (char)rng.Next(offset, 127)).ToArray());
            string bases = new(Enumerable.Range(0, n).Select(_ => "ACGTN"[rng.Next(5)]).ToArray());
            bool valid = qLen == n && q.All(c => c >= offset && c <= 126);

            Func<QualitySequence> make = () => new QualitySequence(bases, q, offset);
            if (valid)
            {
                var qs = make();
                qs.GetQualityString(offset).Should().Be(q);
                qs.Qualities.Should().OnlyContain(b => b <= 126 - offset);
            }
            else
            {
                make.Should().Throw<ArgumentException>().Which.Should().NotBeOfType<ArgumentOutOfRangeException>();
            }

            double p = rng.Next(5) switch
            {
                0 => rng.NextDouble() * 4 - 2,
                1 => double.NaN,
                _ => rng.NextDouble(),
            };
            Func<byte> phred = () => QualitySequence.ErrorProbabilityToPhred(p);
            if (double.IsNaN(p) || p < 0 || p > 1) phred.Should().Throw<ArgumentOutOfRangeException>();
            else phred().Should().BeLessThanOrEqualTo(QualitySequence.MaxSangerPhred);
        }
    }

    /// <summary>
    /// IUPAC code algebra on arbitrary chars: CodesMatch never throws and is symmetric; GetIupacCode on random
    /// base multisets (incl. U, lower case and junk) returns one of the 15 codes and is the inverse of ExpandCode
    /// whenever every symbol is a base.
    /// </summary>
    [Test]
    public void IupacCodeAlgebra_RandomInputs_TotalAndConsistent()
    {
        var rng = new Random(Seed + 5);
        for (int trial = 0; trial < 3000; trial++)
        {
            char a = Garbage(rng, 1)[0], b = Garbage(rng, 1)[0];
            IupacDnaSequence.CodesMatch(a, b).Should().Be(IupacDnaSequence.CodesMatch(b, a));

            string bases = new(Enumerable.Range(0, rng.Next(0, 6)).Select(_ => "ACGTUacgtuX"[rng.Next(11)]).ToArray());
            char code = IupacDnaSequence.GetIupacCode(bases);
            "ACGTRYSWKMBDHVN".Should().Contain(code.ToString());
            var set = bases.Select(c => char.ToUpperInvariant(c) == 'U' ? 'T' : char.ToUpperInvariant(c)).Distinct().OrderBy(c => c).ToList();
            if (set.Count > 0 && set.All(c => "ACGT".Contains(c)))
                IupacDnaSequence.ExpandCode(code).OrderBy(c => c).Should().Equal(set);
            else
                code.Should().Be('N');
        }
    }
}
