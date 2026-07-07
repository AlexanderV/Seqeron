// SEQ-COMPLEX-WINDOW-001 — Windowed Sequence Complexity
// Evidence: docs/Evidence/SEQ-COMPLEX-WINDOW-001-Evidence.md
// TestSpec: tests/TestSpecs/SEQ-COMPLEX-WINDOW-001.md
// Source: Shannon CE (1948) A Mathematical Theory of Communication; Troyanskaya et al. (2002) Bioinformatics 18(5):679-688;
//         Gabrielian & Bolshoy (1999) Comput Chem 23(3-4):263-274; Wikipedia "Linguistic sequence complexity".
//
// Spec: per-window Shannon entropy H = -Σ p·log₂p (bits, 4-base alphabet) and linguistic complexity
// LC = (Σ Vᵢ)/(Σ min(4^i, w-i+1)) with maxWordLength = min(6, w); windows with start i where i+w ≤ L,
// stepped by s. Expected values are derived independently from these formulas, NOT from the implementation.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class SequenceComplexity_CalculateWindowedComplexity_Tests
{
    // ACGTACGT (w=8, maxWord=6): distinct subwords by length 1..6 = 4,4,4,4,4,3 (sum 23);
    // maxima min(4^i,8-i+1) = 4,7,6,5,4,3 (sum 29) ⇒ LC = 23/29 (Wikipedia LC; Evidence dataset).
    private const double UniformWindowLc = 23.0 / 29.0;        // 0.7931034482758621
    // AAAAAAAA (w=8): 1 distinct subword per length (sum 6) over the same maxima (sum 29) ⇒ LC = 6/29.
    private const double HomopolymerWindowLc = 6.0 / 29.0;     // 0.20689655172413793
    private const double MaxDnaEntropy = 2.0;                  // log₂4 — uniform 4-base distribution (Shannon 1948)

    // 24 bp: window0 = ACGTACGT (uniform), window1 = AAAAAAAA (homopolymer), window2 = ACGTACGT (uniform).
    private const string Profiled24 = "ACGTACGTAAAAAAAAACGTACGT";

    #region CalculateWindowedComplexity(DnaSequence, windowSize, stepSize) — enumeration & coordinates

    // M1 — INV-01: L=24, w=8, s=8 ⇒ floor((24-8)/8)+1 = 3 fully-contained windows.
    [Test]
    public void CalculateWindowedComplexity_NonOverlappingWindows_ReturnsExactCount()
    {
        var seq = new DnaSequence(Profiled24);

        var points = SequenceComplexity.CalculateWindowedComplexity(seq, windowSize: 8, stepSize: 8).ToList();

        Assert.That(points.Count, Is.EqualTo(3),
            "L=24, w=8, s=8 ⇒ floor((24-8)/8)+1 = 3 windows fully contained in the sequence.");
    }

    // M2 — INV-02: 0-based start, inclusive end (i+w-1), center = i + w/2.
    [Test]
    public void CalculateWindowedComplexity_NonOverlappingWindows_CoordinatesAreExact()
    {
        var seq = new DnaSequence(Profiled24);

        var points = SequenceComplexity.CalculateWindowedComplexity(seq, windowSize: 8, stepSize: 8).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(points.Select(p => p.WindowStart), Is.EqualTo(new[] { 0, 8, 16 }),
                "Window starts advance by stepSize=8 from 0: 0,8,16.");
            Assert.That(points.Select(p => p.WindowEnd), Is.EqualTo(new[] { 7, 15, 23 }),
                "WindowEnd is inclusive: start + windowSize - 1 = 7,15,23.");
            Assert.That(points.Select(p => p.Position), Is.EqualTo(new[] { 4, 12, 20 }),
                "Position is the window center start + windowSize/2 = 4,12,20.");
        });
    }

    #endregion

    #region Per-window Shannon entropy

    // M3 — uniform window ACGTACGT (A=C=G=T=2) ⇒ Shannon entropy = log₂4 = 2.0 (Shannon 1948, uniform max).
    [Test]
    public void CalculateWindowedComplexity_UniformWindow_ShannonEntropyIsMaximum()
    {
        var seq = new DnaSequence(Profiled24);

        var points = SequenceComplexity.CalculateWindowedComplexity(seq, windowSize: 8, stepSize: 8).ToList();

        Assert.That(points[0].ShannonEntropy, Is.EqualTo(MaxDnaEntropy).Within(1e-10),
            "Window 0 'ACGTACGT' has equal base counts; uniform Shannon entropy = log₂4 = 2.0 bits.");
    }

    // M4 — homopolymer window AAAAAAAA ⇒ deterministic distribution ⇒ Shannon entropy = 0 (Shannon 1948).
    [Test]
    public void CalculateWindowedComplexity_HomopolymerWindow_ShannonEntropyIsZero()
    {
        var seq = new DnaSequence(Profiled24);

        var points = SequenceComplexity.CalculateWindowedComplexity(seq, windowSize: 8, stepSize: 8).ToList();

        Assert.That(points[1].ShannonEntropy, Is.EqualTo(0.0).Within(1e-10),
            "Window 1 'AAAAAAAA' is a single base (p=1); a deterministic distribution has entropy 0.");
    }

    #endregion

    #region Per-window linguistic complexity

    // M5 — uniform window ACGTACGT: LC = (4+4+4+4+4+3)/(4+7+6+5+4+3) = 23/29 (Wikipedia LC, Vmax=min(4^i,8-i+1)).
    [Test]
    public void CalculateWindowedComplexity_UniformWindow_LinguisticComplexityIsExact()
    {
        var seq = new DnaSequence(Profiled24);

        var points = SequenceComplexity.CalculateWindowedComplexity(seq, windowSize: 8, stepSize: 8).ToList();

        Assert.That(points[0].LinguisticComplexity, Is.EqualTo(UniformWindowLc).Within(1e-10),
            "Window 0 'ACGTACGT', maxWord=6: ΣVᵢ=23, ΣVmax=29 ⇒ LC=23/29=0.7931034482758621.");
    }

    // M6 — homopolymer window AAAAAAAA: LC = 6/29 (1 distinct subword per length 1..6 over Σmax=29).
    [Test]
    public void CalculateWindowedComplexity_HomopolymerWindow_LinguisticComplexityIsExact()
    {
        var seq = new DnaSequence(Profiled24);

        var points = SequenceComplexity.CalculateWindowedComplexity(seq, windowSize: 8, stepSize: 8).ToList();

        Assert.That(points[1].LinguisticComplexity, Is.EqualTo(HomopolymerWindowLc).Within(1e-10),
            "Window 1 'AAAAAAAA': one distinct subword per length 1..6 ⇒ ΣVᵢ=6, ΣVmax=29 ⇒ LC=6/29=0.20689655172413793.");
    }

    #endregion

    #region Boundaries and step variants

    // M7 — L < windowSize ⇒ no fully-contained window ⇒ empty profile (corner case: no partial trailing window).
    [Test]
    public void CalculateWindowedComplexity_SequenceShorterThanWindow_ReturnsEmpty()
    {
        var seq = new DnaSequence("ACGTA"); // L=5 < w=8

        var points = SequenceComplexity.CalculateWindowedComplexity(seq, windowSize: 8, stepSize: 8).ToList();

        Assert.That(points, Is.Empty,
            "L=5 < windowSize=8: no window fits entirely in the sequence, so the profile is empty.");
    }

    // S1 — overlapping step s<w: L=24, w=8, s=4 ⇒ floor((24-8)/4)+1 = 5 windows, starts 0,4,8,12,16.
    [Test]
    public void CalculateWindowedComplexity_OverlappingStep_ReturnsExactCountAndStarts()
    {
        var seq = new DnaSequence(Profiled24);

        var points = SequenceComplexity.CalculateWindowedComplexity(seq, windowSize: 8, stepSize: 4).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(points.Count, Is.EqualTo(5),
                "L=24, w=8, s=4 ⇒ floor((24-8)/4)+1 = 5 overlapping windows.");
            Assert.That(points.Select(p => p.WindowStart), Is.EqualTo(new[] { 0, 4, 8, 12, 16 }),
                "Overlapping window starts advance by stepSize=4: 0,4,8,12,16.");
        });
    }

    // S2 — exact-fit single window: L=w=8 ⇒ exactly 1 window spanning the whole sequence.
    [Test]
    public void CalculateWindowedComplexity_SequenceEqualsWindow_ReturnsSingleWindow()
    {
        var seq = new DnaSequence("ACGTACGT"); // L=8=w

        var points = SequenceComplexity.CalculateWindowedComplexity(seq, windowSize: 8, stepSize: 8).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(points.Count, Is.EqualTo(1), "L=windowSize=8 ⇒ exactly one fully-contained window.");
            Assert.That(points[0].WindowStart, Is.EqualTo(0), "The single window starts at index 0.");
            Assert.That(points[0].WindowEnd, Is.EqualTo(7), "The single window ends (inclusive) at index 7.");
        });
    }

    #endregion

    #region Failure modes

    // M8 — null DnaSequence ⇒ ArgumentNullException (contract). Enumeration forced with ToList.
    [Test]
    public void CalculateWindowedComplexity_NullSequence_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => SequenceComplexity.CalculateWindowedComplexity((DnaSequence)null!).ToList(),
            "A null DnaSequence must be rejected with ArgumentNullException.");
    }

    // M9 — windowSize < 1 ⇒ ArgumentOutOfRangeException (contract).
    [Test]
    public void CalculateWindowedComplexity_ZeroWindowSize_Throws()
    {
        var seq = new DnaSequence("ACGTACGT");

        Assert.Throws<ArgumentOutOfRangeException>(
            () => SequenceComplexity.CalculateWindowedComplexity(seq, windowSize: 0).ToList(),
            "windowSize must be ≥ 1; 0 is an invalid window length.");
    }

    // M10 — stepSize < 1 ⇒ ArgumentOutOfRangeException (contract).
    [Test]
    public void CalculateWindowedComplexity_ZeroStepSize_Throws()
    {
        var seq = new DnaSequence("ACGTACGT");

        Assert.Throws<ArgumentOutOfRangeException>(
            () => SequenceComplexity.CalculateWindowedComplexity(seq, windowSize: 4, stepSize: 0).ToList(),
            "stepSize must be ≥ 1; 0 would never advance the window.");
    }

    #endregion

    #region Invariant (C1)

    // C1 — INV-03/INV-04: every window satisfies 0 ≤ H ≤ log₂4 and 0 < LC ≤ 1 (Shannon/LC bounds).
    [TestCase("ACGTACGTGGCCTTAAACGTACGT", 8, 4)]
    [TestCase("AAAAAAAAAAAAAAAA", 8, 8)]
    [TestCase("ACGTACGTACGTACGT", 4, 2)]
    public void CalculateWindowedComplexity_BoundsInvariant_AllWindowsWithinRange(string sequence, int w, int s)
    {
        var seq = new DnaSequence(sequence);

        var points = SequenceComplexity.CalculateWindowedComplexity(seq, w, s).ToList();

        Assert.Multiple(() =>
        {
            foreach (var p in points)
            {
                Assert.That(p.ShannonEntropy, Is.GreaterThanOrEqualTo(0.0).And.LessThanOrEqualTo(MaxDnaEntropy + 1e-10),
                    "Per-window Shannon entropy is bounded by [0, log₂4] for a 4-base alphabet.");
                Assert.That(p.LinguisticComplexity, Is.GreaterThan(0.0).And.LessThanOrEqualTo(1.0 + 1e-10),
                    "Per-window linguistic complexity lies in (0, 1] for DNA windows.");
            }
        });
    }

    #endregion
}
