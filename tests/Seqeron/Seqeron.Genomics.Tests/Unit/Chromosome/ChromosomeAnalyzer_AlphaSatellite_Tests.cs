// CHROM-CENT-001 — Alpha-satellite-specific detection (171 bp monomer + CENP-B box)
// Evidence: docs/Evidence/CHROM-CENT-001-Evidence.md
// TestSpec: tests/TestSpecs/CHROM-CENT-001.md
// Sources:
//   Monomer 171 bp: Willard HF (1985); Waye JS, Willard HF (1987); review Hartley & O'Neill (2019),
//     "Alpha satellite DNA biology" (PMC6121732): "fundamental 171bp monomeric repeat units".
//   CENP-B box 17 bp consensus 5'-YTTCGTTGGAARCGGGA-3': Masumoto H et al. (1989), J Cell Biol 109(4):1963-1973;
//     consensus reported in PMC6121732 and PMC4843215.

using NUnit.Framework;
using Seqeron.Genomics;
using System;
using System.Linq;

namespace Seqeron.Genomics.Tests.Unit.Chromosome;

/// <summary>
/// Tests for the alpha-satellite-SPECIFIC detection added to CHROM-CENT-001:
/// <see cref="ChromosomeAnalyzer.DetectAlphaSatellite"/> and
/// <see cref="ChromosomeAnalyzer.FindCenpBBoxes"/>.
/// </summary>
[TestFixture]
public class ChromosomeAnalyzer_AlphaSatellite_Tests
{
    // A synthetic 171-bp AT-rich monomer. Constructed (not from any reference record) purely as a
    // test fixture: composition is hand-controlled so expected values are derivable. It is NOT used
    // by the implementation (the detector embeds no consensus monomer); it only seeds the test array.
    // Composition: 100 A/T bases + 71 C/G bases over 171 positions  =>  AT content = 100/171.
    private static string BuildAtRichMonomer()
    {
        // 60 'A', 40 'T', 36 'C', 35 'G' = 171 bases; AT = 100, GC = 71.
        var chars = new char[ChromosomeAnalyzer.AlphaSatelliteMonomerLength];
        int idx = 0;
        for (int i = 0; i < 60; i++) chars[idx++] = 'A';
        for (int i = 0; i < 40; i++) chars[idx++] = 'T';
        for (int i = 0; i < 36; i++) chars[idx++] = 'C';
        for (int i = 0; i < 35; i++) chars[idx++] = 'G';
        return new string(chars);
    }

    private static string TandemArray(string monomer, int copies) =>
        string.Concat(Enumerable.Repeat(monomer, copies));

    private static string GenerateRandomSequence(int seed, int length)
    {
        var random = new Random(seed);
        var bases = new[] { 'A', 'C', 'G', 'T' };
        var sequence = new char[length];
        for (int i = 0; i < length; i++)
            sequence[i] = bases[random.Next(4)];
        return new string(sequence);
    }

    #region DetectAlphaSatellite

    [Test]
    public void DetectAlphaSatellite_PerfectTandem171bpAtRichArray_IsDetected()
    {
        // M-ALPHA-1: A clean tandem array of a 171-bp AT-rich monomer is genuine alpha-satellite signal.
        // Expected (hand-derived):
        //   - period 171 perfect tandem => every base equals the base 171 positions back => periodicity = 1.0
        //   - AtContent = 100/171 (exact, since monomer = 100 AT + 71 GC and the array is exact copies)
        //   - IsAlphaSatellite = true (periodicity 1.0 >= 0.50 AND AT 100/171 ≈ 0.585 > 0.50)
        string monomer = BuildAtRichMonomer();
        string array = TandemArray(monomer, 20); // 20 monomers = 3420 bp

        var result = ChromosomeAnalyzer.DetectAlphaSatellite(array);

        Assert.Multiple(() =>
        {
            Assert.That(result.IsAlphaSatellite, Is.True,
                "Tandem 171-bp AT-rich array must be called alpha-satellite");
            Assert.That(result.PeriodicityScore, Is.EqualTo(1.0).Within(1e-10),
                "Perfect tandem array is identical to itself shifted by the monomer period");
            Assert.That(result.BestPeriod, Is.EqualTo(171),
                "Best period must be the 171-bp monomer length");
            Assert.That(result.AtContent, Is.EqualTo(100.0 / 171.0).Within(1e-10),
                "Monomer is 100 AT bases out of 171; the array preserves the ratio exactly");
        });
    }

    [Test]
    public void DetectAlphaSatellite_RandomSequence_IsNotDetected()
    {
        // M-ALPHA-2: Random DNA has no 171-bp tandem periodicity and ~balanced composition,
        // so it must NOT be called alpha-satellite. Periodicity for random 4-letter DNA ≈ 0.25.
        string random = GenerateRandomSequence(seed: 7, length: 3420);

        var result = ChromosomeAnalyzer.DetectAlphaSatellite(random);

        Assert.Multiple(() =>
        {
            Assert.That(result.IsAlphaSatellite, Is.False,
                "Random sequence is not a 171-bp tandem array and must not be alpha-satellite");
            Assert.That(result.PeriodicityScore, Is.LessThan(0.50),
                "Random DNA self-similarity at any period is far below the 0.50 tandem threshold");
        });
    }

    [Test]
    public void DetectAlphaSatellite_AtRichButNonRepetitive_IsNotDetected()
    {
        // M-ALPHA-3: AT-richness ALONE is not sufficient — without 171-bp periodicity it must fail.
        // A non-repetitive but AT-biased sequence over the full 4-letter alphabet (70% A/T, 30% G/C,
        // randomly ordered). Expected self-similarity at any period ≈ Σ p_i^2 = 0.35^2*2 + 0.15^2*2 ≈ 0.29,
        // which is below the 0.50 tandem threshold even though AT content (~0.70) exceeds 0.50.
        var random = new Random(99);
        var chars = new char[3420];
        for (int i = 0; i < chars.Length; i++)
        {
            int r = random.Next(100);
            chars[i] = r < 35 ? 'A' : r < 70 ? 'T' : r < 85 ? 'G' : 'C';
        }
        string atRich = new string(chars);

        var result = ChromosomeAnalyzer.DetectAlphaSatellite(atRich);

        Assert.Multiple(() =>
        {
            Assert.That(result.AtContent, Is.GreaterThan(0.50),
                "Sequence is AT-biased so AT content exceeds the 0.50 alpha-satellite threshold");
            Assert.That(result.PeriodicityScore, Is.LessThan(0.50),
                "A non-repetitive AT-biased sequence has self-similarity ≈ 0.29, below the tandem threshold");
            Assert.That(result.IsAlphaSatellite, Is.False,
                "AT-richness without 171-bp tandem periodicity must not be called alpha-satellite");
        });
    }

    [Test]
    public void DetectAlphaSatellite_GcRichTandem16bpRepeat_IsNotDetected()
    {
        // M-ALPHA-4: A highly repetitive but NON-alphoid array (GC-rich, 16-bp period, not 171 bp)
        // must NOT be flagged. This is the false positive the generic heuristic suffers from.
        string repeat16 = string.Concat(Enumerable.Repeat("GCGCGCGCGCGCGCGC", 200)); // 3200 bp, period 16, GC-rich

        var result = ChromosomeAnalyzer.DetectAlphaSatellite(repeat16);

        Assert.Multiple(() =>
        {
            Assert.That(result.AtContent, Is.EqualTo(0.0).Within(1e-10),
                "Pure GC sequence has AT content 0");
            Assert.That(result.IsAlphaSatellite, Is.False,
                "A GC-rich 16-bp tandem repeat is not AT-rich, so it is not alpha-satellite");
        });
    }

    [Test]
    public void DetectAlphaSatellite_ArrayCarryingCenpBBoxes_CountsThem()
    {
        // M-ALPHA-5: A tandem array whose monomer carries one CENP-B box must report one box per monomer.
        // Build a 171-bp AT-rich monomer with a single CENP-B box instance embedded at a fixed offset.
        string box = "TTTCGTTGGAAGCGGGA"; // Y->T, R->G instance of YTTCGTTGGAARCGGGA (17 bp)
        Assert.That(box.Length, Is.EqualTo(17), "fixture sanity: box is 17 bp");

        // Monomer = 77 'A' + box(17) + 77 'T' = 171 bp; AT content high, one box per monomer.
        string monomer = new string('A', 77) + box + new string('T', 77);
        Assert.That(monomer.Length, Is.EqualTo(171), "fixture sanity: monomer is 171 bp");

        const int copies = 10;
        string array = TandemArray(monomer, copies);

        var result = ChromosomeAnalyzer.DetectAlphaSatellite(array);

        Assert.Multiple(() =>
        {
            Assert.That(result.IsAlphaSatellite, Is.True,
                "171-bp tandem AT-rich array with CENP-B boxes is alpha-satellite");
            Assert.That(result.BestPeriod, Is.EqualTo(171));
            Assert.That(result.CenpBBoxCount, Is.EqualTo(copies),
                "Exactly one CENP-B box per monomer => 10 boxes in 10 tandem copies");
        });
    }

    [Test]
    public void DetectAlphaSatellite_EmptySequence_ReturnsNoSignal()
    {
        var result = ChromosomeAnalyzer.DetectAlphaSatellite("");

        Assert.Multiple(() =>
        {
            Assert.That(result.IsAlphaSatellite, Is.False);
            Assert.That(result.PeriodicityScore, Is.EqualTo(0));
            Assert.That(result.BestPeriod, Is.EqualTo(0));
            Assert.That(result.AtContent, Is.EqualTo(0));
            Assert.That(result.CenpBBoxCount, Is.EqualTo(0));
        });
    }

    [Test]
    public void DetectAlphaSatellite_NullSequence_ReturnsNoSignal()
    {
        var result = ChromosomeAnalyzer.DetectAlphaSatellite(null!);

        Assert.That(result.IsAlphaSatellite, Is.False);
    }

    [Test]
    public void DetectAlphaSatellite_TooShortToMeasurePeriod_ReturnsNoSignal()
    {
        // Below 171 + tolerance + 1 bp the monomer period cannot be measured => no-signal result.
        string shortSeq = new string('A', 100);

        var result = ChromosomeAnalyzer.DetectAlphaSatellite(shortSeq);

        Assert.Multiple(() =>
        {
            Assert.That(result.IsAlphaSatellite, Is.False);
            Assert.That(result.BestPeriod, Is.EqualTo(0));
            Assert.That(result.PeriodicityScore, Is.EqualTo(0));
        });
    }

    [Test]
    public void DetectAlphaSatellite_NonAcgtBases_ExcludedFromAtContentDenominator()
    {
        // M-ALPHA-6: Non-ACGT symbols (e.g. 'N') must not pollute the AT fraction. AT content is
        // defined over ACGT bases only, so N is excluded from BOTH numerator and denominator.
        // Hand-derived: monomer = 60 A + 40 T + 36 C + 30 G + 5 N = 171 bp.
        //   ACGT bases = 166 per monomer; AT = 100 per monomer.
        //   AtContent = 100/166 ≈ 0.6024 (NOT 100/171).
        //   Periodicity stays 1.0 because position i and i-171 are the identical monomer base
        //   (including the N positions, which compare N == N).
        string monomer = new string('A', 60) + new string('T', 40)
                       + new string('C', 36) + new string('G', 30) + new string('N', 5);
        Assert.That(monomer.Length, Is.EqualTo(171), "fixture sanity: monomer is 171 bp");
        string array = TandemArray(monomer, 20);

        var result = ChromosomeAnalyzer.DetectAlphaSatellite(array);

        Assert.Multiple(() =>
        {
            Assert.That(result.AtContent, Is.EqualTo(100.0 / 166.0).Within(1e-10),
                "N bases are excluded from the ACGT denominator: 100 AT / 166 ACGT");
            Assert.That(result.PeriodicityScore, Is.EqualTo(1.0).Within(1e-10),
                "Perfect tandem: every base (including N) equals the base 171 positions upstream");
            Assert.That(result.BestPeriod, Is.EqualTo(171));
            Assert.That(result.IsAlphaSatellite, Is.True,
                "periodicity 1.0 >= 0.50 AND AT 100/166 > 0.50");
        });
    }

    [Test]
    public void DetectAlphaSatellite_MixedCaseInput_MatchesUppercase()
    {
        string monomer = BuildAtRichMonomer();
        string array = TandemArray(monomer, 15);

        var upper = ChromosomeAnalyzer.DetectAlphaSatellite(array);
        var lower = ChromosomeAnalyzer.DetectAlphaSatellite(array.ToLowerInvariant());

        Assert.Multiple(() =>
        {
            Assert.That(lower.IsAlphaSatellite, Is.EqualTo(upper.IsAlphaSatellite));
            Assert.That(lower.PeriodicityScore, Is.EqualTo(upper.PeriodicityScore).Within(1e-10));
            Assert.That(lower.BestPeriod, Is.EqualTo(upper.BestPeriod));
            Assert.That(lower.AtContent, Is.EqualTo(upper.AtContent).Within(1e-10));
            Assert.That(lower.CenpBBoxCount, Is.EqualTo(upper.CenpBBoxCount));
        });
    }

    #endregion

    #region FindCenpBBoxes

    [Test]
    public void FindCenpBBoxes_CanonicalBox_FoundAtPositionZero()
    {
        // C-ALPHA-1: The exact consensus instance with Y->C, R->A must match at index 0.
        // YTTCGTTGGAARCGGGA  with Y=C, R=A  => CTTCGTTGGAAACGGGA (17 bp).
        string box = "CTTCGTTGGAAACGGGA";
        Assert.That(box.Length, Is.EqualTo(17));

        var hits = ChromosomeAnalyzer.FindCenpBBoxes(box);

        Assert.Multiple(() =>
        {
            Assert.That(hits, Has.Count.EqualTo(1), "single 17-bp box yields exactly one hit");
            Assert.That(hits[0], Is.EqualTo(0), "match starts at index 0");
        });
    }

    [Test]
    public void FindCenpBBoxes_BothAmbiguityResolutions_Match()
    {
        // C-ALPHA-2: Y matches both C and T; R matches both A and G. All four corner instances match.
        Assert.Multiple(() =>
        {
            Assert.That(ChromosomeAnalyzer.FindCenpBBoxes("CTTCGTTGGAAACGGGA"), Has.Count.EqualTo(1), "Y=C,R=A");
            Assert.That(ChromosomeAnalyzer.FindCenpBBoxes("TTTCGTTGGAAACGGGA"), Has.Count.EqualTo(1), "Y=T,R=A");
            Assert.That(ChromosomeAnalyzer.FindCenpBBoxes("CTTCGTTGGAAGCGGGA"), Has.Count.EqualTo(1), "Y=C,R=G");
            Assert.That(ChromosomeAnalyzer.FindCenpBBoxes("TTTCGTTGGAAGCGGGA"), Has.Count.EqualTo(1), "Y=T,R=G");
        });
    }

    [Test]
    public void FindCenpBBoxes_NonAmbiguousPositionViolated_NoMatch()
    {
        // C-ALPHA-3: A fixed (non-ambiguous) consensus position must match exactly.
        // Mutate position 2 (consensus 'T') to 'A'  =>  no match.
        string broken = "CTACGTTGGAAACGGGA"; // pos2 T->A
        Assert.That(broken.Length, Is.EqualTo(17));

        var hits = ChromosomeAnalyzer.FindCenpBBoxes(broken);

        Assert.That(hits, Is.Empty, "violating a fixed consensus base must prevent a match");
    }

    [Test]
    public void FindCenpBBoxes_WrongAmbiguityBase_NoMatch()
    {
        // C-ALPHA-4: Y must NOT match A or G. Putting 'A' in the leading Y position => no match.
        string box = "ATTCGTTGGAAACGGGA"; // first base A, but consensus is Y = C/T
        var hits = ChromosomeAnalyzer.FindCenpBBoxes(box);
        Assert.That(hits, Is.Empty, "Y position must reject A");
    }

    [Test]
    public void FindCenpBBoxes_BoxInsideLongerSequence_ReportsCorrectOffset()
    {
        // C-ALPHA-5: Position reporting is 0-based offset into the sequence.
        string flankLeft = new string('A', 50);
        string box = "TTTCGTTGGAAGCGGGA";
        string flankRight = new string('A', 30);
        string seq = flankLeft + box + flankRight;

        var hits = ChromosomeAnalyzer.FindCenpBBoxes(seq);

        Assert.Multiple(() =>
        {
            Assert.That(hits, Has.Count.EqualTo(1));
            Assert.That(hits[0], Is.EqualTo(50), "box begins right after the 50-bp left flank");
        });
    }

    [Test]
    public void FindCenpBBoxes_NoBox_ReturnsEmpty()
    {
        var hits = ChromosomeAnalyzer.FindCenpBBoxes(new string('A', 200));
        Assert.That(hits, Is.Empty);
    }

    [Test]
    public void FindCenpBBoxes_EmptyOrShort_ReturnsEmpty()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ChromosomeAnalyzer.FindCenpBBoxes(""), Is.Empty);
            Assert.That(ChromosomeAnalyzer.FindCenpBBoxes(null!), Is.Empty);
            Assert.That(ChromosomeAnalyzer.FindCenpBBoxes("CTTCGTTGGAAACGGG"), Is.Empty,
                "16-bp input is shorter than the 17-bp box");
        });
    }

    #endregion

    #region Sourced Constants

    [Test]
    public void AlphaSatelliteMonomerLength_Is171()
    {
        // Willard 1985 / Waye & Willard 1987 / PMC6121732: fundamental 171-bp monomer.
        Assert.That(ChromosomeAnalyzer.AlphaSatelliteMonomerLength, Is.EqualTo(171));
    }

    [Test]
    public void CenpBBoxConsensus_Is17bpCanonicalMotif()
    {
        // Masumoto et al. 1989: 17-bp consensus 5'-YTTCGTTGGAARCGGGA-3'.
        Assert.Multiple(() =>
        {
            Assert.That(ChromosomeAnalyzer.CenpBBoxConsensus, Is.EqualTo("YTTCGTTGGAARCGGGA"));
            Assert.That(ChromosomeAnalyzer.CenpBBoxConsensus.Length, Is.EqualTo(17));
        });
    }

    #endregion
}
