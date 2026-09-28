// SEQ-REPLICATION-001 — Replication Origin Prediction (cumulative GC-skew minimum)
// Evidence: docs/Evidence/SEQ-REPLICATION-001-Evidence.md
// TestSpec: tests/TestSpecs/SEQ-REPLICATION-001.md
// Source: Grigoriev A (1998) Nucleic Acids Res 26(10):2286-2290; Rosalind BA1F "Minimum Skew Problem".

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GcSkewCalculator_PredictReplicationOrigin_Tests
{
    // Rosalind BA1F sample genome (length 100). Sample output: minimizing positions "53 97".
    private const string Ba1fSample =
        "CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG";

    #region PredictReplicationOrigin(DnaSequence)

    // M1 — Rosalind BA1F sample: origin is the FIRST prefix index minimizing the cumulative skew.
    // Published sample output is "53 97"; min skew value is -4. A windowed/off-by-one model would not give 53/-4.
    [Test]
    public void PredictReplicationOrigin_Ba1fSample_ReturnsFirstMinimizingPositionAndSkew()
    {
        var sequence = new DnaSequence(Ba1fSample);

        var prediction = GcSkewCalculator.PredictReplicationOrigin(sequence);

        Assert.Multiple(() =>
        {
            Assert.That(prediction.PredictedOrigin, Is.EqualTo(53),
                "BA1F sample minimizes skew at prefix indices 53 and 97; first index reported is 53");
            Assert.That(prediction.OriginSkew, Is.EqualTo(-4.0).Within(1e-10),
                "Global minimum of the BA1F cumulative skew diagram is -4");
        });
    }

    // M2 — Per-nucleotide increments G:+1 C:-1 A/T:0 with Skew_0=0.
    // CCGGGG diagram = 0,-1,-2,-1,0,+1,+2 => min -2 @ prefix 2 (origin), max +2 @ prefix 6 (terminus).
    [Test]
    public void PredictReplicationOrigin_SmallGcSequence_MatchesPerNucleotideDiagram()
    {
        var sequence = new DnaSequence("CCGGGG");

        var prediction = GcSkewCalculator.PredictReplicationOrigin(sequence);

        Assert.Multiple(() =>
        {
            Assert.That(prediction.PredictedOrigin, Is.EqualTo(2), "min -2 first reached at prefix index 2");
            Assert.That(prediction.OriginSkew, Is.EqualTo(-2.0).Within(1e-10), "two leading C bases give skew -2");
            Assert.That(prediction.PredictedTerminus, Is.EqualTo(6), "max +2 reached at prefix index 6 (end)");
            Assert.That(prediction.TerminusSkew, Is.EqualTo(2.0).Within(1e-10), "net (G-C) over whole string is +2");
        });
    }

    // M3 — Terminus is the global MAXIMUM of the diagram (Grigoriev 1998; Wikipedia max=terminus).
    // GGGCCC diagram = 0,+1,+2,+3,+2,+1,0 => max +3 @ prefix 3; min 0 @ prefix 0.
    [Test]
    public void PredictReplicationOrigin_GggcccSequence_TerminusAtGlobalMaximum()
    {
        var sequence = new DnaSequence("GGGCCC");

        var prediction = GcSkewCalculator.PredictReplicationOrigin(sequence);

        Assert.Multiple(() =>
        {
            Assert.That(prediction.PredictedTerminus, Is.EqualTo(3), "skew peaks at +3 after the three G bases");
            Assert.That(prediction.TerminusSkew, Is.EqualTo(3.0).Within(1e-10), "global maximum is +3");
            Assert.That(prediction.PredictedOrigin, Is.EqualTo(0), "diagram never drops below Skew_0=0, so origin is prefix 0");
            Assert.That(prediction.OriginSkew, Is.EqualTo(0.0).Within(1e-10), "global minimum is the starting value 0");
        });
    }

    // M4 — Tie-break: when several positions share the extreme value, the FIRST (smallest) index is reported.
    // CCGGCC diagram = 0,-1,-2,-1,0,-1,-2 => min -2 at prefix indices 2 AND 6; expect 2.
    [Test]
    public void PredictReplicationOrigin_RepeatedMinimum_ReturnsFirstIndex()
    {
        var sequence = new DnaSequence("CCGGCC");

        var prediction = GcSkewCalculator.PredictReplicationOrigin(sequence);

        Assert.Multiple(() =>
        {
            Assert.That(prediction.PredictedOrigin, Is.EqualTo(2),
                "min -2 occurs at prefix 2 and prefix 6; the first (2) is reported, not 6");
            Assert.That(prediction.OriginSkew, Is.EqualTo(-2.0).Within(1e-10), "tied minimum value is -2");
        });
    }

    // M5 — INV-6: A and T bases leave the cumulative diagram unchanged (only G/C count).
    // GAATTG diagram = 0,+1,+1,+1,+1,+1,+2 (the A/T runs are flat) => terminus +2 @ 6, origin 0 @ 0.
    [Test]
    public void PredictReplicationOrigin_AtBasesBetweenGc_DoNotChangeDiagram()
    {
        var sequence = new DnaSequence("GAATTG");

        var prediction = GcSkewCalculator.PredictReplicationOrigin(sequence);

        Assert.Multiple(() =>
        {
            Assert.That(prediction.TerminusSkew, Is.EqualTo(2.0).Within(1e-10), "two G, no C => max skew +2 despite A/T");
            Assert.That(prediction.PredictedTerminus, Is.EqualTo(6), "maximum reached only after the final G at prefix 6");
            Assert.That(prediction.OriginSkew, Is.EqualTo(0.0).Within(1e-10), "no C => diagram never goes below 0");
        });
    }

    // M6 — INV-3: Skew_0 = 0 is always in range, so OriginSkew <= 0 <= TerminusSkew.
    [Test]
    public void PredictReplicationOrigin_Ba1fSample_OriginSkewLeZeroLeTerminusSkew()
    {
        var sequence = new DnaSequence(Ba1fSample);

        var prediction = GcSkewCalculator.PredictReplicationOrigin(sequence);

        Assert.Multiple(() =>
        {
            Assert.That(prediction.OriginSkew, Is.LessThanOrEqualTo(0.0), "min includes Skew_0=0 so cannot exceed 0");
            Assert.That(prediction.TerminusSkew, Is.GreaterThanOrEqualTo(0.0), "max includes Skew_0=0 so cannot be below 0");
        });
    }

    // S1 — Flat diagram (no G/C asymmetry): origin = terminus = 0, skews 0, not significant.
    [Test]
    public void PredictReplicationOrigin_NoGcContent_ReturnsZeroPredictionNotSignificant()
    {
        var sequence = new DnaSequence("AAAATTTT");

        var prediction = GcSkewCalculator.PredictReplicationOrigin(sequence);

        Assert.Multiple(() =>
        {
            Assert.That(prediction.PredictedOrigin, Is.EqualTo(0), "flat diagram: min is Skew_0 at prefix 0");
            Assert.That(prediction.PredictedTerminus, Is.EqualTo(0), "flat diagram: max is also Skew_0 at prefix 0");
            Assert.That(prediction.OriginSkew, Is.EqualTo(0.0).Within(1e-10), "no G/C => skew stays 0");
            Assert.That(prediction.TerminusSkew, Is.EqualTo(0.0).Within(1e-10), "no G/C => skew stays 0");
            Assert.That(prediction.IsSignificant, Is.False, "amplitude 0 (max == min) => no detectable asymmetry");
        });
    }

    // S2 — IsSignificant true when the diagram has non-zero amplitude (max > min).
    [Test]
    public void PredictReplicationOrigin_Ba1fSample_IsSignificantTrue()
    {
        var sequence = new DnaSequence(Ba1fSample);

        var prediction = GcSkewCalculator.PredictReplicationOrigin(sequence);

        Assert.That(prediction.IsSignificant, Is.True,
            "BA1F sample has a real skew amplitude (max > min) => detectable strand asymmetry");
    }

    // S4 — INV-4: positions lie within [0, n].
    [Test]
    public void PredictReplicationOrigin_PositionsWithinSequenceBounds()
    {
        var sequence = new DnaSequence("GGGCCC");

        var prediction = GcSkewCalculator.PredictReplicationOrigin(sequence);

        Assert.Multiple(() =>
        {
            Assert.That(prediction.PredictedOrigin, Is.InRange(0, sequence.Length), "origin must be a valid prefix index");
            Assert.That(prediction.PredictedTerminus, Is.InRange(0, sequence.Length), "terminus must be a valid prefix index");
        });
    }

    // C1 — Null DnaSequence throws.
    [Test]
    public void PredictReplicationOrigin_NullDnaSequence_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            GcSkewCalculator.PredictReplicationOrigin((DnaSequence)null!));
    }

    // C3 — Single base boundary: "G" diagram = 0,+1 => origin 0/skew 0, terminus 1/skew +1.
    [Test]
    public void PredictReplicationOrigin_SingleGuanine_OriginZeroTerminusOne()
    {
        var sequence = new DnaSequence("G");

        var prediction = GcSkewCalculator.PredictReplicationOrigin(sequence);

        Assert.Multiple(() =>
        {
            Assert.That(prediction.PredictedOrigin, Is.EqualTo(0), "Skew_0=0 is the minimum");
            Assert.That(prediction.PredictedTerminus, Is.EqualTo(1), "skew rises to +1 after the single G");
            Assert.That(prediction.TerminusSkew, Is.EqualTo(1.0).Within(1e-10), "one G, no C => +1");
        });
    }

    #endregion

    #region PredictReplicationOrigin(string)

    // M1 (delegate) — string overload reproduces the BA1F sample identically to the DnaSequence overload.
    [Test]
    public void PredictReplicationOrigin_StringOverload_Ba1fSample_Matches()
    {
        var prediction = GcSkewCalculator.PredictReplicationOrigin(Ba1fSample);

        Assert.Multiple(() =>
        {
            Assert.That(prediction.PredictedOrigin, Is.EqualTo(53), "string overload delegates to the same core");
            Assert.That(prediction.OriginSkew, Is.EqualTo(-4.0).Within(1e-10), "same minimum skew -4");
        });
    }

    // S3 — Case-insensitive: lowercase input yields the same result as uppercase.
    [Test]
    public void PredictReplicationOrigin_StringOverload_CaseInsensitive()
    {
        var upper = GcSkewCalculator.PredictReplicationOrigin("CCGGGG");
        var lower = GcSkewCalculator.PredictReplicationOrigin("ccgggg");

        Assert.Multiple(() =>
        {
            Assert.That(lower.PredictedOrigin, Is.EqualTo(upper.PredictedOrigin), "lowercase must give same origin");
            Assert.That(lower.PredictedTerminus, Is.EqualTo(upper.PredictedTerminus), "lowercase must give same terminus");
            Assert.That(lower.OriginSkew, Is.EqualTo(upper.OriginSkew).Within(1e-10), "lowercase must give same skew");
        });
    }

    // C2 — Null string returns the zero prediction (not significant).
    [Test]
    public void PredictReplicationOrigin_NullString_ReturnsZeroPrediction()
    {
        var prediction = GcSkewCalculator.PredictReplicationOrigin((string)null!);

        Assert.Multiple(() =>
        {
            Assert.That(prediction.PredictedOrigin, Is.EqualTo(0), "null => zero origin");
            Assert.That(prediction.PredictedTerminus, Is.EqualTo(0), "null => zero terminus");
            Assert.That(prediction.IsSignificant, Is.False, "null => no signal");
        });
    }

    // C2 — Empty string returns the zero prediction (not significant).
    [Test]
    public void PredictReplicationOrigin_EmptyString_ReturnsZeroPrediction()
    {
        var prediction = GcSkewCalculator.PredictReplicationOrigin(string.Empty);

        Assert.Multiple(() =>
        {
            Assert.That(prediction.PredictedOrigin, Is.EqualTo(0), "empty => zero origin");
            Assert.That(prediction.IsSignificant, Is.False, "empty => no signal");
        });
    }

    #endregion

    #region Reference cross-check (Biopython / python BA1F re-implementation)

    // Synthetic circular genome with known ori/ter (100 kb): C-rich lagging-strand blocks
    // [0,30000) and [80000,100000) ("ACGTC": net −1 per unit), G-rich leading strand [30000,80000)
    // ("AGGTC": net +1 per unit). Leading strand G-rich => cumulative minimum at ori, maximum at ter
    // (Grigoriev 1998; Lobry 1996).
    private static string SyntheticGenome() =>
        string.Concat(Enumerable.Repeat("ACGTC", 6000))
        + string.Concat(Enumerable.Repeat("AGGTC", 10000))
        + string.Concat(Enumerable.Repeat("ACGTC", 4000));

    // R1 — Python BA1F re-implementation (per-nucleotide #G−#C, all minimizers enumerated) gives
    // minimizers [29997, 30000, 30001, ...] with value −6000 and maximizers [79998, 79999] with
    // value +4001; the first of each is reported.
    [Test]
    public void PredictReplicationOrigin_SyntheticGenome_MatchesPythonReference()
    {
        var prediction = GcSkewCalculator.PredictReplicationOrigin(new DnaSequence(SyntheticGenome()));

        Assert.Multiple(() =>
        {
            Assert.That(prediction.PredictedOrigin, Is.EqualTo(29997));
            Assert.That(prediction.OriginSkew, Is.EqualTo(-6000.0));
            Assert.That(prediction.PredictedTerminus, Is.EqualTo(79998));
            Assert.That(prediction.TerminusSkew, Is.EqualTo(4001.0));
            Assert.That(prediction.IsSignificant, Is.True);
        });
    }

    // R2 — Circular genome: rotating the start to 50000 (Grigoriev's "arbitrary start") moves the
    // extrema to prefix 79997 (−4000) and 29998 (+6001); mapped back ((p + 50000) mod 100000) they
    // are the same ori 29997 / ter 79998 (python reference). This genome is G/C-balanced (Skew_n = 0),
    // so rotation introduces no wrap-around step.
    [Test]
    public void PredictReplicationOrigin_RotatedCircularGenome_SameLociAfterMappingBack()
    {
        string g = SyntheticGenome();
        string rotated = g[50000..] + g[..50000];

        var prediction = GcSkewCalculator.PredictReplicationOrigin(rotated);

        Assert.Multiple(() =>
        {
            Assert.That(prediction.PredictedOrigin, Is.EqualTo(79997));
            Assert.That(prediction.OriginSkew, Is.EqualTo(-4000.0));
            Assert.That(prediction.PredictedTerminus, Is.EqualTo(29998));
            Assert.That(prediction.TerminusSkew, Is.EqualTo(6001.0));
            Assert.That((prediction.PredictedOrigin + 50000) % g.Length, Is.EqualTo(29997));
            Assert.That((prediction.PredictedTerminus + 50000) % g.Length, Is.EqualTo(79998));
        });
    }

    // R3 — Grigoriev's windowed cumulative skew (numpy.cumsum of Biopython GC_skew, window 1000)
    // on the same genome has its minimum −10 in the window ending at 30000 and its maximum 20/3 in
    // the window ending at 80000: the per-base prediction lies within one window of both.
    [Test]
    public void PredictReplicationOrigin_SyntheticGenome_AgreesWithWindowedGrigorievDiagram()
    {
        string g = SyntheticGenome();
        var cumulative = GcSkewCalculator.CalculateCumulativeGcSkew(g, 1000).ToList();
        var min = cumulative.MinBy(p => p.CumulativeGcSkew);
        var max = cumulative.MaxBy(p => p.CumulativeGcSkew);
        var prediction = GcSkewCalculator.PredictReplicationOrigin(g);

        Assert.Multiple(() =>
        {
            Assert.That(min.CumulativeGcSkew, Is.EqualTo(-10.0).Within(1e-9));
            Assert.That(min.Position, Is.EqualTo(29500), "window [29000,30000) centre");
            Assert.That(max.CumulativeGcSkew, Is.EqualTo(20.0 / 3.0).Within(1e-9));
            Assert.That(max.Position, Is.EqualTo(79500), "window [79000,80000) centre");
            Assert.That(Math.Abs(prediction.PredictedOrigin - 30000), Is.LessThanOrEqualTo(1000));
            Assert.That(Math.Abs(prediction.PredictedTerminus - 80000), Is.LessThanOrEqualTo(1000));
        });
    }

    // R4 — The prediction is the extremum of the canonical cumulative skew at a one-base window
    // (Biopython GC_skew(seq, 1) cumsum on BA1F: min −4 first at prefix 53, max +2 first at 16).
    [Test]
    public void PredictReplicationOrigin_Ba1fSample_EqualsCanonicalCumulativeSkewAtWindowOne()
    {
        var cumulative = GcSkewCalculator.CalculateCumulativeGcSkew(Ba1fSample, 1).ToList();
        var prediction = GcSkewCalculator.PredictReplicationOrigin(Ba1fSample);

        Assert.Multiple(() =>
        {
            Assert.That(cumulative.Min(p => p.CumulativeGcSkew), Is.EqualTo(prediction.OriginSkew));
            Assert.That(cumulative.Max(p => p.CumulativeGcSkew), Is.EqualTo(prediction.TerminusSkew));
            Assert.That(prediction.PredictedTerminus, Is.EqualTo(16));
            Assert.That(prediction.TerminusSkew, Is.EqualTo(2.0));
        });
    }

    #endregion
}
