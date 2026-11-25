// ANNOT-CODING-001 — Fickett (1982) TESTCODE statistic (F12)
// Evidence: docs/Evidence/ANNOT-CODING-001-Evidence.md (§ Fickett TESTCODE, oracle numbers)
// TestSpec: tests/TestSpecs/ANNOT-CODING-001.md
// Sources: Fickett JW (1982) NAR 10(17):5303-5318, https://doi.org/10.1093/nar/10.17.5303;
//          EMBOSS 6.6.0 tcode.c / tcode.acd / Etcode.dat (kimrutherford/EMBOSS, run as the binary);
//          CPAT 3.0.5 cpmodule/fickett.py `fickett_value` (liguowang/cpat).
// Expected values: Fickett1982 = exact rational TESTCODE (Python fractions) = EMBOSS tcode to its 3 printed
// decimals; Cpat = CPAT fickett_value output (bit-exact repr).

namespace Seqeron.Genomics.Tests.Unit.Annotation;

using System;
using System.Linq;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.GenomeAnnotator;

[TestFixture]
public class GenomeAnnotator_FickettScore_Tests
{
    // Transcript t1 of the CPAT fixture (388 nt, seeded generator, see Evidence).
    private const string T1 =
        "TCATTGGCTATCCTAACCCGACCCTAGGAGCGGTTGGCGTATGGTATGCCGTGAATTTTCTCATTTCCGCACATCGTTCTGCCTATATCTGGACAACATCCCGGCGACTGCGACCCACAGAATCGTCCCTTCTAACGTAGTTCGCATAGTTCCCGTCCGTAGCCGGACTATTCGAACACCCAGTATTCGATTAACTCGGGCTCGTATTAGAGGCGTTAGTGTGTAACAGGTAAGATACGCCAACGGAATTAACCTCTATGGTGACACTCCGCGGAGCCTTCGGACATATAAGTGATCGGGTCTACGTTTGTACTGACGTCTGTGAGTTGGGTCTAATAAATCGCCTACACGTGGAGTCTAACGGGTAGAAGCGTCGAATCCTGATACA";

    // t3 (50 nt): G content 8/50 = 0.16 lies below 0.17, where CPAT's 9-threshold content table
    // (0.19 missing) reads one bin too early.
    private const string T3 = "TTATATTTTCCACCCAACATCGTGCACATTAGGCTCACCGCGCCCTACGG";

    // Human window (CPAT test_files, Human_test_coding_mRNA) whose exact TESTCODE is 0.95.
    private const string Boundary095 =
        "TGCAACTTATTCCCAGATCTTTTCTGAAAAAAAAAGCAGAGAGACTGTTTGCCAGACAGCTGGACCAGAATGAGGTAAAAGAAGATAAGGCCCAACCGAGAAGTCAGAAAAGCTTCACTTTTCATTAGCTCCACATTCCAAAGTAAGAAATCAATGTGCTTAGCTGCAATTCAAAGCAATTAGTCTTTGTTCTCACACAC";

    // Human window whose exact TESTCODE is 0.74.
    private const string Boundary074 =
        "GAGCAATTGTAAATGGTATTGTGCTTTGAATTTTGGTTTTCACATGTTCATTGTTAGTATACAGAAATGTGATTGATTTTTATATGTTGATCTTGTATCCTGCAAACTTGTTGCACTACTTATTAATCCTAGAAGTTTCATTGTGGATTTCTTATGATTTTCTATAATTATCCACCTGTGAATAGGGACAGTTTTATTTC";

    [TestCase(T1, 0.7505, 0.7505)]
    [TestCase(T3, 0.8787, 0.8847)]
    [TestCase("ACGTNNACGTTTGCAANNNNGGCATCGATCGANACGT", 1.3275, 1.2721)]
    [TestCase("gcgcgcgcgcgcatgcatgcatgcaaaaaaaaaaaa", 0.4319, 0.42210000000000003)]
    [TestCase("AC", 0.407, 0.4032)]
    public void CalculateFickettScore_MatchesExactTestcodeAndCpatFickettValue(string sequence, double fickett1982, double cpat)
    {
        Assert.Multiple(() =>
        {
            Assert.That(CalculateFickettScore(sequence), Is.EqualTo(fickett1982),
                "Fickett (1982) tables, A/C/G/T composition, exact 4-decimal value");
            Assert.That(CalculateFickettScore(sequence, FickettTableSet.Fickett1982), Is.EqualTo(fickett1982));
            Assert.That(CalculateFickettScore(sequence, FickettTableSet.Cpat), Is.EqualTo(cpat),
                "CPAT fickett_value, bit-exact");
        });
    }

    [Test]
    public void CalculateFickettScore_ContentBelow017_CpatReadsShiftedBin()
    {
        // T3: G fraction 0.16. Fickett 1982 → last bin (≥ 0.00), G content probability 0.29;
        // CPAT's list lacks 0.19, so "≥ 0" is its 9th bin → 0.33. Difference (0.33 − 0.29) · 0.15 = 0.006.
        double diff = CalculateFickettScore(T3, FickettTableSet.Cpat) - CalculateFickettScore(T3);
        Assert.That(diff, Is.EqualTo(0.006).Within(1e-12));
    }

    [Test]
    public void CalculateFickettScore_CpatShortSequence_ReturnsZero_1982ScoresIt()
    {
        Assert.Multiple(() =>
        {
            Assert.That(CalculateFickettScore("A", FickettTableSet.Cpat), Is.EqualTo(0.0), "CPAT: len < 2 → 0");
            Assert.That(CalculateFickettScore("", FickettTableSet.Cpat), Is.EqualTo(0.0));
            Assert.That(CalculateFickettScore("A"), Is.EqualTo(0.3458), "1982: exact TESTCODE of a single A");
        });
    }

    [Test]
    public void CalculateFickettScore_NoAcgt_1982ReturnsZero_CpatCountsLength()
    {
        Assert.Multiple(() =>
        {
            Assert.That(CalculateFickettScore("NNNN"), Is.EqualTo(0.0), "no A/C/G/T: no composition information");
            Assert.That(CalculateFickettScore("N"), Is.EqualTo(0.0));
            Assert.That(CalculateFickettScore(""), Is.EqualTo(0.0));
            Assert.That(CalculateFickettScore("NNNN", FickettTableSet.Cpat), Is.EqualTo(0.40990000000000004),
                "CPAT fickett_value('NNNN'): fractions 0/4, position 0/(0+1)");
        });
    }

    [Test]
    public void CalculateFickettScore_IsCaseInsensitive()
    {
        Assert.That(CalculateFickettScore(T1.ToLowerInvariant()), Is.EqualTo(CalculateFickettScore(T1)));
        Assert.That(CalculateFickettScore(T1.ToLowerInvariant(), FickettTableSet.Cpat),
            Is.EqualTo(CalculateFickettScore(T1, FickettTableSet.Cpat)));
    }

    [Test]
    public void CalculateFickettScore_InvalidArguments_Throw()
    {
        Assert.Throws<ArgumentNullException>(() => CalculateFickettScore(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => CalculateFickettScore("ACGT", (FickettTableSet)7));
    }

    [TestCase(0.95, FickettCodingCall.Coding)]
    [TestCase(0.9499, FickettCodingCall.NoOpinion)]
    [TestCase(0.74, FickettCodingCall.NoOpinion)]
    [TestCase(0.7399, FickettCodingCall.NonCoding)]
    [TestCase(1.2, FickettCodingCall.Coding)]
    [TestCase(0.0, FickettCodingCall.NonCoding)]
    public void ClassifyFickettScore_AppliesFickettLimits(double score, FickettCodingCall expected)
        => Assert.That(ClassifyFickettScore(score), Is.EqualTo(expected));

    [Test]
    public void CalculateFickettScore_ExactBoundaryWindows_ClassifiedByPublishedRule()
    {
        // EMBOSS tcode sums in float32 and reports these as "No opinion" / "Non-coding"
        // (0.950 / 0.740 printed); the exact statistic equals the limits, so the published
        // rule (≥ 0.95 coding, < 0.74 non-coding) gives Coding / No opinion.
        Assert.Multiple(() =>
        {
            Assert.That(CalculateFickettScore(Boundary095), Is.EqualTo(0.95));
            Assert.That(ClassifyFickettScore(CalculateFickettScore(Boundary095)), Is.EqualTo(FickettCodingCall.Coding));
            Assert.That(CalculateFickettScore(Boundary074), Is.EqualTo(0.74));
            Assert.That(ClassifyFickettScore(CalculateFickettScore(Boundary074)), Is.EqualTo(FickettCodingCall.NoOpinion));
        });
    }

    [Test]
    public void ScanFickettScore_DefaultWindow_MatchesEmbossTcode()
    {
        // EMBOSS: tcode -window 200 -step 3 on T1 (388 nt) → 63 windows, 1-200 … 187-386
        // (loop `pos + window < len`).
        var windows = ScanFickettScore(T1);

        Assert.Multiple(() =>
        {
            Assert.That(windows, Has.Count.EqualTo(63));
            Assert.That(windows[0], Is.EqualTo(new FickettWindow(0, 200, 0.8343, FickettCodingCall.NoOpinion)));
            Assert.That(windows[2], Is.EqualTo(new FickettWindow(6, 206, 0.9794, FickettCodingCall.Coding)));
            Assert.That(windows[60], Is.EqualTo(new FickettWindow(180, 380, 0.7264, FickettCodingCall.NonCoding)));
            Assert.That(windows[62], Is.EqualTo(new FickettWindow(186, 386, 0.7781, FickettCodingCall.NoOpinion)));
            Assert.That(windows.Select(w => w.Start), Is.EqualTo(Enumerable.Range(0, 63).Select(i => 3 * i)));
            Assert.That(windows.All(w => w.Score == CalculateFickettScore(T1.Substring(w.Start, 200))), Is.True,
                "each window is the whole-sequence statistic of its substring");
        });
    }

    [Test]
    public void ScanFickettScore_WindowEqualToLength_YieldsNothing_AsEmboss()
    {
        Assert.That(ScanFickettScore(T1.Substring(0, 200)), Is.Empty);
        Assert.That(ScanFickettScore(T1.Substring(0, 201)), Has.Count.EqualTo(1));
    }

    [Test]
    public void ScanFickettScore_InvalidArguments_Throw()
    {
        Assert.Throws<ArgumentNullException>(() => ScanFickettScore(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => ScanFickettScore(T1, windowSize: 199));
        Assert.Throws<ArgumentOutOfRangeException>(() => ScanFickettScore(T1, stepSize: 0));
    }
}
