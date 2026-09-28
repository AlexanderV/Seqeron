// SEQ-PI-001 — Isoelectric Point (pI) Calculation
// Evidence: docs/Evidence/SEQ-PI-001-Evidence.md
// TestSpec: tests/TestSpecs/SEQ-PI-001.md
// Sources (reference implementations run in the 2026-09 review):
//   • EMBOSS 6.6.0 `iep` binary (nucleus/embiep.c + data/Epk.dat, Amino 7.5 / Carboxyl 3.6),
//     https://raw.githubusercontent.com/kimrutherford/EMBOSS/master/emboss/data/Epk.dat — every
//     EMBOSS expected value below is the `iep -auto` "Isoelectric Point" output (4 dp) rounded to 2 dp.
//   • Biopython 1.88 Bio.SeqUtils.IsoelectricPoint (Bjellqvist 1993/1994 pK set with N/C-terminal
//     residue-specific pKs); Bjellqvist expected values are the exact root of Biopython's
//     charge_at_pH (scipy brentq over [0,14]), which equals IsoelectricPoint.pi() whenever the root
//     lies inside Biopython's [4.05, 12] bisection window.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class SequenceStatistics_CalculateIsoelectricPoint_Tests
{
    // The method returns the correctly rounded 2-dp pI (root located to 1e-9 pH).
    private const double Exact = 1e-9;

    // Swiss-Prot LACI_ECOLI (P03023), 360 aa — the EMBOSS iep documentation usage example 1
    // (taken from the EMBOSS test database test/swiss/seq.dat).
    private const string LacIEcoli =
        "MKPVTLYDVAEYAGVSYQTVSRVVNQASHVSAKTREKVEAAMAELNYIPNRVAQQLAGKQSLLIGVATSSLALHAPSQIVAAIKSRADQLG" +
        "ASVVVSMVERSGVEACKAAVHNLLAQRVSGLIINYPLDDQDAIAVEAACTNVPALFLDVSDQTPINSIIFSHEDGTRLGVEHLVALGHQQIA" +
        "LLAGPLSSVSARLRLAGWHKYLTRNQIQPIAEREGDWSAMSGFQQTMQMLNEGIVPTAMLVANDQMALGAMRAITESGLRVGADISVVGYDD" +
        "TEDSSCYIPPLTTIKQDFRLLGQTSVDRLLQLSQGQAVKGNQLLPVSLVKRKTTLAPNTQTASPRALADSLMQLARQVSRLESGQ";

    #region EMBOSS scale (default overload)

    // EMBOSS iep documentation usage example 1: LACI_ECOLI pI = 6.8385 (reproduced by the
    // EMBOSS 6.6.0 binary; with the stale Amino 8.6 listing it would be 6.8820).
    [Test]
    public void CalculateIsoelectricPoint_LacIEcoli_MatchesEmbossDocExample()
    {
        Assert.That(SequenceStatistics.CalculateIsoelectricPoint(LacIEcoli), Is.EqualTo(6.84).Within(Exact),
            "EMBOSS iep LACI_ECOLI Isoelectric Point = 6.8385");
    }

    // Termini-only: pI = (Amino 7.5 + Carboxyl 3.6)/2 = 5.55; EMBOSS iep "A" and "AG" = 5.5500.
    [TestCase("A")]
    [TestCase("AG")]
    public void CalculateIsoelectricPoint_TerminiOnly_ReturnsEmbossMidpoint555(string seq)
    {
        Assert.That(SequenceStatistics.CalculateIsoelectricPoint(seq), Is.EqualTo(5.55).Within(Exact),
            "INV-04: termini-only pI = (7.5 + 3.6)/2 = 5.55 (EMBOSS iep 5.5500)");
    }

    // EMBOSS 6.6.0 iep outputs (4 dp) → correctly rounded 2 dp.
    [TestCase("D", 3.75)]                         // iep 3.7498
    [TestCase("E", 3.85)]                         // iep 3.8497
    [TestCase("K", 9.15)]                         // iep 9.1500
    [TestCase("R", 10.00)]                        // iep 10.0000
    [TestCase("H", 7.00)]                         // iep 7.0004
    [TestCase("C", 5.53)]                         // iep 5.5289
    [TestCase("Y", 5.55)]                         // iep 5.5494
    [TestCase("DDDD", 3.23)]                      // iep 3.2279
    [TestCase("KKKK", 11.28)]                     // iep 11.2772
    [TestCase("DDDDDDDD", 2.95)]                  // iep 2.9549
    [TestCase("ACDEFGHIKLMNPQRSTVWY", 6.97)]      // iep 6.9681
    [TestCase("FLPVLAGLTPSIVPKLVCLLTKKC", 9.57)]  // iep 9.5678
    [TestCase("DKDK", 5.90)]                      // iep 5.9023
    [TestCase("PETER", 4.26)]                     // iep 4.2577
    [TestCase("MAEGEITTFT", 3.61)]                // iep 3.6135
    public void CalculateIsoelectricPoint_Emboss_MatchesIepBinary(string seq, double expected)
    {
        Assert.That(SequenceStatistics.CalculateIsoelectricPoint(seq), Is.EqualTo(expected).Within(Exact),
            $"EMBOSS 6.6.0 iep pI of {seq}");
        Assert.That(SequenceStatistics.CalculateIsoelectricPoint(seq, SequenceStatistics.PkaScale.Emboss),
            Is.EqualTo(expected).Within(Exact), "explicit Emboss scale equals the default overload");
    }

    // RRRRRRRR: iep 13.3450 sits on a rounding boundary at 4 dp, so check within 0.01; INV-01 bounds.
    [Test]
    public void CalculateIsoelectricPoint_PolyArginine_WithinBoundsAndMatchesIep()
    {
        double pi = SequenceStatistics.CalculateIsoelectricPoint("RRRRRRRR");
        Assert.Multiple(() =>
        {
            Assert.That(pi, Is.EqualTo(13.345).Within(0.01), "EMBOSS iep 13.3450");
            Assert.That(pi, Is.InRange(0.0, 14.0), "INV-01");
        });
    }

    // EMBOSS embIepCompC splits B → D/N (5.5:4.3) and Z → E/Q (6.0:3.9) by Dayhoff frequency,
    // (int)(0.5 + n·f). iep: B 3.7498, Z 3.8497, BBBB 3.4918 (2 D), ZZZZ 3.6135 (2 E).
    [TestCase("B", 3.75)]
    [TestCase("Z", 3.85)]
    [TestCase("BBBB", 3.49)]
    [TestCase("ZZZZ", 3.61)]
    [TestCase("AB", 3.75)]
    [TestCase("AZ", 3.85)]
    public void CalculateIsoelectricPoint_Emboss_AmbiguityCodesSplitByDayhoff(string seq, double expected)
    {
        Assert.That(SequenceStatistics.CalculateIsoelectricPoint(seq), Is.EqualTo(expected).Within(Exact),
            $"EMBOSS iep pI of {seq} (B/Z Dayhoff split)");
    }

    // Non-ionizable characters (whitespace, punctuation, X) are ignored and never throw.
    // "A B!G" ≡ "AB" (iep 3.7498); "XZ" ≡ "AZ" (iep 3.8497); "X-X" ≡ termini only (5.55).
    [Test]
    public void CalculateIsoelectricPoint_NonIonizableCharacters_Ignored()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceStatistics.CalculateIsoelectricPoint("A B!G"), Is.EqualTo(3.75).Within(Exact));
            Assert.That(SequenceStatistics.CalculateIsoelectricPoint("XZ"), Is.EqualTo(3.85).Within(Exact));
            Assert.That(SequenceStatistics.CalculateIsoelectricPoint("X-X"), Is.EqualTo(5.55).Within(Exact));
        });
    }

    // Net charge vs the EMBOSS iep charge table (printed to 2 dp).
    [TestCase("FLPVLAGLTPSIVPKLVCLLTKKC", 5.0, 3.03)]
    [TestCase("FLPVLAGLTPSIVPKLVCLLTKKC", 7.0, 2.70)]
    [TestCase("FLPVLAGLTPSIVPKLVCLLTKKC", 9.0, 0.46)]
    public void CalculateNetCharge_Emboss_MatchesIepChargeTable(string seq, double pH, double expected)
    {
        Assert.That(SequenceStatistics.CalculateNetCharge(seq, pH), Is.EqualTo(expected).Within(0.005),
            $"EMBOSS iep charge of {seq} at pH {pH}");
    }

    [Test]
    public void CalculateNetCharge_LacIEcoli_MatchesEmbossDocTable()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceStatistics.CalculateNetCharge(LacIEcoli, 5.0), Is.EqualTo(7.75).Within(0.005));
            Assert.That(SequenceStatistics.CalculateNetCharge(LacIEcoli, 7.0), Is.EqualTo(-0.63).Within(0.005));
            Assert.That(SequenceStatistics.CalculateNetCharge(LacIEcoli, 9.0), Is.EqualTo(-5.99).Within(0.005));
        });
    }

    [Test]
    public void CalculateIsoelectricPoint_EmptyString_ReturnsNeutralSeven()
    {
        Assert.That(SequenceStatistics.CalculateIsoelectricPoint(""), Is.EqualTo(7.0),
            "Empty input has no defined pI; the documented input-guard sentinel is 7.0");
    }

    [Test]
    public void CalculateIsoelectricPoint_Null_ReturnsNeutralSeven()
    {
        Assert.That(SequenceStatistics.CalculateIsoelectricPoint(null!), Is.EqualTo(7.0),
            "Null input has no defined pI; the documented input-guard sentinel is 7.0");
    }

    [Test]
    public void CalculateIsoelectricPoint_LowercaseInput_MatchesUppercase()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceStatistics.CalculateIsoelectricPoint("dddd"),
                Is.EqualTo(SequenceStatistics.CalculateIsoelectricPoint("DDDD")));
            Assert.That(SequenceStatistics.CalculateIsoelectricPoint("peter", SequenceStatistics.PkaScale.Bjellqvist),
                Is.EqualTo(SequenceStatistics.CalculateIsoelectricPoint("PETER", SequenceStatistics.PkaScale.Bjellqvist)));
        });
    }

    // INV-02 (EMBOSS only): composition-only model — permutations are equal (iep DKDK = KDDK = 5.9023).
    [Test]
    public void CalculateIsoelectricPoint_Emboss_PermutedSequence_HasIdenticalPi()
    {
        Assert.That(SequenceStatistics.CalculateIsoelectricPoint("KDDK"),
            Is.EqualTo(SequenceStatistics.CalculateIsoelectricPoint("DKDK")),
            "INV-02: EMBOSS pI depends only on composition");
    }

    #endregion

    #region Bjellqvist scale (ExPASy / Biopython)

    // Biopython IsoelectricPoint(seq).pi() (Bjellqvist, terminal-residue-specific pKs).
    [TestCase("ACDEFGHIKLMNPQRSTVWY", 6.78)]      // Biopython 6.784552; seqinr computePI doc 6.78454
    [TestCase("INGAR", 9.75)]                     // Biopython doctest 9.75
    [TestCase("PETER", 4.53)]                     // Biopython doctest 4.53 (N-term P 8.36, C-term R)
    [TestCase("FLPVLAGLTPSIVPKLVCLLTKKC", 9.39)]  // Biopython 9.3902
    [TestCase("A", 5.57)]                         // (N-term A 7.59 + C-term 3.55)/2 = 5.57
    [TestCase("D", 4.30)]                         // Biopython 4.2994 (C-term D 4.55)
    [TestCase("K", 8.75)]                         // Biopython 8.7501
    [TestCase("E", 4.60)]                         // Biopython 4.5993 (N-term E 7.7, C-term E 4.75)
    [TestCase("AKD", 6.13)]                       // Biopython 6.1315
    [TestCase("SKE", 5.94)]                       // Biopython 5.9377
    [TestCase("PKD", 6.51)]                       // Biopython 6.5108
    [TestCase("KKKK", 10.48)]                     // Biopython 10.4777
    public void CalculateIsoelectricPoint_Bjellqvist_MatchesBiopython(string seq, double expected)
    {
        Assert.That(SequenceStatistics.CalculateIsoelectricPoint(seq, SequenceStatistics.PkaScale.Bjellqvist),
            Is.EqualTo(expected).Within(Exact), $"Biopython Bjellqvist pI of {seq}");
    }

    [Test]
    public void CalculateIsoelectricPoint_Bjellqvist_LacIEcoli_MatchesBiopython()
    {
        Assert.That(SequenceStatistics.CalculateIsoelectricPoint(LacIEcoli, SequenceStatistics.PkaScale.Bjellqvist),
            Is.EqualTo(6.39).Within(Exact), "Biopython IsoelectricPoint(LACI_ECOLI).pi() = 6.3901");
    }

    // Roots outside Biopython's [4.05, 12] bisection window: Biopython returns the window edge
    // (4.05 / 12.0); the true root of Biopython's own charge_at_pH (brentq) is asserted instead.
    [TestCase("MAEGEITTFT", 3.79)]  // brentq 3.794634 (Biopython pi() clamps to 4.05)
    [TestCase("DDDD", 3.52)]        // brentq 3.521692
    [TestCase("RRRRRRRR", 12.85)]   // brentq 12.845100 (Biopython pi() clamps to 12.0)
    public void CalculateIsoelectricPoint_Bjellqvist_ExtremeRoots_AreExact(string seq, double expected)
    {
        Assert.That(SequenceStatistics.CalculateIsoelectricPoint(seq, SequenceStatistics.PkaScale.Bjellqvist),
            Is.EqualTo(expected).Within(Exact), $"root of Biopython charge_at_pH for {seq}");
    }

    // Terminal-residue-specific pKs make the Bjellqvist pI order-dependent: AKD ≠ PKD ≠ DKA.
    [Test]
    public void CalculateIsoelectricPoint_Bjellqvist_DependsOnTerminalResidues()
    {
        double akd = SequenceStatistics.CalculateIsoelectricPoint("AKD", SequenceStatistics.PkaScale.Bjellqvist);
        double pkd = SequenceStatistics.CalculateIsoelectricPoint("PKD", SequenceStatistics.PkaScale.Bjellqvist);
        Assert.That(pkd, Is.Not.EqualTo(akd), "N-terminal P (8.36) vs A (7.59) shifts pI (6.51 vs 6.13)");
    }

    // Net charge vs Biopython IsoelectricPoint.charge_at_pH (6 dp).
    [TestCase("FLPVLAGLTPSIVPKLVCLLTKKC", 5.0, 3.030883)]
    [TestCase("FLPVLAGLTPSIVPKLVCLLTKKC", 7.0, 2.737303)]
    [TestCase("FLPVLAGLTPSIVPKLVCLLTKKC", 9.0, 0.757930)]
    [TestCase("INGAR", 7.0, 0.760092)]
    [TestCase("PETER", 7.0, -1.035860)]
    public void CalculateNetCharge_Bjellqvist_MatchesBiopython(string seq, double pH, double expected)
    {
        Assert.That(SequenceStatistics.CalculateNetCharge(seq, pH, SequenceStatistics.PkaScale.Bjellqvist),
            Is.EqualTo(expected).Within(5e-7), $"Biopython charge_at_pH({pH}) of {seq}");
    }

    [Test]
    public void CalculateIsoelectricPoint_Bjellqvist_EmptyOrNull_ReturnsNeutralSeven()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceStatistics.CalculateIsoelectricPoint("", SequenceStatistics.PkaScale.Bjellqvist), Is.EqualTo(7.0));
            Assert.That(SequenceStatistics.CalculateIsoelectricPoint(null!, SequenceStatistics.PkaScale.Bjellqvist), Is.EqualTo(7.0));
            Assert.That(SequenceStatistics.CalculateNetCharge("", 7.0, SequenceStatistics.PkaScale.Bjellqvist), Is.EqualTo(0.0));
        });
    }

    #endregion
}
