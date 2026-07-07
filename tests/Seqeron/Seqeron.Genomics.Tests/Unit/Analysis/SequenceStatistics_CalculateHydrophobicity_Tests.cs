// SEQ-HYDRO-001 — Hydrophobicity Analysis (Kyte-Doolittle GRAVY + sliding-window profile)
// Evidence: docs/Evidence/SEQ-HYDRO-001-Evidence.md
// TestSpec: tests/TestSpecs/SEQ-HYDRO-001.md
// Source: Kyte J, Doolittle RF (1982). J Mol Biol 157(1):105-132.
//         Biopython Bio.SeqUtils.ProtParamData.kd (scale); ProtParam.gravy / protein_scale.
//         Expasy ProtParam doc: GRAVY = sum(hydropathy)/residue count.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class SequenceStatistics_CalculateHydrophobicity_Tests
{
    // Expected values are exact derivations from the Kyte-Doolittle kd scale.
    private const double Tolerance = 1e-10;

    #region CalculateHydrophobicity (GRAVY)

    // M1 — single residue A: GRAVY = 1.8 / 1
    // Evidence: kd['A']=1.8 (Biopython); GRAVY = sum/length (Expasy).
    [Test]
    public void CalculateHydrophobicity_SingleResidueA_ReturnsKdValue()
    {
        double gravy = SequenceStatistics.CalculateHydrophobicity("A");

        Assert.That(gravy, Is.EqualTo(1.8).Within(Tolerance),
            "GRAVY of one residue equals its kd value (1.8 for Ala) (INV-01)");
    }

    // M2 — hydrophobic FLIV: (2.8+3.8+4.5+4.2)/4 = 15.3/4 = 3.825
    // Evidence: kd F=2.8, L=3.8, I=4.5, V=4.2 (Biopython); GRAVY=sum/length (Expasy).
    [Test]
    public void CalculateHydrophobicity_HydrophobicPeptideFLIV_ReturnsExactGravy()
    {
        double gravy = SequenceStatistics.CalculateHydrophobicity("FLIV");

        Assert.That(gravy, Is.EqualTo(3.825).Within(Tolerance),
            "Sum of kd values 15.3 over 4 residues = 3.825 (positive => hydrophobic)");
    }

    // M3 — hydrophilic RKDE: (-4.5-3.9-3.5-3.5)/4 = -15.4/4 = -3.85
    // Evidence: kd R=-4.5, K=-3.9, D=-3.5, E=-3.5 (Biopython).
    [Test]
    public void CalculateHydrophobicity_HydrophilicPeptideRKDE_ReturnsExactGravy()
    {
        double gravy = SequenceStatistics.CalculateHydrophobicity("RKDE");

        Assert.That(gravy, Is.EqualTo(-3.85).Within(Tolerance),
            "Sum of kd values -15.4 over 4 residues = -3.85 (negative => hydrophilic)");
    }

    // M4 — case-insensitive: "fliv" equals "FLIV"
    // Evidence: scale defined on uppercase; impl uppercases input (INV-04).
    [Test]
    public void CalculateHydrophobicity_LowercaseInput_MatchesUppercase()
    {
        double gravy = SequenceStatistics.CalculateHydrophobicity("fliv");

        Assert.That(gravy, Is.EqualTo(3.825).Within(Tolerance),
            "GRAVY is case-insensitive; lowercase 'fliv' equals 'FLIV' = 3.825 (INV-04)");
    }

    // M5 — empty string -> 0
    [Test]
    public void CalculateHydrophobicity_EmptyString_ReturnsZero()
    {
        double gravy = SequenceStatistics.CalculateHydrophobicity("");

        Assert.That(gravy, Is.EqualTo(0),
            "Empty sequence has no residues to average, GRAVY = 0 (INV-05)");
    }

    // M6 — null -> 0
    [Test]
    public void CalculateHydrophobicity_Null_ReturnsZero()
    {
        double gravy = SequenceStatistics.CalculateHydrophobicity(null!);

        Assert.That(gravy, Is.EqualTo(0),
            "Null input returns 0 without throwing (INV-05)");
    }

    // S1 — unknown residues are skipped: "AX" divides by recognized count (1) -> 1.8
    // Evidence: scale defines only the 20 standard residues; impl skips unknowns (deviation 5.4).
    [Test]
    public void CalculateHydrophobicity_UnknownResidueSkipped_DividesByRecognizedCount()
    {
        double gravy = SequenceStatistics.CalculateHydrophobicity("AX");

        Assert.That(gravy, Is.EqualTo(1.8).Within(Tolerance),
            "'X' is not in the kd scale and is skipped; GRAVY = 1.8/1 = 1.8 (deviation 5.4)");
    }

    // S1b — all non-standard residues: recognized count is 0 -> GRAVY 0 (not NaN)
    // Evidence: scale defines only the 20 standard residues (Kyte-Doolittle 1982 / Expasy);
    //           impl returns 0 when no residue is recognized (contract 6.1).
    [Test]
    public void CalculateHydrophobicity_AllUnknownResidues_ReturnsZero()
    {
        double gravy = SequenceStatistics.CalculateHydrophobicity("XXZB");

        Assert.That(gravy, Is.EqualTo(0),
            "No residue is in the kd scale; recognized count is 0, GRAVY = 0 (not NaN) (contract 6.1)");
    }

    #endregion

    #region CalculateHydrophobicityProfile (sliding window)

    // M7 — FLIV, window 3: two windows, exact means
    // Evidence: protein_scale yields n-W+1 unweighted means (Biopython, edge=1.0).
    [Test]
    public void CalculateHydrophobicityProfile_FLIV_Window3_ReturnsExactWindowMeans()
    {
        var profile = SequenceStatistics.CalculateHydrophobicityProfile("FLIV", 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(profile, Has.Count.EqualTo(2),
                "n - W + 1 = 4 - 3 + 1 = 2 windows (INV-02)");
            Assert.That(profile[0], Is.EqualTo((2.8 + 3.8 + 4.5) / 3).Within(Tolerance),
                "Window 1 mean of F,L,I = 11.1/3 = 3.7 (unweighted, INV-03)");
            Assert.That(profile[1], Is.EqualTo((3.8 + 4.5 + 4.2) / 3).Within(Tolerance),
                "Window 2 mean of L,I,V = 12.5/3 = 4.16666666667 (unweighted, INV-03)");
        });
    }

    // M8 — window larger than sequence -> empty
    // Evidence: range(n-W+1) yields 0 iterations when W>n (Biopython).
    [Test]
    public void CalculateHydrophobicityProfile_WindowLargerThanSequence_ReturnsEmpty()
    {
        var profile = SequenceStatistics.CalculateHydrophobicityProfile("AG", 3).ToList();

        Assert.That(profile, Is.Empty,
            "Window 3 > length 2 produces no windows (INV-02)");
    }

    // M9 — empty and null -> empty profile
    [Test]
    public void CalculateHydrophobicityProfile_EmptyOrNull_ReturnsEmpty()
    {
        var fromEmpty = SequenceStatistics.CalculateHydrophobicityProfile("", 3).ToList();
        var fromNull = SequenceStatistics.CalculateHydrophobicityProfile(null!, 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(fromEmpty, Is.Empty, "Empty sequence yields an empty profile (INV-05)");
            Assert.That(fromNull, Is.Empty, "Null sequence yields an empty profile (INV-05)");
        });
    }

    // S2 — profile length is n-W+1 for a longer sequence
    // Evidence: protein_scale loop bound (Biopython).
    [Test]
    public void CalculateHydrophobicityProfile_LongerSequence_HasNMinusWPlus1Values()
    {
        var profile = SequenceStatistics.CalculateHydrophobicityProfile("FLIVAG", 3).ToList();

        Assert.That(profile, Has.Count.EqualTo(4),
            "n - W + 1 = 6 - 3 + 1 = 4 windows (INV-02)");
    }

    // S3 — profile divides by window size (W), so a non-standard residue contributes 0
    // to its window's sum (still divided by W, not by recognized count).
    // Evidence: contract 5.2/6.1 deviation; window mean is sum/W with unknowns adding 0.
    [Test]
    public void CalculateHydrophobicityProfile_UnknownResidueInWindow_ContributesZero()
    {
        var profile = SequenceStatistics.CalculateHydrophobicityProfile("FXIV", 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(profile, Has.Count.EqualTo(2),
                "n - W + 1 = 4 - 3 + 1 = 2 windows (INV-02)");
            Assert.That(profile[0], Is.EqualTo((2.8 + 0.0 + 4.5) / 3).Within(Tolerance),
                "Window F,X,I: X is non-standard => 0; sum 7.3 divided by W=3 (deviation 5.4)");
            Assert.That(profile[1], Is.EqualTo((0.0 + 4.5 + 4.2) / 3).Within(Tolerance),
                "Window X,I,V: X is non-standard => 0; sum 8.7 divided by W=3 (deviation 5.4)");
        });
    }

    // C1 — transmembrane-style window (W=19) over a hydrophobic stretch exceeds 1.6
    // Evidence: kd I=4.5; GCAT/Kyte-Doolittle: transmembrane peaks > 1.6 at window 19.
    [Test]
    public void CalculateHydrophobicityProfile_HydrophobicStretchWindow19_ExceedsTransmembraneThreshold()
    {
        const string isoleucine19 = "IIIIIIIIIIIIIIIIIII"; // 19 residues, all kd 4.5

        var profile = SequenceStatistics.CalculateHydrophobicityProfile(isoleucine19, 19).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(profile, Has.Count.EqualTo(1),
                "n - W + 1 = 19 - 19 + 1 = 1 window (INV-02)");
            Assert.That(profile[0], Is.EqualTo(4.5).Within(Tolerance),
                "Uniform window of Ile (kd 4.5) averages to 4.5, above the 1.6 transmembrane threshold");
        });
    }

    #endregion
}
