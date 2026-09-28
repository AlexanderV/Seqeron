namespace Seqeron.Genomics.Tests.Unit.Core;

/// <summary>
/// Review 2026-09, B02 duplication sweep (F28–F30): <see cref="ProteinSequence"/> molecular weight,
/// isoelectric point and GRAVY delegate to the canonical <see cref="ProteinPhysicochemistry"/>
/// (the implementation of SequenceStatistics SEQ-MW-001 / SEQ-PI-001 / SEQ-HYDRO-001), and
/// <see cref="ProteinSequence.FindMotif"/> uses the suffix-tree exact matcher.
/// Reference values: Biopython 1.88 (molecular_weight, ProteinAnalysis.gravy, IsoelectricPoint)
/// and the EMBOSS iep documentation example (LACI_ECOLI pI 6.8385).
/// </summary>
[TestFixture]
public class ProteinSequence_CanonicalProperties_Tests
{
    private const double Tolerance = 1e-9;

    // Swiss-Prot LACI_ECOLI (P03023), 360 aa — EMBOSS iep documentation usage example 1.
    private const string LacIEcoli =
        "MKPVTLYDVAEYAGVSYQTVSRVVNQASHVSAKTREKVEAAMAELNYIPNRVAQQLAGKQSLLIGVATSSLALHAPSQIVAAIKSRADQLG" +
        "ASVVVSMVERSGVEACKAAVHNLLAQRVSGLIINYPLDDQDAIAVEAACTNVPALFLDVSDQTPINSIIFSHEDGTRLGVEHLVALGHQQIA" +
        "LLAGPLSSVSARLRLAGWHKYLTRNQIQPIAEREGDWSAMSGFQQTMQMLNEGIVPTAMLVANDQMALGAMRAITESGLRVGADISVVGYDD" +
        "TEDSSCYIPPLTTIKQDFRLLGQTSVDRLLQLSQGQAVKGNQLLPVSLVKRKTTLAPNTQTASPRALADSLMQLARQVSRLESGQ";

    // Human insulin A chain.
    private const string InsulinA = "GIVEQCCTSICSLYQLENYCN";

    #region F28 — molecular weight (Biopython masses, water 18.0153, unrounded)

    [TestCase("MK", 277.3836)]                       // Biopython 1.88: 277.3836 (old code 277.38)
    [TestCase(InsulinA, 2383.6961000000006)]         // Biopython 1.88 (old code 2383.72)
    [TestCase(LacIEcoli, 38589.71700000004)]         // Biopython 1.88 (old code 38589.84)
    [TestCase("GGGAAAVVV", 699.7959000000001)]       // Biopython 1.88 (old code 699.81)
    public void MolecularWeight_MatchesBiopython(string sequence, double expected)
    {
        Assert.That(new ProteinSequence(sequence).MolecularWeight(), Is.EqualTo(expected).Within(1e-6));
    }

    [Test]
    public void MolecularWeight_AmbiguousResiduesCarryNoMass_AsCanonical()
    {
        // B/Z/J/X/'*' have no mass (skipped), exactly as SequenceStatistics.CalculateMolecularWeight.
        var protein = new ProteinSequence("MBZJX*K");
        Assert.That(protein.MolecularWeight(), Is.EqualTo(277.3836).Within(1e-9));
    }

    [Test]
    public void Properties_MolecularWeight_AreBiopythonProteinWeights()
    {
        // Biopython Bio/Data/IUPACData.py protein_weights.
        Assert.Multiple(() =>
        {
            Assert.That(ProteinSequence.Properties['A'].MolecularWeight, Is.EqualTo(89.0932));
            Assert.That(ProteinSequence.Properties['W'].MolecularWeight, Is.EqualTo(204.2252));
            Assert.That(ProteinSequence.Properties['R'].MolecularWeight, Is.EqualTo(174.201));
            foreach (var (aa, props) in ProteinSequence.Properties)
                Assert.That(props.MolecularWeight, Is.EqualTo(ProteinPhysicochemistry.AverageAminoAcidMasses[aa]));
        });
    }

    #endregion

    #region F29 — isoelectric point (EMBOSS iep pK set)

    [Test]
    public void IsoelectricPoint_LacIEcoli_MatchesEmbossIepExample()
    {
        // EMBOSS iep documentation example 1: LACI_ECOLI Isoelectric Point = 6.8385.
        // The previous unsourced pK mix (N-term 9.69, C-term 2.34, H 6.0, K 10.5, C 8.3) gave 6.42.
        Assert.That(new ProteinSequence(LacIEcoli).IsoelectricPoint(), Is.EqualTo(6.84).Within(Tolerance));
    }

    [TestCase("A", 5.55)]         // termini only: (7.5 + 3.6)/2 — EMBOSS iep "A" = 5.5500 (old code 6.01)
    [TestCase("GGGAAAVVV", 5.55)] // no ionizable side chain (old code 6.01)
    [TestCase("MK", 9.15)]        // (old code 10.10)
    [TestCase("DDDDEEEE", 3.03)]  // (old code 2.63)
    public void IsoelectricPoint_EmbossScale(string sequence, double expected)
    {
        Assert.That(new ProteinSequence(sequence).IsoelectricPoint(), Is.EqualTo(expected).Within(Tolerance));
    }

    #endregion

    #region F30 — GRAVY unrounded (Biopython ProteinAnalysis.gravy)

    [TestCase(InsulinA, 0.21428571428571427)]     // Biopython 1.88 (old code rounded to 0.214)
    [TestCase(LacIEcoli, 0.03166666666666661)]    // Biopython 1.88 (old code 0.032)
    [TestCase("GGGAAAVVV", 1.8666666666666667)]   // Biopython 1.88 (old code 1.867)
    public void Gravy_MatchesBiopythonUnrounded(string sequence, double expected)
    {
        Assert.That(new ProteinSequence(sequence).Gravy(), Is.EqualTo(expected).Within(1e-12));
    }

    #endregion

    #region Single implementation: ProteinSequence == ProteinPhysicochemistry == SequenceStatistics

    [Test]
    public void ProteinProperties_EqualSequenceStatistics_OnRandomProteins()
    {
        const string alphabet = "ACDEFGHIKLMNPQRSTVWYBZJX*";
        var random = new Random(20260928);
        for (int n = 0; n < 500; n++)
        {
            var chars = new char[random.Next(1, 120)];
            for (int i = 0; i < chars.Length; i++)
                chars[i] = alphabet[random.Next(alphabet.Length)];
            string seq = new(chars);
            var protein = new ProteinSequence(seq);

            Assert.That(protein.MolecularWeight(), Is.EqualTo(SequenceStatistics.CalculateMolecularWeight(seq)), seq);
            Assert.That(protein.IsoelectricPoint(), Is.EqualTo(SequenceStatistics.CalculateIsoelectricPoint(seq)), seq);
            Assert.That(protein.Gravy(), Is.EqualTo(SequenceStatistics.CalculateHydrophobicity(seq)), seq);
            Assert.That(ProteinPhysicochemistry.IsoelectricPoint(seq, ProteinPkaScale.Bjellqvist),
                Is.EqualTo(SequenceStatistics.CalculateIsoelectricPoint(seq, SequenceStatistics.PkaScale.Bjellqvist)), seq);
            Assert.That(ProteinPhysicochemistry.NetCharge(seq, 7.0),
                Is.EqualTo(SequenceStatistics.CalculateNetCharge(seq, 7.0)), seq);
        }
    }

    [Test]
    public void ProteinPhysicochemistry_Bjellqvist_LacIEcoli_MatchesBiopython()
    {
        // Biopython 1.88 IsoelectricPoint(LACI_ECOLI).pi() = 6.3901.
        Assert.That(ProteinPhysicochemistry.IsoelectricPoint(LacIEcoli, ProteinPkaScale.Bjellqvist),
            Is.EqualTo(6.39).Within(Tolerance));
    }

    [Test]
    public void ProteinPhysicochemistry_EmptyInput_Sentinels()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ProteinPhysicochemistry.MolecularWeight(""), Is.EqualTo(0));
            Assert.That(ProteinPhysicochemistry.Gravy(null), Is.EqualTo(0));
            Assert.That(ProteinPhysicochemistry.IsoelectricPoint(""), Is.EqualTo(7.0));
            Assert.That(new ProteinSequence("").IsoelectricPoint(), Is.EqualTo(0)); // ProteinSequence contract
        });
    }

    #endregion

    #region FindMotif — suffix-tree exact matching

    [Test]
    public void FindMotif_OverlappingOccurrences_AscendingOrder()
    {
        var protein = new ProteinSequence("AAAAKAAA");
        Assert.That(protein.FindMotif("aa").ToList(), Is.EqualTo(new[] { 0, 1, 2, 5, 6 }));
    }

    [Test]
    public void FindMotif_PatternLongerThanSequence_ReturnsEmpty()
    {
        Assert.Multiple(() =>
        {
            Assert.That(new ProteinSequence("MK").FindMotif("MKV"), Is.Empty);
            Assert.That(new ProteinSequence("").FindMotif("M"), Is.Empty);
        });
    }

    #endregion
}
