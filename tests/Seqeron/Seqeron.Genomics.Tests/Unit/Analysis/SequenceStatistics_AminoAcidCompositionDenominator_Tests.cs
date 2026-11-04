// SEQ-AA-COMP-001 — amino-acid composition ratio denominator (review 2026-09, FIN-B03 WP10, audit row A17).
// Evidence: docs/Validation/review-2026-09/B03.md, F35.
// Reference runs (2026-10-14):
//   Biopython 1.88 ProteinAnalysis (self.length = len(seq.upper()); amino_acids_percent = count*100/length;
//   aromaticity = Σ percent[Y,W,F]/100):
//     "MKW"  -> M,K,W 33.333333333333336 %, aromaticity 0.33333333333333337
//     "MKW*" -> M,K,W 25.0 %, aromaticity 0.25
//     "MKW-A"-> M,K,W,A 20.0 %, aromaticity 0.2
//     "mkwfy**" -> M,K,W,F,Y 14.285714285714286 %, aromaticity 0.4285714285714286
//     "FWY12K"  -> 16.666666666666668 % each, aromaticity 0.5
//     "ACDEFGHIKLMNPQRSTVWYX" -> 4.761904761904762 % each, aromaticity 0.14285714285714285
//   EMBOSS 6.6.0 pepstats Mole% (denominator = residues): "MKW*" -> Residues = 3, K 33.333;
//     "MKW-A" -> Residues = 4, K 25.000; "MKWX" -> Residues = 4, K 25.000; "MKWB" -> D(Asx row) 25.000.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class SequenceStatistics_AminoAcidCompositionDenominator_Tests
{
    private const double Tolerance = 1e-12;
    private const SequenceStatistics.AminoAcidCompositionDenominator Residues =
        SequenceStatistics.AminoAcidCompositionDenominator.Residues;
    private const SequenceStatistics.AminoAcidCompositionDenominator SeqLength =
        SequenceStatistics.AminoAcidCompositionDenominator.SequenceLength;

    [TestCase("MKW", 0.33333333333333337)]
    [TestCase("MKW*", 0.25)]
    [TestCase("MKW-A", 0.2)]
    [TestCase("mkwfy**", 0.4285714285714286)]
    [TestCase("FWY12K", 0.5)]
    [TestCase("ACDEFGHIKLMNPQRSTVWYX", 0.14285714285714285)]
    public void Aromaticity_SequenceLengthDenominator_MatchesBiopython(string sequence, double expected)
    {
        var comp = SequenceStatistics.CalculateAminoAcidComposition(sequence, SeqLength);
        Assert.That(comp.AromaticResidueRatio, Is.EqualTo(expected).Within(Tolerance));
    }

    [TestCase("MKW*", 'K', 25.0)]
    [TestCase("MKW-A", 'A', 20.0)]
    [TestCase("mkwfy**", 'F', 14.285714285714286)]
    [TestCase("FWY12K", 'K', 16.666666666666668)]
    [TestCase("ACDEFGHIKLMNPQRSTVWYX", 'Q', 4.761904761904762)]
    public void Fractions_SequenceLengthDenominator_MatchBiopythonAminoAcidsPercent(string sequence, char residue, double percent)
    {
        var fractions = SequenceStatistics.CalculateAminoAcidFractions(sequence, SeqLength);
        Assert.That(fractions[residue] * 100, Is.EqualTo(percent).Within(Tolerance));
    }

    [TestCase("MKW*", 'K', 33.333)]
    [TestCase("MKW-A", 'K', 25.000)]
    [TestCase("MKWX", 'K', 25.000)]
    public void Fractions_ResidueDenominator_MatchEmbossPepstatsMolePercent(string sequence, char residue, double molePercent)
    {
        var fractions = SequenceStatistics.CalculateAminoAcidFractions(sequence);
        Assert.That(Math.Round(fractions[residue] * 100, 3), Is.EqualTo(molePercent));
    }

    [Test]
    public void Fractions_KeysAreResidueLetters_AmbiguityCodesIncluded()
    {
        var fractions = SequenceStatistics.CalculateAminoAcidFractions("mKwB*");
        Assert.Multiple(() =>
        {
            Assert.That(fractions.Keys, Is.EquivalentTo(new[] { 'M', 'K', 'W', 'B' }));
            Assert.That(fractions.Values.Sum(), Is.EqualTo(1.0).Within(Tolerance));
            Assert.That(SequenceStatistics.CalculateAminoAcidFractions("mKwB*", SeqLength)['B'], Is.EqualTo(0.2).Within(Tolerance));
            Assert.That(SequenceStatistics.CalculateAminoAcidFractions(""), Is.Empty);
            Assert.That(SequenceStatistics.CalculateAminoAcidFractions(null!, SeqLength), Is.Empty);
        });
    }

    [Test]
    public void Composition_SequenceLengthDenominator_OnlyRatiosChange()
    {
        var residues = SequenceStatistics.CalculateAminoAcidComposition("MKW*", Residues);
        var full = SequenceStatistics.CalculateAminoAcidComposition("MKW*", SeqLength);
        Assert.Multiple(() =>
        {
            Assert.That(residues.AromaticResidueRatio, Is.EqualTo(1.0 / 3).Within(Tolerance));
            Assert.That(residues.ChargedResidueRatio, Is.EqualTo(1.0 / 3).Within(Tolerance));
            Assert.That(full.ChargedResidueRatio, Is.EqualTo(0.25).Within(Tolerance));
            Assert.That(full.Length, Is.EqualTo(3), "Length stays the residue count");
            Assert.That(full.Counts, Is.EquivalentTo(residues.Counts));
            Assert.That(full.MolecularWeight, Is.EqualTo(residues.MolecularWeight));
            Assert.That(full.IsoelectricPoint, Is.EqualTo(residues.IsoelectricPoint));
            Assert.That(full.Hydrophobicity, Is.EqualTo(residues.Hydrophobicity));
            Assert.That(SequenceStatistics.CalculateAminoAcidComposition("**", SeqLength).AromaticResidueRatio, Is.EqualTo(0));
        });
    }

    [Test]
    public void Composition_DefaultDenominator_EqualsOneArgumentOverload_Random()
    {
        var rng = new Random(20261014);
        const string alphabet = "ACDEFGHIKLMNPQRSTVWYacdefwyBXZ*-1 ſ";
        for (int iteration = 0; iteration < 500; iteration++)
        {
            var chars = new char[rng.Next(0, 40)];
            for (int i = 0; i < chars.Length; i++) chars[i] = alphabet[rng.Next(alphabet.Length)];
            string seq = new(chars);

            var a = SequenceStatistics.CalculateAminoAcidComposition(seq);
            var b = SequenceStatistics.CalculateAminoAcidComposition(seq, Residues);
            Assert.That(b.Length, Is.EqualTo(a.Length), seq);
            Assert.That(b.Counts, Is.EquivalentTo(a.Counts), seq);
            Assert.That(b.ChargedResidueRatio, Is.EqualTo(a.ChargedResidueRatio), seq);
            Assert.That(b.AromaticResidueRatio, Is.EqualTo(a.AromaticResidueRatio), seq);

            var fractions = SequenceStatistics.CalculateAminoAcidFractions(seq, SeqLength);
            foreach (var (residue, count) in a.Counts)
                Assert.That(fractions[residue], Is.EqualTo((double)count / seq.Length).Within(Tolerance), seq);
        }
    }

    [Test]
    public void UndefinedDenominator_Throws()
    {
        var bad = (SequenceStatistics.AminoAcidCompositionDenominator)7;
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceStatistics.CalculateAminoAcidComposition("MKW", bad));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceStatistics.CalculateAminoAcidFractions("MKW", bad));
        });
    }
}
