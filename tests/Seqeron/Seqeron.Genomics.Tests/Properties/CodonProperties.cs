using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for codon tables, translation, and codon usage.
/// Verifies invariants of the genetic code and translation process.
///
/// Test Units: TRANS-CODON-001, TRANS-PROT-001, CODON-USAGE-001, CODON-CAI-001 (Property Extensions)
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Translation")]
public class CodonProperties
{
    private static Arbitrary<string> CodingDnaArbitrary() =>
        Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length >= 3)
            .Select(a => new string(a, 0, a.Length - a.Length % 3)) // trim to codon boundary
            .Where(s => s.Length >= 3)
            .ToArbitrary();

    /// <summary>
    /// The standard genetic code must have exactly 64 codons (4³).
    /// </summary>
    [Test]
    [Category("Property")]
    public void StandardCode_Has64Codons()
    {
        var code = GeneticCode.Standard;
        Assert.That(code.CodonTable.Count, Is.EqualTo(64),
            "Standard genetic code must map all 64 codons");
    }

    /// <summary>
    /// Every codon in the table must be exactly 3 characters.
    /// </summary>
    [Test]
    [Category("Property")]
    public void AllCodons_HaveLength3()
    {
        var code = GeneticCode.Standard;
        Assert.That(code.CodonTable.Keys.All(c => c.Length == 3), Is.True,
            "All codons must be exactly 3 nucleotides");
    }

    /// <summary>
    /// Start and stop codon sets must be non-empty and disjoint.
    /// </summary>
    [Test]
    [Category("Property")]
    public void StartAndStopCodons_AreDisjoint()
    {
        var code = GeneticCode.Standard;
        Assert.That(code.StartCodons.Count, Is.GreaterThan(0));
        Assert.That(code.StopCodons.Count, Is.GreaterThan(0));
        Assert.That(code.StartCodons.Intersect(code.StopCodons).Count(), Is.EqualTo(0),
            "Start and stop codon sets must be disjoint");
    }

    /// <summary>
    /// ATG must be a start codon in the standard code.
    /// </summary>
    [Test]
    [Category("Property")]
    public void ATG_IsStartCodon()
    {
        Assert.That(GeneticCode.Standard.IsStartCodon("ATG"), Is.True);
    }

    /// <summary>
    /// TAA, TAG, TGA must be stop codons in the standard code.
    /// </summary>
    [TestCase("TAA")]
    [TestCase("TAG")]
    [TestCase("TGA")]
    [Category("Property")]
    public void KnownStopCodons_AreStopCodons(string codon)
    {
        Assert.That(GeneticCode.Standard.IsStopCodon(codon), Is.True);
    }

    /// <summary>
    /// Translation preserves length: protein length == DNA length / 3.
    /// (Excluding stop codons which produce '*')
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Translation_OutputLength_MatchesInput()
    {
        return Prop.ForAll(CodingDnaArbitrary(), seq =>
        {
            if (seq.Length < 3) return true.ToProperty();
            var protein = Translator.Translate(seq);
            int expectedLen = seq.Length / 3;
            return (protein.Sequence.Length == expectedLen)
                .Label($"Protein len={protein.Sequence.Length}, expected={expectedLen}");
        });
    }

    /// <summary>
    /// Translation output only contains valid amino acid characters or '*' (stop).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Translation_ProducesValidAminoAcids()
    {
        const string validAa = "ACDEFGHIKLMNPQRSTVWY*";
        return Prop.ForAll(CodingDnaArbitrary(), seq =>
        {
            if (seq.Length < 3) return true.ToProperty();
            var protein = Translator.Translate(seq);
            return protein.Sequence.All(c => validAa.Contains(c))
                .Label($"Invalid amino acid in: {protein.Sequence[..Math.Min(20, protein.Sequence.Length)]}");
        });
    }

    /// <summary>
    /// Codon usage counts sum to total number of codons in the sequence.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CodonUsage_CountsSumToTotal()
    {
        return Prop.ForAll(CodingDnaArbitrary(), seq =>
        {
            if (seq.Length < 3) return true.ToProperty();
            var usage = CodonOptimizer.CalculateCodonUsage(seq);
            int total = usage.Values.Sum();
            int expected = seq.Length / 3;
            return (total == expected)
                .Label($"Usage sum={total}, expected={expected}");
        });
    }

    /// <summary>
    /// CAI is always in (0, 1] for valid coding sequences.
    /// Evidence: CAI = geometric mean of relative adaptiveness values, each in (0, 1].
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CAI_InRange()
    {
        // Generate coding sequences that are multiples of 3, contain only ACGT, 
        // and have at least one non-stop codon
        var codingGen = Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length >= 3)
            .Select(a => new string(a, 0, a.Length - a.Length % 3))
            .Where(s => s.Length >= 3)
            // Filter out sequences that are all stop codons
            .Where(s =>
            {
                var code = GeneticCode.Standard;
                for (int i = 0; i < s.Length - 2; i += 3)
                {
                    string codon = s.Substring(i, 3);
                    if (!code.IsStopCodon(codon)) return true;
                }
                return false;
            })
            .ToArbitrary();

        return Prop.ForAll(codingGen, seq =>
        {
            double cai = CodonOptimizer.CalculateCAI(seq, CodonOptimizer.EColiK12);
            return (cai >= 0.0 && cai <= 1.0 + 0.0001)
                .Label($"CAI={cai:F4}, must be in [0, 1]");
        });
    }

    /// <summary>
    /// GetCodonsForAminoAcid returns codons that all translate back to the same amino acid.
    /// </summary>
    [TestCase('M')]
    [TestCase('W')]
    [TestCase('L')]
    [TestCase('S')]
    [Category("Property")]
    public void Codons_TranslateBackToSameAminoAcid(char aminoAcid)
    {
        var code = GeneticCode.Standard;
        var codons = code.GetCodonsForAminoAcid(aminoAcid).ToList();
        Assert.That(codons, Is.Not.Empty, $"Must have codons for '{aminoAcid}'");
        foreach (var codon in codons)
            Assert.That(code.Translate(codon), Is.EqualTo(aminoAcid),
                $"Codon {codon} should translate to {aminoAcid}");
    }
}
