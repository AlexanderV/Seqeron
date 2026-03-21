using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for codon tables, translation, and codon usage.
/// Verifies invariants of the genetic code and translation process.
///
/// Test Units: TRANS-CODON-001, TRANS-PROT-001, CODON-USAGE-001, CODON-CAI-001, CODON-OPT-001
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

    #region CODON-OPT-001: P: optimized translates to same protein; R: only valid codons; D: deterministic

    /// <summary>
    /// INV-1: Optimized sequence translates to the same protein as the original.
    /// Evidence: Codon optimization replaces synonymous codons only — the amino acid
    /// sequence encoded by the DNA must be preserved exactly.
    /// This is the fundamental invariant of codon optimization.
    /// </summary>
    [Test]
    [Category("Property")]
    public void CodonOptimization_PreservesProteinSequence()
    {
        // A known coding sequence (GFP fragment)
        string codingSeq = "AUGGCUAGCAAAGGA";
        var result = CodonOptimizer.OptimizeSequence(codingSeq, CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        // Translate both original and optimized
        string originalProtein = TranslateRna(result.OriginalSequence);
        string optimizedProtein = TranslateRna(result.OptimizedSequence);

        Assert.That(optimizedProtein, Is.EqualTo(originalProtein),
            $"Optimized protein '{optimizedProtein}' ≠ original '{originalProtein}'");
    }

    /// <summary>
    /// INV-2: Optimized sequence translates to same protein for arbitrary coding DNA.
    /// Evidence: Synonymous codon substitution preserves translation (by definition of the genetic code).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CodonOptimization_PreservesProtein_Property()
    {
        return Prop.ForAll(CodingDnaArbitrary(), seq =>
        {
            if (seq.Length < 3) return true.ToProperty();
            var result = CodonOptimizer.OptimizeSequence(seq, CodonOptimizer.EColiK12);
            if (result.OptimizedSequence.Length == 0) return true.ToProperty();

            string origProtein = TranslateRna(result.OriginalSequence);
            string optProtein = TranslateRna(result.OptimizedSequence);
            return (origProtein == optProtein)
                .Label($"Protein mismatch: orig='{origProtein[..Math.Min(10, origProtein.Length)]}' ≠ opt='{optProtein[..Math.Min(10, optProtein.Length)]}'");
        });
    }

    /// <summary>
    /// INV-3: Optimized sequence contains only valid RNA codons (A, C, G, U).
    /// Evidence: CodonOptimizer works with RNA codons internally.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CodonOptimization_OnlyValidCodons()
    {
        return Prop.ForAll(CodingDnaArbitrary(), seq =>
        {
            if (seq.Length < 3) return true.ToProperty();
            var result = CodonOptimizer.OptimizeSequence(seq, CodonOptimizer.EColiK12);
            if (result.OptimizedSequence.Length == 0) return true.ToProperty();

            bool valid = result.OptimizedSequence.All(c => "ACGU".Contains(c));
            return valid.Label($"Invalid chars in optimized: {result.OptimizedSequence[..Math.Min(20, result.OptimizedSequence.Length)]}");
        });
    }

    /// <summary>
    /// INV-4: Optimized sequence length equals original sequence length.
    /// Evidence: Codon optimization replaces codons 1:1, preserving total length.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CodonOptimization_PreservesLength()
    {
        return Prop.ForAll(CodingDnaArbitrary(), seq =>
        {
            if (seq.Length < 3) return true.ToProperty();
            var result = CodonOptimizer.OptimizeSequence(seq, CodonOptimizer.EColiK12);
            return (result.OptimizedSequence.Length == result.OriginalSequence.Length)
                .Label($"Length changed: orig={result.OriginalSequence.Length}, opt={result.OptimizedSequence.Length}");
        });
    }

    /// <summary>
    /// INV-5: Codon optimization is deterministic.
    /// Evidence: OptimizeSequence is a pure function.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CodonOptimization_IsDeterministic()
    {
        return Prop.ForAll(CodingDnaArbitrary(), seq =>
        {
            if (seq.Length < 3) return true.ToProperty();
            var r1 = CodonOptimizer.OptimizeSequence(seq, CodonOptimizer.EColiK12);
            var r2 = CodonOptimizer.OptimizeSequence(seq, CodonOptimizer.EColiK12);
            return (r1.OptimizedSequence == r2.OptimizedSequence &&
                    r1.OptimizedCAI == r2.OptimizedCAI)
                .Label("OptimizeSequence must be deterministic");
        });
    }

    #endregion

    private static string TranslateRna(string rna)
    {
        var sb = new System.Text.StringBuilder();
        for (int i = 0; i + 2 < rna.Length; i += 3)
        {
            string codon = rna.Substring(i, 3);
            // Convert to DNA for GeneticCode
            string dnaCodon = codon.Replace('U', 'T');
            sb.Append(GeneticCode.Standard.Translate(dnaCodon));
        }
        return sb.ToString();
    }
}
