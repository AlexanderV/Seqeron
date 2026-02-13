using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for sequence composition: validation, complexity, entropy.
///
/// Test Units: SEQ-VALID-001, SEQ-COMPLEX-001, SEQ-ENTROPY-001 (Property Extensions)
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Composition")]
public class SequenceCompositionProperties
{
    private static Arbitrary<string> DnaStringArbitrary() =>
        Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length > 0)
            .Select(a => new string(a))
            .ToArbitrary();

    // -- SEQ-VALID-001 --

    /// <summary>
    /// Pure ACGT sequences are always valid DNA.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PureAcgt_IsValidDna()
    {
        return Prop.ForAll(
            DnaStringArbitrary(),
            seq => seq.AsSpan().IsValidDna().Label("Pure ACGT must be valid DNA"));
    }

    /// <summary>
    /// Pure ACGU sequences are always valid RNA.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PureAcgu_IsValidRna()
    {
        return Prop.ForAll(
            Gen.Elements('A', 'C', 'G', 'U')
                .ArrayOf()
                .Where(a => a.Length > 0)
                .Select(a => new string(a))
                .ToArbitrary(),
            seq => seq.AsSpan().IsValidRna().Label("Pure ACGU must be valid RNA"));
    }

    /// <summary>
    /// Sequences containing non-nucleotide chars (digits, Z) are invalid DNA.
    /// </summary>
    [TestCase("ACGT1ACGT")]
    [TestCase("ZZZZZ")]
    [TestCase("ACG TACGT")]
    [Category("Property")]
    public void InvalidChars_NotValidDna(string seq)
    {
        Assert.That(seq.AsSpan().IsValidDna(), Is.False);
    }

    // -- SEQ-COMPLEX-001 --

    /// <summary>
    /// Linguistic complexity of a homopolymer is minimal (near 0).
    /// </summary>
    [TestCase("AAAAAAAAAA")]
    [TestCase("CCCCCCCCCC")]
    [Category("Property")]
    public void Homopolymer_HasLowComplexity(string seq)
    {
        var dna = new DnaSequence(seq);
        double complexity = SequenceComplexity.CalculateLinguisticComplexity(dna);
        Assert.That(complexity, Is.LessThan(0.25), $"Complexity={complexity:F4} for homopolymer");
    }

    /// <summary>
    /// Complexity of a random-looking sequence is high.
    /// </summary>
    [Test]
    [Category("Property")]
    public void ComplexSequence_HasHighComplexity()
    {
        var dna = new DnaSequence("ACGTGCTAGCATGCATGCATGCTAGCTAGCATG");
        double complexity = SequenceComplexity.CalculateLinguisticComplexity(dna);
        Assert.That(complexity, Is.GreaterThan(0.5), $"Complexity={complexity:F4}");
    }

    /// <summary>
    /// Complexity is always in [0, 1].
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Complexity_InRange()
    {
        return Prop.ForAll(
            Gen.Elements('A', 'C', 'G', 'T')
                .ArrayOf()
                .Where(a => a.Length >= 4)
                .Select(a => new DnaSequence(new string(a)))
                .ToArbitrary(),
            dna =>
            {
                double c = SequenceComplexity.CalculateLinguisticComplexity(dna);
                return (c >= 0.0 && c <= 1.0 + 0.0001)
                    .Label($"Complexity={c:F4} must be in [0,1]");
            });
    }

    // -- SEQ-ENTROPY-001 --

    /// <summary>
    /// Shannon entropy of a homopolymer is 0.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Homopolymer_HasZeroEntropy()
    {
        var dna = new DnaSequence("AAAAAAAAAA");
        double entropy = SequenceComplexity.CalculateShannonEntropy(dna);
        Assert.That(entropy, Is.EqualTo(0.0).Within(0.0001));
    }

    /// <summary>
    /// Shannon entropy of equifrequent bases is maximal (2.0 bits for 4 symbols).
    /// </summary>
    [Test]
    [Category("Property")]
    public void Equifrequent_HasMaximalEntropy()
    {
        var dna = new DnaSequence("ACGTACGTACGTACGTACGT");
        double entropy = SequenceComplexity.CalculateShannonEntropy(dna);
        Assert.That(entropy, Is.EqualTo(2.0).Within(0.01));
    }

    /// <summary>
    /// Shannon entropy is always non-negative.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Entropy_IsNonNegative()
    {
        return Prop.ForAll(
            Gen.Elements('A', 'C', 'G', 'T')
                .ArrayOf()
                .Where(a => a.Length >= 4)
                .Select(a => new DnaSequence(new string(a)))
                .ToArbitrary(),
            dna =>
            {
                double e = SequenceComplexity.CalculateShannonEntropy(dna);
                return (e >= -0.0001).Label($"Entropy={e:F4} must be ≥ 0");
            });
    }
}
