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

    private static Arbitrary<DnaSequence> DnaSequenceArbitrary(int minLen = 4) =>
        Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length >= minLen)
            .Select(a => new DnaSequence(new string(a)))
            .ToArbitrary();

    #region SEQ-VALID-001: R: result ∈ {true,false}; P: case-insensitive; D: deterministic; P: valid DNA ⊂ valid IUPAC

    /// <summary>
    /// INV-1: Pure ACGT sequences are always valid DNA.
    /// Evidence: IsValidDna accepts only {A,C,G,T} characters.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PureAcgt_IsValidDna()
    {
        return Prop.ForAll(
            DnaStringArbitrary(),
            seq => seq.AsSpan().IsValidDna().Label("Pure ACGT must be valid DNA"));
    }

    /// <summary>
    /// INV-2: Pure ACGU sequences are always valid RNA.
    /// Evidence: IsValidRna accepts only {A,C,G,U} characters.
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
    /// INV-3: DNA validation is case-insensitive.
    /// Evidence: IsValidDna uses char.ToUpperInvariant internally.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Validation_IsCaseInsensitive()
    {
        return Prop.ForAll(
            DnaStringArbitrary(),
            seq =>
            {
                bool upper = seq.ToUpperInvariant().AsSpan().IsValidDna();
                bool lower = seq.ToLowerInvariant().AsSpan().IsValidDna();
                bool mixed = seq.AsSpan().IsValidDna();
                return (upper == lower && lower == mixed)
                    .Label($"Case sensitivity mismatch: upper={upper}, lower={lower}, mixed={mixed}");
            });
    }

    /// <summary>
    /// INV-4: DNA validation is deterministic — same input always yields same result.
    /// Evidence: IsValidDna is a pure function with no side effects.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Validation_IsDeterministic()
    {
        return Prop.ForAll(
            DnaStringArbitrary(),
            seq =>
            {
                bool result1 = seq.AsSpan().IsValidDna();
                bool result2 = seq.AsSpan().IsValidDna();
                return (result1 == result2)
                    .Label("IsValidDna must be deterministic");
            });
    }

    /// <summary>
    /// INV-5: Valid DNA ⊂ valid IUPAC — every DNA base (A,C,G,T) is a valid IUPAC code.
    /// Evidence: IUPAC-IUB standard includes A,C,G,T as the four standard nucleotides.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ValidDna_IsSubsetOfIupac()
    {
        return Prop.ForAll(
            DnaStringArbitrary(),
            seq =>
            {
                bool allIupac = seq.All(c => IupacHelper.MatchesIupac(char.ToUpperInvariant(c), char.ToUpperInvariant(c)));
                return allIupac.Label("Every valid DNA base must be a valid IUPAC code");
            });
    }

    /// <summary>
    /// INV-6: Non-DNA characters fail validation.
    /// Evidence: Characters outside {A,C,G,T} are rejected.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property InvalidCharacters_FailDnaValidation()
    {
        return Prop.ForAll(
            Gen.Elements('X', 'Z', '1', '@', ' ', 'N')
                .Select(c => c.ToString())
                .ToArbitrary(),
            seq => (!seq.AsSpan().IsValidDna())
                .Label($"'{seq}' must not be valid DNA"));
    }

    #endregion

    #region SEQ-COMPLEX-001: R: complexity ∈ [0,1]; M: homopolymer → min; M: random long → high; P: permutation invariance

    /// <summary>
    /// INV-1: Linguistic complexity is always in [0, 1].
    /// Evidence: LC = observed_subwords / possible_subwords, both positive with observed ≤ possible.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Complexity_InRange()
    {
        return Prop.ForAll(
            DnaSequenceArbitrary(),
            dna =>
            {
                double c = SequenceComplexity.CalculateLinguisticComplexity(dna);
                return (c >= 0.0 && c <= 1.0 + 0.0001)
                    .Label($"Complexity={c:F4} must be in [0,1]");
            });
    }

    /// <summary>
    /// INV-2: Homopolymer sequences have minimal complexity.
    /// Evidence: Homopolymers produce only 1 distinct subword per length → LC close to 0.
    /// </summary>
    [TestCase("AAAAAAAAAA")]
    [TestCase("CCCCCCCCCC")]
    [TestCase("GGGGGGGGGG")]
    [TestCase("TTTTTTTTTT")]
    [Category("Property")]
    public void Homopolymer_HasLowComplexity(string seq)
    {
        var dna = new DnaSequence(seq);
        double complexity = SequenceComplexity.CalculateLinguisticComplexity(dna);
        Assert.That(complexity, Is.LessThan(0.25), $"Complexity={complexity:F4} for homopolymer");
    }

    /// <summary>
    /// INV-3: Complex sequences with high diversity have high complexity.
    /// Evidence: Many distinct subwords → LC approaches 1.
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
    /// INV-4: Complexity is invariant under bijective base substitution.
    /// Evidence: Replacing A↔C and G↔T preserves the number of distinct subwords.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Complexity_BaseSubstitutionInvariant()
    {
        return Prop.ForAll(
            DnaSequenceArbitrary(6),
            dna =>
            {
                double original = SequenceComplexity.CalculateLinguisticComplexity(dna);
                string substituted = new(dna.Sequence.Select(c => c switch
                {
                    'A' => 'C',
                    'C' => 'A',
                    'G' => 'T',
                    'T' => 'G',
                    _ => c
                }).ToArray());
                double after = SequenceComplexity.CalculateLinguisticComplexity(new DnaSequence(substituted));
                return (Math.Abs(original - after) < 0.0001)
                    .Label($"Original={original:F4}, Substituted={after:F4}");
            });
    }

    #endregion

    #region SEQ-ENTROPY-001: R: entropy ≥ 0; P: permutation invariance; M: uniform dist → max entropy; D: deterministic

    /// <summary>
    /// INV-1: Shannon entropy is always non-negative.
    /// Evidence: H = -Σ p·log₂(p) where p ∈ (0,1] → each term ≥ 0.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Entropy_IsNonNegative()
    {
        return Prop.ForAll(
            DnaSequenceArbitrary(),
            dna =>
            {
                double e = SequenceComplexity.CalculateShannonEntropy(dna);
                return (e >= -0.0001).Label($"Entropy={e:F4} must be ≥ 0");
            });
    }

    /// <summary>
    /// INV-2: Homopolymer has zero entropy.
    /// Evidence: Single symbol → p=1 → H = -1·log₂(1) = 0.
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
    /// INV-3: Equifrequent bases yield maximum entropy (2.0 bits for 4-symbol alphabet).
    /// Evidence: H_max = log₂(4) = 2.0 when all p_i = 1/4.
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
    /// INV-4: Shannon entropy is invariant under sequence permutation (shuffle).
    /// Evidence: Entropy depends only on character frequencies, not ordering.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Entropy_PermutationInvariant()
    {
        return Prop.ForAll(
            DnaSequenceArbitrary(8),
            dna =>
            {
                double original = SequenceComplexity.CalculateShannonEntropy(dna);
                // Reverse the sequence — same character frequencies, different order
                string reversed = new(dna.Sequence.Reverse().ToArray());
                double reversed_e = SequenceComplexity.CalculateShannonEntropy(new DnaSequence(reversed));
                return (Math.Abs(original - reversed_e) < 0.0001)
                    .Label($"Original={original:F4}, Reversed={reversed_e:F4}");
            });
    }

    /// <summary>
    /// INV-5: Shannon entropy is deterministic — same input always yields same result.
    /// Evidence: CalculateShannonEntropy is a pure function.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Entropy_IsDeterministic()
    {
        return Prop.ForAll(
            DnaSequenceArbitrary(),
            dna =>
            {
                double e1 = SequenceComplexity.CalculateShannonEntropy(dna);
                double e2 = SequenceComplexity.CalculateShannonEntropy(dna);
                return (Math.Abs(e1 - e2) < 0.0001)
                    .Label($"Entropy must be deterministic: {e1:F4} vs {e2:F4}");
            });
    }

    #endregion
}
