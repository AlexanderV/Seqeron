using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for core sequence operations (complement, reverse complement).
/// Verifies algebraic invariants that must hold for ALL valid DNA sequences.
///
/// Test Unit: SEQ-COMP-001, SEQ-REVCOMP-001 (Property Extensions)
/// Evidence: DNA complementarity rules (Watson-Crick base pairing).
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Core")]
public class SequenceProperties
{
    private static Arbitrary<string> DnaArbitrary() =>
        Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length > 0)
            .Select(a => new string(a))
            .ToArbitrary();

    /// <summary>
    /// Complement is an involution: complement(complement(x)) == x.
    /// Evidence: Watson-Crick pairing — A↔T, G↔C applied twice returns original.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Complement_IsInvolution()
    {
        return Prop.ForAll(DnaArbitrary(), seq =>
        {
            var dna = new DnaSequence(seq);
            var doubleComp = dna.Complement().Complement();
            return (doubleComp.Sequence == dna.Sequence)
                .Label($"complement(complement(\"{seq}\")) should equal original");
        });
    }

    /// <summary>
    /// Reverse complement is an involution: revcomp(revcomp(x)) == x.
    /// Evidence: Reverse of reverse is identity; complement of complement is identity.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ReverseComplement_IsInvolution()
    {
        return Prop.ForAll(DnaArbitrary(), seq =>
        {
            var dna = new DnaSequence(seq);
            var doubleRc = dna.ReverseComplement().ReverseComplement();
            return (doubleRc.Sequence == dna.Sequence)
                .Label($"revcomp(revcomp(\"{seq}\")) should equal original");
        });
    }

    /// <summary>
    /// Complement preserves length: |complement(x)| == |x|.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Complement_PreservesLength()
    {
        return Prop.ForAll(DnaArbitrary(), seq =>
        {
            var dna = new DnaSequence(seq);
            return (dna.Complement().Length == dna.Length)
                .Label($"len={dna.Length}, complement len={dna.Complement().Length}");
        });
    }

    /// <summary>
    /// Reverse complement preserves length: |revcomp(x)| == |x|.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ReverseComplement_PreservesLength()
    {
        return Prop.ForAll(DnaArbitrary(), seq =>
        {
            var dna = new DnaSequence(seq);
            return (dna.ReverseComplement().Length == dna.Length)
                .Label($"len={dna.Length}, revcomp len={dna.ReverseComplement().Length}");
        });
    }

    /// <summary>
    /// Complement + Reverse == Reverse + Complement (operations commute).
    /// Evidence: Complement is position-independent, reverse is value-independent.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Complement_And_Reverse_Commute()
    {
        return Prop.ForAll(DnaArbitrary(), seq =>
        {
            var dna = new DnaSequence(seq);
            // Path 1: complement then reverse
            string compThenRev = new string(dna.Complement().Sequence.Reverse().ToArray());
            // Path 2: reverse then complement (= ReverseComplement)
            string revThenComp = dna.ReverseComplement().Sequence;

            return (compThenRev == revThenComp)
                .Label("complement→reverse should equal reverse→complement");
        });
    }

    /// <summary>
    /// DnaSequence construction round-trip: new DnaSequence(seq).Sequence == seq (uppercase).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Construction_PreservesSequenceContent()
    {
        return Prop.ForAll(DnaArbitrary(), seq =>
        {
            var dna = new DnaSequence(seq);
            return (dna.Sequence == seq.ToUpperInvariant())
                .Label($"Expected \"{seq.ToUpperInvariant()}\", got \"{dna.Sequence}\"");
        });
    }

    /// <summary>
    /// Subsequence length matches requested length.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Subsequence_HasCorrectLength()
    {
        return Prop.ForAll(DnaArbitrary(), seq =>
        {
            var dna = new DnaSequence(seq);
            if (dna.Length < 2) return true.ToProperty();

            int start = 0;
            int length = Math.Min(dna.Length / 2, dna.Length - start);
            var sub = dna.Subsequence(start, length);
            return (sub.Length == length)
                .Label($"Expected length {length}, got {sub.Length}");
        });
    }
}
