// MOTIF-CONS-001 — Consensus Sequence from a Multiple Alignment
// Evidence: docs/Evidence/MOTIF-CONS-001-Evidence.md
// TestSpec: tests/TestSpecs/MOTIF-CONS-001.md
// Source: Rosalind, Consensus and Profile (CONS), https://rosalind.info/problems/cons/ (sample dataset).
//         Wikipedia, Consensus sequence (most-frequent residue per column).
//         Los Alamos HIV Database, Advanced Consensus (alphabetical tie-break).

using System;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Canonical test class for MOTIF-CONS-001: column-wise most-frequent consensus from an
/// aligned set of DNA sequences. Verifies
/// <see cref="MotifFinder.CreateConsensusFromAlignment(System.Collections.Generic.IEnumerable{string})"/>
/// against the Rosalind CONS worked example and the alphabetical tie-break rule.
/// </summary>
[TestFixture]
public class MotifFinder_CreateConsensusFromAlignment_Tests
{
    #region CreateConsensusFromAlignment — MUST

    // M1 — Rosalind CONS sample: 7 strings of length 8. Published profile column maxima
    // (A,T,G,C,A,A,C,T) give consensus "ATGCAACT". Source: https://rosalind.info/problems/cons/.
    [Test]
    public void CreateConsensusFromAlignment_RosalindSample_ReturnsPublishedConsensus()
    {
        var aligned = new[]
        {
            "ATCCAGCT", "GGGCAACT", "ATGGATCT", "AAGCAACC",
            "TTGGAACT", "ATGCCATT", "ATGGCACT"
        };

        string consensus = MotifFinder.CreateConsensusFromAlignment(aligned);

        Assert.That(consensus, Is.EqualTo("ATGCAACT"),
            "Rosalind CONS sample: taking the maximum-count base per column yields ATGCAACT.");
    }

    // M2 — Identical sequences: every column is unanimous, so the consensus equals the input (INV-04).
    [Test]
    public void CreateConsensusFromAlignment_IdenticalSequences_ReturnsThatSequence()
    {
        var aligned = new[] { "ACGT", "ACGT", "ACGT" };

        string consensus = MotifFinder.CreateConsensusFromAlignment(aligned);

        Assert.That(consensus, Is.EqualTo("ACGT"),
            "Unanimous columns must reproduce the input sequence exactly (INV-04).");
    }

    // M3 — Alphabetical tie-break: column 1 has A and G tied (count 1 each); the
    // alphabetically-earlier base A is chosen (INV-03). Column 2 is unanimous T.
    [Test]
    public void CreateConsensusFromAlignment_TiedColumn_ChoosesAlphabeticallyEarliestBase()
    {
        var aligned = new[] { "AT", "GT" };

        string consensus = MotifFinder.CreateConsensusFromAlignment(aligned);

        Assert.That(consensus, Is.EqualTo("AT"),
            "A and G tie in column 1; alphabetical tie-break (A<C<G<T) selects A (INV-03).");
    }

    // M4 — Single sequence: each column's only base is its maximum, so the input is returned.
    [Test]
    public void CreateConsensusFromAlignment_SingleSequence_ReturnsItUnchanged()
    {
        var aligned = new[] { "GATTACA" };

        string consensus = MotifFinder.CreateConsensusFromAlignment(aligned);

        Assert.That(consensus, Is.EqualTo("GATTACA"),
            "With one sequence every column has a single base, which is the consensus.");
    }

    // M5 — Case-insensitive: lowercase input is normalised to uppercase before counting.
    [Test]
    public void CreateConsensusFromAlignment_LowercaseInput_NormalisedToUppercase()
    {
        var aligned = new[] { "acgt", "ACGT" };

        string consensus = MotifFinder.CreateConsensusFromAlignment(aligned);

        Assert.That(consensus, Is.EqualTo("ACGT"),
            "Lowercase and uppercase of the same base must count together; output is uppercase.");
    }

    #endregion

    #region CreateConsensusFromAlignment — SHOULD (edge cases)

    // S1 — Null collection throws ArgumentNullException.
    [Test]
    public void CreateConsensusFromAlignment_NullInput_ThrowsArgumentNullException()
    {
        Assert.That(() => MotifFinder.CreateConsensusFromAlignment(null!),
            NUnit.Framework.Throws.ArgumentNullException,
            "A null sequence collection is invalid input.");
    }

    // S2 — Empty collection returns empty string (nothing to summarise).
    [Test]
    public void CreateConsensusFromAlignment_EmptyCollection_ReturnsEmptyString()
    {
        string consensus = MotifFinder.CreateConsensusFromAlignment(Array.Empty<string>());

        Assert.That(consensus, Is.EqualTo(string.Empty),
            "An empty alignment has no columns, so the consensus is the empty string.");
    }

    // S3 — Unequal lengths violate the equal-length precondition (Rosalind CONS) → ArgumentException.
    [Test]
    public void CreateConsensusFromAlignment_UnequalLengths_ThrowsArgumentException()
    {
        var aligned = new[] { "AC", "ACG" };

        Assert.That(() => MotifFinder.CreateConsensusFromAlignment(aligned),
            NUnit.Framework.Throws.ArgumentException,
            "Consensus is defined only over equal-length aligned strings (Rosalind CONS).");
    }

    #endregion

    #region CreateConsensusFromAlignment — COULD

    // C1 — Non-ACGT character is rejected (alphabet validation, as in CreatePwm).
    [Test]
    public void CreateConsensusFromAlignment_InvalidCharacter_ThrowsArgumentException()
    {
        var aligned = new[] { "AX" };

        Assert.That(() => MotifFinder.CreateConsensusFromAlignment(aligned),
            NUnit.Framework.Throws.ArgumentException,
            "Only A, C, G, T are valid; 'X' must be rejected.");
    }

    // C2 — Pure majority (no tie): a column with A,A,C has A as the strict maximum.
    [Test]
    public void CreateConsensusFromAlignment_MajorityColumn_ChoosesMostFrequentBase()
    {
        var aligned = new[] { "AA", "AA", "CC" };
        // Column 1: A,A,C -> A (count 2 > 1). Column 2: A,A,C -> A (count 2 > 1).

        string consensus = MotifFinder.CreateConsensusFromAlignment(aligned);

        Assert.That(consensus, Is.EqualTo("AA"),
            "A appears twice vs C once in each column; the strict majority A is the consensus.");
    }

    #endregion
}
