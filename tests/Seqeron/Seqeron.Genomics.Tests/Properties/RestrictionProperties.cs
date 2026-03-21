using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for restriction enzyme analysis: site finding and digestion.
///
/// Test Units: RESTR-FIND-001, RESTR-DIGEST-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("MolTools")]
public class RestrictionProperties
{
    // Sequence containing known EcoRI (GAATTC) and BamHI (GGATCC) sites
    private const string TestSequence =
        "AAAGAATTCAAAGGATCCAAAGAATTCAAAGGATCCAAA";

    private static Arbitrary<string> DnaArbitrary(int minLen = 20) =>
        Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length >= minLen)
            .Select(a => new string(a))
            .ToArbitrary();

    #region RESTR-FIND-001: R: positions valid; P: site sequence matches enzyme recognition; D: deterministic

    /// <summary>
    /// INV-1: All site positions are within sequence bounds.
    /// Evidence: A recognition site at position p with length L requires 0 ≤ p and p + L ≤ seqLen.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindSites_Positions_WithinBounds()
    {
        var dna = new DnaSequence(TestSequence);
        var sites = RestrictionAnalyzer.FindSites(dna, "EcoRI").ToList();

        foreach (var site in sites)
        {
            Assert.That(site.Position, Is.GreaterThanOrEqualTo(0));
            Assert.That(site.Position, Is.LessThan(TestSequence.Length));
        }
    }

    /// <summary>
    /// INV-2: Found restriction sites contain the recognized sequence matching the enzyme.
    /// Evidence: The enzyme's recognition site must be present at the reported position.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindSites_RecognizedSequence_MatchesEnzyme()
    {
        var dna = new DnaSequence(TestSequence);
        var enzyme = RestrictionAnalyzer.GetEnzyme("EcoRI");
        Assert.That(enzyme, Is.Not.Null, "EcoRI not found in enzyme database");

        var sites = RestrictionAnalyzer.FindSites(dna, "EcoRI").ToList();

        foreach (var site in sites)
            Assert.That(site.RecognizedSequence, Does.Contain(enzyme!.RecognitionSequence).IgnoreCase,
                $"Site at {site.Position} recognized '{site.RecognizedSequence}' doesn't match '{enzyme.RecognitionSequence}'");
    }

    /// <summary>
    /// INV-3: The recognition sequence at the reported position matches the enzyme pattern.
    /// Evidence: Extracting the substring at the position must yield the recognition site.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindSites_SubstringAtPosition_MatchesRecognition()
    {
        var dna = new DnaSequence(TestSequence);
        var enzyme = RestrictionAnalyzer.GetEnzyme("EcoRI");
        Assert.That(enzyme, Is.Not.Null);

        var sites = RestrictionAnalyzer.FindSites(dna, "EcoRI").ToList();

        foreach (var site in sites)
        {
            string actual = TestSequence.Substring(site.Position, enzyme!.RecognitionLength);
            Assert.That(actual, Is.EqualTo(enzyme.RecognitionSequence).IgnoreCase,
                $"Substring at position {site.Position} = '{actual}', expected '{enzyme.RecognitionSequence}'");
        }
    }

    /// <summary>
    /// INV-4: FindSites is deterministic — same input always yields same results.
    /// Evidence: FindSites is a pure function with no side effects.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindSites_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(30), seq =>
        {
            var dna = new DnaSequence(seq);
            var sites1 = RestrictionAnalyzer.FindSites(dna, "EcoRI").ToList();
            var sites2 = RestrictionAnalyzer.FindSites(dna, "EcoRI").ToList();
            bool equal = sites1.Count == sites2.Count &&
                         sites1.Zip(sites2).All(pair =>
                             pair.First.Position == pair.Second.Position);
            return equal.Label("FindSites must be deterministic");
        });
    }

    /// <summary>
    /// INV-5: FindSites positions are valid for random DNA sequences.
    /// Evidence: All positions must be within [0, seqLen).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindSites_Positions_ValidForRandomDna()
    {
        return Prop.ForAll(DnaArbitrary(30), seq =>
        {
            var dna = new DnaSequence(seq);
            var sites = RestrictionAnalyzer.FindSites(dna, "EcoRI").ToList();
            return sites.All(s => s.Position >= 0 && s.Position < seq.Length)
                .Label("All positions must be in [0, seqLen)");
        });
    }

    /// <summary>
    /// FindSites for non-existing recognition site returns empty.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindSites_NoSite_ReturnsEmpty()
    {
        var dna = new DnaSequence("ACGTACGTACGTACGT");
        var sites = RestrictionAnalyzer.FindSites(dna, "HindIII").ToList();
        Assert.That(sites, Is.Empty);
    }

    /// <summary>
    /// Enzymes dictionary is not empty and all entries have recognition sequences.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Enzymes_AllHaveRecognitionSequence()
    {
        var enzymes = RestrictionAnalyzer.Enzymes;
        Assert.That(enzymes.Count, Is.GreaterThan(0), "Enzyme database should not be empty");

        foreach (var kvp in enzymes)
        {
            Assert.That(kvp.Value.RecognitionSequence, Is.Not.Null.And.Not.Empty,
                $"Enzyme {kvp.Key} has no recognition sequence");
            Assert.That(kvp.Value.RecognitionLength, Is.GreaterThan(0),
                $"Enzyme {kvp.Key} has zero recognition length");
        }
    }

    #endregion

    #region RESTR-DIGEST-001: P: sum(fragment lengths) = seq length; R: fragments ≥ 1; D: deterministic

    /// <summary>
    /// INV-1: Digest fragment lengths sum to original sequence length.
    /// Evidence: Digestion partitions the sequence without loss or overlap.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Digest_FragmentLengths_SumToSequenceLength()
    {
        var dna = new DnaSequence(TestSequence);
        var fragments = RestrictionAnalyzer.Digest(dna, "EcoRI").ToList();

        int totalLength = fragments.Sum(f => f.Length);
        Assert.That(totalLength, Is.EqualTo(TestSequence.Length),
            $"Fragment lengths sum {totalLength} ≠ sequence length {TestSequence.Length}");
    }

    /// <summary>
    /// INV-2: Fragment lengths sum to sequence length for random DNA with embedded sites.
    /// Evidence: Universal partition invariant — no bases are lost or duplicated.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Digest_FragmentLengthSum_EqualsSequenceLength_Property()
    {
        return Prop.ForAll(DnaArbitrary(30), seq =>
        {
            var dna = new DnaSequence(seq);
            var fragments = RestrictionAnalyzer.Digest(dna, "EcoRI").ToList();
            if (fragments.Count == 0) return true.ToProperty();
            int totalLen = fragments.Sum(f => f.Length);
            return (totalLen == seq.Length)
                .Label($"Fragment sum={totalLen}, seqLen={seq.Length}");
        });
    }

    /// <summary>
    /// INV-3: Number of digest fragments is ≥ 1 (at least the whole sequence).
    /// Evidence: Even without cut sites, the entire sequence is returned as one fragment.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Digest_FragmentCount_AtLeastOne()
    {
        return Prop.ForAll(DnaArbitrary(30), seq =>
        {
            var dna = new DnaSequence(seq);
            var fragments = RestrictionAnalyzer.Digest(dna, "EcoRI").ToList();
            return (fragments.Count >= 1)
                .Label($"Expected ≥1 fragment, got {fragments.Count}");
        });
    }

    /// <summary>
    /// INV-4: When cut sites exist, number of fragments ≥ 2.
    /// Evidence: At least one cut produces at least two pieces.
    /// Note: fragment count may differ from total FindSites count because
    /// palindromic enzymes report both strand orientations.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Digest_WithSites_AtLeastTwoFragments()
    {
        var dna = new DnaSequence(TestSequence);
        var sites = RestrictionAnalyzer.FindSites(dna, "EcoRI").ToList();
        var fragments = RestrictionAnalyzer.Digest(dna, "EcoRI").ToList();

        Assert.That(sites.Count, Is.GreaterThan(0), "Should find at least one EcoRI site");
        Assert.That(fragments.Count, Is.GreaterThanOrEqualTo(2),
            $"Expected ≥ 2 fragments when sites exist, got {fragments.Count}");
    }

    /// <summary>
    /// INV-5: Digest is deterministic — same input always yields same fragments.
    /// Evidence: Digest is a pure function with no side effects.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Digest_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(30), seq =>
        {
            var dna = new DnaSequence(seq);
            var f1 = RestrictionAnalyzer.Digest(dna, "EcoRI").ToList();
            var f2 = RestrictionAnalyzer.Digest(dna, "EcoRI").ToList();
            bool equal = f1.Count == f2.Count &&
                         f1.Zip(f2).All(p => p.First.Length == p.Second.Length &&
                                              p.First.StartPosition == p.Second.StartPosition);
            return equal.Label("Digest must be deterministic");
        });
    }

    /// <summary>
    /// INV-6: All fragment lengths are positive.
    /// Evidence: A zero-length fragment indicates a degenerate digestion.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Digest_AllFragments_HavePositiveLength()
    {
        var dna = new DnaSequence(TestSequence);
        var fragments = RestrictionAnalyzer.Digest(dna, "EcoRI", "BamHI").ToList();

        foreach (var f in fragments)
            Assert.That(f.Length, Is.GreaterThan(0), $"Fragment {f.FragmentNumber} has non-positive length {f.Length}");
    }

    /// <summary>
    /// INV-7: Multi-enzyme digest produces more or equal fragments than single enzyme.
    /// Evidence: Additional cut sites can only increase the number of fragments.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Digest_MultiEnzyme_MoreOrEqualFragments()
    {
        var dna = new DnaSequence(TestSequence);
        var singleFragments = RestrictionAnalyzer.Digest(dna, "EcoRI").ToList();
        var multiFragments = RestrictionAnalyzer.Digest(dna, "EcoRI", "BamHI").ToList();

        Assert.That(multiFragments.Count, Is.GreaterThanOrEqualTo(singleFragments.Count),
            "Multi-enzyme digest should produce ≥ fragments than single enzyme");
    }

    /// <summary>
    /// INV-8: DigestSummary total matches fragment count.
    /// Evidence: Summary is derived from the same digest operation.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DigestSummary_TotalFragments_MatchesDigest()
    {
        var dna = new DnaSequence(TestSequence);
        var summary = RestrictionAnalyzer.GetDigestSummary(dna, "EcoRI");
        var fragments = RestrictionAnalyzer.Digest(dna, "EcoRI").ToList();

        Assert.That(summary.TotalFragments, Is.EqualTo(fragments.Count));
    }

    #endregion
}
