// VARIANT-SNP-001 — SNP Detection
// Evidence: docs/Evidence/VARIANT-SNP-001-Evidence.md
// TestSpec: tests/TestSpecs/VARIANT-SNP-001.md
// Source: VCFv4.3 spec (samtools/hts-specs); Hamming Distance (PMC5410656);
//         Wikipedia Transversion (Futuyma 2013); Collins & Jukes (1994) Genomics 20(3):386-396.

namespace Seqeron.Genomics.Tests.Unit.Annotation;

[TestFixture]
public class VariantCaller_FindSnps_Tests
{
    #region FindSnpsDirect

    // M1 — Source 3 (Hamming distance of equal strings is 0); INV-01.
    [Test]
    public void FindSnpsDirect_IdenticalSequences_ReturnsNoSnps()
    {
        var snps = VariantCaller.FindSnpsDirect("ATGC", "ATGC").ToList();

        Assert.That(snps, Is.Empty,
            "Identical equal-length sequences have Hamming distance 0, so no substitutions and zero SNPs (PMC5410656).");
    }

    // M2 — Source 1 (simple SNP = single-base substitution) + Source 3 (mismatch position); INV-04.
    [Test]
    public void FindSnpsDirect_SingleSubstitution_ReturnsExactSnp()
    {
        // "ATGC" vs "ATTC": the only differing index is 2 (G vs T).
        var snps = VariantCaller.FindSnpsDirect("ATGC", "ATTC").ToList();

        Assert.That(snps, Has.Count.EqualTo(1), "Exactly one position differs, so exactly one SNP.");
        Assert.Multiple(() =>
        {
            Assert.That(snps[0].Type, Is.EqualTo(VariantType.SNP), "A single-base substitution is a SNP (VCFv4.3 §1.1).");
            Assert.That(snps[0].Position, Is.EqualTo(2), "The substituted base is at 0-based position 2.");
            Assert.That(snps[0].ReferenceAllele, Is.EqualTo("G"), "Reference base at position 2 is G.");
            Assert.That(snps[0].AlternateAllele, Is.EqualTo("T"), "Query base at position 2 is T.");
            Assert.That(snps[0].QueryPosition, Is.EqualTo(2), "Positional comparison: query index equals reference index.");
        });
    }

    // M3 — Source 3 (each mismatch is one substitution); INV-04.
    [Test]
    public void FindSnpsDirect_MultipleSubstitutions_ReturnsSnpsAtExactPositions()
    {
        // "AAAA" vs "TGTA": indices 0,1,2 differ (A!=T, A!=G, A!=T); index 3 matches (A==A).
        var snps = VariantCaller.FindSnpsDirect("AAAA", "TGTA").ToList();

        Assert.That(snps, Has.Count.EqualTo(3), "Three differing positions (0,1,2) yield three SNPs; position 3 matches.");
        Assert.Multiple(() =>
        {
            Assert.That(snps.Select(s => s.Position), Is.EqualTo(new[] { 0, 1, 2 }), "Exact 0-based mismatch positions.");
            Assert.That(snps.Select(s => s.ReferenceAllele), Is.EqualTo(new[] { "A", "A", "A" }), "All reference bases are A.");
            Assert.That(snps.Select(s => s.AlternateAllele), Is.EqualTo(new[] { "T", "G", "T" }), "Query bases at the mismatches.");
        });
    }

    // M4 — Source 1 (a SNP is a substitution: type SNP, REF != ALT); INV-02, INV-03.
    [Test]
    public void FindSnpsDirect_AllResults_AreSnpsWithDistinctAlleles()
    {
        var snps = VariantCaller.FindSnpsDirect("AAAA", "TGTA").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(snps.All(s => s.Type == VariantType.SNP), "Every reported variant is a SNP (INV-02).");
            Assert.That(snps.All(s => s.ReferenceAllele != s.AlternateAllele),
                "Every SNP is a substitution, so reference and alternate alleles differ (INV-03).");
        });
    }

    // M5 — Source 3 (Hamming defined for equal-length only → common prefix); INV-06, ASSUMPTION-1.
    [Test]
    public void FindSnpsDirect_UnequalLengths_ComparesCommonPrefixOnly()
    {
        // ref "ATGCAA" (len 6) vs query "ATTC" (len 4): compares indices 0..3 only.
        // index 0 A==A, 1 T==T, 2 G!=T (SNP), 3 C==C; trailing "AA" of the reference is ignored.
        var snps = VariantCaller.FindSnpsDirect("ATGCAA", "ATTC").ToList();

        Assert.That(snps, Has.Count.EqualTo(1),
            "Only the common prefix (length 4) is compared; one mismatch within it (Hamming defined for equal length, PMC5410656).");
        Assert.Multiple(() =>
        {
            Assert.That(snps[0].Position, Is.EqualTo(2), "The single in-prefix mismatch is at position 2.");
            Assert.That(snps[0].ReferenceAllele, Is.EqualTo("G"), "Reference base at position 2 is G.");
            Assert.That(snps[0].AlternateAllele, Is.EqualTo("T"), "Query base at position 2 is T.");
        });
    }

    // S1 — documented input contract (empty input).
    [Test]
    public void FindSnpsDirect_EmptyInput_ReturnsEmpty()
    {
        var snps = VariantCaller.FindSnpsDirect("", "").ToList();

        Assert.That(snps, Is.Empty, "Empty inputs have no positions to compare and therefore no SNPs.");
    }

    // S2 — documented input contract (one empty operand).
    [Test]
    public void FindSnpsDirect_OneEmptyOperand_ReturnsEmpty()
    {
        var snps = VariantCaller.FindSnpsDirect("ATGC", "").ToList();

        Assert.That(snps, Is.Empty, "With one empty operand the common prefix is empty, so no SNPs.");
    }

    // C1 — property (Source 3): SNP count equals Hamming distance over an equal-length pair; INV-05.
    [Test]
    public void FindSnpsDirect_AnyEqualLengthPair_CountEqualsHammingDistance()
    {
        // Deterministic equal-length pair; differences at indices 0,2,5,7.
        const string reference = "ACGTACGT";
        const string query     = "TCATACGA";
        int hammingDistance = reference.Where((c, i) => c != query[i]).Count();

        var snps = VariantCaller.FindSnpsDirect(reference, query).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(snps, Has.Count.EqualTo(hammingDistance),
                "The number of SNPs from FindSnpsDirect equals the Hamming distance of the two equal-length sequences (PMC5410656, INV-05).");
            Assert.That(snps.All(s => s.Type == VariantType.SNP), "Every reported variant is a SNP.");
            Assert.That(snps.All(s => s.ReferenceAllele != s.AlternateAllele), "Every SNP has distinct alleles (INV-03).");
        });
    }

    #endregion

    #region FindSnps (alignment-based delegate)

    // M6 — Source 1 (SNP-only filter excludes indels); INV-02.
    [Test]
    public void FindSnps_SubstitutionOnlyInput_ReturnsSnpsOnly()
    {
        // Equal-length single-substitution input: position 3 C->A; no indels expected.
        var reference = new DnaSequence("ATGCATGC");
        var query = new DnaSequence("ATGAATGC");

        var snps = VariantCaller.FindSnps(reference, query).ToList();

        Assert.That(snps, Has.Count.EqualTo(1), "A single substitution yields exactly one SNP and no indels.");
        Assert.Multiple(() =>
        {
            Assert.That(snps.All(s => s.Type == VariantType.SNP), "FindSnps returns SNPs only (insertions/deletions filtered out).");
            Assert.That(snps[0].Position, Is.EqualTo(3), "The substituted base is at 0-based reference position 3.");
            Assert.That(snps[0].ReferenceAllele, Is.EqualTo("C"), "Reference base at position 3 is C.");
            Assert.That(snps[0].AlternateAllele, Is.EqualTo("A"), "Query base at position 3 is A.");
        });
    }

    // S5 — delegation: identical sequences have no differences, so no SNPs; INV-01.
    [Test]
    public void FindSnps_IdenticalSequences_ReturnsNoSnps()
    {
        var reference = new DnaSequence("ATGCATGC");
        var query = new DnaSequence("ATGCATGC");

        var snps = VariantCaller.FindSnps(reference, query).ToList();

        Assert.That(snps, Is.Empty, "Identical sequences have no substitutions, so FindSnps returns no SNPs.");
    }

    // S3 — input validation propagated from CallVariants.
    [Test]
    public void FindSnps_NullReference_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(
            () => VariantCaller.FindSnps(null!, new DnaSequence("ATGC")).ToList(),
            "A null reference is invalid input and must throw ArgumentNullException.");
    }

    // S4 — input validation propagated from CallVariants.
    [Test]
    public void FindSnps_NullQuery_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(
            () => VariantCaller.FindSnps(new DnaSequence("ATGC"), null!).ToList(),
            "A null query is invalid input and must throw ArgumentNullException.");
    }

    #endregion
}
