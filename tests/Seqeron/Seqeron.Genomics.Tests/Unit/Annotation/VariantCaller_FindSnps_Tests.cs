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

    // M5 — Source 3 (Hamming mismatch set defined only for equal-length strings); INV-06.
    // Reference cross-check: scipy 1.17.1 spatial.distance.hamming(list("ATGCAA"), list("ATTC")) and
    // scikit-bio 0.7.4 DNA.mismatches both raise ValueError ("must have equal lengths"). The previous
    // contract silently compared the common prefix only and hid the unmatched tail.
    [Test]
    public void FindSnpsDirect_UnequalLengths_ThrowsArgumentException()
    {
        Assert.Throws<ArgumentException>(
            () => VariantCaller.FindSnpsDirect("ATGCAA", "ATTC").ToList(),
            "Positional (Hamming) comparison is undefined for unequal lengths; the tail must not be silently ignored.");
        Assert.Throws<ArgumentException>(
            () => VariantCaller.FindSnpsDirect("ATGC", "ATGCAAAA").ToList(),
            "A longer query is equally undefined for positional SNP detection.");
    }

    // M7 — Source 1 (VCFv4.3 l.339/350: REF/ALT bases "A,C,G,T,N (case insensitive)").
    // Reference cross-check: bcftools norm -c e (pysam 0.24.1) rejects REF=C ALT=c as "Duplicate alleles".
    [Test]
    public void FindSnpsDirect_CaseOnlyDifference_IsNotASnp()
    {
        Assert.That(VariantCaller.FindSnpsDirect("acgt", "ACGT").ToList(), Is.Empty,
            "Soft-masked lowercase bases equal to the uppercase bases are matches, not SNPs.");
    }

    // M8 — Source 1 (case-insensitive comparison; alleles reported as given).
    [Test]
    public void FindSnpsDirect_LowercaseSubstitution_ReturnsSnpWithGivenAlleles()
    {
        var snps = VariantCaller.FindSnpsDirect("acGt", "ACtT").ToList();

        Assert.That(snps, Has.Count.EqualTo(1), "Only index 2 (G vs t) is a real substitution.");
        Assert.Multiple(() =>
        {
            Assert.That(snps[0].Position, Is.EqualTo(2));
            Assert.That(snps[0].ReferenceAllele, Is.EqualTo("G"));
            Assert.That(snps[0].AlternateAllele, Is.EqualTo("t"));
            Assert.That(VariantCaller.ClassifyMutation(snps[0]), Is.EqualTo(MutationType.Transversion),
                "G->T is a purine->pyrimidine transversion regardless of case.");
        });
    }

    // M9 — Source 1 (VCFv4.3 l.350: ALT is a base string over A,C,G,T,N or '*'; a gap is not a base,
    // a gap column is an insertion/deletion, not a substitution).
    [Test]
    public void FindSnpsDirect_GapColumns_AreNotSnps_PositionsAreUngapped()
    {
        // Aligned: ref "AC-TA" / query "ACGTC": column 2 is an insertion (not a SNP); column 4 A->C is a SNP
        // at ungapped reference index 3 and ungapped query index 4.
        var snps = VariantCaller.FindSnpsDirect("AC-TA", "ACGTC").ToList();

        Assert.That(snps, Has.Count.EqualTo(1), "The gap column is an indel and must not be reported as a SNP.");
        Assert.Multiple(() =>
        {
            Assert.That(snps[0].Position, Is.EqualTo(3), "Reference coordinate skips the gap column.");
            Assert.That(snps[0].QueryPosition, Is.EqualTo(4));
            Assert.That(snps[0].ReferenceAllele, Is.EqualTo("A"));
            Assert.That(snps[0].AlternateAllele, Is.EqualTo("C"));
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

    #region CalculateStatistics (uncovered public method, reviewed with VARIANT-SNP-001)

    // Reference cross-check: bcftools stats (pysam 0.24.1) on the 3 SNP records of this pair
    // (c:1 A>G, c:6 C>T, c:11 G>T) reports "number of SNPs: 3", "number of indels: 0", TSTV ts=2 tv=1 ts/tv=2.00.
    [Test]
    public void CalculateStatistics_SubstitutionOnlyPair_MatchesBcftoolsStats()
    {
        var reference = new DnaSequence("ACGTACGTACGTACGTACGT");
        var query = new DnaSequence("GCGTATGTACTTACGTACGT");

        var stats = VariantCaller.CalculateStatistics(reference, query);

        Assert.Multiple(() =>
        {
            Assert.That(stats.Snps, Is.EqualTo(3));
            Assert.That(stats.Insertions + stats.Deletions, Is.EqualTo(0));
            Assert.That(stats.TotalVariants, Is.EqualTo(3));
            Assert.That(stats.TiTvRatio, Is.EqualTo(2.0).Within(1e-12), "bcftools stats ts/tv = 2.00");
            Assert.That(stats.VariantDensity, Is.EqualTo(150.0).Within(1e-12), "3 variants / 20 bp * 1000");
            Assert.That(VariantCaller.CalculateTiTvRatio(
                    VariantCaller.FindSnpsDirect(reference.Sequence, query.Sequence)),
                Is.EqualTo(2.0).Within(1e-12), "Direct and alignment-based SNP sets agree on Ti/Tv.");
        });
    }


    // Transition/transversion is defined only between purines A,G and pyrimidines C,T (Evidence: Transversion
    // definition). Reference cross-check: bcftools stats (pysam 0.24.1) on records A>N and C>T reports
    // "number of SNPs: 2" but TSTV ts=1 tv=0 ts/tv=0.00 — the N change is neither Ti nor Tv.
    [Test]
    public void ClassifyMutation_AmbiguousBase_IsOther_AndExcludedFromTiTv()
    {
        var aToN = new Variant(0, "A", "N", VariantType.SNP, 0);
        var cToT = new Variant(1, "C", "T", VariantType.SNP, 1);

        Assert.Multiple(() =>
        {
            Assert.That(VariantCaller.ClassifyMutation(aToN), Is.EqualTo(MutationType.Other));
            Assert.That(VariantCaller.ClassifyMutation(new Variant(0, "n", "g", VariantType.SNP, 0)), Is.EqualTo(MutationType.Other));
            Assert.That(VariantCaller.CalculateTiTvRatio(new[] { aToN, cToT }), Is.EqualTo(0.0),
                "bcftools stats: ts=1, tv=0 -> ts/tv 0.00");
        });
    }

    #endregion
}
