// SEQ-STATS-001 — Sequence Composition Statistics
// Evidence: docs/Evidence/SEQ-STATS-001-Evidence.md
// TestSpec: tests/TestSpecs/SEQ-STATS-001.md
// Source: Lobry J.R. (1996). Mol Biol Evol 13(5):660-665; Biopython Bio.SeqUtils (gc_fraction, GC_skew); Wikipedia "GC skew".

using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class SequenceStatistics_CalculateNucleotideComposition_Tests
{
    private const double Tolerance = 1e-10;

    #region CalculateNucleotideComposition

    // M1 — exact per-base counts and length partition
    // Evidence: nucleotide composition counts each base; AAUUGGCC = A2 T0 G2 C2 U2
    [Test]
    public void CalculateNucleotideComposition_RnaSequence_ReturnsExactCounts()
    {
        var comp = SequenceStatistics.CalculateNucleotideComposition("AAUUGGCC");

        Assert.Multiple(() =>
        {
            Assert.That(comp.Length, Is.EqualTo(8), "Length equals input character count");
            Assert.That(comp.CountA, Is.EqualTo(2), "Two A");
            Assert.That(comp.CountT, Is.EqualTo(0), "No T (RNA)");
            Assert.That(comp.CountU, Is.EqualTo(2), "Two U");
            Assert.That(comp.CountG, Is.EqualTo(2), "Two G");
            Assert.That(comp.CountC, Is.EqualTo(2), "Two C");
            Assert.That(comp.CountN, Is.EqualTo(0), "No N");
            Assert.That(comp.CountOther, Is.EqualTo(0), "No other symbols");
        });
    }

    // M2 — GC content = (G+C)/total
    // Evidence: Biopython gc_fraction; ATGC -> (1+1)/4 = 0.5
    [Test]
    public void CalculateNucleotideComposition_BalancedDna_ReturnsGcContentHalf()
    {
        var comp = SequenceStatistics.CalculateNucleotideComposition("ATGC");

        Assert.That(comp.GcContent, Is.EqualTo(0.5).Within(Tolerance),
            "GC content (G+C)/total = 2/4 = 0.5 per Biopython gc_fraction");
    }

    // M3 — all-GC sequence -> GC content 1.0
    // Evidence: Biopython gc_fraction; GGGC -> 4/4 = 1.0
    [Test]
    public void CalculateNucleotideComposition_AllGc_ReturnsGcContentOne()
    {
        var comp = SequenceStatistics.CalculateNucleotideComposition("GGGC");

        Assert.That(comp.GcContent, Is.EqualTo(1.0).Within(Tolerance),
            "All bases are G/C so GC content = 1.0");
    }

    // M4 — AT content = (A+T+U)/total
    // Evidence: ATGC -> (1+1)/4 = 0.5
    [Test]
    public void CalculateNucleotideComposition_BalancedDna_ReturnsAtContentHalf()
    {
        var comp = SequenceStatistics.CalculateNucleotideComposition("ATGC");

        Assert.That(comp.AtContent, Is.EqualTo(0.5).Within(Tolerance),
            "AT content (A+T)/total = 2/4 = 0.5");
    }

    // M4b — AT content includes U for RNA: (A+T+U)/total
    // Evidence: SEQ-STATS-001 Evidence worked example AAUUGGCC -> A2 T0 U2 G2 C2,
    //           AT content (A+T+U)/total = (2+0+2)/8 = 0.5. Locks the documented U-inclusion branch.
    [Test]
    public void CalculateNucleotideComposition_RnaSequence_AtContentIncludesUracil()
    {
        var comp = SequenceStatistics.CalculateNucleotideComposition("AAUUGGCC");

        Assert.That(comp.AtContent, Is.EqualTo(0.5).Within(Tolerance),
            "AT content (A+T+U)/total = (2+0+2)/8 = 0.5 (U counted with A/T for RNA)");
    }

    // M7b — AT skew uses DNA formula (A-T)/(A+T) without U; AAUUGGCC -> (2-0)/(2+0) = 1.0
    // Evidence: SEQ-STATS-001 Evidence worked example AAUUGGCC AT skew = 1.0 (Wikipedia "GC skew" formula).
    [Test]
    public void CalculateNucleotideComposition_RnaSequence_AtSkewExcludesUracil()
    {
        var comp = SequenceStatistics.CalculateNucleotideComposition("AAUUGGCC");

        Assert.That(comp.AtSkew, Is.EqualTo(1.0).Within(Tolerance),
            "AT skew (A-T)/(A+T) = (2-0)/(2+0) = 1.0; U excluded per Lobry/Wikipedia DNA-specific formula");
    }

    // M5 — GC skew positive = (G-C)/(G+C)
    // Evidence: Wikipedia/Biopython GC_skew; GGGC -> (3-1)/4 = 0.5
    [Test]
    public void CalculateNucleotideComposition_GRich_ReturnsPositiveGcSkew()
    {
        var comp = SequenceStatistics.CalculateNucleotideComposition("GGGC");

        Assert.That(comp.GcSkew, Is.EqualTo(0.5).Within(Tolerance),
            "GC skew (G-C)/(G+C) = (3-1)/4 = 0.5 per Lobry/Biopython");
    }

    // M6 — GC skew negative
    // Evidence: GCCC -> (1-3)/4 = -0.5
    [Test]
    public void CalculateNucleotideComposition_CRich_ReturnsNegativeGcSkew()
    {
        var comp = SequenceStatistics.CalculateNucleotideComposition("GCCC");

        Assert.That(comp.GcSkew, Is.EqualTo(-0.5).Within(Tolerance),
            "GC skew (G-C)/(G+C) = (1-3)/4 = -0.5; negative = C-rich");
    }

    // M7 — AT skew = (A-T)/(A+T)
    // Evidence: Wikipedia "GC skew"; AAAT -> (3-1)/4 = 0.5
    [Test]
    public void CalculateNucleotideComposition_ARich_ReturnsPositiveAtSkew()
    {
        var comp = SequenceStatistics.CalculateNucleotideComposition("AAAT");

        Assert.That(comp.AtSkew, Is.EqualTo(0.5).Within(Tolerance),
            "AT skew (A-T)/(A+T) = (3-1)/4 = 0.5 per Wikipedia 'GC skew'");
    }

    // M8 — empty sequence -> all-zero composition
    // Evidence: Biopython gc_fraction returns 0 for empty sequence
    [Test]
    public void CalculateNucleotideComposition_EmptyString_ReturnsAllZero()
    {
        var comp = SequenceStatistics.CalculateNucleotideComposition("");

        Assert.Multiple(() =>
        {
            Assert.That(comp.Length, Is.EqualTo(0), "Empty length");
            Assert.That(comp.CountA + comp.CountT + comp.CountG + comp.CountC + comp.CountU, Is.EqualTo(0), "No counts");
            Assert.That(comp.GcContent, Is.EqualTo(0.0).Within(Tolerance), "GC content 0 for empty per Biopython");
            Assert.That(comp.GcSkew, Is.EqualTo(0.0).Within(Tolerance), "GC skew 0 for empty");
            Assert.That(comp.AtSkew, Is.EqualTo(0.0).Within(Tolerance), "AT skew 0 for empty");
        });
    }

    // M9 — null behaves like empty (null-safe contract)
    [Test]
    public void CalculateNucleotideComposition_Null_ReturnsAllZeroWithoutThrowing()
    {
        var comp = SequenceStatistics.CalculateNucleotideComposition(null!);

        Assert.Multiple(() =>
        {
            Assert.That(comp.Length, Is.EqualTo(0), "Null treated as empty");
            Assert.That(comp.GcContent, Is.EqualTo(0.0).Within(Tolerance), "GC content 0 for null");
        });
    }

    // S1 — case-insensitive counting
    // Evidence: Biopython counts lowercase ("CGScgs", s.count("g"))
    [Test]
    public void CalculateNucleotideComposition_MixedCase_MatchesUpperCase()
    {
        var upper = SequenceStatistics.CalculateNucleotideComposition("ATGC");
        var lower = SequenceStatistics.CalculateNucleotideComposition("atgc");
        var mixed = SequenceStatistics.CalculateNucleotideComposition("AtGc");

        Assert.Multiple(() =>
        {
            Assert.That(lower.GcContent, Is.EqualTo(upper.GcContent).Within(Tolerance), "Lowercase GC content equals uppercase");
            Assert.That(mixed.GcContent, Is.EqualTo(upper.GcContent).Within(Tolerance), "Mixed-case GC content equals uppercase");
            Assert.That(lower.CountA, Is.EqualTo(upper.CountA), "Lowercase A count equals uppercase");
            Assert.That(lower.CountG, Is.EqualTo(upper.CountG), "Lowercase G count equals uppercase");
        });
    }

    // S2 — no G/C -> GC skew 0 (zero-denominator guard)
    // Evidence: Biopython GC_skew catches ZeroDivisionError -> 0.0
    [Test]
    public void CalculateNucleotideComposition_NoGc_ReturnsGcSkewZero()
    {
        var comp = SequenceStatistics.CalculateNucleotideComposition("AAAT");

        Assert.That(comp.GcSkew, Is.EqualTo(0.0).Within(Tolerance),
            "GC skew 0 when G+C=0 per Biopython zero-division handling");
    }

    // S3 — no A/T -> AT skew 0 (zero-denominator guard)
    [Test]
    public void CalculateNucleotideComposition_NoAt_ReturnsAtSkewZero()
    {
        var comp = SequenceStatistics.CalculateNucleotideComposition("GGGC");

        Assert.That(comp.AtSkew, Is.EqualTo(0.0).Within(Tolerance),
            "AT skew 0 when A+T=0");
    }

    // S4 — N and non-alphabet symbols counted separately
    // Evidence: standard alphabet {A,T,G,C,U}; others -> N/Other
    [Test]
    public void CalculateNucleotideComposition_WithNsAndOther_PartitionsCorrectly()
    {
        var comp = SequenceStatistics.CalculateNucleotideComposition("ATGCNNXX");

        Assert.Multiple(() =>
        {
            Assert.That(comp.CountN, Is.EqualTo(2), "Two N");
            Assert.That(comp.CountOther, Is.EqualTo(2), "Two other (X)");
            Assert.That(comp.Length, Is.EqualTo(8), "Length includes N and other");
            Assert.That(comp.GcContent, Is.EqualTo(0.5).Within(Tolerance), "GC over ACGTU only: 2/4 = 0.5");
        });
    }

    // S5 — INV-03: counts partition the sequence
    [Test]
    public void CalculateNucleotideComposition_ArbitrarySequence_CountsSumToLength()
    {
        var comp = SequenceStatistics.CalculateNucleotideComposition("ATGCUNXatgc");

        int sum = comp.CountA + comp.CountT + comp.CountG + comp.CountC +
                  comp.CountU + comp.CountN + comp.CountOther;

        Assert.That(sum, Is.EqualTo(comp.Length),
            "INV-03: per-base counts partition every character of the sequence");
    }

    #endregion

    #region SummarizeNucleotideSequence (Delegate)

    // C1 — summary delegates composition values
    [Test]
    public void SummarizeNucleotideSequence_BalancedDna_AggregatesComposition()
    {
        var summary = SequenceStatistics.SummarizeNucleotideSequence("ATGC");

        Assert.Multiple(() =>
        {
            Assert.That(summary.Length, Is.EqualTo(4), "Summary length matches composition");
            Assert.That(summary.GcContent, Is.EqualTo(0.5).Within(Tolerance), "Summary GC content delegates to composition (0.5)");
        });
    }

    #endregion

    #region CalculateAminoAcidComposition (Delegate)

    // C2 — protein residue counts and length
    [Test]
    public void CalculateAminoAcidComposition_Protein_ReturnsExactResidueCounts()
    {
        var comp = SequenceStatistics.CalculateAminoAcidComposition("MKVLWA");

        Assert.Multiple(() =>
        {
            Assert.That(comp.Length, Is.EqualTo(6), "Six residues");
            Assert.That(comp.Counts['M'], Is.EqualTo(1), "One M");
            Assert.That(comp.Counts['K'], Is.EqualTo(1), "One K");
            Assert.That(comp.Counts['V'], Is.EqualTo(1), "One V");
            Assert.That(comp.Counts['L'], Is.EqualTo(1), "One L");
            Assert.That(comp.Counts['W'], Is.EqualTo(1), "One W");
            Assert.That(comp.Counts['A'], Is.EqualTo(1), "One A");
        });
    }

    // C2b — empty protein -> zero length
    [Test]
    public void CalculateAminoAcidComposition_EmptyString_ReturnsZeroLength()
    {
        var comp = SequenceStatistics.CalculateAminoAcidComposition("");

        Assert.That(comp.Length, Is.EqualTo(0), "Empty protein has zero length");
    }

    #endregion
}
