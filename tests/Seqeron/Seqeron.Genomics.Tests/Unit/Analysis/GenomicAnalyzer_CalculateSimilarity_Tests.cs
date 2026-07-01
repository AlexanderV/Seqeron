// GENOMIC-SIMILARITY-001 — Sequence Similarity (k-mer Jaccard index)
// Evidence: docs/Evidence/GENOMIC-SIMILARITY-001-Evidence.md
// TestSpec: tests/TestSpecs/GENOMIC-SIMILARITY-001.md
// Source: Jaccard P. (1901). Bull. Soc. Vaudoise Sci. Nat. 37(142):547–579 (J=|A∩B|/|A∪B|);
//         Ondov et al. (2016). Genome Biology 17:132 (k-mer-set Jaccard). DOI: 10.1186/s13059-016-0997-x.

using System;
using NUnit.Framework;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GenomicAnalyzer_CalculateSimilarity_Tests
{
    // The method reports J×100 (percentage). Jaccard ∈ [0,1] (Jaccard 1901).
    private const double Tolerance = 1e-10;

    #region CalculateSimilarity

    // M1 — Partial overlap, exact fraction (Jaccard 1901; Ondov 2016).
    // k=3 distinct sets: ACGTACGT -> {ACG,CGT,GTA,TAC}; ACGTACGA -> {ACG,CGT,GTA,TAC,CGA}.
    // |A∩B|=4, |A∪B|=5 -> J=4/5=0.8 -> 80.0. (Identity impl would give 100; off-by-one in k would differ.)
    [Test]
    public void CalculateSimilarity_PartialOverlap_ReturnsExactJaccardPercent()
    {
        var seq1 = new DnaSequence("ACGTACGT");
        var seq2 = new DnaSequence("ACGTACGA");

        double result = GenomicAnalyzer.CalculateSimilarity(seq1, seq2, kmerSize: 3);

        Assert.That(result, Is.EqualTo(80.0).Within(Tolerance),
            "4 shared distinct 3-mers over 5 distinct 3-mers in union = 4/5 = 0.8 -> 80.0 (Jaccard 1901).");
    }

    // M2 — Identical sequences -> J=1 -> 100.0 (Jaccard 1901: A=B => J=1).
    [Test]
    public void CalculateSimilarity_IdenticalSequences_Returns100()
    {
        var seq1 = new DnaSequence("ACGTACGT");
        var seq2 = new DnaSequence("ACGTACGT");

        double result = GenomicAnalyzer.CalculateSimilarity(seq1, seq2, kmerSize: 3);

        Assert.That(result, Is.EqualTo(100.0).Within(Tolerance),
            "Identical sequences share all k-mers: A=B so |A∩B|=|A∪B|, J=1 -> 100.0.");
    }

    // M3 — Disjoint k-mer sets -> J=0 -> 0.0 (Jaccard 1901: disjoint => J=0).
    // AAAAA -> {AAA}; CCCCC -> {CCC}; intersection empty, union {AAA,CCC}.
    [Test]
    public void CalculateSimilarity_DisjointKmerSets_Returns0()
    {
        var seq1 = new DnaSequence("AAAAA");
        var seq2 = new DnaSequence("CCCCC");

        double result = GenomicAnalyzer.CalculateSimilarity(seq1, seq2, kmerSize: 3);

        Assert.That(result, Is.EqualTo(0.0).Within(Tolerance),
            "{AAA} and {CCC} are disjoint: |A∩B|=0, |A∪B|=2, J=0 -> 0.0.");
    }

    // M4 — Non-integer fraction (Jaccard formula). ACGT -> {ACG,CGT}; ACGA -> {ACG,CGA}.
    // |A∩B|=1 (ACG), |A∪B|=3 -> J=1/3 -> 100/3 = 33.333... Exact value, not a rounded one.
    [Test]
    public void CalculateSimilarity_NonIntegerFraction_ReturnsExactThird()
    {
        var seq1 = new DnaSequence("ACGT");
        var seq2 = new DnaSequence("ACGA");

        double result = GenomicAnalyzer.CalculateSimilarity(seq1, seq2, kmerSize: 3);

        Assert.That(result, Is.EqualTo(100.0 / 3.0).Within(Tolerance),
            "One shared 3-mer (ACG) over union {ACG,CGT,CGA}=3 -> 1/3 -> 100/3 = 33.333...");
    }

    // M5 — Distinct-set semantics (Ondov 2016: distinct k-mers; within-sequence repeats counted once).
    // AAAAAA -> {AAA}; AAAA -> {AAA}. Both sets equal => J=1 even though k-mer counts differ (4 vs 2).
    // A multiset/bag implementation would NOT return 100 here, so this discriminates set vs bag.
    [Test]
    public void CalculateSimilarity_RepeatedKmers_TreatedAsSet_Returns100()
    {
        var seq1 = new DnaSequence("AAAAAA");
        var seq2 = new DnaSequence("AAAA");

        double result = GenomicAnalyzer.CalculateSimilarity(seq1, seq2, kmerSize: 3);

        Assert.That(result, Is.EqualTo(100.0).Within(Tolerance),
            "Both reduce to the distinct set {AAA}; set Jaccard J=1 -> 100.0 despite differing k-mer counts.");
    }

    // M6 — Symmetry (Jaccard 1901: ∩,∪ commutative). f(a,b,k) == f(b,a,k).
    [Test]
    public void CalculateSimilarity_Symmetric_OrderIndependent()
    {
        var seq1 = new DnaSequence("ACGTACGT");
        var seq2 = new DnaSequence("ACGTACGA");

        double ab = GenomicAnalyzer.CalculateSimilarity(seq1, seq2, kmerSize: 3);
        double ba = GenomicAnalyzer.CalculateSimilarity(seq2, seq1, kmerSize: 3);

        Assert.Multiple(() =>
        {
            Assert.That(ab, Is.EqualTo(80.0).Within(Tolerance), "Forward order is 80.0 (4/5).");
            Assert.That(ba, Is.EqualTo(ab).Within(Tolerance), "Jaccard is symmetric: f(b,a)=f(a,b).");
        });
    }

    // S1 — Both empty: empty union, Jaccard undefined (Jaccard 1901: non-empty sets only).
    // Implementation contract (ASM-1): returns 0.0.
    [Test]
    public void CalculateSimilarity_BothEmpty_Returns0()
    {
        var seq1 = new DnaSequence("");
        var seq2 = new DnaSequence("");

        double result = GenomicAnalyzer.CalculateSimilarity(seq1, seq2, kmerSize: 3);

        Assert.That(result, Is.EqualTo(0.0).Within(Tolerance),
            "Both empty -> empty k-mer sets -> empty union (Jaccard undefined); impl returns 0.0 (ASM-1).");
    }

    // S2 — Both shorter than k: no k-mers extractable, empty union -> 0.0 (ASM-1).
    [Test]
    public void CalculateSimilarity_BothShorterThanK_Returns0()
    {
        var seq1 = new DnaSequence("AC");
        var seq2 = new DnaSequence("GT");

        double result = GenomicAnalyzer.CalculateSimilarity(seq1, seq2, kmerSize: 3);

        Assert.That(result, Is.EqualTo(0.0).Within(Tolerance),
            "Sequences shorter than k yield empty k-mer sets -> empty union -> 0.0 (ASM-1).");
    }

    // S3 — Null first sequence -> ArgumentNullException (documented failure mode).
    [Test]
    public void CalculateSimilarity_NullSequence1_Throws()
    {
        var seq2 = new DnaSequence("ACGT");

        Assert.Throws<ArgumentNullException>(
            () => GenomicAnalyzer.CalculateSimilarity(null!, seq2, kmerSize: 3),
            "A null first sequence must raise ArgumentNullException.");
    }

    // S4 — Null second sequence -> ArgumentNullException.
    [Test]
    public void CalculateSimilarity_NullSequence2_Throws()
    {
        var seq1 = new DnaSequence("ACGT");

        Assert.Throws<ArgumentNullException>(
            () => GenomicAnalyzer.CalculateSimilarity(seq1, null!, kmerSize: 3),
            "A null second sequence must raise ArgumentNullException.");
    }

    // S5 — kmerSize below 1 is meaningless (a k-mer needs length >= 1) -> ArgumentOutOfRangeException.
    [Test]
    public void CalculateSimilarity_InvalidKmerSize_Throws()
    {
        var seq1 = new DnaSequence("ACGT");
        var seq2 = new DnaSequence("ACGA");

        Assert.Throws<ArgumentOutOfRangeException>(
            () => GenomicAnalyzer.CalculateSimilarity(seq1, seq2, kmerSize: 0),
            "kmerSize must be at least 1; 0 must raise ArgumentOutOfRangeException.");
    }

    // C1 — Range invariant INV-01: every result in [0,100] (Jaccard ∈ [0,1] scaled ×100).
    [Test]
    public void CalculateSimilarity_VariedInputs_AlwaysInRange()
    {
        var pairs = new[]
        {
            ("ACGTACGT", "ACGTACGA"),
            ("AAAAA", "CCCCC"),
            ("ACGTGGTACC", "TTACGTGGAA"),
            ("ACGTACGT", "ACGTACGT"),
        };

        Assert.Multiple(() =>
        {
            foreach (var (a, b) in pairs)
            {
                double r = GenomicAnalyzer.CalculateSimilarity(new DnaSequence(a), new DnaSequence(b), kmerSize: 3);
                Assert.That(r, Is.InRange(0.0, 100.0),
                    $"Jaccard×100 must lie in [0,100] for ({a},{b}); got {r}.");
            }
        });
    }

    // C2 — One side empty, other non-empty: empty intersection over non-empty union -> J=0 -> 0.0.
    [Test]
    public void CalculateSimilarity_OneSideEmpty_Returns0()
    {
        var seq1 = new DnaSequence("ACGTAC");
        var seq2 = new DnaSequence("");

        double result = GenomicAnalyzer.CalculateSimilarity(seq1, seq2, kmerSize: 3);

        Assert.That(result, Is.EqualTo(0.0).Within(Tolerance),
            "One empty k-mer set: |A∩B|=0 over non-empty union -> J=0 -> 0.0.");
    }

    #endregion
}
