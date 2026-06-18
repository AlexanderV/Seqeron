// ONCO-SIG-001 — SBS-96 Trinucleotide Context Catalog (pyrimidine-strand folding)
// Evidence: docs/Evidence/ONCO-SIG-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-SIG-001.md
// Source: Alexandrov L.B. et al. (2013). Nature 500(7463):415-421. https://www.nature.com/articles/nature12477
//         COSMIC SBS96. https://cancer.sanger.ac.uk/signatures/sbs/sbs96/
//         Bergstrom E.N. et al. (2019). SigProfilerMatrixGenerator. BMC Genomics 20:685. https://pmc.ncbi.nlm.nih.gov/articles/PMC6717374/
//         Complementarity (molecular biology), Watson-Crick A<->T, C<->G. https://en.wikipedia.org/wiki/Complementarity_(molecular_biology)

using System.Collections.Generic;
using System.Linq;
using Seqeron.Genomics.Oncology;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class OncologyAnalyzer_ClassifySbsContext_Tests
{
    #region ClassifySbsContext (pyrimidine substitutions, unchanged)

    // M1 — Pyrimidine C>A: already on the pyrimidine strand, kept as-is (COSMIC SBS96).
    [Test]
    public void ClassifySbsContext_PyrimidineCtoA_KeepsContext()
    {
        Assert.That(OncologyAnalyzer.ClassifySbsContext('A', 'C', 'A', 'A'), Is.EqualTo("A[C>A]A"),
            "A C>A mutation with A/A flanks is a canonical pyrimidine SBS and must be left unfolded.");
    }

    // M2 — Pyrimidine C>T with asymmetric flanks: order 5'..3' must be preserved (COSMIC SBS96).
    [Test]
    public void ClassifySbsContext_PyrimidineCtoT_PreservesFlankOrder()
    {
        Assert.That(OncologyAnalyzer.ClassifySbsContext('T', 'C', 'T', 'G'), Is.EqualTo("T[C>T]G"),
            "5'=T and 3'=G must be reported in order; a swapped/off-by-one impl would give G[C>T]T.");
    }

    // M3 — Pyrimidine T>C: the second pyrimidine reference base (COSMIC SBS96).
    [Test]
    public void ClassifySbsContext_PyrimidineTtoC_KeepsContext()
    {
        Assert.That(OncologyAnalyzer.ClassifySbsContext('G', 'T', 'C', 'A'), Is.EqualTo("G[T>C]A"),
            "T>C is a pyrimidine substitution and must not be reverse-complemented.");
    }

    #endregion

    #region ClassifySbsContext (purine reference, reverse-complement folding)

    // M4 — Purine G>T at 5'-T G A-3' folds to T[C>A]A (SigProfiler revcomp rule; worked example).
    [Test]
    public void ClassifySbsContext_PurineGtoT_FoldsToReverseComplement()
    {
        Assert.That(OncologyAnalyzer.ClassifySbsContext('T', 'G', 'T', 'A'), Is.EqualTo("T[C>A]A"),
            "G is a purine: revcomp(TGA)=TCA, G>T->C>A => T[C>A]A (Bergstrom 2019 folding rule).");
    }

    // M5 — Purine A>G at 5'-C A T-3' folds to A[T>C]G (SigProfiler revcomp rule; worked example).
    [Test]
    public void ClassifySbsContext_PurineAtoG_FoldsToReverseComplement()
    {
        Assert.That(OncologyAnalyzer.ClassifySbsContext('C', 'A', 'G', 'T'), Is.EqualTo("A[T>C]G"),
            "A is a purine: revcomp(CAT)=ATG, A>G->T>C => A[T>C]G.");
    }

    // M6 — Purine G>C at 5'-G G C-3' folds to G[C>G]C (SigProfiler revcomp rule; worked example).
    [Test]
    public void ClassifySbsContext_PurineGtoC_FoldsToReverseComplement()
    {
        Assert.That(OncologyAnalyzer.ClassifySbsContext('G', 'G', 'C', 'C'), Is.EqualTo("G[C>G]C"),
            "G is a purine: revcomp(GGC)=GCC, G>C->C>G => G[C>G]C.");
    }

    // M7 — Purine A>T at 5'-A A A-3' folds to T[T>A]T (SigProfiler revcomp rule; worked example).
    [Test]
    public void ClassifySbsContext_PurineAtoT_FoldsToReverseComplement()
    {
        Assert.That(OncologyAnalyzer.ClassifySbsContext('A', 'A', 'T', 'A'), Is.EqualTo("T[T>A]T"),
            "A is a purine: revcomp(AAA)=TTT, A>T->T>A => T[T>A]T.");
    }

    // M4b — Purine G>A at 5'-A G C-3' folds to G[C>T]T (covers the G>A->C>T fold; Wikipedia/SigProfiler).
    [Test]
    public void ClassifySbsContext_PurineGtoA_FoldsToReverseComplement()
    {
        Assert.That(OncologyAnalyzer.ClassifySbsContext('A', 'G', 'A', 'C'), Is.EqualTo("G[C>T]T"),
            "G is a purine: revcomp(AGC)=GCT, G>A->C>T => G[C>T]T (Wikipedia: G>A counted as its C>T equivalent).");
    }

    // M4c — Purine A>C at 5'-G A T-3' folds to A[T>G]C (covers the A>C->T>G fold; Wikipedia/SigProfiler).
    [Test]
    public void ClassifySbsContext_PurineAtoC_FoldsToReverseComplement()
    {
        Assert.That(OncologyAnalyzer.ClassifySbsContext('G', 'A', 'C', 'T'), Is.EqualTo("A[T>G]C"),
            "A is a purine: revcomp(GAT)=ATC, A>C->T>G => A[T>G]C (Wikipedia: A>C counted as its T>G equivalent).");
    }

    // C1 — Lower-case input classifies identically (bases are upper-cased; robustness).
    [Test]
    public void ClassifySbsContext_LowerCaseBases_ClassifySameAsUpperCase()
    {
        Assert.That(OncologyAnalyzer.ClassifySbsContext('a', 'c', 'a', 'a'), Is.EqualTo("A[C>A]A"),
            "Lower-case bases must be normalised to upper-case before classification.");
    }

    #endregion

    #region ClassifySbsContext (validation)

    // S3 — Non-ACGT context base has no defined trinucleotide context.
    [Test]
    public void ClassifySbsContext_NonAcgtContextBase_Throws()
    {
        Assert.Throws<System.ArgumentException>(
            () => OncologyAnalyzer.ClassifySbsContext('N', 'C', 'A', 'A'),
            "A flanking base that is not A/C/G/T cannot define a trinucleotide context.");
    }

    // S4 — ref == alt is not a substitution.
    [Test]
    public void ClassifySbsContext_ReferenceEqualsAlternate_Throws()
    {
        Assert.Throws<System.ArgumentException>(
            () => OncologyAnalyzer.ClassifySbsContext('A', 'C', 'C', 'A'),
            "reference == alternate is not a mutation and must be rejected.");
    }

    // S5 — Invalid (non-ACGT) reference/alternate base.
    [Test]
    public void ClassifySbsContext_InvalidBase_Throws()
    {
        Assert.Throws<System.ArgumentException>(
            () => OncologyAnalyzer.ClassifySbsContext('A', 'X', 'A', 'A'),
            "A reference base that is not A/C/G/T must be rejected.");
    }

    // S5b — Invalid (non-ACGT) alternate base is rejected.
    [Test]
    public void ClassifySbsContext_InvalidAlternateBase_Throws()
    {
        Assert.Throws<System.ArgumentException>(
            () => OncologyAnalyzer.ClassifySbsContext('A', 'C', 'Z', 'A'),
            "An alternate base that is not A/C/G/T must be rejected.");
    }

    // S5c — Invalid (non-ACGT) 3' flanking base is rejected.
    [Test]
    public void ClassifySbsContext_Invalid3PrimeBase_Throws()
    {
        Assert.Throws<System.ArgumentException>(
            () => OncologyAnalyzer.ClassifySbsContext('A', 'C', 'A', 'N'),
            "A 3' flanking base that is not A/C/G/T must be rejected.");
    }

    #endregion

    #region EnumerateSbs96Channels

    // M8 — The channel space is exactly the 96 canonical pyrimidine labels (6 x 4 x 4; INV-02).
    [Test]
    public void EnumerateSbs96Channels_Returns96DistinctPyrimidineChannels()
    {
        var channels = OncologyAnalyzer.EnumerateSbs96Channels();

        Assert.Multiple(() =>
        {
            Assert.That(channels.Count, Is.EqualTo(96),
                "SBS-96 has exactly 6 substitutions x 4 5'-bases x 4 3'-bases = 96 channels (Alexandrov 2013).");
            Assert.That(channels.Distinct().Count(), Is.EqualTo(96),
                "All 96 channel labels must be distinct.");
            // INV-01: every reference (mutated) base must be a pyrimidine C or T.
            Assert.That(channels.All(c => c[2] is 'C' or 'T'), Is.True,
                "Every channel's reference base (position 2 in 5'[REF>ALT]3') must be a pyrimidine.");
            // Exactly the six pyrimidine substitutions, 16 contexts each.
            var subs = channels.Select(c => c.Substring(2, 3)).Distinct().OrderBy(s => s).ToArray();
            Assert.That(subs, Is.EqualTo(new[] { "C>A", "C>G", "C>T", "T>A", "T>C", "T>G" }),
                "The six substitution subtypes must be exactly the COSMIC pyrimidine set.");
            Assert.That(channels.Contains("A[C>A]A") && channels.Contains("T[T>A]T"), Is.True,
                "Spot-check that representative folded/unfolded labels are members of the channel set.");
        });
    }

    #endregion

    #region Build96ContextCatalog

    // M9 — Catalog tallies match a hand computed multiset; sum equals the variant count (INV-03).
    [Test]
    public void Build96ContextCatalog_TalliesVariants_MatchesHandCount()
    {
        var variants = new (char, char, char, char)[]
        {
            ('A', 'C', 'A', 'A'),  // A[C>A]A
            ('A', 'C', 'A', 'A'),  // A[C>A]A  (duplicate)
            ('T', 'C', 'T', 'G'),  // T[C>T]G
            ('T', 'G', 'T', 'A'),  // purine G>T -> folds to T[C>A]A
        };

        var catalog = OncologyAnalyzer.Build96ContextCatalog(variants);

        Assert.Multiple(() =>
        {
            Assert.That(catalog.Count, Is.EqualTo(96), "All 96 channels must be present in the spectrum.");
            Assert.That(catalog["A[C>A]A"], Is.EqualTo(2), "Two variants fall in A[C>A]A.");
            Assert.That(catalog["T[C>T]G"], Is.EqualTo(1), "One variant falls in T[C>T]G.");
            Assert.That(catalog["T[C>A]A"], Is.EqualTo(1), "The folded purine G>T falls in T[C>A]A.");
            Assert.That(catalog.Values.Sum(), Is.EqualTo(variants.Length),
                "INV-03: the sum of channel counts equals the number of classifiable variants.");
            Assert.That(catalog.Values.Count(v => v != 0), Is.EqualTo(3),
                "Exactly three distinct channels are populated by this input.");
        });
    }

    // M10 — A purine variant and its pyrimidine-strand equivalent are co-counted in one channel (INV-04).
    [Test]
    public void Build96ContextCatalog_PurineAndPyrimidineForms_CoCountInSameChannel()
    {
        var variants = new (char, char, char, char)[]
        {
            ('T', 'G', 'T', 'A'),  // purine G>T -> T[C>A]A
            ('T', 'C', 'A', 'A'),  // pyrimidine form, already T[C>A]A
        };

        var catalog = OncologyAnalyzer.Build96ContextCatalog(variants);

        Assert.Multiple(() =>
        {
            Assert.That(catalog["T[C>A]A"], Is.EqualTo(2),
                "INV-04: G>T and its reverse-complement C>A form map to the same channel and must co-count.");
            Assert.That(catalog.Values.Sum(), Is.EqualTo(2), "Both variants are counted exactly once.");
        });
    }

    // S1 — Null input.
    [Test]
    public void Build96ContextCatalog_NullVariants_Throws()
    {
        Assert.Throws<System.ArgumentNullException>(
            () => OncologyAnalyzer.Build96ContextCatalog(null!),
            "A null variant collection must be rejected.");
    }

    // S2 — Empty input yields all 96 channels with count 0.
    [Test]
    public void Build96ContextCatalog_EmptyVariants_AllChannelsZero()
    {
        var catalog = OncologyAnalyzer.Build96ContextCatalog(
            new List<(char, char, char, char)>());

        Assert.Multiple(() =>
        {
            Assert.That(catalog.Count, Is.EqualTo(96),
                "The empty spectrum still has the fixed 96-channel shape.");
            Assert.That(catalog.Values.Sum(), Is.EqualTo(0), "An empty input produces zero counts everywhere.");
        });
    }

    #endregion
}
