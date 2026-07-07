// KMER-BOTH-001 — Both-Strand K-mer Counting (forward + reverse complement)
// Evidence: docs/Evidence/KMER-BOTH-001-Evidence.md
// TestSpec: tests/TestSpecs/KMER-BOTH-001.md
// Source: Anvar SY et al. (2014). Genome Biology 15:555 (kPAL "balance" = sum of each k-mer
//         and its reverse complement); Shporer S et al. (2016), BMC Genomics — inversion
//         symmetry (count[w] = forward[w] + forward[RC(w)]); Marçais G & Kingsford C (2011),
//         Bioinformatics 27(6):764–770 (k-mer window count, total = 2·(L−k+1)).

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// Tests for KMER-BOTH-001: KmerAnalyzer.CountKmersBothStrands.
///
/// Expected values are derived from the sources, not from the implementation:
/// the both-strand count of a k-mer w is forward[w] + forward[RC(w)] (kPAL balance,
/// Anvar et al. 2014; inversion symmetry, Shporer et al. 2016); the grand total over
/// all keys is 2·(L − k + 1) (Marçais & Kingsford 2011). Worked-example dictionaries
/// were computed by hand from these definitions in the Evidence artifact.
/// </summary>
[TestFixture]
public class KmerAnalyzer_CountKmersBothStrands_Tests
{
    #region CountKmersBothStrands(string) — MUST

    // M1 — Worked example ATGGC, k=2. Forward {AT,TG,GG,GC}; RC(ATGGC)=GCCAT → {GC,CC,CA,AT}.
    // Summed per key: AT:2, TG:1, GG:1, GC:2, CC:1, CA:1. (Evidence §Test Datasets.)
    [Test]
    public void CountKmersBothStrands_WorkedExampleAtggcK2_ReturnsExactDictionary()
    {
        var counts = KmerAnalyzer.CountKmersBothStrands("ATGGC", 2);

        var expected = new Dictionary<string, int>
        {
            ["AT"] = 2, ["TG"] = 1, ["GG"] = 1, ["GC"] = 2, ["CC"] = 1, ["CA"] = 1,
        };
        Assert.That(counts, Is.EquivalentTo(expected),
            "ATGGC k=2 both-strand counts must equal forward[w]+forward[RC(w)] per key " +
            "(kPAL balance / inversion symmetry); exactly these 6 keys, no others.");
    }

    // M2 — ACGT, k=2: every 2-mer here is its own reverse complement (RC(AC)=GT? no).
    // RC(ACGT)=ACGT (palindrome of the whole word) → RC-strand 2-mers {AC,CG,GT} identical to
    // forward, so each key doubles: AC:2, CG:2, GT:2. (Evidence §Test Datasets.)
    [Test]
    public void CountKmersBothStrands_PalindromicAcgtK2_DoublesEachKey()
    {
        var counts = KmerAnalyzer.CountKmersBothStrands("ACGT", 2);

        var expected = new Dictionary<string, int> { ["AC"] = 2, ["CG"] = 2, ["GT"] = 2 };
        Assert.That(counts, Is.EquivalentTo(expected),
            "ACGT k=2: RC(ACGT)=ACGT so the reverse strand yields the same 2-mers; " +
            "each count doubles to 2 and there are exactly 3 keys.");
    }

    // M3 — AAA, k=2: forward {AA:2}; RC(AAA)=TTT → {TT:2}. Combined AA:2, TT:2.
    [Test]
    public void CountKmersBothStrands_NonPalindromicAaaK2_AddsReverseComplementKey()
    {
        var counts = KmerAnalyzer.CountKmersBothStrands("AAA", 2);

        var expected = new Dictionary<string, int> { ["AA"] = 2, ["TT"] = 2 };
        Assert.That(counts, Is.EquivalentTo(expected),
            "AAA k=2: forward AA:2; reverse complement TTT contributes TT:2 " +
            "(count[w]=forward[w]+forward[RC(w)]).");
    }

    // M4 — INV-02: grand total over all keys = 2·(L − k + 1) = 2·(5−2+1) = 8 for ATGGC.
    [Test]
    public void CountKmersBothStrands_AtggcK2_GrandTotalEqualsTwiceWindowCount()
    {
        var counts = KmerAnalyzer.CountKmersBothStrands("ATGGC", 2);

        const int length = 5, k = 2;
        int expectedTotal = 2 * (length - k + 1); // both strands each contribute L−k+1 windows
        Assert.That(counts.Values.Sum(), Is.EqualTo(expectedTotal),
            "Both-strand grand total must be 2·(L−k+1)=8 (Marçais & Kingsford 2011 window count).");
    }

    // M5 — INV-01: for every key w, count[w] must equal forward[w] + forward[RC(w)],
    // computed independently from single-strand CountKmers.
    [Test]
    public void CountKmersBothStrands_AtggcK2_EachKeyEqualsForwardPlusReverseComplement()
    {
        const string seq = "ATGGC";
        const int k = 2;
        var both = KmerAnalyzer.CountKmersBothStrands(seq, k);
        var forward = KmerAnalyzer.CountKmers(seq, k);

        Assert.Multiple(() =>
        {
            foreach (var kvp in both)
            {
                string rc = DnaSequence.GetReverseComplementString(kvp.Key);
                int expected = forward.GetValueOrDefault(kvp.Key) + forward.GetValueOrDefault(rc);
                Assert.That(kvp.Value, Is.EqualTo(expected),
                    $"count[{kvp.Key}] must equal forward[{kvp.Key}]+forward[{rc}] (inversion symmetry).");
            }
        });
    }

    // M6 — INV-03: the both-strand profile is strand-symmetric: count[w] == count[RC(w)].
    [Test]
    public void CountKmersBothStrands_AtggcK2_IsStrandSymmetric()
    {
        var counts = KmerAnalyzer.CountKmersBothStrands("ATGGC", 2);

        Assert.Multiple(() =>
        {
            foreach (var kvp in counts)
            {
                string rc = DnaSequence.GetReverseComplementString(kvp.Key);
                Assert.That(counts.GetValueOrDefault(rc), Is.EqualTo(kvp.Value),
                    $"count[{kvp.Key}] must equal count[{rc}] (strand-symmetric profile, kPAL balance).");
            }
        });
    }

    #endregion

    #region CountKmersBothStrands(string) — SHOULD

    // S1 — Case-insensitivity: lowercase input yields the same dictionary as uppercase.
    [Test]
    public void CountKmersBothStrands_LowercaseInput_EqualsUppercase()
    {
        var lower = KmerAnalyzer.CountKmersBothStrands("atggc", 2);
        var upper = KmerAnalyzer.CountKmersBothStrands("ATGGC", 2);

        Assert.That(lower, Is.EquivalentTo(upper),
            "Counting is case-insensitive (input upper-cased internally), so atggc == ATGGC.");
    }

    // S2 — k = L: one window per strand. ATGC (L=4) forward {ATGC:1}; RC(ATGC)=GCAT → {GCAT:1}.
    [Test]
    public void CountKmersBothStrands_KEqualsLength_OneWindowPerStrand()
    {
        var counts = KmerAnalyzer.CountKmersBothStrands("ATGC", 4);

        var expected = new Dictionary<string, int> { ["ATGC"] = 1, ["GCAT"] = 1 };
        Assert.That(counts, Is.EquivalentTo(expected),
            "k=L gives one window per strand: forward ATGC and reverse complement GCAT, each count 1.");
    }

    // S3 — DnaSequence overload delegates to the string overload (smoke).
    [Test]
    public void CountKmersBothStrands_DnaSequenceOverload_MatchesStringOverload()
    {
        var fromString = KmerAnalyzer.CountKmersBothStrands("ATGGC", 2);
        var fromDna = KmerAnalyzer.CountKmersBothStrands(new DnaSequence("ATGGC"), 2);

        Assert.That(fromDna, Is.EquivalentTo(fromString),
            "The DnaSequence overload must delegate to the string overload and produce identical counts.");
    }

    #endregion

    #region CountKmersBothStrands(string) — COULD (edge / failure modes)

    // C1 — Empty sequence ⇒ empty dictionary (no windows).
    [Test]
    public void CountKmersBothStrands_EmptySequence_ReturnsEmpty()
    {
        var counts = KmerAnalyzer.CountKmersBothStrands(string.Empty, 2);
        Assert.That(counts, Is.Empty, "Empty sequence has no k-mer windows on either strand.");
    }

    // C2 — Null sequence ⇒ empty dictionary (null-safe).
    [Test]
    public void CountKmersBothStrands_NullSequence_ReturnsEmpty()
    {
        var counts = KmerAnalyzer.CountKmersBothStrands((string?)null!, 2);
        Assert.That(counts, Is.Empty, "Null sequence is treated as empty: no k-mers on either strand.");
    }

    // C3 — k > L ⇒ empty dictionary (L − k + 1 ≤ 0).
    [Test]
    public void CountKmersBothStrands_KGreaterThanLength_ReturnsEmpty()
    {
        var counts = KmerAnalyzer.CountKmersBothStrands("AC", 5);
        Assert.That(counts, Is.Empty, "k>L gives no windows (L−k+1≤0) on either strand.");
    }

    // C4 — k ≤ 0 ⇒ ArgumentOutOfRangeException (API contract, sibling CountKmers).
    [Test]
    public void CountKmersBothStrands_NonPositiveK_Throws()
    {
        Assert.That(() => KmerAnalyzer.CountKmersBothStrands("ACGT", 0),
            NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>(),
            "k must be positive; k=0 must throw ArgumentOutOfRangeException.");
    }

    #endregion
}
