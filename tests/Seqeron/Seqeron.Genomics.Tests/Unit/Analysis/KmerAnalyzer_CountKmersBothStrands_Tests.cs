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

    // M2 — ACGT, k=2: RC(ACGT)=ACGT (the whole word is a reverse-complement palindrome;
    // RC(AC)=GT and RC(CG)=CG) → RC-strand 2-mers {AC,CG,GT} identical to
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

    #region Reference cross-check — kPAL Profile.balance (LUMC/kPAL klib.py, run 2026-09-28)

    // Reference sequence with the EcoRI (GAATTC) and BamHI (GGATCC) reverse-complement palindromes.
    private const string KpalReferenceSequence = "GAATTCACGTTGCAGGATCCATGC"; // L = 24

    // R1 — odd k (no reverse-complement palindromes can exist for odd k). Expected dictionary is the
    // verbatim output of kPAL's own code: Profile.from_sequences([s], 3) followed by Profile.balance()
    // (Anvar et al. 2014; https://raw.githubusercontent.com/LUMC/kPAL/master/kpal/klib.py), and is
    // identical to Counter(s) + Counter(Bio.Seq.reverse_complement(s)) (Biopython). Σ = 2·(24−3+1) = 44.
    [Test]
    public void CountKmersBothStrands_KpalReference_OddK3_MatchesKpalBalance()
    {
        var counts = KmerAnalyzer.CountKmersBothStrands(KpalReferenceSequence, 3);

        var expected = new Dictionary<string, int>
        {
            ["AAC"] = 1, ["AAT"] = 2, ["ACG"] = 2, ["AGG"] = 1, ["ATC"] = 2, ["ATG"] = 2, ["ATT"] = 2,
            ["CAA"] = 1, ["CAC"] = 1, ["CAG"] = 1, ["CAT"] = 2, ["CCA"] = 1, ["CCT"] = 1, ["CGT"] = 2,
            ["CTG"] = 1, ["GAA"] = 2, ["GAT"] = 2, ["GCA"] = 3, ["GGA"] = 2, ["GTG"] = 1, ["GTT"] = 1,
            ["TCA"] = 1, ["TCC"] = 2, ["TGA"] = 1, ["TGC"] = 3, ["TGG"] = 1, ["TTC"] = 2, ["TTG"] = 1,
        };
        Assert.Multiple(() =>
        {
            Assert.That(counts, Is.EquivalentTo(expected), "Must equal kPAL balanced profile (k=3, 28 keys).");
            Assert.That(counts.Values.Sum(), Is.EqualTo(44), "Σ = 2·(L−k+1) = 44.");
        });
    }

    // R2 — even k with five reverse-complement palindromes (AATT, ACGT, CATG, GATC, TGCA), each
    // occurring once on the forward strand. kPAL balance doubles a palindrome (`elif i == i_rc:
    // counts[i] += counts[i]`), so each is 2. Expected dictionary = kPAL Profile.balance output (k=4).
    // Contrast (Jellyfish -C / canonical counting, Marçais & Kingsford 2011): the same palindromes
    // count 1 each, because every occurrence is tallied once under its canonical key — additive
    // both-strand counting is exactly 2× canonical for palindromes and equal to canonical otherwise.
    [Test]
    public void CountKmersBothStrands_KpalReference_EvenK4_MatchesKpalBalanceAndDoublesPalindromes()
    {
        var counts = KmerAnalyzer.CountKmersBothStrands(KpalReferenceSequence, 4);

        var expected = new Dictionary<string, int>
        {
            ["AACG"] = 1, ["AATT"] = 2, ["ACGT"] = 2, ["AGGA"] = 1, ["ATCC"] = 2, ["ATGC"] = 1, ["ATGG"] = 1,
            ["ATTC"] = 2, ["CAAC"] = 1, ["CACG"] = 1, ["CAGG"] = 1, ["CATG"] = 2, ["CCAT"] = 1, ["CCTG"] = 1,
            ["CGTG"] = 1, ["CGTT"] = 1, ["CTGC"] = 1, ["GAAT"] = 2, ["GATC"] = 2, ["GCAA"] = 1, ["GCAG"] = 1,
            ["GCAT"] = 1, ["GGAT"] = 2, ["GTGA"] = 1, ["GTTG"] = 1, ["TCAC"] = 1, ["TCCA"] = 1, ["TCCT"] = 1,
            ["TGAA"] = 1, ["TGCA"] = 2, ["TGGA"] = 1, ["TTCA"] = 1, ["TTGC"] = 1,
        };
        Assert.Multiple(() =>
        {
            Assert.That(counts, Is.EquivalentTo(expected), "Must equal kPAL balanced profile (k=4, 33 keys).");
            Assert.That(counts.Values.Sum(), Is.EqualTo(42), "Σ = 2·(L−k+1) = 42.");
            foreach (var palindrome in new[] { "AATT", "ACGT", "CATG", "GATC", "TGCA" })
                Assert.That(counts[palindrome], Is.EqualTo(2),
                    $"{palindrome} is its own reverse complement: kPAL balance doubles it (Jellyfish -C would give 1).");
        });
    }

    // R3 — k=6: the restriction-site palindromes GAATTC and GGATCC each occur once and are doubled
    // (kPAL balance output: GAATTC:2, GGATCC:2; 36 keys, Σ = 2·(24−6+1) = 38).
    [Test]
    public void CountKmersBothStrands_KpalReference_K6_RestrictionSitePalindromesDoubled()
    {
        var counts = KmerAnalyzer.CountKmersBothStrands(KpalReferenceSequence, 6);

        Assert.Multiple(() =>
        {
            Assert.That(counts["GAATTC"], Is.EqualTo(2));
            Assert.That(counts["GGATCC"], Is.EqualTo(2));
            Assert.That(counts, Has.Count.EqualTo(36));
            Assert.That(counts.Values.Sum(), Is.EqualTo(38));
            Assert.That(counts.Where(kv => kv.Key is not ("GAATTC" or "GGATCC")).All(kv => kv.Value == 1), Is.True,
                "All other 6-mers occur once across both strands (kPAL balance output).");
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

    // S4 — Repository convention (not a kPAL behaviour: kPAL splits on non-ACGT and skips such
    // k-mers): ambiguity codes are kept as literal keys (inherited from CountKmers, KMER-COUNT-001)
    // and complemented by the canonical IUPAC GetComplementBase (N↔N, R↔Y). "AAN" → forward {AA, AN};
    // RC = "NTT" → {NT, TT}. "ARC" (R=A|G) → forward {AR, RC}; RC = "GYT" → {GY, YT}.
    [Test]
    public void CountKmersBothStrands_IupacCodes_KeptAsLiteralKeysWithIupacComplement()
    {
        Assert.Multiple(() =>
        {
            Assert.That(KmerAnalyzer.CountKmersBothStrands("AAN", 2), Is.EquivalentTo(
                new Dictionary<string, int> { ["AA"] = 1, ["AN"] = 1, ["NT"] = 1, ["TT"] = 1 }));
            Assert.That(KmerAnalyzer.CountKmersBothStrands("arc", 2), Is.EquivalentTo(
                new Dictionary<string, int> { ["AR"] = 1, ["RC"] = 1, ["GY"] = 1, ["YT"] = 1 }));
        });
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
