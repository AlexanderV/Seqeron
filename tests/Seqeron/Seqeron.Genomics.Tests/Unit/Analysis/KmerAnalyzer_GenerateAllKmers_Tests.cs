// KMER-GENERATE-001 — K-mer Generation (all possible k-mers over an alphabet)
// Evidence: docs/Evidence/KMER-GENERATE-001-Evidence.md
// TestSpec: tests/TestSpecs/KMER-GENERATE-001.md
// Source: Wikipedia — K-mer (https://en.wikipedia.org/wiki/K-mer);
//         Clavijo BJ (2018), BioInfoLogics — k-mer counting, part I
//         (https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/);
//         Python Std Library — itertools.product (https://docs.python.org/3/library/itertools.html);
//         Rosalind LEXF — Enumerating k-mers Lexicographically (https://rosalind.info/problems/lexf/)

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// Tests for KMER-GENERATE-001: KmerAnalyzer.GenerateAllKmers.
///
/// Expected values are derived from the sources, not from the implementation:
/// the k-mer universe over an n-letter alphabet has n^k members (Wikipedia/BioInfoLogics),
/// enumerated as the k-fold Cartesian product (itertools.product); for a sorted alphabet
/// the output is lexicographic with the rightmost position advancing fastest.
/// </summary>
[TestFixture]
public class KmerAnalyzer_GenerateAllKmers_Tests
{
    // Default DNA alphabet {A,C,G,T} is already in sorted order -> lexicographic output.
    private const string Dna = "ACGT";

    // 20 standard amino acids, sorted; |alphabet| = 20 -> universe size 20^k.
    private const string Protein = "ACDEFGHIKLMNPQRSTVWY";

    #region GenerateAllKmers — MUST

    // M1 — n^k with n=4, k=1: the four DNA monomers A,C,G,T in lexicographic order.
    [Test]
    public void GenerateAllKmers_K1Dna_ReturnsFourMonomersInOrder()
    {
        var result = KmerAnalyzer.GenerateAllKmers(1).ToList();

        var expected = new[] { "A", "C", "G", "T" };
        Assert.That(result, Is.EqualTo(expected),
            "k=1 over {A,C,G,T} must be exactly the 4 monomers (4^1=4) in lexicographic order.");
    }

    // M2 — n^k with n=4, k=2: all 16 (4^2) 2-mers, lexicographic odometer order AA..TT.
    [Test]
    public void GenerateAllKmers_K2Dna_ReturnsSixteenTwoMersLexicographic()
    {
        var result = KmerAnalyzer.GenerateAllKmers(2).ToList();

        // Cartesian product {A,C,G,T} x {A,C,G,T}, rightmost position fastest (itertools.product).
        var expected = new[]
        {
            "AA", "AC", "AG", "AT",
            "CA", "CC", "CG", "CT",
            "GA", "GC", "GG", "GT",
            "TA", "TC", "TG", "TT"
        };
        Assert.That(result, Is.EqualTo(expected),
            "k=2 must be all 16 (4^2) 2-mers in lexicographic order per the Cartesian-product odometer ordering.");
    }

    // M3 — k=3: count is 4^3=64; boundaries first=AAA, second=AAC, last=TTT (odometer order).
    [Test]
    public void GenerateAllKmers_K3Dna_HasSixtyFourWithCorrectBoundaries()
    {
        var result = KmerAnalyzer.GenerateAllKmers(3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(result.Count, Is.EqualTo(64),
                "k=3 universe size is 4^3 = 64 (Wikipedia n^k).");
            Assert.That(result[0], Is.EqualTo("AAA"),
                "First k-mer in lexicographic order is AAA.");
            Assert.That(result[1], Is.EqualTo("AAC"),
                "Second is AAC: rightmost position advances first (odometer ordering).");
            Assert.That(result[^1], Is.EqualTo("TTT"),
                "Last k-mer in lexicographic order is TTT.");
        });
    }

    // M4 — universe size equals 4^k for k = 1..6 (default DNA alphabet).
    [Test]
    public void GenerateAllKmers_DnaVariousK_CountEqualsFourToTheK()
    {
        Assert.Multiple(() =>
        {
            for (int k = 1; k <= 6; k++)
            {
                int expected = (int)Math.Pow(4, k); // 4,16,64,256,1024,4096
                Assert.That(KmerAnalyzer.GenerateAllKmers(k).Count(), Is.EqualTo(expected),
                    $"k={k}: number of all possible DNA k-mers must be 4^{k} = {expected}.");
            }
        });
    }

    // M5 — n^k generalises to any alphabet: 20 amino acids, k=2 -> 20^2 = 400.
    [Test]
    public void GenerateAllKmers_ProteinAlphabetK2_ReturnsFourHundred()
    {
        var result = KmerAnalyzer.GenerateAllKmers(2, Protein).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(result.Count, Is.EqualTo(400),
                "Universe size for a 20-letter alphabet at k=2 is 20^2 = 400 (n^k).");
            Assert.That(result.Distinct().Count(), Is.EqualTo(400),
                "All 400 protein 2-mers must be distinct.");
        });
    }

    // M6 — output is exactly the Cartesian-product set: no duplicates, distinct count = 4^k.
    [Test]
    public void GenerateAllKmers_K4Dna_NoDuplicatesEqualsCartesianSet()
    {
        var result = KmerAnalyzer.GenerateAllKmers(4).ToList();

        // Independent reference: build the 4-fold Cartesian product of {A,C,G,T}.
        var reference =
            from a in Dna
            from b in Dna
            from c in Dna
            from d in Dna
            select string.Concat(a, b, c, d);
        var referenceSet = reference.ToHashSet();

        Assert.Multiple(() =>
        {
            Assert.That(result.Count, Is.EqualTo(256),
                "k=4 universe size is 4^4 = 256.");
            Assert.That(result.Distinct().Count(), Is.EqualTo(256),
                "Output must contain no duplicate k-mers (each is a unique length-k tuple).");
            Assert.That(result.ToHashSet(), Is.EquivalentTo(referenceSet),
                "Output set must equal the independently built 4-fold Cartesian product of {A,C,G,T}.");
        });
    }

    #endregion

    #region GenerateAllKmers — SHOULD (edge cases)

    // S1 — single-letter alphabet: 1^k = 1, the homopolymer only.
    [Test]
    public void GenerateAllKmers_SingleLetterAlphabet_ReturnsOnlyHomopolymer()
    {
        var result = KmerAnalyzer.GenerateAllKmers(4, "A").ToList();

        Assert.That(result, Is.EqualTo(new[] { "AAAA" }),
            "A 1-letter alphabet yields exactly one k-mer (1^4 = 1): the homopolymer AAAA.");
    }

    // S2 — INV-03: every k-mer has length exactly k and uses only alphabet characters.
    [Test]
    public void GenerateAllKmers_K2Dna_EveryKmerHasLengthKFromAlphabet()
    {
        var alphabetSet = Dna.ToHashSet();

        var result = KmerAnalyzer.GenerateAllKmers(2).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(result.All(km => km.Length == 2), Is.True,
                "Every generated k-mer must have length exactly k=2.");
            Assert.That(result.All(km => km.All(alphabetSet.Contains)), Is.True,
                "Every character of every k-mer must come from the alphabet {A,C,G,T}.");
        });
    }

    // S3 — k must be positive: k=0 and k<0 throw ArgumentOutOfRangeException.
    [Test]
    public void GenerateAllKmers_NonPositiveK_ThrowsArgumentOutOfRange()
    {
        Assert.Multiple(() =>
        {
            // Enumerate to trigger the validation (deferred-execution method).
            Assert.Throws<ArgumentOutOfRangeException>(
                () => KmerAnalyzer.GenerateAllKmers(0).ToList(),
                "k=0 is not a valid k-mer length.");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => KmerAnalyzer.GenerateAllKmers(-1).ToList(),
                "Negative k is not a valid k-mer length.");
        });
    }

    // S4 — null/empty alphabet: no symbols means no k-mers can be formed -> ArgumentException.
    [Test]
    public void GenerateAllKmers_NullOrEmptyAlphabet_ThrowsArgumentException()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(
                () => KmerAnalyzer.GenerateAllKmers(2, "").ToList(),
                "An empty alphabet cannot form any k-mer.");
            Assert.Throws<ArgumentException>(
                () => KmerAnalyzer.GenerateAllKmers(2, null!).ToList(),
                "A null alphabet cannot form any k-mer.");
        });
    }

    #endregion

    #region GenerateAllKmers — COULD

    // C1 — ordering follows the supplied alphabet, not a forced sort: "TGCA", k=1 -> T,G,C,A.
    [Test]
    public void GenerateAllKmers_UnsortedAlphabetK1_FollowsAlphabetOrder()
    {
        var result = KmerAnalyzer.GenerateAllKmers(1, "TGCA").ToList();

        Assert.That(result, Is.EqualTo(new[] { "T", "G", "C", "A" }),
            "Output order follows the alphabet's own order; lexicographic order holds only for a sorted alphabet.");
    }

    #endregion

    #region GenerateAllKmers — Rosalind LEXF / itertools.product cross-checks (review 2026-09)

    // C2 — Rosalind LEXF sample dataset: ordered alphabet "T A G C", n=2. Lexicographic order is
    // defined by the alphabet's own order (T < A < G < C). Sample output (Rosalind LEXF, as archived
    // in github.com/mtarbit/Rosalind-Problems e015-lexf.py); identical to
    // [''.join(p) for p in itertools.product("TAGC", repeat=2)].
    [Test]
    public void GenerateAllKmers_RosalindLexfSample_TagcN2_MatchesSampleOutput()
    {
        var result = KmerAnalyzer.GenerateAllKmers(2, "TAGC").ToList();

        var expected = new[]
        {
            "TT", "TA", "TG", "TC", "AT", "AA", "AG", "AC",
            "GT", "GA", "GG", "GC", "CT", "CA", "CG", "CC",
        };
        Assert.That(result, Is.EqualTo(expected),
            "Rosalind LEXF sample (T A G C, n=2) must be reproduced exactly, in the alphabet's order.");
    }

    // C3 — LEXF order definition for an unsorted alphabet at k=3: s <Lex t iff the first mismatching
    // symbol of s precedes that of t in the ALPHABET order. Checked pairwise against that definition.
    [Test]
    public void GenerateAllKmers_UnsortedAlphabetK3_IsLexicographicInAlphabetOrder()
    {
        const string alphabet = "GTCA";
        var rank = alphabet.Select((c, i) => (c, i)).ToDictionary(t => t.c, t => t.i);
        var result = KmerAnalyzer.GenerateAllKmers(3, alphabet).ToList();

        bool LexLess(string s, string t)
        {
            for (int j = 0; j < s.Length; j++)
                if (s[j] != t[j])
                    return rank[s[j]] < rank[t[j]];
            return false;
        }

        Assert.Multiple(() =>
        {
            Assert.That(result.Count, Is.EqualTo(64), "4^3 = 64.");
            Assert.That(result[0], Is.EqualTo("GGG"));
            Assert.That(result[1], Is.EqualTo("GGT"));
            Assert.That(result[^1], Is.EqualTo("AAA"));
            for (int i = 1; i < result.Count; i++)
                Assert.That(LexLess(result[i - 1], result[i]), Is.True,
                    $"{result[i - 1]} must precede {result[i]} in LEXF order over (G,T,C,A).");
        });
    }

    // C4 — itertools.product semantics for a repeated symbol: product("AAC", repeat=2) =
    // AA AA AC AA AA AC CA CA CC (Python 3 reference output). Count stays |alphabet|^k = 9;
    // output is duplicate-free only when the alphabet symbols are distinct.
    [Test]
    public void GenerateAllKmers_RepeatedSymbolAlphabet_MatchesItertoolsProduct()
    {
        var result = KmerAnalyzer.GenerateAllKmers(2, "AAC").ToList();

        Assert.That(result, Is.EqualTo(new[] { "AA", "AA", "AC", "AA", "AA", "AC", "CA", "CA", "CC" }));
    }

    // C5 — O(k) working space (doc §4.3): a very long homopolymer k-mer is produced directly.
    // The former recursive prefix-extension kept all k prefixes alive (O(k^2) chars, ~20 GB here).
    [Test]
    public void GenerateAllKmers_SingleLetterHugeK_ProducesHomopolymerInLinearSpace()
    {
        const int k = 100_000;
        var result = KmerAnalyzer.GenerateAllKmers(k, "A").ToList();

        Assert.That(result, Has.Count.EqualTo(1));
        Assert.That(result[0], Is.EqualTo(new string('A', k)));
    }

    // C6 — validation is eager (not deferred to enumeration).
    [Test]
    public void GenerateAllKmers_InvalidArguments_ThrowEagerlyWithoutEnumeration()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.GenerateAllKmers(0));
            Assert.Throws<ArgumentException>(() => KmerAnalyzer.GenerateAllKmers(2, ""));
        });
    }

    #endregion
}
