using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// B09 WP78 R14-18 (F360): the bit-parallel engine of <see cref="ApproximateMatcher.FindWithMismatches(string, string, int, bool, CancellationToken)"/> (2 bits per base, XOR + popcount for A/C/G/T patterns of up to 32 symbols;
/// one-hot nibbles and a shared-bit test for IUPAC patterns of up to 16 symbols) against the character loop that it replaces: the same windows, the same mismatch counts and positions, the same order, on random inputs
/// (lower case, N / IUPAC symbols and U in the sequence, literal and IUPAC patterns, k = 0 ... 6).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("PAT-APPROX-001")]
public class ApproximateMatcher_PackedHamming_Tests
{
    // the previous implementation, verbatim in its semantics: every window is compared symbol by symbol (literal equality or IupacHelper.MatchesIupac)
    private static List<(int Position, string Window, int Distance, int[] Positions)> Reference(string sequence, string pattern, int k, bool iupac)
    {
        var result = new List<(int, string, int, int[])>();
        if (string.IsNullOrEmpty(sequence) || string.IsNullOrEmpty(pattern))
            return result;
        string seq = sequence.ToUpperInvariant(), pat = pattern.ToUpperInvariant();
        for (int i = 0; i + pat.Length <= seq.Length; i++)
        {
            var mismatches = new List<int>();
            for (int j = 0; j < pat.Length; j++)
            {
                bool miss = iupac ? !IupacHelper.MatchesIupac(seq[i + j], pat[j]) : seq[i + j] != pat[j];
                if (miss)
                    mismatches.Add(j);
            }

            if (mismatches.Count <= k)
                result.Add((i, seq.Substring(i, pat.Length), mismatches.Count, mismatches.ToArray()));
        }

        return result;
    }

    private static void AssertEqual(string sequence, string pattern, int k, bool iupac)
    {
        var expected = Reference(sequence, pattern, k, iupac);
        var actual = ApproximateMatcher.FindWithMismatches(sequence, pattern, k, iupac, CancellationToken.None)
            .Select(m => (m.Position, m.MatchedSequence, m.Distance, m.MismatchPositions.ToArray())).ToList();
        Assert.That(actual.Count, Is.EqualTo(expected.Count), $"seq={sequence} pat={pattern} k={k} iupac={iupac}");
        for (int n = 0; n < expected.Count; n++)
        {
            Assert.That(actual[n].Position, Is.EqualTo(expected[n].Position), $"hit {n}: seq={sequence} pat={pattern} k={k} iupac={iupac}");
            Assert.That(actual[n].MatchedSequence, Is.EqualTo(expected[n].Window));
            Assert.That(actual[n].Distance, Is.EqualTo(expected[n].Distance));
            Assert.That(actual[n].Item4, Is.EqualTo(expected[n].Positions));
        }
    }

    private static string Random(Random rng, string alphabet, int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = alphabet[rng.Next(alphabet.Length)];
        return new string(chars);
    }

    [Test]
    public void LiteralPatternOfBases_EqualsTheCharacterLoop_OnRandomInputs()
    {
        var rng = new Random(7801);
        for (int round = 0; round < 1500; round++)
        {
            // the sequence side: lower case, N, IUPAC letters, U, a gap: none of them is one of A/C/G/T, so every window that holds one mismatches there
            string alphabet = rng.Next(3) == 0 ? "ACGT" : rng.Next(2) == 0 ? "ACGTacgtN" : "ACGTRYKMSWNUn-";
            string seq = Random(rng, alphabet, rng.Next(0, 300));
            int m = rng.Next(1, 36);
            string pat = rng.Next(4) == 0 && seq.Length >= m ? seq.Substring(rng.Next(0, seq.Length - m + 1), m).ToUpperInvariant().Replace('N', 'A').Replace('U', 'T').Replace('R', 'G').Replace('Y', 'C').Replace('K', 'G').Replace('M', 'A').Replace('S', 'C').Replace('W', 'A').Replace('-', 'T') : Random(rng, "ACGT", m);
            AssertEqual(seq, pat, rng.Next(0, 7), iupac: false);
            AssertEqual(seq, pat, rng.Next(0, 7), iupac: true);
        }
    }

    [Test]
    public void LiteralPatternWithOtherSymbols_KeepsTheCharacterLoop_SoNMatchesN()
    {
        var rng = new Random(7802);
        for (int round = 0; round < 600; round++)
        {
            string seq = Random(rng, "ACGTNn", rng.Next(0, 200));
            string pat = Random(rng, "ACGTN", rng.Next(1, 14));
            AssertEqual(seq, pat, rng.Next(0, 4), iupac: false);
        }
    }

    [Test]
    public void IupacPattern_EqualsTheCharacterLoop_OnRandomInputs()
    {
        var rng = new Random(7803);
        for (int round = 0; round < 2000; round++)
        {
            string seq = Random(rng, rng.Next(2) == 0 ? "ACGT" : "ACGTacgtNRYU-", rng.Next(0, 300));
            string pat = Random(rng, "ACGTNRYSWKMBDHV", rng.Next(1, 22));   // 17 ... 21 symbols use the character loop
            AssertEqual(seq, pat, rng.Next(0, 6), iupac: true);
        }
    }

    [Test]
    public void BoundaryLengths_OfTheWord_AreExact()
    {
        var rng = new Random(7804);
        foreach (int m in new[] { 1, 2, 15, 16, 17, 31, 32, 33, 40 })
        {
            string seq = Random(rng, "ACGT", 400);
            string exact = seq.Substring(120, m);
            AssertEqual(seq, exact, 0, iupac: false);
            AssertEqual(seq, exact, 3, iupac: false);
            AssertEqual(seq, exact, 3, iupac: true);
            if (m <= 24)
            {
                string iupacPat = string.Concat(exact.Select((c, i) => i % 4 == 1 ? 'N' : i % 5 == 2 ? 'R' : c));
                AssertEqual(seq, iupacPat, 2, iupac: true);
            }
        }
    }

    [Test]
    public void Cancellation_IsObservedInThePackedScan()
    {
        string seq = new string('A', 5000);
        using var cts = new CancellationTokenSource();
        cts.Cancel();
        Assert.Throws<OperationCanceledException>(() => ApproximateMatcher.FindWithMismatches(seq, "AAAA", 1, false, cts.Token).ToList());
        Assert.Throws<OperationCanceledException>(() => ApproximateMatcher.FindWithMismatches(seq, "AANA", 1, true, cts.Token).ToList());
    }

    [Test]
    public void KnownWindow_IsReportedWithItsMismatchPositions()
    {
        // ACGTACGT against ACGAACGT: one mismatch at pattern position 3
        var hits = ApproximateMatcher.FindWithMismatches("ACGTACGT", "ACGAACGT", 1).ToList();
        Assert.That(hits, Has.Count.EqualTo(1));
        Assert.That(hits[0].Position, Is.EqualTo(0));
        Assert.That(hits[0].Distance, Is.EqualTo(1));
        Assert.That(hits[0].MismatchPositions, Is.EqualTo(new[] { 3 }));
        Assert.That(hits[0].MatchedSequence, Is.EqualTo("ACGTACGT"));
    }
}
