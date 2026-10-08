namespace Seqeron.Genomics.Tests.Metamorphic;

/// <summary>
/// Heavy-tier metamorphic relations for the behaviours introduced by review batch B04
/// (docs/Validation/review-2026-09/B04.md): LZ76 (F1), DUST divisor (F2), SDUST symmetry (F3),
/// ACGT-only palindromes (F12), TRF approximate repeats (F14–F18), entropy kernel (F19) and the
/// linguistic-complexity alphabet (F20). Fixed-seed random inputs; every relation is exact.
///
/// Test Units: SEQ-COMPLEX-COMPRESS-001, SEQ-COMPLEX-DUST-001, SEQ-COMPLEX-001, SEQ-ENTROPY-001,
/// REP-PALIN-001, REP-APPROX-001.
/// </summary>
[TestFixture]
[Category("Metamorphic")]
[Category("Complexity")]
public class B04ComplexityRepeatsMetamorphicTests
{
    private const int Seed = 20260930;

    private static string Random(Random rng, string alphabet, int length) =>
        new(Enumerable.Range(0, length).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());

    // Repeat-rich DNA: a random 1–4 letter sub-alphabet, sometimes with a planted short tandem block.
    private static string RepeatRich(Random rng, int length)
    {
        string[] alphabets = ["A", "AC", "AT", "CG", "ACG", "ACGT"];
        string s = Random(rng, alphabets[rng.Next(alphabets.Length)], length);
        if (length > 10 && rng.Next(2) == 0)
        {
            string unit = Random(rng, "ACGT", rng.Next(1, 4));
            string block = string.Concat(Enumerable.Repeat(unit, rng.Next(3, 8)));
            int at = rng.Next(0, length);
            s = (s[..at] + block + s[at..])[..length];
        }
        return s;
    }

    private static string Reverse(string s) => new(s.Reverse().ToArray());

    private static string ReverseComplementAcgt(string s) =>
        new(s.Reverse().Select(c => c switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', _ => 'C' }).ToArray());

    /// <summary>
    /// MR (F1): LZ76 complexity depends only on the equality pattern of symbols, so any bijective relabelling
    /// of the alphabet (here A→C→G→T→A and the complement A↔T, C↔G) leaves c(S) unchanged.
    /// </summary>
    [Test]
    public void LempelZiv_BijectiveRelabelling_PreservesComplexity()
    {
        var rng = new Random(Seed);
        for (int trial = 0; trial < 300; trial++)
        {
            string s = RepeatRich(rng, rng.Next(0, 120));
            string rotated = new(s.Select(c => c switch { 'A' => 'C', 'C' => 'G', 'G' => 'T', _ => 'A' }).ToArray());
            string complement = new(s.Select(c => c switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', _ => 'C' }).ToArray());
            int c0 = SequenceComplexity.CalculateLempelZivComplexity(s);
            Assert.That(SequenceComplexity.CalculateLempelZivComplexity(rotated), Is.EqualTo(c0), $"rotated '{s}'");
            Assert.That(SequenceComplexity.CalculateLempelZivComplexity(complement), Is.EqualTo(c0), $"complement '{s}'");
        }
    }

    /// <summary>
    /// MR (F2): reversing the sequence reverses every triplet, which permutes triplet identities but keeps the
    /// multiset of counts and ℓ, so the DUST score Σc(c−1)/2 / (ℓ−1) is unchanged; the same holds for the
    /// reverse complement.
    /// </summary>
    [Test]
    public void Dust_ReverseAndReverseComplement_PreserveScore()
    {
        var rng = new Random(Seed + 1);
        for (int trial = 0; trial < 300; trial++)
        {
            string s = RepeatRich(rng, rng.Next(0, 150));
            double d = SequenceComplexity.CalculateDustScore(s);
            Assert.That(SequenceComplexity.CalculateDustScore(Reverse(s)), Is.EqualTo(d), $"reverse '{s}'");
            Assert.That(SequenceComplexity.CalculateDustScore(ReverseComplementAcgt(s)), Is.EqualTo(d), $"revcomp '{s}'");
        }
    }

    /// <summary>
    /// MR (F3): SDUST is symmetric (Morgulis et al. 2006): masking the reverse (or the reverse complement) of S
    /// masks exactly the mirrored positions of the mask of S.
    /// </summary>
    [Test]
    public void Sdust_ReverseAndReverseComplement_MirrorMask()
    {
        var rng = new Random(Seed + 2);
        for (int trial = 0; trial < 400; trial++)
        {
            string s = RepeatRich(rng, rng.Next(0, 200));
            int w = rng.Next(3, 70);
            double t = new[] { 0.5, 1.0, 2.0, 2.5, 4.0 }[rng.Next(5)];
            string mask = new(SequenceComplexity.MaskLowComplexity(new DnaSequence(s), w, t, 'x')
                .Select(c => c == 'x' ? '1' : '0').ToArray());
            string mirrored = Reverse(mask);
            foreach (string transformed in new[] { Reverse(s), ReverseComplementAcgt(s) })
            {
                string m2 = new(SequenceComplexity.MaskLowComplexity(new DnaSequence(transformed), w, t, 'x')
                    .Select(c => c == 'x' ? '1' : '0').ToArray());
                Assert.That(m2, Is.EqualTo(mirrored), $"'{s}' W={w} T={t}");
            }
        }
    }

    /// <summary>
    /// MR (F12): a reverse-complement palindrome of S at (p, L) is a palindrome of revcomp(S) at (n − p − L, L)
    /// with the same text, so the palindrome set of revcomp(S) is the mirror image of that of S.
    /// </summary>
    [Test]
    public void Palindromes_ReverseComplementInput_MirrorsPalindromeSet()
    {
        var rng = new Random(Seed + 3);
        for (int trial = 0; trial < 300; trial++)
        {
            string s = Random(rng, rng.Next(2) == 0 ? "AT" : "ACGT", rng.Next(0, 120));
            int n = s.Length;
            var fwd = RepeatFinder.FindPalindromes(s, 4, 16)
                .Select(p => (Pos: n - p.Position - p.Length, p.Length, p.Sequence))
                .OrderBy(p => p.Pos).ThenBy(p => p.Length).ToList();
            var rc = RepeatFinder.FindPalindromes(ReverseComplementAcgt(s), 4, 16)
                .Select(p => (Pos: p.Position, p.Length, p.Sequence)).ToList();
            Assert.That(rc, Is.EqualTo(fwd), $"'{s}'");
        }
    }

    /// <summary>
    /// MR (F14–F18): the string overload is case-insensitive — lowercasing the input (or any part of it) yields
    /// the identical TRF result list; and inserting a non-ACGT symbol run that cannot match (N) far from a repeat
    /// never creates a result made only of N (TRF scores A/C/G/T only).
    /// </summary>
    [Test]
    public void ApproximateRepeats_Lowercasing_IdenticalResults()
    {
        var rng = new Random(Seed + 4);
        string[] units = ["CAG", "AC", "GATA", "TTAGGG", "A", "ACGTC"];
        for (int trial = 0; trial < 60; trial++)
        {
            string unit = units[rng.Next(units.Length)];
            var core = string.Concat(Enumerable.Repeat(unit, rng.Next(4, 12))).ToCharArray();
            for (int m = rng.Next(0, 3); m > 0; m--) core[rng.Next(core.Length)] = "ACGT"[rng.Next(4)];
            string s = Random(rng, "ACGT", rng.Next(0, 20)) + new string(core) + Random(rng, "ACGT", rng.Next(0, 20));

            var upper = RepeatFinder.FindApproximateTandemRepeats(s, 1, 12, 20).ToList();
            var lower = RepeatFinder.FindApproximateTandemRepeats(s.ToLowerInvariant(), 1, 12, 20).ToList();
            var mixed = RepeatFinder.FindApproximateTandemRepeats(
                new string(s.Select((c, i) => i % 2 == 0 ? char.ToLowerInvariant(c) : c).ToArray()), 1, 12, 20).ToList();
            Assert.That(lower, Is.EqualTo(upper), $"lowercase '{s}'");
            Assert.That(mixed, Is.EqualTo(upper), $"mixed case '{s}'");

            string withN = s + new string('N', 30);
            var nResults = RepeatFinder.FindApproximateTandemRepeats(withN, 1, 12, 20).ToList();
            Assert.That(nResults.Any(r => r.Start >= s.Length), Is.False, $"N-only repeat reported for '{withN}'");
        }
    }

    /// <summary>
    /// MR (F19): Shannon entropy depends only on base counts, so every permutation (Fisher–Yates shuffle) of S
    /// gives the bit-identical value; and the RNA counterpart (T→U) has the same entropy.
    /// </summary>
    [Test]
    public void ShannonEntropy_ShuffleAndRnaCounterpart_Identical()
    {
        var rng = new Random(Seed + 5);
        for (int trial = 0; trial < 300; trial++)
        {
            string s = RepeatRich(rng, rng.Next(0, 100));
            double h = SequenceComplexity.CalculateShannonEntropy(s);
            var a = s.ToCharArray();
            for (int i = a.Length - 1; i > 0; i--)
            {
                int j = rng.Next(i + 1);
                (a[i], a[j]) = (a[j], a[i]);
            }
            Assert.That(SequenceComplexity.CalculateShannonEntropy(new string(a)), Is.EqualTo(h), $"shuffle '{s}'");
            Assert.That(SequenceComplexity.CalculateShannonEntropy(s.Replace('T', 'U')), Is.EqualTo(h), $"RNA '{s}'");
        }
    }

    /// <summary>
    /// MR (F20): the RNA counterpart (T→U) of a DNA string has the same linguistic complexity (alphabet
    /// {A,C,G,U}, a = 4, identical V_i) on both the hash (m ≤ 12) and suffix-tree (m &gt; 12) paths; lowercasing
    /// is also invariant.
    /// </summary>
    [Test]
    public void LinguisticComplexity_RnaCounterpartAndCase_Identical()
    {
        var rng = new Random(Seed + 6);
        for (int trial = 0; trial < 200; trial++)
        {
            string s = RepeatRich(rng, rng.Next(1, 80));
            foreach (int m in new[] { 1, 4, 12, 13, int.MaxValue })
            {
                double lc = SequenceComplexity.CalculateLinguisticComplexity(s, m);
                Assert.That(SequenceComplexity.CalculateLinguisticComplexity(s.Replace('T', 'U'), m), Is.EqualTo(lc), $"RNA '{s}' m={m}");
                Assert.That(SequenceComplexity.CalculateLinguisticComplexity(s.ToLowerInvariant(), m), Is.EqualTo(lc), $"lower '{s}' m={m}");
            }
        }
    }
}
