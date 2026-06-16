// SEQ-ENTROPY-PROFILE-001 — Shannon Entropy Profile (sliding-window Shannon entropy)
// Evidence: docs/Evidence/SEQ-ENTROPY-PROFILE-001-Evidence.md
// TestSpec: tests/TestSpecs/SEQ-ENTROPY-PROFILE-001.md
// Source: Shannon C. E. (1948). A Mathematical Theory of Communication.
//         Bell System Technical Journal 27(3):379-423.
//         Wikipedia, Entropy (information theory) (citing Shannon 1948).
//         Entropy-Based Biological Sequence Study, IntechOpen (Eq. 3).

using System.Linq;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class SequenceStatistics_CalculateEntropyProfile_Tests
{
    // Expected entropy values are derived by hand from H = -Σ pᵢ log₂ pᵢ (base 2 → bits),
    // computed independently of the implementation, per the Evidence datasets.
    private const double Tolerance = 1e-10;

    // -(3/4·log2(3/4) + 1/4·log2(1/4)) — the skewed 3:1 window value used across cases.
    private const double H_3to1 = 0.8112781244591328;

    #region CalculateShannonEntropy (per-window kernel)

    // M1 — uniform 4-symbol window: H = log2(4) = 2 bits (max for DNA alphabet).
    // Evidence: Wikipedia max=log2(n); IntechOpen 2-bit DNA maximum.
    [Test]
    public void CalculateShannonEntropy_UniformFourSymbols_ReturnsTwoBits()
    {
        double h = SequenceStatistics.CalculateShannonEntropy("ATGC");

        Assert.That(h, Is.EqualTo(2.0).Within(Tolerance),
            "Equal A,T,G,C → H = log2(4) = 2 bits, the DNA maximum (INV-04/INV-02)");
    }

    // M2 — two equally-frequent symbols: H = 1 bit.
    // Evidence: Shannon 1948 formula; -2·(1/2·log2(1/2)) = 1.
    [Test]
    public void CalculateShannonEntropy_TwoEqualSymbols_ReturnsOneBit()
    {
        double h = SequenceStatistics.CalculateShannonEntropy("AATT");

        Assert.That(h, Is.EqualTo(1.0).Within(Tolerance),
            "A=2,T=2 → H = -2·(1/2·log2 1/2) = 1 bit");
    }

    // M3 — skewed 3:1 window: H = 0.8112781244591328 (not the trivial 0, 1, or 2).
    // Evidence: Shannon 1948 formula -(3/4·log2 3/4 + 1/4·log2 1/4).
    [Test]
    public void CalculateShannonEntropy_Skewed3To1_ReturnsExactValue()
    {
        double h = SequenceStatistics.CalculateShannonEntropy("AAAT");

        Assert.That(h, Is.EqualTo(H_3to1).Within(Tolerance),
            "A=3,T=1 → H = -(3/4 log2 3/4 + 1/4 log2 1/4) = 0.8112781244591328");
    }

    // M4 — three-symbol 3:2:1 distribution: H = 1.4591479170272448.
    // Evidence: Shannon 1948 formula -(1/2 log2 1/2 + 1/3 log2 1/3 + 1/6 log2 1/6).
    [Test]
    public void CalculateShannonEntropy_ThreeSymbols3To2To1_ReturnsExactValue()
    {
        double h = SequenceStatistics.CalculateShannonEntropy("AAATTC");

        Assert.That(h, Is.EqualTo(1.4591479170272448).Within(Tolerance),
            "A=3,T=2,C=1 → H = -(1/2 log2 1/2 + 1/3 log2 1/3 + 1/6 log2 1/6) = 1.4591479170272448");
    }

    // M5 — homopolymer window: H = 0 (zero-probability convention).
    // Evidence: Shannon 1948 / Wikipedia; single symbol p=1, log2 1 = 0 (INV-03).
    [Test]
    public void CalculateShannonEntropy_Homopolymer_ReturnsZero()
    {
        double h = SequenceStatistics.CalculateShannonEntropy("AAAA");

        Assert.That(h, Is.EqualTo(0.0).Within(Tolerance),
            "Single symbol p=1 → H = -1·log2 1 = 0 bits (INV-03)");
    }

    #endregion

    #region CalculateEntropyProfile (sliding window)

    // M6 — full profile, window 4, step 1 over "AAATGC": three windows AAAT, AATG, ATGC.
    // Evidence: IntechOpen sliding window + Shannon formula per window.
    [Test]
    public void CalculateEntropyProfile_Window4Step1_ReturnsExactProfile()
    {
        var profile = SequenceStatistics.CalculateEntropyProfile("AAATGC", 4, 1).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(profile.Length, Is.EqualTo(3),
                "n=6,w=4,step=1 → ⌊(6-4)/1⌋+1 = 3 windows (INV-05)");
            Assert.That(profile[0], Is.EqualTo(H_3to1).Within(Tolerance),
                "Window AAAT → 0.8112781244591328");
            Assert.That(profile[1], Is.EqualTo(1.5).Within(Tolerance),
                "Window AATG (A=2,T=1,G=1) → -(1/2 log2 1/2 + 2·1/4 log2 1/4) = 1.5");
            Assert.That(profile[2], Is.EqualTo(2.0).Within(Tolerance),
                "Window ATGC → log2 4 = 2.0");
        });
    }

    // M7 — full profile, window 4, step 2 over "AAATGCAA": windows AAAT(0), ATGC(2), GCAA(4).
    // Evidence: IntechOpen sliding window with step W>1 + Shannon formula per window.
    [Test]
    public void CalculateEntropyProfile_Window4Step2_ReturnsExactProfile()
    {
        var profile = SequenceStatistics.CalculateEntropyProfile("AAATGCAA", 4, 2).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(profile.Length, Is.EqualTo(3),
                "n=8,w=4,step=2 → offsets 0,2,4 = 3 windows (INV-05)");
            Assert.That(profile[0], Is.EqualTo(H_3to1).Within(Tolerance),
                "Window AAAT → 0.8112781244591328");
            Assert.That(profile[1], Is.EqualTo(2.0).Within(Tolerance),
                "Window ATGC → log2 4 = 2.0");
            Assert.That(profile[2], Is.EqualTo(1.5).Within(Tolerance),
                "Window GCAA (G=1,C=1,A=2) → 1.5");
        });
    }

    // M8 — window count obeys INV-05 for several configurations.
    // Evidence: IntechOpen sliding window of width W (INV-05).
    [Test]
    public void CalculateEntropyProfile_WindowCount_MatchesInvariant()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceStatistics.CalculateEntropyProfile("AAATGCAA", 4, 1).Count(),
                Is.EqualTo(5), "n=8,w=4,step=1 → 5 windows (INV-05)");
            Assert.That(SequenceStatistics.CalculateEntropyProfile("AAATGCAA", 4, 2).Count(),
                Is.EqualTo(3), "n=8,w=4,step=2 → 3 windows (INV-05)");
            Assert.That(SequenceStatistics.CalculateEntropyProfile("AAATGCAA", 4, 3).Count(),
                Is.EqualTo(2), "n=8,w=4,step=3 → offsets 0,3 = 2 windows (INV-05)");
        });
    }

    #endregion

    #region Edge cases and invariants

    // S1 — windowSize greater than sequence length yields an empty profile.
    // Evidence: INV-05 (no full window exists when W > n).
    [Test]
    public void CalculateEntropyProfile_WindowLargerThanSequence_ReturnsEmpty()
    {
        var profile = SequenceStatistics.CalculateEntropyProfile("AAA", 4, 1).ToArray();

        Assert.That(profile, Is.Empty,
            "W=4 > n=3 → no full window → empty profile (INV-05)");
    }

    // S2 — windowSize equal to length yields exactly one window equal to the whole-sequence entropy.
    // Evidence: INV-05; AATT → 1 bit.
    [Test]
    public void CalculateEntropyProfile_WindowEqualsLength_ReturnsSingleValue()
    {
        var profile = SequenceStatistics.CalculateEntropyProfile("AATT", 4, 1).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(profile.Length, Is.EqualTo(1),
                "W == n → exactly one window (INV-05)");
            Assert.That(profile[0], Is.EqualTo(1.0).Within(Tolerance),
                "AATT → 1 bit");
        });
    }

    // S3 — null and empty input yield empty profiles.
    // Evidence: guarded input (§3.3).
    [Test]
    public void CalculateEntropyProfile_NullOrEmpty_ReturnsEmpty()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceStatistics.CalculateEntropyProfile(null!, 4, 1), Is.Empty,
                "null sequence → empty profile");
            Assert.That(SequenceStatistics.CalculateEntropyProfile("", 4, 1), Is.Empty,
                "empty sequence → empty profile");
        });
    }

    // S4 — INV-02: no DNA-window entropy exceeds log2(4) = 2 bits.
    // Evidence: Wikipedia max=log2(n); IntechOpen 2-bit DNA maximum.
    [Test]
    public void CalculateEntropyProfile_DnaWindows_NeverExceedTwoBits()
    {
        var profile = SequenceStatistics.CalculateEntropyProfile("ATGCATGCATGCGGGGAAAACCCC", 6, 1).ToArray();

        Assert.That(profile, Is.Not.Empty, "profile must contain windows");
        Assert.That(profile.All(h => h <= 2.0 + Tolerance), Is.True,
            "every DNA-window entropy ≤ log2(4) = 2 bits (INV-02)");
    }

    // S5 — INV-01: entropy is never negative.
    // Evidence: Wikipedia (entropy non-negative).
    [Test]
    public void CalculateEntropyProfile_AllWindows_NonNegative()
    {
        var profile = SequenceStatistics.CalculateEntropyProfile("ATGCATGCATGCGGGGAAAACCCC", 6, 1).ToArray();

        Assert.That(profile.All(h => h >= -Tolerance), Is.True,
            "every window entropy ≥ 0 (INV-01)");
    }

    // C1 — case-insensitivity: lowercase input produces the same profile as uppercase.
    // Evidence: implementation case-folds before counting (§3.3).
    [Test]
    public void CalculateEntropyProfile_LowercaseInput_MatchesUppercase()
    {
        var upper = SequenceStatistics.CalculateEntropyProfile("AAATGC", 4, 1).ToArray();
        var lower = SequenceStatistics.CalculateEntropyProfile("aaatgc", 4, 1).ToArray();

        Assert.That(lower, Is.EqualTo(upper).Within(Tolerance),
            "case-folded counting → lowercase profile equals uppercase profile");
    }

    // C2 — degenerate N is counted as its own symbol (no special-casing).
    // ACGN has 4 distinct equally-frequent letters → H = log2(4) = 2.0.
    // Expected value derived from H = -Σ pᵢ log₂ pᵢ (Shannon 1948), independent of code;
    // documents the §3.3 behavior "degenerate/N symbols are counted as their own symbol".
    [Test]
    public void CalculateShannonEntropy_NCountedAsOwnSymbol_ReturnsLog2OfDistinctCount()
    {
        double h = SequenceStatistics.CalculateShannonEntropy("ACGN");

        Assert.That(h, Is.EqualTo(2.0).Within(Tolerance),
            "A,C,G,N each once → 4 distinct equal symbols → H = log2(4) = 2 bits");
    }

    // C3 — no T↔U normalization: T and U are distinct symbols.
    // ATUG has 4 distinct equally-frequent letters → H = log2(4) = 2.0 (not 3-symbol log2(3)).
    // Expected value from the Shannon formula; documents §3.3 "no T↔U normalization".
    [Test]
    public void CalculateShannonEntropy_TAndUAreDistinctSymbols_NoNormalization()
    {
        double h = SequenceStatistics.CalculateShannonEntropy("ATUG");

        Assert.That(h, Is.EqualTo(2.0).Within(Tolerance),
            "A,T,U,G distinct (no T↔U merge) → 4 equal symbols → H = log2(4) = 2 bits");
    }

    // C4 — non-letter characters are ignored (only char.IsLetter symbols counted).
    // "A-C-G-T" contributes 4 distinct equal letters → H = log2(4) = 2.0, dashes ignored.
    // Expected value from the Shannon formula; documents §3.3 "non-letters ignored".
    [Test]
    public void CalculateShannonEntropy_NonLetterCharacters_AreIgnored()
    {
        double h = SequenceStatistics.CalculateShannonEntropy("A-C-G-T");

        Assert.That(h, Is.EqualTo(2.0).Within(Tolerance),
            "dashes ignored; A,C,G,T each once → H = log2(4) = 2 bits");
    }

    // INV-2 (general k>4) — entropy is alphabet-sensitive; a protein-style window with
    // 8 distinct equal residues yields log2(8) = 3 bits, exceeding the 2-bit DNA ceiling.
    // Expected value from H = log₂ k at the uniform distribution (Wikipedia max=log₂ n);
    // documents §6.2 "protein windows can exceed 2 bits since k > 4".
    [Test]
    public void CalculateShannonEntropy_EightDistinctSymbols_ExceedsTwoBits()
    {
        double h = SequenceStatistics.CalculateShannonEntropy("ACDEFGHI");

        Assert.That(h, Is.EqualTo(3.0).Within(Tolerance),
            "8 distinct equal residues → H = log2(8) = 3 bits (> 2-bit DNA max, INV-02 with k=8)");
    }

    #endregion
}
