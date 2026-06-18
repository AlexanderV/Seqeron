using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for miRNA analysis: pre-miRNA hairpins and seed sequence analysis.
///
/// Test Units: MIRNA-PRECURSOR-001, MIRNA-SEED-001
/// MIRNA-SEED-001 property tests removed — consolidated into canonical MiRnaAnalyzer_SeedAnalysis_Tests.cs
/// (3 duplicates of M-003/M-007/M-009; 1 weak: IsSubstringOfMiRna can't distinguish extraction position)
/// MIRNA-TARGET-001 property tests removed — consolidated into canonical MiRnaAnalyzer_TargetPrediction_Tests.cs
/// (FindTargetSites_Positions_WithinBounds → S-002; FindTargetSites_Score_InRange → M-009)
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Annotation")]
public class MiRnaProperties
{
    // -- MIRNA-PRECURSOR-001 --

    /// <summary>
    /// Pre-miRNA hairpin has structure with balanced parentheses.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindPreMiRnaHairpins_Structure_BalancedParentheses()
    {
        // A sequence long enough to contain potential hairpins
        string sequence = string.Concat(Enumerable.Repeat("GCGCUUUUGCGC", 20));
        var hairpins = MiRnaAnalyzer.FindPreMiRnaHairpins(sequence, minHairpinLength: 30).ToList();

        foreach (var hp in hairpins)
        {
            int opens = hp.Structure.Count(c => c == '(');
            int closes = hp.Structure.Count(c => c == ')');
            Assert.That(opens, Is.EqualTo(closes),
                $"Unbalanced structure: {opens} opens vs {closes} closes");
        }
    }

    /// <summary>
    /// Pre-miRNA start/end are within sequence bounds.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindPreMiRnaHairpins_Bounds_WithinSequence()
    {
        string sequence = string.Concat(Enumerable.Repeat("GCGCUUUUGCGC", 20));
        var hairpins = MiRnaAnalyzer.FindPreMiRnaHairpins(sequence, minHairpinLength: 30).ToList();

        foreach (var hp in hairpins)
        {
            Assert.That(hp.Start, Is.GreaterThanOrEqualTo(0));
            Assert.That(hp.End, Is.LessThanOrEqualTo(sequence.Length));
            Assert.That(hp.End, Is.GreaterThan(hp.Start));
        }
    }

    /// <summary>
    /// Pre-miRNA structure length equals sequence length.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindPreMiRnaHairpins_StructureLength_EqualsSequenceLength()
    {
        string sequence = string.Concat(Enumerable.Repeat("GCGCUUUUGCGC", 20));
        var hairpins = MiRnaAnalyzer.FindPreMiRnaHairpins(sequence, minHairpinLength: 30).ToList();

        foreach (var hp in hairpins)
            Assert.That(hp.Structure.Length, Is.EqualTo(hp.Sequence.Length),
                $"Structure length {hp.Structure.Length} ≠ sequence length {hp.Sequence.Length}");
    }

    #region MIRNA-SEED-001

    // ---------------------------------------------------------------------
    // MIRNA-SEED-001 — Seed Sequence Analysis
    //
    // Reference: docs/algorithms/MiRNA/Seed_Sequence_Analysis.md
    //   §2.2 Core model: seed(m) = m_2 m_3 ... m_8 (1-based positions 2..8),
    //        i.e. zero-based indices 1..7 — EXACTLY 7 nt, the canonical
    //        TargetScan seed length used throughout this repository.
    //   §2.4 Invariants INV-01..INV-05.
    //   §5.2 Two distinct normalization paths:
    //        - GetSeedSequence uppercases ONLY (keeps DNA 'T'),
    //        - CreateMiRna normalizes 'T'->'U' AND uppercases, stores
    //          SeedStart=1, SeedEnd=7.
    //
    // NOTE on the checklist phrasing "seed len = 6-8": that loose range is
    // NOT the contract of THIS implementation. The implementation produces a
    // seed that is EITHER the empty string (input shorter than 8 nt) OR a
    // sequence of length EXACTLY 7 (the canonical TargetScan 7-mer seed,
    // positions 2..8). The tests below assert that exact 7-or-empty contract.
    //
    // All expected values are derived INDEPENDENTLY from the doc, never read
    // back from the implementation under test.
    // ---------------------------------------------------------------------

    private const int CanonicalSeedLength = 7;

    /// <summary>RNA/DNA alphabet over which random miRNA sequences are drawn.</summary>
    private static readonly char[] NucleotideAlphabet = { 'A', 'C', 'G', 'U', 'T', 'a', 'c', 'g', 'u', 't' };

    /// <summary>
    /// Generates random nucleotide strings (mixed case, mixed DNA/RNA alphabet)
    /// of length in <paramref name="minLen"/>..<paramref name="maxLen"/> inclusive.
    /// </summary>
    private static Gen<string> SequenceGen(int minLen, int maxLen) =>
        from len in Gen.Choose(minLen, maxLen)
        from chars in Gen.Elements(NucleotideAlphabet).ArrayOf(len)
        select new string(chars);

    /// <summary>Sequences guaranteed long enough to yield a canonical 7-nt seed (length ≥ 8).</summary>
    private static Arbitrary<string> SeedableSequenceArbitrary() =>
        SequenceGen(8, 40).ToArbitrary();

    /// <summary>Sequences of any length, including those shorter than 8 nt and the empty string.</summary>
    private static Arbitrary<string> AnyLengthSequenceArbitrary() =>
        SequenceGen(0, 40).ToArbitrary();

    /// <summary>Short sequences (length 0..7) that can never produce a seed.</summary>
    private static Arbitrary<string> TooShortSequenceArbitrary() =>
        SequenceGen(0, 7).ToArbitrary();

    /// <summary>
    /// Generates a pair of seedable sequences that share an IDENTICAL seed window
    /// (zero-based indices 1..7) under <c>CreateMiRna</c> normalization, by reusing
    /// the same 7-nt RNA core while letting the surrounding bases vary freely.
    /// </summary>
    private static Arbitrary<(string a, string b)> SameSeedPairArbitrary() =>
        (from core in Gen.Elements('A', 'C', 'G', 'U').ArrayOf(CanonicalSeedLength)
         from lead1 in Gen.Elements('A', 'C', 'G', 'U')
         from lead2 in Gen.Elements('A', 'C', 'G', 'U')
         from tail1 in Gen.Elements('A', 'C', 'G', 'U').ArrayOf()
         from tail2 in Gen.Elements('A', 'C', 'G', 'U').ArrayOf()
         let seed = new string(core)
         select (lead1 + seed + new string(tail1), lead2 + seed + new string(tail2)))
        .ToArbitrary();

    /// <summary>
    /// Generates a pair of seedable sequences whose seed windows differ in EXACTLY
    /// one position. The two cores are identical except at a chosen index, where the
    /// second core uses the cyclic-next base — guaranteeing a single mismatch.
    /// </summary>
    private static Arbitrary<(string a, string b)> OneMismatchSeedPairArbitrary() =>
        (from core in Gen.Elements('A', 'C', 'G', 'U').ArrayOf(CanonicalSeedLength)
         from idx in Gen.Choose(0, CanonicalSeedLength - 1)
         from lead1 in Gen.Elements('A', 'C', 'G', 'U')
         from lead2 in Gen.Elements('A', 'C', 'G', 'U')
         let core2 = MutateOnePosition(core, idx)
         select (lead1 + new string(core), lead2 + new string(core2)))
        .ToArbitrary();

    /// <summary>Returns a copy of <paramref name="core"/> with the base at <paramref name="idx"/> changed to a different base.</summary>
    private static char[] MutateOnePosition(char[] core, int idx)
    {
        const string bases = "ACGU";
        char[] mutated = (char[])core.Clone();
        char original = mutated[idx];
        int next = (bases.IndexOf(original) + 1) % bases.Length;
        mutated[idx] = bases[next];
        return mutated;
    }

    /// <summary>
    /// Independent oracle for <c>GetSeedSequence</c>: empty for inputs shorter than 8 nt,
    /// otherwise the uppercased characters at zero-based indices 1..7 (positions 2..8).
    /// DNA 'T' is preserved (uppercase only — no T→U normalization).
    /// </summary>
    private static string ExpectedDirectSeed(string s)
    {
        if (string.IsNullOrEmpty(s) || s.Length < 8)
            return "";

        char[] window = new char[CanonicalSeedLength];
        for (int i = 0; i < CanonicalSeedLength; i++)
            window[i] = char.ToUpperInvariant(s[i + 1]); // zero-based index 1..7
        return new string(window);
    }

    /// <summary>Independent oracle for the normalized sequence stored by <c>CreateMiRna</c>: uppercase with every T→U.</summary>
    private static string ExpectedNormalizedSequence(string s)
    {
        char[] result = new char[s.Length];
        for (int i = 0; i < s.Length; i++)
        {
            char u = char.ToUpperInvariant(s[i]);
            result[i] = u == 'T' ? 'U' : u;
        }
        return new string(result);
    }

    /// <summary>Independent positional match count over two equal-length seeds.</summary>
    private static int CountMatches(string a, string b)
    {
        int matches = 0;
        for (int i = 0; i < a.Length; i++)
            if (a[i] == b[i]) matches++;
        return matches;
    }

    /// <summary>
    /// INV-01 + seed-length truth: GetSeedSequence is empty for inputs shorter than 8 nt
    /// and otherwise equals s.Substring(1,7).ToUpperInvariant() with length EXACTLY 7.
    /// (Asserts the canonical 7-or-empty contract; the checklist "6-8" is not this contract.)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GetSeedSequence_MatchesDirectOracle_And_IsSevenOrEmpty()
    {
        return Prop.ForAll(AnyLengthSequenceArbitrary(), seq =>
        {
            string actual = MiRnaAnalyzer.GetSeedSequence(seq);
            string expected = ExpectedDirectSeed(seq);
            bool lengthOk = actual.Length == 0 || actual.Length == CanonicalSeedLength;
            return (actual == expected && lengthOk)
                .Label($"INV-01: seed='{actual}' expected='{expected}' len={actual.Length} for '{seq}' (len {seq?.Length})");
        });
    }

    /// <summary>
    /// P (positions 2-8 / 5' end): for length ≥ 8 the seed equals, char-by-char, the
    /// uppercased characters at zero-based indices 1..7 (1-based positions 2..8).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GetSeedSequence_IsCharsAtIndicesOneThroughSeven()
    {
        return Prop.ForAll(SeedableSequenceArbitrary(), seq =>
        {
            string seed = MiRnaAnalyzer.GetSeedSequence(seq);
            if (seed.Length != CanonicalSeedLength)
                return false.Label($"P: expected 7-nt seed for length-{seq.Length} '{seq}', got '{seed}'");

            for (int i = 0; i < CanonicalSeedLength; i++)
            {
                char expectedChar = char.ToUpperInvariant(seq[i + 1]);
                if (seed[i] != expectedChar)
                    return false.Label($"P: seed[{i}]='{seed[i]}' ≠ upper(seq[{i + 1}])='{expectedChar}' for '{seq}'");
            }
            return true.ToProperty();
        });
    }

    /// <summary>
    /// INV-02: CreateMiRna(name,seq).SeedSequence == GetSeedSequence(CreateMiRna(name,seq).Sequence),
    /// and the stored Sequence equals seq uppercased with every T→U (independent oracle).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CreateMiRna_NormalizesSequence_And_SeedDerivesFromIt()
    {
        return Prop.ForAll(AnyLengthSequenceArbitrary(), seq =>
        {
            var mirna = MiRnaAnalyzer.CreateMiRna("test", seq);

            bool seqNormalized = mirna.Sequence == ExpectedNormalizedSequence(seq);
            bool seedFromNormalized = mirna.SeedSequence == MiRnaAnalyzer.GetSeedSequence(mirna.Sequence);
            bool seedMatchesOracle = mirna.SeedSequence == ExpectedDirectSeed(ExpectedNormalizedSequence(seq));

            return (seqNormalized && seedFromNormalized && seedMatchesOracle)
                .Label($"INV-02: stored='{mirna.Sequence}' seed='{mirna.SeedSequence}' for input '{seq}'");
        });
    }

    /// <summary>
    /// T-vs-U distinction (§5.2): when the seed window contains DNA 'T', the DIRECT
    /// GetSeedSequence path keeps 'T', while the CreateMiRna path converts those
    /// positions to 'U'. Proves the two normalization paths genuinely differ.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DirectSeed_KeepsT_While_CreateMiRna_ConvertsToU()
    {
        // Build sequences whose seed window (indices 1..7) is guaranteed to contain at least one T.
        var withTinSeedArbitrary =
            (from lead in Gen.Elements('A', 'C', 'G')
             from seedCore in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(CanonicalSeedLength).Where(a => a.Contains('T'))
             from tail in Gen.Elements('A', 'C', 'G', 'T').ArrayOf()
             select lead + new string(seedCore) + new string(tail)).ToArbitrary();

        return Prop.ForAll(withTinSeedArbitrary, seq =>
        {
            string directSeed = MiRnaAnalyzer.GetSeedSequence(seq);
            string createdSeed = MiRnaAnalyzer.CreateMiRna("test", seq).SeedSequence;

            bool directHasT = directSeed.Contains('T');
            bool createdHasNoT = !createdSeed.Contains('T');
            bool createdIsDirectWithUForT = createdSeed == directSeed.Replace('T', 'U');

            return (directHasT && createdHasNoT && createdIsDirectWithUForT)
                .Label($"T-vs-U: direct='{directSeed}' created='{createdSeed}' for '{seq}'");
        });
    }

    /// <summary>
    /// INV-03: CreateMiRna stores SeedStart == 1 and SeedEnd == 7 for ANY input
    /// (zero-based coordinates of the canonical seed window).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CreateMiRna_StoresFixedSeedBounds()
    {
        return Prop.ForAll(AnyLengthSequenceArbitrary(), seq =>
        {
            var mirna = MiRnaAnalyzer.CreateMiRna("test", seq);
            return (mirna.SeedStart == 1 && mirna.SeedEnd == 7)
                .Label($"INV-03: SeedStart={mirna.SeedStart} SeedEnd={mirna.SeedEnd} for '{seq}'");
        });
    }

    /// <summary>
    /// INV-04: CompareSeedRegions(a,b).IsSameFamily is true iff the two stored seeds are
    /// exactly equal — AND false whenever either stored seed is empty.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CompareSeedRegions_IsSameFamily_IffSeedsExactlyEqual()
    {
        // Mix seedable and too-short sequences so empty-seed cases are exercised too.
        var pairArbitrary =
            (from a in SequenceGen(0, 30)
             from b in SequenceGen(0, 30)
             select (a, b)).ToArbitrary();

        return Prop.ForAll(pairArbitrary, pair =>
        {
            var (a, b) = pair;
            var m1 = MiRnaAnalyzer.CreateMiRna("a", a);
            var m2 = MiRnaAnalyzer.CreateMiRna("b", b);
            var cmp = MiRnaAnalyzer.CompareSeedRegions(m1, m2);

            bool eitherEmpty = m1.SeedSequence.Length == 0 || m2.SeedSequence.Length == 0;
            bool expectedFamily = !eitherEmpty && m1.SeedSequence == m2.SeedSequence;

            return (cmp.IsSameFamily == expectedFamily)
                .Label($"INV-04: family={cmp.IsSameFamily} expected={expectedFamily} seeds='{m1.SeedSequence}'/'{m2.SeedSequence}'");
        });
    }

    /// <summary>
    /// INV-05: for two present (non-empty, canonical 7-nt) seeds, Matches + Mismatches == 7,
    /// Matches equals the independently counted positional equality, and Mismatches == 7 - Matches.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CompareSeedRegions_MatchesPlusMismatches_EqualsSeven()
    {
        var pairArbitrary =
            (from a in SequenceGen(8, 30)
             from b in SequenceGen(8, 30)
             select (a, b)).ToArbitrary();

        return Prop.ForAll(pairArbitrary, pair =>
        {
            var (a, b) = pair;
            var m1 = MiRnaAnalyzer.CreateMiRna("a", a);
            var m2 = MiRnaAnalyzer.CreateMiRna("b", b);
            var cmp = MiRnaAnalyzer.CompareSeedRegions(m1, m2);

            int expectedMatches = CountMatches(m1.SeedSequence, m2.SeedSequence);

            bool ok = cmp.Matches == expectedMatches
                      && cmp.Mismatches == CanonicalSeedLength - expectedMatches
                      && cmp.Matches + cmp.Mismatches == CanonicalSeedLength;

            return ok.Label($"INV-05: matches={cmp.Matches} mismatches={cmp.Mismatches} expectedMatches={expectedMatches} seeds='{m1.SeedSequence}'/'{m2.SeedSequence}'");
        });
    }

    /// <summary>
    /// INV-05 (empty branch): when either stored seed is empty, the comparison is
    /// fully zeroed — Matches == 0, Mismatches == 0, IsSameFamily == false.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CompareSeedRegions_EmptySeed_ProducesZeroedResult()
    {
        return Prop.ForAll(TooShortSequenceArbitrary(), SeedableSequenceArbitrary(), (shortSeq, longSeq) =>
        {
            var empty = MiRnaAnalyzer.CreateMiRna("empty", shortSeq);
            var present = MiRnaAnalyzer.CreateMiRna("present", longSeq);

            var cmp1 = MiRnaAnalyzer.CompareSeedRegions(empty, present);
            var cmp2 = MiRnaAnalyzer.CompareSeedRegions(present, empty);

            bool zeroed =
                cmp1.Matches == 0 && cmp1.Mismatches == 0 && !cmp1.IsSameFamily &&
                cmp2.Matches == 0 && cmp2.Mismatches == 0 && !cmp2.IsSameFamily;

            return zeroed.Label($"INV-05 empty: emptySeed='{empty.SeedSequence}' present='{present.SeedSequence}'");
        });
    }

    /// <summary>
    /// SAME-SEED family property: sequences sharing an identical 7-nt seed core are always
    /// reported as the same family with Matches == 7 and Mismatches == 0.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property SameSeedPair_IsFamily_WithSevenMatches()
    {
        return Prop.ForAll(SameSeedPairArbitrary(), pair =>
        {
            var (a, b) = pair;
            var m1 = MiRnaAnalyzer.CreateMiRna("a", a);
            var m2 = MiRnaAnalyzer.CreateMiRna("b", b);
            var cmp = MiRnaAnalyzer.CompareSeedRegions(m1, m2);

            bool ok = m1.SeedSequence == m2.SeedSequence
                      && cmp.IsSameFamily
                      && cmp.Matches == CanonicalSeedLength
                      && cmp.Mismatches == 0;

            return ok.Label($"same-seed: '{m1.SeedSequence}'='{m2.SeedSequence}' matches={cmp.Matches} family={cmp.IsSameFamily}");
        });
    }

    /// <summary>
    /// ONE-MISMATCH property: sequences whose seeds differ in exactly one position are
    /// never the same family and report Matches == 6, Mismatches == 1.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property OneMismatchSeedPair_NotFamily_WithSixMatches()
    {
        return Prop.ForAll(OneMismatchSeedPairArbitrary(), pair =>
        {
            var (a, b) = pair;
            var m1 = MiRnaAnalyzer.CreateMiRna("a", a);
            var m2 = MiRnaAnalyzer.CreateMiRna("b", b);
            var cmp = MiRnaAnalyzer.CompareSeedRegions(m1, m2);

            bool ok = m1.SeedSequence != m2.SeedSequence
                      && !cmp.IsSameFamily
                      && cmp.Matches == CanonicalSeedLength - 1
                      && cmp.Mismatches == 1;

            return ok.Label($"one-mismatch: '{m1.SeedSequence}' vs '{m2.SeedSequence}' matches={cmp.Matches} mismatches={cmp.Mismatches}");
        });
    }

    /// <summary>
    /// D (determinism): GetSeedSequence and CompareSeedRegions return identical results
    /// across repeated calls on the same inputs.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property SeedAnalysis_IsDeterministic()
    {
        var pairArbitrary =
            (from a in SequenceGen(0, 30)
             from b in SequenceGen(0, 30)
             select (a, b)).ToArbitrary();

        return Prop.ForAll(pairArbitrary, pair =>
        {
            var (a, b) = pair;
            string seedA1 = MiRnaAnalyzer.GetSeedSequence(a);
            string seedA2 = MiRnaAnalyzer.GetSeedSequence(a);

            var m1 = MiRnaAnalyzer.CreateMiRna("a", a);
            var m2 = MiRnaAnalyzer.CreateMiRna("b", b);
            var cmp1 = MiRnaAnalyzer.CompareSeedRegions(m1, m2);
            var cmp2 = MiRnaAnalyzer.CompareSeedRegions(m1, m2);

            return (seedA1 == seedA2 && cmp1.Equals(cmp2))
                .Label($"D: seed '{seedA1}'/'{seedA2}', cmp {cmp1} vs {cmp2}");
        });
    }

    // -- Anchored examples (deterministic, doc-derived) --

    /// <summary>
    /// Worked example (§7.1): hsa-let-7a-5p seed is "GAGGUAG" and is the same family as let-7b-5p.
    /// </summary>
    [Test]
    public void Anchor_WorkedExample_Let7a_SeedAndFamily()
    {
        var let7a = MiRnaAnalyzer.CreateMiRna("hsa-let-7a-5p", "UGAGGUAGUAGGUUGUAUAGUU");
        var let7b = MiRnaAnalyzer.CreateMiRna("hsa-let-7b-5p", "UGAGGUAGUAGGUUGUGUGGUU");

        Assert.That(MiRnaAnalyzer.GetSeedSequence(let7a.Sequence), Is.EqualTo("GAGGUAG"));
        Assert.That(let7a.SeedSequence, Is.EqualTo("GAGGUAG"));

        var cmp = MiRnaAnalyzer.CompareSeedRegions(let7a, let7b);
        Assert.That(cmp.IsSameFamily, Is.True);
        Assert.That(cmp.Matches, Is.EqualTo(7));
        Assert.That(cmp.Mismatches, Is.EqualTo(0));
    }

    /// <summary>Anchor: an input shorter than 8 nt yields the empty seed (INV-01 / §6.1).</summary>
    [TestCase("")]
    [TestCase("A")]
    [TestCase("UGAGGUA")] // exactly 7 nt — still too short
    public void Anchor_TooShortInput_YieldsEmptySeed(string seq)
    {
        Assert.That(MiRnaAnalyzer.GetSeedSequence(seq), Is.EqualTo(""));
    }

    /// <summary>Anchor: null input yields the empty seed (§6.1).</summary>
    [Test]
    public void Anchor_NullInput_YieldsEmptySeed()
    {
        Assert.That(MiRnaAnalyzer.GetSeedSequence(null!), Is.EqualTo(""));
    }

    /// <summary>Anchor: a same-seed pair scores Matches 7 / Mismatches 0 / family true.</summary>
    [Test]
    public void Anchor_SameSeedPair_SevenMatchesAndFamily()
    {
        var a = MiRnaAnalyzer.CreateMiRna("a", "AGAGGUAGCCCC"); // seed window indices 1..7 = "GAGGUAG"
        var b = MiRnaAnalyzer.CreateMiRna("b", "UGAGGUAGUUUU"); // seed window indices 1..7 = "GAGGUAG"

        Assert.That(a.SeedSequence, Is.EqualTo("GAGGUAG"));
        Assert.That(b.SeedSequence, Is.EqualTo("GAGGUAG"));

        var cmp = MiRnaAnalyzer.CompareSeedRegions(a, b);
        Assert.That(cmp.Matches, Is.EqualTo(7));
        Assert.That(cmp.Mismatches, Is.EqualTo(0));
        Assert.That(cmp.IsSameFamily, Is.True);
    }

    /// <summary>Anchor: a one-mismatch pair scores Matches 6 / Mismatches 1 / family false.</summary>
    [Test]
    public void Anchor_OneMismatchPair_SixMatchesNoFamily()
    {
        var a = MiRnaAnalyzer.CreateMiRna("a", "AGAGGUAG"); // seed = "GAGGUAG"
        var b = MiRnaAnalyzer.CreateMiRna("b", "AGAGGUAC"); // seed = "GAGGUAC" (last position differs)

        Assert.That(a.SeedSequence, Is.EqualTo("GAGGUAG"));
        Assert.That(b.SeedSequence, Is.EqualTo("GAGGUAC"));

        var cmp = MiRnaAnalyzer.CompareSeedRegions(a, b);
        Assert.That(cmp.Matches, Is.EqualTo(6));
        Assert.That(cmp.Mismatches, Is.EqualTo(1));
        Assert.That(cmp.IsSameFamily, Is.False);
    }

    #endregion
}
