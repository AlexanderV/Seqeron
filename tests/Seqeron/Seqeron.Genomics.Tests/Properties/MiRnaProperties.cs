using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for miRNA analysis: pre-miRNA hairpins and seed sequence analysis.
///
/// Test Units: MIRNA-PRECURSOR-001, MIRNA-SEED-001, MIRNA-PAIR-001
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

    #region MIRNA-TARGET-001

    // ---------------------------------------------------------------------
    // MIRNA-TARGET-001 — Target Site Prediction (MiRnaAnalyzer.FindTargetSites)
    //
    // Reference: docs/algorithms/MiRNA/Target_Site_Prediction.md
    //   §2.2 seedRC = revcomp(seed); 6mer core = seedRC[2..7] (0-based [1..6]),
    //        offset 6mer = seedRC[1..6] (0-based [0..5]); canonical site table.
    //   §2.4 INV-01 (0<=Score<=1), INV-02 (SeedMatchLength∈{6,7,8}),
    //        INV-03 (zero-based inclusive End=Start+len-1), INV-04 (canonical
    //        seed sites suppress overlapping offset-6mer), INV-05 (antiparallel).
    //   §5.2 BASE SCORES: 8mer=1.0, 7mer-m8=0.52, 7mer-A1=0.32, 6mer=0.15,
    //        offset-6mer=0.10; +0.05 when duplex Matches>10; -0.01 per mismatch;
    //        clamp to [0,1].
    //   §6.1 edge cases: empty/null mRNA, empty miRNA, seed<7nt ⇒ no sites;
    //        minScore>1.0 ⇒ no sites.
    //   §7.1 worked example: let-7a / "GGGGGCUACCUCAGGGGG" / minScore 0.1
    //        ⇒ sites[0] is Seed8mer with SeedMatchLength 8.
    //
    // ALL expected values are derived INDEPENDENTLY here: an independent
    // reverse-complement, an independent antiparallel WC + G:U duplex match/
    // mismatch counter, and an independent re-statement of the §5.2 score
    // formula. Nothing is read back from the implementation under test.
    // ---------------------------------------------------------------------

    /// <summary>Pure RNA alphabet for constructing valid (T→U normalized) miRNAs and mRNAs.</summary>
    private static readonly char[] RnaAlphabet = { 'A', 'C', 'G', 'U' };

    /// <summary>Independent reverse complement over the RNA alphabet (A↔U, G↔C); other → 'N'.</summary>
    private static string OracleRevComp(string s)
    {
        char[] r = new char[s.Length];
        for (int i = 0; i < s.Length; i++)
        {
            char c = char.ToUpperInvariant(s[s.Length - 1 - i]);
            r[i] = c switch { 'A' => 'U', 'U' => 'A', 'T' => 'A', 'G' => 'C', 'C' => 'G', _ => 'N' };
        }
        return new string(r);
    }

    /// <summary>Independent base-pairing predicate: Watson-Crick (A-U,U-A,G-C,C-G) or G:U wobble; DNA T≡U.</summary>
    private static bool OracleCanPair(char a, char b)
    {
        char x = a == 'T' ? 'U' : char.ToUpperInvariant(a);
        char y = b == 'T' ? 'U' : char.ToUpperInvariant(b);
        return (x == 'A' && y == 'U') || (x == 'U' && y == 'A') ||
               (x == 'G' && y == 'C') || (x == 'C' && y == 'G') ||
               (x == 'G' && y == 'U') || (x == 'U' && y == 'G');
    }

    /// <summary>True only for the two G:U wobble pairs (G-U, U-G).</summary>
    private static bool OracleIsWobble(char a, char b)
    {
        char x = a == 'T' ? 'U' : char.ToUpperInvariant(a);
        char y = b == 'T' ? 'U' : char.ToUpperInvariant(b);
        return (x == 'G' && y == 'U') || (x == 'U' && y == 'G');
    }

    /// <summary>
    /// Independent antiparallel duplex match/mismatch count over the first
    /// min(|miRNA|,|target|) positions: miRNA index i pairs with target[len-1-i]
    /// (INV-05). Mirrors <c>AlignMiRnaToTarget</c> EXACTLY: a Watson-Crick pair counts
    /// as a match, a G:U wobble counts as NEITHER match nor mismatch (it is tracked
    /// separately as a wobble in the implementation), and a non-pairing position is a
    /// mismatch. This distinction matters because the §5.2 score's "+0.05 when matches &gt; 10"
    /// bonus is driven by the Watson-Crick match count only.
    /// </summary>
    private static (int matches, int mismatches) OracleDuplex(string mirna, string target)
    {
        string m = mirna.ToUpperInvariant().Replace('T', 'U');
        string t = target.ToUpperInvariant().Replace('T', 'U');
        int len = Math.Min(m.Length, t.Length);
        int matches = 0, mismatches = 0;
        for (int i = 0; i < len; i++)
        {
            char a = m[i], b = t[t.Length - 1 - i];
            if (!OracleCanPair(a, b))
                mismatches++;
            else if (!OracleIsWobble(a, b))
                matches++;
            // G:U wobble: counted as neither match nor mismatch (matches the implementation).
        }
        return (matches, mismatches);
    }

    /// <summary>Independent §5.2 base score per emitted site class.</summary>
    private static double OracleBaseScore(MiRnaAnalyzer.TargetSiteType type) => type switch
    {
        MiRnaAnalyzer.TargetSiteType.Seed8mer => 1.0,
        MiRnaAnalyzer.TargetSiteType.Seed7merM8 => 0.52,
        MiRnaAnalyzer.TargetSiteType.Seed7merA1 => 0.32,
        MiRnaAnalyzer.TargetSiteType.Seed6mer => 0.15,
        MiRnaAnalyzer.TargetSiteType.Offset6mer => 0.10,
        _ => 0.01
    };

    /// <summary>
    /// Independent §5.2 score formula: base(class) + (Matches>10 ? +0.05) - 0.01·Mismatches,
    /// clamped to [0,1]. The duplex is the full extended target window (miRNA-length or tail).
    /// </summary>
    private static double OracleScore(MiRnaAnalyzer.TargetSiteType type, string mirna, string extendedTarget)
    {
        var (matches, mismatches) = OracleDuplex(mirna, extendedTarget);
        double s = OracleBaseScore(type);
        if (matches > 10) s += 0.05;
        s -= mismatches * 0.01;
        return Math.Max(0.0, Math.Min(1.0, s));
    }

    /// <summary>The extended target window the finder scores: from <paramref name="pos"/>, min(|miRNA|, tail).</summary>
    private static string ExtendedWindow(string mrna, int pos, int mirnaLen)
    {
        string normalized = mrna.ToUpperInvariant().Replace('T', 'U');
        int len = Math.Min(mirnaLen, normalized.Length - pos);
        return normalized.Substring(pos, len);
    }

    /// <summary>Generates random pure-RNA strings of length in [minLen,maxLen].</summary>
    private static Gen<string> RnaGen(int minLen, int maxLen) =>
        from len in Gen.Choose(minLen, maxLen)
        from chars in Gen.Elements(RnaAlphabet).ArrayOf(len)
        select new string(chars);

    /// <summary>A valid miRNA (length ≥ 8 ⇒ canonical 7-nt seed) and a random RNA mRNA to scan.</summary>
    private static Arbitrary<(MiRnaAnalyzer.MiRna miRna, string mrna)> MiRnaAndMrnaArbitrary() =>
        (from mirnaSeq in RnaGen(8, 26)
         from mrna in RnaGen(0, 60)
         select (MiRnaAnalyzer.CreateMiRna("rnd", mirnaSeq), mrna)).ToArbitrary();

    /// <summary>
    /// Builds an mRNA that embeds a PERFECT canonical 8mer for a miRNA that begins with 'U'.
    /// The target's 8-nt site window is exactly revcomp(miRNA[0..7]) placed at the mRNA tail, so
    /// the antiparallel duplex over those 8 positions is fully Watson-Crick paired: the
    /// position-8 base = comp(miRNA[7]) (= seedRC[0]) sits at the site start, the 6mer core
    /// follows, and the A1 slot = comp(miRNA[0]) = 'A' (since miRNA[0]=='U'). With the window
    /// length pinned to 8 (no tail beyond), there are zero duplex mismatches and no &gt;10-match
    /// bonus, so §5.2 gives an exact score of 1.0. The site starts at <c>leftPad.Length</c>.
    /// </summary>
    private static (string mrna, int siteStart) BuildClean8mer(string mirnaSeq, string leftPad)
    {
        // Window must equal revcomp(miRNA[0..7]) and be the whole tail (extended length == 8).
        string window = OracleRevComp(mirnaSeq.Substring(0, 8));
        return (leftPad + window, leftPad.Length);
    }

    /// <summary>
    /// Builds an mRNA with a CLEAN 6mer-only site: the 6mer core (revcomp of miRNA pos 2-7)
    /// is placed at target index 0 so no upstream position-8 base can exist, and the base
    /// immediately downstream is forced to a non-'A' so the A1 rule fails — yielding Seed6mer.
    /// The remainder of the window is the reverse complement of the miRNA so the duplex is
    /// fully paired (zero mismatches), making the §5.2 score exactly predictable.
    /// </summary>
    private static string BuildClean6mer(string mirnaSeq)
    {
        // Full antiparallel-perfect target, then trim everything left of the 6mer core so the
        // core starts at index 0 (kills position-8). revcomp(miRNA) places the core at index
        // |miRNA|-7; the base at index |miRNA|-1 (the A1 slot) is comp(miRNA[0]).
        string full = OracleRevComp(mirnaSeq);          // length = |miRNA|
        int coreIdx = mirnaSeq.Length - 7;              // documented core offset within revcomp
        string fromCore = full.Substring(coreIdx);      // 7 chars: core(6) + A1 slot(1)
        // fromCore[6] is the A1 slot = comp(miRNA[0]); force it to a guaranteed non-'A' base.
        char a1 = fromCore[6];
        char forcedNonA = a1 == 'C' ? 'G' : 'C';        // never 'A'
        return string.Concat(fromCore.AsSpan(0, 6), forcedNonA.ToString());
    }

    // -- INV-01 / INV-02 / INV-03 / P : invariants over arbitrary inputs --

    /// <summary>
    /// INV-01 (R: score ∈ [0,1]): every emitted site has 0.0 ≤ Score ≤ 1.0 for arbitrary
    /// (mRNA, miRNA). minScore is left at a permissive 0.0 so the widest set of sites is checked.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindTargetSites_Score_InUnitInterval()
    {
        return Prop.ForAll(MiRnaAndMrnaArbitrary(), input =>
        {
            var (miRna, mrna) = input;
            var sites = MiRnaAnalyzer.FindTargetSites(mrna, miRna, minScore: 0.0).ToList();
            bool ok = sites.All(s => s.Score >= 0.0 && s.Score <= 1.0);
            return ok.Label($"INV-01: scores={string.Join(",", sites.Select(s => s.Score))} mrna='{mrna}' seed='{miRna.SeedSequence}'");
        });
    }

    /// <summary>INV-02: every emitted site has SeedMatchLength ∈ {6,7,8}.</summary>
    [FsCheck.NUnit.Property]
    public Property FindTargetSites_SeedMatchLength_Is6Or7Or8()
    {
        return Prop.ForAll(MiRnaAndMrnaArbitrary(), input =>
        {
            var (miRna, mrna) = input;
            var sites = MiRnaAnalyzer.FindTargetSites(mrna, miRna, minScore: 0.0).ToList();
            bool ok = sites.All(s => s.SeedMatchLength is 6 or 7 or 8);
            return ok.Label($"INV-02: lengths={string.Join(",", sites.Select(s => s.SeedMatchLength))}");
        });
    }

    /// <summary>
    /// INV-03 (coordinates): 0 ≤ Start ≤ End &lt; |mRNA|, and the canonical-window length is
    /// consistent with End = Start + length - 1, where length is fixed by the site class
    /// (8mer→8, 7mer→7, 6mer/offset→6).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindTargetSites_Coordinates_AreConsistentAndInBounds()
    {
        return Prop.ForAll(MiRnaAndMrnaArbitrary(), input =>
        {
            var (miRna, mrna) = input;
            var sites = MiRnaAnalyzer.FindTargetSites(mrna, miRna, minScore: 0.0).ToList();

            foreach (var s in sites)
            {
                int classLen = s.Type switch
                {
                    MiRnaAnalyzer.TargetSiteType.Seed8mer => 8,
                    MiRnaAnalyzer.TargetSiteType.Seed7merM8 => 7,
                    MiRnaAnalyzer.TargetSiteType.Seed7merA1 => 7,
                    _ => 6 // Seed6mer / Offset6mer
                };
                bool ok = s.Start >= 0 && s.Start <= s.End && s.End < mrna.Length
                          && s.End == s.Start + classLen - 1;
                if (!ok)
                    return false.Label($"INV-03: Start={s.Start} End={s.End} type={s.Type} len={mrna.Length}");
            }
            return true.ToProperty();
        });
    }

    /// <summary>
    /// P (no spurious sites): every emitted site's target really contains the miRNA seed
    /// reverse-complement core. Independently recompute the 6mer core = revcomp(seed)[1..6]
    /// and the offset pattern = revcomp(seed)[0..5], then assert the mRNA carries that exact
    /// 6-mer at the documented offset of the reported site.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindTargetSites_EverySite_CarriesSeedReverseComplementCore()
    {
        return Prop.ForAll(MiRnaAndMrnaArbitrary(), input =>
        {
            var (miRna, mrna) = input;
            string normalized = mrna.ToUpperInvariant().Replace('T', 'U');
            string seedRC = OracleRevComp(miRna.SeedSequence);
            if (seedRC.Length < 7) return true.ToProperty();

            string sixmerCore = seedRC.Substring(1, 6); // RC of miRNA positions 2-7
            string offsetPat = seedRC.Substring(0, 6);  // RC of miRNA positions 3-8

            var sites = MiRnaAnalyzer.FindTargetSites(mrna, miRna, minScore: 0.0).ToList();
            foreach (var s in sites)
            {
                // Canonical classes: core starts one past Start for 8mer/7mer-m8 (which include
                // the upstream position-8 base), at Start for 7mer-A1/6mer. Offset uses Start.
                int coreOffset = s.Type switch
                {
                    MiRnaAnalyzer.TargetSiteType.Seed8mer => s.Start + 1,
                    MiRnaAnalyzer.TargetSiteType.Seed7merM8 => s.Start + 1,
                    MiRnaAnalyzer.TargetSiteType.Seed7merA1 => s.Start,
                    MiRnaAnalyzer.TargetSiteType.Seed6mer => s.Start,
                    _ => s.Start // Offset6mer
                };
                string expectedCore = s.Type == MiRnaAnalyzer.TargetSiteType.Offset6mer ? offsetPat : sixmerCore;
                string actual = normalized.Substring(coreOffset, 6);
                if (actual != expectedCore)
                    return false.Label($"P: site {s.Type}@{s.Start} core='{actual}' expected='{expectedCore}' in '{normalized}'");
            }
            return true.ToProperty();
        });
    }

    /// <summary>
    /// Score == independent §5.2 oracle: for every emitted site, the implementation's Score
    /// equals base(class) + (Matches>10?+0.05) - 0.01·Mismatches clamped to [0,1], where
    /// matches/mismatches are computed by the independent antiparallel duplex oracle over the
    /// extended target window. (Proves the documented scoring formula, not the code's output.)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindTargetSites_Score_MatchesDocFormula()
    {
        return Prop.ForAll(MiRnaAndMrnaArbitrary(), input =>
        {
            var (miRna, mrna) = input;
            var sites = MiRnaAnalyzer.FindTargetSites(mrna, miRna, minScore: 0.0).ToList();
            foreach (var s in sites)
            {
                string window = ExtendedWindow(mrna, s.Start, miRna.Sequence.Length);
                double expected = OracleScore(s.Type, miRna.Sequence, window);
                if (Math.Abs(s.Score - expected) > 1e-9)
                    return false.Label($"score: {s.Type}@{s.Start} got={s.Score} expected={expected} window='{window}'");
            }
            return true.ToProperty();
        });
    }

    // -- M : perfect seed match → higher score (central business invariant) --

    /// <summary>Generates a valid miRNA that BEGINS WITH 'U' (so revcomp(miRNA) yields a canonical 8mer).</summary>
    private static Arbitrary<MiRnaAnalyzer.MiRna> UStartMiRnaArbitrary() =>
        (from rest in Gen.Elements(RnaAlphabet).ArrayOf(11) // total length 12 ⇒ canonical 7-nt seed
         select MiRnaAnalyzer.CreateMiRna("u", "U" + new string(rest))).ToArbitrary();

    /// <summary>
    /// M (perfect seed match → higher score): for the SAME miRNA, a constructed perfect canonical
    /// 8mer always scores exactly 1.0 (§5.2, clamped) and strictly higher than a constructed
    /// clean 6mer-only site. Both site classes are built independently from the seed reverse
    /// complement; the 8mer/6mer classification and ordering are asserted directly.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PerfectSeedMatch_ScoresHigherThanSixmer()
    {
        return Prop.ForAll(UStartMiRnaArbitrary(), miRna =>
        {
            var (mrna8, start8) = BuildClean8mer(miRna.Sequence, "GG");
            string mrna6 = BuildClean6mer(miRna.Sequence);

            var site8 = MiRnaAnalyzer.FindTargetSites(mrna8, miRna, minScore: 0.0)
                .FirstOrDefault(s => s.Start == start8 && s.Type == MiRnaAnalyzer.TargetSiteType.Seed8mer);
            var site6 = MiRnaAnalyzer.FindTargetSites(mrna6, miRna, minScore: 0.0)
                .FirstOrDefault(s => s.Type == MiRnaAnalyzer.TargetSiteType.Seed6mer);

            bool eight = site8.Type == MiRnaAnalyzer.TargetSiteType.Seed8mer
                         && site8.SeedMatchLength == 8
                         && Math.Abs(site8.Score - 1.0) < 1e-9;
            bool six = site6.Type == MiRnaAnalyzer.TargetSiteType.Seed6mer
                       && site6.SeedMatchLength == 6;
            bool higher = site8.Score > site6.Score;

            return (eight && six && higher)
                .Label($"M: 8mer score={site8.Score} ({site8.Type}) vs 6mer score={site6.Score} ({site6.Type}) for {miRna.Sequence}");
        });
    }

    // -- minScore filtering: monotonic superset & threshold edges --

    /// <summary>
    /// minScore monotonic superset: lowering the threshold never removes a site. The site set
    /// at the higher threshold is a subset of the set at the lower threshold for the same
    /// (mRNA, miRNA). Compared on the full record so identity is exact.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindTargetSites_LowerMinScore_IsSuperset()
    {
        var arb = (from mirnaSeq in RnaGen(8, 26)
                   from mrna in RnaGen(0, 60)
                   from t1 in Gen.Choose(0, 100)
                   from t2 in Gen.Choose(0, 100)
                   select (MiRnaAnalyzer.CreateMiRna("rnd", mirnaSeq), mrna, t1 / 100.0, t2 / 100.0)).ToArbitrary();

        return Prop.ForAll(arb, input =>
        {
            var (miRna, mrna, ta, tb) = input;
            double hi = Math.Max(ta, tb);
            double lo = Math.Min(ta, tb);

            var atHi = MiRnaAnalyzer.FindTargetSites(mrna, miRna, hi).ToList();
            var atLo = MiRnaAnalyzer.FindTargetSites(mrna, miRna, lo).ToList();

            bool subset = atHi.All(s => atLo.Contains(s));
            return subset.Label($"superset: |hi({hi})|={atHi.Count} ⊄ |lo({lo})|={atLo.Count}");
        });
    }

    /// <summary>minScore &gt; 1.0 ⇒ no sites (all scores are clamped to ≤ 1.0; §6.1).</summary>
    [FsCheck.NUnit.Property]
    public Property FindTargetSites_MinScoreAboveOne_YieldsNoSites()
    {
        var arb = (from mirnaSeq in RnaGen(8, 26)
                   from mrna in RnaGen(0, 60)
                   select (MiRnaAnalyzer.CreateMiRna("rnd", mirnaSeq), mrna)).ToArbitrary();

        return Prop.ForAll(arb, input =>
        {
            var (miRna, mrna) = input;
            var sites = MiRnaAnalyzer.FindTargetSites(mrna, miRna, minScore: 1.0001).ToList();
            return (sites.Count == 0).Label($"minScore>1: got {sites.Count} sites");
        });
    }

    // -- INV-04 : canonical sites suppress overlapping offset-6mer --

    /// <summary>
    /// INV-04 (priority): a clean canonical 8mer covers positions that would otherwise also be
    /// read as an offset-6mer; no Offset6mer site may be emitted overlapping the canonical
    /// coordinates. Constructed on a miRNA beginning with U so the 8mer is guaranteed.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindTargetSites_CanonicalSite_SuppressesOverlappingOffset6mer()
    {
        return Prop.ForAll(UStartMiRnaArbitrary(), miRna =>
        {
            var (mrna, start8) = BuildClean8mer(miRna.Sequence, "GG");
            var sites = MiRnaAnalyzer.FindTargetSites(mrna, miRna, minScore: 0.0).ToList();

            var canonical = sites.FirstOrDefault(s => s.Start == start8 && s.Type == MiRnaAnalyzer.TargetSiteType.Seed8mer);
            if (canonical.Type != MiRnaAnalyzer.TargetSiteType.Seed8mer)
                return false.Label($"INV-04 setup: expected 8mer at {start8} in '{mrna}'");

            bool noOverlappingOffset = !sites.Any(s =>
                s.Type == MiRnaAnalyzer.TargetSiteType.Offset6mer &&
                s.Start <= canonical.End && s.End >= canonical.Start);

            return noOverlappingOffset
                .Label($"INV-04: offset-6mer overlaps canonical [{canonical.Start},{canonical.End}] in '{mrna}'");
        });
    }

    // -- D : determinism --

    /// <summary>
    /// D (determinism): identical (mRNA, miRNA, minScore) ⇒ identical site list, field-for-field
    /// and in the same order.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindTargetSites_IsDeterministic()
    {
        return Prop.ForAll(MiRnaAndMrnaArbitrary(), input =>
        {
            var (miRna, mrna) = input;
            var run1 = MiRnaAnalyzer.FindTargetSites(mrna, miRna, 0.0).ToList();
            var run2 = MiRnaAnalyzer.FindTargetSites(mrna, miRna, 0.0).ToList();
            return run1.SequenceEqual(run2).Label($"D: {run1.Count} vs {run2.Count} sites for '{mrna}'");
        });
    }

    // -- Edge cases (§6.1) --

    /// <summary>Edge: empty or null mRNA ⇒ no sites (§6.1).</summary>
    [TestCase("")]
    [TestCase(null)]
    public void Anchor_EmptyOrNullMrna_YieldsNoSites(string? mrna)
    {
        var let7a = MiRnaAnalyzer.CreateMiRna("let-7a", "UGAGGUAGUAGGUUGUAUAGUU");
        Assert.That(MiRnaAnalyzer.FindTargetSites(mrna!, let7a, 0.0).ToList(), Is.Empty);
    }

    /// <summary>Edge: empty miRNA sequence ⇒ no sites (§6.1).</summary>
    [Test]
    public void Anchor_EmptyMiRnaSequence_YieldsNoSites()
    {
        var empty = MiRnaAnalyzer.CreateMiRna("empty", "");
        Assert.That(MiRnaAnalyzer.FindTargetSites("GGGGGCUACCUCAGGGGG", empty, 0.0).ToList(), Is.Empty);
    }

    /// <summary>Edge: miRNA shorter than 8 nt has an empty seed (seedRC &lt; 7) ⇒ no sites (§6.1).</summary>
    [TestCase("U")]
    [TestCase("UGAGGUA")] // 7 nt — still too short to yield a canonical seed
    public void Anchor_SeedTooShort_YieldsNoSites(string shortMiRna)
    {
        var m = MiRnaAnalyzer.CreateMiRna("short", shortMiRna);
        Assert.That(m.SeedSequence, Is.Empty); // independent precondition: no canonical seed
        Assert.That(MiRnaAnalyzer.FindTargetSites("GGGGGCUACCUCAGGGGG", m, 0.0).ToList(), Is.Empty);
    }

    // -- Anchored examples (deterministic, doc-derived) --

    /// <summary>
    /// Worked example (§7.1): let-7a vs "GGGGGCUACCUCAGGGGG" with minScore 0.1 yields a first
    /// site of type Seed8mer with SeedMatchLength 8 (the documented assertion). Start/End are the
    /// independently derived seed-RC offsets (core "UACCUC" at index 6 ⇒ 8mer Start 5, End 12),
    /// and the Score equals the independent §5.2 oracle over the same extended window.
    /// </summary>
    [Test]
    public void Anchor_WorkedExample_Let7a_Yields8mer()
    {
        var let7a = MiRnaAnalyzer.CreateMiRna("let-7a", "UGAGGUAGUAGGUUGUAUAGUU");
        var sites = MiRnaAnalyzer.FindTargetSites("GGGGGCUACCUCAGGGGG", let7a, minScore: 0.1).ToList();

        Assert.That(sites, Is.Not.Empty);
        Assert.That(sites[0].Type, Is.EqualTo(MiRnaAnalyzer.TargetSiteType.Seed8mer));
        Assert.That(sites[0].SeedMatchLength, Is.EqualTo(8));
        Assert.That(sites[0].Start, Is.EqualTo(5));   // §2.2 worked offsets: core at 6 ⇒ 8mer start 5
        Assert.That(sites[0].End, Is.EqualTo(12));     // End = Start + 8 - 1

        string window = ExtendedWindow("GGGGGCUACCUCAGGGGG", 5, let7a.Sequence.Length);
        double expected = OracleScore(MiRnaAnalyzer.TargetSiteType.Seed8mer, let7a.Sequence, window);
        Assert.That(sites[0].Score, Is.EqualTo(expected).Within(1e-9));
    }

    /// <summary>
    /// §5.2 base-score anchor — perfect 8mer ⇒ exactly 1.0. Constructed so the antiparallel
    /// duplex is fully Watson-Crick paired (zero mismatches); the 8mer base 1.0 plus the
    /// &gt;10-match bonus both clamp to 1.0.
    /// </summary>
    [Test]
    public void Anchor_Clean8mer_ScoresExactlyOne()
    {
        var miRna = MiRnaAnalyzer.CreateMiRna("u8", "UACGUACGUACG"); // begins with U ⇒ 8mer
        var (mrna, start) = BuildClean8mer(miRna.Sequence, "GG");

        var site = MiRnaAnalyzer.FindTargetSites(mrna, miRna, 0.0)
            .Single(s => s.Start == start && s.Type == MiRnaAnalyzer.TargetSiteType.Seed8mer);

        Assert.That(site.Score, Is.EqualTo(1.0).Within(1e-9));
    }

    /// <summary>
    /// §5.2 base-score anchor — clean 6mer ⇒ base 0.15 minus 0.01 per duplex mismatch, with no
    /// &gt;10-match bonus (short window). The expected value is computed by the independent
    /// duplex oracle over the same window, then asserted against the implementation.
    /// </summary>
    [Test]
    public void Anchor_Clean6mer_ScoresPerDocFormula()
    {
        var miRna = MiRnaAnalyzer.CreateMiRna("m6", "GACGUCAUCGUACG");
        string mrna = BuildClean6mer(miRna.Sequence);

        var site = MiRnaAnalyzer.FindTargetSites(mrna, miRna, 0.0)
            .Single(s => s.Type == MiRnaAnalyzer.TargetSiteType.Seed6mer);

        string window = ExtendedWindow(mrna, site.Start, miRna.Sequence.Length);
        double expected = OracleScore(MiRnaAnalyzer.TargetSiteType.Seed6mer, miRna.Sequence, window);

        Assert.That(expected, Is.LessThanOrEqualTo(0.15)); // base 0.15, no >10-match bonus possible
        Assert.That(site.Score, Is.EqualTo(expected).Within(1e-9));
        Assert.That(site.Score, Is.LessThan(1.0)); // strictly below a perfect 8mer
    }

    #endregion

    #region MIRNA-PRECURSOR-001

    // ---------------------------------------------------------------------
    // MIRNA-PRECURSOR-001 — Pre-miRNA Hairpin Detection
    //   (MiRnaAnalyzer.FindPreMiRnaHairpins / private AnalyzeHairpin)
    //
    // Reference: docs/algorithms/MiRNA/Pre_miRNA_Detection.md
    //   §2.2 Core model: a candidate S of length n is accepted only if it
    //        decomposes as S = 5'stem + loop + 3'stem with UNINTERRUPTED
    //        mirrored pairing from the ends inward under the pairing set
    //        {A-U, U-A, G-C, C-G, G-U, U-G}; stem length ≥ 18; loop length
    //        = n − 2·stem with 3 ≤ loop ≤ 25; window length default 55–120.
    //   §2.4 INV-01 balanced dot-bracket, |Structure| = |Sequence|.
    //        INV-02 stem ≥ 18.   INV-03 loop ∈ [3,25].
    //        INV-04 zero-based inclusive coordinates, End = Start + len − 1.
    //        INV-05 |Mature| = |Star| = min(matureLength, stem).
    //   §6.1 edge cases: null/empty/too-short ⇒ none; stem < 18 ⇒ reject;
    //        loop < 3 or > 25 ⇒ reject.
    //   §7.1 worked example anchor.
    //
    // THEORY ENCODING (independent of the code under test):
    //   We CONSTRUCT perfect hairpins  Sequence = arm5 + ('A'·L) + rnaRevComp(arm5)
    //   where arm5 is `stem` random A/C/G/U bases and rnaRevComp is a strict
    //   Watson-Crick reverse complement (A↔U, G↔C). By construction every
    //   mirrored position i ∈ [0,stem) is a Watson-Crick pair, and the all-'A'
    //   loop boundary (A:A) does NOT pair, so the uninterrupted stem stops at
    //   EXACTLY `stem`. The expected dot-bracket, stem, loop, mature, star and
    //   coordinates are all derived HERE from the construction, never read back
    //   from the implementation.
    //
    // IMPLEMENTATION-DERIVED CONSTRAINT (honest, from AnalyzeHairpin source):
    //   The detector caps the scanned stem at maxStem = min(n/2 − 5, 35). For
    //   the constructed full-length window to be detected with the FULL intended
    //   stem we therefore require  stem ≤ n/2 − 5, i.e. loop ≥ 10 (and stem ≤ 35).
    //   The generator below honours this so the construction's stem is the one
    //   the detector reports. This is a true property of the source, not a
    //   weakening: the dot-bracket / stem / loop oracle is still derived from
    //   the perfect-hairpin theory.
    // ---------------------------------------------------------------------

    private const int MinStem = 18;          // §2.2 stem ≥ 18 bp
    private const int MaxStemCap = 35;        // AnalyzeHairpin hard cap on scanned stem
    private const int MinLoop = 3;            // §2.2 loop ≥ 3
    private const int MaxLoop = 25;           // §2.2 loop ≤ 25
    private const int MinWindow = 55;         // default minHairpinLength
    private const int MaxWindow = 120;        // default maxHairpinLength
    private const int DefaultMatureLength = 22;

    /// <summary>Pure RNA alphabet used to build hairpin arms.</summary>
    private static readonly char[] HairpinArmAlphabet = { 'A', 'C', 'G', 'U' };

    /// <summary>The allowed pre-miRNA pairing set {A-U, U-A, G-C, C-G, G-U, U-G} (§2.2) — same rule as <see cref="OracleCanPair"/>.</summary>
    private static bool IsAllowedHairpinPair(char a, char b) => OracleCanPair(a, b);

    /// <summary>
    /// Independent STRICT Watson-Crick reverse complement (A↔U, G↔C only — no wobble),
    /// so that <c>arm5 + rnaRevComp(arm5)</c> pairs every mirrored position perfectly.
    /// </summary>
    private static string RnaRevComp(string arm)
    {
        char[] r = new char[arm.Length];
        for (int i = 0; i < arm.Length; i++)
        {
            char c = arm[arm.Length - 1 - i];
            r[i] = c switch { 'A' => 'U', 'U' => 'A', 'G' => 'C', 'C' => 'G', _ => 'N' };
        }
        return new string(r);
    }

    /// <summary>The expected dot-bracket for a perfect hairpin: stem '(' + loop '.' + stem ')'.</summary>
    private static string ExpectedStructure(int stem, int loop) =>
        new string('(', stem) + new string('.', loop) + new string(')', stem);

    /// <summary>
    /// Builds a perfect hairpin <c>arm5 + ('A'·loop) + rnaRevComp(arm5)</c> from explicit
    /// stem and loop parameters using <paramref name="seed"/> for reproducible arm bases.
    /// Returns the sequence together with the theory-derived stem and loop it encodes.
    /// </summary>
    private static (string sequence, int stem, int loop) BuildPerfectHairpin(int stem, int loop, int seed)
    {
        var rng = new Random(seed);
        char[] arm = new char[stem];
        for (int i = 0; i < stem; i++)
            arm[i] = HairpinArmAlphabet[rng.Next(HairpinArmAlphabet.Length)];
        string arm5 = new string(arm);
        string sequence = arm5 + new string('A', loop) + RnaRevComp(arm5);
        return (sequence, stem, loop);
    }

    /// <summary>
    /// Generator of perfect hairpins whose FULL-length window is accepted with the EXACT
    /// constructed stem. Constraints: loop ∈ [10,25] (so stem ≤ n/2 − 5 holds and the all-'A'
    /// loop boundary cleanly terminates the stem), stem ∈ [18, min(35,(120−loop)/2)], and total
    /// length ∈ [55,120]. The seed drives reproducible arm bases.
    /// </summary>
    private static Arbitrary<(string sequence, int stem, int loop)> PerfectHairpinArbitrary() =>
        (from loop in Gen.Choose(10, MaxLoop)
         let hiStem = Math.Min(MaxStemCap, (MaxWindow - loop) / 2)
         from stem in Gen.Choose(MinStem, hiStem)
         from seed in Gen.Choose(1, 1_000_000)
         where 2 * stem + loop >= MinWindow && 2 * stem + loop <= MaxWindow
         select BuildPerfectHairpin(stem, loop, seed)).ToArbitrary();

    /// <summary>Arbitrary random RNA strings (mostly non-hairpins) of length 0..130.</summary>
    private static Arbitrary<string> ArbitraryRnaArbitrary() =>
        (from len in Gen.Choose(0, 130)
         from chars in Gen.Elements(HairpinArmAlphabet).ArrayOf(len)
         select new string(chars)).ToArbitrary();

    /// <summary>
    /// Independent re-statement of the full INV-01..INV-05 + R contract for ONE emitted
    /// candidate against the original scanned input. All expected values are recomputed here.
    /// </summary>
    private static (bool ok, string detail) ValidateCandidate(MiRnaAnalyzer.PreMiRna hp, string input, int matureLength)
    {
        // INV-04: coordinates zero-based inclusive, within the input, End = Start + len − 1.
        if (hp.Start < 0 || hp.End < hp.Start || hp.End >= input.Length)
            return (false, $"INV-04 bounds Start={hp.Start} End={hp.End} |input|={input.Length}");
        if (hp.End != hp.Start + hp.Sequence.Length - 1)
            return (false, $"INV-04 End {hp.End} ≠ Start+len-1 {hp.Start + hp.Sequence.Length - 1}");

        // The emitted Sequence must be the corresponding normalized (T→U, uppercase) substring.
        string expectedSub = input.ToUpperInvariant().Replace('T', 'U').Substring(hp.Start, hp.End - hp.Start + 1);
        if (hp.Sequence != expectedSub)
            return (false, $"INV-04 Sequence '{hp.Sequence}' ≠ input substring '{expectedSub}'");

        // INV-01: |Structure| = |Sequence|, balanced, of form (^k .^m )^k.
        if (hp.Structure.Length != hp.Sequence.Length)
            return (false, $"INV-01 |Structure|={hp.Structure.Length} ≠ |Sequence|={hp.Sequence.Length}");
        int opens = hp.Structure.Count(c => c == '(');
        int closes = hp.Structure.Count(c => c == ')');
        int dots = hp.Structure.Count(c => c == '.');
        if (opens != closes)
            return (false, $"INV-01 unbalanced {opens} '(' vs {closes} ')'");
        if (opens + closes + dots != hp.Structure.Length)
            return (false, $"INV-01 structure has unexpected characters: '{hp.Structure}'");
        int stem = opens; // leading '(' count
        int loop = dots;
        if (hp.Structure != ExpectedStructure(stem, loop))
            return (false, $"INV-01 structure '{hp.Structure}' ≠ canonical stem/loop form");

        // INV-02 / INV-03: stem ≥ 18, loop ∈ [3,25], and loop = n − 2·stem.
        if (stem < MinStem)
            return (false, $"INV-02 stem {stem} < {MinStem}");
        if (loop < MinLoop || loop > MaxLoop)
            return (false, $"INV-03 loop {loop} ∉ [{MinLoop},{MaxLoop}]");
        if (loop != hp.Sequence.Length - 2 * stem)
            return (false, $"loop {loop} ≠ n − 2·stem {hp.Sequence.Length - 2 * stem}");

        // The reported stem must be a REAL uninterrupted run of allowed mirrored pairs.
        for (int i = 0; i < stem; i++)
        {
            if (!IsAllowedHairpinPair(hp.Sequence[i], hp.Sequence[hp.Sequence.Length - 1 - i]))
                return (false, $"stem position {i} ({hp.Sequence[i]}:{hp.Sequence[hp.Sequence.Length - 1 - i]}) is not an allowed pair");
        }

        // INV-05: |Mature| = |Star| = min(matureLength, stem); Mature is the 5' arm prefix.
        int expectedArmLen = Math.Min(matureLength, stem);
        if (hp.MatureSequence.Length != expectedArmLen || hp.StarSequence.Length != expectedArmLen)
            return (false, $"INV-05 |Mature|={hp.MatureSequence.Length} |Star|={hp.StarSequence.Length} expected {expectedArmLen}");
        if (hp.MatureSequence != hp.Sequence.Substring(0, expectedArmLen))
            return (false, $"INV-05 Mature '{hp.MatureSequence}' ≠ 5' arm '{hp.Sequence.Substring(0, expectedArmLen)}'");
        if (hp.StarSequence != hp.Sequence.Substring(hp.Sequence.Length - expectedArmLen, expectedArmLen))
            return (false, $"INV-05 Star '{hp.StarSequence}' ≠ 3' arm");

        // R: precursor length strictly greater than mature length.
        if (hp.Sequence.Length <= hp.MatureSequence.Length)
            return (false, $"R precursor len {hp.Sequence.Length} ≤ mature len {hp.MatureSequence.Length}");

        return (true, "ok");
    }

    /// <summary>
    /// P (hairpin present) + INV-01..INV-05 + R on CONSTRUCTED perfect hairpins.
    /// For each built hairpin, the finder must emit a candidate whose Sequence equals the
    /// full constructed sequence, with the theory-derived dot-bracket (stem '(' + loop '.' +
    /// stem ')'), and that candidate satisfies the full independently-recomputed contract.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindPreMiRnaHairpins_PerfectHairpin_AcceptedWithExpectedStructure()
    {
        return Prop.ForAll(PerfectHairpinArbitrary(), built =>
        {
            var (sequence, stem, loop) = built;
            string expectedStructure = ExpectedStructure(stem, loop);

            var candidates = MiRnaAnalyzer.FindPreMiRnaHairpins(sequence).ToList();
            var full = candidates.Where(c => c.Sequence == sequence).ToList();

            if (full.Count == 0)
                return false.Label($"P: no full-length candidate for stem={stem} loop={loop} seq='{sequence}'");

            foreach (var hp in full)
            {
                if (hp.Structure != expectedStructure)
                    return false.Label($"INV-01: structure '{hp.Structure}' ≠ expected '{expectedStructure}'");
                if (hp.Structure.Length != sequence.Length)
                    return false.Label($"INV-01: |structure| {hp.Structure.Length} ≠ |seq| {sequence.Length}");
                if (hp.Structure.Count(c => c == '(') != hp.Structure.Count(c => c == ')'))
                    return false.Label("INV-01: unbalanced parentheses");

                // Independently verify the hairpin is real (mirrored allowed pairing across the stem).
                for (int i = 0; i < stem; i++)
                {
                    if (!IsAllowedHairpinPair(sequence[i], sequence[sequence.Length - 1 - i]))
                        return false.Label($"verify: stem pos {i} not an allowed pair");
                }

                var (ok, detail) = ValidateCandidate(hp, sequence, DefaultMatureLength);
                if (!ok)
                    return false.Label($"contract: {detail}");
            }

            return true.ToProperty();
        });
    }

    /// <summary>
    /// INV-05 + R focus on constructed hairpins: |Mature| = |Star| = min(22, stem), Mature is
    /// the 5' arm prefix, Star is the 3' arm suffix, and precursor length &gt; mature length.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindPreMiRnaHairpins_PerfectHairpin_MatureStarLengthsAndArms()
    {
        return Prop.ForAll(PerfectHairpinArbitrary(), built =>
        {
            var (sequence, stem, loop) = built;
            int expectedArmLen = Math.Min(DefaultMatureLength, stem);

            var full = MiRnaAnalyzer.FindPreMiRnaHairpins(sequence)
                .Where(c => c.Sequence == sequence).ToList();
            if (full.Count == 0)
                return false.Label($"setup: no full-length candidate stem={stem} loop={loop}");

            foreach (var hp in full)
            {
                bool ok = hp.MatureSequence.Length == expectedArmLen
                          && hp.StarSequence.Length == expectedArmLen
                          && hp.MatureSequence == sequence.Substring(0, expectedArmLen)
                          && hp.StarSequence == sequence.Substring(sequence.Length - expectedArmLen, expectedArmLen)
                          && hp.Sequence.Length > hp.MatureSequence.Length;
                if (!ok)
                    return false.Label($"INV-05/R: mature='{hp.MatureSequence}' star='{hp.StarSequence}' expectedLen={expectedArmLen} seqLen={hp.Sequence.Length}");
            }
            return true.ToProperty();
        });
    }

    /// <summary>
    /// INV-04 (coordinates) on constructed hairpins: the full-length candidate sits at Start=0,
    /// End=|seq|−1, and its Sequence is exactly the input window.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindPreMiRnaHairpins_PerfectHairpin_CoordinatesExact()
    {
        return Prop.ForAll(PerfectHairpinArbitrary(), built =>
        {
            var (sequence, _, _) = built;
            var hp = MiRnaAnalyzer.FindPreMiRnaHairpins(sequence)
                .FirstOrDefault(c => c.Sequence == sequence);
            if (hp.Sequence != sequence)
                return false.Label($"setup: no full-length candidate for '{sequence}'");

            bool ok = hp.Start == 0
                      && hp.End == sequence.Length - 1
                      && hp.Sequence == sequence.Substring(hp.Start, hp.End - hp.Start + 1);
            return ok.Label($"INV-04: Start={hp.Start} End={hp.End} |seq|={sequence.Length}");
        });
    }

    /// <summary>
    /// Robustness over ARBITRARY input: whatever candidates the finder emits on random RNA
    /// (mostly non-hairpins) all satisfy INV-01..INV-05 + R. Most random input yields none —
    /// that is acceptable; the assertion is "no malformed output".
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindPreMiRnaHairpins_ArbitraryInput_AllCandidatesWellFormed()
    {
        return Prop.ForAll(ArbitraryRnaArbitrary(), seq =>
        {
            var candidates = MiRnaAnalyzer.FindPreMiRnaHairpins(seq).ToList();
            foreach (var hp in candidates)
            {
                var (ok, detail) = ValidateCandidate(hp, seq, DefaultMatureLength);
                if (!ok)
                    return false.Label($"malformed candidate: {detail} (input '{seq}')");
            }
            return true.ToProperty();
        });
    }

    /// <summary>
    /// D (determinism): identical input ⇒ identical candidate list, field-for-field and in the
    /// same order, on both constructed hairpins and arbitrary RNA.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindPreMiRnaHairpins_IsDeterministic()
    {
        var arb = Gen.OneOf(
                PerfectHairpinArbitrary().Generator.Select(b => b.sequence),
                ArbitraryRnaArbitrary().Generator)
            .ToArbitrary();

        return Prop.ForAll(arb, seq =>
        {
            var run1 = MiRnaAnalyzer.FindPreMiRnaHairpins(seq).ToList();
            var run2 = MiRnaAnalyzer.FindPreMiRnaHairpins(seq).ToList();
            return run1.SequenceEqual(run2).Label($"D: {run1.Count} vs {run2.Count} candidates for '{seq}'");
        });
    }

    // -- Edge cases (§6.1) and worked-example anchor (§7.1) --

    /// <summary>Edge (§6.1): null or empty input yields no candidates.</summary>
    [TestCase(null)]
    [TestCase("")]
    public void Anchor_NullOrEmptyInput_YieldsNoCandidates(string? seq)
    {
        Assert.That(MiRnaAnalyzer.FindPreMiRnaHairpins(seq!).ToList(), Is.Empty);
    }

    /// <summary>Edge (§6.1): input shorter than minHairpinLength (55) yields no candidates.</summary>
    [Test]
    public void Anchor_InputShorterThanMinWindow_YieldsNoCandidates()
    {
        // A would-be perfect hairpin of total length 54 (< 55) — too short to be scanned.
        var (shortSeq, _, _) = BuildPerfectHairpin(stem: 22, loop: 10, seed: 7); // length 54
        Assert.That(shortSeq.Length, Is.EqualTo(54));
        Assert.That(MiRnaAnalyzer.FindPreMiRnaHairpins(shortSeq).ToList(), Is.Empty);
    }

    /// <summary>
    /// Edge (§6.1): a hairpin with stem 17 (just below the 18 bp minimum) is rejected — no
    /// candidate equal to that constructed sequence is emitted. Loop padded so total ≥ 55.
    /// </summary>
    [Test]
    public void Anchor_StemBelowMinimum_IsRejected()
    {
        // stem 17, loop 21 ⇒ length 55. Detected uninterrupted stem = 17 < 18 ⇒ reject.
        var (seq, _, _) = BuildPerfectHairpin(stem: 17, loop: 21, seed: 3);
        Assert.That(seq.Length, Is.EqualTo(55));
        var full = MiRnaAnalyzer.FindPreMiRnaHairpins(seq).Where(c => c.Sequence == seq).ToList();
        Assert.That(full, Is.Empty, "stem-17 hairpin must be rejected (INV-02)");
    }

    /// <summary>
    /// Edge (§6.1, INV-03 upper bound): a hairpin whose loop is 26 (above the 25-nt maximum) is
    /// rejected — the constructed full-length sequence must not appear among emitted candidates.
    /// (stem 22, loop 26 ⇒ length 70; the detected uninterrupted stem is 22, leaving loop 26.)
    /// </summary>
    [Test]
    public void Anchor_LoopAboveMaximum_IsRejected()
    {
        var (seq, _, _) = BuildPerfectHairpin(stem: 22, loop: 26, seed: 11);
        Assert.That(seq.Length, Is.EqualTo(70));
        var full = MiRnaAnalyzer.FindPreMiRnaHairpins(seq).Where(c => c.Sequence == seq).ToList();
        Assert.That(full, Is.Empty, "loop-26 hairpin must be rejected (INV-03)");
    }

    /// <summary>
    /// INV-03 lower bound — honest reachability note. The detector caps the scanned stem at
    /// maxStem = n/2 − 5, so the reported loop is ALWAYS ≥ 10 at any accepted window; a loop
    /// below 3 is therefore unreachable from the scan. A perfect hairpin constructed with a tiny
    /// loop (here 2) is consequently NOT rejected: it is re-decomposed into a stem capped at
    /// maxStem with a valid loop ≥ 3. This test asserts that genuine behaviour (any emitted
    /// candidate still honours INV-03's lower bound), proving the bound holds rather than
    /// fabricating an unreachable loop&lt;3 acceptance.
    /// </summary>
    [Test]
    public void Anchor_TinyConstructedLoop_ReDecomposedWithinLoopBounds()
    {
        // stem 27, loop 2 ⇒ length 56; maxStem = 56/2 − 5 = 23 caps the stem, loop becomes 10.
        var (seq, _, _) = BuildPerfectHairpin(stem: 27, loop: 2, seed: 11);
        Assert.That(seq.Length, Is.EqualTo(56));

        var candidates = MiRnaAnalyzer.FindPreMiRnaHairpins(seq).ToList();
        foreach (var hp in candidates)
        {
            int loop = hp.Structure.Count(c => c == '.');
            Assert.That(loop, Is.InRange(MinLoop, MaxLoop), "every emitted loop must satisfy INV-03");
            var (ok, detail) = ValidateCandidate(hp, seq, DefaultMatureLength);
            Assert.That(ok, Is.True, detail);
        }
    }

    /// <summary>
    /// Worked-example anchor (§7.1): the documented 57-nt sequence yields at least one candidate,
    /// the full-length window is accepted, and that candidate satisfies the full INV contract.
    /// Independently: its uninterrupted stem is 23, loop 11, dot-bracket = (^23 .^11 )^23.
    /// </summary>
    [Test]
    public void Anchor_WorkedExample_YieldsValidHairpin()
    {
        string sequence =
            "GCAUAGCUAGCUAGCUAGCUAGCUA" +
            "GAAAUUU" +
            "UAGCUAGCUAGCUAGCUAGCUAUGC";

        var candidates = MiRnaAnalyzer.FindPreMiRnaHairpins(sequence).ToList();
        Assert.That(candidates, Is.Not.Empty);

        var full = candidates.Single(c => c.Sequence == sequence.ToUpperInvariant().Replace('T', 'U'));
        Assert.That(full.Structure, Is.EqualTo(ExpectedStructure(23, 11)));
        Assert.That(full.MatureSequence, Has.Length.EqualTo(Math.Min(DefaultMatureLength, 23)));
        Assert.That(full.MatureSequence, Is.EqualTo(full.Sequence.Substring(0, full.MatureSequence.Length)));
        Assert.That(full.Sequence.Length, Is.GreaterThan(full.MatureSequence.Length));

        var (ok, detail) = ValidateCandidate(full, sequence, DefaultMatureLength);
        Assert.That(ok, Is.True, detail);
    }

    #endregion

    #region MIRNA-PAIR-001: P: seed region paired; R: alignment counts ≥ 0; D: deterministic

    // AlignMiRnaToTarget aligns a miRNA against an mRNA target and reports match/mismatch/wobble/gap
    // counts. A target that is the reverse complement of the miRNA pairs fully (including the seed).

    /// <summary>Generates an RNA miRNA of length 18..24 over {A,C,G,U}.</summary>
    private static Arbitrary<string> MiRnaRnaArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            int len = 18 + rng.Next(7);
            var c = new char[len];
            for (int i = 0; i < len; i++) c[i] = "ACGU"[rng.Next(4)];
            return new string(c);
        }).ToArbitrary();

    /// <summary>
    /// INV-1 (R): all duplex counts (matches, mismatches, wobbles, gaps) are non-negative for any
    /// miRNA/target pair.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MiRnaPair_CountsAreNonNegative()
    {
        return Prop.ForAll(MiRnaRnaArbitrary(), MiRnaRnaArbitrary(), (mirna, target) =>
        {
            var d = MiRnaAnalyzer.AlignMiRnaToTarget(mirna, target);
            return (d.Matches >= 0 && d.Mismatches >= 0 && d.GUWobbles >= 0 && d.Gaps >= 0)
                .Label($"negative duplex count: M={d.Matches} X={d.Mismatches} GU={d.GUWobbles} gaps={d.Gaps}");
        });
    }

    /// <summary>
    /// INV-2 (P): a target that is the reverse complement of the miRNA pairs across the whole length —
    /// in particular the 7-nt seed region pairs — with no mismatches.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MiRnaPair_ReverseComplementTarget_PairsSeed()
    {
        return Prop.ForAll(MiRnaRnaArbitrary(), mirna =>
        {
            string target = MiRnaAnalyzer.GetReverseComplement(mirna);
            var d = MiRnaAnalyzer.AlignMiRnaToTarget(mirna, target);
            // Watson-Crick pairs are counted as matches; the seed (positions 2-8) is within these.
            return (d.Matches >= 7 && d.Mismatches == 0)
                .Label($"reverse-complement target did not fully pair the seed: M={d.Matches} X={d.Mismatches}");
        });
    }

    /// <summary>
    /// INV-3 (D): Duplex alignment is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MiRnaPair_IsDeterministic()
    {
        return Prop.ForAll(MiRnaRnaArbitrary(), MiRnaRnaArbitrary(), (mirna, target) =>
        {
            var a = MiRnaAnalyzer.AlignMiRnaToTarget(mirna, target);
            var b = MiRnaAnalyzer.AlignMiRnaToTarget(mirna, target);
            return (a.Matches == b.Matches && a.Mismatches == b.Mismatches && a.AlignmentString == b.AlignmentString)
                .Label("AlignMiRnaToTarget must be deterministic");
        });
    }

    #endregion

    #region MIRNA-CONTEXT-001: R: context++ score ≤ 0 (more negative = stronger); D: deterministic

    // ScoreTargetSiteContextPlusPlus — TargetScan context++ (Agarwal et al. 2015): an additive log-fold-change
    // repression model whose (partial) context score is ≤ 0 for canonical seed sites (more negative = stronger).

    /// <summary>
    /// INV-1 (R): the partial context++ score of a perfect canonical 8mer site is ≤ 0 (repression is a
    /// non-positive log fold-change; a stronger site is more negative).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ContextPlusPlus_8merScore_IsNonPositive()
    {
        return Prop.ForAll(UStartMiRnaArbitrary(), miRna =>
        {
            var (mrna8, start8) = BuildClean8mer(miRna.Sequence, "GG");
            var site = MiRnaAnalyzer.FindTargetSites(mrna8, miRna, minScore: 0.0)
                .FirstOrDefault(s => s.Start == start8 && s.Type == MiRnaAnalyzer.TargetSiteType.Seed8mer);
            if (site.Type != MiRnaAnalyzer.TargetSiteType.Seed8mer) return true.ToProperty();
            var ctx = MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus(mrna8, miRna, site);
            return (ctx.ContextScorePartial <= 1e-9).Label($"context++ partial {ctx.ContextScorePartial} must be ≤ 0");
        });
    }

    /// <summary>INV-2 (D): context++ scoring is deterministic.</summary>
    [FsCheck.NUnit.Property]
    public Property ContextPlusPlus_IsDeterministic()
    {
        return Prop.ForAll(UStartMiRnaArbitrary(), miRna =>
        {
            var (mrna8, start8) = BuildClean8mer(miRna.Sequence, "GG");
            var site = MiRnaAnalyzer.FindTargetSites(mrna8, miRna, minScore: 0.0)
                .FirstOrDefault(s => s.Start == start8 && s.Type == MiRnaAnalyzer.TargetSiteType.Seed8mer);
            if (site.Type != MiRnaAnalyzer.TargetSiteType.Seed8mer) return true.ToProperty();
            return (MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus(mrna8, miRna, site).ContextScorePartial
                    == MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus(mrna8, miRna, site).ContextScorePartial)
                .Label("ScoreTargetSiteContextPlusPlus must be deterministic");
        });
    }

    #endregion

    #region MIRNA-PCT-001: R: PCT ∈ [0,1]; M: higher branch length → higher PCT; D: deterministic

    // PctFromBranchLength — TargetScan PCT sigmoid (Friedman et al. 2009): PCT = B0 + B1/(1+e^(−B2·BL+B3)),
    // clamped to ≥ 0. Generated over the published parameter regime (B0+B1 ≤ 1, B2 > 0) where PCT ∈ [0,1].

    private static Arbitrary<(double bl, MiRnaAnalyzer.PctSigmoidParameters p)> PctArbitrary() =>
        (from blCenti in Gen.Choose(0, 500)
         from b0c in Gen.Choose(0, 10)
         from b1c in Gen.Choose(5, 90)
         from b2c in Gen.Choose(1, 200)
         from b3c in Gen.Choose(-300, 300)
         select (blCenti / 100.0,
                 new MiRnaAnalyzer.PctSigmoidParameters(b0c / 100.0, b1c / 100.0, b2c / 100.0, b3c / 100.0)))
        .ToArbitrary();

    /// <summary>INV-1 (R): PCT lies in [0,1] over the published parameter regime (B0+B1 ≤ 1).</summary>
    [FsCheck.NUnit.Property]
    public Property Pct_InUnitInterval()
    {
        return Prop.ForAll(PctArbitrary(), t =>
        {
            double pct = MiRnaAnalyzer.PctFromBranchLength(t.bl, t.p);
            return (pct >= -1e-9 && pct <= 1.0 + 1e-9).Label($"PCT {pct} outside [0,1]");
        });
    }

    /// <summary>INV-2 (M): with B2 &gt; 0 and B1 &gt; 0, a longer branch length gives a higher PCT (the sigmoid is increasing).</summary>
    [FsCheck.NUnit.Property]
    public Property Pct_MonotoneInBranchLength()
    {
        var gen = (from t in PctArbitrary().Generator
                   from extraCenti in Gen.Choose(1, 300)
                   select (t.p, lo: t.bl, hi: t.bl + extraCenti / 100.0)).ToArbitrary();
        return Prop.ForAll(gen, t =>
        {
            double pctLo = MiRnaAnalyzer.PctFromBranchLength(t.lo, t.p);
            double pctHi = MiRnaAnalyzer.PctFromBranchLength(t.hi, t.p);
            return (pctHi >= pctLo - 1e-9).Label($"PCT not monotone: BL {t.lo}→{pctLo}, {t.hi}→{pctHi}");
        });
    }

    /// <summary>INV-3 (D): PCT is deterministic.</summary>
    [FsCheck.NUnit.Property]
    public Property Pct_IsDeterministic()
    {
        return Prop.ForAll(PctArbitrary(), t =>
            (MiRnaAnalyzer.PctFromBranchLength(t.bl, t.p) == MiRnaAnalyzer.PctFromBranchLength(t.bl, t.p))
                .Label("PctFromBranchLength must be deterministic"));
    }

    #endregion

    #region MIRNA-CLASSIFY-001: R: probability ∈ [0,1]; D: deterministic; threshold split positive/negative

    // ClassifyPreMiRna — bundled logistic-regression P(natural) ∈ [0,1]; IsNatural = P ≥ threshold.

    /// <summary>
    /// INV-1 (R + threshold): for any predicted pre-miRNA hairpin the natural probability is in [0,1] and the
    /// positive/negative call equals (probability ≥ threshold).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ClassifyPreMiRna_Probability_InUnitInterval_AndThresholdSplit()
    {
        return Prop.ForAll(PerfectHairpinArbitrary(), built =>
        {
            const double threshold = 0.5;
            var c = MiRnaAnalyzer.ClassifyPreMiRna(built.sequence, threshold);
            if (c is null) return true.ToProperty();
            var v = c.Value;
            bool ok = v.NaturalProbability is >= 0.0 and <= 1.0
                      && v.IsNatural == (v.NaturalProbability >= threshold);
            return ok.Label($"P={v.NaturalProbability}, IsNatural={v.IsNatural}");
        });
    }

    /// <summary>INV-2 (D): pre-miRNA classification is deterministic.</summary>
    [FsCheck.NUnit.Property]
    public Property ClassifyPreMiRna_IsDeterministic()
    {
        return Prop.ForAll(PerfectHairpinArbitrary(), built =>
        {
            var a = MiRnaAnalyzer.ClassifyPreMiRna(built.sequence);
            var b = MiRnaAnalyzer.ClassifyPreMiRna(built.sequence);
            return ((a is null && b is null)
                    || (a is not null && b is not null
                        && a.Value.NaturalProbability == b.Value.NaturalProbability
                        && a.Value.IsNatural == b.Value.IsNatural))
                .Label("ClassifyPreMiRna must be deterministic");
        });
    }

    #endregion

    #region MIRNA-CLEAVAGE-001: R: cleavage positions within precursor; R: 2-nt 3' overhang; D: deterministic

    // PredictDroshaDicerCleavage — Drosha basal-junction ruler (~11 bp, Han 2006) + Dicer 5'-counting (~22 nt,
    // Park 2011); the RNase III staggered cut leaves a 2-nt 3' overhang. Returns null when the geometry runs
    // off the sequence.

    private static Arbitrary<(string seq, int basal)> CleavageArbitrary() =>
        (from len in Gen.Choose(40, 130)
         from chars in Gen.Elements('A', 'C', 'G', 'U').ArrayOf(len)
         from basal in Gen.Choose(0, len / 3)
         select (new string(chars), basal)).ToArbitrary();

    /// <summary>
    /// INV-1 (R): whenever a cleavage is predicted, every reported cut/span coordinate lies within the
    /// precursor and the 3' overhang is exactly 2 nt (the RNase III staggered-cut signature).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DroshaDicerCleavage_PositionsValid_And2ntOverhang()
    {
        return Prop.ForAll(CleavageArbitrary(), t =>
        {
            var c = MiRnaAnalyzer.PredictDroshaDicerCleavage(t.seq, t.basal);
            if (c is null) return true.ToProperty();
            var v = c.Value;
            int n = t.seq.Length;
            bool inBounds = new[] { v.BasalJunction, v.DroshaCut5Prime, v.DroshaCut3Prime,
                                    v.MatureStart, v.MatureEnd, v.StarStart, v.StarEnd }
                .All(p => p >= 0 && p <= n);
            return (inBounds && v.ThreePrimeOverhang == 2)
                .Label($"positions in-bounds={inBounds}, overhang={v.ThreePrimeOverhang}");
        });
    }

    /// <summary>INV-2 (D): Drosha/Dicer cleavage prediction is deterministic.</summary>
    [FsCheck.NUnit.Property]
    public Property DroshaDicerCleavage_IsDeterministic()
    {
        return Prop.ForAll(CleavageArbitrary(), t =>
            (MiRnaAnalyzer.PredictDroshaDicerCleavage(t.seq, t.basal)
                == MiRnaAnalyzer.PredictDroshaDicerCleavage(t.seq, t.basal))
                .Label("PredictDroshaDicerCleavage must be deterministic"));
    }

    #endregion
}
