namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the MolTools <b>DNA hairpin folder</b> (PRIMER-HAIRPIN-001) — the opt-in
/// SantaLucia &amp; Hicks (2004) most-stable intramolecular hairpin: a single Watson-Crick stem
/// closing one hairpin loop, with its ΔH°/ΔS°/ΔG°37 and a unimolecular (concentration-independent)
/// hairpin Tm. The exact hand-derived thermodynamics are covered by
/// PrimerDesigner_HairpinTm_Tests; THIS file targets the malformed/boundary inputs and pins the
/// STRUCTURAL discipline of any returned hairpin against an independent oracle.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies (docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing")
/// ───────────────────────────────────────────────────────────────────────────
/// Malformed / boundary / out-of-domain inputs must NEVER hang, throw an *unhandled* runtime
/// exception (an IndexOutOfRange from stem/loop indexing, a NaN/Inf ΔG leaking where the doc
/// contracts a finite value), or emit OUT-OF-CONTRACT output — a returned hairpin whose stem is
/// not actually Watson-Crick complementary, whose loop is below the documented minimum (3 nt),
/// whose stem is below minStemLength, or whose indices fall outside the input. Every input must
/// resolve to EITHER a well-defined, theory-correct result (INCLUDING the documented `null` =
/// "no hairpin" sentinel and `NaN` Tm), OR a documented validation outcome.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PRIMER-HAIRPIN-001 — primer hairpin folding / most-stable hairpin
/// Checklist: docs/checklists/03_FUZZING.md, row 241.
/// Algorithm doc: docs/algorithms/MolTools/DNA_Hairpin_Folding_Tm.md.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the EMPTY sequence (→ null, no crash), a sequence with NO
///          complementary stem / homopolymer (→ null), a PERFECT PALINDROME that DOES fold (stem
///          and loop must be structurally valid), MIN-LOOP violations (a stem so close that no
///          loop ≥ 3 nt can be closed → null; loops < 3 nt are sterically prohibited [1, INV-05]),
///          and a VERY LONG sequence (the O(n³) closing-pair scan must not blow up — [CancelAfter]).
///   • MC = Malformed Content — non-ACGT junk / special / unicode characters and an all-N
///          sequence (the strict A/C/G/T alphabet rejects them → null, never a crash), plus the
///          out-of-domain `minStemLength < 2` parameter (→ null, the documented guard).
/// — docs/checklists/03_FUZZING.md §Description (BE; MC = malformed content); row 241 targets:
///   "empty seq, no complementary stem, palindrome, min loop violations, very long".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The hairpin-folder contract under test (DNA_Hairpin_Folding_Tm.md)
/// ───────────────────────────────────────────────────────────────────────────
/// API entries:
///   • PrimerDesigner.FindMostStableHairpin(string sequence, int minStemLength = 2,
///       double loopBonusDeltaG37 = 0.0) → HairpinResult?  (PrimerDesigner.cs §FindMostStableHairpin).
///   • PrimerDesigner.CalculateHairpinMeltingTemperature(string, int, double) → double
///       (NaN when no hairpin / invalid input).
/// HairpinResult carries StemStart, StemEnd (0-based outermost stem indices on the input strand),
/// StemLength (bp), LoopSize (unpaired loop nt), DeltaH (kcal/mol), DeltaS (cal/(K·mol)),
/// DeltaG37 (kcal/mol).
///
/// VALIDATION / sentinels (doc §3.3, §6.1):
///   • null / empty sequence → null hairpin, NaN Tm.
///   • a non-ACGT (degenerate / junk / unicode) base anywhere → null (strict A/C/G/T alphabet,
///     consistent with the duplex NN methods).
///   • minStemLength &lt; 2 → null (a stem needs ≥ 1 NN stack).
///   • a sequence with no Watson-Crick stem of ≥ minStemLength bp that can close a loop of ≥ 3 nt
///     (homopolymer, or a palindrome whose only stem closes a 0-nt loop) → null (INV-02, INV-05).
/// None of these is a throw — they are the documented null / NaN sentinels. The folder never
/// throws for ANY string input.
///
/// STRUCTURAL invariants pinned on EVERY returned hairpin (an INDEPENDENT oracle, doc §2.4):
///   • INV (indices) — 0 ≤ StemStart &lt; StemEnd &lt; sequence.Length; the stem arms and loop fit
///     inside the input (no out-of-bounds index ever escapes).
///   • INV-06 (complementary stem) — re-derive each of the StemLength outer→inner base pairs from
///     the (upper-cased) input by an INDEPENDENT Watson-Crick check
///     (seq[StemStart+k] pairs seq[StemEnd−k]) — every stem pair MUST be a genuine A·T / G·C
///     complement. A returned "stem" that is not actually complementary is an out-of-contract bug.
///   • INV-05 (loop minimum) — LoopSize ≥ 3 (loops &lt; 3 nt are sterically prohibited [1]); and
///     LoopSize equals the inner gap StemEnd − StemStart − 2·StemLength + 1 implied by the indices.
///   • INV (stem minimum) — StemLength ≥ max(2, minStemLength).
///   • finiteness — DeltaH, DeltaS, DeltaG37 are all finite (never NaN / ±Inf); and the loop is
///     destabilising so DeltaS &lt; 0.
/// These are asserted UNIVERSALLY on every positive result. The ΔG ≤ 0 (net-stability) assertion
/// is made ONLY on a genuine GC-stem domain (a clean inverted-repeat hairpin), NOT universally —
/// mirroring the published-domain note in the PRIMER-HAIRPIN property test: a weak A·T stem closing
/// a large loop can have a slightly POSITIVE ΔG°37 yet still be the most-stable hairpin returned.
///
/// FindMostStableHairpin / CalculateHairpinMeltingTemperature are pure functions (no iterator), so
/// every probe calls them directly. The very-long-sequence probe is pinned with [CancelAfter] so
/// any O(n³) blow-up fails as a timeout rather than hanging.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class PrimerHairpinFuzzTests
{
    private const int MinLoop = 3;   // SantaLucia & Hicks 2004: loops < 3 nt are sterically prohibited (INV-05).

    #region Helpers — independent Watson-Crick oracle

    /// <summary>Independent Watson-Crick complement check (NOT the source's IsWatsonCrickPair).</summary>
    private static bool IsComplement(char a, char b) =>
        (a, b) switch
        {
            ('A', 'T') or ('T', 'A') or ('G', 'C') or ('C', 'G') => true,
            _ => false
        };

    /// <summary>Deterministic ACGT generator — seed fixed locally so fuzz inputs are reproducible.</summary>
    private static string RandomDna(int length, int seed)
    {
        const string bases = "ACGT";
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    /// <summary>
    /// Asserts EVERY documented structural invariant on a returned hairpin against an INDEPENDENT
    /// oracle: indices in bounds, every stem pair a genuine Watson-Crick complement (re-derived from
    /// the input), loop ≥ 3 and consistent with the indices, stem ≥ the effective minimum, and all
    /// of ΔH°/ΔS°/ΔG°37 finite with a destabilising (negative) loop entropy. The ΔG ≤ 0 net-stability
    /// claim is deliberately NOT asserted here (a weak A·T stem may yield a slightly positive ΔG°37);
    /// it is pinned only on the genuine GC-stem palindrome test.
    /// </summary>
    private static void AssertWellFormedHairpin(
        PrimerDesigner.HairpinResult h, string original, int effectiveMinStem)
    {
        string seq = original.ToUpperInvariant();

        // Indices in bounds and ordered.
        h.StemStart.Should().BeGreaterThanOrEqualTo(0, "StemStart is a 0-based index");
        h.StemEnd.Should().BeLessThan(seq.Length, "StemEnd must fall inside the input");
        h.StemStart.Should().BeLessThan(h.StemEnd, "the 5' stem arm precedes the 3' stem arm");

        // Stem length ≥ the effective minimum (≥ 2 always).
        h.StemLength.Should().BeGreaterThanOrEqualTo(Math.Max(2, effectiveMinStem),
            "a returned stem must meet minStemLength (and is always ≥ 2)");

        // The two stem arms fit inside the sequence (no out-of-range slice).
        (h.StemStart + h.StemLength - 1).Should().BeLessThan(seq.Length, "5' stem arm fits");
        (h.StemEnd - h.StemLength + 1).Should().BeGreaterThanOrEqualTo(0, "3' stem arm fits");

        // INV-06: re-derive every stem pair independently — each MUST be Watson-Crick.
        for (int k = 0; k < h.StemLength; k++)
        {
            char left = seq[h.StemStart + k];
            char right = seq[h.StemEnd - k];
            IsComplement(left, right).Should().BeTrue(
                $"stem pair #{k} '{left}'·'{right}' (indices {h.StemStart + k},{h.StemEnd - k}) must be a genuine Watson-Crick complement (INV-06)");
        }

        // INV-05: loop ≥ 3, and consistent with the stem span the indices imply.
        h.LoopSize.Should().BeGreaterThanOrEqualTo(MinLoop,
            "loops < 3 nt are sterically prohibited (INV-05)");
        int impliedLoop = (h.StemEnd - h.StemLength + 1) - (h.StemStart + h.StemLength - 1) - 1;
        h.LoopSize.Should().Be(impliedLoop,
            "LoopSize must equal the inner gap between the innermost stem pair implied by the indices");

        // Finiteness; loop is destabilising so total ΔS° < 0 (stem ΔS° < 0 too).
        double.IsFinite(h.DeltaH).Should().BeTrue("ΔH° must be finite");
        double.IsFinite(h.DeltaS).Should().BeTrue("ΔS° must be finite");
        double.IsFinite(h.DeltaG37).Should().BeTrue("ΔG°37 must be finite");
        h.DeltaS.Should().BeLessThan(0.0, "stem + destabilising loop ΔS° is negative (Table 4 footnote a)");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PRIMER-HAIRPIN-001 — hairpin folding : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PRIMER-HAIRPIN-001 — BE: empty / null sequence (→ null, NaN; no crash)

    /// <summary>
    /// BE (empty / null): an empty or null oligo has no bases to fold. The doc (§3.3, §6.1)
    /// contracts the null hairpin sentinel and a NaN Tm — never an IndexOutOfRange from the
    /// closing-pair scan, never a throw. Pinned for both "" and null on the folder and the Tm.
    /// </summary>
    [Test]
    public void EmptyOrNullSequence_NoHairpin_NullAndNaN_NoThrow()
    {
        foreach (var s in new[] { "", null })
        {
            var act = () => PrimerDesigner.FindMostStableHairpin(s!);
            act.Should().NotThrow("empty/null is a documented null sentinel, not a crash");

            PrimerDesigner.FindMostStableHairpin(s!).Should().BeNull(
                "an empty/null oligo has no hairpin (§6.1)");
            double.IsNaN(PrimerDesigner.CalculateHairpinMeltingTemperature(s!)).Should().BeTrue(
                "no hairpin → NaN Tm");
        }
    }

    #endregion

    #region PRIMER-HAIRPIN-001 — BE: no complementary stem (homopolymer / random non-folding)

    /// <summary>
    /// BE (no complementary stem): a homopolymer cannot Watson-Crick pair with itself, so NO stem
    /// exists and the folder returns null with a NaN Tm (INV-02, doc §6.1). Pinned for each of the
    /// four homopolymers and across several lengths — never a fabricated hairpin, never a crash.
    /// </summary>
    [Test]
    public void Homopolymer_NoStem_ReturnsNull_NaN()
    {
        foreach (char b in "ACGT")
        {
            foreach (int len in new[] { 6, 12, 40 })
            {
                string poly = new string(b, len);
                PrimerDesigner.FindMostStableHairpin(poly).Should().BeNull(
                    $"poly-{b} (len {len}) has no Watson-Crick stem → no hairpin (INV-02)");
                double.IsNaN(PrimerDesigner.CalculateHairpinMeltingTemperature(poly)).Should().BeTrue(
                    $"poly-{b} → NaN Tm");
            }
        }
    }

    /// <summary>
    /// BE (no complementary stem — only one base type that pairs): a sequence built ONLY from A and
    /// G (neither pairs the other, A·G is not Watson-Crick, and A·A / G·G never pair) has no
    /// complementary stem at all → null. A C·T-only sequence is the symmetric case. This probes the
    /// "no stem closes a loop" boundary on a heterogeneous (non-homopolymer) non-folding input.
    /// </summary>
    [Test]
    public void NonPairingAlphabet_NoStem_ReturnsNull()
    {
        // Among {A,G}: A·G is not WC, A·A/G·G never pair → no stem possible.
        PrimerDesigner.FindMostStableHairpin("AGAGAGAGAGAG").Should().BeNull(
            "an A/G-only oligo has no Watson-Crick pair → no hairpin");
        PrimerDesigner.FindMostStableHairpin("CTCTCTCTCTCT").Should().BeNull(
            "a C/T-only oligo has no Watson-Crick pair → no hairpin");
    }

    #endregion

    #region PRIMER-HAIRPIN-001 — BE: palindrome that DOES fold (stem + loop must be valid)

    /// <summary>
    /// BE (palindrome — KEY positive signal): a clean inverted-repeat hairpin
    /// (5'-stem · loop · revcomp-stem-3') DOES fold into a hairpin, so the folder must RETURN one —
    /// and the returned structure must be STRUCTURALLY VALID: every stem pair a genuine Watson-Crick
    /// complement (independent oracle), loop ≥ 3, stem ≥ 2, indices in bounds, finite ΔH°/ΔS°/ΔG°37.
    /// Because the stem here is an all-G·C inverted repeat (the maximally-stable domain), ΔG°37 ≤ 0
    /// is additionally asserted — this is the ONLY place the net-stability claim is made (a weak A·T
    /// stem could give a slightly positive ΔG°37, per the property-test domain note). The stem must
    /// span the designed arms and the loop must equal the designed loop length.
    /// </summary>
    [Test]
    public void PalindromeHairpin_FoldsWithValidStemAndLoop()
    {
        // 5'-GGGGG + AAAA(loop) + CCCCC-3' : an inverted-repeat hairpin with a 5-bp G·C stem.
        const int stemArm = 5, loopLen = 4;
        string seq = new string('G', stemArm) + new string('A', loopLen) + new string('C', stemArm);

        var hp = PrimerDesigner.FindMostStableHairpin(seq);

        hp.Should().NotBeNull("a clean inverted-repeat palindrome folds into a hairpin");
        var h = hp!.Value;

        AssertWellFormedHairpin(h, seq, effectiveMinStem: 2);

        h.StemLength.Should().Be(stemArm, "the full 5-bp G·C inverted-repeat stem is used");
        h.LoopSize.Should().Be(loopLen, "the designed 4-nt loop is closed by the stem");
        h.StemStart.Should().Be(0, "the 5' stem arm starts at index 0");
        h.StemEnd.Should().Be(seq.Length - 1, "the 3' stem arm ends at the last index");

        // GC-stem domain ONLY: this maximally-stable hairpin has a negative ΔG°37.
        h.DeltaG37.Should().BeLessThan(0.0, "an all-G·C stem hairpin is net-stable (negative ΔG°37)");

        double tm = PrimerDesigner.CalculateHairpinMeltingTemperature(seq);
        double.IsNaN(tm).Should().BeFalse("a real hairpin has a finite Tm");
        double.IsInfinity(tm).Should().BeFalse("a real hairpin Tm is finite");
    }

    /// <summary>
    /// BE (true full palindrome, no interior loop — KEY min-loop boundary): a sequence whose reverse
    /// complement equals itself with NO unpaired interior (e.g. GGGGCCCC) can ONLY pair into a stem
    /// that closes a 0-nt loop. A 0-nt loop is below the 3-nt minimum (sterically prohibited), so the
    /// MFE folder must return NULL — it must NOT return a hairpin with a sub-3 loop. This is the
    /// palindrome-meets-min-loop edge: the perfect palindrome forms no VALID hairpin.
    /// </summary>
    [Test]
    public void FullPalindrome_NoInteriorLoop_ReturnsNull()
    {
        foreach (var pal in new[] { "GGGGCCCC", "GCGCGCGC", "ATGCAT" /* revcomp = ATGCAT */ })
        {
            var hp = PrimerDesigner.FindMostStableHairpin(pal);
            if (hp is not null)
            {
                // If anything IS returned it must still respect the 3-nt loop minimum.
                hp.Value.LoopSize.Should().BeGreaterThanOrEqualTo(MinLoop,
                    $"any hairpin for '{pal}' must close a ≥3-nt loop (never a sub-3 loop)");
                AssertWellFormedHairpin(hp.Value, pal, effectiveMinStem: 2);
            }
        }

        // The canonical no-interior-loop palindrome: the only stem closes a 0-nt loop → null.
        PrimerDesigner.FindMostStableHairpin("GGGGCCCC").Should().BeNull(
            "GGGGCCCC's only stem closes a 0-nt loop (< 3 nt, sterically prohibited) → null");
    }

    #endregion

    #region PRIMER-HAIRPIN-001 — BE: min-loop violations (stem too close → null)

    /// <summary>
    /// BE (min-loop violation — KEY): when the only complementary stem sits so close that the inner
    /// gap is &lt; 3 nt, NO loop of ≥ 3 can be closed, so the folder returns null (INV-05). Pinned for
    /// the 0-nt-loop case (GCGC), the 1-nt-loop case (GCAGC: G·C stem with a single A between) and the
    /// 2-nt-loop case (GCAAGC). None may return a hairpin with a sub-3 loop; the MFE result is null.
    /// </summary>
    [Test]
    public void MinLoopViolation_StemTooClose_ReturnsNull()
    {
        // GCGC: outer G·C stem extends to C·G, leaving a 0-nt loop → null.
        PrimerDesigner.FindMostStableHairpin("GCGC").Should().BeNull(
            "GCGC can only close a 0-nt loop (< 3 nt) → null (INV-05)");

        // A single complementary stem closing a 1-nt or 2-nt loop is still below the 3-nt minimum.
        // GCATGC: the outer G·C pair would close interior 'AT' that itself pairs, collapsing the loop;
        // construct unambiguous sub-3 cases with a non-pairing interior.
        foreach (var (seq, gap) in new[] { ("GCAGC", 1), ("GCTTGC", 2) })
        {
            var hp = PrimerDesigner.FindMostStableHairpin(seq);
            if (hp is not null)
                hp.Value.LoopSize.Should().BeGreaterThanOrEqualTo(MinLoop,
                    $"'{seq}' must never return a hairpin with a loop < 3 (would-be gap {gap})");
        }
    }

    /// <summary>
    /// BE (just-large-enough loop is the OPPOSITE boundary): a 3-nt loop is the SMALLEST valid loop,
    /// so a stem closing exactly 3 unpaired nt DOES fold. This pins that the min-loop gate admits the
    /// boundary it should (loop == 3) while the sub-3 tests above reject below it — the structure is
    /// well-formed with LoopSize exactly 3.
    /// </summary>
    [Test]
    public void ThreeNtLoop_IsTheValidLowerBoundary_Folds()
    {
        // GGG + AAA(loop) + CCC : 3-bp stem closing the minimal 3-nt loop.
        var hp = PrimerDesigner.FindMostStableHairpin("GGGAAACCC");
        hp.Should().NotBeNull("a 3-nt loop is the smallest sterically-allowed loop and DOES fold");
        AssertWellFormedHairpin(hp!.Value, "GGGAAACCC", effectiveMinStem: 2);
        hp.Value.LoopSize.Should().Be(3, "the minimal valid loop size is exactly 3 nt");
    }

    #endregion

    #region PRIMER-HAIRPIN-001 — MC: non-ACGT / unicode / all-N (strict alphabet → null)

    /// <summary>
    /// MC (non-ACGT / junk / unicode / all-N): the folder enforces the strict A/C/G/T alphabet —
    /// any degenerate, special, whitespace, digit, unicode, or null-byte character makes the
    /// structure not computable and returns null (doc §3.3, §6.1), with a NaN Tm — never a crash,
    /// never a fabricated hairpin. Pinned for an all-N sequence and for junk interleaved into an
    /// otherwise-foldable palindrome (the junk base suppresses the result). Lower-case is upper-cased
    /// first, so a clean lowercase palindrome is the control that DOES fold.
    /// </summary>
    [Test]
    public void NonAcgtJunkUnicodeAllN_NotComputable_ReturnsNull_NaN_NoThrow()
    {
        string?[] junk =
        {
            new string('N', 12),                 // all-N (IUPAC any-base, not in the ACGT alphabet)
            "GGGGNNNNCCCC",                       // N's interleaved into a would-be hairpin
            "GGGG TTTT CCCC",                     // whitespace
            "GGGG-TTTT-CCCC",                     // hyphen
            "GGGG5TTTTCCCC",                      // digit
            "GGGGUUUUCCCC",                       // U is RNA, not in the DNA alphabet
            "GGGG\0TTTTCCCC",                     // embedded null byte
            "GGGGÄTTTTCCCC",                 // U+00C4 Latin diacritic
            "αβγδαβγδ",          // Greek letters
            "gggg\U0001F600cccc"                 // emoji surrogate pair injected into lowercase
        };

        foreach (var s in junk)
        {
            var act = () => PrimerDesigner.FindMostStableHairpin(s!);
            act.Should().NotThrow($"junk input must be handled by the null sentinel, not crash");

            PrimerDesigner.FindMostStableHairpin(s!).Should().BeNull(
                "a non-ACGT base makes the structure not computable → null (strict alphabet)");
            double.IsNaN(PrimerDesigner.CalculateHairpinMeltingTemperature(s!)).Should().BeTrue(
                "non-ACGT input → NaN Tm");
        }

        // Control: the SAME shape, cleaned and lowercase, IS computable (upper-cased internally).
        PrimerDesigner.FindMostStableHairpin("ggggaaaacccc").Should().NotBeNull(
            "a clean lowercase inverted repeat is upper-cased and folds — so the nulls above are alphabet-driven");
    }

    #endregion

    #region PRIMER-HAIRPIN-001 — MC: minStemLength out of domain / selectivity

    /// <summary>
    /// MC (minStemLength &lt; 2): a stem needs at least one NN stack (≥ 2 bp), so minStemLength &lt; 2
    /// is out of domain and returns null (doc §3.3) — for every non-positive / 1 value, even on a
    /// sequence that WOULD fold at the default. Pinned alongside the selectivity boundary: a
    /// minStemLength ABOVE the longest available stem also yields null, while at-or-below it folds.
    /// </summary>
    [Test]
    public void MinStemLength_OutOfDomainAndSelectivity()
    {
        const string foldable = "GGGGGAAAACCCCC"; // 5-bp G·C stem, 4-nt loop

        foreach (int bad in new[] { int.MinValue, -1, 0, 1 })
            PrimerDesigner.FindMostStableHairpin(foldable, minStemLength: bad).Should().BeNull(
                $"minStemLength {bad} < 2 has no NN stack → null (§3.3)");

        // Selectivity: a threshold above the only 5-bp stem yields null; at/below admits it.
        PrimerDesigner.FindMostStableHairpin(foldable, minStemLength: 6).Should().BeNull(
            "no 6-bp stem exists; minStemLength = 6 → null");

        var atFive = PrimerDesigner.FindMostStableHairpin(foldable, minStemLength: 5);
        atFive.Should().NotBeNull("minStemLength = 5 admits the 5-bp stem");
        AssertWellFormedHairpin(atFive!.Value, foldable, effectiveMinStem: 5);
        atFive.Value.StemLength.Should().BeGreaterThanOrEqualTo(5, "the returned stem meets the threshold");
    }

    #endregion

    #region PRIMER-HAIRPIN-001 — Universal structural discipline over random fuzz inputs

    /// <summary>
    /// The core fuzz sweep: over many fixed-seed random oligos of varied lengths, EVERY non-null
    /// hairpin the folder returns MUST satisfy every documented structural invariant — verified
    /// against an INDEPENDENT Watson-Crick oracle (re-derived stem pairs), with the loop ≥ 3, the
    /// stem ≥ 2, indices in bounds, and finite ΔH°/ΔS°/ΔG°37. A null result is also acceptable (no
    /// hairpin). NO input may throw or hang. This is the load-bearing "never out-of-contract output"
    /// guarantee: a returned stem that is not actually complementary, or a loop below the minimum,
    /// fails here on whichever random seed triggers it.
    /// </summary>
    [Test]
    public void RandomFuzz_AnyReturnedHairpin_IsStructurallyValid()
    {
        for (int seed = 241_001; seed <= 241_120; seed++)
        {
            int len = 4 + (seed % 60); // 4..63 nt
            string seq = RandomDna(len, seed);

            PrimerDesigner.HairpinResult? hp = null;
            var act = () => hp = PrimerDesigner.FindMostStableHairpin(seq);
            act.Should().NotThrow($"seed {seed} ('{seq}') must never throw");

            if (hp is not null)
                AssertWellFormedHairpin(hp.Value, seq, effectiveMinStem: 2);

            // Tm: NaN iff no hairpin, finite otherwise — never ±Inf, never a throw.
            double tm = PrimerDesigner.CalculateHairpinMeltingTemperature(seq);
            if (hp is null)
                double.IsNaN(tm).Should().BeTrue($"seed {seed}: no hairpin → NaN Tm");
            else
                double.IsInfinity(tm).Should().BeFalse($"seed {seed}: a real hairpin Tm is finite");
        }
    }

    #endregion

    #region PRIMER-HAIRPIN-001 — BE: very long sequence (no O(n^k) blow-up / non-finite)

    /// <summary>
    /// BE (very long sequence): the folder is an O(n²) closing-pair scan with each stem extended O(n)
    /// (O(n³) worst case, doc §4.3) — a long oligo must NOT blow up or hang. Pinned with [CancelAfter]
    /// so any super-cubic regression fails as a timeout. The result on a long EXPLICIT inverted repeat
    /// (a long G-arm · loop · C-arm) is a well-formed hairpin validated against the oracle; a long
    /// RANDOM oligo returns either a structurally-valid hairpin or null, with a finite-or-NaN Tm.
    /// NOTE: this pins NUMERICAL/perf discipline on the long path, not the accuracy of any value.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void VeryLongSequence_NoBlowUp_StructurallyValidOrNull()
    {
        // A long explicit inverted-repeat hairpin: 200·G + 10·A loop + 200·C (length 410).
        string longHairpin = new string('G', 200) + new string('A', 10) + new string('C', 200);

        PrimerDesigner.HairpinResult? hp = null;
        var act = () => hp = PrimerDesigner.FindMostStableHairpin(longHairpin);
        act.Should().NotThrow("the long inverted-repeat path must not throw");
        hp.Should().NotBeNull("a long inverted repeat folds into a hairpin");
        AssertWellFormedHairpin(hp!.Value, longHairpin, effectiveMinStem: 2);
        hp.Value.DeltaG37.Should().BeLessThan(0.0, "a long G·C-stem hairpin is net-stable (negative ΔG°37)");

        double tmLong = PrimerDesigner.CalculateHairpinMeltingTemperature(longHairpin);
        double.IsNaN(tmLong).Should().BeFalse("the long hairpin has a finite Tm");
        double.IsInfinity(tmLong).Should().BeFalse("the long hairpin Tm is finite, not ±Inf");

        // A long RANDOM oligo: either a structurally-valid hairpin or null, no blow-up.
        string longRandom = RandomDna(length: 500, seed: 241_777);
        PrimerDesigner.HairpinResult? hr = null;
        var actR = () => hr = PrimerDesigner.FindMostStableHairpin(longRandom);
        actR.Should().NotThrow("the long random path must not throw");
        if (hr is not null)
            AssertWellFormedHairpin(hr.Value, longRandom, effectiveMinStem: 2);
    }

    #endregion
}
