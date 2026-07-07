using System.Text;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Translation area, unit TRANS-SIXFRAME-001 — six-frame
/// translation of a double-stranded DNA sequence, exposed by
/// <see cref="Translator.TranslateSixFrames(DnaSequence, GeneticCode?)"/> in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Core/Translator.cs.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, invalid and boundary inputs to a unit and asserts that
/// the code NEVER fails in an undisciplined way: no hang, no NaN/Infinity, no
/// state corruption, and no *unhandled* runtime exception
/// (KeyNotFoundException, NullReferenceException, IndexOutOfRangeException, …).
/// Every input must produce EITHER a well-defined, theory-correct value OR a
/// *documented, intentional* validation exception (ArgumentNullException /
/// ArgumentException). A raw runtime exception or a hang on garbage input is a
/// bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: TRANS-SIXFRAME-001 — six-frame translation (Translation)
/// Checklist: docs/checklists/03_FUZZING.md, row 223.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — length NOT a multiple of 3 (trailing partial
///     codon dropped per frame), the EMPTY sequence (defined result: six empty
///     frames, no crash), and the SINGLE-BASE sequence (too short for any codon
///     ⇒ all six frames empty).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes); row 223.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The six-frame contract under test (from Six_Frame_Translation.md)
/// ───────────────────────────────────────────────────────────────────────────
/// A duplex sequence has exactly six reading frames: three forward frames at
/// offsets 0,1,2 of the given strand, and three reverse frames at offsets 0,1,2
/// of its reverse complement (Six_Frame_Translation.md §2.1). The exact
/// behaviour asserted here is:
///   • EXACTLY SIX frames keyed {+1,+2,+3,−1,−2,−3}, no key 0
///     (INV-01; §2.4, §3.2).
///   • Forward frame +f = translation at offset f−1 of the INPUT strand
///     (INV-02; §2.2).
///   • Reverse frame −f = translation at offset f−1 of the REVERSE COMPLEMENT
///     (Biopython frames[-(i+1)] = translate(anti[i:]); INV-03; §2.2, §5.4).
///   • Each frame length = ⌊(len − offset)/3⌋, with the trailing partial codon
///     ignored (INV-04; §6.1 "Length not multiple of 3").
///   • Empty sequence → six EMPTY frames, no ORFs, no crash (§3.3, §6.1).
///   • Single base (and length 2) → too short for any codon ⇒ all six frames
///     empty (INV-04 with len < 3 ⇒ ⌊(len−offset)/3⌋ = 0; §6.1).
///   • Internal stop codons are rendered as '*' (translation NOT terminated
///     early), unlike Translate(..., toFirstStop: true) (§5.2).
///   • Null dna → ArgumentNullException (intentional; §3.3, §6.1) — never a bare
///     NullReferenceException.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Independence of expected values (anti-tautology discipline)
/// ───────────────────────────────────────────────────────────────────────────
/// Expected frame strings are NEVER read back from the implementation's arrays.
/// They are derived from PRIMARY sources independently of Translator.cs:
///   • <see cref="StandardCodonTable"/> is the canonical NCBI Standard Code
///     (transl_table=1) built from the textbook TCAG-ordered amino-acid string
///     [NCBI; Six_Frame_Translation.md §4.2 ref [3]] — NOT GeneticCode's map.
///   • <see cref="ReverseComplement"/> recomputes the antiparallel strand
///     independently (A↔T, C↔G, reversed).
///   • <see cref="ReferenceFrame"/> reproduces the offset codon-reading model of
///     §2.2 (consume only complete codons), so a wrong implementation would FAIL
///     these tests rather than silently agree with itself.
/// The hand-checked worked example from §7.1 is pinned verbatim and each of its
/// six frames is also re-derived from the canonical table below.
///
/// Determinism note: each randomized test uses a LOCALLY-seeded
/// new Random(seed); there is NO shared static Rng. Hang-sensitive / sweep
/// tests carry [CancelAfter] as a tripwire.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class TranslationSixFrameFuzzTests
{
    // ═══════════════════════════════════════════════════════════════════════
    //  Independent reference oracle — built from the PRIMARY genetic-code table,
    //  not from GeneticCode / Translator internals.
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>The six frame keys the contract requires (INV-01): +1..+3, −1..−3.</summary>
    private static readonly int[] FrameKeys = { 1, 2, 3, -1, -2, -3 };

    /// <summary>
    /// Canonical NCBI Standard Code (transl_table=1). The amino-acid string is
    /// the textbook one-letter ordering for codons enumerated with the bases in
    /// the order T,C,A,G for each position (NCBI "The Genetic Codes";
    /// Six_Frame_Translation.md §4.2 ref [3]). Built here independently so the
    /// expected values cannot be a copy of the code under test.
    /// </summary>
    private static readonly IReadOnlyDictionary<string, char> StandardCodonTable = BuildStandardTable();

    private static Dictionary<string, char> BuildStandardTable()
    {
        const string bases = "TCAG";
        // Amino acids for the 64 codons in TCAG×TCAG×TCAG order (NCBI table 1).
        const string aminoAcids =
            "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

        var table = new Dictionary<string, char>(64);
        int i = 0;
        foreach (var b1 in bases)
            foreach (var b2 in bases)
                foreach (var b3 in bases)
                    table[$"{b1}{b2}{b3}"] = aminoAcids[i++];
        return table;
    }

    /// <summary>
    /// Independent reverse complement over the ACGT alphabet (the only alphabet
    /// <see cref="DnaSequence"/> accepts): complement each base (A↔T, C↔G) and
    /// reverse. Mirrors the antiparallel-strand definition (§2.1) without calling
    /// the production ReverseComplement.
    /// </summary>
    private static string ReverseComplement(string dna)
    {
        var chars = new char[dna.Length];
        for (int i = 0; i < dna.Length; i++)
        {
            char c = dna[dna.Length - 1 - i];
            chars[i] = c switch
            {
                'A' => 'T',
                'T' => 'A',
                'C' => 'G',
                'G' => 'C',
                _ => throw new InvalidOperationException($"unexpected base '{c}'"),
            };
        }
        return new string(chars);
    }

    /// <summary>
    /// Independent single-frame translation reproducing the §2.2 model: read
    /// codons at the given offset, consume only COMPLETE triplets, map each via
    /// the canonical Standard table, mark stops as '*'. Inputs are DNA (ACGT);
    /// the table keys are DNA codons, so no T→U step is needed here.
    /// </summary>
    private static string ReferenceFrame(string dna, int offset)
    {
        var sb = new StringBuilder();
        for (int i = offset; i + 3 <= dna.Length; i += 3)
            sb.Append(StandardCodonTable[dna.Substring(i, 3)]);
        return sb.ToString();
    }

    /// <summary>
    /// Independent reference for ALL six frames, keyed exactly as the contract
    /// requires: +f = forward offset f−1, −f = reverse-complement offset f−1.
    /// </summary>
    private static Dictionary<int, string> ReferenceSixFrames(string dna)
    {
        string rc = ReverseComplement(dna);
        var expected = new Dictionary<int, string>();
        for (int offset = 0; offset < 3; offset++)
        {
            expected[offset + 1] = ReferenceFrame(dna, offset);
            expected[-(offset + 1)] = ReferenceFrame(rc, offset);
        }
        return expected;
    }

    /// <summary>Materialises the production six-frame result as plain strings.</summary>
    private static Dictionary<int, string> ActualSixFrames(string dna)
    {
        var frames = Translator.TranslateSixFrames(new DnaSequence(dna));
        return frames.ToDictionary(kv => kv.Key, kv => kv.Value.Sequence);
    }

    // ═══════════════════════════════════════════════════════════════════════
    //  TRANS-SIXFRAME-001 — six-frame translation : fuzz targets
    // ═══════════════════════════════════════════════════════════════════════

    #region TRANS-SIXFRAME-001 — six-frame translation

    #region Positive sanity — the hand-checked worked example (all six frames)

    /// <summary>
    /// Positive sanity (worked example, §7.1): for the documented sequence
    /// "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG" the doc pins
    /// frames[+1] == "MAIVMGR*KGAR*" and frames[−1] == "LSGTLSAAHYNGH". Both are
    /// re-derived here from the canonical NCBI Standard table and the independent
    /// reverse complement (NOT from the implementation), proving the fuzz targets
    /// below are measured against a functioning, theory-correct translator.
    /// The remaining four frames are pinned against the same independent oracle.
    /// — Six_Frame_Translation.md §7.1, §2.2; INV-02, INV-03.
    /// </summary>
    [Test]
    public void SixFrames_WorkedExample_MatchesDocAndIndependentReference()
    {
        const string seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"; // 39 bases

        var actual = ActualSixFrames(seq);

        // Verbatim from the algorithm doc §7.1 (independently re-derived below).
        actual[1].Should().Be("MAIVMGR*KGAR*", "doc §7.1 pins frame +1");
        actual[-1].Should().Be("LSGTLSAAHYNGH", "doc §7.1 pins frame -1");

        // Independent oracle agrees with the doc on +1/-1 (proves the oracle is
        // correct), then pins every frame against that same oracle.
        var expected = ReferenceSixFrames(seq);
        expected[1].Should().Be("MAIVMGR*KGAR*");
        expected[-1].Should().Be("LSGTLSAAHYNGH");

        foreach (var key in FrameKeys)
            actual[key].Should().Be(expected[key],
                "frame {0} must equal the independent canonical-table translation", key);

        // Internal stop codons are rendered as '*', not truncated (§5.2).
        actual[1].Should().Contain("*");
    }

    #endregion

    #region Contract — exactly six frames, keyed {+1,+2,+3,-1,-2,-3}, no key 0

    /// <summary>
    /// Contract INV-01: <c>TranslateSixFrames</c> returns EXACTLY six entries
    /// keyed {+1,+2,+3,−1,−2,−3}; key 0 never appears and no extra keys leak in.
    /// Asserted across a spread of representative lengths (including non-%3 and
    /// boundary lengths) so the key set is independent of input length.
    /// — Six_Frame_Translation.md §3.2, §2.4 INV-01.
    /// </summary>
    [Test]
    public void SixFrames_AlwaysExactlySixKeys_NoKeyZero()
    {
        foreach (var seq in new[] { "", "A", "AT", "ACG", "ACGT", "ATGGCCTGA", "ATGGCCTGAA" })
        {
            var frames = Translator.TranslateSixFrames(new DnaSequence(seq));

            frames.Should().HaveCount(6, "a duplex always has six reading frames (INV-01)");
            frames.Keys.Should().BeEquivalentTo(FrameKeys);
            frames.ContainsKey(0).Should().BeFalse("frame 0 is not a reading frame");
        }
    }

    #endregion

    #region BE: empty input → six empty frames, no crash

    /// <summary>
    /// Fuzz target "empty" (BE): an empty sequence is a DEFINED result, not a
    /// crash — it yields six frames, every one of which is the empty protein
    /// (no complete codons exist). Both <c>new DnaSequence("")</c> and the
    /// implicit-empty case are exercised.
    /// — Six_Frame_Translation.md §3.3, §6.1 "Empty sequence".
    /// </summary>
    [Test]
    public void SixFrames_EmptySequence_YieldsSixEmptyFrames_NoCrash()
    {
        Action act = () => Translator.TranslateSixFrames(new DnaSequence(""));
        act.Should().NotThrow("an empty sequence is a defined input, not an error");

        var actual = ActualSixFrames("");
        actual.Should().HaveCount(6);
        foreach (var key in FrameKeys)
            actual[key].Should().BeEmpty("frame {0} of an empty sequence has no codons", key);
    }

    #endregion

    #region BE: single base (and length 2) → no codon in any frame → all six empty

    /// <summary>
    /// Fuzz target "single base" (BE): a one-base sequence is too short for any
    /// codon in ANY of the six frames, so every frame is empty — and crucially
    /// there is NO IndexOutOfRangeException from attempting to read a partial
    /// triplet. Length-2 inputs are swept too (still &lt; 3, still all empty),
    /// pinning ⌊(len−offset)/3⌋ = 0 for len &lt; 3 (INV-04).
    /// — Six_Frame_Translation.md §6.1; INV-04.
    /// </summary>
    [Test]
    public void SixFrames_SingleBaseOrLengthTwo_AllFramesEmpty_NeverIndexOutOfRange()
    {
        foreach (var single in new[] { "A", "C", "G", "T" })
        {
            Action act = () => Translator.TranslateSixFrames(new DnaSequence(single));
            act.Should().NotThrow("single base \"{0}\" has no complete codon", single);

            var actual = ActualSixFrames(single);
            actual.Should().HaveCount(6);
            foreach (var key in FrameKeys)
                actual[key].Should().BeEmpty(
                    "single-base input \"{0}\" yields no codon in frame {1}", single, key);
        }

        foreach (var pair in new[] { "AT", "GC", "TA", "CG", "AA", "GG" })
        {
            var actual = ActualSixFrames(pair);
            foreach (var key in FrameKeys)
                actual[key].Should().BeEmpty(
                    "length-2 input \"{0}\" still has no complete codon (frame {1})", pair, key);
        }
    }

    #endregion

    #region BE: length not a multiple of 3 → trailing partial codon dropped per frame

    /// <summary>
    /// Fuzz target "length not %3" (BE): when the input length is not a multiple
    /// of three, each frame consumes only COMPLETE codons and the trailing 1–2
    /// dangling bases are silently dropped, so every frame length is exactly
    /// ⌊(len − offset)/3⌋ (INV-04) — never an IndexOutOfRange on the partial tail.
    /// "ATGGCC" (a clean 2-codon prefix) is extended by 1 and 2 dangling bases;
    /// forward frame +1 must stay "MA" in all three cases, and every frame's
    /// length must match the floor formula and the independent oracle.
    /// — Six_Frame_Translation.md §6.1 "Length not multiple of 3"; INV-04.
    /// </summary>
    [Test]
    public void SixFrames_LengthNotMultipleOfThree_DropsTrailingPartialCodon()
    {
        foreach (var dangling in new[] { "", "A", "AT" }) // len 6, 7, 8
        {
            string seq = "ATGGCC" + dangling;
            var actual = ActualSixFrames(seq);
            var expected = ReferenceSixFrames(seq);

            // Forward frame +1 reads ATG GCC = "MA" regardless of the dropped tail.
            actual[1].Should().Be("MA",
                "trailing partial codon \"{0}\" must not change forward frame +1", dangling);

            foreach (var key in FrameKeys)
            {
                int offset = Math.Abs(key) - 1;
                int expectedLen = (seq.Length - offset) / 3; // ⌊(len−offset)/3⌋

                actual[key].Length.Should().Be(expectedLen,
                    "frame {0} length must be floor((len-offset)/3) (INV-04)", key);
                actual[key].Should().Be(expected[key],
                    "frame {0} must equal the independent canonical translation", key);
            }
        }
    }

    #endregion

    #region BE: randomized boundary sweep — no crash/hang AND full contract pinned

    /// <summary>
    /// Randomized boundary sweep (BE): a locally-seeded Random generates ACGT
    /// sequences whose lengths cluster around the codon boundaries (0,1,2 and
    /// around multiples of 3, plus a few large lengths). For EVERY generated
    /// input the full contract must hold simultaneously: no crash/hang, exactly
    /// six keys {+1..+3,−1..−3}, each frame equal to the independent
    /// canonical-table translation, and each frame length = ⌊(len−offset)/3⌋.
    /// A wrong frame ordering, a swapped strand, a wrong offset, or a mishandled
    /// partial codon would all FAIL here. [CancelAfter] is a hang tripwire.
    /// — Six_Frame_Translation.md §2.2, §6.1; INV-01..INV-04.
    /// </summary>
    [Test]
    [CancelAfter(15000)]
    public void SixFrames_RandomBoundaryLengths_PinFullContract_NoCrashNoHang(
        CancellationToken token)
    {
        const string alphabet = "ACGT";
        var rng = new Random(20260621);

        // Lengths concentrated on the BE-relevant region: 0,1,2, around %3, plus big.
        var lengths = new List<int> { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
        for (int k = 0; k < 40; k++)
            lengths.Add(rng.Next(0, 200));
        lengths.Add(1000);
        lengths.Add(1001);
        lengths.Add(1002);

        foreach (var len in lengths)
        {
            token.ThrowIfCancellationRequested();

            var sb = new StringBuilder(len);
            for (int i = 0; i < len; i++)
                sb.Append(alphabet[rng.Next(alphabet.Length)]);
            string seq = sb.ToString();

            Dictionary<int, string> actual = null!;
            Action act = () => actual = ActualSixFrames(seq);
            act.Should().NotThrow("six-frame translation must never crash (len {0})", len);

            actual.Should().HaveCount(6, "INV-01: exactly six frames (len {0})", len);
            actual.Keys.Should().BeEquivalentTo(FrameKeys);

            var expected = ReferenceSixFrames(seq);
            foreach (var key in FrameKeys)
            {
                int offset = Math.Abs(key) - 1;
                int expectedLen = Math.Max(0, (seq.Length - offset)) / 3;

                actual[key].Should().Be(expected[key],
                    "frame {0} (len {1}) must equal the independent translation", key, len);
                actual[key].Length.Should().Be(expectedLen,
                    "frame {0} length = floor((len-offset)/3) (INV-04, len {1})", key, len);
                // No corrupt residues: every char is a valid AA letter or stop '*'.
                actual[key].Should().MatchRegex("^[A-Z*]*$",
                    "frame {0} must contain only amino-acid letters or '*'", key);
            }
        }
    }

    #endregion

    #region Null input → ArgumentNullException, never NullReference

    /// <summary>
    /// Boundary "null dna" (BE): a null sequence is rejected with the documented,
    /// intentional ArgumentNullException — NOT a bare NullReferenceException.
    /// — Six_Frame_Translation.md §3.3, §6.1 "Null dna".
    /// </summary>
    [Test]
    public void SixFrames_NullDna_ThrowsArgumentNullException_NotNullReference()
    {
        Action act = () => Translator.TranslateSixFrames(null!);
        act.Should().Throw<ArgumentNullException>("null dna is rejected intentionally")
            .And.ParamName.Should().Be("dna");
    }

    #endregion

    #endregion
}
