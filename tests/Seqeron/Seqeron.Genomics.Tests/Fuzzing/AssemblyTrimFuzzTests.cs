using System.Text;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Assembly area — Quality Trimming (ASSEMBLY-TRIM-001), the
/// BWA / cutadapt running-sum quality trimmer
/// <see cref="SequenceAssembler.QualityTrimReads(IReadOnlyList{ValueTuple{string,string}}, int, int)"/>.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and
/// asserts the code NEVER fails in an undisciplined way: no hang or infinite loop
/// (the two linear end-passes per read must always terminate — Quality_Trimming.md
/// §4.3 O(n·r)), no state corruption, no nonsense output (a trimmed read that is
/// NOT a contiguous substring of its input, a trimmed length &gt; the original, a
/// negative-length substring when the whole read is consumed, a read/quality length
/// desync, or a non-deterministic result), and no *unhandled* runtime exception —
/// in particular NO IndexOutOfRangeException / ArgumentOutOfRangeException from a
/// negative-length <c>Substring</c> on a fully-trimmed read and NO DivideByZero /
/// NullReference on the EMPTY boundaries of this row. Every input must resolve to
/// EITHER a well-defined, theory-correct result OR a *documented, intentional*
/// validation exception (ArgumentNullException for null <c>reads</c> — §3.3, §6.1).
/// A raw runtime exception, a hang, a non-substring result, a length-growing trim,
/// or an order-dependent result is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ASSEMBLY-TRIM-001 — Quality Trimming (BWA / cutadapt running-sum)
/// Checklist: docs/checklists/03_FUZZING.md, row 148.
/// Algorithm doc: docs/algorithms/Assembly/Quality_Trimming.md
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row "all-low-quality, all-high-quality, empty, quality cutoff 0":
///          – ALL-LOW-QUALITY: every Phred score below the cutoff → both end passes
///            consume the read → the good-quality window is empty (cutadapt
///            start>=stop ⇒ (0,0)) → trimmed length 0; with minLength ≥ 1 the read
///            is DROPPED. No negative-length Substring, no IndexOutOfRange (§6.1
///            "all-low ⇒ dropped (length 0)").
///          – ALL-HIGH-QUALITY: every Phred score ≥ the cutoff → partial-sum minimum
///            sits at the read end → read returned UNCHANGED, full length (§6.1
///            "all-high ⇒ unchanged"; INV-02).
///          – EMPTY: empty <c>reads</c> → empty result; an empty read (length-0
///            sequence + length-0 quality) → length 0, dropped, with NO crash, NO
///            DivideByZero, NO negative-length Substring (§6.1).
///          – QUALITY CUTOFF 0: cutoff &lt; 1 disables trimming entirely — every base
///            "passes" and the read is returned unchanged (subject only to the
///            min-length filter). This pins the INCLUSIVE/EXCLUSIVE boundary at
///            cutoff 0: cutoff 0 (and any negative cutoff) is the BWA <c>trim_qual
///            &lt; 1</c> guard, NOT a per-base "q ≥ 0" comparison (§4.1 step 1,
///            INV-03, §6.1 "cutoff ≤ 0 ⇒ unchanged").
/// — docs/checklists/03_FUZZING.md §Description (BE = Boundary Exploitation:
///   граничні значення 0, -1, MaxInt, empty).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test (Quality_Trimming.md §2.2, §4.1, §3, §6.1)
/// ───────────────────────────────────────────────────────────────────────────
/// For a read with Phred+33 qualities and integer cutoff C:
///   • if C &lt; 1, trimming is disabled — keep the full read (INV-03, §4.1 step 1);
///   • else compute s_i = (ASCII−33) − C and run TWO independent passes over the
///     FULL read (cutadapt quality_trim_index):
///       3' pass — accumulate s_i from the 3' end; STOP as soon as the running sum
///                 becomes positive (the BWA/cutadapt s&lt;0 early break, sign-flipped);
///                 set `end` to the index of the minimal partial sum reached;
///       5' pass — same from the 5' end; set `start` to one past the minimal index;
///   • if start ≥ end the good-quality segment is empty (cutadapt start>=stop ⇒ (0,0)) —
///     drop the read;
///   • else emit sequence[start..end) iff (end − start) ≥ minLength, else drop it.
/// Output: surviving trimmed sequences, in input order (INV-01 substring of input,
/// INV-02 length ≤ original). Phred decoded as ASCII−33 (INV-05). Null reads →
/// ArgumentNullException (§3.3). Defaults: minQuality = 20, minLength = 50.
///   SequenceAssembler.QualityTrimReads(
///       IReadOnlyList&lt;(string sequence, string quality)&gt; reads,
///       int minQuality = 20, int minLength = 50) → IReadOnlyList&lt;string&gt;
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class AssemblyTrimFuzzTests
{
    // Sanger/Phred+33 ASCII offset (Quality_Trimming.md §2.1, INV-05).
    private const int PhredAsciiOffset = 33;

    private static readonly char[] DnaUpper = { 'A', 'C', 'G', 'T' };

    #region Helpers

    /// <summary>
    /// Asserts every surviving trimmed read is WELL-FORMED with respect to the
    /// ORIGINAL input set, per the documented contract:
    ///   §3.2  the result is never null and never longer than the input set;
    ///   INV-04 every survivor has length ≥ minLength (none of them is empty when
    ///          minLength ≥ 1, so no negative/empty substring leaked through);
    ///   INV-01 every survivor is a CONTIGUOUS SUBSTRING of some input sequence
    ///          (the trim is a single [start,end) cut per end — read/quality stay in
    ///          sync because the same [start,end) span indexes both);
    ///   INV-02 every survivor's length ≤ the matching input sequence's length.
    /// The substring/length check is the read/quality length-sync guard: the method
    /// returns sequence.Substring(start, end−start), and a desync or negative length
    /// would surface here (or as an exception) rather than as a silent pass.
    /// </summary>
    private static void AssertWellFormed(
        IReadOnlyList<string> survivors,
        IReadOnlyList<(string sequence, string quality)> input,
        int minLength)
    {
        survivors.Should().NotBeNull("the trimmer never returns null (§3.2)");
        survivors.Count.Should().BeLessThanOrEqualTo(input.Count,
            "trimming only drops reads, never invents them (§3.2)");

        // Multiset of all contiguous substrings available in the input, so that each
        // survivor can be matched to an input source span.
        var inputSequences = input.Select(r => r.sequence ?? string.Empty).ToList();

        foreach (string survivor in survivors)
        {
            survivor.Should().NotBeNull("INV-01: each survivor is a real string");
            survivor.Length.Should().BeGreaterThanOrEqualTo(Math.Max(0, minLength) == 0 ? 0 : minLength,
                "INV-04: every survivor has length ≥ minLength (no empty/negative-length leak)");

            bool isSubstringOfSome = inputSequences.Any(seq =>
                survivor.Length <= seq.Length && seq.Contains(survivor, StringComparison.Ordinal));
            isSubstringOfSome.Should().BeTrue(
                "INV-01/INV-02: every trimmed read is a contiguous substring of an input sequence " +
                "no longer than the original — read and quality were cut with one [start,end) span");
        }
    }

    /// <summary>
    /// Independent reference implementation of the documented cutadapt
    /// <c>quality_trim_index</c> running-sum trim (Quality_Trimming.md §2.2/§4.1):
    /// two independent full-read passes with the <c>s &lt; 0</c> early break and the
    /// <c>start &gt;= stop ⇒ (0,0)</c> drop rule. Returns the (start, end) span the
    /// production code MUST agree with for one read. This is a DIFFERENT formulation
    /// than the unit under test (kept here only to fuzz-oracle the production code).
    /// </summary>
    private static (int start, int end) ReferenceTrim(string quality, int cutoff)
    {
        int n = quality.Length;
        if (cutoff < 1)
        {
            return (0, n); // BWA trim_qual < 1 guard: no trimming.
        }

        // 3' end pass → end (exclusive upper bound).
        int sum = 0, min = 0, end = n;
        for (int i = n - 1; i >= 0; i--)
        {
            sum += (quality[i] - PhredAsciiOffset) - cutoff;
            if (sum > 0) break;
            if (sum < min) { min = sum; end = i; }
        }

        // 5' end pass → start (inclusive lower bound).
        sum = 0; min = 0; int start = 0;
        for (int i = 0; i < n; i++)
        {
            sum += (quality[i] - PhredAsciiOffset) - cutoff;
            if (sum > 0) break;
            if (sum < min) { min = sum; start = i + 1; }
        }

        if (start >= end) return (0, 0); // cutadapt empty good-quality segment.
        return (start, end);
    }

    /// <summary>The surviving trimmed sequence for ONE read per the reference oracle, or null if dropped.</summary>
    private static string? ReferenceSurvivor(string sequence, string quality, int cutoff, int minLength)
    {
        var (start, end) = ReferenceTrim(quality, cutoff);
        int len = end - start;
        return len >= minLength ? sequence.Substring(start, len) : null;
    }

    /// <summary>A random DNA sequence of the given length over upper-case ACGT.</summary>
    private static string RandomDna(Random rng, int length)
    {
        var sb = new StringBuilder(length);
        for (int i = 0; i < length; i++)
            sb.Append(DnaUpper[rng.Next(DnaUpper.Length)]);
        return sb.ToString();
    }

    /// <summary>A random Phred+33 quality string of the given length spanning a chosen Phred range.</summary>
    private static string RandomQuality(Random rng, int length, int minPhred, int maxPhred)
    {
        var sb = new StringBuilder(length);
        for (int i = 0; i < length; i++)
            sb.Append((char)(PhredAsciiOffset + rng.Next(minPhred, maxPhred + 1)));
        return sb.ToString();
    }

    /// <summary>A Phred+33 quality string in which every base has the same Phred score.</summary>
    private static string ConstantQuality(int length, int phred)
        => new string((char)(PhredAsciiOffset + phred), length);

    #endregion

    #region ASSEMBLY-TRIM-001 — Quality Trimming (BE: all-low, all-high, empty, cutoff 0)

    #region Positive sanity — documented core span is trimmed correctly, ends in sync

    // Documented cutadapt worked example (§7.1): qualities 42,40,26,27,8,7,11,4,2,3,
    // cutoff 10 → keep the FIRST FOUR bases. The high-quality core "ACGT" is retained
    // with correct start (0) and end (4); read and quality were cut with one span.
    [Test]
    public void QualityTrimReads_HighQualityCore_TrimmedToDocumentedSpan()
    {
        // "KI;<)(,%#$" decodes to Phred 42,40,26,27,8,7,11,4,2,3 (Quality_Trimming.md §7.1).
        var reads = new List<(string, string)> { ("ACGTACGTAC", "KI;<)(,%#$") };

        IReadOnlyList<string> result = SequenceAssembler.QualityTrimReads(reads, minQuality: 10, minLength: 1);

        result.Should().ContainSingle().Which.Should().Be("ACGT",
            "the running-sum minimum is at index 4 → keep the first four bases (§7.1)");
        AssertWellFormed(result, reads, 1);
    }

    // Low-quality ends, high-quality core: both ends are trimmed to the exact core span,
    // start and end correct, read/quality in sync (positive both-end sanity).
    [Test]
    public void QualityTrimReads_LowQualityBothEnds_TrimmedToCore()
    {
        // Phred: 2,2,2, 40,40,40,40, 2,2  ("###IIII##") over "AAACCCCGG".
        var reads = new List<(string, string)> { ("AAACCCCGG", "###IIII##") };

        IReadOnlyList<string> result = SequenceAssembler.QualityTrimReads(reads, minQuality: 20, minLength: 1);

        result.Should().ContainSingle().Which.Should().Be("CCCC",
            "both low-quality tails are removed, leaving the high-quality core in sync");
        AssertWellFormed(result, reads, 1);
    }

    #endregion

    #region BE — Boundary: all-high-quality (every base ≥ cutoff ⇒ unchanged, full length)

    // Every base well above the cutoff: the partial-sum minimum sits at the read end →
    // read returned UNCHANGED at full length (§6.1 "all-high ⇒ unchanged", INV-02).
    [Test]
    public void QualityTrimReads_AllHighQuality_ReturnsUnchanged()
    {
        var reads = new List<(string, string)> { ("ACGTACGTAC", ConstantQuality(10, 40)) };

        IReadOnlyList<string> result = SequenceAssembler.QualityTrimReads(reads, minQuality: 20, minLength: 1);

        result.Should().ContainSingle().Which.Should().Be("ACGTACGTAC",
            "all bases ≥ cutoff ⇒ nothing trimmed, full length retained (§6.1, INV-02)");
        AssertWellFormed(result, reads, 1);
    }

    // Phred EXACTLY at the cutoff (q == C ⇒ s_i = 0): a zero contributes nothing, never
    // makes the running sum positive, and never produces a new strict minimum → the read
    // is kept unchanged. This pins the inclusive boundary q ≥ C (not q > C).
    [Test]
    public void QualityTrimReads_QualityExactlyAtCutoff_NotTrimmed()
    {
        var reads = new List<(string, string)> { ("ACGTAC", ConstantQuality(6, 20)) };

        IReadOnlyList<string> result = SequenceAssembler.QualityTrimReads(reads, minQuality: 20, minLength: 1);

        result.Should().ContainSingle().Which.Should().Be("ACGTAC",
            "q == cutoff ⇒ s_i = 0 ⇒ no strict minimum ⇒ kept (inclusive boundary q ≥ C)");
        AssertWellFormed(result, reads, 1);
    }

    // Fuzz: random reads whose every base is strictly above the cutoff are ALWAYS returned
    // unchanged at full length, for many cutoffs/lengths — never over-trimmed, never thrown.
    [Test]
    [CancelAfter(30_000)]
    public void QualityTrimReads_AllHighQuality_AlwaysUnchanged()
    {
        var rng = new Random(148_001);
        for (int trial = 0; trial < 800; trial++)
        {
            int cutoff = rng.Next(1, 41);
            int len = rng.Next(1, 60);
            // Every Phred strictly above the cutoff (and within the valid 0..93 range).
            int maxPhred = Math.Min(93, cutoff + 1 + rng.Next(0, 30));
            string seq = RandomDna(rng, len);
            string qual = RandomQuality(rng, len, cutoff + 1, maxPhred);
            var reads = new List<(string, string)> { (seq, qual) };

            IReadOnlyList<string> result = SequenceAssembler.QualityTrimReads(reads, cutoff, minLength: 1);

            result.Should().ContainSingle().Which.Should().Be(seq,
                "all bases > cutoff ⇒ unchanged at full length (§6.1, INV-02)");
            AssertWellFormed(result, reads, 1);
        }
    }

    #endregion

    #region BE — Boundary: all-low-quality (every base < cutoff ⇒ window empty, dropped, no IOoR)

    // Every base below the cutoff: both end passes consume the read, start >= stop ⇒ (0,0),
    // trimmed length 0 < minLength ⇒ DROPPED. No negative-length Substring, no IndexOutOfRange
    // (§6.1 "all-low ⇒ dropped (length 0)").
    [Test]
    public void QualityTrimReads_AllLowQuality_DroppedNoNegativeSubstring()
    {
        var reads = new List<(string, string)> { ("ACGTAC", ConstantQuality(6, 2)) };

        IReadOnlyList<string> result = SequenceAssembler.QualityTrimReads(reads, minQuality: 20, minLength: 1);

        result.Should().BeEmpty(
            "every base < cutoff ⇒ empty good-quality window ⇒ dropped, no negative-length Substring (§6.1)");
    }

    // Fuzz: random reads whose every base is strictly below the cutoff must ALWAYS be dropped
    // (with minLength ≥ 1) and must NEVER throw a negative-length Substring / IndexOutOfRange.
    [Test]
    [CancelAfter(30_000)]
    public void QualityTrimReads_AllLowQuality_AlwaysDropped_NeverThrows()
    {
        var rng = new Random(148_002);
        for (int trial = 0; trial < 800; trial++)
        {
            int cutoff = rng.Next(1, 41);
            int len = rng.Next(1, 60);
            // Every Phred strictly below the cutoff (and ≥ 0).
            int minPhred = 0;
            int maxPhred = cutoff - 1;
            string seq = RandomDna(rng, len);
            string qual = RandomQuality(rng, len, minPhred, maxPhred);
            var reads = new List<(string, string)> { (seq, qual) };

            IReadOnlyList<string> result = SequenceAssembler.QualityTrimReads(reads, cutoff, minLength: 1);

            result.Should().BeEmpty(
                "all bases < cutoff ⇒ empty window ⇒ dropped; no IndexOutOfRange on full trim (§6.1)");
        }
    }

    #endregion

    #region BE — Boundary: quality cutoff 0 / negative (BWA trim_qual < 1 guard ⇒ unchanged)

    // Cutoff 0: the BWA trim_qual < 1 guard disables trimming entirely. Even all-low-quality
    // bases pass through UNCHANGED — this is NOT a per-base "q ≥ 0" comparison. Pins the
    // cutoff-0 boundary (§4.1 step 1, INV-03, §6.1 "cutoff ≤ 0 ⇒ unchanged").
    [Test]
    public void QualityTrimReads_CutoffZero_DisablesTrimming_ReturnsUnchanged()
    {
        // All Phred 0 ("!!!!!!"): with any cutoff ≥ 1 this would be fully trimmed; cutoff 0 keeps it.
        var reads = new List<(string, string)> { ("ACGTAC", ConstantQuality(6, 0)) };

        IReadOnlyList<string> result = SequenceAssembler.QualityTrimReads(reads, minQuality: 0, minLength: 1);

        result.Should().ContainSingle().Which.Should().Be("ACGTAC",
            "cutoff 0 < 1 disables trimming entirely — even Phred-0 bases pass (§4.1 step 1, INV-03)");
        AssertWellFormed(result, reads, 1);
    }

    // Negative cutoff is also below the trim_qual < 1 guard ⇒ trimming disabled. Confirms the
    // boundary is "cutoff < 1" not "cutoff ≤ 0" alone, and guards against off-by-one at 0.
    [Test]
    public void QualityTrimReads_NegativeCutoff_DisablesTrimming()
    {
        var reads = new List<(string, string)> { ("ACGTAC", ConstantQuality(6, 0)) };

        foreach (int cutoff in new[] { -1, -1000, int.MinValue })
        {
            IReadOnlyList<string> result = SequenceAssembler.QualityTrimReads(reads, cutoff, minLength: 1);

            result.Should().ContainSingle().Which.Should().Be("ACGTAC",
                "cutoff {0} < 1 ⇒ trimming disabled (BWA trim_qual < 1 guard)", cutoff);
        }
    }

    // Cutoff 1 (the first ENABLED cutoff): trimming IS active — an all-Phred-0 read is now
    // dropped, in contrast with cutoff 0. This is the exclusive side of the cutoff-0 boundary.
    [Test]
    public void QualityTrimReads_CutoffOne_EnablesTrimming_AllZeroDropped()
    {
        var reads = new List<(string, string)> { ("ACGTAC", ConstantQuality(6, 0)) };

        IReadOnlyList<string> result = SequenceAssembler.QualityTrimReads(reads, minQuality: 1, minLength: 1);

        result.Should().BeEmpty("cutoff 1 ≥ 1 enables trimming; all-Phred-0 ⇒ empty window ⇒ dropped");
    }

    #endregion

    #region BE — Boundary: empty (empty list / empty read ⇒ guarded, no crash, no /0)

    // No reads → documented empty result; nothing to process, no NullReference (§6.1).
    [Test]
    public void QualityTrimReads_EmptyList_ReturnsEmpty()
    {
        IReadOnlyList<string> result =
            SequenceAssembler.QualityTrimReads(new List<(string, string)>(), minQuality: 20, minLength: 1);

        result.Should().BeEmpty("an empty read list ⇒ empty trimmed result (§6.1)");
    }

    // A single empty read (length-0 sequence + length-0 quality): length 0 < minLength ⇒ dropped,
    // with NO DivideByZero and NO negative-length Substring on the empty scan (§6.1).
    [Test]
    public void QualityTrimReads_SingleEmptyRead_DroppedNoCrash()
    {
        var reads = new List<(string, string)> { ("", "") };

        IReadOnlyList<string> result = SequenceAssembler.QualityTrimReads(reads, minQuality: 20, minLength: 1);

        result.Should().BeEmpty("empty read ⇒ length 0 < minLength ⇒ dropped, no crash (§6.1)");
    }

    // An empty read survives ONLY when minLength is 0 (length 0 ≥ 0): it is emitted as the
    // empty string with NO negative-length Substring. Guards the end−start == 0 path.
    [Test]
    public void QualityTrimReads_EmptyReadWithMinLengthZero_EmittedEmptyNoCrash()
    {
        var reads = new List<(string, string)> { ("", "") };

        Action act = () =>
        {
            IReadOnlyList<string> result = SequenceAssembler.QualityTrimReads(reads, minQuality: 20, minLength: 0);
            result.Should().ContainSingle().Which.Should().BeEmpty(
                "length 0 ≥ minLength 0 ⇒ empty survivor, no negative-length Substring");
        };

        act.Should().NotThrow("the end−start == 0 substring path must be safe");
    }

    // Mix of empty, all-low and all-high reads at the default-ish cutoff: empties/all-low drop,
    // all-high survives, order preserved, no crash on the empty element.
    [Test]
    public void QualityTrimReads_MixedEmptyLowHigh_NoCrash_OrderPreserved()
    {
        var reads = new List<(string, string)>
        {
            ("", ""),                                    // empty   → dropped
            ("AAAAAA", ConstantQuality(6, 2)),           // all-low → dropped
            ("CCCCCC", ConstantQuality(6, 40)),          // all-high → "CCCCCC"
        };

        IReadOnlyList<string> result = SequenceAssembler.QualityTrimReads(reads, minQuality: 20, minLength: 1);

        result.Should().ContainSingle().Which.Should().Be("CCCCCC",
            "only the all-high read survives; empty/all-low dropped, order preserved");
        AssertWellFormed(result, reads, 1);
    }

    #endregion

    #region Validation contract (§3.3)

    // Null reads → documented ArgumentNullException (§3.3, §6.1).
    [Test]
    public void QualityTrimReads_NullReads_Throws()
    {
        Action act = () => SequenceAssembler.QualityTrimReads(null!);

        act.Should().Throw<ArgumentNullException>("null reads is the documented validation contract (§3.3)");
    }

    #endregion

    #region BE — Broad fuzz: random reads/cutoffs match the cutadapt oracle, never throw

    // Sweeps read length, Phred range, cutoff (including 0/negative) and minLength over many
    // random reads and asserts the production trimmer:
    //   • NEVER throws an unexpected runtime exception (no negative-length Substring on full
    //     trim, no IndexOutOfRange at the cutoff boundary);
    //   • EXACTLY matches the independent cutadapt quality_trim_index reference oracle
    //     (same survivor / same drop) — the strongest contract check;
    //   • produces only well-formed survivors (contiguous substrings, length ≤ original,
    //     read/quality in sync) and is DETERMINISTIC.
    [Test]
    [CancelAfter(60_000)]
    public void QualityTrimReads_RandomInputs_MatchOracle_NeverThrows_Deterministic()
    {
        var rng = new Random(148_003);
        for (int trial = 0; trial < 3000; trial++)
        {
            int len = rng.Next(0, 60);
            string seq = RandomDna(rng, len);
            // Wide Phred range so windows cross, isolated good bases appear, etc.
            int lo = rng.Next(0, 45);
            int hi = Math.Min(93, lo + rng.Next(0, 45));
            string qual = RandomQuality(rng, len, lo, hi);
            int cutoff = rng.Next(-3, 45);   // include 0 and negatives (trimming-disabled guard)
            int minLength = rng.Next(0, 10);

            var reads = new List<(string, string)> { (seq, qual) };

            IReadOnlyList<string> result = SequenceAssembler.QualityTrimReads(reads, cutoff, minLength);
            IReadOnlyList<string> again = SequenceAssembler.QualityTrimReads(reads, cutoff, minLength);

            string? expected = ReferenceSurvivor(seq, qual, cutoff, minLength);

            if (expected is null)
            {
                result.Should().BeEmpty(
                    "oracle drops this read (empty window or below minLength); seq='{0}' qual='{1}' C={2} minLen={3}",
                    seq, qual, cutoff, minLength);
            }
            else
            {
                result.Should().ContainSingle().Which.Should().Be(expected,
                    "production must match cutadapt quality_trim_index; seq='{0}' qual='{1}' C={2} minLen={3}",
                    seq, qual, cutoff, minLength);
            }

            result.Should().Equal(again, "trimming is deterministic for identical input");
            AssertWellFormed(result, reads, minLength);
        }
    }

    // Multi-read fuzz: a batch of mixed reads is trimmed read-by-read independently and the
    // surviving sequences appear in input order, exactly matching the per-read oracle.
    [Test]
    [CancelAfter(60_000)]
    public void QualityTrimReads_RandomBatches_PerReadIndependent_OrderPreserved()
    {
        var rng = new Random(148_004);
        for (int trial = 0; trial < 600; trial++)
        {
            int count = rng.Next(0, 8);
            var reads = new List<(string, string)>(count);
            for (int r = 0; r < count; r++)
            {
                int len = rng.Next(0, 40);
                string seq = RandomDna(rng, len);
                int lo = rng.Next(0, 45);
                int hi = Math.Min(93, lo + rng.Next(0, 40));
                reads.Add((seq, RandomQuality(rng, len, lo, hi)));
            }
            int cutoff = rng.Next(-2, 40);
            int minLength = rng.Next(0, 6);

            IReadOnlyList<string> result = SequenceAssembler.QualityTrimReads(reads, cutoff, minLength);

            var expected = reads
                .Select(rd => ReferenceSurvivor(rd.Item1, rd.Item2, cutoff, minLength))
                .Where(s => s is not null)
                .Select(s => s!)
                .ToList();

            result.Should().Equal(expected,
                "each read trimmed independently, survivors in input order, matching the oracle");
            AssertWellFormed(result, reads, minLength);
        }
    }

    #endregion

    #endregion
}
