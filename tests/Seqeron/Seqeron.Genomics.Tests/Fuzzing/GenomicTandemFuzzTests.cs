using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Analysis area — exact tandem-repeat detection (GENOMIC-TANDEM-001), the brute-force
/// detector <see cref="GenomicAnalyzer.FindTandemRepeats(DnaSequence, int, int)"/>: for a DNA sequence it
/// reports every EXACT tandem block — a unit (period) U repeated k ≥ minRepetitions contiguous copies —
/// as a <see cref="TandemRepeat"/> carrying the unit string, the 0-based start <c>Position</c>, the copy
/// count <c>Repetitions</c>, the derived <c>TotalLength = Unit.Length × Repetitions</c> and the rebuilt
/// <c>FullSequence</c>. Per Benson (1999): "a tandem repeat in DNA is two or more contiguous … copies of
/// a pattern of nucleotides", i.e. S[p .. p + k|U|) = U^k with k ≥ 2.
/// — docs/algorithms/Genomic_Analysis/Tandem_Repeat_Detection.md §1, §2.2, §2.4, §3, §4.1, §6.1.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate and boundary sequences to the unit and asserts the code NEVER fails in an
/// undisciplined way: no crash / IndexOutOfRange on the empty sequence, a single base, a sequence with NO
/// tandem (all-distinct / non-periodic content), or a sequence that is ENTIRELY one tandem (homopolymer /
/// "ATATAT…"), and NO hang / quadratic blow-up on dense full-tandem input (guarded with [CancelAfter]).
/// It also asserts the code NEVER produces nonsense: never a copy count below max(2, minRepetitions)
/// (INV-01 — a single would-be unit with k = 1 is NOT a tandem), never a block running out of bounds
/// (INV-03), never TotalLength ≠ Unit.Length × Repetitions (INV-02), never a reported block whose text is
/// not the unit repeated exactly that many times contiguously, never a non-deterministic result. The
/// argument guards (minUnitLength ≥ 1, minRepetitions ≥ 2) are exercised at their boundaries: a 0-length
/// unit never terminates the scan and a k of 0 would divide by zero in the unit-length bound, so both are
/// rejected with <see cref="ArgumentOutOfRangeException"/> rather than hanging / crashing. The subject is
/// a <see cref="DnaSequence"/> (uppercased and validated to {A,C,G,T} at construction), so out-of-domain
/// residues / null bytes / unicode cannot reach the scan. A raw runtime exception, a hang, a k &lt; 2
/// "repeat", an out-of-bounds block, a fabricated/missed tandem, or a non-deterministic output is a bug.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: GENOMIC-TANDEM-001 — exact tandem-repeat detection (FindTandemRepeats), the Analysis-area
/// GenomicAnalyzer tandem detector (rows 175–180). Checklist: docs/checklists/03_FUZZING.md row 180
/// (BE — "no tandem, full tandem, single unit"). Algorithm doc:
/// docs/algorithms/Genomic_Analysis/Tandem_Repeat_Detection.md (Test Unit ID GENOMIC-TANDEM-001).
/// Distinct from row 14 (REP-TANDEM-001): that is the canonical consolidation of the SAME method under the
/// Repeats area; THIS row fuzzes the Analysis-area facade. (The two share the implementation by design.)
///
/// Fuzz strategy exercised for THIS unit (BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt,
/// empty — docs/checklists/03_FUZZING.md §Description), mapped to the row's three targets:
///   • NO TANDEM — a non-periodic subject ("ACGT", "ACGTGCA", random distinct runs): NO tandem reaches the
///     copy threshold, empty enumeration, NO false positive, NO crash (§6.1 "no tandem present → empty").
///   • FULL TANDEM — a subject that is ENTIRELY one repeated unit: a homopolymer ("AAAA…", unit "A") or a
///     period-p tile ("ATATAT…", unit "AT") — the documented maximal block is reported with the correct
///     unit, copy count and span (§6.1 "entire sequence is one tandem → one result spanning the region").
///     No quadratic hang on the dense case (kept modest, guarded by [CancelAfter]; §4.3 O(n²·m)).
///   • SINGLE UNIT — a subject containing exactly ONE copy of a would-be unit (no actual repetition): NOT
///     reported, because a tandem needs k ≥ 2 (INV-01). The k = 1 / k = 2 frontier is the BE boundary; a
///     unit-length-1 (homopolymer) and the argument-guard boundaries (minUnitLength {0,1}, minRepetitions
///     {1,2}) are the minimal / off-by-one cases.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test (Tandem_Repeat_Detection.md §2.4 INV-01..03, §3, §4.1, §6.1)
/// ───────────────────────────────────────────────────────────────────────────
///   INV-01 every reported repeat has Repetitions ≥ max(2, minRepetitions) (k ≥ 2);
///   INV-02 TotalLength == Unit.Length × Repetitions;
///   INV-03 Position + Unit.Length × Repetitions ≤ sequence.Length (block is in-bounds, contiguous);
///   Unit.Length ≥ minUnitLength; the text at Position spells Unit repeated exactly Repetitions times.
///   Empty / non-periodic subject → empty enumeration. minUnitLength &lt; 1 or minRepetitions &lt; 2 throw
///   ArgumentOutOfRangeException (eager, before enumeration). Null sequence → ArgumentNullException.
///   GenomicAnalyzer.FindTandemRepeats(DnaSequence, int minUnitLength, int minRepetitions)
///       → IEnumerable&lt;TandemRepeat&gt;
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class GenomicTandemFuzzTests
{
    private static readonly char[] Alphabet = { 'A', 'C', 'G', 'T' };

    #region Helpers

    /// <summary>A random ACGT string of the given length (length 0 ⇒ empty string).</summary>
    private static string RandomDna(Random rng, int length)
    {
        var sb = new StringBuilder(length);
        for (int i = 0; i < length; i++)
            sb.Append(Alphabet[rng.Next(Alphabet.Length)]);
        return sb.ToString();
    }

    /// <summary>
    /// Independent oracle for FindTandemRepeats, built directly from the documented algorithm (§4.1),
    /// NOT from the unit. Mirrors the greedy scan: for each candidate unit length from minUnitLength up to
    /// ⌊n / minRepetitions⌋, scan start positions left-to-right; at each start extract the candidate unit,
    /// count consecutive exact copies, and if the count reaches minRepetitions yield (unit, start, count)
    /// and advance the start cursor to the end of that block within the current unit-length pass. This is
    /// the exact U^k, k ≥ minRepetitions definition (§2.2) realised with the documented forward-skip.
    /// </summary>
    private static List<(string Unit, int Position, int Repetitions)> Oracle(string seq, int minUnitLength, int minRepetitions)
    {
        var hits = new List<(string, int, int)>();
        for (int unitLen = minUnitLength; unitLen <= seq.Length / minRepetitions; unitLen++)
        {
            for (int start = 0; start <= seq.Length - unitLen * minRepetitions; start++)
            {
                string unit = seq.Substring(start, unitLen);
                int reps = 1;
                int pos = start + unitLen;
                while (pos + unitLen <= seq.Length && seq.Substring(pos, unitLen) == unit)
                {
                    reps++;
                    pos += unitLen;
                }
                if (reps >= minRepetitions)
                {
                    hits.Add((unit, start, reps));
                    start = pos - unitLen; // forward-skip; for() then does start++
                }
            }
        }
        return hits;
    }

    /// <summary>
    /// Asserts a single <see cref="TandemRepeat"/> is WELL-FORMED against the documented invariants and the
    /// text it was found in: unit length ≥ minUnitLength, copy count ≥ max(2, minRepetitions) (INV-01),
    /// TotalLength == Unit.Length × Repetitions (INV-02), the block is in-bounds (INV-03), and the text at
    /// Position spells exactly the unit repeated Repetitions times contiguously (FullSequence agrees).
    /// </summary>
    private static void AssertTandemWellFormed(TandemRepeat r, string seq, int minUnitLength, int minRepetitions)
    {
        r.Unit.Should().NotBeNullOrEmpty("a tandem unit is a non-empty pattern");
        r.Unit.Length.Should().BeGreaterThanOrEqualTo(minUnitLength, "unit length ≥ minUnitLength");

        r.Repetitions.Should().BeGreaterThanOrEqualTo(Math.Max(2, minRepetitions),
            "a tandem has ≥ max(2, minRepetitions) contiguous copies; k ≥ 2 (INV-01)");

        r.TotalLength.Should().Be(r.Unit.Length * r.Repetitions, "TotalLength = period × copy number (INV-02)");

        r.Position.Should().BeGreaterThanOrEqualTo(0, "Position is a 0-based offset");
        (r.Position + r.TotalLength).Should().BeLessThanOrEqualTo(seq.Length,
            "the whole block is in-bounds (no read past the end) (INV-03)");

        // The text at Position is exactly Unit^Repetitions, contiguously (the heart of the definition).
        string spanned = seq.Substring(r.Position, r.TotalLength);
        spanned.Should().Be(r.FullSequence, "the spanned text is the unit repeated Repetitions times (U^k)");
        r.FullSequence.Should().Be(string.Concat(Enumerable.Repeat(r.Unit, r.Repetitions)),
            "FullSequence is the unit concatenated Repetitions times");
        for (int c = 0; c < r.Repetitions; c++)
            seq.Substring(r.Position + c * r.Unit.Length, r.Unit.Length).Should()
                .Be(r.Unit, "every contiguous copy equals the unit (exact tandem)");
    }

    /// <summary>
    /// Runs FindTandemRepeats, asserts every yielded tandem is well-formed, asserts the reported sequence
    /// EXACTLY equals the independent oracle (unit, position, copy count — none missed, none fabricated, no
    /// off-by-one in copy count), and asserts determinism (two enumerations agree). Returns the result.
    /// </summary>
    private static List<TandemRepeat> AssertWellFormed(string text, int minUnitLength, int minRepetitions)
    {
        var dna = new DnaSequence(text);
        string seq = dna.Sequence;

        var result = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength, minRepetitions).ToList();

        foreach (var r in result)
            AssertTandemWellFormed(r, seq, minUnitLength, minRepetitions);

        var actual = result.Select(r => (r.Unit, r.Position, r.Repetitions)).ToList();
        var expected = Oracle(seq, minUnitLength, minRepetitions);
        actual.Should().Equal(expected,
            "the tandem sequence must exactly match the spec oracle (unit, position, copy count, order)");

        var again = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength, minRepetitions)
            .Select(r => (r.Unit, r.Position, r.Repetitions)).ToList();
        again.Should().Equal(actual, "enumeration is deterministic");

        return result;
    }

    #endregion

    #region GENOMIC-TANDEM-001 — tandem repeat detection (BE: no tandem, full tandem, single unit)

    // ─── Positive sanity: a known tandem is found with the documented unit / count / span ───────────

    [Test]
    public void PositiveSanity_DocWorkedExample_FindsUnitWithCorrectCount()
    {
        // Doc §7.1 Wikipedia worked example: "ATTCGATTCGATTCG" = "ATTCG" × 3 (period 5, copy number 3).
        var dna = new DnaSequence("ATTCGATTCGATTCG");
        var hits = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 5, minRepetitions: 2).ToList();

        hits.Should().ContainSingle("exactly one length-5 tandem block (doc §7.1)");
        hits[0].Unit.Should().Be("ATTCG");
        hits[0].Position.Should().Be(0);
        hits[0].Repetitions.Should().Be(3);
        hits[0].TotalLength.Should().Be(15);
        hits[0].FullSequence.Should().Be("ATTCGATTCGATTCG");

        AssertWellFormed("ATTCGATTCGATTCG", minUnitLength: 5, minRepetitions: 2);
    }

    [Test]
    public void PositiveSanity_CagRepeatInFlankingSequence_FoundAtCorrectCoordinates()
    {
        // A "CAG" × 5 disease-style microsatellite embedded in non-repeating flanks: GCT | CAGCAGCAGCAGCAG | TAC.
        // The detector must report unit "CAG", Position 3, Repetitions 5, span 15 — correct copy count + span.
        const string flankLeft = "GCT";   // length 3
        const string flankRight = "TAC";  // length 3
        string text = flankLeft + string.Concat(Enumerable.Repeat("CAG", 5)) + flankRight;

        var dna = new DnaSequence(text);
        var hits = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 3, minRepetitions: 3).ToList();

        var cag = hits.Where(h => h.Unit == "CAG" && h.Repetitions == 5).ToList();
        cag.Should().ContainSingle("the embedded CAG × 5 tandem is reported once at length-3 period");
        cag[0].Position.Should().Be(flankLeft.Length, "block starts immediately after the left flank");
        cag[0].TotalLength.Should().Be(15, "span = 3 × 5 (INV-02)");
        cag[0].FullSequence.Should().Be("CAGCAGCAGCAGCAG");

        AssertWellFormed(text, minUnitLength: 3, minRepetitions: 3);
    }

    // ─── NO TANDEM: non-periodic content ⇒ empty enumeration, no false positive ─────────────────────

    [Test]
    public void NoTandem_NonPeriodicSequence_EmptyNoFalsePositive()
    {
        // "ACGT" / "ACGTGCA": no unit (period ≥ 2) recurs twice contiguously ⇒ no tandem at the defaults.
        AssertWellFormed("ACGT", minUnitLength: 2, minRepetitions: 2)
            .Should().BeEmpty("ACGT is non-periodic; no length-≥2 unit repeats (§6.1)");
        AssertWellFormed("ACGTGCA", minUnitLength: 2, minRepetitions: 2)
            .Should().BeEmpty("ACGTGCA has no contiguous length-≥2 tandem (§6.1)");
        // Even allowing mononucleotide units, distinct adjacent bases give no homopolymer run (and the
        // whole string is not periodic at ANY period ≤ n/2): "ACGTGCAT" has no contiguous repeat.
        AssertWellFormed("ACGTGCAT", minUnitLength: 1, minRepetitions: 2)
            .Should().BeEmpty("ACGTGCAT is aperiodic; no unit (any length) repeats contiguously");
    }

    [Test]
    public void NoTandem_Empty_NoCrash()
    {
        // Empty subject: the scan loops never execute ⇒ empty result, no IndexOutOfRange (§6.1).
        AssertWellFormed(string.Empty, minUnitLength: 1, minRepetitions: 2).Should().BeEmpty("ε has no tandem");
        AssertWellFormed(string.Empty, minUnitLength: 2, minRepetitions: 5).Should().BeEmpty("ε has no tandem");
    }

    [Test]
    public void NoTandem_SingleBase_NoCrash()
    {
        // A single character cannot host two contiguous copies of any unit ⇒ empty, no crash.
        foreach (var s in new[] { "A", "C", "G", "T" })
        {
            AssertWellFormed(s, minUnitLength: 1, minRepetitions: 2).Should()
                .BeEmpty($"'{s}' is one base; no unit recurs (k < 2)");
            AssertWellFormed(s, minUnitLength: 2, minRepetitions: 2).Should()
                .BeEmpty($"'{s}' is shorter than even one length-2 unit");
        }
    }

    [Test]
    public void NoTandem_RandomNonPeriodicWholeStringUnit_NeverTandem()
    {
        // For minUnitLength = n and minRepetitions = 2, no length-n unit can be followed by a second copy
        // (the string is only n long) ⇒ always empty — a clean "no tandem" boundary, no off-by-one.
        var rng = new Random(180_001);
        for (int t = 0; t < 300; t++)
        {
            string s = RandomDna(rng, rng.Next(1, 40));
            AssertWellFormed(s, minUnitLength: s.Length, minRepetitions: 2).Should()
                .BeEmpty("a length-n unit cannot have a second contiguous copy in a length-n string");
        }
    }

    // ─── FULL TANDEM: the whole sequence is one repeated unit (homopolymer / period-p tile) ─────────

    [Test]
    public void FullTandem_Homopolymer_MaximalUnitLengthOneBlock()
    {
        // "AAAAAAAA" (n=8) is entirely the unit "A". At minUnitLength=1, minRepetitions=2 the documented
        // maximal block is "A" × 8 spanning the whole sequence at Position 0 (§6.1 "entire sequence is one
        // tandem → one result"). The greedy forward-skip means it is reported exactly once at unitLen=1.
        var dna = new DnaSequence("AAAAAAAA");
        var hits = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 1, minRepetitions: 2).ToList();

        var len1 = hits.Where(h => h.Unit.Length == 1).ToList();
        len1.Should().ContainSingle("the homopolymer yields one maximal length-1 tandem block");
        len1[0].Unit.Should().Be("A");
        len1[0].Position.Should().Be(0);
        len1[0].Repetitions.Should().Be(8, "all 8 copies form one block (correct copy count, no off-by-one)");
        len1[0].TotalLength.Should().Be(8, "the block spans the whole sequence (INV-02/03)");

        AssertWellFormed("AAAAAAAA", minUnitLength: 1, minRepetitions: 2);
    }

    [Test]
    public void FullTandem_Period2Tile_MaximalUnitReported()
    {
        // "ATATATAT" (n=8) is entirely "AT" × 4. At minUnitLength=2, minRepetitions=2 the maximal length-2
        // block is "AT" × 4 at Position 0 spanning the whole sequence (documented full-tandem case §6.1).
        var dna = new DnaSequence("ATATATAT");
        var hits = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 2, minRepetitions: 2).ToList();

        var at = hits.Where(h => h.Unit == "AT").ToList();
        at.Should().ContainSingle("the period-2 tile yields one maximal 'AT' block");
        at[0].Position.Should().Be(0);
        at[0].Repetitions.Should().Be(4, "AT repeats exactly 4 times (correct copy count)");
        at[0].TotalLength.Should().Be(8, "block spans the whole sequence");

        AssertWellFormed("ATATATAT", minUnitLength: 2, minRepetitions: 2);
    }

    [Test]
    [CancelAfter(60_000)]
    public void FullTandem_DenseHomopolymer_NoQuadraticHang()
    {
        // Dense full-tandem input is the worst case for the O(n²·m) scan (§4.3). Keep modest and guard with
        // [CancelAfter]: a long homopolymer must complete, report the maximal length-1 block, and not hang.
        foreach (int n in new[] { 60, 150, 250 })
        {
            string s = new string('A', n);
            var result = AssertWellFormed(s, minUnitLength: 1, minRepetitions: 2);
            var len1 = result.Where(h => h.Unit.Length == 1).ToList();
            len1.Should().ContainSingle($"length-{n} homopolymer has one maximal length-1 block");
            len1[0].Repetitions.Should().Be(n, "all n bases form one block");
        }
    }

    [Test]
    [CancelAfter(60_000)]
    public void FullTandem_DensePeriod3Tile_NoQuadraticHang()
    {
        // A long "CAG"-tiled sequence (Huntington-style) is dense in tandems; verify it completes, matches
        // the oracle, and reports the maximal "CAG" block with the full copy count — no hang, no off-by-one.
        foreach (int copies in new[] { 30, 80 })
        {
            string s = string.Concat(Enumerable.Repeat("CAG", copies));
            var result = AssertWellFormed(s, minUnitLength: 3, minRepetitions: 2);
            var cag = result.Where(h => h.Unit == "CAG").ToList();
            cag.Should().ContainSingle($"CAG × {copies} yields one maximal length-3 block");
            cag[0].Repetitions.Should().Be(copies, "the full copy count is reported");
            cag[0].Position.Should().Be(0);
        }
    }

    // ─── SINGLE UNIT: exactly one copy (k = 1) ⇒ NOT a tandem (the k = 1 / k = 2 frontier) ─────────

    [Test]
    public void SingleUnit_OneCopyInFlanks_NotReported()
    {
        // "GGG CAG TTT": the candidate unit "CAG" appears exactly ONCE (k = 1). A tandem needs k ≥ 2
        // (INV-01), so it must NOT be reported — the core false-positive guard for the "single unit" target.
        string text = "GGG" + "CAG" + "TTT";
        var result = AssertWellFormed(text, minUnitLength: 3, minRepetitions: 2);
        result.Should().NotContain(h => h.Unit == "CAG",
            "a unit appearing once (k = 1) is not a tandem (INV-01)");
    }

    [Test]
    public void SingleUnit_Versus_TwoCopies_TheBoundary()
    {
        // The k = 1 → k = 2 frontier on the SAME unit: "AT" once (in non-periodic context) is no tandem;
        // "AT" twice ("ATAT") is the minimal tandem. Exercises the off-by-one boundary on the copy count.
        AssertWellFormed("ATGC", minUnitLength: 2, minRepetitions: 2)
            .Should().NotContain(h => h.Unit == "AT", "'AT' occurs once in ATGC (k = 1) ⇒ no tandem");

        var twice = AssertWellFormed("ATAT", minUnitLength: 2, minRepetitions: 2);
        var at = twice.Where(h => h.Unit == "AT").ToList();
        at.Should().ContainSingle("'ATAT' is exactly 'AT' × 2 — the minimal tandem (k = 2)");
        at[0].Repetitions.Should().Be(2);
        at[0].Position.Should().Be(0);
        at[0].TotalLength.Should().Be(4);
    }

    [Test]
    public void SingleUnit_MinRepetitionsThreeRejectsTwoCopies_TheThreshold()
    {
        // With minRepetitions = 3, two copies ("AT" × 2) is BELOW threshold and must NOT be reported, while
        // three copies ("AT" × 3) must be. The k = minReps frontier — never an off-by-one under threshold.
        AssertWellFormed("ATAT", minUnitLength: 2, minRepetitions: 3)
            .Should().BeEmpty("'AT' × 2 is below minRepetitions = 3 (k must reach the threshold)");

        var three = AssertWellFormed("ATATAT", minUnitLength: 2, minRepetitions: 3);
        var at = three.Where(h => h.Unit == "AT").ToList();
        at.Should().ContainSingle("'AT' × 3 meets minRepetitions = 3");
        at[0].Repetitions.Should().Be(3);
    }

    // ─── Argument-guard boundaries: minUnitLength {0} and minRepetitions {1} (BE: 0 / -1) ──────────

    [Test]
    public void Guard_MinUnitLengthBelowOne_Throws_NoHang()
    {
        // minUnitLength < 1 is rejected eagerly (a 0-length unit never advances the scan ⇒ would hang).
        var dna = new DnaSequence("AAAA");
        foreach (int bad in new[] { 0, -1, int.MinValue })
        {
            Action act = () => GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: bad, minRepetitions: 2);
            act.Should().Throw<ArgumentOutOfRangeException>("a 0-length unit never terminates the scan")
                .Which.ParamName.Should().Be("minUnitLength");
        }
    }

    [Test]
    public void Guard_MinRepetitionsBelowTwo_Throws_NoDivideByZero()
    {
        // minRepetitions < 2 is rejected eagerly: k must be ≥ 2 (definition), and k = 0 would divide by
        // zero in the seq.Length / minRepetitions unit-length bound. Boundary values 1, 0, -1.
        var dna = new DnaSequence("AAAA");
        foreach (int bad in new[] { 1, 0, -1, int.MinValue })
        {
            Action act = () => GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 1, minRepetitions: bad);
            act.Should().Throw<ArgumentOutOfRangeException>("k must be ≥ 2; k = 0 would divide by zero")
                .Which.ParamName.Should().Be("minRepetitions");
        }
    }

    [Test]
    public void Guard_ThrowsEagerly_BeforeEnumeration()
    {
        // The guards live in the public method (not the iterator body), so they throw at call time — not
        // deferred to first MoveNext. Constructing the enumerable with bad args must itself throw.
        var dna = new DnaSequence("AAAA");
        Action a1 = () => { _ = GenomicAnalyzer.FindTandemRepeats(dna, 0, 2); };
        Action a2 = () => { _ = GenomicAnalyzer.FindTandemRepeats(dna, 1, 1); };
        a1.Should().Throw<ArgumentOutOfRangeException>("guard is eager, not deferred to MoveNext");
        a2.Should().Throw<ArgumentOutOfRangeException>("guard is eager, not deferred to MoveNext");
    }

    [Test]
    public void Guard_NullSequence_ThrowsArgumentNull()
    {
        Action act = () => GenomicAnalyzer.FindTandemRepeats(null!, 2, 2);
        act.Should().Throw<ArgumentNullException>("a null sequence is rejected up front")
            .Which.ParamName.Should().Be("sequence");
    }

    // ─── Broad random fuzz: well-formedness + oracle equivalence over many shapes ───────────────────

    [Test]
    [CancelAfter(90_000)]
    public void RandomSequences_WellFormed_AndMatchOracle()
    {
        var rng = new Random(180_002);
        for (int t = 0; t < 1500; t++)
        {
            string text = RandomDna(rng, rng.Next(0, 50));
            int minUnit = rng.Next(1, 5);
            int minReps = rng.Next(2, 5);
            AssertWellFormed(text, minUnit, minReps);
        }
    }

    [Test]
    [CancelAfter(90_000)]
    public void RandomSeededTandems_AlwaysFound_AndWellFormed()
    {
        // Embed a guaranteed tandem (unit repeated k ≥ minReps copies) inside random non-tandem flanks,
        // then verify full well-formedness + oracle equivalence and that the seeded unit/count is present.
        var rng = new Random(180_003);
        for (int t = 0; t < 700; t++)
        {
            string unit = RandomDna(rng, rng.Next(1, 5));
            int copies = rng.Next(2, 7);
            string left = RandomDna(rng, rng.Next(0, 6));
            string right = RandomDna(rng, rng.Next(0, 6));
            string text = left + string.Concat(Enumerable.Repeat(unit, copies)) + right;

            var result = AssertWellFormed(text, minUnitLength: unit.Length, minRepetitions: 2);
            // Some unit of length == unit.Length must repeat ≥ copies times somewhere (the embedded block,
            // possibly extended by adjacent equal flank tiles). A non-empty result is guaranteed.
            result.Should().Contain(h => h.Unit.Length == unit.Length && h.Repetitions >= copies,
                "the embedded length-|unit| tandem (≥ copies copies) is reported (INV-01)");
        }
    }

    #endregion
}
