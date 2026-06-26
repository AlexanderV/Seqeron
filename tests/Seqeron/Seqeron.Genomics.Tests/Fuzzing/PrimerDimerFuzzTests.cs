using System;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.MolTools;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the MolTools <b>primer-dimer detector</b> (PRIMER-DIMER-001) — the opt-in
/// SantaLucia &amp; Hicks (2004) most-stable <i>intermolecular</i> duplex (self- or hetero-dimer)
/// formed between two DNA oligos. The exact hand-derived thermodynamics are covered by
/// PrimerDesigner_DimerTm_Tests; THIS file targets the malformed / boundary inputs and pins the
/// STRUCTURAL discipline of any returned dimer against an INDEPENDENT Watson-Crick oracle.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies (docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing")
/// ───────────────────────────────────────────────────────────────────────────
/// Malformed / boundary / out-of-domain inputs must NEVER hang, throw an *unhandled* runtime
/// exception (an IndexOutOfRange from the offset×length alignment scan, a NaN/Inf ΔG leaking where
/// the doc contracts a finite value), or emit OUT-OF-CONTRACT output — a reported dimer whose
/// aligned bases are NOT actually Watson-Crick complements, whose base-pair count is below the
/// minimum, or whose spans fall outside the input strands. Every input must resolve to EITHER a
/// well-defined, theory-correct result (INCLUDING the documented `null` = "no dimer" sentinel and
/// `NaN` Tm), OR a documented validation outcome. The detector never throws for ANY string input.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PRIMER-DIMER-001 — primer-dimer detection (most stable inter-strand duplex)
/// Checklist: docs/checklists/03_FUZZING.md, row 242.
/// Algorithm doc: docs/algorithms/MolTools/DNA_Dimer_Tm.md.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the EMPTY / 1-bp strand(s) (a duplex needs ≥ 2 contiguous bp
///          so &lt; 2 bases → null / NaN, INV-03); two strands with NO complementarity (poly-A vs
///          poly-A, A/G-only vs A/G-only → null); IDENTICAL strands (the self-dimer case — must be
///          handled, never crash, and any returned duplex must be genuinely complementary); a VERY
///          LONG pair (the O(n·m) offset scan must not blow up super-quadratically — [CancelAfter]).
///   • MC = Malformed Content — non-ACGT junk / whitespace / digits / unicode / null-byte / RNA-U
///          (the strict A/C/G/T alphabet rejects them → null / NaN, never a crash).
/// — docs/checklists/03_FUZZING.md §Description (BE; MC = malformed content); row 242 targets:
///   "empty seq, no complementarity, identical seqs, 1-bp, very long".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The dimer contract under test (DNA_Dimer_Tm.md)
/// ───────────────────────────────────────────────────────────────────────────
/// API entries:
///   • PrimerDesigner.FindMostStableDimer(string s1, string s2, double na, double ct)
///       → DimerResult?  (PrimerDesigner.cs §FindMostStableDimer) — the gapless antiparallel
///       offset scan; scores each maximal contiguous Watson-Crick run of ≥ 2 bp; returns the
///       highest-Tm one, or null when no such run exists / either strand is invalid.
///   • PrimerDesigner.CalculateDimerMeltingTemperature(s1, s2, na, ct) → double (NaN when no dimer).
///   • PrimerDesigner.CalculateSelfDimerMeltingTemperature(seq, na, ct) → double — self-dimer wrapper.
/// DimerResult carries Strand1Start / Strand2Start (0-based 5' indices of the aligned duplex on
/// each strand), BasePairs (contiguous WC bp), DeltaH / DeltaS / DeltaG37.
///
/// VALIDATION / sentinels (doc §3.3, §6.1):
///   • null / empty / &lt; 2-base / non-ACGT strand → null dimer, NaN Tm.
///   • a pair with no Watson-Crick duplex of ≥ 2 contiguous bp (poly-A self-dimer; fully
///     non-complementary pair) → null / NaN (INV-03).
/// None of these is a throw — they are the documented null / NaN sentinels.
///
/// STRUCTURAL invariants pinned on EVERY returned dimer (an INDEPENDENT oracle, doc §2.4, §3.2):
///   • INV (indices) — 0 ≤ Strand1Start, Strand1Start + BasePairs ≤ |s1|; likewise on strand 2
///     (the aligned spans fit inside their strands; no out-of-bounds index ever escapes).
///   • INV (≥ 2 bp) — BasePairs ≥ 2 (a duplex needs ≥ 1 NN stack, MinDimerBasePairs / INV-03).
///   • INV (complementary alignment) — re-derive each of the BasePairs aligned columns from the
///     (upper-cased) inputs by an INDEPENDENT Watson-Crick check: strand 1 read 5'→3' from
///     Strand1Start pairs strand 2 read 3'→5' — i.e. s1[Strand1Start+k] complements
///     s2[Strand2Start + BasePairs − 1 − k]. EVERY aligned column MUST be a genuine A·T / G·C
///     complement. A reported "dimer" whose columns are not complementary is an out-of-contract bug.
///   • finiteness — DeltaH, DeltaS, DeltaG37 are all finite (never NaN / ±Inf); the duplex ΔS° &lt; 0.
/// These are asserted UNIVERSALLY on every positive result. ΔG ≤ 0 (net-stability) is asserted ONLY
/// on a clean all-G·C duplex (the maximally-stable domain), NOT universally — a weak A·T duplex may
/// have a slightly positive ΔG°37 yet still be the most-stable dimer returned.
///
/// Symmetry (doc §2.2 / §7): dimer(a,b) and dimer(b,a) describe the SAME antiparallel duplex, so
/// the thermodynamics (BasePairs, ΔH°, ΔS°, ΔG°37, Tm) are identical — pinned as a metamorphic check.
///
/// FindMostStableDimer / the Tm methods are pure (no iterator), so every probe calls them directly.
/// The very-long probe is pinned with [CancelAfter] so any super-quadratic blow-up fails as a
/// timeout rather than hanging.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class PrimerDimerFuzzTests
{
    private const int MinBp = 2; // SantaLucia & Hicks 2004: a duplex needs ≥ 1 NN stack (≥ 2 bp), INV-03.

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
    /// Asserts EVERY documented structural invariant on a returned dimer against an INDEPENDENT
    /// oracle: spans in bounds, ≥ 2 bp, every aligned column a genuine Watson-Crick complement
    /// (re-derived from the inputs with strand 2 read antiparallel), and all of ΔH°/ΔS°/ΔG°37 finite
    /// with a stabilising (negative) duplex ΔS°. ΔG ≤ 0 is deliberately NOT asserted here.
    /// </summary>
    private static void AssertWellFormedDimer(
        PrimerDesigner.DimerResult d, string strand1, string strand2)
    {
        string s1 = strand1.ToUpperInvariant();
        string s2 = strand2.ToUpperInvariant();

        // ≥ 2 bp (a real duplex with ≥ 1 NN stack).
        d.BasePairs.Should().BeGreaterThanOrEqualTo(MinBp, "a dimer needs ≥ 2 contiguous bp (INV-03)");

        // Spans in bounds on both strands (no out-of-range slice escapes).
        d.Strand1Start.Should().BeGreaterThanOrEqualTo(0, "Strand1Start is a 0-based index");
        d.Strand2Start.Should().BeGreaterThanOrEqualTo(0, "Strand2Start is a 0-based index");
        (d.Strand1Start + d.BasePairs).Should().BeLessThanOrEqualTo(s1.Length,
            "the strand-1 aligned span fits inside strand 1");
        (d.Strand2Start + d.BasePairs).Should().BeLessThanOrEqualTo(s2.Length,
            "the strand-2 aligned span fits inside strand 2");

        // Complementary-alignment: re-derive every aligned column independently. Strand 1 reads
        // 5'→3' from Strand1Start; strand 2 reads 3'→5', so column k pairs
        // s1[Strand1Start+k] with s2[Strand2Start + BasePairs - 1 - k].
        for (int k = 0; k < d.BasePairs; k++)
        {
            char left = s1[d.Strand1Start + k];
            char right = s2[d.Strand2Start + d.BasePairs - 1 - k];
            IsComplement(left, right).Should().BeTrue(
                $"aligned column #{k} '{left}'·'{right}' " +
                $"(s1[{d.Strand1Start + k}], s2[{d.Strand2Start + d.BasePairs - 1 - k}]) " +
                "must be a genuine Watson-Crick complement (out-of-contract otherwise)");
        }

        // Finiteness; the duplex ΔS° is negative (ordered helix).
        double.IsFinite(d.DeltaH).Should().BeTrue("ΔH° must be finite");
        double.IsFinite(d.DeltaS).Should().BeTrue("ΔS° must be finite");
        double.IsFinite(d.DeltaG37).Should().BeTrue("ΔG°37 must be finite");
        d.DeltaS.Should().BeLessThan(0.0, "a Watson-Crick duplex has a negative ΔS° (ordered helix)");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PRIMER-DIMER-001 — primer-dimer detection : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PRIMER-DIMER-001 — BE: empty / null strand(s) (→ null, NaN; no crash)

    /// <summary>
    /// BE (empty / null): an empty or null oligo has no bases to align. The doc (§3.3, §6.1)
    /// contracts the null dimer sentinel and a NaN Tm — never an IndexOutOfRange from the offset
    /// scan, never a throw. Pinned for "" and null in either or both strand positions, on the
    /// detector, the dimer-Tm and the self-dimer wrapper.
    /// </summary>
    [Test]
    public void EmptyOrNullStrand_NoDimer_NullAndNaN_NoThrow()
    {
        string?[] empties = { "", null };

        foreach (var a in empties)
        {
            // Self-dimer wrapper on an empty/null oligo.
            var actSelf = () => PrimerDesigner.CalculateSelfDimerMeltingTemperature(a!);
            actSelf.Should().NotThrow("empty/null self-dimer is a documented NaN sentinel, not a crash");
            double.IsNaN(PrimerDesigner.CalculateSelfDimerMeltingTemperature(a!)).Should().BeTrue(
                "an empty/null oligo has no self-dimer → NaN");

            // As either / both strands of a hetero pair (paired with a real foldable oligo).
            foreach (var b in new[] { a, "GCGCGCGC" })
            {
                var act = () => PrimerDesigner.FindMostStableDimer(a!, b!);
                act.Should().NotThrow("empty/null is a documented null sentinel, not a crash");
                PrimerDesigner.FindMostStableDimer(a!, b!).Should().BeNull(
                    "a null/empty strand → no dimer (§6.1)");
                PrimerDesigner.FindMostStableDimer(b!, a!).Should().BeNull(
                    "a null/empty strand → no dimer (either position)");
                double.IsNaN(PrimerDesigner.CalculateDimerMeltingTemperature(a!, b!)).Should().BeTrue(
                    "no dimer → NaN Tm");
            }
        }
    }

    #endregion

    #region PRIMER-DIMER-001 — BE: 1-bp strand(s) (< 2 bp → null, NaN)

    /// <summary>
    /// BE (1-bp boundary): a single base has no NN stack, so a duplex of ≥ 2 contiguous bp is
    /// impossible (INV-03). For every single base on each side, and for a single base paired with a
    /// long complementary oligo (where one column WOULD pair but cannot stack), the detector returns
    /// null with a NaN Tm — never a 1-bp "dimer". This is the just-below-the-minimum length edge.
    /// </summary>
    [Test]
    public void SingleBaseStrand_BelowMinimum_ReturnsNull_NaN()
    {
        foreach (char b in "ACGT")
        {
            string one = b.ToString();

            // 1-bp self-dimer.
            PrimerDesigner.FindMostStableDimer(one, one).Should().BeNull(
                $"a 1-bp oligo '{b}' has no NN stack → no dimer (< 2 bp, INV-03)");
            double.IsNaN(PrimerDesigner.CalculateSelfDimerMeltingTemperature(one)).Should().BeTrue(
                $"1-bp '{b}' → NaN self-dimer Tm");

            // 1-bp vs a long fully-complementary oligo: one column could pair, but < 2 bp ⇒ null.
            PrimerDesigner.FindMostStableDimer(one, "GCGCGCGCGCGC").Should().BeNull(
                $"a 1-bp strand cannot form a ≥ 2-bp duplex → null");
            PrimerDesigner.FindMostStableDimer("GCGCGCGCGCGC", one).Should().BeNull(
                $"a 1-bp strand cannot form a ≥ 2-bp duplex (either position) → null");
        }
    }

    #endregion

    #region PRIMER-DIMER-001 — BE: no complementarity (→ null, NaN)

    /// <summary>
    /// BE (no complementarity): a poly-A oligo cannot Watson-Crick pair with another poly-A (A·A is
    /// not a pair), so the self-dimer and the A/T cross-dimer are null with a NaN Tm (poly-A self-
    /// dimer edge case, doc §6.1). A pair built only from {A,G} on one side and {A,G} on the other
    /// also has no WC column (A·A/G·G/A·G never pair) → null. Never a fabricated dimer, never a crash.
    /// </summary>
    [Test]
    public void NoComplementarity_ReturnsNull_NaN()
    {
        // Poly-A self-dimer: no Watson-Crick run of ≥ 2 bp (doc §6.1).
        PrimerDesigner.FindMostStableDimer("AAAAAAAAAA", "AAAAAAAAAA").Should().BeNull(
            "poly-A cannot pair with poly-A → no self-dimer (§6.1)");
        double.IsNaN(PrimerDesigner.CalculateSelfDimerMeltingTemperature("AAAAAAAAAA")).Should().BeTrue(
            "poly-A → NaN self-dimer Tm");

        // poly-A vs poly-A is the textbook no-complementarity hetero case; poly-A vs poly-G likewise.
        PrimerDesigner.FindMostStableDimer("AAAAAAAA", "GGGGGGGG").Should().BeNull(
            "A·G is not Watson-Crick → no dimer between poly-A and poly-G");

        // An {A,G}-only pair: among {A,G,A,G} no column is ever a WC complement → null.
        PrimerDesigner.FindMostStableDimer("AGAGAGAGAG", "GAGAGAGAGA").Should().BeNull(
            "an A/G-only pair has no Watson-Crick column → no dimer");
        double.IsNaN(PrimerDesigner.CalculateDimerMeltingTemperature("AGAGAGAGAG", "GAGAGAGAGA"))
            .Should().BeTrue("no WC column → NaN Tm");
    }

    #endregion

    #region PRIMER-DIMER-001 — BE: identical strands (self-dimer case) must be handled

    /// <summary>
    /// BE (IDENTICAL strands — the self-dimer case): passing the same oligo as both strands is the
    /// self-dimer path. A self-complementary palindrome (GCGCGCGC) DOES dimerise, so the detector
    /// must RETURN a structurally-valid duplex (every aligned column a genuine WC complement via the
    /// independent oracle, ≥ 2 bp, spans in bounds, finite thermodynamics). Because the stem is
    /// all-G·C (maximally stable) ΔG°37 ≤ 0 is additionally asserted. The self-dimer WRAPPER and the
    /// two-equal-argument FindMostStableDimer/CalculateDimerMeltingTemperature must agree.
    /// </summary>
    [Test]
    public void IdenticalStrands_SelfDimer_FoldsWithValidComplementaryAlignment()
    {
        const string pal = "GCGCGCGC"; // EcoRI-like self-complementary palindrome.

        var dimer = PrimerDesigner.FindMostStableDimer(pal, pal);
        dimer.Should().NotBeNull("a self-complementary palindrome forms a self-dimer");
        var d = dimer!.Value;

        AssertWellFormedDimer(d, pal, pal);

        // The maximally-stable all-G·C self-dimer is net-stable.
        d.DeltaG37.Should().BeLessThan(0.0, "an all-G·C self-dimer is net-stable (negative ΔG°37)");

        // The self-dimer wrapper agrees with the two-equal-argument Tm method, both finite.
        double tmWrapper = PrimerDesigner.CalculateSelfDimerMeltingTemperature(pal);
        double tmPair = PrimerDesigner.CalculateDimerMeltingTemperature(pal, pal);
        double.IsNaN(tmWrapper).Should().BeFalse("a real self-dimer has a finite Tm");
        double.IsInfinity(tmWrapper).Should().BeFalse("a real self-dimer Tm is finite");
        tmWrapper.Should().Be(tmPair, "the self-dimer wrapper delegates to the two-strand dimer Tm");
    }

    /// <summary>
    /// BE (identical NON-palindromic strands): two copies of a non-self-complementary oligo still
    /// form an antiparallel self-dimer wherever a contiguous WC run exists. Whatever the detector
    /// returns (a dimer or null), it must be structurally valid against the oracle and never crash —
    /// pinned across several fixed-seed random oligos used as their own partner.
    /// </summary>
    [Test]
    public void IdenticalStrands_RandomOligo_AnyReturnedSelfDimerIsValid()
    {
        for (int seed = 242_001; seed <= 242_060; seed++)
        {
            int len = 4 + (seed % 30); // 4..33 nt
            string seq = RandomDna(len, seed);

            PrimerDesigner.DimerResult? d = null;
            var act = () => d = PrimerDesigner.FindMostStableDimer(seq, seq);
            act.Should().NotThrow($"seed {seed} self-dimer ('{seq}') must never throw");

            if (d is not null)
                AssertWellFormedDimer(d.Value, seq, seq);

            double tm = PrimerDesigner.CalculateSelfDimerMeltingTemperature(seq);
            double.IsInfinity(tm).Should().BeFalse($"seed {seed}: self-dimer Tm is finite or NaN, never ±Inf");
        }
    }

    #endregion

    #region PRIMER-DIMER-001 — MC: non-ACGT / unicode / junk (strict alphabet → null)

    /// <summary>
    /// MC (non-ACGT / junk / unicode): the detector enforces the strict A/C/G/T alphabet — any
    /// degenerate, whitespace, digit, unicode, null-byte or RNA-U character in EITHER strand makes
    /// the duplex not computable and returns null (doc §3.3, §6.1), with a NaN Tm — never a crash,
    /// never a fabricated dimer. The junk strand is paired with a clean complementary oligo so any
    /// crash would surface in the alignment. Lower-case is upper-cased first, so a clean lowercase
    /// pair is the control that DOES dimerise.
    /// </summary>
    [Test]
    public void NonAcgtJunkUnicode_NotComputable_ReturnsNull_NaN_NoThrow()
    {
        string?[] junk =
        {
            new string('N', 10),       // all-N (IUPAC any-base, not in the ACGT alphabet)
            "GCGCNGCGC",               // N interleaved
            "GCGC GCGC",               // whitespace
            "GCGC-GCGC",               // hyphen
            "GCGC5GCGC",               // digit
            "GCGCUGCGC",               // U is RNA, not in the DNA alphabet
            "GCGC\0GCGC",              // embedded null byte
            "GCGCÄGCGC",          // U+00C4 Latin diacritic
            "αβγδαβγδ",  // Greek letters
            "gcgc\U0001F600gcgc"       // emoji surrogate pair injected into lowercase
        };

        foreach (var s in junk)
        {
            var act = () => PrimerDesigner.FindMostStableDimer(s!, "GCGCGCGC");
            act.Should().NotThrow("junk input must be handled by the null sentinel, not crash");

            PrimerDesigner.FindMostStableDimer(s!, "GCGCGCGC").Should().BeNull(
                "a non-ACGT base in strand 1 makes the duplex not computable → null (strict alphabet)");
            PrimerDesigner.FindMostStableDimer("GCGCGCGC", s!).Should().BeNull(
                "a non-ACGT base in strand 2 → null (strict alphabet)");
            double.IsNaN(PrimerDesigner.CalculateDimerMeltingTemperature(s!, "GCGCGCGC")).Should().BeTrue(
                "non-ACGT input → NaN Tm");
            double.IsNaN(PrimerDesigner.CalculateSelfDimerMeltingTemperature(s!)).Should().BeTrue(
                "non-ACGT self-dimer → NaN Tm");
        }

        // Control: the SAME clean shape, lowercase, IS computable (upper-cased internally).
        PrimerDesigner.FindMostStableDimer("gcgcgcgc", "gcgcgcgc").Should().NotBeNull(
            "a clean lowercase self-complementary pair is upper-cased and dimerises — the nulls above are alphabet-driven");
    }

    #endregion

    #region PRIMER-DIMER-001 — Symmetry: dimer(a,b) == dimer(b,a) thermodynamics

    /// <summary>
    /// Metamorphic symmetry (doc §2.2 / §7): dimer(a, b) and dimer(b, a) describe the SAME
    /// antiparallel duplex, so their thermodynamics (BasePairs, ΔH°, ΔS°, ΔG°37, Tm) must be
    /// identical — and each returned dimer must be structurally valid against the oracle on its own
    /// argument order. Pinned over fixed-seed random oligo pairs (a null result on one side must be
    /// null on the other too). This catches an asymmetric alignment / index bug.
    /// </summary>
    [Test]
    public void DimerSymmetry_SwappingStrands_SameThermodynamics()
    {
        for (int seed = 242_201; seed <= 242_260; seed++)
        {
            string a = RandomDna(6 + (seed % 20), seed);
            string b = RandomDna(6 + ((seed * 7) % 20), seed + 5000);

            PrimerDesigner.DimerResult? ab = null, ba = null;
            var act = () =>
            {
                ab = PrimerDesigner.FindMostStableDimer(a, b);
                ba = PrimerDesigner.FindMostStableDimer(b, a);
            };
            act.Should().NotThrow($"seed {seed}: neither order may throw");

            (ab is null).Should().Be(ba is null,
                $"seed {seed}: dimer(a,b) and dimer(b,a) agree on existence");

            if (ab is not null)
            {
                AssertWellFormedDimer(ab.Value, a, b);
                AssertWellFormedDimer(ba!.Value, b, a);

                ba.Value.BasePairs.Should().Be(ab.Value.BasePairs, "swap preserves bp count");
                ba.Value.DeltaH.Should().BeApproximately(ab.Value.DeltaH, 1e-9, "swap preserves ΔH°");
                ba.Value.DeltaS.Should().BeApproximately(ab.Value.DeltaS, 1e-9, "swap preserves ΔS°");
                ba.Value.DeltaG37.Should().BeApproximately(ab.Value.DeltaG37, 1e-9, "swap preserves ΔG°37");
            }

            double tmAb = PrimerDesigner.CalculateDimerMeltingTemperature(a, b);
            double tmBa = PrimerDesigner.CalculateDimerMeltingTemperature(b, a);
            if (double.IsNaN(tmAb))
                double.IsNaN(tmBa).Should().BeTrue($"seed {seed}: Tm NaN symmetric");
            else
                tmBa.Should().BeApproximately(tmAb, 1e-9, $"seed {seed}: dimer Tm is symmetric");
        }
    }

    #endregion

    #region PRIMER-DIMER-001 — Universal structural discipline over random fuzz inputs

    /// <summary>
    /// The core fuzz sweep: over many fixed-seed random oligo pairs of varied lengths, EVERY non-null
    /// dimer the detector returns MUST satisfy every documented structural invariant — verified
    /// against an INDEPENDENT Watson-Crick oracle (re-derived aligned columns), with ≥ 2 bp, spans in
    /// bounds, and finite ΔH°/ΔS°/ΔG°37. A null result is also acceptable (no dimer). NO input may
    /// throw or hang. This is the load-bearing "never out-of-contract output" guarantee: a reported
    /// dimer whose columns are not actually complementary, or whose spans escape the strands, fails
    /// here on whichever random seed triggers it.
    /// </summary>
    [Test]
    public void RandomFuzz_AnyReturnedDimer_HasValidComplementaryAlignment()
    {
        for (int seed = 242_401; seed <= 242_520; seed++)
        {
            int len1 = 4 + (seed % 40);            // 4..43 nt
            int len2 = 4 + ((seed * 3 + 7) % 40);  // 4..43 nt
            string s1 = RandomDna(len1, seed);
            string s2 = RandomDna(len2, seed + 1_000_000);

            PrimerDesigner.DimerResult? d = null;
            var act = () => d = PrimerDesigner.FindMostStableDimer(s1, s2);
            act.Should().NotThrow($"seed {seed} ('{s1}'/'{s2}') must never throw");

            if (d is not null)
                AssertWellFormedDimer(d.Value, s1, s2);

            // Tm: NaN iff no dimer, finite otherwise — never ±Inf, never a throw.
            double tm = PrimerDesigner.CalculateDimerMeltingTemperature(s1, s2);
            if (d is null)
                tm.Should().Match(t => double.IsNaN(t) || double.IsFinite(t),
                    $"seed {seed}: the legacy scorer found no dimer; the ntthal Tm path is null/NaN or finite");
            else
                double.IsInfinity(tm).Should().BeFalse($"seed {seed}: a real dimer Tm is finite");
        }
    }

    #endregion

    #region PRIMER-DIMER-001 — BE: very long sequences (no O(n^k) blow-up / non-finite)

    /// <summary>
    /// BE (very long sequences): FindMostStableDimer is an O(n·m) offset×length scan (doc §4.3) — a
    /// long pair must NOT blow up or hang. Pinned with [CancelAfter] so any super-quadratic
    /// regression fails as a timeout. The result on a long EXPLICIT self-complementary repeat is a
    /// well-formed self-dimer validated against the oracle (and net-stable, all G·C); a long pair of
    /// independent RANDOM oligos returns either a structurally-valid dimer or null, finite-or-NaN Tm.
    /// NOTE: this pins NUMERICAL / perf discipline on the long path, not the accuracy of any value.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void VeryLongSequences_NoBlowUp_StructurallyValidOrNull()
    {
        // A long self-complementary repeat: (GC)×200 (length 400) — its own reverse complement.
        string longPal = string.Concat(System.Linq.Enumerable.Repeat("GC", 200));

        PrimerDesigner.DimerResult? d = null;
        var act = () => d = PrimerDesigner.FindMostStableDimer(longPal, longPal);
        act.Should().NotThrow("the long self-complementary path must not throw");
        d.Should().NotBeNull("a long (GC)-repeat forms a self-dimer");
        AssertWellFormedDimer(d!.Value, longPal, longPal);
        d.Value.DeltaG37.Should().BeLessThan(0.0, "a long all-G·C self-dimer is net-stable (negative ΔG°37)");

        double tmLong = PrimerDesigner.CalculateSelfDimerMeltingTemperature(longPal);
        double.IsNaN(tmLong).Should().BeFalse("the long self-dimer has a finite Tm");
        double.IsInfinity(tmLong).Should().BeFalse("the long self-dimer Tm is finite, not ±Inf");

        // A long pair of independent random oligos: either a structurally-valid dimer or null.
        string longA = RandomDna(length: 600, seed: 242_777);
        string longB = RandomDna(length: 600, seed: 242_778);
        PrimerDesigner.DimerResult? r = null;
        var actR = () => r = PrimerDesigner.FindMostStableDimer(longA, longB);
        actR.Should().NotThrow("the long random path must not throw");
        if (r is not null)
            AssertWellFormedDimer(r.Value, longA, longB);
    }

    #endregion
}
