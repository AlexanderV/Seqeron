using System;
using System.Collections.Generic;
using System.Text;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Statistics-area molecular-weight units.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, invalid and boundary inputs to a unit and asserts that
/// the code NEVER fails in an undisciplined way: no hang, no state corruption,
/// and no *unhandled* runtime exception (KeyNotFoundException, OverflowException,
/// NullReferenceException, NaN/Infinity result, …). Every input must result in
/// EITHER a well-defined, theory-correct value, OR a *documented, intentional*
/// validation exception. A raw runtime exception, a NaN, a negative mass, or a
/// hang on garbage input is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SEQ-MW-001 — Molecular Weight Calculation (Statistics)
/// Checklist: docs/checklists/03_FUZZING.md, row 124.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — empty / null (n=0), single monomer (n=1, the
///          critical (n−1)·water = 0 boundary that must NOT go negative), and a
///          very long sequence (linear growth, no overflow, terminates).
///   • MC = Malformed Content — characters absent from the mass table (IUPAC
///          ambiguity codes, gaps, stop, digits, unicode/control junk): skipped,
///          no KeyNotFound, no invented mass.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes BE, MC).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The molecular-weight contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// API entries (SequenceStatistics.cs, src/.../Seqeron.Genomics.Analysis):
///   • SequenceStatistics.CalculateMolecularWeight(string)            — protein Mw (Da)
///   • SequenceStatistics.CalculateNucleotideMolecularWeight(string, bool isDna) — DNA/RNA Mw (Da)
///
/// Documented behaviour (Molecular_Weight_Calculation.md, Test Unit ID SEQ-MW-001):
///   • §2.2 / §4.1: average-isotopic Mw = Σ monomer_mass(recognized_i) − (n − 1)·W,
///     i.e. one water (W = 18.0153 Da) removed per condensation bond; n recognized
///     monomers form (n − 1) bonds.
///   • §3.3 / §5.4: unknown symbols are SKIPPED (no mass, no bond) via TryGetValue —
///     NO KeyNotFound. A sequence with no recognized monomer returns 0.
///   • §2.4 INV-01 / §6.1: MW(null) = MW(empty) = 0 — no exception, NO negative mass
///     from a naive (0 − 1)·W term (guarded: residues == 0 ⇒ return 0).
///   • §2.4 INV-03 / §6.1: a single recognized monomer ⇒ its FREE monomer mass with
///     ZERO bonds ((n−1) = 0), NOT mass − W — the n=1 boundary must not subtract water.
///   • §2.4 INV-02: MW > 0 for any non-empty recognized sequence (each monomer mass > W).
///   • §2.4 INV-04 / §6.1: case-insensitive (input ToUpperInvariant'd).
///
/// Average-mass tables pinned below are independent oracle copies of the documented
/// constants (Molecular_Weight_Calculation.md §2.2, Biopython IUPACData [5]).
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class SequenceMolecularWeightFuzzTests
{
    #region Helpers — oracle mass tables (independent copy of the documented constants)

    private const double Tolerance = 1e-6;

    /// <summary>Average isotopic mass of water (Da) — the per-bond loss term (§2.2).</summary>
    private const double Water = 18.0153;

    /// <summary>The 20 standard one-letter residues — the ONLY characters with a protein mass entry.</summary>
    private const string StandardResidues = "ARNDCEQGHILKMFPSTWYV";

    /// <summary>Average free-amino-acid masses (Da) — oracle copy of §2.2.</summary>
    private static readonly IReadOnlyDictionary<char, double> Aa = new Dictionary<char, double>
    {
        { 'A', 89.0932 },  { 'C', 121.1582 }, { 'D', 133.1027 }, { 'E', 147.1293 },
        { 'F', 165.1891 }, { 'G', 75.0666 },  { 'H', 155.1546 }, { 'I', 131.1729 },
        { 'K', 146.1876 }, { 'L', 131.1729 }, { 'M', 149.2113 }, { 'N', 132.1179 },
        { 'P', 115.1305 }, { 'Q', 146.1445 }, { 'R', 174.201 },  { 'S', 105.0926 },
        { 'T', 119.1192 }, { 'V', 117.1463 }, { 'W', 204.2252 }, { 'Y', 181.1885 }
    };

    /// <summary>Average DNA 5'-monophosphate masses (Da) — oracle copy of §2.2.</summary>
    private static readonly IReadOnlyDictionary<char, double> Dna = new Dictionary<char, double>
    {
        { 'A', 331.2218 }, { 'C', 307.1971 }, { 'G', 347.2212 }, { 'T', 322.2085 }
    };

    /// <summary>Average RNA 5'-monophosphate masses (Da) — oracle copy of §2.2.</summary>
    private static readonly IReadOnlyDictionary<char, double> Rna = new Dictionary<char, double>
    {
        { 'A', 347.2212 }, { 'C', 323.1965 }, { 'G', 363.2206 }, { 'U', 324.1813 }
    };

    /// <summary>Smallest free monomer mass across all tables — the lightest possible single monomer (Glycine,
    /// 75.0666 Da). No recognized non-empty sequence may weigh less than this (INV-02).</summary>
    private const double LightestMonomer = 75.0666;

    /// <summary>Independent Mw oracle that mirrors the documented contract exactly:
    /// Σ recognized monomer masses − (recognizedCount − 1)·W, 0 when none recognized.</summary>
    private static double MwOracle(string seq, IReadOnlyDictionary<char, double> table)
    {
        if (string.IsNullOrEmpty(seq)) return 0.0;
        double sum = 0;
        int n = 0;
        foreach (char c in seq.ToUpperInvariant())
            if (table.TryGetValue(c, out double m)) { sum += m; n++; }
        return n == 0 ? 0.0 : sum - (n - 1) * Water;
    }

    /// <summary>The universal well-formedness contract for ANY input: Mw is a finite number
    /// (never NaN / ±Infinity) and is NON-NEGATIVE — there is no input for which a real
    /// biopolymer mass is negative, and the (n−1)·W water term must never drive the result
    /// below zero (the n=0/n=1 boundary guard). 0 (empty / no-recognized-monomer sentinel)
    /// is admissible.</summary>
    private static void AssertWellFormed(double mw)
    {
        double.IsFinite(mw).Should().BeTrue(
            "Mw must be finite — no NaN, no ±Infinity, no overflow (§8 fuzzing contract)");
        mw.Should().BeGreaterThanOrEqualTo(0.0,
            "Mw is a physical mass and the (n−1)·W water term must never make it negative (INV-01/INV-03)");
    }

    /// <summary>Random string of arbitrary BMP code points (control chars, the null byte,
    /// lone surrogate halves, unicode letters/digits) — fuzz fodder with no mass entry.</summary>
    private static string RandomBmpChars(Random rng, int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = (char)rng.Next(0x0000, 0x10000);
        return new string(chars);
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-MW-001 — molecular weight : fuzz targets (BE, MC)
    // ═══════════════════════════════════════════════════════════════════

    #region Positive sanity — hand-computed exact Mw

    /// <summary>
    /// Positive baseline (not a boundary): the documented worked examples must reproduce
    /// EXACTLY, including the (n−1)·W water-loss term.
    ///   • protein "AGC" = 89.0932 + 75.0666 + 121.1582 − 2·18.0153 = 249.2874 Da
    ///   • DNA "AGC"     = 331.2218 + 347.2212 + 307.1971 − 2·18.0153 = 949.6095 Da
    ///   • RNA "AGC"     = 347.2212 + 363.2206 + 323.1965 − 2·18.0153 = 997.6077 Da
    /// Confirms the suite asserts the real BUSINESS formula (Σ mass − (n−1)·W), not merely
    /// non-throwing. — Molecular_Weight_Calculation.md §7.1.
    /// </summary>
    [Test]
    public void Mw_DocumentedWorkedExamples_MatchHandComputedExactly()
    {
        SequenceStatistics.CalculateMolecularWeight("AGC")
            .Should().BeApproximately(249.2874, Tolerance,
                "89.0932+75.0666+121.1582 − 2·18.0153 = 249.2874 (§7.1)");

        SequenceStatistics.CalculateNucleotideMolecularWeight("AGC", isDna: true)
            .Should().BeApproximately(949.6095, Tolerance,
                "331.2218+347.2212+307.1971 − 2·18.0153 = 949.6095 (§7.1)");

        SequenceStatistics.CalculateNucleotideMolecularWeight("AGC", isDna: false)
            .Should().BeApproximately(997.6077, Tolerance,
                "347.2212+363.2206+323.1965 − 2·18.0153 = 997.6077 (§2.2 table)");
    }

    /// <summary>
    /// Positive baseline: a dipeptide is the INV-03 boundary where exactly ONE water is lost.
    /// "GG" = 2·75.0666 − 1·18.0153 = 132.1179 Da; the difference between MW(GG) and 2·MW(G)
    /// must equal exactly one water. Pins the per-bond water correction directionality.
    /// — Molecular_Weight_Calculation.md §2.4 INV-03.
    /// </summary>
    [Test]
    public void Mw_Dipeptide_RemovesExactlyOneWater()
    {
        double single = SequenceStatistics.CalculateMolecularWeight("G");
        double di = SequenceStatistics.CalculateMolecularWeight("GG");

        di.Should().BeApproximately(132.1179, Tolerance, "2·75.0666 − 18.0153 = 132.1179 (INV-03)");
        (2 * single - di).Should().BeApproximately(Water, Tolerance,
            "a dipeptide is exactly one water lighter than two free residues (INV-03)");
    }

    /// <summary>
    /// Positive baseline: Mw is case-insensitive (INV-04) — the lowercase form must equal
    /// the uppercase Mw exactly. Guards against a missing ToUpperInvariant that would skip
    /// lowercase monomers and corrupt both the mass sum and the bond count.
    /// — Molecular_Weight_Calculation.md §2.4 INV-04.
    /// </summary>
    [Test]
    public void Mw_LowercaseInput_EqualsUppercaseMw()
    {
        SequenceStatistics.CalculateMolecularWeight("acdefghiklmnpqrstvwy")
            .Should().BeApproximately(
                SequenceStatistics.CalculateMolecularWeight("ACDEFGHIKLMNPQRSTVWY"), Tolerance,
                "INV-04: protein Mw is case-insensitive");

        SequenceStatistics.CalculateNucleotideMolecularWeight("acgt", isDna: true)
            .Should().BeApproximately(
                SequenceStatistics.CalculateNucleotideMolecularWeight("ACGT", isDna: true), Tolerance,
                "INV-04: DNA Mw is case-insensitive");

        SequenceStatistics.CalculateNucleotideMolecularWeight("acgu", isDna: false)
            .Should().BeApproximately(
                SequenceStatistics.CalculateNucleotideMolecularWeight("ACGU", isDna: false), Tolerance,
                "INV-04: RNA Mw is case-insensitive");
    }

    #endregion

    #region BE — Boundary: empty / null (Mw 0, no negative mass)

    /// <summary>
    /// BE: the empty string is the n=0 lower size boundary. The documented result is 0,
    /// reached via the `residues == 0 ⇒ return 0` guard — crucially NOT a naive
    /// (0 − 1)·W = −18.0153 Da NEGATIVE mass. Verified for protein, DNA and RNA.
    /// — Molecular_Weight_Calculation.md §6.1 / INV-01.
    /// </summary>
    [Test]
    public void Mw_EmptyString_IsZero_NotNegative()
    {
        SequenceStatistics.CalculateMolecularWeight(string.Empty)
            .Should().Be(0.0, "no residues ⇒ Mw 0, not (0−1)·W (INV-01)");
        SequenceStatistics.CalculateNucleotideMolecularWeight(string.Empty, isDna: true)
            .Should().Be(0.0, "no monomers ⇒ DNA Mw 0 (INV-01)");
        SequenceStatistics.CalculateNucleotideMolecularWeight(string.Empty, isDna: false)
            .Should().Be(0.0, "no monomers ⇒ RNA Mw 0 (INV-01)");

        AssertWellFormed(SequenceStatistics.CalculateMolecularWeight(string.Empty));
    }

    /// <summary>
    /// BE: null is treated identically to empty (IsNullOrEmpty short-circuit) — Mw 0, no
    /// NullReferenceException. — Molecular_Weight_Calculation.md §3.3 (null/empty → 0).
    /// </summary>
    [Test]
    public void Mw_Null_IsZero_NoThrow()
    {
        var protein = () => SequenceStatistics.CalculateMolecularWeight(null!);
        var dna = () => SequenceStatistics.CalculateNucleotideMolecularWeight(null!, isDna: true);

        protein.Should().NotThrow("null is documented as 'no sequence', not an error");
        dna.Should().NotThrow("null is documented as 'no sequence', not an error");
        protein().Should().Be(0.0);
        dna().Should().Be(0.0);
    }

    #endregion

    #region BE — Boundary: single monomer (free mass, NO water subtraction)

    /// <summary>
    /// BE/INV-03: a single recognized monomer is the n=1 boundary where (n−1) = 0, so its
    /// Mw must equal the FREE monomer mass EXACTLY — with NO water subtracted (mass − 0·W).
    /// A negative or water-subtracted single value would prove the (n−1) term mishandled.
    /// Verified for all 20 residues, all 4 DNA bases and all 4 RNA bases against the oracle.
    /// — Molecular_Weight_Calculation.md §6.1 (single monomer ⇒ free mass, zero bonds).
    /// </summary>
    [Test]
    public void Mw_SingleRecognizedMonomer_EqualsFreeMass_NoWaterLoss()
    {
        foreach (var (residue, mass) in Aa)
            SequenceStatistics.CalculateMolecularWeight(residue.ToString())
                .Should().BeApproximately(mass, Tolerance,
                    $"single residue '{residue}' ⇒ free mass {mass}, no peptide bond (INV-03)");

        foreach (var (b, mass) in Dna)
            SequenceStatistics.CalculateNucleotideMolecularWeight(b.ToString(), isDna: true)
                .Should().BeApproximately(mass, Tolerance,
                    $"single DNA base '{b}' ⇒ free monophosphate mass {mass}, no bond");

        foreach (var (b, mass) in Rna)
            SequenceStatistics.CalculateNucleotideMolecularWeight(b.ToString(), isDna: false)
                .Should().BeApproximately(mass, Tolerance,
                    $"single RNA base '{b}' ⇒ free monophosphate mass {mass}, no bond");
    }

    /// <summary>
    /// BE/MC: a single UNRECOGNIZED character is the n=1 boundary where the only monomer is
    /// unknown — recognized count 0, so the documented result is the 0 sentinel, with NO
    /// KeyNotFound and NO negative (0−1)·W mass.
    /// — Molecular_Weight_Calculation.md §5.4 (unknown skipped) + §6.1 (none recognized ⇒ 0).
    /// </summary>
    [TestCase("X")]
    [TestCase("B")]
    [TestCase("Z")]
    [TestCase("-")]
    [TestCase("*")]
    [TestCase("1")]
    [TestCase("?")]
    public void Mw_SingleUnrecognizedChar_IsZero_NoThrow(string seq)
    {
        var protein = () => SequenceStatistics.CalculateMolecularWeight(seq);
        protein.Should().NotThrow($"'{seq}' has no mass entry but must not throw (TryGetValue)");
        protein().Should().Be(0.0, "the only monomer is unrecognized ⇒ count 0 ⇒ Mw 0");

        SequenceStatistics.CalculateNucleotideMolecularWeight(seq, isDna: true)
            .Should().Be(0.0, "unrecognized DNA char ⇒ Mw 0");
    }

    #endregion

    #region MC — Malformed Content: unknown residues skipped (no KeyNotFound)

    /// <summary>
    /// MC: unknown symbols are skipped from BOTH the mass sum and the bond count, so a
    /// recognized sequence PADDED with junk anywhere yields the SAME Mw as the bare sequence.
    /// Guards against a length-based bond count that would over-subtract water for skipped
    /// chars, and confirms no KeyNotFound on ambiguity/gap/stop/digit codes.
    /// — Molecular_Weight_Calculation.md §3.3 / §5.4 (unknown ⇒ no mass, no bond).
    /// </summary>
    [Test]
    public void Mw_JunkPadding_ExcludedFromMassAndBondCount()
    {
        // Protein: "GG" and a junk-padded "G-X*G1" must both reduce to the two G residues.
        SequenceStatistics.CalculateMolecularWeight("G-X*G1")
            .Should().BeApproximately(SequenceStatistics.CalculateMolecularWeight("GG"), Tolerance,
                "junk is skipped — only the two recognized G residues form a single bond");

        // DNA: interior junk must not add mass nor an extra water loss.
        SequenceStatistics.CalculateNucleotideMolecularWeight("A-C*G", isDna: true)
            .Should().BeApproximately(SequenceStatistics.CalculateNucleotideMolecularWeight("ACG", isDna: true), Tolerance,
                "skipped chars form no phosphodiester bond");
    }

    /// <summary>
    /// MC: a string of ALL the common non-standard codes (ambiguity B/Z/X/J, gap, stop,
    /// digits, punctuation) has zero recognized monomers, so Mw is the 0 sentinel with no
    /// KeyNotFound and no negative mass — for protein, DNA and RNA.
    /// — Molecular_Weight_Calculation.md §3.3 (any non-standard char skipped, no exception).
    /// </summary>
    [Test]
    public void Mw_AllNonStandardCodes_IsZero_NoThrow()
    {
        const string junk = "BZXJ-*0123.,;:!@#";

        var protein = () => SequenceStatistics.CalculateMolecularWeight(junk);
        var dna = () => SequenceStatistics.CalculateNucleotideMolecularWeight(junk, isDna: true);
        var rna = () => SequenceStatistics.CalculateNucleotideMolecularWeight(junk, isDna: false);

        protein.Should().NotThrow("no non-standard char has a protein mass entry");
        dna.Should().NotThrow("no non-standard char has a DNA mass entry");
        rna.Should().NotThrow("no non-standard char has an RNA mass entry");

        protein().Should().Be(0.0, "zero recognized residues ⇒ Mw 0");
        dna().Should().Be(0.0);
        rna().Should().Be(0.0);
    }

    /// <summary>
    /// MC: protein-table characters fed to the DNA path (e.g. 'W', 'K' — letters that ARE
    /// amino acids but are NOT DNA bases) must be skipped, not throw. Confirms the per-path
    /// table boundary: only A/C/G/T contribute to DNA Mw.
    /// — Molecular_Weight_Calculation.md §4.2 (DNA/RNA table selected by isDna).
    /// </summary>
    [Test]
    public void Mw_AminoAcidLettersInNucleotidePath_AreSkipped()
    {
        // 'W','K','S','Y' etc. are valid amino acids but not A/C/G/T; only A and G count.
        SequenceStatistics.CalculateNucleotideMolecularWeight("AWKGSY", isDna: true)
            .Should().BeApproximately(SequenceStatistics.CalculateNucleotideMolecularWeight("AG", isDna: true), Tolerance,
                "non-base letters are skipped; only A and G contribute (one bond)");

        // 'U' is RNA-only: in the DNA path it is unrecognized and skipped.
        SequenceStatistics.CalculateNucleotideMolecularWeight("AUG", isDna: true)
            .Should().BeApproximately(SequenceStatistics.CalculateNucleotideMolecularWeight("AG", isDna: true), Tolerance,
                "U has no DNA mass entry ⇒ skipped");
    }

    #endregion

    #region BE — Boundary: very long input (linear growth, no overflow, terminates)

    /// <summary>
    /// BE: a very long homopolymer must remain finite and grow LINEARLY with length, with no
    /// overflow (accumulation in double) and bounded time. For N identical residues of mass m,
    /// Mw = N·m − (N−1)·W; the exact closed form is asserted against the oracle for N up to
    /// 1,000,000. — Molecular_Weight_Calculation.md §4.3 (O(n), double accumulation).
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Mw_VeryLongHomopolymer_GrowsLinearly_NoOverflow()
    {
        foreach (int n in new[] { 1000, 100_000, 1_000_000 })
        {
            string protein = new string('G', n);
            double mw = SequenceStatistics.CalculateMolecularWeight(protein);

            double expected = n * Aa['G'] - (n - 1) * Water;
            mw.Should().BeApproximately(expected, expected * 1e-9,
                $"N·m − (N−1)·W for N={n} all-Glycine (linear, no overflow)");
            AssertWellFormed(mw);

            string dna = new string('A', n);
            double dnaMw = SequenceStatistics.CalculateNucleotideMolecularWeight(dna, isDna: true);
            dnaMw.Should().BeApproximately(n * Dna['A'] - (n - 1) * Water, (n * Dna['A']) * 1e-9,
                $"DNA poly-A Mw grows linearly for N={n}");
            AssertWellFormed(dnaMw);
        }
    }

    /// <summary>
    /// BE/MC: a very long string saturated with unrecognized junk must terminate quickly and
    /// return the 0 sentinel — no KeyNotFound on any char, no negative mass, no hang.
    /// </summary>
    [Test]
    [CancelAfter(15000)]
    public void Mw_VeryLongAllJunk_IsZero_NoCrash()
    {
        string junk = new string('X', 500_000);
        SequenceStatistics.CalculateMolecularWeight(junk)
            .Should().Be(0.0, "no recognized residue at scale ⇒ Mw 0");
        SequenceStatistics.CalculateNucleotideMolecularWeight(junk, isDna: false)
            .Should().Be(0.0, "no recognized RNA monomer at scale ⇒ Mw 0");
    }

    #endregion

    #region MC / BE — Random garbage: never throws, always well-formed, matches oracle

    /// <summary>
    /// MC/BE: a large batch of arbitrary BMP strings (control chars, the null byte, lone
    /// surrogate halves, unicode letters/digits, occasionally seeded with real monomers) must
    /// NEVER throw and ALWAYS yield a well-formed Mw (finite, non-negative) that matches the
    /// independent oracle. Core fuzz guarantee: no KeyNotFound, no NaN, no overflow, no
    /// negative mass on garbage of any shape or length (incl. 0).
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Mw_RandomGarbageStrings_NeverThrow_MatchOracle()
    {
        var rng = new Random(124_001);

        for (int iteration = 0; iteration < 3000; iteration++)
        {
            int len = rng.Next(0, 200);
            string input = RandomBmpChars(rng, len);

            double protein = 0, dna = 0, rnaMw = 0;
            var actP = () => protein = SequenceStatistics.CalculateMolecularWeight(input);
            var actD = () => dna = SequenceStatistics.CalculateNucleotideMolecularWeight(input, isDna: true);
            var actR = () => rnaMw = SequenceStatistics.CalculateNucleotideMolecularWeight(input, isDna: false);

            actP.Should().NotThrow($"garbage (len {len}) must never crash protein Mw");
            actD.Should().NotThrow($"garbage (len {len}) must never crash DNA Mw");
            actR.Should().NotThrow($"garbage (len {len}) must never crash RNA Mw");

            protein.Should().BeApproximately(MwOracle(input, Aa), Tolerance, "protein Mw matches oracle");
            dna.Should().BeApproximately(MwOracle(input, Dna), Tolerance, "DNA Mw matches oracle");
            rnaMw.Should().BeApproximately(MwOracle(input, Rna), Tolerance, "RNA Mw matches oracle");

            AssertWellFormed(protein);
            AssertWellFormed(dna);
            AssertWellFormed(rnaMw);
        }
    }

    /// <summary>
    /// MC/BE: randomly built sequences over the real alphabets (n ≥ 1, mixed with random junk)
    /// must equal the independent oracle (Σ recognized mass − (recognized−1)·W) exactly over
    /// many shapes, and any non-empty recognized sequence must satisfy INV-02 (Mw ≥ the
    /// lightest single monomer mass).
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Mw_RandomSequencesWithJunk_MatchOracle_AndPositive()
    {
        var rng = new Random(124_002);
        const string junk = "XBZJ-*0123 .?\t\n";

        for (int iteration = 0; iteration < 2000; iteration++)
        {
            int monomers = rng.Next(1, 200);
            var sb = new StringBuilder();
            int recognized = 0;
            for (int i = 0; i < monomers; i++)
            {
                sb.Append(StandardResidues[rng.Next(StandardResidues.Length)]);
                recognized++;
                int pad = rng.Next(0, 3);
                for (int p = 0; p < pad; p++) sb.Append(junk[rng.Next(junk.Length)]);
            }
            string seq = sb.ToString();

            double mw = SequenceStatistics.CalculateMolecularWeight(seq);

            mw.Should().BeApproximately(MwOracle(seq, Aa), Tolerance,
                "protein Mw = Σ recognized mass − (recognized−1)·W");
            mw.Should().BeGreaterThanOrEqualTo(LightestMonomer - Tolerance,
                "INV-02: a non-empty recognized sequence weighs at least one (lightest) monomer");
            AssertWellFormed(mw);
        }
    }

    #endregion
}
