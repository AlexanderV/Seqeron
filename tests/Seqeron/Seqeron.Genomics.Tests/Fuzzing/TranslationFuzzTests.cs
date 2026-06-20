using System;
using System.Collections.Generic;
using System.Threading;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Translation area — the codon table / genetic code lookup
/// (TRANS-CODON-001), exposed by <see cref="GeneticCode"/> in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Core/GeneticCode.cs.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, invalid and boundary inputs to a unit and asserts that
/// the code NEVER fails in an undisciplined way: no hang, no state corruption,
/// and no *unhandled* runtime exception (KeyNotFoundException,
/// NullReferenceException, IndexOutOfRangeException, …). Every input must result
/// in EITHER a well-defined, theory-correct value, OR a *documented, intentional*
/// validation exception (ArgumentException / ArgumentNullException). A raw runtime
/// exception or a hang on garbage input is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: TRANS-CODON-001 — codon table / genetic code (Translation)
/// Checklist: docs/checklists/03_FUZZING.md, row 62.
/// Fuzz strategy exercised for THIS unit:
///   • MC = Malformed Content — an invalid / unsupported NCBI table ID, a null
///     codon, and the verification that a non-standard table is actually applied
///     (and is not silently the standard code).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The genetic-code lookup contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// A genetic code maps a nucleotide triplet (codon) to a single-letter amino-acid
/// code; the NCBI numbers the alternative codes as translation tables (1=Standard,
/// 2=Vertebrate Mitochondrial, 3=Yeast Mitochondrial, 11=Bacterial/Plastid).
/// [NCBI, The Genetic Codes]
///   — docs/algorithms/Translation/Codon_Translation.md §1, §2.1.
///
/// API entry points (GeneticCode.cs):
///   • GeneticCode.GetByTableNumber(int tableNumber) — table is identified by its
///     NCBI NUMERIC id. The repository supports ONLY tables 1, 2, 3 and 11. Any
///     other id (0, negative, 9999, etc.) is rejected with an INTENTIONAL
///     ArgumentException — NOT an unhandled KeyNotFoundException
///     (GeneticCode.cs lines 154–161; Codon_Translation.md §2.2 row
///     "[GetByTableNumber] tableNumber", §4.4).
///   • char GeneticCode.Translate(string codon) — throws ArgumentException when the
///     input is null, empty, or not exactly three characters; returns the mapped
///     amino acid; returns 'X' for an ambiguous-but-valid IUPAC codon (e.g. NNN);
///     and throws ArgumentException for a non-IUPAC codon (e.g. "XYZ")
///     (GeneticCode.cs lines 63–79; Codon_Translation.md §3.1, §4.4). Note the
///     contract uses ArgumentException for the null codon (it is the `codon`
///     argument being invalid), NOT a bare NullReferenceException — the key
///     no-crash guarantee for the "null code" fuzz target.
///   • IsStartCodon / IsStopCodon — return FALSE (never throw) on null / empty /
///     wrong-length input (GeneticCode.cs lines 94–113; Codon_Translation.md §4.4).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Proving the non-standard table is actually applied (not ignored)
/// ───────────────────────────────────────────────────────────────────────────
/// Each alternative table is built from the standard base table with table-specific
/// overrides (GeneticCode.cs lines 174–225). To prove the override is real and not
/// silently discarded, this suite pins KNOWN codon reassignments that DIFFER from
/// the standard code (Codon_Translation.md §2.1 "Key AA Differences"):
///   • Table 2 (Vertebrate Mitochondrial) vs Table 1 (Standard):
///       UGA  Stop ('*') → Tryptophan ('W')   (GeneticCode.cs line 181)
///       AGA  Arginine ('R') → Stop ('*')      (GeneticCode.cs line 182)
///       AGG  Arginine ('R') → Stop ('*')      (GeneticCode.cs line 183)
///       AUA  Isoleucine ('I') → Methionine ('M') (GeneticCode.cs line 180)
///   • Table 3 (Yeast Mitochondrial) vs Table 1: CUN Leucine ('L') → Threonine ('T')
///       (GeneticCode.cs lines 199–202).
///   • Table 11 (Bacterial/Plastid) shares the standard CODON table but adds extra
///       START codons (e.g. AUU, AUC, AUA, GUG), so its difference from the standard
///       is in IsStartCodon, not in Translate (GeneticCode.cs lines 214–224).
/// If any of these reassignments did NOT hold, the alternative table would be the
/// standard code in disguise — exactly the "silently used the wrong genetic code"
/// failure this fuzz unit is designed to catch.
///
/// Determinism note: every test uses FIXED, hand-chosen inputs (the genetic-code
/// lookup is a pure deterministic mapping, so no Rng is needed). A local
/// new Random(seed) drives the only randomized scan (the invalid-table-id sweep),
/// with NO shared static Rng.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class TranslationFuzzTests
{
    /// <summary>The four NCBI table numbers this repository supports.</summary>
    private static readonly int[] SupportedTables = { 1, 2, 3, 11 };

    // ═══════════════════════════════════════════════════════════════════
    //  TRANS-CODON-001 — codon table / genetic code : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region TRANS-CODON-001 — codon table / genetic code

    #region Positive sanity — the standard table maps correctly, a non-standard one differs

    /// <summary>
    /// Positive sanity: the standard genetic code (NCBI Table 1) maps representative
    /// codons exactly per theory — AUG/ATG → 'M' (Methionine), UAA/UGA/UAG → '*'
    /// (stop), GGU → 'G' — and a valid-but-ambiguous IUPAC codon (NNN) → 'X'. This
    /// is the working happy path that proves the fuzz targets below are measured
    /// against a functioning lookup, not a uniformly-broken one.
    /// — Codon_Translation.md §6 (worked examples), §3.1.
    /// </summary>
    [Test]
    public void StandardTable_MapsRepresentativeCodons_PerTheory()
    {
        var standard = GeneticCode.GetByTableNumber(1);
        standard.Should().BeSameAs(GeneticCode.Standard);
        standard.TableNumber.Should().Be(1);

        // ATG (DNA) and AUG (RNA) both normalize to Methionine — case-insensitive.
        standard.Translate("ATG").Should().Be('M');
        standard.Translate("AUG").Should().Be('M');
        standard.Translate("aug").Should().Be('M');

        // The three standard stop codons map to the stop sentinel '*'.
        standard.Translate("UAA").Should().Be('*');
        standard.Translate("UAG").Should().Be('*');
        standard.Translate("UGA").Should().Be('*');

        // A plain coding codon.
        standard.Translate("GGU").Should().Be('G');

        // Ambiguous-but-valid IUPAC codon is untranslatable, not invalid → 'X'.
        standard.Translate("NNN").Should().Be('X');

        standard.IsStartCodon("AUG").Should().BeTrue();
        standard.IsStopCodon("UAA").Should().BeTrue();
        standard.IsStopCodon("AUG").Should().BeFalse();
    }

    /// <summary>
    /// Positive sanity (non-standard table is genuinely applied): NCBI Table 2
    /// (Vertebrate Mitochondrial) reassigns documented codons relative to the
    /// standard code. We pin the KNOWN differences so a "silently used the standard
    /// table" regression is caught: UGA *→W, AGA/AGG R→*, AUA I→M. We also assert
    /// these codons DIFFER from the standard table at the same codon, proving the
    /// override (not just an accidental match) is in force.
    /// — Codon_Translation.md §2.1; GeneticCode.cs lines 174–190.
    /// </summary>
    [Test]
    public void NonStandardTable2_AppliesDocumentedReassignments_DifferingFromStandard()
    {
        var standard = GeneticCode.GetByTableNumber(1);
        var mito = GeneticCode.GetByTableNumber(2);
        mito.Should().BeSameAs(GeneticCode.VertebrateMitochondrial);
        mito.TableNumber.Should().Be(2);

        // Reassignments (table 2 value) — these are the load-bearing assertions.
        mito.Translate("UGA").Should().Be('W');   // Stop → Tryptophan
        mito.Translate("AGA").Should().Be('*');   // Arginine → Stop
        mito.Translate("AGG").Should().Be('*');   // Arginine → Stop
        mito.Translate("AUA").Should().Be('M');   // Isoleucine → Methionine

        // ... and each one genuinely DIFFERS from the standard code at that codon.
        standard.Translate("UGA").Should().Be('*');
        standard.Translate("AGA").Should().Be('R');
        standard.Translate("AGG").Should().Be('R');
        standard.Translate("AUA").Should().Be('I');

        mito.Translate("UGA").Should().NotBe(standard.Translate("UGA"));
        mito.Translate("AGA").Should().NotBe(standard.Translate("AGA"));
        mito.Translate("AUA").Should().NotBe(standard.Translate("AUA"));

        // A non-reassigned codon stays identical across both tables (AUG → M).
        mito.Translate("AUG").Should().Be(standard.Translate("AUG"));

        // The reassigned start/stop sets are applied too (UGA is no longer a stop;
        // AGA/AGG become stops; AUA becomes a start).
        mito.IsStopCodon("UGA").Should().BeFalse();
        standard.IsStopCodon("UGA").Should().BeTrue();
        mito.IsStopCodon("AGA").Should().BeTrue();
        mito.IsStartCodon("AUA").Should().BeTrue();
    }

    /// <summary>
    /// Positive sanity (a second non-standard table is applied): NCBI Table 3
    /// (Yeast Mitochondrial) reassigns the four CUN Leucine codons to Threonine,
    /// and NCBI Table 11 (Bacterial/Plastid) shares the standard codon table but
    /// adds extra START codons. These guard against the alternative tables being
    /// the standard code in disguise.
    /// — Codon_Translation.md §2.1; GeneticCode.cs lines 193–224.
    /// </summary>
    [Test]
    public void NonStandardTables3And11_ApplyDocumentedDifferences()
    {
        var standard = GeneticCode.GetByTableNumber(1);
        var yeast = GeneticCode.GetByTableNumber(3);
        var bacterial = GeneticCode.GetByTableNumber(11);

        yeast.TableNumber.Should().Be(3);
        bacterial.TableNumber.Should().Be(11);

        // Table 3: CUN Leucine → Threonine (a real, documented reassignment).
        foreach (var codon in new[] { "CUU", "CUC", "CUA", "CUG" })
        {
            yeast.Translate(codon).Should().Be('T');
            standard.Translate(codon).Should().Be('L');
        }

        // Table 11: codon→AA mapping equals the standard table for every codon...
        foreach (var codon in standard.CodonTable.Keys)
            bacterial.Translate(codon).Should().Be(standard.Translate(codon));

        // ...but it adds START codons absent from the standard set (e.g. AUU, GUG).
        bacterial.IsStartCodon("AUU").Should().BeTrue();
        standard.IsStartCodon("AUU").Should().BeFalse();
        bacterial.IsStartCodon("GUG").Should().BeTrue();
        standard.IsStartCodon("GUG").Should().BeFalse();
    }

    #endregion

    #region Fuzz target: invalid table ID → intentional ArgumentException, never KeyNotFound

    /// <summary>
    /// Fuzz target "invalid table ID" (MC): an unsupported NCBI table number — 0,
    /// a negative id, or a large bogus id like 9999 — must be REJECTED with the
    /// documented, intentional ArgumentException, NOT an unhandled
    /// KeyNotFoundException or any other raw runtime exception. The switch in
    /// GetByTableNumber has an explicit `_ => throw new ArgumentException(...)`
    /// default arm precisely so a lookup miss cannot escape as a crash.
    /// — GeneticCode.cs lines 154–161; Codon_Translation.md §4.4.
    /// </summary>
    [Test]
    public void GetByTableNumber_UnsupportedId_ThrowsArgumentException_NotKeyNotFound()
    {
        int[] invalidIds =
        {
            0, -1, -42, 4, 5, 6, 7, 8, 9, 10, 12, 13, 25, 9999,
            int.MaxValue, int.MinValue,
        };

        foreach (var id in invalidIds)
        {
            Action act = () => GeneticCode.GetByTableNumber(id);
            act.Should().Throw<ArgumentException>(
                    "table id {0} is unsupported and must be rejected intentionally", id)
                .And.ParamName.Should().Be("tableNumber");
        }
    }

    /// <summary>
    /// Fuzz target "invalid table ID" — randomized sweep (MC): a local fixed-seed
    /// Random generates many integers; each one must either be a supported table
    /// (1/2/3/11, returning a non-null GeneticCode with the matching TableNumber) or
    /// throw ArgumentException. No other exception type may ever escape, and no input
    /// may return null — proving the lookup is total over the entire int domain.
    /// </summary>
    [Test]
    public void GetByTableNumber_RandomIds_AlwaysSupportedOrArgumentException()
    {
        var rng = new Random(20260620);

        for (int i = 0; i < 5000; i++)
        {
            int id = rng.Next(int.MinValue, int.MaxValue);

            GeneticCode? result = null;
            ArgumentException? thrown = null;
            try
            {
                result = GeneticCode.GetByTableNumber(id);
            }
            catch (ArgumentException ex)
            {
                thrown = ex;
            }

            if (Array.IndexOf(SupportedTables, id) >= 0)
            {
                thrown.Should().BeNull("table id {0} is supported", id);
                result.Should().NotBeNull();
                result!.TableNumber.Should().Be(id);
            }
            else
            {
                result.Should().BeNull("table id {0} is unsupported", id);
                thrown.Should().NotBeNull(
                    "unsupported table id {0} must throw ArgumentException, not crash", id);
            }
        }
    }

    #endregion

    #region Fuzz target: null / malformed codon → ArgumentException, never NullReference

    /// <summary>
    /// Fuzz target "null code" (MC): Translate(null) must throw the documented
    /// ArgumentException (the `codon` argument is invalid), NEVER a bare
    /// NullReferenceException. Empty and wrong-length codons are also malformed
    /// content and follow the same intentional-rejection contract.
    /// — GeneticCode.cs lines 63–66; Codon_Translation.md §3.1, §4.4.
    /// </summary>
    [Test]
    public void Translate_NullEmptyOrWrongLengthCodon_ThrowsArgumentException_NotNullReference()
    {
        var standard = GeneticCode.GetByTableNumber(1);

        foreach (var bad in new[] { null, "", "A", "AU", "AUGG", "AUGAUG" })
        {
            Action act = () => standard.Translate(bad!);
            act.Should().Throw<ArgumentException>(
                    "malformed codon {0} must be rejected intentionally",
                    bad is null ? "null" : $"\"{bad}\"")
                .And.ParamName.Should().Be("codon");
        }
    }

    /// <summary>
    /// Fuzz target "null code" companion (MC): IsStartCodon / IsStopCodon are the
    /// non-throwing classifiers — on null, empty or wrong-length input they must
    /// return FALSE rather than throw a NullReferenceException. This is the
    /// documented divergence from Translate (which throws on the same input) and
    /// must not crash on malformed content.
    /// — GeneticCode.cs lines 94–113; Codon_Translation.md §4.4.
    /// </summary>
    [Test]
    public void IsStartStopCodon_NullOrMalformedCodon_ReturnsFalse_NeverThrows()
    {
        var standard = GeneticCode.GetByTableNumber(1);

        foreach (var bad in new[] { null, "", "A", "AU", "AUGG", "ZZZ", "12G" })
        {
            standard.IsStartCodon(bad!).Should().BeFalse();
            standard.IsStopCodon(bad!).Should().BeFalse();
        }
    }

    /// <summary>
    /// Fuzz target "malformed codon content" (MC): a 3-character codon whose symbols
    /// are NOT valid IUPAC nucleotides (e.g. "XYZ", "12G", "@#$") is invalid content
    /// and must throw ArgumentException — distinct from a valid-but-ambiguous IUPAC
    /// codon ("NNN", "ANR") which is untranslatable and returns 'X'. Neither path may
    /// surface a KeyNotFoundException from the underlying dictionary lookup.
    /// The scan is bounded by [CancelAfter] as a hang tripwire.
    /// — GeneticCode.cs lines 70–78; Codon_Translation.md §3.1, §6.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void Translate_NonIupacCodon_ThrowsArgumentException_AmbiguousIupacReturnsX(
        CancellationToken token)
    {
        foreach (var table in SupportedTables)
        {
            var code = GeneticCode.GetByTableNumber(table);

            // Non-IUPAC 3-char content → intentional ArgumentException (not KeyNotFound).
            foreach (var bad in new[] { "XYZ", "12G", "@#$", "  G", "...", "JOZ" })
            {
                token.ThrowIfCancellationRequested();
                Action act = () => code.Translate(bad);
                act.Should().Throw<ArgumentException>(
                        "non-IUPAC codon \"{0}\" is invalid content on table {1}", bad, table)
                    .And.ParamName.Should().Be("codon");
            }

            // Valid-but-ambiguous IUPAC codons → untranslatable sentinel 'X', no throw.
            foreach (var ambiguous in new[] { "NNN", "ANR", "YYY", "NNG" })
            {
                token.ThrowIfCancellationRequested();
                code.Translate(ambiguous).Should().Be('X',
                    "ambiguous IUPAC codon \"{0}\" is untranslatable, not invalid", ambiguous);
            }
        }
    }

    #endregion

    #endregion
}
