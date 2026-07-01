// CHROM-CENT-001 — Suprachromosomal-family (SF) assignment (bundled CC0 reference)
// Evidence: docs/Evidence/CHROM-CENT-001-Evidence.md
// TestSpec: tests/TestSpecs/CHROM-CENT-001.md
// Sources:
//   SF taxonomy + A/B-box rule: McNulty SM, Sullivan BA (2018), "Alpha satellite DNA biology"
//     (Chromosome Res; PMC6121732): SF1/SF2 dimeric, SF3 pentameric, SF4 monomeric, SF5 irregular;
//     "A-type monomers include J1, D2, W4, W5, M1, and R2 monomers, while B-type consist of J2, D1,
//     W1-W3, and R1 monomers. B-type monomers contain CENP-B boxes; A-type contain pJalpha boxes."
//   SF classification origin: Shepelev et al. (2009), PLOS Genet 5:e1000641 (CC BY).
//   Bundled reference (CC0): Dfam ALR (DF000000029, A), ALRa (DF000000014, A), ALRb (DF000000015, B);
//     retrieved 2026-06-25 from the Dfam REST API. ALRb carries the 17-bp CENP-B box at consensus
//     position 126; ALR/ALRa do not. (Masumoto et al. 1989 for the CENP-B box.)

using NUnit.Framework;
using Seqeron.Genomics.Chromosome;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests.Unit.Chromosome;

/// <summary>
/// Tests for the opt-in suprachromosomal-family (SF) assignment added to CHROM-CENT-001:
/// <see cref="ChromosomeAnalyzer.AssignSuprachromosomalFamily(string, IReadOnlyList{ChromosomeAnalyzer.AlphaSatelliteReferenceMonomer}?)"/>
/// and <see cref="ChromosomeAnalyzer.LoadBundledAlphaSatelliteReference"/>.
/// <para>
/// Period-dependent arrays are built from the REAL bundled Dfam reference monomers. For the
/// pentameric (period-5) case, five DISTINCT monomers are needed, so mild point-variants of the real
/// reference strings are used (each stays well above the alpha-satellite identity gate and keeps its
/// parent's A/B-box type). All A/B and identity assertions trace to the retrieved Dfam sequences and
/// the sourced A/B rule, not to code echoes.
/// </para>
/// </summary>
[TestFixture]
public class ChromosomeAnalyzer_SuprachromosomalFamily_Tests
{
    private static IReadOnlyList<ChromosomeAnalyzer.AlphaSatelliteReferenceMonomer> Ref()
        => ChromosomeAnalyzer.LoadBundledAlphaSatelliteReference();

    private static string A() => Ref().First(r => r.Name == "ALRa").Sequence;  // A-type, 172 bp
    private static string A171() => Ref().First(r => r.Name == "ALR").Sequence; // A-type, 171 bp
    private static string B() => Ref().First(r => r.Name == "ALRb").Sequence;  // B-type, 169 bp

    private static string Rep(string unit, int n) => string.Concat(Enumerable.Repeat(unit, n));

    // Deterministic mild point mutation: flip a few positions (A<->C, G<->T cycle) to make a DISTINCT
    // but still-alpha-satellite monomer. Used only to synthesise the 5 distinct pentamer monomers.
    private static string Mut(string s, params int[] positions)
    {
        var c = s.ToCharArray();
        foreach (int p in positions)
            c[p] = c[p] switch { 'A' => 'C', 'C' => 'G', 'G' => 'T', _ => 'A' };
        return new string(c);
    }

    #region LoadBundledAlphaSatelliteReference

    // M-SF-1 — bundled CC0 reference loads as the three Dfam consensus monomers with correct A/B type.
    [Test]
    public void LoadBundledAlphaSatelliteReference_ReturnsThreeDfamConsensusMonomers()
    {
        var reference = Ref();

        Assert.That(reference.Count, Is.EqualTo(3),
            "The bundled CC0 reference is exactly the three Dfam alpha-satellite consensus monomers.");

        var alr = reference.First(r => r.Name == "ALR");
        var alra = reference.First(r => r.Name == "ALRa");
        var alrb = reference.First(r => r.Name == "ALRb");

        Assert.Multiple(() =>
        {
            Assert.That(alr.Accession, Is.EqualTo("DF000000029"), "ALR Dfam accession.");
            Assert.That(alr.Sequence.Length, Is.EqualTo(171), "ALR consensus length is 171 bp (Dfam).");
            Assert.That(alr.BoxType, Is.EqualTo(ChromosomeAnalyzer.AlphaSatelliteBoxType.A),
                "ALR has no CENP-B box -> A-type (PMC6121732 A/B rule).");

            Assert.That(alra.Accession, Is.EqualTo("DF000000014"), "ALRa Dfam accession.");
            Assert.That(alra.Sequence.Length, Is.EqualTo(172), "ALRa consensus length is 172 bp (Dfam).");
            Assert.That(alra.BoxType, Is.EqualTo(ChromosomeAnalyzer.AlphaSatelliteBoxType.A),
                "ALRa has no CENP-B box -> A-type.");

            Assert.That(alrb.Accession, Is.EqualTo("DF000000015"), "ALRb Dfam accession.");
            Assert.That(alrb.Sequence.Length, Is.EqualTo(169), "ALRb consensus length is 169 bp (Dfam).");
            Assert.That(alrb.BoxType, Is.EqualTo(ChromosomeAnalyzer.AlphaSatelliteBoxType.B),
                "ALRb carries the CENP-B box -> B-type.");
        });
    }

    // M-SF-2 — the sequence-defined A/B distinction: ALRb carries the 17-bp CENP-B box at consensus
    // position 126; ALR and ALRa do not. This is the basis of the A/B typing.
    [Test]
    public void BundledReference_CenpBBoxPresentOnlyInBTypeAlrb()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ChromosomeAnalyzer.FindCenpBBoxes(B()), Is.EqualTo(new[] { 126 }),
                "ALRb (B-type) carries the 17-bp CENP-B box at consensus position 126 (Masumoto 1989; Dfam).");
            Assert.That(ChromosomeAnalyzer.FindCenpBBoxes(A()), Is.Empty,
                "ALRa (A-type) has no CENP-B box.");
            Assert.That(ChromosomeAnalyzer.FindCenpBBoxes(A171()), Is.Empty,
                "ALR (A-type) has no CENP-B box.");
        });
    }

    #endregion

    #region AssignSuprachromosomalFamily — SF rules

    // M-SF-3 — monomeric A-type array (M1 is A-type, monomeric) -> SF4.
    [Test]
    public void AssignSuprachromosomalFamily_MonomericATypeArray_ReturnsSf4()
    {
        var result = ChromosomeAnalyzer.AssignSuprachromosomalFamily(Rep(A(), 8));

        Assert.Multiple(() =>
        {
            Assert.That(result.IsAlphaSatellite, Is.True, "A real Dfam alpha-satellite monomer matches the reference.");
            Assert.That(result.Family, Is.EqualTo(ChromosomeAnalyzer.SuprachromosomalFamily.Sf4),
                "Monomeric (period 1) A-type only is the SF4 signature (M1; PMC6121732).");
            Assert.That(result.MonomersPerUnit, Is.EqualTo(1), "Homogeneous monomeric array -> period 1.");
            Assert.That(result.BoxTypePattern, Is.EqualTo(new[] { ChromosomeAnalyzer.AlphaSatelliteBoxType.A }),
                "Single A-type monomer in the unit.");
            Assert.That(result.BestReferenceName, Is.EqualTo("ALRa"), "Best match is the source A-type reference.");
        });
    }

    // M-SF-4 — dimeric A.B array (J1.J2 / D1.D2 are A+B per unit) -> {SF1, SF2}.
    [Test]
    public void AssignSuprachromosomalFamily_DimericABArray_ReturnsSf1OrSf2Dimeric()
    {
        var result = ChromosomeAnalyzer.AssignSuprachromosomalFamily(Rep(A() + B(), 6));

        Assert.Multiple(() =>
        {
            Assert.That(result.IsAlphaSatellite, Is.True, "Dimer of real alpha-satellite monomers.");
            Assert.That(result.Family, Is.EqualTo(ChromosomeAnalyzer.SuprachromosomalFamily.Sf1OrSf2Dimeric),
                "Dimeric (period 2) A+B is the SF1/SF2 signature; the CC0 reference cannot separate SF1 from SF2.");
            Assert.That(result.MonomersPerUnit, Is.EqualTo(2), "A-B dimer -> period 2.");
            Assert.That(result.BoxTypePattern, Is.EqualTo(new[]
                {
                    ChromosomeAnalyzer.AlphaSatelliteBoxType.A,
                    ChromosomeAnalyzer.AlphaSatelliteBoxType.B
                }),
                "One A-type then one B-type monomer per HOR unit.");
        });
    }

    // M-SF-5 — pentameric 3 B-type + 2 A-type unit (W1-W3 = B, W4-W5 = A) -> SF3.
    [Test]
    public void AssignSuprachromosomalFamily_PentamericThreeBTwoA_ReturnsSf3()
    {
        // Five DISTINCT monomers so the HOR detector resolves a true period of 5. W1-W3 B-type (from
        // ALRb), W4-W5 A-type (from ALRa); the mutated copies stay > the 60% identity gate.
        string w1 = B();
        string w2 = Mut(B(), 10, 40, 70);
        string w3 = Mut(B(), 15, 45, 75, 100);
        string w4 = A();
        string w5 = Mut(A(), 12, 42, 72, 102);

        var result = ChromosomeAnalyzer.AssignSuprachromosomalFamily(Rep(w1 + w2 + w3 + w4 + w5, 6));

        Assert.Multiple(() =>
        {
            Assert.That(result.IsAlphaSatellite, Is.True, "Pentamer of real-derived alpha-satellite monomers.");
            Assert.That(result.MonomersPerUnit, Is.EqualTo(5), "Five distinct monomers per repeated unit -> period 5.");
            Assert.That(result.Family, Is.EqualTo(ChromosomeAnalyzer.SuprachromosomalFamily.Sf3),
                "Pentameric (period multiple of 5) is the SF3 signature (W1-W5; PMC6121732).");
            Assert.That(result.BoxTypePattern, Is.EqualTo(new[]
                {
                    ChromosomeAnalyzer.AlphaSatelliteBoxType.B,
                    ChromosomeAnalyzer.AlphaSatelliteBoxType.B,
                    ChromosomeAnalyzer.AlphaSatelliteBoxType.B,
                    ChromosomeAnalyzer.AlphaSatelliteBoxType.A,
                    ChromosomeAnalyzer.AlphaSatelliteBoxType.A
                }),
                "W1-W3 are B-type, W4-W5 are A-type (PMC6121732).");
        });
    }

    // M-SF-6 — irregular A/B mix with no regular HOR period (R1.R2-like) -> SF5.
    [Test]
    public void AssignSuprachromosomalFamily_IrregularABMix_ReturnsSf5()
    {
        // Distinct A- and B-type monomers in an irregular order so no regular period is found, but both
        // box types are present -> SF5 (R1 = B-type, R2 = A-type, alternating irregularly).
        string a1 = Mut(A171(), 20, 60, 110);
        string a2 = A();
        string b1 = Mut(B(), 22, 62, 112);
        string b2 = B();

        var result = ChromosomeAnalyzer.AssignSuprachromosomalFamily(a1 + b1 + a2 + b2 + b1 + a1);

        Assert.Multiple(() =>
        {
            Assert.That(result.IsAlphaSatellite, Is.True, "Irregular array of real-derived alpha-satellite monomers.");
            Assert.That(result.Family, Is.EqualTo(ChromosomeAnalyzer.SuprachromosomalFamily.Sf5),
                "Irregular A/B mix with no regular HOR period is the SF5 signature (R1.R2; PMC6121732, Shepelev 2009).");
        });
    }

    // M-SF-7 — random non-alpha-satellite sequence -> not alpha-satellite, Unknown.
    [Test]
    public void AssignSuprachromosomalFamily_RandomSequence_NotAlphaSatelliteUnknown()
    {
        var rng = new Random(1234); // fixed seed -> deterministic
        string random = string.Concat(Enumerable.Range(0, 400).Select(_ => "ACGT"[rng.Next(4)]));

        var result = ChromosomeAnalyzer.AssignSuprachromosomalFamily(random);

        Assert.Multiple(() =>
        {
            Assert.That(result.IsAlphaSatellite, Is.False,
                "Random DNA does not clear the 60% identity gate to the alpha-satellite reference.");
            Assert.That(result.Family, Is.EqualTo(ChromosomeAnalyzer.SuprachromosomalFamily.Unknown),
                "No SF is assigned to a non-alpha-satellite sequence.");
            Assert.That(result.BestReferenceName, Is.Null, "No reference monomer matched.");
        });
    }

    // M-SF-8 — opt-in addition is byte-additive: the existing detectors are unchanged on a fixed array.
    [Test]
    public void AssignSuprachromosomalFamily_DoesNotChangeExistingDetectors()
    {
        string array = Rep(A() + B(), 6);

        var alphaBefore = ChromosomeAnalyzer.DetectAlphaSatellite(array);
        var horBefore = ChromosomeAnalyzer.DetectHigherOrderRepeat(array);

        // Invoke the new opt-in method.
        _ = ChromosomeAnalyzer.AssignSuprachromosomalFamily(array);

        var alphaAfter = ChromosomeAnalyzer.DetectAlphaSatellite(array);
        var horAfter = ChromosomeAnalyzer.DetectHigherOrderRepeat(array);

        Assert.Multiple(() =>
        {
            Assert.That(alphaAfter, Is.EqualTo(alphaBefore), "DetectAlphaSatellite result is unchanged (additive contract).");
            Assert.That(horAfter, Is.EqualTo(horBefore), "DetectHigherOrderRepeat result is unchanged (additive contract).");
        });
    }

    #endregion

    #region Caller-supplied reference / case / sanity (Should / Could)

    // S-SF-1 — a caller-supplied reference overrides the bundled set.
    [Test]
    public void AssignSuprachromosomalFamily_CallerSuppliedReference_IsUsed()
    {
        // A single A-type custom reference: every monomer types A, so a dimer still resolves period 2
        // but with an all-A pattern -> {SF1, SF2} bucket by period (composition aside).
        var custom = new List<ChromosomeAnalyzer.AlphaSatelliteReferenceMonomer>
        {
            new("CUSTOM_A", "NA", A(), ChromosomeAnalyzer.AlphaSatelliteBoxType.A)
        };

        var result = ChromosomeAnalyzer.AssignSuprachromosomalFamily(Rep(A(), 6), custom);

        Assert.Multiple(() =>
        {
            Assert.That(result.IsAlphaSatellite, Is.True, "Monomers match the caller-supplied reference.");
            Assert.That(result.BestReferenceName, Is.EqualTo("CUSTOM_A"),
                "The best match comes from the caller-supplied reference, not the bundled set.");
            Assert.That(result.Family, Is.EqualTo(ChromosomeAnalyzer.SuprachromosomalFamily.Sf4),
                "Monomeric A-type with the custom A reference -> SF4.");
        });
    }

    // S-SF-2 — case insensitivity.
    [Test]
    public void AssignSuprachromosomalFamily_IsCaseInsensitive()
    {
        string array = Rep(A() + B(), 6);

        var upper = ChromosomeAnalyzer.AssignSuprachromosomalFamily(array.ToUpperInvariant());
        var lower = ChromosomeAnalyzer.AssignSuprachromosomalFamily(array.ToLowerInvariant());

        // The record holds an IReadOnlyList (reference-equality in the synthesized record equality),
        // so compare the meaningful fields explicitly.
        Assert.Multiple(() =>
        {
            Assert.That(lower.IsAlphaSatellite, Is.EqualTo(upper.IsAlphaSatellite), "Same alpha-satellite call.");
            Assert.That(lower.Family, Is.EqualTo(upper.Family), "Same family.");
            Assert.That(lower.MonomersPerUnit, Is.EqualTo(upper.MonomersPerUnit), "Same period.");
            Assert.That(lower.BestReferenceName, Is.EqualTo(upper.BestReferenceName), "Same best reference.");
            Assert.That(lower.MeanReferenceIdentity, Is.EqualTo(upper.MeanReferenceIdentity).Within(1e-10), "Same identity.");
            Assert.That(lower.BoxTypePattern, Is.EqualTo(upper.BoxTypePattern), "Same A/B pattern.");
        });
    }

    // C-SF-1 — mean reference identity is a valid percentage and higher for real alpha than random.
    [Test]
    public void AssignSuprachromosomalFamily_MeanReferenceIdentity_HigherForAlphaThanRandom()
    {
        var rng = new Random(99);
        string random = string.Concat(Enumerable.Range(0, 400).Select(_ => "ACGT"[rng.Next(4)]));

        double alphaId = ChromosomeAnalyzer.AssignSuprachromosomalFamily(Rep(A(), 4)).MeanReferenceIdentity;
        double randomId = ChromosomeAnalyzer.AssignSuprachromosomalFamily(random).MeanReferenceIdentity;

        Assert.Multiple(() =>
        {
            Assert.That(alphaId, Is.InRange(0.0, 100.0), "Identity is a percentage in [0,100].");
            Assert.That(randomId, Is.InRange(0.0, 100.0), "Identity is a percentage in [0,100].");
            Assert.That(alphaId, Is.GreaterThan(randomId),
                "A real alpha-satellite monomer is closer to the reference than random DNA.");
        });
    }

    #endregion

    #region Edge cases

    // Empty / null / too-short / empty-reference.
    [Test]
    public void AssignSuprachromosomalFamily_EmptyOrNullOrShort_NotAlphaSatelliteUnknown()
    {
        var empty = ChromosomeAnalyzer.AssignSuprachromosomalFamily("");
        var nul = ChromosomeAnalyzer.AssignSuprachromosomalFamily(null!);
        var shortSeq = ChromosomeAnalyzer.AssignSuprachromosomalFamily(new string('A', 100)); // < 171

        Assert.Multiple(() =>
        {
            Assert.That(empty.IsAlphaSatellite, Is.False, "Empty -> not alpha-satellite.");
            Assert.That(empty.Family, Is.EqualTo(ChromosomeAnalyzer.SuprachromosomalFamily.Unknown), "Empty -> Unknown.");
            Assert.That(empty.MonomersPerUnit, Is.EqualTo(0), "Empty -> no monomers.");
            Assert.That(empty.BoxTypePattern, Is.Empty, "Empty -> empty pattern.");

            Assert.That(nul.IsAlphaSatellite, Is.False, "Null -> not alpha-satellite.");
            Assert.That(nul.Family, Is.EqualTo(ChromosomeAnalyzer.SuprachromosomalFamily.Unknown), "Null -> Unknown.");

            Assert.That(shortSeq.IsAlphaSatellite, Is.False, "Shorter than one monomer -> not alpha-satellite.");
            Assert.That(shortSeq.Family, Is.EqualTo(ChromosomeAnalyzer.SuprachromosomalFamily.Unknown), "Short -> Unknown.");
        });
    }

    [Test]
    public void AssignSuprachromosomalFamily_EmptyReference_Throws()
    {
        var emptyRef = new List<ChromosomeAnalyzer.AlphaSatelliteReferenceMonomer>();

        Assert.Throws<ArgumentException>(
            () => ChromosomeAnalyzer.AssignSuprachromosomalFamily(Rep(A(), 4), emptyRef),
            "An empty reference set is rejected.");
    }

    #endregion
}
