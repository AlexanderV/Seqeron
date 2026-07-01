// SEQ-RNACOMP-001 — RNA-specific Complement
// Evidence: docs/Evidence/SEQ-RNACOMP-001-Evidence.md
// TestSpec: tests/TestSpecs/SEQ-RNACOMP-001.md
// Source: Biopython Bio/Data/IUPACData.py (ambiguous_rna_complement); Bio/Seq.py (complement_rna, "T treated as U");
//         Bio.Seq 1.79 docs; bioinformatics.org IUPAC table; NC-IUB 1984 (Cornish-Bowden, NAR 13(9):3021).

using System;
using NUnit.Framework;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests.Unit.Core;

/// <summary>
/// Tests for <see cref="SequenceExtensions.GetRnaComplementBase(char)"/> (SEQ-RNACOMP-001).
/// Expected values are taken from the IUPAC RNA complement table per Biopython
/// `ambiguous_rna_complement` and the NC-IUB 1984 nomenclature, with the repository
/// uppercase-normalization convention (see Evidence Assumption #1).
/// </summary>
[TestFixture]
public class SequenceExtensions_GetRnaComplementBase_Tests
{
    #region GetRnaComplementBase

    // M1 — Standard RNA pairing A↔U, C↔G (Biopython ambiguous_rna_complement; bioinformatics.org IUPAC table).
    [Test]
    public void GetRnaComplementBase_StandardBasesUppercase_ReturnsRnaComplement()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceExtensions.GetRnaComplementBase('A'), Is.EqualTo('U'), "A pairs with U in RNA (not T)");
            Assert.That(SequenceExtensions.GetRnaComplementBase('U'), Is.EqualTo('A'), "U pairs with A");
            Assert.That(SequenceExtensions.GetRnaComplementBase('C'), Is.EqualTo('G'), "C pairs with G");
            Assert.That(SequenceExtensions.GetRnaComplementBase('G'), Is.EqualTo('C'), "G pairs with C");
        });
    }

    // M2 — Thymine treated as uracil → T maps to A (Biopython Seq.py: ambiguous_rna_complement["T"]=...["U"]; docs "Any T is treated as a U").
    [Test]
    public void GetRnaComplementBase_Thymine_ReturnsAdenine()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceExtensions.GetRnaComplementBase('T'), Is.EqualTo('A'), "T is treated as U in RNA context; complement of U is A");
            Assert.That(SequenceExtensions.GetRnaComplementBase('t'), Is.EqualTo('A'), "lowercase t is also treated as U");
        });
    }

    // M3 — All eleven IUPAC ambiguity codes (Biopython ambiguous_rna_complement; bioinformatics.org IUPAC table).
    [Test]
    public void GetRnaComplementBase_IupacAmbiguityCodesUppercase_ReturnsRnaComplements()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceExtensions.GetRnaComplementBase('R'), Is.EqualTo('Y'), "R (A|G) -> Y (C|U)");
            Assert.That(SequenceExtensions.GetRnaComplementBase('Y'), Is.EqualTo('R'), "Y (C|U) -> R (A|G)");
            Assert.That(SequenceExtensions.GetRnaComplementBase('S'), Is.EqualTo('S'), "S (G|C) is self-complementary");
            Assert.That(SequenceExtensions.GetRnaComplementBase('W'), Is.EqualTo('W'), "W (A|U) is self-complementary");
            Assert.That(SequenceExtensions.GetRnaComplementBase('K'), Is.EqualTo('M'), "K (G|U) -> M (A|C)");
            Assert.That(SequenceExtensions.GetRnaComplementBase('M'), Is.EqualTo('K'), "M (A|C) -> K (G|U)");
            Assert.That(SequenceExtensions.GetRnaComplementBase('B'), Is.EqualTo('V'), "B (C|G|U) -> V (A|C|G)");
            Assert.That(SequenceExtensions.GetRnaComplementBase('V'), Is.EqualTo('B'), "V (A|C|G) -> B (C|G|U)");
            Assert.That(SequenceExtensions.GetRnaComplementBase('D'), Is.EqualTo('H'), "D (A|G|U) -> H (A|C|U)");
            Assert.That(SequenceExtensions.GetRnaComplementBase('H'), Is.EqualTo('D'), "H (A|C|U) -> D (A|G|U)");
            Assert.That(SequenceExtensions.GetRnaComplementBase('N'), Is.EqualTo('N'), "N (any) is self-complementary");
        });
    }

    // M4 — Lowercase recognized input returns uppercase complement (repo convention; mapping per ambiguous_rna_complement).
    [Test]
    public void GetRnaComplementBase_LowercaseRecognized_ReturnsUppercase()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceExtensions.GetRnaComplementBase('a'), Is.EqualTo('U'), "a -> U (uppercased)");
            Assert.That(SequenceExtensions.GetRnaComplementBase('u'), Is.EqualTo('A'), "u -> A");
            Assert.That(SequenceExtensions.GetRnaComplementBase('c'), Is.EqualTo('G'), "c -> G");
            Assert.That(SequenceExtensions.GetRnaComplementBase('g'), Is.EqualTo('C'), "g -> C");
            Assert.That(SequenceExtensions.GetRnaComplementBase('r'), Is.EqualTo('Y'), "r -> Y");
            Assert.That(SequenceExtensions.GetRnaComplementBase('y'), Is.EqualTo('R'), "y -> R");
            Assert.That(SequenceExtensions.GetRnaComplementBase('n'), Is.EqualTo('N'), "n -> N");
        });
    }

    // M5 — Full-alphabet worked example from Biopython reverse_complement_rna("ACGTUacgtuXYZxyz") = "zrxZRXaacguAACGU".
    // Un-reversed forward complement (Biopython, case-preserved) = "UGCAAugcaaXRZxrz". Under the repo uppercase
    // convention, recognized bases (incl. Y/y->R) are uppercased; non-IUPAC chars (X, Z, x, z) pass through verbatim.
    [Test]
    public void GetRnaComplementBase_BiopythonAlphabetExample_MatchesPerBase()
    {
        const string input = "ACGTUacgtuXYZxyz";
        const string expected = "UGCAAUGCAAXRZxRz"; // repo convention applied to Biopython "UGCAAugcaaXRZxrz"

        var actual = new char[input.Length];
        for (int i = 0; i < input.Length; i++)
            actual[i] = SequenceExtensions.GetRnaComplementBase(input[i]);

        Assert.That(new string(actual), Is.EqualTo(expected),
            "Per-base RNA complement of the Biopython worked-example alphabet under repo uppercase convention");
    }

    // M6 — RNA-specific: A -> U, distinct from the DNA complement A -> T (complement_rna("ACG")="UGC" vs complement="TGC").
    [Test]
    public void GetRnaComplementBase_Adenine_DiffersFromDnaComplement()
    {
        char rna = SequenceExtensions.GetRnaComplementBase('A');
        char dna = SequenceExtensions.GetComplementBase('A');

        Assert.Multiple(() =>
        {
            Assert.That(rna, Is.EqualTo('U'), "RNA complement of A is U");
            Assert.That(dna, Is.EqualTo('T'), "DNA complement of A is T");
            Assert.That(rna, Is.Not.EqualTo(dna), "RNA and DNA complements of A must differ (U vs T)");
        });
    }

    // S1 — Non-IUPAC characters pass through unchanged (not in the nucleotide table; gaps map to gaps).
    [Test]
    public void GetRnaComplementBase_NonIupacCharacters_PassThroughUnchanged()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceExtensions.GetRnaComplementBase('-'), Is.EqualTo('-'), "gap '-' passes through");
            Assert.That(SequenceExtensions.GetRnaComplementBase('.'), Is.EqualTo('.'), "gap '.' passes through");
            Assert.That(SequenceExtensions.GetRnaComplementBase('5'), Is.EqualTo('5'), "digit passes through");
            Assert.That(SequenceExtensions.GetRnaComplementBase('Z'), Is.EqualTo('Z'), "non-nucleotide letter Z passes through");
            Assert.That(SequenceExtensions.GetRnaComplementBase('z'), Is.EqualTo('z'), "lowercase z passes through with original case");
            Assert.That(SequenceExtensions.GetRnaComplementBase(' '), Is.EqualTo(' '), "space passes through");
        });
    }

    // S2 — Involution on the canonical RNA bases and ambiguity codes: complement(complement(x)) == x (INV-06).
    [Test]
    public void GetRnaComplementBase_Involution_AllRnaBasesAndCodes()
    {
        const string alphabet = "AUCGRYSWKMBDHVN";
        Assert.Multiple(() =>
        {
            foreach (char x in alphabet)
            {
                char once = SequenceExtensions.GetRnaComplementBase(x);
                char twice = SequenceExtensions.GetRnaComplementBase(once);
                Assert.That(twice, Is.EqualTo(x), $"complement is an involution for '{x}' (got '{once}' then '{twice}')");
            }
        });
    }

    // C1 — RNA alphabet: no recognized base or ambiguity code ever emits 'T' (INV-02).
    [Test]
    public void GetRnaComplementBase_RecognizedBases_NeverEmitThymine()
    {
        const string recognized = "AUCGTRYSWKMBDHVN";
        Assert.Multiple(() =>
        {
            foreach (char x in recognized)
                Assert.That(SequenceExtensions.GetRnaComplementBase(x), Is.Not.EqualTo('T'),
                    $"RNA complement must use U, never T (input '{x}')");
        });
    }

    #endregion
}
