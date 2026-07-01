// RNA-PAIR-001 — RNA Base Pairing
// Evidence: docs/Evidence/RNA-PAIR-001-Evidence.md
// TestSpec: tests/TestSpecs/RNA-PAIR-001.md
// Source: Crick FHC (1966) J Mol Biol 19(2):548-555 (G-U wobble); Wikipedia Base pair (A-U/G-C
//         canonical Watson-Crick); IUPAC-IUB (1970) Biochemistry 9(20):4022-4027 + Biopython
//         complement_rna (RNA complement A->U, U->A, G->C, C->G, T->A).

using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Analysis.RnaSecondaryStructure;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class RnaSecondaryStructure_CanPair_Tests
{
    #region CanPair

    // M1 — Canonical Watson-Crick pairs. Evidence: Wikipedia Base pair — canonical RNA pairs are A•U and G•C.
    [Test]
    public void CanPair_WatsonCrickPairs_ReturnsTrue()
    {
        Assert.Multiple(() =>
        {
            Assert.That(CanPair('A', 'U'), Is.True, "A-U is a canonical Watson-Crick pair (2 H-bonds)");
            Assert.That(CanPair('U', 'A'), Is.True, "U-A is the same Watson-Crick pair as A-U");
            Assert.That(CanPair('G', 'C'), Is.True, "G-C is a canonical Watson-Crick pair (3 H-bonds)");
            Assert.That(CanPair('C', 'G'), Is.True, "C-G is the same Watson-Crick pair as G-C");
        });
    }

    // M2 — G:U wobble pairs. Evidence: Crick (1966) wobble hypothesis; G-U is the standard RNA wobble pair.
    [Test]
    public void CanPair_WobblePairs_ReturnsTrue()
    {
        Assert.Multiple(() =>
        {
            Assert.That(CanPair('G', 'U'), Is.True, "G-U is the standard RNA wobble pair (Crick 1966)");
            Assert.That(CanPair('U', 'G'), Is.True, "U-G is the same wobble pair as G-U (Crick 1966)");
        });
    }

    // M3 — Non-pairing combinations. Evidence: Crick (1966)/pairing rules — A pairs only with U, C only with G.
    [Test]
    public void CanPair_NonPairs_ReturnsFalse()
    {
        Assert.Multiple(() =>
        {
            Assert.That(CanPair('A', 'A'), Is.False, "A-A does not pair");
            Assert.That(CanPair('A', 'G'), Is.False, "A pairs only with U, not G");
            Assert.That(CanPair('A', 'C'), Is.False, "A pairs only with U, not C");
            Assert.That(CanPair('C', 'U'), Is.False, "C pairs only with G, not U");
            Assert.That(CanPair('G', 'G'), Is.False, "G-G does not pair");
            Assert.That(CanPair('C', 'C'), Is.False, "C-C does not pair");
        });
    }

    // S1 — Case-insensitivity. Evidence: normalization contract; lower/upper denote the same nucleotide.
    [Test]
    public void CanPair_LowercaseInput_SameAsUppercase()
    {
        Assert.Multiple(() =>
        {
            Assert.That(CanPair('a', 'u'), Is.True, "lowercase a-u must pair like A-U");
            Assert.That(CanPair('g', 'c'), Is.True, "lowercase g-c must pair like G-C");
            Assert.That(CanPair('g', 'u'), Is.True, "lowercase g-u must pair like G-U wobble");
            Assert.That(CanPair('a', 'g'), Is.False, "lowercase a-g must not pair");
        });
    }

    // S2 — RNA alphabet only: T is not an RNA base, so CanPair does not pair it.
    // Evidence: Crick (1966)/Watson-Crick define pairing over the RNA alphabet {A,C,G,U}; the
    // sources do not define a DNA base T in RNA pairing, so CanPair returns false for T inputs.
    [Test]
    public void CanPair_DnaT_NotAnRnaBase_ReturnsFalse()
    {
        Assert.Multiple(() =>
        {
            Assert.That(CanPair('T', 'A'), Is.False, "T is not an RNA base; RNA pairing is defined over {A,C,G,U}");
            Assert.That(CanPair('A', 'T'), Is.False, "A pairs with U (RNA), not T");
            Assert.That(CanPair('G', 'T'), Is.False, "G wobble-pairs with U (RNA), not T");
        });
    }

    // C1 — Out-of-range characters. Evidence: robustness — bounds-checked lookup returns false, no exception.
    [Test]
    public void CanPair_OutOfRangeChar_ReturnsFalse()
    {
        Assert.Multiple(() =>
        {
            Assert.That(CanPair('A', 'Ѐ'), Is.False, "non-ASCII base does not pair and must not throw");
            Assert.That(CanPair('5', 'U'), Is.False, "digit is not a base and does not pair");
        });
    }

    #endregion

    #region GetBasePairType

    // M4 — Watson-Crick type. Evidence: Wikipedia Base pair — A•U and G•C are the canonical WC pairs.
    [Test]
    public void GetBasePairType_WatsonCrick_ReturnsWatsonCrick()
    {
        Assert.Multiple(() =>
        {
            Assert.That(GetBasePairType('A', 'U'), Is.EqualTo(BasePairType.WatsonCrick), "A-U is Watson-Crick");
            Assert.That(GetBasePairType('U', 'A'), Is.EqualTo(BasePairType.WatsonCrick), "U-A is Watson-Crick");
            Assert.That(GetBasePairType('G', 'C'), Is.EqualTo(BasePairType.WatsonCrick), "G-C is Watson-Crick");
            Assert.That(GetBasePairType('C', 'G'), Is.EqualTo(BasePairType.WatsonCrick), "C-G is Watson-Crick");
        });
    }

    // M5 — Wobble type (INV-04). Evidence: Crick (1966); wobble does NOT follow Watson-Crick rules.
    [Test]
    public void GetBasePairType_Wobble_ReturnsWobble()
    {
        Assert.Multiple(() =>
        {
            Assert.That(GetBasePairType('G', 'U'), Is.EqualTo(BasePairType.Wobble), "G-U is Wobble, not Watson-Crick (Crick 1966)");
            Assert.That(GetBasePairType('U', 'G'), Is.EqualTo(BasePairType.Wobble), "U-G is Wobble, not Watson-Crick (Crick 1966)");
        });
    }

    // M6 — Non-pairs return null. Evidence: pairing rules — A only with U, C only with G.
    [Test]
    public void GetBasePairType_NonPairs_ReturnsNull()
    {
        Assert.Multiple(() =>
        {
            Assert.That(GetBasePairType('A', 'A'), Is.Null, "A-A is not a pair");
            Assert.That(GetBasePairType('A', 'G'), Is.Null, "A-G is not a pair");
            Assert.That(GetBasePairType('C', 'U'), Is.Null, "C-U is not a pair");
        });
    }

    #endregion

    #region GetComplement

    // M7 — RNA complement of standard bases. Evidence: IUPAC-IUB (1970) + Biopython complement_rna
    //      (CGAUT -> GCUAA, i.e. C->G, G->C, A->U, U->A, T->A).
    [Test]
    public void GetComplement_StandardBases_ReturnsRnaComplement()
    {
        Assert.Multiple(() =>
        {
            Assert.That(GetComplement('A'), Is.EqualTo('U'), "complement of A is U in RNA");
            Assert.That(GetComplement('U'), Is.EqualTo('A'), "complement of U is A");
            Assert.That(GetComplement('G'), Is.EqualTo('C'), "complement of G is C");
            Assert.That(GetComplement('C'), Is.EqualTo('G'), "complement of C is G");
            Assert.That(GetComplement('T'), Is.EqualTo('A'), "T treated as U; complement is A (Biopython complement_rna)");
        });
    }

    // S3 — IUPAC degenerate complements. Evidence: IUPAC-IUB (1970) Biochemistry 9(20):4022-4027,
    //      via Wikipedia "Nucleic acid notation" full complement table:
    //      W<->W, S<->S, M<->K, K<->M, R<->Y, Y<->R, B<->V, V<->B, D<->H, H<->D, N<->N.
    //      All 11 degenerate branches of GetRnaComplementBase are exercised against the sourced table.
    [Test]
    public void GetComplement_IupacDegenerate_ReturnsExpected()
    {
        Assert.Multiple(() =>
        {
            Assert.That(GetComplement('N'), Is.EqualTo('N'), "N (any) complements to N");
            Assert.That(GetComplement('R'), Is.EqualTo('Y'), "R (A|G purine) complements to Y (C|U pyrimidine)");
            Assert.That(GetComplement('Y'), Is.EqualTo('R'), "Y (C|U) complements to R (A|G)");
            Assert.That(GetComplement('W'), Is.EqualTo('W'), "W (weak, A|U) is self-complementary");
            Assert.That(GetComplement('S'), Is.EqualTo('S'), "S (strong, G|C) is self-complementary");
            Assert.That(GetComplement('M'), Is.EqualTo('K'), "M (amino, A|C) complements to K (keto, G|U)");
            Assert.That(GetComplement('K'), Is.EqualTo('M'), "K (keto, G|U) complements to M (amino, A|C)");
            Assert.That(GetComplement('B'), Is.EqualTo('V'), "B (not A, C|G|U) complements to V (not U, A|C|G)");
            Assert.That(GetComplement('V'), Is.EqualTo('B'), "V (not U, A|C|G) complements to B (not A, C|G|U)");
            Assert.That(GetComplement('D'), Is.EqualTo('H'), "D (not C, A|G|U) complements to H (not G, A|C|U)");
            Assert.That(GetComplement('H'), Is.EqualTo('D'), "H (not G, A|C|U) complements to D (not C, A|G|U)");
        });
    }

    // S3b — Lowercase IUPAC degenerate complements preserve case mapping (normalization contract).
    // Evidence: same IUPAC-IUB table; GetRnaComplementBase accepts lowercase via 'X' or 'x' arms.
    [Test]
    public void GetComplement_LowercaseIupacDegenerate_ReturnsExpected()
    {
        Assert.Multiple(() =>
        {
            Assert.That(GetComplement('r'), Is.EqualTo('Y'), "lowercase r complements to Y");
            Assert.That(GetComplement('n'), Is.EqualTo('N'), "lowercase n complements to N");
            Assert.That(GetComplement('k'), Is.EqualTo('M'), "lowercase k complements to M");
        });
    }

    // S3c — Non-IUPAC characters pass through unchanged (Core helper contract).
    // Evidence: doc 6.1 "Non-IUPAC char in GetComplement -> passed through unchanged".
    [Test]
    public void GetComplement_NonIupacChar_PassesThroughUnchanged()
    {
        Assert.Multiple(() =>
        {
            Assert.That(GetComplement('-'), Is.EqualTo('-'), "gap character passes through unchanged");
            Assert.That(GetComplement('X'), Is.EqualTo('X'), "non-IUPAC 'X' passes through unchanged");
        });
    }

    #endregion

    #region Invariants (property-based)

    // M8 — Symmetry (INV-01, INV-02): pairing is reciprocal. Evidence: Wikipedia Base pair — a pair is between
    //      two complementary bases; A•U ≡ U•A. Verified over the full RNA alphabet cross-product.
    [Test]
    public void CanPair_And_Type_AreSymmetric()
    {
        const string alphabet = "ACGU";
        Assert.Multiple(() =>
        {
            foreach (char x in alphabet)
            {
                foreach (char y in alphabet)
                {
                    Assert.That(CanPair(x, y), Is.EqualTo(CanPair(y, x)),
                        $"CanPair must be symmetric for ({x},{y})");
                    Assert.That(GetBasePairType(x, y), Is.EqualTo(GetBasePairType(y, x)),
                        $"GetBasePairType must be symmetric for ({x},{y})");
                }
            }
        });
    }

    // M9 — Consistency (INV-03): CanPair is true iff GetBasePairType is non-null.
    //      Evidence: a pair exists iff it has a defined type (shared classification). Full cross-product.
    [Test]
    public void CanPair_AgreesWith_GetBasePairType()
    {
        const string alphabet = "ACGU";
        Assert.Multiple(() =>
        {
            foreach (char x in alphabet)
            {
                foreach (char y in alphabet)
                {
                    bool canPair = CanPair(x, y);
                    bool hasType = GetBasePairType(x, y) is not null;
                    Assert.That(canPair, Is.EqualTo(hasType),
                        $"CanPair and (type != null) must agree for ({x},{y})");
                }
            }
        });
    }

    #endregion
}
