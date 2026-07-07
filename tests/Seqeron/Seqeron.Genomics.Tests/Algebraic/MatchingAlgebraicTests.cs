using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Algebraic;

/// <summary>
/// Algebraic-law tests for the Matching area (approximate-matching metrics).
///
/// Algebraic testing pins the formal laws an operation must obey for EVERY input.
/// For string-distance functions the relevant algebra is the *metric-space axioms*:
/// identity of indiscernibles d(x,x)=0, symmetry d(a,b)=d(b,a), and the triangle
/// inequality d(a,c) ≤ d(a,b)+d(b,c). A function that violates any axiom is not a
/// metric and silently breaks every downstream algorithm (clustering, nearest-
/// neighbour, phylogeny) that assumes one.
/// — docs/checklists/06_ALGEBRAIC_TESTING.md §Description, rows 9–10;
///   docs/ADVANCED_TESTING_CHECKLIST.md §4 "Algebraic Testing".
/// </summary>
[TestFixture]
[Category("Algebraic")]
[Category("Matching")]
public class MatchingAlgebraicTests
{
    /// <summary>
    /// Generates three DNA strings of the SAME length, so equal-length metrics
    /// (Hamming) are defined on every triple and the triangle inequality is
    /// meaningful.
    /// </summary>
    private static Arbitrary<(string A, string B, string C)> EqualLengthTriple() =>
        (from n in Gen.Choose(1, 40)
         from a in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(n)
         from b in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(n)
         from c in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(n)
         select (new string(a), new string(b), new string(c)))
        .ToArbitrary();

    /// <summary>
    /// Generates three DNA strings of INDEPENDENT lengths (including empty), the
    /// natural domain of edit distance, which is defined for unequal lengths.
    /// </summary>
    private static Arbitrary<(string A, string B, string C)> FreeLengthTriple() =>
        (from a in Gen.Elements('A', 'C', 'G', 'T').ArrayOf()
         from b in Gen.Elements('A', 'C', 'G', 'T').ArrayOf()
         from c in Gen.Elements('A', 'C', 'G', 'T').ArrayOf()
         select (new string(a), new string(b), new string(c)))
        .ToArbitrary();

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PAT-APPROX-001 — Hamming distance (Matching)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 9.
    //
    // Model: Hamming distance = number of positions at which two equal-length
    //        strings differ (case-insensitive). It is a metric on the space of
    //        fixed-length strings over the alphabet.
    //   — ApproximateMatcher.HammingDistance; Wikipedia "Hamming distance".
    //
    // Laws under test (checklist row 9):
    //   • ID   — d(x, x) = 0 (identity of indiscernibles; reflexive zero).
    //   • COMM — d(a, b) = d(b, a) (symmetry; counting mismatches is order-free).
    //   • TRI  — d(a, c) ≤ d(a, b) + d(b, c) (triangle inequality; a position
    //            that differs between a and c must differ from b in at least one
    //            of the two comparisons).
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>
    /// ID: d(x, x) = 0 — a sequence is at zero Hamming distance from itself.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Hamming_Identity_SelfDistanceIsZero()
    {
        return Prop.ForAll(EqualLengthTriple(), t =>
            (ApproximateMatcher.HammingDistance(t.A, t.A) == 0)
                .Label($"d(x,x) != 0 for \"{t.A}\""));
    }

    /// <summary>
    /// COMM: d(a, b) = d(b, a) — Hamming distance is symmetric.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Hamming_Commutative_Symmetric()
    {
        return Prop.ForAll(EqualLengthTriple(), t =>
        {
            int ab = ApproximateMatcher.HammingDistance(t.A, t.B);
            int ba = ApproximateMatcher.HammingDistance(t.B, t.A);
            return (ab == ba).Label($"d(a,b)={ab} != d(b,a)={ba}");
        });
    }

    /// <summary>
    /// TRI: d(a, c) ≤ d(a, b) + d(b, c) — the triangle inequality.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Hamming_TriangleInequality()
    {
        return Prop.ForAll(EqualLengthTriple(), t =>
        {
            int ac = ApproximateMatcher.HammingDistance(t.A, t.C);
            int ab = ApproximateMatcher.HammingDistance(t.A, t.B);
            int bc = ApproximateMatcher.HammingDistance(t.B, t.C);
            return (ac <= ab + bc).Label($"d(a,c)={ac} > d(a,b)+d(b,c)={ab}+{bc}");
        });
    }

    /// <summary>
    /// Identity of indiscernibles (strict): d(a, b) = 0 ⟺ a = b.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Hamming_ZeroIffEqual()
    {
        return Prop.ForAll(EqualLengthTriple(), t =>
        {
            int ab = ApproximateMatcher.HammingDistance(t.A, t.B);
            bool equal = string.Equals(t.A, t.B, StringComparison.OrdinalIgnoreCase);
            return ((ab == 0) == equal).Label($"d={ab}, equal={equal} for a=\"{t.A}\" b=\"{t.B}\"");
        });
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PAT-APPROX-002 — Edit (Levenshtein) distance (Matching)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 10.
    //
    // Model: minimum number of single-character insertions, deletions and
    //        substitutions transforming one string into the other. A metric on
    //        the free monoid over the alphabet (unequal lengths allowed).
    //   — ApproximateMatcher.EditDistance; Wikipedia "Levenshtein distance".
    //
    // Laws under test (checklist row 10):
    //   • ID   — d(x, x) = 0 (no edits needed; and d=0 ⟺ x=y).
    //   • COMM — d(a, b) = d(b, a) (each edit operation has an inverse of equal
    //            unit cost, so the optimal cost is direction-independent).
    //   • TRI  — d(a, c) ≤ d(a, b) + d(b, c) (concatenating the edit script
    //            a→b with b→c yields a valid, not-necessarily-optimal, a→c script).
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>
    /// ID: d(x, x) = 0, and more strictly d(a, b) = 0 ⟺ a = b.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Edit_Identity_ZeroIffEqual()
    {
        return Prop.ForAll(FreeLengthTriple(), t =>
        {
            int self = ApproximateMatcher.EditDistance(t.A, t.A);
            int ab = ApproximateMatcher.EditDistance(t.A, t.B);
            bool equal = t.A == t.B;
            return (self == 0 && (ab == 0) == equal)
                .Label($"self={self}, d(a,b)={ab}, equal={equal}");
        });
    }

    /// <summary>
    /// COMM: d(a, b) = d(b, a) — edit distance is symmetric.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Edit_Commutative_Symmetric()
    {
        return Prop.ForAll(FreeLengthTriple(), t =>
        {
            int ab = ApproximateMatcher.EditDistance(t.A, t.B);
            int ba = ApproximateMatcher.EditDistance(t.B, t.A);
            return (ab == ba).Label($"d(a,b)={ab} != d(b,a)={ba}");
        });
    }

    /// <summary>
    /// TRI: d(a, c) ≤ d(a, b) + d(b, c) — the triangle inequality.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Edit_TriangleInequality()
    {
        return Prop.ForAll(FreeLengthTriple(), t =>
        {
            int ac = ApproximateMatcher.EditDistance(t.A, t.C);
            int ab = ApproximateMatcher.EditDistance(t.A, t.B);
            int bc = ApproximateMatcher.EditDistance(t.B, t.C);
            return (ac <= ab + bc).Label($"d(a,c)={ac} > d(a,b)+d(b,c)={ab}+{bc}");
        });
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PAT-IUPAC-001 — IUPAC ambiguity-code matching (Matching)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 11.
    //
    // Model: each IUPAC code denotes a SET of concrete bases; matching is set
    //        membership: MatchesIupac(b, code) ⟺ b ∈ set(code). The codes form a
    //        bounded lattice under set inclusion with N (= {A,C,G,T}) as the top.
    //   — docs IUPAC-IUB 1970 / NC-IUB 1984; IupacHelper.MatchesIupac.
    //
    // Laws under test (checklist row 11):
    //   • ID   — exact base match is reflexive: each concrete base used as its own
    //            code matches only itself (the singleton atoms of the lattice).
    //   • DIST — N is the top element: every base matches N, and set(N) ⊇ set(code)
    //            for every code, i.e. anything any code matches, N matches too.
    //   • DIST — complementary two-way codes partition the alphabet: R⊎Y, S⊎W,
    //            K⊎M each cover {A,C,G,T} exactly once (conservation); the
    //            three-way codes B,D,H,V are exactly "all bases but one".
    // ═══════════════════════════════════════════════════════════════════════

    private static readonly char[] Bases = { 'A', 'C', 'G', 'T' };

    /// <summary>
    /// ID: exact base match is reflexive and discriminating — base b used as a
    /// code matches b and rejects every other base (the singleton atoms).
    /// </summary>
    [Test]
    public void Iupac_Identity_ExactMatchIsReflexiveAtom()
    {
        foreach (char code in Bases)
        {
            foreach (char b in Bases)
            {
                IupacHelper.MatchesIupac(b, code).Should().Be(b == code,
                    $"base '{b}' vs exact code '{code}'");
            }
        }
    }

    /// <summary>
    /// DIST: N is the top element — it matches every concrete base, and its match
    /// set is a superset of every other code's match set.
    /// </summary>
    [Test]
    public void Iupac_Distributive_NIsTopElement()
    {
        // N matches every base.
        foreach (char b in Bases)
            IupacHelper.MatchesIupac(b, 'N').Should().BeTrue($"N must match '{b}'");

        // set(N) ⊇ set(code) for every standard code.
        const string codes = "ACGTNRYSWKMBDHV";
        foreach (char code in codes)
        {
            foreach (char b in Bases)
            {
                if (IupacHelper.MatchesIupac(b, code))
                    IupacHelper.MatchesIupac(b, 'N').Should().BeTrue(
                        $"'{b}' matched '{code}' but not top element N");
            }
        }
    }

    /// <summary>
    /// DIST conservation: the three complementary two-way codes each partition the
    /// alphabet — every base matches exactly one of {R,Y}, one of {S,W}, one of {K,M}.
    /// </summary>
    [Test]
    public void Iupac_Distributive_TwoWayCodesPartitionAlphabet()
    {
        foreach (char b in Bases)
        {
            (IupacHelper.MatchesIupac(b, 'R') ^ IupacHelper.MatchesIupac(b, 'Y'))
                .Should().BeTrue($"R/Y must partition '{b}'");
            (IupacHelper.MatchesIupac(b, 'S') ^ IupacHelper.MatchesIupac(b, 'W'))
                .Should().BeTrue($"S/W must partition '{b}'");
            (IupacHelper.MatchesIupac(b, 'K') ^ IupacHelper.MatchesIupac(b, 'M'))
                .Should().BeTrue($"K/M must partition '{b}'");
        }
    }

    /// <summary>
    /// DIST: the three-way codes B,D,H,V are exactly the complements of the four
    /// singleton atoms — B = ¬A, D = ¬C, H = ¬G, V = ¬T.
    /// </summary>
    [Test]
    public void Iupac_Distributive_ThreeWayCodesAreSingletonComplements()
    {
        var excludedBy = new (char Code, char Excluded)[] { ('B', 'A'), ('D', 'C'), ('H', 'G'), ('V', 'T') };
        foreach (var (code, excluded) in excludedBy)
        {
            foreach (char b in Bases)
            {
                IupacHelper.MatchesIupac(b, code).Should().Be(b != excluded,
                    $"code '{code}' should match '{b}' iff '{b}' != '{excluded}'");
            }
        }
    }
}
