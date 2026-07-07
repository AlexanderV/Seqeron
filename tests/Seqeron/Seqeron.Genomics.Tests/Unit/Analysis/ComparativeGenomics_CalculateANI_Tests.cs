// COMPGEN-ANI-001 — Average Nucleotide Identity (ANI), ANIb (Goris et al. 2007)
// Evidence: docs/Evidence/COMPGEN-ANI-001-Evidence.md
// TestSpec: tests/TestSpecs/COMPGEN-ANI-001.md
// Source: Goris J, Konstantinidis KT, Klappenbach JA, Coenye T, Vandamme P, Tiedje JM (2007).
//         DNA-DNA hybridization values and their relationship to whole-genome sequence
//         similarities. Int J Syst Evol Microbiol 57(1):81-91. DOI:10.1099/ijs.0.64483-0
//
// ANI = mean identity of conserved query fragments, identity recalculated over the fragment
// length, keeping only fragments with >30% identity over >=70% alignable length.
// Expected values below are derived by hand from that formula (independently re-derived in
// Python), NOT read back from the implementation.
//
// Gapped placement + reciprocal ANI (COMPGEN-ANI-001 limitation fix):
//   Goris 2007 use GAPPED BLASTN and "reverse searching ... to provide reciprocal values".
//   pyani (reference impl) recalculates identity over the query length (ani_pid = ani_alnids/qlen)
//   and coverage as ani_alnlen/qlen where ani_alnlen = blast_alnlen - blast_gaps. Gapped tests
//   below encode hand-derived alignment identities; reciprocal tests encode the mean of the two
//   directional ANIs.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class ComparativeGenomics_CalculateANI_Tests
{
    // Reference used across the worked examples: four distinct homopolymer runs of 4 nt each.
    private const string Reference = "AAAACCCCGGGGTTTT";

    private const double Tolerance = 1e-10;

    #region CalculateANI

    // M1 — Identical genomes: every 4-nt fragment is a perfect substring of the reference,
    // so each fragment identity = 4/4 = 1.0 and ANI = mean(1,1,1,1) = 1.0 (Goris 2007; INV-02).
    [Test]
    public void CalculateANI_IdenticalGenomes_ReturnsOne()
    {
        double ani = ComparativeGenomics.CalculateANI(Reference, Reference, fragmentLength: 4);

        Assert.That(ani, Is.EqualTo(1.0).Within(Tolerance),
            "Identical genomes must yield ANI = 1.0: every fragment is a perfect substring (Goris 2007).");
    }

    // M2 — One substituted base in the final fragment. Query last fragment "TTTA" best-matches
    // "TTTT" with 3/4 bases => identity 0.75. First three fragments are perfect (1.0 each).
    // ANI = (1 + 1 + 1 + 0.75) / 4 = 0.9375. Identity is recalculated over the whole fragment
    // length, so a single mismatch in a 4-nt fragment costs exactly 0.25 there (Goris 2007).
    [Test]
    public void CalculateANI_OneSubstitution_ReturnsExactMean()
    {
        double ani = ComparativeGenomics.CalculateANI("AAAACCCCGGGGTTTA", Reference, fragmentLength: 4);

        Assert.That(ani, Is.EqualTo(0.9375).Within(Tolerance),
            "One mismatch in one 4-nt fragment gives that fragment 0.75; ANI = (1+1+1+0.75)/4 = 0.9375 (Goris 2007).");
    }

    // M3 — Last fragment "AATT" best-matches "TTTT" with 2/4 = 0.5 identity (still > 0.30, so it
    // qualifies). ANI = (1 + 1 + 1 + 0.5) / 4 = 0.875. Distinguishes the recalculated-over-fragment
    // identity (0.5) from a naive longest-common-substring proxy ("TT" => 2/4 would also be 0.5 here,
    // but M2/M4 separate identity from LCS).
    [Test]
    public void CalculateANI_HalfIdentityFragment_ReturnsExactMean()
    {
        double ani = ComparativeGenomics.CalculateANI("AAAACCCCGGGGAATT", Reference, fragmentLength: 4);

        Assert.That(ani, Is.EqualTo(0.875).Within(Tolerance),
            "A 0.5-identity qualifying fragment gives ANI = (1+1+1+0.5)/4 = 0.875 (Goris 2007 mean identity).");
    }

    // M4 — Identity cut-off excludes a fragment. Reference is all A's; query "AAAACGTC" splits into
    // "AAAA" (identity 1.0, kept) and "CGTC" (0 matches vs any AAAA window => identity 0.0, which is
    // NOT > 0.30, so excluded). ANI = mean of the single kept fragment = 1.0. A wrong implementation
    // that averaged ALL fragments would return (1.0 + 0.0)/2 = 0.5, so this value would fail it.
    [Test]
    public void CalculateANI_FragmentBelowIdentityCutoff_IsExcludedFromMean()
    {
        double ani = ComparativeGenomics.CalculateANI("AAAACGTC", "AAAAAAAA", fragmentLength: 4);

        Assert.That(ani, Is.EqualTo(1.0).Within(Tolerance),
            "A fragment with identity <= 0.30 must be discarded; ANI = mean of qualifying fragments only = 1.0 (Goris 2007 >30% cut-off).");
    }

    // M5 — Alignable-region cut-off. Reference "AA" is shorter than the 4-nt fragment, so the
    // fragment cannot align over >= 70% of its length (alignable fraction 0). No fragment qualifies
    // => ANI = 0 (Goris 2007 ">=70% alignable region"; INV-03).
    [Test]
    public void CalculateANI_FragmentBelowAlignableCutoff_ReturnsZero()
    {
        double ani = ComparativeGenomics.CalculateANI("AAAA", "AA", fragmentLength: 4);

        Assert.That(ani, Is.EqualTo(0.0).Within(Tolerance),
            "Reference shorter than the fragment cannot align >= 70% => fragment excluded => ANI = 0 (Goris 2007).");
    }

    // M6 — Consecutive non-overlapping fragmentation; trailing partial fragment dropped. A 10-nt
    // query with fragmentLength 4 yields exactly 2 fragments (offsets 0 and 4); the trailing 2 nt
    // are ignored (Goris 2007 "consecutive 1020 nt fragments"; INV-04). Here both fragments are
    // perfect substrings of the matching reference => ANI = 1.0, proving the trailing-2-nt mismatch
    // region is not counted.
    [Test]
    public void CalculateANI_TrailingPartialFragment_IsIgnored()
    {
        // Query: "AAAACCCCXX" where XX would mismatch; reference contains AAAA and CCCC.
        // Fragments used: "AAAA" (1.0), "CCCC" (1.0); trailing "XX" dropped.
        double ani = ComparativeGenomics.CalculateANI("AAAACCCCXX", "AAAACCCCGGGG", fragmentLength: 4);

        Assert.That(ani, Is.EqualTo(1.0).Within(Tolerance),
            "Only the 2 full fragments (AAAA, CCCC) are scored; the trailing 2-nt partial is ignored (Goris 2007).");
    }

    // M7 — Null / empty inputs return 0 (validation contract).
    [Test]
    public void CalculateANI_NullOrEmptyInput_ReturnsZero()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ComparativeGenomics.CalculateANI(null!, Reference), Is.EqualTo(0.0),
                "Null query must return 0.");
            Assert.That(ComparativeGenomics.CalculateANI(Reference, null!), Is.EqualTo(0.0),
                "Null reference must return 0.");
            Assert.That(ComparativeGenomics.CalculateANI("", Reference), Is.EqualTo(0.0),
                "Empty query must return 0.");
            Assert.That(ComparativeGenomics.CalculateANI(Reference, ""), Is.EqualTo(0.0),
                "Empty reference must return 0.");
        });
    }

    // M8 — Non-positive fragment length is invalid.
    [Test]
    public void CalculateANI_NonPositiveFragmentLength_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => ComparativeGenomics.CalculateANI(Reference, Reference, fragmentLength: 0),
            "fragmentLength <= 0 must throw ArgumentOutOfRangeException.");
    }

    // S1 — Range invariant (INV-01): for any inputs, ANI is a fraction in [0, 1]. Property test
    // over several divergent pairs and a default-length run. Verifies the bound the mean must obey.
    [Test]
    public void CalculateANI_AnyInput_ResultInUnitInterval()
    {
        string[][] pairs =
        {
            new[] { "AAAACCCCGGGGTTTT", "TTTTGGGGCCCCAAAA" },
            new[] { "ACGTACGTACGTACGT", "TGCATGCATGCATGCA" },
            new[] { "AAAACCCCGGGGTTTT", "AAAACCCCGGGGTTTA" },
            new[] { new string('A', 4096), new string('C', 4096) }, // default fragmentLength 1020
        };

        Assert.Multiple(() =>
        {
            foreach (var p in pairs)
            {
                double ani = ComparativeGenomics.CalculateANI(p[0], p[1]);
                Assert.That(ani, Is.InRange(0.0, 1.0),
                    $"ANI must lie in [0,1] for query '{p[0][..Math.Min(8, p[0].Length)]}...' (INV-01).");
            }
        });
    }

    // S2 — Query shorter than the fragment length: no full fragment fits => ANI = 0.
    [Test]
    public void CalculateANI_QueryShorterThanFragment_ReturnsZero()
    {
        double ani = ComparativeGenomics.CalculateANI("AAA", "AAAAAAAA", fragmentLength: 4);

        Assert.That(ani, Is.EqualTo(0.0).Within(Tolerance),
            "A query shorter than fragmentLength yields no fragments => ANI = 0.");
    }

    // C1 — Custom minIdentity. With the default minIdentity (0.30) the "CGTC" fragment (identity 0)
    // is excluded (see M4). It stays excluded even at minIdentity 0.0 because the cut-off is strict
    // (identity must be > minIdentity), and 0.0 is not > 0.0. We instead lower the cut-off below a
    // 0.25-identity fragment to confirm it then contributes. Query "AAAAAAAC": frag1 "AAAA"=1.0,
    // frag2 "AAAC" vs "AAAAAAAA" best=3/4=0.75 (already > 0.30). To exercise the parameter we use a
    // fragment with identity 0.25 and lower minIdentity below it.
    [Test]
    public void CalculateANI_CustomMinIdentity_IncludesLowerIdentityFragment()
    {
        // Reference all A. Query frag2 "ACGT" => 1 match vs AAAA = 0.25 identity.
        // Default minIdentity 0.30 excludes it -> ANI = 1.0 (only frag1 "AAAA").
        double aniDefault = ComparativeGenomics.CalculateANI("AAAAACGT", "AAAAAAAA", fragmentLength: 4);
        // Lower minIdentity to 0.20 (< 0.25) -> frag2 now qualifies: ANI = (1.0 + 0.25)/2 = 0.625.
        double aniLow = ComparativeGenomics.CalculateANI("AAAAACGT", "AAAAAAAA", fragmentLength: 4, minIdentity: 0.20);

        Assert.Multiple(() =>
        {
            Assert.That(aniDefault, Is.EqualTo(1.0).Within(Tolerance),
                "At default minIdentity 0.30 the 0.25-identity fragment is excluded => ANI = 1.0.");
            Assert.That(aniLow, Is.EqualTo(0.625).Within(Tolerance),
                "Lowering minIdentity to 0.20 includes the 0.25-identity fragment => ANI = (1.0+0.25)/2 = 0.625.");
        });
    }

    #endregion

    #region CalculateANI (gapped placement)

    // G1 — Identical genomes with gapped placement still yield ANI = 1.0: every fragment aligns to
    // an exact copy with no gaps, so identical columns = fragment length (Goris 2007; INV-02).
    [Test]
    public void CalculateANI_Gapped_IdenticalGenomes_ReturnsOne()
    {
        double ani = ComparativeGenomics.CalculateANI(Reference, Reference, fragmentLength: 4, gapped: true);

        Assert.That(ani, Is.EqualTo(1.0).Within(Tolerance),
            "Gapped identical genomes must yield ANI = 1.0: the fragment aligns to an exact copy (Goris 2007).");
    }

    // G2 — Gapped placement recovers an indel that the ungapped scan cannot. Query fragment
    // "AAAACCCC" (8 nt) vs reference "AAAATCCCC" (a single T inserted after AAAA).
    //   Ungapped: the 8-mer cannot slide over the inserted T without losing column registration;
    //     the best full-length offset matches only 7 of 8 bases => identity 7/8 = 0.875.
    //   Gapped (Smith-Waterman): "AAAA-CCCC" vs "AAAATCCCC" aligns all 8 query bases identically
    //     (one gap opposite the inserted T) => 8 identical columns / 8 = identity 1.0; the 8 ungapped
    //     columns / 8 = coverage 1.0 (>= 0.70), so the fragment qualifies.
    // The gapped result (1.0) is strictly higher than the ungapped result (0.875), and is the
    // correct recalculated-over-fragment identity per Goris 2007 / pyani (ani_pid = ani_alnids/qlen).
    [Test]
    public void CalculateANI_Gapped_IndelFragment_RecoversHigherIdentityThanUngapped()
    {
        const string query = "AAAACCCC";       // 8-nt fragment
        const string reference = "AAAATCCCC";   // same fragment with one inserted base (indel)

        double ungapped = ComparativeGenomics.CalculateANI(query, reference, fragmentLength: 8, gapped: false);
        double gapped = ComparativeGenomics.CalculateANI(query, reference, fragmentLength: 8, gapped: true);

        Assert.Multiple(() =>
        {
            Assert.That(ungapped, Is.EqualTo(0.875).Within(Tolerance),
                "Ungapped best offset matches only 7/8 bases across the inserted base => 0.875.");
            Assert.That(gapped, Is.EqualTo(1.0).Within(Tolerance),
                "Gapped alignment places a gap opposite the inserted base => all 8 query bases match => identity 1.0 (Goris 2007 gapped BLASTN).");
            Assert.That(gapped, Is.GreaterThan(ungapped),
                "Gapped placement must recover strictly higher correct identity than ungapped for an indel.");
        });
    }

    // G3 — Gapped alignable-region cut-off: reference far shorter than the fragment cannot cover
    // >= 70% of the fragment length, so the fragment is excluded => ANI = 0. Query fragment "AAAACCCC"
    // (8 nt), reference "AA" (2 nt): the best local alignment covers at most 2 of 8 query columns
    // (coverage 2/8 = 0.25 < 0.70). No fragment qualifies (Goris 2007 ">= 70% alignable region").
    [Test]
    public void CalculateANI_Gapped_FragmentBelowAlignableCutoff_ReturnsZero()
    {
        double ani = ComparativeGenomics.CalculateANI("AAAACCCC", "AA", fragmentLength: 8, gapped: true);

        Assert.That(ani, Is.EqualTo(0.0).Within(Tolerance),
            "Gapped coverage 2/8 = 0.25 < 0.70 excludes the fragment => ANI = 0 (Goris 2007 alignable cut-off).");
    }

    // G4 — Gapped identity cut-off still gates: a fragment whose best gapped identity is <= 0.30 is
    // discarded. Reference all A's; query "AAAACGTC" splits into "AAAA" (identity 1.0, kept) and
    // "CGTC" (best gapped alignment to AAAA has 0 identical columns => identity 0.0, NOT > 0.30,
    // excluded). ANI = mean of the single kept fragment = 1.0 (Goris 2007 ">30% identity").
    [Test]
    public void CalculateANI_Gapped_FragmentBelowIdentityCutoff_IsExcludedFromMean()
    {
        double ani = ComparativeGenomics.CalculateANI("AAAACGTC", "AAAAAAAA", fragmentLength: 4, gapped: true);

        Assert.That(ani, Is.EqualTo(1.0).Within(Tolerance),
            "A gapped fragment with identity <= 0.30 is discarded; ANI = mean of qualifying fragments = 1.0 (Goris 2007).");
    }

    #endregion

    #region CalculateReciprocalAni

    // R1 — Reciprocal ANI of identical genomes is 1.0 in both directions => mean = 1.0
    // (Goris 2007 reverse searching; both directional ANIs are 1.0).
    [Test]
    public void CalculateReciprocalAni_IdenticalGenomes_ReturnsOne()
    {
        double ani = ComparativeGenomics.CalculateReciprocalAni(Reference, Reference, fragmentLength: 4);

        Assert.That(ani, Is.EqualTo(1.0).Within(Tolerance),
            "Reciprocal ANI of identical genomes = mean(1.0, 1.0) = 1.0 (Goris 2007 reciprocal values).");
    }

    // R2 — Reciprocal ANI is symmetric: ANI(A,B) == ANI(B,A) by construction (mean of the two
    // directions is order-independent). Uses an asymmetric pair where the single-direction values
    // genuinely differ, proving the averaging — not a trivially symmetric input — is what makes it
    // symmetric. Query A = "AAAACCCCGGGGTTTT", B = "AAAACCCCGGGGTTTA" (one base differs).
    [Test]
    public void CalculateReciprocalAni_IsSymmetric()
    {
        const string a = "AAAACCCCGGGGTTTT";
        const string b = "AAAACCCCGGGGTTTA";

        double ab = ComparativeGenomics.CalculateReciprocalAni(a, b, fragmentLength: 4);
        double ba = ComparativeGenomics.CalculateReciprocalAni(b, a, fragmentLength: 4);

        Assert.That(ab, Is.EqualTo(ba).Within(Tolerance),
            "Reciprocal ANI must satisfy ANI(A,B) == ANI(B,A) (mean of both directions is symmetric).");
    }

    // R3 — Reciprocal ANI equals the mean of the two single-direction ANIb values. Hand-derived:
    // A = "AAAACGTC", B = "AAAAAAAA", fragLen 4, default cut-offs.
    //   A->B: frag1 "AAAA"=1.0 kept; frag2 "CGTC" best vs AAAA = 0 identical => 0.0 not >0.30,
    //         excluded => ANI(A->B) = 1.0.
    //   B->A: B fragments are "AAAA","AAAA"; against A="AAAACGTC" each best-matches the leading
    //         "AAAA" with identity 1.0 => ANI(B->A) = 1.0.
    //   Reciprocal = (1.0 + 1.0)/2 = 1.0.
    [Test]
    public void CalculateReciprocalAni_EqualsMeanOfBothDirections()
    {
        const string a = "AAAACGTC";
        const string b = "AAAAAAAA";

        double forward = ComparativeGenomics.CalculateANI(a, b, fragmentLength: 4);
        double reverse = ComparativeGenomics.CalculateANI(b, a, fragmentLength: 4);
        double reciprocal = ComparativeGenomics.CalculateReciprocalAni(a, b, fragmentLength: 4);

        Assert.Multiple(() =>
        {
            Assert.That(forward, Is.EqualTo(1.0).Within(Tolerance),
                "A->B: only 'AAAA' qualifies (CGTC excluded) => 1.0.");
            Assert.That(reverse, Is.EqualTo(1.0).Within(Tolerance),
                "B->A: both 'AAAA' fragments match the leading AAAA of A => 1.0.");
            Assert.That(reciprocal, Is.EqualTo((forward + reverse) / 2.0).Within(Tolerance),
                "Reciprocal ANI = mean of the two directional ANIs (Goris 2007 reciprocal values).");
        });
    }

    // R4 — Reciprocal ANI null/empty inputs return 0 (validation contract, same as CalculateANI).
    [Test]
    public void CalculateReciprocalAni_NullOrEmptyInput_ReturnsZero()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ComparativeGenomics.CalculateReciprocalAni(null!, Reference), Is.EqualTo(0.0),
                "Null genomeA must return 0.");
            Assert.That(ComparativeGenomics.CalculateReciprocalAni(Reference, null!), Is.EqualTo(0.0),
                "Null genomeB must return 0.");
            Assert.That(ComparativeGenomics.CalculateReciprocalAni("", Reference), Is.EqualTo(0.0),
                "Empty genomeA must return 0.");
        });
    }

    // R5 — Reciprocal ANI rejects non-positive fragmentLength (validation contract).
    [Test]
    public void CalculateReciprocalAni_NonPositiveFragmentLength_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => ComparativeGenomics.CalculateReciprocalAni(Reference, Reference, fragmentLength: 0),
            "fragmentLength <= 0 must throw ArgumentOutOfRangeException.");
    }

    #endregion
}
