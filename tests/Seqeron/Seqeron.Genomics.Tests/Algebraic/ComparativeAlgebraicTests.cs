using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Algebraic;

using Gene = ComparativeGenomics.Gene;

/// <summary>
/// Algebraic-law tests for the Comparative area (ANI, genome comparison,
/// orthologs, reciprocal best hits, reversal distance).
///
/// Algebraic testing pins the metric axioms of the reversal (breakpoint) distance,
/// the symmetry/determinism of the comparison and ortholog relations, and the
/// self-identity of ANI.
/// — docs/checklists/06_ALGEBRAIC_TESTING.md §Description, rows 131, 133, 135, 136, 138.
/// </summary>
[TestFixture]
[Category("Algebraic")]
[Category("Comparative")]
public class ComparativeAlgebraicTests
{
    // ═══════════════════════════════════════════════════════════════════════
    // Unit: COMPGEN-REVERSAL-001 — Reversal (breakpoint) distance (Comparative)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 138.
    //
    // Model: the unsigned breakpoint reversal-distance lower bound ⌈b/2⌉ where b is
    //        the number of breakpoints between two gene-order permutations. The
    //        breakpoint distance is a metric on permutations, and ⌈b/2⌉ preserves
    //        the metric axioms.
    //   — docs/algorithms/Comparative_Genomics; ComparativeGenomics.CalculateReversalDistance.
    //
    // Laws under test (checklist row 138):
    //   • ID       — d(A, A) = 0.
    //   • COMM     — d(A, B) = d(B, A).
    //   • TRIANGLE — d(A, C) ≤ d(A, B) + d(B, C).
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>Three permutations of the same marker set {0..n-1}.</summary>
    private static Arbitrary<(int[] A, int[] B, int[] C)> PermutationTriple() =>
        (from n in Gen.Choose(1, 12)
         from ka in Gen.Choose(0, 1_000_000).ArrayOf(n)
         from kb in Gen.Choose(0, 1_000_000).ArrayOf(n)
         from kc in Gen.Choose(0, 1_000_000).ArrayOf(n)
         select (Perm(ka), Perm(kb), Perm(kc)))
        .ToArbitrary();

    // Order the markers 0..n-1 by independent random keys → a permutation.
    private static int[] Perm(int[] keys) =>
        Enumerable.Range(0, keys.Length).OrderBy(i => keys[i]).ThenBy(i => i).ToArray();

    [FsCheck.NUnit.Property]
    public Property Reversal_MetricAxioms()
    {
        return Prop.ForAll(PermutationTriple(), t =>
        {
            int aa = ComparativeGenomics.CalculateReversalDistance(t.A, t.A);
            int ab = ComparativeGenomics.CalculateReversalDistance(t.A, t.B);
            int ba = ComparativeGenomics.CalculateReversalDistance(t.B, t.A);
            int ac = ComparativeGenomics.CalculateReversalDistance(t.A, t.C);
            int bc = ComparativeGenomics.CalculateReversalDistance(t.B, t.C);
            bool id = aa == 0;
            bool comm = ab == ba;
            bool tri = ac <= ab + bc;
            return (id && comm && tri).Label($"id={id} comm={comm} tri={tri} (ac={ac},ab={ab},bc={bc})");
        });
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: COMPGEN-ANI-001 — Average Nucleotide Identity (Comparative)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 131.
    //
    // Model: ANIb (Goris 2007) is the mean identity of conserved genome fragments,
    //        a fraction in [0,1] (×100 for the percentage the checklist quotes). A
    //        genome compared to itself is perfectly conserved (ANI = 1.0).
    //   — docs/algorithms/Comparative_Genomics; ComparativeGenomics.CalculateANI.
    //
    // Laws (row 131): ID — ANI(A, A) = 1.0 (perfect self-identity).
    //                 COMM — ANI symmetric on shared content.
    // ═══════════════════════════════════════════════════════════════════════

    private static Arbitrary<string> DnaArbitrary(int minLen) =>
        Gen.Elements('A', 'C', 'G', 'T').ArrayOf().Where(a => a.Length >= minLen)
            .Select(a => new string(a)).ToArbitrary();

    /// <summary>ID: a genome compared to itself has ANI = 1.0 (max self-identity).</summary>
    [FsCheck.NUnit.Property]
    public Property Ani_Identity_SelfIsOne()
    {
        return Prop.ForAll(DnaArbitrary(minLen: 40), seq =>
        {
            double ani = ComparativeGenomics.CalculateANI(seq, seq, fragmentLength: 20);
            return (Math.Abs(ani - 1.0) < 1e-9).Label($"ANI(self)={ani} for len={seq.Length}");
        });
    }

    /// <summary>
    /// COMM: ANI is symmetric on shared content — two genomes built from the same
    /// 20-mer blocks in swapped order are each other's perfect match, so
    /// ANI(A,B) = ANI(B,A).
    /// </summary>
    [Test]
    public void Ani_Commutative_Symmetric()
    {
        string x = "ACGTACGTACGTACGTACGT";
        string y = "TTTTGGGGCCCCAAAATTTT";
        string a = x + y;
        string b = y + x;
        double ab = ComparativeGenomics.CalculateANI(a, b, fragmentLength: 20);
        double ba = ComparativeGenomics.CalculateANI(b, a, fragmentLength: 20);
        ab.Should().BeApproximately(ba, 1e-9);
        ab.Should().BeApproximately(1.0, 1e-9);
    }

    // ── Gene / genome fixtures for RBH / ortholog / compare laws ──────────────

    private const string SeqA = "ATGCGTACGTTAGCCGATCGATCGTAGCTA";
    private const string SeqB = "TTTTGGGGCCCCAAAATTTTGGGGCCCCAA";

    private static List<Gene> Genome(string genomeId) => new()
    {
        new Gene($"{genomeId}_A", genomeId, 0, 30, '+', SeqA),
        new Gene($"{genomeId}_B", genomeId, 40, 70, '+', SeqB),
    };

    private static HashSet<(string, string)> UnorderedPairs(IEnumerable<ComparativeGenomics.OrthologPair> pairs) =>
        pairs.Select(p => string.CompareOrdinal(p.Gene1Id, p.Gene2Id) <= 0
                ? (p.Gene1Id, p.Gene2Id) : (p.Gene2Id, p.Gene1Id)).ToHashSet();

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: COMPGEN-RBH-001 — Reciprocal best hits (Comparative), row 136.
    // COMM — RBH(A,B) and RBH(B,A) yield the same unordered ortholog pairs.
    // IDEMP — deterministic (re-running gives the identical result).
    //   — ComparativeGenomics.FindReciprocalBestHits.
    // ═══════════════════════════════════════════════════════════════════════

    [Test]
    public void Rbh_Commutative_Symmetric()
    {
        var ab = ComparativeGenomics.FindReciprocalBestHits(Genome("g1"), Genome("g2")).ToList();
        var ba = ComparativeGenomics.FindReciprocalBestHits(Genome("g2"), Genome("g1")).ToList();
        ab.Should().NotBeEmpty();
        UnorderedPairs(ab).Should().BeEquivalentTo(UnorderedPairs(ba));
    }

    [Test]
    public void Rbh_Idempotent_Deterministic()
    {
        var first = ComparativeGenomics.FindReciprocalBestHits(Genome("g1"), Genome("g2")).ToList();
        var second = ComparativeGenomics.FindReciprocalBestHits(Genome("g1"), Genome("g2")).ToList();
        second.Should().BeEquivalentTo(first);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: COMPGEN-ORTHO-001 — Orthologs (Comparative), row 135.
    // COMM — the ortholog relation is symmetric.  IDEMP — deterministic.
    //   — ComparativeGenomics.FindOrthologs (RBH-based).
    // ═══════════════════════════════════════════════════════════════════════

    [Test]
    public void Orthologs_Commutative_Symmetric()
    {
        var ab = ComparativeGenomics.FindOrthologs(Genome("g1"), Genome("g2")).ToList();
        var ba = ComparativeGenomics.FindOrthologs(Genome("g2"), Genome("g1")).ToList();
        ab.Should().NotBeEmpty();
        UnorderedPairs(ab).Should().BeEquivalentTo(UnorderedPairs(ba));
    }

    [Test]
    public void Orthologs_Idempotent_Deterministic()
    {
        var first = ComparativeGenomics.FindOrthologs(Genome("g1"), Genome("g2")).ToList();
        var second = ComparativeGenomics.FindOrthologs(Genome("g1"), Genome("g2")).ToList();
        second.Should().BeEquivalentTo(first);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: COMPGEN-COMPARE-001 — Genome comparison (Comparative), row 133.
    // COMM — conserved-gene count is symmetric; genome-specific counts swap.
    // IDEMP — deterministic.
    //   — ComparativeGenomics.CompareGenomes.
    // ═══════════════════════════════════════════════════════════════════════

    [Test]
    public void Compare_Commutative_Symmetric()
    {
        var ab = ComparativeGenomics.CompareGenomes(Genome("g1"), Genome("g2"));
        var ba = ComparativeGenomics.CompareGenomes(Genome("g2"), Genome("g1"));
        ab.ConservedGenes.Should().Be(ba.ConservedGenes);
        ab.GenomeSpecificGenes1.Should().Be(ba.GenomeSpecificGenes2);
        ab.GenomeSpecificGenes2.Should().Be(ba.GenomeSpecificGenes1);
    }

    [Test]
    public void Compare_Idempotent_Deterministic()
    {
        var first = ComparativeGenomics.CompareGenomes(Genome("g1"), Genome("g2"));
        var second = ComparativeGenomics.CompareGenomes(Genome("g1"), Genome("g2"));
        second.ConservedGenes.Should().Be(first.ConservedGenes);
        second.OverallSynteny.Should().Be(first.OverallSynteny);
    }
}
