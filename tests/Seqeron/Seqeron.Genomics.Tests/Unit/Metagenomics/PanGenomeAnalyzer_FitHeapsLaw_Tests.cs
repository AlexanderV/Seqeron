// PANGEN-HEAP-001 — Pan-Genome Growth Model (Heaps' Law)
// Evidence: docs/Evidence/PANGEN-HEAP-001-Evidence.md
// TestSpec: tests/TestSpecs/PANGEN-HEAP-001.md
// Source: Tettelin H et al. (2005) PNAS 102(39):13950-13955; Tettelin H et al. (2008)
//         Curr Opin Microbiol 11(5):472-477; micropan heaps() (R/powerlaw.R).

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Metagenomics;

namespace Seqeron.Genomics.Tests.Unit.Metagenomics;

[TestFixture]
public class PanGenomeAnalyzer_FitHeapsLaw_Tests
{
    // Exact tolerance for floating-point comparisons against analytic optima.
    private const double Tol = 1e-9;

    // Closed-form (Evidence §Test Datasets): new-gene curve x=[2,3], y=[8,4] lies on a
    // single power curve y = K * x^(-alpha): alpha = ln2/ln(3/2), K = 8 * 2^alpha.
    private static readonly double ExpectedAlphaClosed = Math.Log(2.0) / Math.Log(1.5);
    private static readonly double ExpectedKClosed = 8.0 * Math.Pow(2.0, ExpectedAlphaClosed);

    /// <summary>
    /// Builds presence rows in a fixed order from explicit cluster-id sets per genome.
    /// Cluster columns are derived deterministically; counts are order-invariant in column.
    /// </summary>
    private static List<PanGenomeAnalyzer.GenePresenceRow> Rows(
        params (string GenomeId, string[] Clusters)[] genomes)
    {
        var allClusters = genomes.SelectMany(g => g.Clusters).Distinct().ToList();
        var rows = new List<PanGenomeAnalyzer.GenePresenceRow>();
        foreach (var (genomeId, clusters) in genomes)
        {
            var set = new HashSet<string>(clusters);
            var presence = allClusters.ToDictionary(c => c, c => set.Contains(c));
            rows.Add(new PanGenomeAnalyzer.GenePresenceRow(
                GenomeId: genomeId,
                GenePresence: presence,
                TotalGenes: allClusters.Count,
                PresentGenes: set.Count));
        }
        return rows;
    }

    #region FitHeapsLaw (canonical)

    // M1 — Exact closed-form power-curve fit (Evidence derived dataset; micropan model).
    // Fixed order G1,G2,G3 with new-gene curve y=[8,4] -> alpha=ln2/ln1.5>1 -> closed.
    [Test]
    public void FitHeapsLaw_ClosedPowerCurve_RecoversExactParameters()
    {
        var rows = Rows(
            ("g1", new[] { "core" }),
            ("g2", new[] { "core", "n1", "n2", "n3", "n4", "n5", "n6", "n7", "n8" }),
            ("g3", new[] { "core", "m1", "m2", "m3", "m4" }));

        var fit = PanGenomeAnalyzer.FitHeapsLaw(rows, permutations: 1);

        Assert.Multiple(() =>
        {
            Assert.That(fit.Alpha, Is.EqualTo(ExpectedAlphaClosed).Within(Tol),
                "alpha must equal ln2/ln(3/2) for the exact power curve y=8*x^-alpha through (2,8),(3,4)");
            Assert.That(fit.Intercept, Is.EqualTo(ExpectedKClosed).Within(1e-6),
                "Intercept must equal 8*2^alpha (the analytic K for that curve)");
            Assert.That(fit.IsOpen, Is.False,
                "alpha > 1 means a closed pan-genome (micropan rule alpha>1 => closed)");
        });
    }

    // M2 — Constant new-gene curve -> best power fit alpha=0, K=mean=1 -> open.
    [Test]
    public void FitHeapsLaw_ConstantNewGeneCurve_ReturnsAlphaZeroOpen()
    {
        // G1{core}, G2 adds exactly 1 new, G3 adds exactly 1 new -> y=[1,1].
        var rows = Rows(
            ("g1", new[] { "core" }),
            ("g2", new[] { "core", "a" }),
            ("g3", new[] { "core", "b" }));

        var fit = PanGenomeAnalyzer.FitHeapsLaw(rows, permutations: 1);

        Assert.Multiple(() =>
        {
            Assert.That(fit.Alpha, Is.EqualTo(0.0).Within(Tol),
                "constant new-gene count is best fit by alpha=0 (no decay)");
            Assert.That(fit.Intercept, Is.EqualTo(1.0).Within(Tol),
                "K equals the constant new-gene count (mean = 1)");
            Assert.That(fit.IsOpen, Is.True,
                "alpha < 1 means an open pan-genome (micropan rule alpha<1 => open)");
        });
    }

    // M4 — New-gene count is FIRST APPEARANCE (micropan (cm==1)[i] & (cm==0)[i-1]), not a
    // recount of all-present clusters. Distinguishing matrix (fixed order g1,g2,g3):
    //   g1={core,a,b}, g2={core,a,b,c}, g3={core,a,b,c,d}
    // First-appearance curve: g2 introduces {c}=1, g3 introduces {d}=1 -> y=[1,1]
    // (a,b first appear at g1, never recounted; core never recounted) -> alpha=0, K=1, OPEN.
    // A WRONG impl that recounted non-core present clusters per genome would get
    // y=[3,4] (present-minus-nothing) -> a different fit. The exact (alpha=0,K=1,open)
    // result can only arise from the first-appearance rule on this matrix.
    [Test]
    public void FitHeapsLaw_CountsNewGenesByFirstAppearance()
    {
        var rows = Rows(
            ("g1", new[] { "core", "a", "b" }),
            ("g2", new[] { "core", "a", "b", "c" }),
            ("g3", new[] { "core", "a", "b", "c", "d" }));

        var fit = PanGenomeAnalyzer.FitHeapsLaw(rows, permutations: 1);

        Assert.Multiple(() =>
        {
            // First-appearance: exactly 1 new cluster at g2 (c) and 1 at g3 (d).
            Assert.That(fit.PredictNewGenes(2), Is.EqualTo(1.0).Within(1e-6),
                "g2 introduces exactly one new cluster (c); a,b,core already seen at g1");
            Assert.That(fit.PredictNewGenes(3), Is.EqualTo(1.0).Within(1e-6),
                "g3 introduces exactly one new cluster (d); a,b,c,core already seen");
            // Constant new-gene curve y=[1,1] -> alpha=0, K=1, open (Evidence derived dataset).
            Assert.That(fit.Alpha, Is.EqualTo(0.0).Within(Tol),
                "constant first-appearance curve y=[1,1] -> alpha=0");
            Assert.That(fit.Intercept, Is.EqualTo(1.0).Within(Tol),
                "constant first-appearance curve y=[1,1] -> K=1");
            Assert.That(fit.IsOpen, Is.True, "alpha < 1 -> open");
        });
    }

    // M4b — Distinguishing closed curve: a cluster present in g1, ABSENT in g2, then present
    // again in g3 must NOT be recounted as new at g3 (it first appeared at g1). Fixed order:
    //   g1={core,x}, g2={core,a,b}, g3={core,x,c}
    // First-appearance: g2 introduces {a,b}=2; g3 introduces {c}=1 (x already seen at g1) ->
    // y=[2,1] on the exact power curve K*N^-alpha with alpha=ln2/ln(3/2), K=2*2^alpha.
    // A WRONG impl comparing only to the immediately-previous genome would count x as new at
    // g3 (present g3, absent g2) -> y=[2,2] -> alpha=0/open, a different verdict.
    [Test]
    public void FitHeapsLaw_FirstAppearance_NotImmediatePredecessor()
    {
        var rows = Rows(
            ("g1", new[] { "core", "x" }),
            ("g2", new[] { "core", "a", "b" }),
            ("g3", new[] { "core", "x", "c" }));

        var fit = PanGenomeAnalyzer.FitHeapsLaw(rows, permutations: 1);

        // Exact curve y=[2,1]: alpha = ln2/ln(3/2), K = 2*2^alpha (Evidence model + hand calc).
        double expectedAlpha = Math.Log(2.0) / Math.Log(1.5);
        double expectedK = 2.0 * Math.Pow(2.0, expectedAlpha);
        Assert.Multiple(() =>
        {
            Assert.That(fit.PredictNewGenes(2), Is.EqualTo(2.0).Within(1e-6),
                "g2 introduces a,b = 2 new clusters");
            Assert.That(fit.PredictNewGenes(3), Is.EqualTo(1.0).Within(1e-6),
                "g3 introduces only c = 1 new cluster; x first appeared at g1, not recounted");
            Assert.That(fit.Alpha, Is.EqualTo(expectedAlpha).Within(Tol),
                "first-appearance curve y=[2,1] -> alpha=ln2/ln(3/2) (>1, closed)");
            Assert.That(fit.Intercept, Is.EqualTo(expectedK).Within(1e-6),
                "first-appearance curve y=[2,1] -> K=2*2^alpha");
            Assert.That(fit.IsOpen, Is.False, "alpha > 1 -> closed");
        });
    }

    // M6 — Fewer than two genomes -> degenerate fit, no exception.
    [Test]
    public void FitHeapsLaw_FewerThanTwoGenomes_ReturnsDegenerateFit()
    {
        var single = Rows(("g1", new[] { "a", "b", "c" }));
        var none = new List<PanGenomeAnalyzer.GenePresenceRow>();

        var fitSingle = PanGenomeAnalyzer.FitHeapsLaw(single, permutations: 5);
        var fitNone = PanGenomeAnalyzer.FitHeapsLaw(none, permutations: 5);

        Assert.Multiple(() =>
        {
            Assert.That(fitSingle.Intercept, Is.EqualTo(0.0),
                "one genome: new-gene curve N=2..G is empty -> degenerate Intercept 0");
            Assert.That(fitSingle.Alpha, Is.EqualTo(0.0), "one genome: degenerate Alpha 0");
            Assert.That(fitSingle.IsOpen, Is.False, "degenerate fit is reported not-open");
            Assert.That(fitSingle.PredictNewGenes(10), Is.EqualTo(0.0), "degenerate predictor returns 0");
            Assert.That(fitNone.Intercept, Is.EqualTo(0.0), "zero genomes: degenerate Intercept 0");
        });
    }

    // M7 — Null / empty matrix -> degenerate fit, no exception.
    [Test]
    public void FitHeapsLaw_NullOrEmptyMatrix_ReturnsDegenerateFit()
    {
        var fitNull = PanGenomeAnalyzer.FitHeapsLaw((IEnumerable<PanGenomeAnalyzer.GenePresenceRow>)null!);
        var fitEmpty = PanGenomeAnalyzer.FitHeapsLaw(Enumerable.Empty<PanGenomeAnalyzer.GenePresenceRow>());

        Assert.Multiple(() =>
        {
            Assert.That(fitNull.Intercept, Is.EqualTo(0.0), "null matrix -> degenerate Intercept 0");
            Assert.That(fitNull.PredictNewGenes(5), Is.EqualTo(0.0), "null matrix -> predictor 0");
            Assert.That(fitEmpty.Intercept, Is.EqualTo(0.0), "empty matrix -> degenerate Intercept 0");
        });
    }

    // M9 — Dictionary overload clusters then delegates to the matrix fit.
    [Test]
    public void FitHeapsLaw_DictionaryOverload_DelegatesToMatrixFit()
    {
        // Distinct sequences -> distinct singleton clusters; shared sequence 'CORE' clusters.
        var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>
        {
            ["g1"] = new List<(string, string)> { ("a", "CCCCCCCC") },
            ["g2"] = new List<(string, string)> { ("b", "CCCCCCCC"), ("c", "AAAAAAAA"), ("d", "GGGGGGGG") },
            ["g3"] = new List<(string, string)> { ("e", "CCCCCCCC"), ("f", "TTTTTTTT") }
        };

        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes).ToList();
        var matrix = PanGenomeAnalyzer.CreatePresenceAbsenceMatrix(genomes, clusters);

        var viaDict = PanGenomeAnalyzer.FitHeapsLaw(genomes, permutations: 8);
        var viaMatrix = PanGenomeAnalyzer.FitHeapsLaw(matrix, permutations: 8);

        Assert.Multiple(() =>
        {
            Assert.That(viaDict.Intercept, Is.EqualTo(viaMatrix.Intercept).Within(Tol),
                "dictionary overload must produce the same Intercept as the matrix path");
            Assert.That(viaDict.Alpha, Is.EqualTo(viaMatrix.Alpha).Within(Tol),
                "dictionary overload must produce the same Alpha as the matrix path");
            Assert.That(viaDict.IsOpen, Is.EqualTo(viaMatrix.IsOpen),
                "dictionary overload must produce the same open/closed verdict");
        });
    }

    // S1 — predictor n(N)=K*N^-alpha is non-increasing for alpha>0 (INV-06).
    [Test]
    public void FitHeapsLaw_PredictNewGenes_IsNonIncreasing()
    {
        var rows = Rows(
            ("g1", new[] { "core" }),
            ("g2", new[] { "core", "n1", "n2", "n3", "n4", "n5", "n6", "n7", "n8" }),
            ("g3", new[] { "core", "m1", "m2", "m3", "m4" }));

        var fit = PanGenomeAnalyzer.FitHeapsLaw(rows, permutations: 1);

        Assert.That(fit.PredictNewGenes(10), Is.LessThanOrEqualTo(fit.PredictNewGenes(5)),
            "for alpha>0 the new-gene predictor must not increase with genome index");
    }

    // S2 — Fitted parameters respect the micropan box bounds (INV-04/INV-05).
    [Test]
    public void FitHeapsLaw_RespectsParameterBounds()
    {
        var rows = Rows(
            ("g1", new[] { "core" }),
            ("g2", new[] { "core", "n1", "n2", "n3", "n4", "n5", "n6", "n7", "n8" }),
            ("g3", new[] { "core", "m1", "m2", "m3", "m4" }));

        var fit = PanGenomeAnalyzer.FitHeapsLaw(rows, permutations: 1);

        Assert.Multiple(() =>
        {
            Assert.That(fit.Alpha, Is.InRange(0.0, 2.0), "alpha must stay within micropan bounds [0,2]");
            Assert.That(fit.Intercept, Is.InRange(0.0, 10000.0),
                "Intercept must stay within micropan bounds [0,10000]");
        });
    }

    // C1 — Deterministic: same input -> identical fit across calls (fixed seed).
    [Test]
    public void FitHeapsLaw_IsDeterministic()
    {
        var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>
        {
            ["g1"] = new List<(string, string)> { ("a", "CCCCCCCC"), ("b", "AAAAAAAA") },
            ["g2"] = new List<(string, string)> { ("c", "CCCCCCCC"), ("d", "GGGGGGGG") },
            ["g3"] = new List<(string, string)> { ("e", "CCCCCCCC"), ("f", "TTTTTTTT") },
            ["g4"] = new List<(string, string)> { ("g", "CCCCCCCC"), ("h", "AATTAATT") }
        };

        var a = PanGenomeAnalyzer.FitHeapsLaw(genomes, permutations: 20);
        var b = PanGenomeAnalyzer.FitHeapsLaw(genomes, permutations: 20);

        Assert.Multiple(() =>
        {
            Assert.That(a.Intercept, Is.EqualTo(b.Intercept).Within(Tol), "fit must be deterministic (Intercept)");
            Assert.That(a.Alpha, Is.EqualTo(b.Alpha).Within(Tol), "fit must be deterministic (Alpha)");
        });
    }

    #endregion

    #region CreatePresenceAbsenceMatrix (canonical)

    // M8 — Presence/absence flags and counts are exact per genome (INV-02/INV-03).
    [Test]
    public void CreatePresenceAbsenceMatrix_PresenceAbsence_ExactFlags()
    {
        // g1 has both clusters; g2 has only the shared one.
        var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>
        {
            ["g1"] = new List<(string, string)> { ("shared", "ATGC"), ("only1", "GGGG") },
            ["g2"] = new List<(string, string)> { ("shared", "ATGC") }
        };
        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes).ToList();

        var matrix = PanGenomeAnalyzer.CreatePresenceAbsenceMatrix(genomes, clusters).ToList();

        var r1 = matrix.First(r => r.GenomeId == "g1");
        var r2 = matrix.First(r => r.GenomeId == "g2");
        Assert.Multiple(() =>
        {
            Assert.That(matrix.Count, Is.EqualTo(2), "one row per genome");
            Assert.That(r1.PresentGenes, Is.EqualTo(2), "g1 contains both clusters");
            Assert.That(r2.PresentGenes, Is.EqualTo(1), "g2 contains only the shared cluster");
            Assert.That(r1.GenePresence.Values.Count(v => v), Is.EqualTo(2),
                "g1 presence flags: both true");
            Assert.That(r2.GenePresence.Values.Count(v => v), Is.EqualTo(1),
                "g2 presence flags: exactly one true");
            Assert.That(r1.TotalGenes, Is.EqualTo(clusters.Count),
                "TotalGenes equals total cluster count");
        });
    }

    // M5 — Binarization (micropan pan.matrix[pan.matrix>0] <- 1): a cluster whose gene family
    // has MULTIPLE member genes in one genome is present ONCE, not once per member. g1 carries
    // two identical-sequence genes (a,b) -> they cluster into a single cluster (GeneIds=[a,b]);
    // that cluster must contribute PresentGenes = 1, not 2. A wrong impl counting gene
    // occurrences (or cluster members) instead of distinct present clusters would yield 2.
    [Test]
    public void CreatePresenceAbsenceMatrix_DuplicatePresence_CountsOnce()
    {
        var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>
        {
            // Two genes, identical sequence -> CD-HIT identity 1.0 -> ONE cluster with 2 members.
            ["g1"] = new List<(string, string)> { ("a", "ATGCATGC"), ("b", "ATGCATGC") }
        };
        var clusters = PanGenomeAnalyzer.ClusterGenes(genomes).ToList();

        var matrix = PanGenomeAnalyzer.CreatePresenceAbsenceMatrix(genomes, clusters).ToList();
        var row = matrix.Single();

        Assert.Multiple(() =>
        {
            Assert.That(clusters.Count, Is.EqualTo(1),
                "two identical-sequence genes collapse into a single ortholog cluster");
            Assert.That(row.GenePresence.Count, Is.EqualTo(1), "exactly one cluster column");
            Assert.That(row.PresentGenes, Is.EqualTo(1),
                "the multi-member cluster is present ONCE (binary), not once per member gene");
            Assert.That(row.GenePresence.Values.Count(v => v), Is.EqualTo(1),
                "exactly one present flag for the single cluster");
        });
    }

    // M5b — Binarization in the FIT path: duplicate copies of a cluster in a genome do not
    // inflate the new-gene count. Two genomes both carry a 2-member shared cluster plus
    // distinct singletons; the second genome contributes exactly its singleton as new (1),
    // never the shared cluster twice. Verified via the new-gene count at g2 = 1.
    [Test]
    public void FitHeapsLaw_BinarizesMultiMemberPresence()
    {
        var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>
        {
            // shared 2-member cluster ("ATGCATGC" twice) + one unique gene each.
            ["g1"] = new List<(string, string)> { ("s1", "ATGCATGC"), ("s2", "ATGCATGC"), ("u1", "AAAAAAAA") },
            ["g2"] = new List<(string, string)> { ("s3", "ATGCATGC"), ("s4", "ATGCATGC"), ("u2", "GGGGGGGG") }
        };

        var fit = PanGenomeAnalyzer.FitHeapsLaw(genomes, permutations: 1);

        // 2 genomes -> single new-gene point at N=2. g1 seeds {shared, u1}; g2 adds only u2
        // (its two shared-cluster members collapse to the already-seen shared cluster) -> new=1.
        Assert.That(fit.PredictNewGenes(2), Is.EqualTo(1.0).Within(1e-6),
            "g2 adds exactly one new cluster (u2); its duplicate shared-cluster members count once");
    }

    #endregion
}
