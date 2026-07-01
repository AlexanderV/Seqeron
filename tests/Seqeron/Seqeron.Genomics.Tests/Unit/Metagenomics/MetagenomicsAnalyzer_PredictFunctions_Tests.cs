// META-FUNC-001 — Functional Prediction (homology transfer + pathway enrichment)
// Evidence: docs/Evidence/META-FUNC-001-Evidence.md
// TestSpec: tests/TestSpecs/META-FUNC-001.md
// Source: Altschul et al. (1990), J Mol Biol 215(3):403-410, NCBI BLAST tutorial;
//         NCBI blast_stat.c (ungapped BLOSUM62 lambda/K); NCBI BLOSUM62 matrix;
//         PNNL Proteomics Data Analysis tutorial §8.2 (hypergeometric ORA).

using NUnit.Framework;
using Seqeron.Genomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests.Unit.Metagenomics;

[TestFixture]
public class MetagenomicsAnalyzer_PredictFunctions_Tests
{
    // Expected bit-score/E-value LITERALS below are taken from the external sources
    // (NCBI BLAST tutorial formulas + ungapped BLOSUM62 lambda=0.3176, K=0.134, W diagonal=11),
    // NOT recomputed from the implementation's own constants — see docs/Evidence/META-FUNC-001-Evidence.md.

    private static IReadOnlyDictionary<string, (string, string, string)> Db(
        params (string Sig, string Func, string Path, string Ko)[] entries)
        => entries.ToDictionary(e => e.Sig, e => (e.Func, e.Path, e.Ko));

    #region PredictFunctions — bit score / E-value (M1, M2)

    // M1 — Bit score of an exact BLOSUM62 self-match. Signature "WWW" → S = 3*11 = 33,
    // S' = (lambda*S - ln K)/ln 2 = 18.0202932787533. Evidence: Altschul tutorial + NCBI lambda/K + W diagonal=11.
    [Test]
    public void PredictFunctions_ExactMatch_ComputesBlastBitScore()
    {
        var proteins = new[] { ("g1", "WWW") };
        var db = Db(("WWW", "tryptophanase", "Amino acid metabolism", "K01667"));

        var result = MetagenomicsAnalyzer.PredictFunctions(proteins, db).Single();

        // Sourced literal (Evidence §"BLAST bit-score / E-value self-consistency"): S=33 (3*W=3*11),
        // S' = (0.3176*33 - ln 0.134)/ln 2 = 18.0202932787533 (NOT recomputed from the code's constants).
        Assert.That(result.BitScore, Is.EqualTo(18.0202932787533).Within(1e-12),
            "Bit score must equal the sourced 18.0202932787533 for raw score S=33 (WWW, BLOSUM62 W=11).");
    }

    // M2 — E-value of the same hit, m=n=3: E = K*m*n*e^(-lambda*S) = 3.3852730346546e-5.
    // Evidence: Altschul tutorial E = K m n e^(-lambda S) = m n 2^(-S').
    [Test]
    public void PredictFunctions_ExactMatch_ComputesBlastEValue()
    {
        var proteins = new[] { ("g1", "WWW") };
        var db = Db(("WWW", "tryptophanase", "Amino acid metabolism", "K01667"));

        var result = MetagenomicsAnalyzer.PredictFunctions(proteins, db).Single();

        // Sourced literal (Evidence dataset): E = K*m*n*e^(-lambda*S) = m*n*2^(-S') = 3.3852730346546e-5
        // for S=33, m=n=3 (the two cited E-value forms agree to machine precision).
        // Asserted against the published number, not a recomputation of the code's own formula.
        double bitForm = 3.0 * 3.0 * Math.Pow(2.0, -18.0202932787533); // E from the sourced bit score S'
        Assert.Multiple(() =>
        {
            Assert.That(result.EValue, Is.EqualTo(3.3852730346546e-5).Within(1e-18),
                "E-value must equal the sourced 3.3852730346546e-5 (K*m*n*e^(-lambda*S), S=33, m=n=3).");
            // INV-02: the two cited forms agree; bitForm uses the 13-sig-fig bit-score literal, so
            // they coincide only to the rounding of that literal (~1.5e-18 here) — assert that.
            Assert.That(result.EValue, Is.EqualTo(bitForm).Within(1e-17),
                "K*m*n*e^(-lambda*S) must equal m*n*2^(-S') (INV-02), the second cited E-value form.");
        });
    }

    #endregion

    #region PredictFunctions — best hit, fields, matching (M3, M6, S1, S2, S3, S4)

    // M3 — Best-hit selection: gene matches two signatures; transfer the one with the lower E-value.
    // "WW" → S=22 (higher score → lower E) beats "AAAA" → S=16. Evidence: BLAST best-hit (E ranks significance).
    [Test]
    public void PredictFunctions_MultipleMatches_TransfersLowestEValueHit()
    {
        var proteins = new[] { ("g1", "AAAAWW") }; // contains both "AAAA" (S=16) and "WW" (S=22)
        var db = Db(
            ("AAAA", "weakHit", "PathA", "K00001"),
            ("WW", "strongHit", "PathB", "K00002"));

        var result = MetagenomicsAnalyzer.PredictFunctions(proteins, db).Single();

        Assert.Multiple(() =>
        {
            Assert.That(result.Function, Is.EqualTo("strongHit"),
                "Best hit is the higher-scoring 'WW' (S=22 > 16), which has the lower E-value.");
            // WW: m=6 (query len), n=2; S = 2*11 = 22. Sourced E = 0.134*6*2*e^(-0.3176*22)
            // = 0.0014851955539388528 (hand-computed from the cited formula + lambda/K).
            Assert.That(result.EValue, Is.EqualTo(0.0014851955539388528).Within(1e-15),
                "Reported E-value must be the sourced E of the best 'WW' hit (m=6, n=2, S=22).");
        });
    }

    // M6 — Annotation fields are transferred from the matched database entry.
    [Test]
    public void PredictFunctions_Match_TransfersFunctionPathwayAndKo()
    {
        var proteins = new[] { ("geneX", "MMM") };
        var db = Db(("MMM", "methyltransferase", "Methane metabolism", "K00399"));

        var result = MetagenomicsAnalyzer.PredictFunctions(proteins, db).Single();

        Assert.Multiple(() =>
        {
            Assert.That(result.GeneId, Is.EqualTo("geneX"), "GeneId carried through.");
            Assert.That(result.Function, Is.EqualTo("methyltransferase"), "Function transferred from DB.");
            Assert.That(result.Pathway, Is.EqualTo("Methane metabolism"), "Pathway transferred from DB.");
            Assert.That(result.KoNumber, Is.EqualTo("K00399"), "KO number transferred from DB.");
        });
    }

    // S1 — Empty/whitespace protein sequence yields no annotation.
    [Test]
    public void PredictFunctions_EmptySequence_YieldsNoAnnotation()
    {
        var proteins = new[] { ("g1", ""), ("g2", "   ") };
        var db = Db(("WWW", "f", "p", "k"));

        var result = MetagenomicsAnalyzer.PredictFunctions(proteins, db).ToList();

        Assert.That(result, Is.Empty, "Empty or whitespace sequences cannot match and yield no annotation.");
    }

    // S2 — A gene that matches no signature yields no annotation.
    [Test]
    public void PredictFunctions_NoMatchingSignature_YieldsNoAnnotation()
    {
        var proteins = new[] { ("g1", "ACDEFG") };
        var db = Db(("WWWWW", "f", "p", "k"));

        var result = MetagenomicsAnalyzer.PredictFunctions(proteins, db).ToList();

        Assert.That(result, Is.Empty, "Homology transfer requires a hit; no signature occurrence ⇒ no annotation.");
    }

    // S3 — Empty database yields no annotations (and does not throw).
    [Test]
    public void PredictFunctions_EmptyDatabase_YieldsNoAnnotation()
    {
        var proteins = new[] { ("g1", "WWW") };
        var db = Db();

        var result = MetagenomicsAnalyzer.PredictFunctions(proteins, db).ToList();

        Assert.That(result, Is.Empty, "No signatures to match ⇒ no annotations.");
    }

    // S4 — Null arguments throw ArgumentNullException.
    [Test]
    public void PredictFunctions_NullArguments_Throws()
    {
        var db = Db(("WWW", "f", "p", "k"));
        Assert.Multiple(() =>
        {
            Assert.That(() => MetagenomicsAnalyzer.PredictFunctions(null!, db).ToList(),
                NUnit.Framework.Throws.ArgumentNullException, "Null proteins must throw.");
            Assert.That(() => MetagenomicsAnalyzer.PredictFunctions(new[] { ("g", "WWW") }, null!).ToList(),
                NUnit.Framework.Throws.ArgumentNullException, "Null database must throw.");
        });
    }

    #endregion

    #region FindPathwayEnrichment — hypergeometric (M4, M5, S5, S6)

    // M4 — ORA worked example: N=8000, M=400, n=100, x=20 → P(X>=20) = 7.88e-8.
    // Evidence: PNNL ORA §8.2 worked example.
    [Test]
    public void FindPathwayEnrichment_WorkedExample_MatchesPublishedPValue()
    {
        // Background: 8000 genes b0001..b8000. Pathway "S": first 400. Query: 100 genes,
        // 20 of which are in S (b0001..b0020) and 80 outside S (b0401..b0480).
        var background = Enumerable.Range(1, 8000).Select(i => $"b{i:D4}").ToArray();
        var pathwayMembers = Enumerable.Range(1, 400).Select(i => $"b{i:D4}").ToArray();
        var query = Enumerable.Range(1, 20).Select(i => $"b{i:D4}")
            .Concat(Enumerable.Range(401, 80).Select(i => $"b{i:D4}")).ToArray();

        var pathwayDb = new Dictionary<string, IReadOnlyCollection<string>>
        {
            ["S"] = pathwayMembers,
        };

        var result = MetagenomicsAnalyzer.FindPathwayEnrichment(query, pathwayDb, background).Single();

        Assert.Multiple(() =>
        {
            Assert.That(result.Overlap, Is.EqualTo(20), "Overlap x must be 20.");
            Assert.That(result.PathwaySize, Is.EqualTo(400), "Pathway size M must be 400.");
            Assert.That(result.QuerySize, Is.EqualTo(100), "Query size n must be 100.");
            Assert.That(result.BackgroundSize, Is.EqualTo(8000), "Background size N must be 8000.");
            Assert.That(result.PValue, Is.EqualTo(7.88e-8).Within(1e-9),
                "Hypergeometric P(X>=20) must equal the published 7.88e-8.");
        });
    }

    // M5 / INV-05 — No overlap (x=0) ⇒ p-value = 1.
    [Test]
    public void FindPathwayEnrichment_NoOverlap_ReturnsPValueOne()
    {
        var background = new[] { "a", "b", "c", "d" };
        var query = new[] { "a", "b" };
        var pathwayDb = new Dictionary<string, IReadOnlyCollection<string>>
        {
            ["P"] = new[] { "c", "d" }, // disjoint from query
        };

        var result = MetagenomicsAnalyzer.FindPathwayEnrichment(query, pathwayDb, background).Single();

        Assert.That(result.PValue, Is.EqualTo(1.0).Within(1e-12),
            "x=0 ⇒ empty sum ⇒ P(X>=0)=1 (no over-representation).");
    }

    // S5 / INV-05 — Empty pathway (M=0) ⇒ p-value = 1.
    [Test]
    public void FindPathwayEnrichment_EmptyPathway_ReturnsPValueOne()
    {
        var background = new[] { "a", "b", "c" };
        var query = new[] { "a", "b" };
        var pathwayDb = new Dictionary<string, IReadOnlyCollection<string>>
        {
            ["P"] = Array.Empty<string>(),
        };

        var result = MetagenomicsAnalyzer.FindPathwayEnrichment(query, pathwayDb, background).Single();

        Assert.That(result.PValue, Is.EqualTo(1.0).Within(1e-12),
            "M=0 (degenerate margin) ⇒ p-value = 1.");
    }

    // S6 / INV-05 — Empty query (n=0) ⇒ p-value = 1.
    [Test]
    public void FindPathwayEnrichment_EmptyQuery_ReturnsPValueOne()
    {
        var background = new[] { "a", "b", "c" };
        var query = Array.Empty<string>();
        var pathwayDb = new Dictionary<string, IReadOnlyCollection<string>>
        {
            ["P"] = new[] { "a", "b" },
        };

        var result = MetagenomicsAnalyzer.FindPathwayEnrichment(query, pathwayDb, background).Single();

        Assert.That(result.PValue, Is.EqualTo(1.0).Within(1e-12),
            "n=0 (empty query) ⇒ p-value = 1.");
    }

    // S4 (enrichment) — Null arguments throw.
    [Test]
    public void FindPathwayEnrichment_NullArguments_Throws()
    {
        var pathwayDb = new Dictionary<string, IReadOnlyCollection<string>> { ["P"] = new[] { "a" } };
        Assert.Multiple(() =>
        {
            Assert.That(() => MetagenomicsAnalyzer.FindPathwayEnrichment(null!, pathwayDb),
                NUnit.Framework.Throws.ArgumentNullException, "Null query must throw.");
            Assert.That(() => MetagenomicsAnalyzer.FindPathwayEnrichment(new[] { "a" }, null!),
                NUnit.Framework.Throws.ArgumentNullException, "Null pathway database must throw.");
        });
    }

    #endregion

    #region Property-based invariants (C1, C2)

    // C1 / INV-01, INV-03 — Bit score strictly increases and E-value strictly decreases with raw score.
    [Test]
    public void FunctionalBitScore_And_EValue_AreMonotonicInRawScore()
    {
        Assert.Multiple(() =>
        {
            for (int s = 1; s <= 50; s++)
            {
                Assert.That(MetagenomicsAnalyzer.FunctionalBitScore(s + 1),
                    Is.GreaterThan(MetagenomicsAnalyzer.FunctionalBitScore(s)),
                    $"Bit score must strictly increase with raw score (at S={s}).");
                Assert.That(MetagenomicsAnalyzer.ExpectedValue(s + 1, 100, 100),
                    Is.LessThan(MetagenomicsAnalyzer.ExpectedValue(s, 100, 100)),
                    $"E-value must strictly decrease with raw score (at S={s}).");
            }
        });
    }

    // Blosum62SelfScore — direct check against NCBI BLOSUM62 diagonal values
    // (W=11, A=4, C=9, M=5, R=5, N=6 from https://ftp.ncbi.nlm.nih.gov/blast/matrices/BLOSUM62).
    // "WACMRN" self-score = 11+4+9+5+5+6 = 40. Also: empty ⇒ 0; unknown residue contributes 0.
    [Test]
    public void Blosum62SelfScore_SumsDiagonalScores()
    {
        Assert.Multiple(() =>
        {
            Assert.That(MetagenomicsAnalyzer.Blosum62SelfScore("W"), Is.EqualTo(11), "W diagonal = 11.");
            Assert.That(MetagenomicsAnalyzer.Blosum62SelfScore("WACMRN"), Is.EqualTo(40),
                "Self-score is the sum of BLOSUM62 diagonals: 11+4+9+5+5+6 = 40.");
            Assert.That(MetagenomicsAnalyzer.Blosum62SelfScore(""), Is.EqualTo(0), "Empty segment ⇒ 0.");
            Assert.That(MetagenomicsAnalyzer.Blosum62SelfScore("XZ"), Is.EqualTo(0),
                "Unknown residues (not in the 20-AA diagonal) contribute 0.");
        });
    }

    // HypergeometricUpperTail — small EXACT hand-computed case independent of the published example.
    // N=10, M=4, n=5, x=3 ⇒ P(X>=3) = [C(4,3)C(6,2) + C(4,4)C(6,1)] / C(10,5) = (4*15 + 1*6)/252
    //                                = 66/252 = 11/42 = 0.2619047619047619 (exact rational).
    [Test]
    public void HypergeometricUpperTail_SmallExactCase_MatchesHandComputation()
    {
        double p = MetagenomicsAnalyzer.HypergeometricUpperTail(x: 3, bigN: 10, bigM: 4, n: 5);
        Assert.That(p, Is.EqualTo(11.0 / 42.0).Within(1e-12),
            "P(X>=3) for N=10,M=4,n=5 must equal the exact hand-computed 11/42.");
    }

    // FindPathwayEnrichment — null background falls back to the union of all pathway members
    // (then the query is always added to the universe). Documented in §3.1 / §4.B.
    [Test]
    public void FindPathwayEnrichment_NullBackground_UsesUnionOfPathwayMembers()
    {
        var query = new[] { "a", "b" };
        var pathwayDb = new Dictionary<string, IReadOnlyCollection<string>>
        {
            ["P"] = new[] { "a", "b", "c", "d" }, // 4 members; union background = {a,b,c,d}, +query = {a,b,c,d}
        };

        // No background arg ⇒ background = union of members = {a,b,c,d}; N=4, M=4, n=2, x=2.
        // P(X>=2) with M=N: every draw is in the pathway ⇒ p = 1.
        var result = MetagenomicsAnalyzer.FindPathwayEnrichment(query, pathwayDb).Single();

        Assert.Multiple(() =>
        {
            Assert.That(result.BackgroundSize, Is.EqualTo(4), "Background = union of pathway members {a,b,c,d}.");
            Assert.That(result.PathwaySize, Is.EqualTo(4), "M = 4 (all members in background).");
            Assert.That(result.Overlap, Is.EqualTo(2), "Both query genes are in the pathway.");
            Assert.That(result.PValue, Is.EqualTo(1.0).Within(1e-12),
                "With M=N every sampled gene is a success ⇒ P(X>=x)=1.");
        });
    }

    // C2 / INV-04 — Hypergeometric p-value is always within [0,1] over a sweep of valid inputs.
    [Test]
    public void HypergeometricUpperTail_AlwaysWithinUnitInterval()
    {
        Assert.Multiple(() =>
        {
            foreach (int bigN in new[] { 10, 100, 1000 })
                foreach (int bigM in new[] { 1, bigN / 2, bigN })
                    foreach (int n in new[] { 1, bigN / 3, bigN })
                        foreach (int x in new[] { 0, 1, Math.Min(bigM, n) })
                        {
                            double p = MetagenomicsAnalyzer.HypergeometricUpperTail(x, bigN, bigM, n);
                            Assert.That(p, Is.InRange(0.0, 1.0),
                                $"p-value must be in [0,1] for N={bigN}, M={bigM}, n={n}, x={x}.");
                        }
        });
    }

    #endregion
}
