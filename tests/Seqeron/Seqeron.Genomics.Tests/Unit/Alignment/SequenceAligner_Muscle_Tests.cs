using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A1-10c part a, B09 F65): <see cref="SequenceAligner.MultipleAlignMuscle(IEnumerable{string}, MuscleOptions?, MatrixScoring?)"/>, the port of MUSCLE 3.8.1551's
/// progressive alignment (Edgar 2004, NAR 32:1792) for -maxiters 1 and 2. The installed muscle is 5.1 (a different program), so the row-for-row oracle is the MUSCLE 3.8.1551 binary of the
/// Docker hand-off package docs/Validation/review-2026-09/handoff/A1-10c-muscle (420 sets) - AWAITING ITS RETURN. What is locked here is everything that can be derived from the published
/// definitions and the program's own source without running it: the log-expectation score formula and its parameters (VTML200 table, center -0.52), the SPN nucleotide score, the
/// ClustalW weights worked example of clwwt.cpp, the UPGMA test matrix of upgma2.cpp, the Kimura / Dayhoff-table distance (the 75 %/250 PAM comment of msadistkimura.cpp), the k-mer
/// distance definition, and the optimality of the profile-profile dynamic program against an exhaustive enumeration of all alignments under the same position-specific gap model.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_Muscle_Tests
{
    private static readonly string[] Family =
    [
        "MKVLAAGIVGLLLAQPWTTQAEAKP", "MKVLAGIVGLLAQPWTQAEAKPLL", "MRVLAAGIIGLLLQPWTTQAEAK", "KVLAAGIVALLLAQPWTTAEAKPS", "MKVAAGIVGLLLAQPWTTQEAKP",
    ];

    private static char[] MuscleRow(string s) => s.ToCharArray();

    private static SequenceAligner.MuscleEngine Engine(MuscleOptions? o = null, params string[] seqs) =>
        new(o ?? new MuscleOptions { SequenceType = MuscleSequenceType.Protein }, seqs.Length == 0 ? ["ACDEFGHIK", "ACDEFGHIK"] : seqs);

    // ---- scores: the definitions ---------------------------------------------------------------------------------------------------------------------

    /// <summary>LE score (Edgar 2004): for two single-residue columns the score is ln(VTML200 likelihood ratio) - center, center = -0.52 (params.cpp SetDefaultsLE), occupancies 1.</summary>
    [Test]
    public void LogExpectationScore_OfSingleResidueColumns_IsLnOfTheVtml200Entry_MinusCenter()
    {
        var e = Engine();
        var a = e.ProfileFromRows([['A']], [1f])[0];
        var c = e.ProfileFromRows([['C']], [1f])[0];
        Assert.Multiple(() =>
        {
            Assert.That(e.ScoreProfPos2(a, a), Is.EqualTo(MathF.Log(2.25080f) + 0.52f).Within(1e-6));   // VTML_LA[A][A] = 2.25080
            Assert.That(e.ScoreProfPos2(c, c), Is.EqualTo(MathF.Log(15.79469f) + 0.52f).Within(1e-6));  // VTML_LA[C][C] = 15.79469
            Assert.That(e.ScoreProfPos2(a, c), Is.EqualTo(MathF.Log(1.31180f) + 0.52f).Within(1e-6));   // VTML_LA[A][C] = 1.31180 (symmetric table)
            Assert.That(e.ScoreProfPos2(c, a), Is.EqualTo(e.ScoreProfPos2(a, c)));
        });
    }

    /// <summary>The LE score is weighted by the column occupancies (a column that is half gaps counts half) and a profile of an all-wildcard X column scores the mean.</summary>
    [Test]
    public void LogExpectationScore_IsWeightedByOccupancy()
    {
        var e = Engine();
        // two sequences, the second column has one gap: occupancy 0.5 (weights 0.5 each)
        var half = e.ProfileFromRows([['A', 'A'], ['A', '-']], [0.5f, 0.5f]);
        var full = e.ProfileFromRows([['A', 'A']], [1f]);
        float s1 = e.ScoreProfPos2(half[1], full[0]);
        float s0 = e.ScoreProfPos2(half[0], full[0]);
        Assert.That(half[1].Occ, Is.EqualTo(0.5f));
        Assert.That(s1, Is.EqualTo(s0 * 0.5f).Within(1e-5));
    }

    /// <summary>SPN (nucleotides): the BLASTZ matrix with the center 2 x 30 added (params.cpp: center pre-added into the matrix), score = sum f_i g_j m_ij - 0.</summary>
    [Test]
    public void NucleotideScore_IsTheBlastzMatrixPlusCenter()
    {
        var e = Engine(new MuscleOptions { SequenceType = MuscleSequenceType.Dna }, "ACGTACGT", "ACGTACGT");
        var a = e.ProfileFromRows([['A']], [1f])[0];
        var c = e.ProfileFromRows([['C']], [1f])[0];
        var g = e.ProfileFromRows([['G']], [1f])[0];
        var t = e.ProfileFromRows([['T']], [1f])[0];
        Assert.Multiple(() =>
        {
            Assert.That(e.ScoreProfPos2(a, a), Is.EqualTo(91f + 60f));   // A,A = 91 (BLASTZ) + 60
            Assert.That(e.ScoreProfPos2(a, c), Is.EqualTo(-114f + 60f));
            Assert.That(e.ScoreProfPos2(a, g), Is.EqualTo(-31f + 60f));   // transition
            Assert.That(e.ScoreProfPos2(a, t), Is.EqualTo(-123f + 60f));
            Assert.That(e.ScoreProfPos2(c, c), Is.EqualTo(100f + 60f));
            Assert.That(e.ScoreProfPos2(g, g), Is.EqualTo(100f + 60f));
        });
    }

    [Test]
    public void WildcardCounts_FollowMsa2Cpp_IncludingTheNucleotideAxRQuirk()
    {
        var p = Engine();
        var prot = p.ProfileFromRows([['B', 'Z', 'X', 'A']], [1f]);
        var d = Engine(new MuscleOptions { SequenceType = MuscleSequenceType.Dna }, "ACGTACGT", "ACGTACGT");
        var dna = d.ProfileFromRows([['X', 'N', 'R', 'A']], [1f]);
        Assert.Multiple(() =>
        {
            Assert.That(prot[0].Counts[2], Is.EqualTo(0.5f), "B = D or N (D = 2, N = 11 in the A C D E F G H I K L M N P Q R S T V W Y order)");
            Assert.That(prot[0].Counts[11], Is.EqualTo(0.5f));
            Assert.That(prot[1].Counts[3], Is.EqualTo(0.5f), "Z = E or Q");
            Assert.That(prot[1].Counts[13], Is.EqualTo(0.5f));
            Assert.That(prot[2].Counts.Sum(), Is.EqualTo(1f).Within(1e-6), "X = 1/20 each");
            Assert.That(prot[2].Counts[0], Is.EqualTo(0.05f));
            // nucleotides: MSA::GetFractionalWeightedCounts compares the NX letter with AX_R (= 14 = NX_X) and AX_Y (= 19): X counts as G or A, N and R as w/20 each of the four letters
            Assert.That(dna[0].Counts[2], Is.EqualTo(0.5f));
            Assert.That(dna[0].Counts[0], Is.EqualTo(0.5f));
            Assert.That(dna[1].Counts[..4], Is.EqualTo(new[] { 0.05f, 0.05f, 0.05f, 0.05f }));
            Assert.That(dna[2].Counts[..4], Is.EqualTo(new[] { 0.05f, 0.05f, 0.05f, 0.05f }));
            Assert.That(dna[3].Counts[0], Is.EqualTo(1f));
        });
    }

    // ---- distances and trees ------------------------------------------------------------------------------------------------------------------------

    /// <summary>msadistkimura.cpp: -ln(1 - p - p^2/5) below 75 % divergence; the ClustalW (dayhoff.h) table above ("250 PAMs = 80.3% observed d"); 10 above 93 %.</summary>
    [Test]
    public void KimuraDistance_FollowsTheFormulaAndTheDayhoffTable()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceAligner.MuscleEngine.KimuraDist(0.8), Is.EqualTo(-Math.Log(1 - 0.2 - 0.2 * 0.2 / 5)).Within(1e-12));
            Assert.That(SequenceAligner.MuscleEngine.KimuraDist(1.0), Is.EqualTo(0.0).Within(1e-12));
            Assert.That(SequenceAligner.MuscleEngine.KimuraDist(0.25), Is.EqualTo(1.95).Within(1e-12), "p = 0.75: table entry 0 = 195 PAMs");
            Assert.That(SequenceAligner.MuscleEngine.KimuraDist(1 - 0.803), Is.EqualTo(2.50).Within(1e-12), "source comment: 250 PAMs = 80.3 % observed d");
            Assert.That(SequenceAligner.MuscleEngine.KimuraDist(1 - 0.930), Is.EqualTo(9.88).Within(1e-12), "last table entry 988 PAMs = 93.0 %");
            Assert.That(SequenceAligner.MuscleEngine.KimuraDist(0.05), Is.EqualTo(10.0), "over 93 %: 10");
        });
    }

    /// <summary>The test matrix of upgma2.cpp (DistCalcTest): ((0,1),(2,3)),4) with average linkage; edge lengths = height differences, heights = half the cluster distance.</summary>
    [Test]
    public void Upgma_OnTheSourceTestMatrix_BuildsTheExpectedTree()
    {
        var e = Engine(new MuscleOptions { SequenceType = MuscleSequenceType.Protein, Cluster1 = MuscleClusterMethod.Upgma }, "AAAAAA", "CCCCCC", "DDDDDD", "EEEEEE", "FFFFFF");
        float[,] d =
        {
            { 0, 2, 14, 14, 20 }, { 2, 0, 14, 14, 20 }, { 14, 14, 0, 4, 20 }, { 14, 14, 4, 0, 20 }, { 20, 20, 20, 20, 0 },
        };
        var t = e.Upgma(d, MuscleClusterMethod.Upgma, [0, 1, 2, 3, 4]);
        Assert.Multiple(() =>
        {
            Assert.That((t.Left[5], t.Right[5]), Is.EqualTo((0, 1)));
            Assert.That((t.Left[6], t.Right[6]), Is.EqualTo((2, 3)));
            Assert.That((t.Left[7], t.Right[7]), Is.EqualTo((5, 6)));
            Assert.That((t.Left[8], t.Right[8]), Is.EqualTo((7, 4)));
            Assert.That(t.Root, Is.EqualTo(8));
            Assert.That(new[] { t.EdgeLength[0], t.EdgeLength[1], t.EdgeLength[2], t.EdgeLength[3], t.EdgeLength[5], t.EdgeLength[6], t.EdgeLength[7], t.EdgeLength[4] },
                Is.EqualTo(new[] { 1.0, 1.0, 2.0, 2.0, 6.0, 5.0, 3.0, 10.0 }));
            Assert.That(SequenceAligner.MuscleEngine.Newick(t), Is.EqualTo("(((0,1),(2,3)),4);"));
        });
    }

    /// <summary>UPGMB = 0.1 x average + 0.9 x minimum (g_dSUEFF) for the distance of a new cluster, ties broken by the lowest row.</summary>
    [Test]
    public void Upgmb_UsesTheBiasedLinkage()
    {
        var e = Engine(new MuscleOptions { SequenceType = MuscleSequenceType.Protein }, "AAAAAA", "CCCCCC", "DDDDDD", "EEEEEE");
        float[,] d = { { 0, 2, 10, 20 }, { 2, 0, 30, 40 }, { 10, 30, 0, 6 }, { 20, 40, 6, 0 } };
        var t = e.Upgma(d, MuscleClusterMethod.Upgmb, [0, 1, 2, 3]);
        // merge (0,1) at height 1; distance {0,1}-2 = 0.1 * 20 + 0.9 * 10 = 11, {0,1}-3 = 0.1 * 30 + 0.9 * 20 = 21; (2,3) = 6 is the next minimum: height 3
        Assert.Multiple(() =>
        {
            Assert.That((t.Left[4], t.Right[4]), Is.EqualTo((0, 1)));
            Assert.That((t.Left[5], t.Right[5]), Is.EqualTo((2, 3)));
            Assert.That(t.EdgeLength[2], Is.EqualTo(3.0));
            // {0,1}-{2,3}: 0.1 * avg(11, 21) + 0.9 * min(11, 21) = 1.6 + 9.9 = 11.5, height 5.75
            Assert.That(t.EdgeLength[4], Is.EqualTo(5.75 - 1.0).Within(1e-6));
            Assert.That(t.EdgeLength[5], Is.EqualTo(5.75 - 3.0).Within(1e-6));
        });
    }

    /// <summary>clwwt.cpp's worked example: edge strength = length / leaves under the edge, weight = sum along the root path (A 0.2, B 0.2, C 1.0, D 1.1), then normalised to sum 1.</summary>
    [Test]
    public void ClustalWWeights_OnTheSourceExample()
    {
        // x (root) -> A (0.2), y (0.3); y -> B (0.1), z (0.4); z -> C (0.7), D (0.8); leaves 0..3 = A, B, C, D; internal z = 4, y = 5, x = 6
        var t = new SequenceAligner.MuscleTree
        {
            LeafCount = 4,
            Root = 6,
            Left = [-1, -1, -1, -1, 2, 1, 0],
            Right = [-1, -1, -1, -1, 3, 4, 5],
            Parent = [6, 5, 4, 4, 5, 6, -1],
            EdgeLength = [0.2, 0.1, 0.7, 0.8, 0.4, 0.3, 0.0],
            LeafId = [0, 1, 2, 3, 0, 0, 0],
        };
        var w = SequenceAligner.MuscleEngine.CalcClustalWWeights(t);
        Assert.That(w, Is.EqualTo(new[] { 0.2f / 2.5f, 0.2f / 2.5f, 1.0f / 2.5f, 1.1f / 2.5f }).Within(1e-6));
        Assert.That(w.Sum(), Is.EqualTo(1f).Within(1e-6));
    }

    [Test]
    public void KmerDistance_IsZeroForIdenticalAndThreeWhenNoSextetIsShared()
    {
        // 6-mers over the six residue groups (A G P S T / I L M V / N D Q E B Z / R H K / F W Y / C)
        var e = Engine(new MuscleOptions { SequenceType = MuscleSequenceType.Protein }, "AGPSTAGPSTAGPST", "AGPSTAGPSTAGPST", "CCCCCCCCCCCCCCC", "FWYFWYFWYFWYFWY");
        var d = e.ComputeKmerDistances();
        Assert.Multiple(() =>
        {
            Assert.That(d[1, 0], Is.EqualTo(0f));
            Assert.That(d[2, 0], Is.EqualTo(3f), "no common sextet: 3 x (c11 - 0) / c11");
            Assert.That(d[3, 2], Is.EqualTo(3f));
            Assert.That(d[0, 2], Is.EqualTo(d[2, 0]));
        });
    }

    // ---- the dynamic program ------------------------------------------------------------------------------------------------------------------------

    private static float PathScore(SequenceAligner.MuscleEngine e, SequenceAligner.MuscleProfPos[] a, SequenceAligner.MuscleProfPos[] b, List<(char Type, int A, int B)> path, float ext)
    {
        double score = 0;
        int i = 0;
        while (i < path.Count)
        {
            char t = path[i].Type;
            if (t == 'M')
            {
                score += e.ScoreProfPos2(a[path[i].A - 1], b[path[i].B - 1]);
                i++;
                continue;
            }

            int j = i;
            while (j + 1 < path.Count && path[j + 1].Type == t)
                j++;
            var prof = t == 'D' ? a : b;
            int first = (t == 'D' ? path[i].A : path[i].B) - 1, last = (t == 'D' ? path[j].A : path[j].B) - 1;
            score += prof[first].GapOpen + (j - i) * ext;
            if (j + 1 < path.Count)
                score += prof[last].GapClose;   // a gap run that is followed by a match pays the close penalty; a trailing run does not
            i = j + 1;
        }

        return (float)score;
    }

    private static IEnumerable<List<(char Type, int A, int B)>> AllPaths(int la, int lb)
    {
        IEnumerable<List<(char, int, int)>> Rec(int ia, int ib, char prev, List<(char, int, int)> cur)
        {
            if (ia == la && ib == lb)
            {
                yield return new List<(char, int, int)>(cur);
                yield break;
            }

            if (ia < la && ib < lb)
            {
                cur.Add(('M', ia + 1, ib + 1));
                foreach (var p in Rec(ia + 1, ib + 1, 'M', cur)) yield return p;
                cur.RemoveAt(cur.Count - 1);
            }

            if (ia < la && prev != 'I')
            {
                cur.Add(('D', ia + 1, ib));
                foreach (var p in Rec(ia + 1, ib, 'D', cur)) yield return p;
                cur.RemoveAt(cur.Count - 1);
            }

            if (ib < lb && prev != 'D')
            {
                cur.Add(('I', ia, ib + 1));
                foreach (var p in Rec(ia, ib + 1, 'I', cur)) yield return p;
                cur.RemoveAt(cur.Count - 1);
            }
        }

        return Rec(0, 0, 'M', []);
    }

    /// <summary>
    /// nwsmall.cpp NWSmall is a Gotoh dynamic program with position-specific open / close penalties; its path must be a best alignment under that model. Checked against the exhaustive
    /// enumeration of every alignment (no D-I adjacency: a gap run in the other sequence needs a match in between) on random profiles of 2-5 columns (single sequences and 2-3 sequence profiles,
    /// protein LE and nucleotide SPN, extension 0 and non-zero).
    /// </summary>
    [Test]
    public void NwSmall_FindsAnOptimalPath_ComparedWithExhaustiveEnumeration()
    {
        var rnd = new Random(20261128);
        const string aa = "ACDEFGHIKLMNPQRSTVWY";
        int checkedCases = 0;
        foreach (var (type, alphabet, o) in new[]
                 {
                     ("protein", aa, new MuscleOptions { SequenceType = MuscleSequenceType.Protein }),
                     ("protein-ext", aa, new MuscleOptions { SequenceType = MuscleSequenceType.Protein, GapExtend = -0.4 }),
                     ("dna", "ACGT", new MuscleOptions { SequenceType = MuscleSequenceType.Dna }),
                     ("dna-ext", "ACGT", new MuscleOptions { SequenceType = MuscleSequenceType.Dna, GapExtend = -30 }),
                 })
        {
            for (int trial = 0; trial < 120; trial++)
            {
                var e = Engine(o, type.StartsWith("dna") ? ["ACGTACGT", "ACGTACGT"] : ["ACDEFGHIK", "ACDEFGHIK"]);
                int la = rnd.Next(2, 6), lb = rnd.Next(1, 6);
                SequenceAligner.MuscleProfPos[] Make(int len, int seqs)
                {
                    var rows = new char[seqs][];
                    for (int s = 0; s < seqs; s++)
                    {
                        rows[s] = new char[len];
                        for (int c = 0; c < len; c++)
                            rows[s][c] = rnd.NextDouble() < 0.15 ? '-' : alphabet[rnd.Next(alphabet.Length)];
                    }

                    var w = Enumerable.Repeat(1f / seqs, seqs).ToArray();
                    return e.ProfileFromRows(rows, w);
                }

                var pa = Make(la, rnd.Next(1, 4));
                var pb = Make(lb, rnd.Next(1, 4));
                var path = e.NwSmall(pa, la, pb, lb);   // also applies SetTermGaps to both profiles, in place
                float ext = (float)(o.GapExtend ?? 0);
                float got = PathScore(e, pa, pb, path, ext);
                float best = AllPaths(la, lb).Max(p => PathScore(e, pa, pb, p, ext));
                Assert.That(got, Is.EqualTo(best).Within(1e-3 * Math.Max(1f, Math.Abs(best))), $"{type} trial {trial} la={la} lb={lb}");
                checkedCases++;
            }
        }

        Assert.That(checkedCases, Is.EqualTo(480));
    }

    [Test]
    public void NwSmall_AlignsIdenticalProfilesWithoutGaps_AndPlacesASingleDeletionAsOneGapRun()
    {
        var e = Engine();
        var a = e.ProfileFromRows([MuscleRow("ACDEFGHIKLMN")], [1f]);
        var same = e.ProfileFromRows([MuscleRow("ACDEFGHIKLMN")], [1f]);
        Assert.That(e.NwSmall(a, a.Length, same, same.Length).All(x => x.Type == 'M'));

        var shorter = e.ProfileFromRows([MuscleRow("ACDEFGHKLMN")], [1f]);   // I deleted
        var path = e.NwSmall(a, a.Length, shorter, shorter.Length);
        Assert.Multiple(() =>
        {
            Assert.That(path.Count(x => x.Type == 'D'), Is.EqualTo(1));
            Assert.That(path.Count(x => x.Type == 'I'), Is.Zero);
            Assert.That(path.Count, Is.EqualTo(12));
        });
    }

    // ---- the alignment ------------------------------------------------------------------------------------------------------------------------------

    private static string[] Align(IEnumerable<string> seqs, MuscleOptions? o = null) =>
        SequenceAligner.MultipleAlignMuscle(seqs, o ?? new MuscleOptions { SequenceType = MuscleSequenceType.Protein }).AlignedSequences.ToArray();

    [Test]
    public void MultipleAlignMuscle_Rows_SpellTheInputs_HaveEqualLength_AndAreDeterministic()
    {
        foreach (int iters in new[] { 1, 2 })
        {
            var rows = Align(Family, new MuscleOptions { SequenceType = MuscleSequenceType.Protein, MaxIterations = iters });
            Assert.Multiple(() =>
            {
                Assert.That(rows.Select(r => r.Length).Distinct().Count(), Is.EqualTo(1));
                for (int i = 0; i < Family.Length; i++)
                    Assert.That(rows[i].Replace("-", ""), Is.EqualTo(Family[i]), $"row {i}");
                Assert.That(Align(Family, new MuscleOptions { SequenceType = MuscleSequenceType.Protein, MaxIterations = iters }), Is.EqualTo(rows));
            });
        }
    }

    [Test]
    public void MultipleAlignMuscle_ProteinAndDnaOverloads_EqualTheStringOverload()
    {
        var dna = new[] { "ACGTACGTTAGCTAGCTAGGATCGATCGGATC", "ACGTACGTAGCTAGCTAGGATCGATCGATC", "ACGTTCGTTAGCTAGCTAGGATCGGATCGGATC", "ACGTACGTTAGCTTAGCTAGGATCGATCGGATC" };
        var viaDna = SequenceAligner.MultipleAlignMuscle(dna.Select(s => new DnaSequence(s)), new MuscleOptions { MaxIterations = 1 });
        var viaString = SequenceAligner.MultipleAlignMuscle(dna, new MuscleOptions { SequenceType = MuscleSequenceType.Dna, MaxIterations = 1 });
        var viaAuto = SequenceAligner.MultipleAlignMuscle(dna, new MuscleOptions { MaxIterations = 1 });   // GuessAlpha: >= 95 % of ACGTN
        var viaProt = SequenceAligner.MultipleAlignMuscle(Family.Select(s => new ProteinSequence(s)), new MuscleOptions { MaxIterations = 2 });
        Assert.Multiple(() =>
        {
            Assert.That(viaDna.AlignedSequences, Is.EqualTo(viaString.AlignedSequences));
            Assert.That(viaAuto.AlignedSequences, Is.EqualTo(viaString.AlignedSequences));
            Assert.That(viaProt.AlignedSequences, Is.EqualTo(Align(Family)));
        });
    }

    [Test]
    public void MultipleAlignMuscle_TwoSequences_StopAfterTheProgressiveStage_AndOneSequenceIsReturned()
    {
        var two = Family.Take(2).ToArray();
        Assert.That(Align(two, new MuscleOptions { SequenceType = MuscleSequenceType.Protein, MaxIterations = 2 }),
            Is.EqualTo(Align(two, new MuscleOptions { SequenceType = MuscleSequenceType.Protein, MaxIterations = 1 })));
        Assert.That(Align([Family[0]]), Is.EqualTo(new[] { Family[0] }));
        Assert.That(SequenceAligner.MultipleAlignMuscle(Array.Empty<string>()).AlignedSequences, Is.Empty);
    }

    [Test]
    public void MultipleAlignMuscle_MethionineHack_LetterRules_AndTreeOrder()
    {
        // a leading M is aligned as X and restored; J / O / U are replaced by X (MUSCLE: "invalid letters found"); gaps, white space and case are dropped / folded
        var rows = Align(["MKVLAAGIVG", "mkvl-aagivg", "KVLAAJIVG", "M K V L A A O I V G"]);
        Assert.Multiple(() =>
        {
            Assert.That(rows.Select(r => r.Length).Distinct().Count(), Is.EqualTo(1));
            Assert.That(rows[0].Replace("-", ""), Is.EqualTo("MKVLAAGIVG"));
            Assert.That(rows[1].Replace("-", ""), Is.EqualTo("MKVLAAGIVG"));
            Assert.That(rows[2].Replace("-", ""), Is.EqualTo("KVLAAXIVG"));
            Assert.That(rows[3].Replace("-", ""), Is.EqualTo("MKVLAAXIVG"));
            Assert.That(rows[0].TrimStart('-')[0], Is.EqualTo('M'));
        });

        var tree = SequenceAligner.MultipleAlignMuscle(Family, new MuscleOptions { SequenceType = MuscleSequenceType.Protein, TreeOrder = true }).AlignedSequences;
        Assert.That(tree.Select(r => r.Replace("-", "")), Is.EquivalentTo(Family));
    }

    [Test]
    public void MultipleAlignMuscle_Iteration2_RealignsOnlyChangedSubtrees_ButKeepsEveryRowValid()
    {
        var rnd = new Random(7);
        int changed = 0;
        for (int trial = 0; trial < 60; trial++)
        {
            string baseSeq = new(Enumerable.Range(0, 60).Select(_ => "ACDEFGHIKLMNPQRSTVWY"[rnd.Next(20)]).ToArray());
            var seqs = Enumerable.Range(0, 5 + rnd.Next(5)).Select(_ => new string(baseSeq.Where(_ => rnd.NextDouble() > 0.1).ToArray())).ToArray();
            var o1 = new MuscleOptions { SequenceType = MuscleSequenceType.Protein, MaxIterations = 1 };
            var r1 = Align(seqs, o1);
            var r2 = Align(seqs, o1 with { MaxIterations = 2 });
            var none = Align(seqs, o1 with { MaxIterations = 2, MaxTrees = 0 });
            Assert.That(none, Is.EqualTo(r1), "MaxTrees = 0 switches the tree refinement off");
            for (int i = 0; i < seqs.Length; i++)
                Assert.That(r2[i].Replace("-", ""), Is.EqualTo(seqs[i]), $"trial {trial} row {i}");
            if (!r1.SequenceEqual(r2))
                changed++;
        }

        Assert.That(changed, Is.GreaterThan(0), "iteration 2 changes some draft alignments");
    }

    [Test]
    public void MultipleAlignMuscle_Options_AreValidated()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<NotSupportedException>(() => Align(Family, new MuscleOptions { SequenceType = MuscleSequenceType.Protein, MaxIterations = 3 }));
            Assert.Throws<ArgumentException>(() => Align(Family, new MuscleOptions { SequenceType = MuscleSequenceType.Protein, MaxIterations = 0 }));
            Assert.Throws<ArgumentException>(() => Align(Family, new MuscleOptions { SequenceType = MuscleSequenceType.Protein, ScoreModel = MuscleScoreModel.SumOfPairsNucleotide }));
            Assert.Throws<ArgumentException>(() => Align(Family, new MuscleOptions { SequenceType = MuscleSequenceType.Protein, Distance1 = MuscleDistance.Kmer4_6 }));
            Assert.Throws<ArgumentException>(() => Align(["ACGTACGT", "ACGTACGA"], new MuscleOptions { SequenceType = MuscleSequenceType.Dna, ScoreModel = MuscleScoreModel.LogExpectation }));
            Assert.Throws<ArgumentException>(() => Align(["ACD", ""]));
            Assert.Throws<ArgumentException>(() => Align(["ACD", "AC1D"]));
        });
    }

    /// <summary>The SP (PAM200), SV (VTML240) and RNA parameter sets produce valid alignments with their own gap parameters, and a more expensive gap opening never adds gaps.</summary>
    [Test]
    public void MultipleAlignMuscle_ScoreModels_AndGapParameters_ProduceValidAlignments()
    {
        foreach (var model in new[] { MuscleScoreModel.LogExpectation, MuscleScoreModel.SumOfPairs, MuscleScoreModel.SumOfPairsVtml })
        {
            var rows = Align(Family, new MuscleOptions { SequenceType = MuscleSequenceType.Protein, ScoreModel = model });
            for (int i = 0; i < Family.Length; i++)
                Assert.That(rows[i].Replace("-", ""), Is.EqualTo(Family[i]), model.ToString());
        }

        int Gaps(double open) => Align(Family, new MuscleOptions { SequenceType = MuscleSequenceType.Protein, MaxIterations = 1, GapOpen = open }).Sum(r => r.Count(c => c == '-'));
        Assert.That(Gaps(-50), Is.LessThanOrEqualTo(Gaps(-0.5)));

        var rna = SequenceAligner.MultipleAlignMuscle(["ACGUACGUUAGCUAGC", "ACGUACGUAGCUAGC", "ACGUUCGUUAGCUAGC"], new MuscleOptions { SequenceType = MuscleSequenceType.Rna });
        Assert.That(rna.AlignedSequences.Select(r => r.Replace("-", "")), Is.EqualTo(new[] { "ACGUACGUUAGCUAGC", "ACGUACGUAGCUAGC", "ACGUUCGUUAGCUAGC" }));
    }
}
