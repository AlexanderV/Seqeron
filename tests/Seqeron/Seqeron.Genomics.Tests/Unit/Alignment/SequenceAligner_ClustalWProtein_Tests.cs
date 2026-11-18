using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A1-15b part b, B09 F62): ClustalW 2.1's PROTEIN pairwise stages — the full dynamic-programming stage
/// (<see cref="SequenceAligner.ClustalWPairwiseAlign(string, string, ClustalWProteinPairwiseOptions)"/>: GONNET250 / BLOSUM30 / PAM350 / identity matrix,
/// length-dependent gap opening, <c>-PWGAPOPEN -PWGAPEXT</c>) and the fast k-tuple stage
/// (<see cref="SequenceAligner.ClustalWQuickPairwiseAlign(string, string, ClustalWSequenceType, ClustalWQuickOptions?, int?)"/>), and their use as the guide-tree
/// distance of the protein multiple aligners (<see cref="MsaDistanceMode.ClustalWAffine"/> / <see cref="MsaDistanceMode.ClustalWQuick"/>).
///
/// Oracles (TestData/Alignment/clustalw_protein_*.tsv, msa_matrix_clustalw_oracle.tsv, scripts/alignment/generate_clustalw_protein_pairwise_oracle.py):
/// the installed <c>clustalw</c> 2.1 <c>-TYPE=PROTEIN -ALIGN</c> "Score:" lines (4 204 pairs, four matrices, ten gap-penalty variants, B/Z/X residues,
/// low-complexity and repeat pairs; every row also equals a line-by-line Python port with the matrices.h tables), <c>-QUICKTREE -ALIGN</c> (3 774 pairs),
/// the canonical references of the matrix tables (EMBOSS EBLOSUM30, Biopython GONNET1992 × 10) and independent Python MSA rows driven by the ClustalW distances.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_ClustalWProtein_Tests
{
    private static IEnumerable<string[]> ReadTsv(string resource)
    {
        string name = "Seqeron.Genomics.Tests.TestData.Alignment." + resource;
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        string? line;
        while ((line = reader.ReadLine()) is not null)
            if (line.Length > 0)
                yield return line.Split('\t');
    }

    private static readonly Lazy<List<string[]>> FullRows = new(() => ReadTsv("clustalw_protein_pairwise_oracle.tsv").ToList());
    private static readonly Lazy<List<string[]>> QuickRows = new(() => ReadTsv("clustalw_protein_quick_oracle.tsv").ToList());
    private static readonly Lazy<List<string[]>> MsaRows = new(() => ReadTsv("msa_matrix_clustalw_oracle.tsv").ToList());

    private static double D(string s) => double.Parse(s, CultureInfo.InvariantCulture);

    private static ClustalWPairwiseMatrix MatrixOf(string name) => name switch
    {
        "GONNET" => ClustalWPairwiseMatrix.Gonnet,
        "BLOSUM" => ClustalWPairwiseMatrix.Blosum,
        "PAM" => ClustalWPairwiseMatrix.Pam,
        "ID" => ClustalWPairwiseMatrix.Identity,
        _ => throw new ArgumentException(name),
    };

    private static ClustalWProteinPairwiseOptions FullOptions(string[] f) => new() { Matrix = MatrixOf(f[0]), GapOpen = D(f[1]), GapExtend = D(f[2]) };

    #region Full stage

    [Test]
    public void FullStage_PrintedScoreIdentitiesAndDpScore_MatchClustalW21Oracle()
    {
        var rows = FullRows.Value;
        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(4200));
        int bad = 0;
        string? first = null;
        foreach (var f in rows)
        {
            var r = SequenceAligner.ClustalWPairwiseAlign(f[3], f[4], FullOptions(f));
            bool ok = r.ClustalWScore == int.Parse(f[5]) && r.Identities == int.Parse(f[6]) && r.Score == int.Parse(f[7]);
            if (!ok && ++bad == 1)
                first = $"{string.Join('/', f.Take(3))} {f[3]} / {f[4]}: expected score {f[5]} identities {f[6]} dp {f[7]}, got {r.ClustalWScore} / {r.Identities} / {r.Score}";
        }

        Assert.That(bad, Is.Zero, first);
    }

    /// <summary>
    /// The returned local alignment (the Myers-Miller edit script of <c>diff</c>, exactly as ClustalW's <c>tracePath</c> walks it) never re-scores above
    /// the dynamic-programming optimum computed with ClustalW's matrix × intScale and the oracle's integer gap penalties, and it re-scores TO it whenever the
    /// path is gap-free or opens at most one gap run. (With the very large protein gap-opening penalties ClustalW's <c>diff</c> can return a path that
    /// scores below the forward optimum on unrelated pairs; the printed identity score — which is what ClustalW reports and what the oracle fixes — is
    /// reproduced regardless, see <see cref="FullStage_PrintedScoreIdentitiesAndDpScore_MatchClustalW21Oracle"/>.)
    /// </summary>
    [Test]
    public void FullStage_AlignedRows_RescoreBelowOrToTheDpScore()
    {
        var tables = new Dictionary<string, (SubstitutionMatrix Table, int Scale)>
        {
            ["GONNET"] = (SubstitutionMatrix.ClustalWGonnet250, 10),
            ["BLOSUM"] = (SubstitutionMatrix.ClustalWBlosum30, 100),
            ["PAM"] = (SubstitutionMatrix.ClustalWPam350, 100),
            ["ID"] = (SubstitutionMatrix.ClustalWIdentity, 100),
        };
        int checkedRows = 0, equalRows = 0, simple = 0, simpleEqual = 0;
        foreach (var f in FullRows.Value.Where((_, i) => i % 3 == 0))
        {
            var r = SequenceAligner.ClustalWPairwiseAlign(f[3], f[4], FullOptions(f));
            if (r.Score == 0)
                continue;
            var (table, scale) = tables[f[0]];
            int go = int.Parse(f[8]), ge = int.Parse(f[9]);
            int score = 0, runs = 0;
            bool inGapRow1 = false, inGapRow2 = false;
            for (int i = 0; i < r.AlignedSequence1.Length; i++)
            {
                char a = r.AlignedSequence1[i], b = r.AlignedSequence2[i];
                if (a == '-')
                {
                    if (!inGapRow1) runs++;
                    score -= (inGapRow1 ? 0 : go) + ge;
                    inGapRow1 = true;
                    inGapRow2 = false;
                }
                else if (b == '-')
                {
                    if (!inGapRow2) runs++;
                    score -= (inGapRow2 ? 0 : go) + ge;
                    inGapRow2 = true;
                    inGapRow1 = false;
                }
                else
                {
                    score += table.GetScore(a, b) * scale;
                    inGapRow1 = inGapRow2 = false;
                }
            }

            Assert.That(score, Is.LessThanOrEqualTo(r.Score), $"{f[0]} {f[3]} / {f[4]}");
            checkedRows++;
            if (score == r.Score) equalRows++;
            if (runs <= 1)
            {
                simple++;
                Assert.That(score, Is.EqualTo(r.Score), $"{f[0]} {f[3]} / {f[4]} (at most one gap run)");
                simpleEqual++;
            }
        }

        Assert.That(checkedRows, Is.GreaterThan(1000));
        Assert.That(simple, Is.GreaterThan(500), "gap-free / one-gap paths are the common case");
        Assert.That(simpleEqual, Is.EqualTo(simple));
        Assert.That(equalRows, Is.GreaterThanOrEqualTo(simple));
    }

    [Test]
    public void FullStage_OracleCoversMatricesGapsAndResidues()
    {
        var rows = FullRows.Value;
        Assert.Multiple(() =>
        {
            Assert.That(rows.Select(f => f[0]).Distinct().Order(), Is.EqualTo(new[] { "BLOSUM", "GONNET", "ID", "PAM" }));
            Assert.That(rows.Select(f => (f[1], f[2])).Distinct().Count(), Is.GreaterThanOrEqualTo(9));
            Assert.That(rows.Count(f => f[3].Contains('B') || f[3].Contains('Z') || f[3].Contains('X') || f[4].Contains('X')), Is.GreaterThan(100), "ambiguity residues");
            Assert.That(rows.Count(f => f[10] == "0"), Is.GreaterThan(500), "identity matrix: matrixAvgScore = 0 (the 2 * (int)(...) gap-open branch)");
            Assert.That(rows.Count(f => f[10] != "0"), Is.GreaterThan(3000));
            // the matrix changes the printed score on many pairs: equality is not an accident of easy inputs
            int differs = 0;
            foreach (var f in rows.Take(400))
            {
                var other = new ClustalWProteinPairwiseOptions { Matrix = ClustalWPairwiseMatrix.Identity, GapOpen = D(f[1]), GapExtend = D(f[2]) };
                if (f[0] != "ID" && SequenceAligner.ClustalWPairwiseAlign(f[3], f[4], other).ClustalWScore != int.Parse(f[5])) differs++;
            }

            Assert.That(differs, Is.GreaterThan(20), "the identity matrix gives other scores than the substitution matrices");
        });
    }

    /// <summary>Hand run through the binary: `clustalw -TYPE=PROTEIN -ALIGN` on three globin fragments prints 40 / 38 / 31 for every -PWMATRIX of this set.</summary>
    [Test]
    public void FullStage_HandRunExample_PrintsSameScores()
    {
        string a = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF";
        string b = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKV";
        string c = "MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKTEAEMKASEDLKKHG";
        var o = new ClustalWProteinPairwiseOptions();
        Assert.Multiple(() =>
        {
            Assert.That(SequenceAligner.ClustalWPairwiseAlign(a, b, o).ClustalWScore, Is.EqualTo(40));
            Assert.That(SequenceAligner.ClustalWPairwiseAlign(a, c, o).ClustalWScore, Is.EqualTo(38));
            Assert.That(SequenceAligner.ClustalWPairwiseAlign(b, c, o).ClustalWScore, Is.EqualTo(31));
            // distance 1 - identities / MIN(len)
            var r = SequenceAligner.ClustalWPairwiseAlign(a, b, o);
            Assert.That(r.Distance, Is.EqualTo(1.0 - (double)r.Identities / a.Length).Within(1e-12));
            Assert.That(r.PercentIdentity, Is.EqualTo(100.0 * r.Identities / a.Length).Within(1e-9));
        });
    }

    [Test]
    public void FullStage_Contracts()
    {
        var o = new ClustalWProteinPairwiseOptions();
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => SequenceAligner.ClustalWPairwiseAlign(null!, "MKV", o));
            Assert.Throws<ArgumentNullException>(() => SequenceAligner.ClustalWPairwiseAlign("MKV", null!, o));
            Assert.Throws<ArgumentNullException>(() => SequenceAligner.ClustalWPairwiseAlign("MKV", "MKV", (ClustalWProteinPairwiseOptions)null!));
            Assert.Throws<ArgumentException>(() => SequenceAligner.ClustalWPairwiseAlign("MKJV", "MKV", o), "J is dropped by clustalw");
            Assert.Throws<ArgumentException>(() => SequenceAligner.ClustalWPairwiseAlign("MK*V", "MKV", o), "'*' is dropped by clustalw");
            Assert.Throws<ArgumentException>(() => SequenceAligner.ClustalWPairwiseAlign("MK1V", "MKV", o));
            Assert.Throws<ArgumentException>(() => SequenceAligner.ClustalWPairwiseAlign("MKV", "MKV", o with { GapOpen = -1 }));
            Assert.Throws<ArgumentException>(() => SequenceAligner.ClustalWPairwiseAlign("MKV", "MKV", o with { GapExtend = double.NaN }));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.ClustalWPairwiseAlign("MKV", "MKV", o with { Matrix = (ClustalWPairwiseMatrix)9 }));
            // empty: distance 1, no alignment (as FullPairwiseAlign)
            var empty = SequenceAligner.ClustalWPairwiseAlign("", "MKV", o);
            Assert.That(empty.Distance, Is.EqualTo(1.0));
            Assert.That(empty.Identities, Is.Zero);
            // case-insensitive, ProteinSequence overload, O and U are residues without table row
            var lower = SequenceAligner.ClustalWPairwiseAlign("mkvlaagivg", "MKVLSAGIVG", o);
            var upper = SequenceAligner.ClustalWPairwiseAlign("MKVLAAGIVG", "MKVLSAGIVG", o);
            Assert.That(lower, Is.EqualTo(upper));
            Assert.That(SequenceAligner.ClustalWPairwiseAlign(new ProteinSequence("MKVLAAGIVG"), new ProteinSequence("MKVLSAGIVG"), o), Is.EqualTo(upper));
            Assert.DoesNotThrow(() => SequenceAligner.ClustalWPairwiseAlign("MKVLOUAGIV", "MKVLOUAGIV", o));
            // defaults of ClustalW's PROTEIN stage
            Assert.That(o.Matrix, Is.EqualTo(ClustalWPairwiseMatrix.Gonnet));
            Assert.That(o.GapOpen, Is.EqualTo(10.0));
            Assert.That(o.GapExtend, Is.EqualTo(0.1));
            Assert.That(o.Quick, Is.EqualTo(new ClustalWQuickOptions { KTuple = 1, Window = 5, TopDiagonals = 5, PairGap = 3 }));
            Assert.That(SequenceAligner.ClustalWProteinPairwiseGapOpen, Is.EqualTo(10.0));
        });
    }

    #endregion

    #region Matrix tables

    [Test]
    public void MatrixTables_EqualTheirCanonicalReferences()
    {
        var rows = ReadTsv("clustalw_protein_matrices_reference.tsv").ToList();
        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(929), "BLOSUM30 23 x 23 (EMBOSS EBLOSUM30) + GONNET250 20 x 20 (round(10 x Biopython GONNET1992))");
        foreach (var f in rows)
        {
            var m = f[0] == "BLOSUM30" ? SubstitutionMatrix.ClustalWBlosum30 : SubstitutionMatrix.ClustalWGonnet250;
            Assert.That(m.GetScore(f[1][0], f[2][0]), Is.EqualTo(int.Parse(f[3])), $"{f[0]} {f[1]}{f[2]}");
        }
    }

    [Test]
    public void MatrixTables_AreSymmetricWithTheClustalWAlphabet()
    {
        foreach (var m in new[] { SubstitutionMatrix.ClustalWGonnet250, SubstitutionMatrix.ClustalWBlosum30, SubstitutionMatrix.ClustalWPam350, SubstitutionMatrix.ClustalWIdentity })
        {
            Assert.That(m.Alphabet, Is.EqualTo("ABCDEFGHIKLMNPQRSTVWXYZ"), m.Name);
            foreach (char a in m.Alphabet)
                foreach (char b in m.Alphabet)
                    Assert.That(m.GetScore(a, b), Is.EqualTo(m.GetScore(b, a)), $"{m.Name} {a}{b}");
        }

        var id = SubstitutionMatrix.ClustalWIdentity;
        Assert.Multiple(() =>
        {
            Assert.That(id.GetScore('A', 'A'), Is.EqualTo(10));
            Assert.That(id.GetScore('B', 'B'), Is.EqualTo(10));
            Assert.That(id.GetScore('A', 'R'), Is.Zero);
            Assert.That(id.GetScore('B', 'D'), Is.Zero);
            // anchor values of the published tables
            Assert.That(SubstitutionMatrix.ClustalWBlosum30.GetScore('W', 'W'), Is.EqualTo(20));
            Assert.That(SubstitutionMatrix.ClustalWBlosum30.GetScore('C', 'C'), Is.EqualTo(17));
            Assert.That(SubstitutionMatrix.ClustalWGonnet250.GetScore('A', 'A'), Is.EqualTo(24), "Gonnet A/A = 2.4 bits x 10");
            Assert.That(SubstitutionMatrix.ClustalWGonnet250.GetScore('W', 'W'), Is.EqualTo(142));
            Assert.That(SubstitutionMatrix.ClustalWPam350.GetScore('W', 'W'), Is.EqualTo(27));
        });
    }

    #endregion

    #region Fast stage

    private static ClustalWQuickOptions QuickOptions(string[] f) => new()
    {
        KTuple = int.Parse(f[1]),
        Window = int.Parse(f[2]),
        TopDiagonals = int.Parse(f[3]),
        PairGap = int.Parse(f[4]),
        Percent = f[5] == "1",
    };

    [Test]
    public void QuickStage_PrintedScore_MatchesClustalW21Oracle()
    {
        var rows = QuickRows.Value;
        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(3700));
        int bad = 0;
        string? first = null;
        foreach (var f in rows)
        {
            var r = SequenceAligner.ClustalWQuickPairwiseAlign(f[7], f[8], ClustalWSequenceType.Protein, QuickOptions(f), int.Parse(f[6]));
            bool ok;
            if (f[9] == "NA")
                ok = !r.Aligned;
            else
            {
                double expected = D(f[9]);
                ok = r.Aligned && Math.Abs(r.Score - expected) <= 5.1e-6 * Math.Max(1.0, Math.Abs(expected));
            }

            if (!ok && ++bad == 1)
                first = $"{string.Join('/', f.Take(7))} {f[7]} / {f[8]}: expected {f[9]}, got {r.Score} (aligned {r.Aligned}, fragment {r.FragmentScore})";
        }

        Assert.That(bad, Is.Zero, first);
    }

    [Test]
    public void QuickStage_DefaultOptions_MatchClustalWProteinDefaultRuns()
    {
        var rows = QuickRows.Value.Where(f => f[1] == "1" && f[2] == "5" && f[3] == "5" && f[4] == "3" && f[5] == "1").ToList();
        Assert.That(rows.Count, Is.GreaterThan(1500));
        int bad = 0;
        foreach (var f in rows)
        {
            var r = SequenceAligner.ClustalWQuickPairwiseAlign(f[7], f[8], ClustalWSequenceType.Protein, null, int.Parse(f[6]));
            bool ok = f[9] == "NA" ? !r.Aligned : r.Aligned && Math.Abs(r.Score - D(f[9])) <= 5.1e-6 * Math.Max(1.0, r.Score);
            if (!ok) bad++;
        }

        Assert.That(bad, Is.Zero);
    }

    [Test]
    public void QuickStage_OracleCoversTheParameterSpace()
    {
        var rows = QuickRows.Value;
        Assert.Multiple(() =>
        {
            Assert.That(rows.Select(f => f[0]).Distinct().Order(), Is.EqualTo(new[] { "long", "mut", "rand", "rep", "x" }));
            Assert.That(rows.Select(f => f[1]).Distinct().Order(), Is.EqualTo(new[] { "1", "2" }), "k-tuple 1..2 (the protein limit)");
            Assert.That(rows.Count(f => f[5] == "0"), Is.GreaterThan(300), "-SCORE=ABSOLUTE");
            Assert.That(rows.Count(f => f[9] == "NA"), Is.GreaterThan(300), "Not Aligned");
            Assert.That(rows.Count(f => f[4] != (int.Parse(f[1]) + 3).ToString()), Is.GreaterThan(300), "explicit -PAIRGAP");
        });
    }

    /// <summary>Hand run through the binary: `clustalw -TYPE=PROTEIN -QUICKTREE -ALIGN` prints 34.0426 / 38.2979 / 26.2295 for the three globin fragments.</summary>
    [Test]
    public void QuickStage_HandRunExample_PrintsSameScores()
    {
        string a = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF";
        string b = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKV";
        string c = "MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKTEAEMKASEDLKKHG";
        int longest = new[] { a.Length, b.Length, c.Length }.Max();
        Assert.Multiple(() =>
        {
            Assert.That(SequenceAligner.ClustalWQuickPairwiseAlign(a, b, ClustalWSequenceType.Protein, null, longest).Score, Is.EqualTo(34.0426).Within(5e-5));
            Assert.That(SequenceAligner.ClustalWQuickPairwiseAlign(a, c, ClustalWSequenceType.Protein, null, longest).Score, Is.EqualTo(38.2979).Within(5e-5));
            Assert.That(SequenceAligner.ClustalWQuickPairwiseAlign(b, c, ClustalWSequenceType.Protein, null, longest).Score, Is.EqualTo(26.2295).Within(5e-5));
            Assert.That(SequenceAligner.ClustalWQuickPairwiseAlign(new ProteinSequence(a), new ProteinSequence(b), null, longest).Score, Is.EqualTo(34.0426).Within(5e-5));
        });
    }

    [Test]
    public void QuickStage_Contracts_AndDnaDelegation()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => SequenceAligner.ClustalWQuickPairwiseAlign(null!, "MKV", ClustalWSequenceType.Protein));
            Assert.Throws<ArgumentException>(() => SequenceAligner.ClustalWQuickPairwiseAlign("MKJV", "MKV", ClustalWSequenceType.Protein));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.ClustalWQuickPairwiseAlign("MKV", "MKV", ClustalWSequenceType.Protein, new ClustalWQuickOptions { KTuple = 3 }),
                "the protein k-tuple is 1..2");
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.ClustalWQuickPairwiseAlign("MKV", "MKV", (ClustalWSequenceType)5));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.ClustalWQuickPairwiseAlign("MKVLAAGIVG", "MKVLAAGIVG", ClustalWSequenceType.Protein, null, 3));
            Assert.That(SequenceAligner.ClustalWQuickPairwiseAlign("", "MKV", ClustalWSequenceType.Protein).Distance, Is.EqualTo(1.0));
            // DNA type = the DNA stage (k-tuple 3 is valid there)
            var dna = SequenceAligner.ClustalWQuickPairwiseAlign("ACGTACGTTTGACC", "ACGTACGATTGACC", ClustalWSequenceType.Dna, new ClustalWQuickOptions { KTuple = 3, PairGap = 7 });
            Assert.That(dna, Is.EqualTo(SequenceAligner.ClustalWQuickPairwiseAlign("ACGTACGTTTGACC", "ACGTACGATTGACC", new ClustalWQuickOptions { KTuple = 3, PairGap = 7 })));
        });
    }

    #endregion

    #region Guide-tree distance of the protein multiple aligners

    private static ProteinSequence[] Proteins(string joined) => joined.Split(';').Select(s => new ProteinSequence(s)).ToArray();

    private static MatrixScoring Scoring(string[] f) =>
        new(f[1] == "PAM250" ? SubstitutionMatrix.Pam250 : SubstitutionMatrix.Blosum62, D(f[2]), D(f[3]));

    [Test]
    public void ClustalWAffineDistance_ProteinMsaRows_EqualTheIndependentPythonImplementation()
    {
        var rows = MsaRows.Value;
        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(60));
        var failures = new List<string>();
        foreach (var f in rows)
        {
            var seqs = Proteins(f[7]);
            var sc = Scoring(f);
            var options = new MsaOptions
            {
                GapModel = f[0] == "affine" ? MsaGapModel.Affine : MsaGapModel.Linear,
                Distance = MsaDistanceMode.ClustalWAffine,
                ClustalWProtein = new ClustalWProteinPairwiseOptions { Matrix = MatrixOf(f[4]), GapOpen = D(f[5]), GapExtend = D(f[6]) },
            };
            var got = new[]
            {
                SequenceAligner.MultipleAlignProgressive(seqs, sc, options).AlignedSequences,
                SequenceAligner.MultipleAlignIterative(seqs, sc, options).AlignedSequences,
                SequenceAligner.MultipleAlignConsistency(seqs, sc, options).AlignedSequences,
            };
            for (int m = 0; m < 3; m++)
                if (string.Join(';', got[m]) != f[8 + m])
                    failures.Add($"method {m} {f[0]} {f[4]} {f[5]}/{f[6]} {f[7]}");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count * 3}; first: {failures.FirstOrDefault()}");
    }

    /// <summary>The ClustalW distance changes the guide tree (and the rows) on many sets: the oracle equality above is not an accident of easy inputs.</summary>
    [Test]
    public void ClustalWAffineDistance_ChangesTheRowsOfManySets()
    {
        int differs = 0;
        foreach (var f in MsaRows.Value.Where(f => f[0] == "linear").Take(30))
        {
            var seqs = Proteins(f[7]);
            var sc = Scoring(f);
            var linear = SequenceAligner.MultipleAlignProgressive(seqs, sc, new MsaOptions { GapModel = MsaGapModel.Linear });
            if (string.Join(';', linear.AlignedSequences) != f[8])
                differs++;
        }

        Assert.That(differs, Is.GreaterThan(5));
    }

    [Test]
    public void QuickAndNeighbourJoining_ProteinMsa_AreValid_AndUseTheProteinStage()
    {
        var rows = MsaRows.Value.Take(24).ToList();
        int quickDiffersFromLinear = 0, njValid = 0;
        foreach (var f in rows)
        {
            var seqs = Proteins(f[7]);
            var sc = Scoring(f);
            var gap = f[0] == "affine" ? MsaGapModel.Affine : MsaGapModel.Linear;
            var quick = SequenceAligner.MultipleAlignProgressive(seqs, sc, new MsaOptions { GapModel = gap, Distance = MsaDistanceMode.ClustalWQuick });
            var linear = SequenceAligner.MultipleAlignProgressive(seqs, sc, new MsaOptions { GapModel = gap });
            var nj = SequenceAligner.MultipleAlignIterative(seqs, sc,
                new MsaOptions { GapModel = gap, Distance = MsaDistanceMode.ClustalWAffine, GuideTree = MsaGuideTree.NeighborJoining });
            var njQuick = SequenceAligner.MultipleAlignConsistency(seqs, sc,
                new MsaOptions { GapModel = gap, Distance = MsaDistanceMode.ClustalWQuick, GuideTree = MsaGuideTree.NeighborJoining });
            foreach (var r in new[] { quick, nj, njQuick })
            {
                Assert.That(r.AlignedSequences.Select(s => s.Replace("-", "")), Is.EqualTo(seqs.Select(s => s.Sequence)));
                Assert.That(r.AlignedSequences.Select(s => s.Length).Distinct().Count(), Is.EqualTo(1));
            }

            njValid++;
            if (string.Join(';', quick.AlignedSequences) != string.Join(';', linear.AlignedSequences))
                quickDiffersFromLinear++;
        }

        Assert.That(njValid, Is.EqualTo(rows.Count));
        Assert.That(quickDiffersFromLinear, Is.GreaterThan(0), "the fast stage gives other guide trees than the linear distance on some sets");
    }

    [Test]
    public void ProteinDistance_WithResiduesClustalWDoesNotRead_IsRejected_AndNucleotideMatricesUseTheDnaStage()
    {
        var sc = new MatrixScoring(SubstitutionMatrix.Blosum62, -9.5, -0.5);
        var clustalw = new MsaOptions { Distance = MsaDistanceMode.ClustalWAffine };
        // '*' is a BLOSUM62 symbol, but ClustalW's protein stage cannot read it
        Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignProgressive(["MKV*L", "MKVAL", "MKVVL"], sc, clustalw));
        // EDNAFULL = nucleotide matrix: the DNA stage (the existing behaviour for DNA) is used
        var edna = new MatrixScoring(SubstitutionMatrix.EdnaFull, -9.5, -0.5);
        string[] dnaSet = ["ACGTACGTTGCA", "ACGTACGTGCA", "ACTTACGTTGCAA", "ACGTTCGTTGA"];
        var viaMatrix = SequenceAligner.MultipleAlignProgressive(dnaSet, edna, clustalw);
        var viaDna = SequenceAligner.MultipleAlignProgressive(dnaSet.Select(s => new DnaSequence(s)).ToList(), new ScoringMatrix(1, -1, -2, -1), clustalw);
        Assert.That(viaMatrix.AlignedSequences.Select(s => s.Replace("-", "")), Is.EqualTo(viaDna.AlignedSequences.Select(s => s.Replace("-", ""))));
    }

    #endregion
}
