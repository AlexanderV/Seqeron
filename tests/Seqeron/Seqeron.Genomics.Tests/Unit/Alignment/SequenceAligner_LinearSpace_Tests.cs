using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-GLOBAL-001 / ALIGN-LOCAL-001 / ALIGN-SEMI-001 (audit A1-6, B09 F46): exact linear-space alignment
/// (<c>*LinearSpace</c> methods and the automatic switch above <see cref="SequenceAligner.LinearSpaceCellThreshold"/> cells).
///
/// Contract: score, aligned rows, coordinates and alignment type are IDENTICAL to the full-matrix aligner of the same model
/// (ties included — the traceback crossing of the middle row is carried through the forward pass), while the memory is O(m + n).
///
/// References: the full-matrix core (<c>forceFullMatrix</c>) row for row on thousands of random / tie-rich pairs with forced
/// deep recursion and on 2.5-5.3 kb pairs (including the automatic switch); Biopython 1.88 <c>PairwiseAligner</c> scores and
/// co-optimal row sets and parasail 1.3.4 (TestData/Alignment/linear_space_oracle.tsv, scripts/alignment/generate_linear_space_oracle.py);
/// EMBOSS 6.6.0 <c>stretcher</c> (Myers-Miller 1988 linear-space global) and <c>matcher</c> (Huang-Miller SIM local).
/// </summary>
[TestFixture]
public class SequenceAligner_LinearSpace_Tests
{
    // ---------------------------------------------------------------- reference plumbing

    private static readonly ScoringMatrix[] Schemes =
    {
        new(2, -3, 0, -2), new(2, -3, -3, -1), new(1, -1, -2, -1), new(5, -4, -10, -1), new(3, -2, 0, -4),
        new(1, -2, -5, -2), new(2, -3, -2, 0), new(1, 0, -1, -1), new(4, -5, 0, -6), new(1, -1, 0, 0),
        SequenceAligner.BlastDna, SequenceAligner.SimpleDna, SequenceAligner.HighIdentityDna,
    };

    /// <summary>The full-matrix aligner of the family (never the linear-space switch).</summary>
    private static AlignmentResult Full(string a, string b, ScoringMatrix sc, bool affine, bool local, EndGaps ends)
    {
        var ps = new PairScoring(sc);
        if (local)
            return affine
                ? SequenceAligner.AffineAlignCore(a, b, ps, AlignmentType.Local, EndGaps.None, forceFullMatrix: true)
                : SequenceAligner.LocalAlignCore(a, b, ps, forceFullMatrix: true);
        if (affine)
            return SequenceAligner.AffineAlignCore(a, b, ps, ends == EndGaps.None ? AlignmentType.Global : AlignmentType.SemiGlobal, ends, forceFullMatrix: true);
        return SequenceAligner.LinearAlignCore(a, b, ps, ends, default, null, forceFullMatrix: true);
    }

    /// <summary>The public linear-space entry point of the family.</summary>
    private static AlignmentResult Fast(string a, string b, ScoringMatrix sc, bool affine, bool local, EndGaps ends, long baseCells)
    {
        if (local)
            return affine
                ? SequenceAligner.LocalAlignAffineLinearSpace(a, b, sc, baseCells)
                : SequenceAligner.LocalAlignLinearSpace(a, b, sc, baseCells);
        if (ends == EndGaps.None)
            return affine
                ? SequenceAligner.GlobalAlignAffineLinearSpace(a, b, sc, baseCells)
                : SequenceAligner.GlobalAlignLinearSpace(a, b, sc, baseCells);
        return affine
            ? SequenceAligner.SemiGlobalAlignAffineLinearSpace(a, b, sc, ends, baseCells)
            : SequenceAligner.SemiGlobalAlignLinearSpace(a, b, sc, ends, baseCells);
    }

    private static string Describe(AlignmentResult r) =>
        $"{r.AlignmentType} score {r.Score} [{r.StartPosition1},{r.StartPosition2}..{r.EndPosition1},{r.EndPosition2}]\n  {r.AlignedSequence1}\n  {r.AlignedSequence2}";

    private static void AssertSame(AlignmentResult expected, AlignmentResult actual, string context) =>
        Assert.That(actual, Is.EqualTo(expected), () => $"{context}\nfull-matrix: {Describe(expected)}\nlinear-space: {Describe(actual)}");

    private static readonly Random Rnd = new(20261118);

    private static string RandomSeq(Random rnd, string alphabet, int length)
    {
        var sb = new StringBuilder(length);
        for (int i = 0; i < length; i++) sb.Append(alphabet[rnd.Next(alphabet.Length)]);
        return sb.ToString();
    }

    /// <summary>Related sequence: substitutions, single insertions/deletions and (optionally) planted long indels.</summary>
    private static string Mutate(Random rnd, string s, string alphabet, double rate, bool longGaps)
    {
        var sb = new StringBuilder();
        foreach (char c in s)
        {
            double r = rnd.NextDouble();
            if (r < rate * 0.5) sb.Append(alphabet[rnd.Next(alphabet.Length)]);
            else if (r < rate * 0.75) { /* deletion */ }
            else if (r < rate) { sb.Append(c); sb.Append(alphabet[rnd.Next(alphabet.Length)]); }
            else sb.Append(c);
        }
        string t = sb.ToString();
        if (longGaps && t.Length > 30)
        {
            int pos = rnd.Next(t.Length - 20);
            int len = rnd.Next(5, 20);
            t = rnd.Next(2) == 0 ? t.Remove(pos, Math.Min(len, t.Length - pos - 1)) : t.Insert(pos, RandomSeq(rnd, alphabet, len));
        }
        return t;
    }

    // ---------------------------------------------------------------- differential: identical rows, forced recursion

    private static readonly EndGaps[] AllPolicies = Enumerable.Range(0, 16).Select(i => (EndGaps)i).ToArray();

    [Test]
    public void LinearSpace_EqualsFullMatrix_RowForRow_AllFamiliesAndPolicies_ForcedRecursion()
    {
        var rnd = new Random(1701);
        string[] alphabets = { "A", "AC", "ACGT" };
        int checkedRows = 0;
        for (int pair = 0; pair < 45; pair++)
        {
            string alphabet = alphabets[pair % 3];
            string a = RandomSeq(rnd, alphabet, rnd.Next(0, 36));
            string b = pair % 4 == 0
                ? RandomSeq(rnd, alphabet, rnd.Next(0, 36))
                : Mutate(rnd, a, alphabet, rnd.Next(2) == 0 ? 0.15 : 0.5, longGaps: true);
            foreach (bool affine in new[] { false, true })
            {
                var sc = Schemes[rnd.Next(Schemes.Length)];
                long baseCells = new long[] { 1, 3, 9, 25, 70 }[rnd.Next(5)];
                // local (ends irrelevant) + every ends-free policy (0 = global)
                AssertSame(Full(a, b, sc, affine, true, EndGaps.None), Fast(a, b, sc, affine, true, EndGaps.None, baseCells),
                    $"local affine={affine} {sc} base={baseCells} '{a}' '{b}'");
                checkedRows++;
                foreach (var ends in AllPolicies)
                {
                    AssertSame(Full(a, b, sc, affine, false, ends), Fast(a, b, sc, affine, false, ends, baseCells),
                        $"ends={ends} affine={affine} {sc} base={baseCells} '{a}' '{b}'");
                    checkedRows++;
                }
            }
        }
        Assert.That(checkedRows, Is.EqualTo(45 * 2 * 17));
    }

    [Test]
    public void LinearSpace_EqualsFullMatrix_ManyRandomPairs_Global_Local_Fitting_Affine_Linear()
    {
        var rnd = new Random(1702);
        int rows = 0;
        for (int pair = 0; pair < 700; pair++)
        {
            string alphabet = new[] { "A", "AC", "ACGT", "ACGT" }[rnd.Next(4)];
            string a = RandomSeq(rnd, alphabet, rnd.Next(1, 60));
            string b = rnd.Next(5) == 0 ? RandomSeq(rnd, alphabet, rnd.Next(1, 60)) : Mutate(rnd, a, alphabet, 0.3, longGaps: true);
            var sc = Schemes[rnd.Next(Schemes.Length)];
            long baseCells = new long[] { 1, 2, 4, 16, 100 }[rnd.Next(5)];
            bool affine = rnd.Next(2) == 0;
            foreach (var (local, ends) in new[] { (false, EndGaps.None), (true, EndGaps.None), (false, EndGaps.Fitting), (false, EndGaps.All) })
            {
                AssertSame(Full(a, b, sc, affine, local, ends), Fast(a, b, sc, affine, local, ends, baseCells),
                    $"local={local} ends={ends} affine={affine} {sc} base={baseCells} '{a}' '{b}'");
                rows++;
            }
        }
        Assert.That(rows, Is.EqualTo(2800));
    }

    [Test]
    public void LinearSpace_MatrixScoring_EqualsFullMatrix_Protein_Dna_FractionalGaps()
    {
        var rnd = new Random(1703);
        const string protein = "ACDEFGHIKLMNPQRSTVWY";
        var gaps = new (double Open, double Extend)[] { (-9.5, -0.5), (-4.5, -0.5), (-2, -1), (0, -4), (-3.25, -0.25), (-11.5, -0.5), (0, 0) };
        int rows = 0;
        for (int pair = 0; pair < 60; pair++)
        {
            bool prot = pair % 2 == 0;
            string alphabet = prot ? protein : "ACGT";
            var matrix = prot ? SubstitutionMatrix.Blosum62 : SubstitutionMatrix.EdnaFull;
            string a = RandomSeq(rnd, alphabet, rnd.Next(1, 70));
            string b = Mutate(rnd, a, alphabet, 0.3, longGaps: true);
            var (open, extend) = gaps[rnd.Next(gaps.Length)];
            var ms = new MatrixScoring(matrix, open, extend);
            long baseCells = rnd.Next(2) == 0 ? 3 : 40;
            AssertSame(SequenceAligner.GlobalAlign(a, b, ms), SequenceAligner.GlobalAlignLinearSpace(a, b, ms, baseCells), $"global {open}/{extend}");
            AssertSame(SequenceAligner.GlobalAlignAffine(a, b, ms), SequenceAligner.GlobalAlignAffineLinearSpace(a, b, ms, baseCells), $"global affine {open}/{extend}");
            AssertSame(SequenceAligner.LocalAlign(a, b, ms), SequenceAligner.LocalAlignLinearSpace(a, b, ms, baseCells), $"local {open}/{extend}");
            AssertSame(SequenceAligner.LocalAlignAffine(a, b, ms), SequenceAligner.LocalAlignAffineLinearSpace(a, b, ms, baseCells), $"local affine {open}/{extend}");
            AssertSame(SequenceAligner.SemiGlobalAlign(a, b, ms), SequenceAligner.SemiGlobalAlignLinearSpace(a, b, ms, EndGaps.Fitting, baseCells), $"fitting {open}/{extend}");
            AssertSame(SequenceAligner.SemiGlobalAlignAffine(a, b, ms), SequenceAligner.SemiGlobalAlignAffineLinearSpace(a, b, ms, EndGaps.Fitting, baseCells), $"fitting affine {open}/{extend}");
            AssertSame(SequenceAligner.SemiGlobalAlignAffine(a, b, ms, EndGaps.All), SequenceAligner.SemiGlobalAlignAffineLinearSpace(a, b, ms, EndGaps.All, baseCells), $"all-free affine {open}/{extend}");
            rows += 7;
        }
        Assert.That(rows, Is.EqualTo(420));
    }

    [Test]
    public void LinearSpace_Overloads_DnaProteinString_Agree()
    {
        var a = new DnaSequence("ACGTTTTACGTACGGATTACA");
        var b = new DnaSequence("ACGTACGTACGATTTACA");
        var s = SequenceAligner.GlobalAlignAffineLinearSpace(a.Sequence, b.Sequence, SequenceAligner.BlastDna, 4);
        Assert.Multiple(() =>
        {
            Assert.That(SequenceAligner.GlobalAlignAffineLinearSpace(a, b, SequenceAligner.BlastDna, 4), Is.EqualTo(s));
            Assert.That(s, Is.EqualTo(SequenceAligner.GlobalAlignAffine(a, b, SequenceAligner.BlastDna)));
            var ms = new MatrixScoring(SubstitutionMatrix.EdnaFull, -11.5, -0.5);
            Assert.That(SequenceAligner.LocalAlignAffineLinearSpace(a, b, ms, 4), Is.EqualTo(SequenceAligner.LocalAlignAffine(a, b, ms)));
            Assert.That(SequenceAligner.SemiGlobalAlignLinearSpace(a, b, SequenceAligner.SimpleDna, EndGaps.Fitting, 4), Is.EqualTo(SequenceAligner.SemiGlobalAlign(a, b, SequenceAligner.SimpleDna)));
            var p1 = new ProteinSequence("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF");
            var p2 = new ProteinSequence("MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKV");
            var pm = MatrixScoring.Emboss(SubstitutionMatrix.Blosum62);
            Assert.That(SequenceAligner.GlobalAlignAffineLinearSpace(p1, p2, pm, 50), Is.EqualTo(SequenceAligner.GlobalAlignAffine(p1, p2, pm)));
            Assert.That(SequenceAligner.GlobalAlignLinearSpace(p1, p2, pm, 50), Is.EqualTo(SequenceAligner.GlobalAlign(p1, p2, pm)));
        });
    }

    // ---------------------------------------------------------------- large differential + automatic switch

    private static (string A, string B) RelatedPair(int seed, int length, double rate)
    {
        var rnd = new Random(seed);
        string a = RandomSeq(rnd, "ACGT", length);
        string b = Mutate(rnd, a, "ACGT", rate, longGaps: false);
        // two planted long indels (50-200 nt), as in resequenced / rearranged regions
        int p1 = rnd.Next(b.Length / 4, b.Length / 2), l1 = rnd.Next(50, 200);
        b = b.Remove(p1, l1);
        int p2 = rnd.Next(b.Length / 2, 3 * b.Length / 4);
        b = b.Insert(p2, RandomSeq(rnd, "ACGT", rnd.Next(50, 200)));
        return (a, b);
    }

    [TestCase(2600, 0.10, false, false)]
    [TestCase(2600, 0.10, true, false)]
    [TestCase(3000, 0.20, true, true)]
    [TestCase(3000, 0.20, false, true)]
    public void LinearSpace_KilobasePairs_EqualFullMatrix_AtDefaultBase(int length, double rate, bool affine, bool local)
    {
        var (a, b) = RelatedPair(length * 7 + (affine ? 1 : 0), length, rate);
        var sc = affine ? SequenceAligner.BlastDna : SequenceAligner.SimpleDna;
        foreach (var ends in local ? new[] { EndGaps.None } : new[] { EndGaps.None, EndGaps.Fitting, EndGaps.All, EndGaps.OverlapSuffix1Prefix2 })
        {
            var full = Full(a, b, sc, affine, local, ends);
            var fast = Fast(a, b, sc, affine, local, ends, SequenceAligner.LinearSpaceBaseCells);
            AssertSame(full, fast, $"{a.Length}x{b.Length} affine={affine} local={local} ends={ends}");
            // and with a deeper recursion
            AssertSame(full, Fast(a, b, sc, affine, local, ends, 50_000), $"{a.Length}x{b.Length} base 50000 ends={ends}");
        }
    }

    [Test]
    public void AutomaticSwitch_AboveThreshold_ReturnsTheFullMatrixAlignment()
    {
        // 5 300 x 5 200 = 27.6 M cells > LinearSpaceCellThreshold: the public aligners run the linear-space program on their own.
        var (a, b) = RelatedPair(99, 5300, 0.12);
        Assert.That((long)(a.Length + 1) * (b.Length + 1), Is.GreaterThan(SequenceAligner.LinearSpaceCellThreshold));
        var bd = SequenceAligner.BlastDna;
        Assert.Multiple(() =>
        {
            AssertSame(Full(a, b, bd, true, false, EndGaps.None), SequenceAligner.GlobalAlignAffine(a, b, bd), "GlobalAlignAffine");
            AssertSame(Full(a, b, bd, true, true, EndGaps.None), SequenceAligner.LocalAlignAffine(a, b, bd), "LocalAlignAffine");
            AssertSame(Full(a, b, bd, true, false, EndGaps.Fitting), SequenceAligner.SemiGlobalAlignAffine(a, b, bd), "SemiGlobalAlignAffine");
            AssertSame(Full(a, b, bd, true, false, EndGaps.All), SequenceAligner.SemiGlobalAlignAffine(a, b, bd, EndGaps.All), "SemiGlobalAlignAffine(All)");
            AssertSame(Full(a, b, bd, false, false, EndGaps.None), SequenceAligner.GlobalAlign(a, b, bd), "GlobalAlign");
            AssertSame(Full(a, b, bd, false, true, EndGaps.None), SequenceAligner.LocalAlign(a, b, bd), "LocalAlign");
            AssertSame(Full(a, b, bd, false, false, EndGaps.Fitting), SequenceAligner.SemiGlobalAlign(a, b, bd), "SemiGlobalAlign");
        });
    }

    [Test]
    public void AutomaticSwitch_MatrixScoring_ProteinAbove25MCells()
    {
        var rnd = new Random(2468);
        const string protein = "ACDEFGHIKLMNPQRSTVWY";
        string a = RandomSeq(rnd, protein, 5200);
        string b = Mutate(rnd, a, protein, 0.25, longGaps: false);
        var ms = MatrixScoring.Emboss(SubstitutionMatrix.Blosum62);
        var viaApi = SequenceAligner.GlobalAlignAffine(a, b, ms);
        var viaForced = SequenceAligner.GlobalAlignAffineLinearSpace(a, b, ms, 2_000_000);
        Assert.That(viaApi, Is.EqualTo(viaForced));
        Assert.That(viaApi.ExactScore, Is.Not.EqualTo(0.0));
        var local = SequenceAligner.LocalAlignAffine(a, b, ms);
        Assert.That(local, Is.EqualTo(SequenceAligner.LocalAlignAffineLinearSpace(a, b, ms, 1_000_000)));
    }

    [Test]
    public void LinearSpace_Memory_IsLinear_NotQuadratic()
    {
        var (a, b) = RelatedPair(5, 5600, 0.10);
        long cells = (long)(a.Length + 1) * (b.Length + 1);
        long fullBytes = cells * 12; // three int matrices (affine)
        SequenceAligner.GlobalAlignAffineLinearSpace(a.Substring(0, 300), b.Substring(0, 300), SequenceAligner.BlastDna, 1000); // JIT warm-up
        long before = GC.GetAllocatedBytesForCurrentThread();
        var r = SequenceAligner.GlobalAlignAffineLinearSpace(a, b, SequenceAligner.BlastDna, 20_000);
        long allocated = GC.GetAllocatedBytesForCurrentThread() - before;
        Assert.That(cells * 12L, Is.GreaterThan(350_000_000L));
        // O(n) rows per recursion node + O(m+n) output + small full-matrix leaves: a tiny fraction of the full matrices.
        Assert.That(allocated, Is.LessThan(fullBytes / 8), () => $"allocated {allocated:N0} bytes, full matrices {fullBytes:N0}");
        Assert.That(r.AlignedSequence1.Replace("-", ""), Is.EqualTo(a));
        Assert.That(r.AlignedSequence2.Replace("-", ""), Is.EqualTo(b));
    }

    [Test]
    public void AutomaticSwitch_Cancellation_IsHonoured()
    {
        string a = RandomSeq(new Random(1), "ACGT", 8000);
        string b = RandomSeq(new Random(2), "ACGT", 8000); // 64 M cells -> linear-space program
        using var cts = new CancellationTokenSource(TimeSpan.FromMilliseconds(30));
        Assert.Throws<OperationCanceledException>(() => SequenceAligner.GlobalAlign(a, b, SequenceAligner.SimpleDna, cts.Token));
    }

    // ---------------------------------------------------------------- hand-checkable values

    [Test]
    public void GlobalAffine_LongGap_KnownValue_ThroughRecursion()
    {
        // BLAST DNA (2/-3, gap 5 + 2k): the 3-nt gap costs 5 + 3*2 = 11; 8 matches = 16 -> score 5.
        var r = SequenceAligner.GlobalAlignAffineLinearSpace("ACGTTTTACGT", "ACGTACGT", SequenceAligner.BlastDna, 1);
        Assert.That(r.Score, Is.EqualTo(5));
        Assert.That(r.AlignedSequence1.Replace("-", ""), Is.EqualTo("ACGTTTTACGT"));
        Assert.That(r, Is.EqualTo(SequenceAligner.GlobalAlignAffine("ACGTTTTACGT", "ACGTACGT", SequenceAligner.BlastDna)));
    }

    [Test]
    public void LocalAffine_ReferenceValue_Biopython()
    {
        // F2: LocalAlignAffine(ACGTACGTAAAAAAACGTACGT, ACGTACGTACGTACGT, HighIdentityDna) = 64 (Biopython 1.88 / parasail 1.3.4).
        var r = SequenceAligner.LocalAlignAffineLinearSpace("ACGTACGTAAAAAAACGTACGT", "ACGTACGTACGTACGT", SequenceAligner.HighIdentityDna, 4);
        Assert.That(r.Score, Is.EqualTo(64));
        Assert.That(r, Is.EqualTo(SequenceAligner.LocalAlignAffine("ACGTACGTAAAAAAACGTACGT", "ACGTACGTACGTACGT", SequenceAligner.HighIdentityDna)));
    }

    [Test]
    public void EmptyInputs_AndTinySequences_FollowTheFullMatrixContract()
    {
        foreach (bool affine in new[] { false, true })
            foreach (bool local in new[] { false, true })
                foreach (var ends in local ? new[] { EndGaps.None } : AllPolicies)
                    foreach (var (a, b) in new[] { ("", ""), ("", "ACGT"), ("ACGT", ""), ("A", "A"), ("A", "C"), ("AC", "A") })
                        AssertSame(Full(a, b, SequenceAligner.SimpleDna, affine, local, ends), Fast(a, b, SequenceAligner.SimpleDna, affine, local, ends, 1),
                            $"affine={affine} local={local} ends={ends} '{a}' '{b}'");
    }

    [Test]
    public void Validation_NullInputs_And_PositiveLocalGapScores()
    {
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.GlobalAlignAffineLinearSpace((string)null!, "A"));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.LocalAlignLinearSpace("A", (string)null!));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.SemiGlobalAlignAffineLinearSpace((DnaSequence)null!, new DnaSequence("A")));
        // a local alignment never starts with a gap only for non-positive gap scores
        Assert.Throws<ArgumentException>(() => SequenceAligner.LocalAlignLinearSpace("ACGT", "ACGT", new ScoringMatrix(1, -1, 0, 1)));
        Assert.Throws<ArgumentException>(() => SequenceAligner.LocalAlignAffineLinearSpace("ACGT", "ACGT", new ScoringMatrix(1, -1, 1, -1)));
        // matrix alphabet is still validated
        Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlignLinearSpace("AC1", "ACG", new MatrixScoring(SubstitutionMatrix.EdnaFull, -2, -1)));
    }

    [Test]
    public void LowerCaseInput_IsUpperCased_LikeEveryAligner()
    {
        Assert.That(SequenceAligner.GlobalAlignAffineLinearSpace("acgtacgt", "ACGTTACGT", SequenceAligner.BlastDna, 2),
            Is.EqualTo(SequenceAligner.GlobalAlignAffine("ACGTACGT", "ACGTTACGT", SequenceAligner.BlastDna)));
    }

    // ---------------------------------------------------------------- oracle: Biopython / parasail / EMBOSS stretcher / matcher

    private sealed record Row(
        string Kind, string Family, string Model, int Flags, string Matrix, int Match, int Mismatch, double GapOpen, double GapExtend,
        string Seq1, string Seq2, double Score, double? Parasail, string Extra);

    private static readonly Lazy<List<Row>> Oracle = new(ReadOracle);

    private static double D(string s) => double.Parse(s, CultureInfo.InvariantCulture);

    private static List<Row> ReadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.linear_space_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<Row>();
        reader.ReadLine();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0 || line[0] == '#') continue;
            var f = line.Split('\t');
            rows.Add(new Row(f[0], f[1], f[2], int.Parse(f[3], CultureInfo.InvariantCulture), f[4],
                int.Parse(f[5], CultureInfo.InvariantCulture), int.Parse(f[6], CultureInfo.InvariantCulture), D(f[7]), D(f[8]),
                f[9], f[10], D(f[11]), f[12] == "NA" ? null : D(f[12]), f[13]));
        }
        return rows;
    }

    /// <summary>Independent re-scoring of rows: aligned columns by <paramref name="sub"/>, a gap run of length k costs
    /// k·extend (+ open when affine); a gap run that touches a free end under <paramref name="ends"/> costs nothing.</summary>
    private static double Rescore(string r1, string r2, EndGaps ends, Func<char, char, double> sub, double open, double extend, bool affine)
    {
        double total = 0;
        int n = r1.Length;
        int c = 0;
        while (c < n)
        {
            if (r1[c] != '-' && r2[c] != '-')
            {
                total += sub(r1[c], r2[c]);
                c++;
                continue;
            }
            bool gapInSeq2 = r2[c] == '-';
            int start = c;
            while (c < n && (gapInSeq2 ? r2[c] == '-' && r1[c] != '-' : r1[c] == '-' && r2[c] != '-')) c++;
            bool atStart = start == 0, atEnd = c == n;
            bool free = gapInSeq2
                ? (atStart && ends.HasFlag(EndGaps.FreeStart1)) || (atEnd && ends.HasFlag(EndGaps.FreeEnd1))
                : (atStart && ends.HasFlag(EndGaps.FreeStart2)) || (atEnd && ends.HasFlag(EndGaps.FreeEnd2));
            if (!free)
                total += (c - start) * extend + (affine ? open : 0);
        }
        return total;
    }

    private static string Canonical(string r1, string r2)
    {
        var o1 = new StringBuilder();
        var o2 = new StringBuilder();
        int c = 0;
        while (c < r1.Length)
        {
            if (r1[c] != '-' && r2[c] != '-')
            {
                o1.Append(r1[c]); o2.Append(r2[c]); c++;
                continue;
            }
            int s = c;
            while (c < r1.Length && (r1[c] == '-' || r2[c] == '-') && !(r1[c] == '-' && r2[c] == '-')) c++;
            for (int k = s; k < c; k++) if (r2[k] == '-') { o1.Append(r1[k]); o2.Append('-'); }
            for (int k = s; k < c; k++) if (r1[k] == '-') { o1.Append('-'); o2.Append(r2[k]); }
        }
        return o1 + "|" + o2;
    }

    [Test]
    public void Oracle_Biopython_Parasail_ScoresAndCoOptimalRows_ThroughRecursion()
    {
        var rows = Oracle.Value.Where(r => r.Kind == "bio").ToList();
        Assert.That(rows, Has.Count.GreaterThanOrEqualTo(1200));
        int membership = 0, parasailChecked = 0, longGapRows = 0;
        foreach (var r in rows)
        {
            bool affine = r.Model == "affine";
            var sc = new ScoringMatrix(r.Match, r.Mismatch, affine ? (int)r.GapOpen : 0, (int)r.GapExtend);
            var ends = (EndGaps)r.Flags;
            bool local = r.Family == "local";
            long baseCells = 1 + (r.Seq1.Length * 7 + r.Seq2.Length) % 60; // deep, varied recursion
            var res = Fast(r.Seq1, r.Seq2, sc, affine, local, r.Family == "global" ? EndGaps.None : ends, baseCells);
            string ctx = $"{r.Family}/{r.Model}/flags={r.Flags} sc=({r.Match},{r.Mismatch},{r.GapOpen},{r.GapExtend}) {r.Seq1} {r.Seq2}";
            Assert.That(res.Score, Is.EqualTo(r.Score).Within(1e-9), ctx);
            if (r.Parasail is double ps)
            {
                Assert.That(res.Score, Is.EqualTo(ps).Within(1e-9), "parasail " + ctx);
                parasailChecked++;
            }

            // the returned alignment re-scores to its score and is one of Biopython's co-optimal alignments
            double re = Rescore(res.AlignedSequence1, res.AlignedSequence2, local ? EndGaps.None : r.Family == "global" ? EndGaps.None : ends,
                (x, y) => x == y ? r.Match : r.Mismatch, r.GapOpen, r.GapExtend, affine);
            Assert.That(re, Is.EqualTo(r.Score).Within(1e-9), "re-score " + ctx);
            Assert.That(res.AlignedSequence1.Replace("-", ""), Is.EqualTo(local ? r.Seq1.Substring(res.StartPosition1, res.EndPosition1 - res.StartPosition1 + 1) : r.Seq1), ctx);
            if (r.Extra != "NA")
            {
                var set = r.Extra.Split(';').Select(p => { var q = p.Split('|'); return Canonical(q[0], q[1]); }).ToHashSet();
                Assert.That(set, Does.Contain(Canonical(res.AlignedSequence1, res.AlignedSequence2)), "co-optimal " + ctx);
                membership++;
            }
            // the linear-space alignment is the full-matrix alignment
            AssertSame(Full(r.Seq1, r.Seq2, sc, affine, local, r.Family == "global" ? EndGaps.None : ends), res, ctx);
            if (Math.Abs(r.Seq1.Length - r.Seq2.Length) >= 5) longGapRows++;
        }
        Assert.Multiple(() =>
        {
            Assert.That(membership, Is.GreaterThanOrEqualTo(400));
            Assert.That(parasailChecked, Is.GreaterThanOrEqualTo(500));
            Assert.That(longGapRows, Is.GreaterThanOrEqualTo(300));
        });
    }

    [Test]
    public void Oracle_EmbossStretcher_MyersMiller_GlobalScoreAndRows()
    {
        var rows = Oracle.Value.Where(r => r.Kind == "stretcher").ToList();
        Assert.That(rows, Has.Count.EqualTo(160));
        int uniqueRows = 0;
        foreach (var r in rows)
        {
            var matrix = r.Matrix == "BLOSUM62" ? SubstitutionMatrix.Blosum62 : SubstitutionMatrix.EdnaFull;
            var ms = new MatrixScoring(matrix, r.GapOpen, r.GapExtend);
            var res = SequenceAligner.GlobalAlignAffineLinearSpace(r.Seq1, r.Seq2, ms, 1 + (r.Seq1.Length + r.Seq2.Length) % 50);
            string ctx = $"stretcher {r.Matrix} {r.GapOpen}/{r.GapExtend} {r.Seq1} {r.Seq2}";
            Assert.That(res.ExactScore, Is.EqualTo(r.Score).Within(1e-6), ctx);
            Assert.That(res, Is.EqualTo(SequenceAligner.GlobalAlignAffine(r.Seq1, r.Seq2, ms)), ctx);
            double re = Rescore(res.AlignedSequence1, res.AlignedSequence2, EndGaps.None, (x, y) => matrix.GetScore(x, y), r.GapOpen, r.GapExtend, affine: true);
            Assert.That(re, Is.EqualTo(r.Score).Within(1e-6), "re-score " + ctx);
            var q = r.Extra.Split('|');
            bool embossRowsAreOptimal = Math.Abs(Rescore(q[0], q[1], EndGaps.None, (x, y) => matrix.GetScore(x, y), r.GapOpen, r.GapExtend, affine: true) - r.Score) < 1e-6;
            Assert.That(embossRowsAreOptimal, Is.True, "EMBOSS rows re-score to the EMBOSS score " + ctx);
            if (q[2] == "1")
            {
                // Biopython finds exactly one optimal alignment: Myers-Miller, the full matrix and the linear-space program must all return it
                uniqueRows++;
                Assert.That(Canonical(res.AlignedSequence1, res.AlignedSequence2), Is.EqualTo(Canonical(q[0], q[1])), "unique optimum " + ctx);
            }
        }
        Assert.That(uniqueRows, Is.GreaterThanOrEqualTo(40));
    }

    [Test]
    public void Oracle_EmbossMatcher_HuangMillerSim_LocalScoreAndRows()
    {
        var rows = Oracle.Value.Where(r => r.Kind == "matcher").ToList();
        Assert.That(rows, Has.Count.EqualTo(160));
        int uniqueRows = 0;
        foreach (var r in rows)
        {
            var matrix = r.Matrix == "BLOSUM62" ? SubstitutionMatrix.Blosum62 : SubstitutionMatrix.EdnaFull;
            var ms = new MatrixScoring(matrix, r.GapOpen, r.GapExtend);
            var res = SequenceAligner.LocalAlignAffineLinearSpace(r.Seq1, r.Seq2, ms, 1 + (r.Seq1.Length + r.Seq2.Length) % 50);
            string ctx = $"matcher {r.Matrix} {r.GapOpen}/{r.GapExtend} {r.Seq1} {r.Seq2}";
            Assert.That(res.ExactScore, Is.EqualTo(r.Score).Within(1e-6), ctx);
            Assert.That(res, Is.EqualTo(SequenceAligner.LocalAlignAffine(r.Seq1, r.Seq2, ms)), ctx);
            double re = Rescore(res.AlignedSequence1, res.AlignedSequence2, EndGaps.None, (x, y) => matrix.GetScore(x, y), r.GapOpen, r.GapExtend, affine: true);
            Assert.That(re, Is.EqualTo(r.Score).Within(1e-6), "re-score " + ctx);
            var q = r.Extra.Split('|');
            bool embossRowsAreOptimal = Math.Abs(Rescore(q[0], q[1], EndGaps.None, (x, y) => matrix.GetScore(x, y), r.GapOpen, r.GapExtend, affine: true) - r.Score) < 1e-6;
            Assert.That(embossRowsAreOptimal, Is.True, "EMBOSS rows re-score to the EMBOSS score " + ctx);
            if (q[6] == "1")
            {
                // exactly one optimal local alignment (Biopython): rows and 1-based coordinates equal matcher's
                uniqueRows++;
                Assert.That(Canonical(res.AlignedSequence1, res.AlignedSequence2), Is.EqualTo(Canonical(q[0], q[1])), "unique optimum " + ctx);
                Assert.That((res.StartPosition1 + 1, res.EndPosition1 + 1, res.StartPosition2 + 1, res.EndPosition2 + 1),
                    Is.EqualTo((int.Parse(q[2]), int.Parse(q[3]), int.Parse(q[4]), int.Parse(q[5]))), ctx);
            }
        }
        Assert.That(uniqueRows, Is.GreaterThanOrEqualTo(80));
    }
}
