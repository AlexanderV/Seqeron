using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// ALIGN-MULTI-001 (audit A1-13b, B09 F66): <see cref="SequenceAligner.MultipleAlignTCoffeeProgram(IEnumerable{string}, TCoffeeProgramOptions?, MatrixScoring?)"/> =
/// the T-Coffee 13.41 program's OWN pipeline: primary library (<c>proba_pair</c> = the default of <c>t_coffee seqs.fa</c>, <c>slow_pair</c> = Myers-Miller global
/// alignment, <c>lalign_id_pair</c> = SIM top-10 local alignments, combinations merged the way the program merges them), relaxation, k-tuple neighbour-joining
/// guide tree, sparse linked-pair-wise progressive stage.
///
/// Oracles (TestData/Alignment, scripts/alignment/generate_tcoffee_program_oracle.py): the installed <c>t_coffee</c> (Debian package t-coffee
/// 13.45.0.4846264+really13.41.0.28bdc39+dfsg-1) run on random DNA and PROTEIN sets -
/// the primary library <c>-lib_only -out_lib</c> writes (every pair with weight and count), the guide tree (parsed from the program's "Group n" log) and the rows of
/// <c>-output fasta_aln</c>. The library, tree and every row of the program must be reproduced exactly.
/// </summary>
[TestFixture]
public class SequenceAligner_TCoffeeProgram_Tests
{
    private sealed record LibraryCase(string Type, string Methods, string[] Seqs, string Library);

    private sealed record MsaCase(string Type, string Methods, string[] Seqs, string Tree, string[] Rows);

    private static IEnumerable<string[]> ReadTsv(string file)
    {
        string name = "Seqeron.Genomics.Tests.TestData.Alignment." + file;
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        reader.ReadLine(); // header
        string? line;
        while ((line = reader.ReadLine()) is not null)
            if (line.Length > 0) yield return line.Split('\t');
    }

    private static List<LibraryCase> LibraryOracle() =>
        ReadTsv("tcoffee_program_library_oracle.tsv").Select(f => new LibraryCase(f[0], f[1], f[2].Split('|'), f[3])).ToList();

    private static List<MsaCase> MsaOracle(string file = "tcoffee_program_msa_oracle.tsv") =>
        ReadTsv(file).Select(f => new MsaCase(f[0], f[1], f[2].Split('|'), f[3], f[4].Split('|'))).ToList();

    private static TCoffeePairMethod[] Methods(string text) =>
        (text == "plain" ? "proba_pair" : text).Split(',').Select(m => m switch
        {
            "slow_pair" => TCoffeePairMethod.SlowPair,
            "lalign_id_pair" => TCoffeePairMethod.LalignIdPair,
            "proba_pair" => TCoffeePairMethod.ProbaPair,
            _ => throw new ArgumentException(m),
        }).ToArray();

    private static string TreeText(TCoffeeTreeNode n) => n.Left is null ? n.LeafIndex.ToString(CultureInfo.InvariantCulture) : "(" + TreeText(n.Left) + "," + TreeText(n.Right!) + ")";

    // ----- the oracle files -----------------------------------------------------------------------------------------------

    [Test]
    public void OracleFiles_AreComplete()
    {
        var lib = LibraryOracle();
        Assert.That(lib.Count, Is.GreaterThanOrEqualTo(800));
        foreach (var m in new[] { "slow_pair", "lalign_id_pair", "proba_pair", "slow_pair,lalign_id_pair" })
            foreach (var t in new[] { "dna", "protein" })
                Assert.That(lib.Count(c => c.Type == t && c.Methods == m), Is.GreaterThanOrEqualTo(60), $"{t} {m}");

        var msa = MsaOracle();
        foreach (var (t, m) in new[] { ("dna", "plain"), ("protein", "plain"), ("dna", "slow_pair,lalign_id_pair"), ("protein", "slow_pair,lalign_id_pair") })
            Assert.That(msa.Count(c => c.Type == t && c.Methods == m), Is.GreaterThanOrEqualTo(200), $"{t} {m}");
        Assert.That(msa.Select(c => c.Seqs.Length).Distinct().Order(), Is.EqualTo(Enumerable.Range(2, 7)));
        Assert.That(msa.Where(c => c.Type == "protein").SelectMany(c => c.Seqs).SelectMany(s => s).Distinct().Count(), Is.GreaterThan(15), "a real protein alphabet");
    }

    // ----- stage a, b, d: the primary library ---------------------------------------------------------------------------------

    /// <summary>
    /// <c>t_coffee -seq s.fa -method M -lib_only -out_lib x.tc_lib</c> = <see cref="SequenceAligner.TCoffeePrimaryLibraryEntries"/>: slow_pair (Myers-Miller with dna_idmat 10 / 0, gap -10 -1, or
    /// BLOSUM62 with (int)(10 x mean mismatch) -1), lalign_id_pair (SIM top 10, DNA 10 / 0 and -10 -1, protein BLOSUM50 -10 -4), proba_pair (ProbCons pair-HMM posteriors), their
    /// merge (the program's forked merge loses and duplicates entries), weights 10 x percent identity and the pair counts.
    /// </summary>
    [Test]
    public void PrimaryLibrary_EveryPairEqualsTheProgram()
    {
        var failures = new List<string>();
        var cases = LibraryOracle();
        foreach (var c in cases)
        {
            string got = string.Join(';', SequenceAligner.TCoffeePrimaryLibraryEntries(c.Seqs, Methods(c.Methods))
                .Select(e => $"{e.SeqA}:{e.PosA + 1}:{e.SeqB}:{e.PosB + 1}:{e.Weight}:{e.Count}"));
            if (got != c.Library) failures.Add($"{c.Type} {c.Methods} {string.Join('|', c.Seqs)}");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {cases.Count} libraries differ from t_coffee; first: {failures.FirstOrDefault()}");
    }

    /// <summary>
    /// The merge of two methods is NOT the union of their libraries: the second method's child returns only the list tail beyond the recorded lengths
    /// (<c>freeze_constraint_list</c>), so an entry inserted before an older one of the same residue is lost and older entries are returned again. Values written by
    /// <c>t_coffee -seq s.fa -method slow_pair,lalign_id_pair -lalign_n_top 5 -lib_only -out_lib</c> for GGATCACAG / CTACACTG (the union has 14 pairs).
    /// </summary>
    [Test]
    public void PrimaryLibrary_TwoMethods_IsNotTheUnion_WorkedExample()
    {
        string[] seqs = { "GGATCACAG", "CTACACTG" };
        string Lib(params TCoffeePairMethod[] m) =>
            string.Join(' ', SequenceAligner.TCoffeePrimaryLibraryEntries(seqs, m, 5).Select(e => $"({e.PosA + 1},{e.PosB + 1},{e.Weight},{e.Count})"));

        Assert.That(Lib(TCoffeePairMethod.SlowPair, TCoffeePairMethod.LalignIdPair),
            Is.EqualTo("(2,1,500,2) (3,2,500,2) (3,3,1000,2) (4,3,500,3) (5,4,500,3) (6,5,800,4) (7,6,800,4) (8,7,800,4) (9,8,500,2)"));
        Assert.That(Lib(TCoffeePairMethod.LalignIdPair, TCoffeePairMethod.SlowPair),
            Is.EqualTo("(2,1,500,1) (3,3,1000,6) (4,2,1000,2) (4,3,500,1) (5,1,800,2) (5,4,500,1) (6,3,1000,2) (6,5,800,2) (7,4,1000,2) (7,6,800,2) (8,5,1000,2) (8,7,800,2) (9,8,800,2)"));
        Assert.That(Lib(TCoffeePairMethod.LalignIdPair),
            Is.EqualTo("(3,3,1000,2) (4,2,1000,1) (5,1,800,1) (6,3,1000,1) (6,5,800,1) (7,4,1000,1) (7,6,800,1) (8,5,1000,1) (8,7,800,1) (9,8,800,1)"));
    }

    /// <summary>
    /// The Myers-Miller routine of the program treats a leading gap in the SECOND sequence as an ordinary (opened) gap when the first sequence has one residue left in a
    /// sub-problem, so for TTTCAACTC against TTTTCAACTC the optimum under its own model (a leading gap, -10) is not what it returns: the library of
    /// <c>t_coffee -seq s.fa -method slow_pair -lib_only</c> holds the pairs 1-1, 2-3, 3-4, ... 9-10, i.e. the alignment <c>T-TTCAACTC</c> / <c>TTTTCAACTC</c>. The same
    /// routine (gap open -4) maps the segments of SIM alignments onto the full sequences (<c>fix_seq_seq</c>), which is why the top-10 library holds pairs the
    /// local alignments never contained.
    /// </summary>
    [Test]
    public void SlowPair_MyersMiller_PlacesTheGapInside_WorkedExample()
    {
        // slow_pair = the same routine with gap open -10 / extend -1: the optimum under the model is a leading gap, the routine's answer an internal one.
        var (a, b) = SequenceAligner.TCoffeeSlowPairAlign("TTTCAACTC", "TTTTCAACTC");
        Assert.That(a, Is.EqualTo("T-TTCAACTC"));
        Assert.That(b, Is.EqualTo("TTTTCAACTC"));
    }

    [Test]
    public void LalignTopAlignments_AreNonIntersecting_AndInScoreOrder()
    {
        var alns = SequenceAligner.TCoffeeLalignTopAlignments("TTTTCGAACTTC", "TTTTCAACTC", 4);
        Assert.That(alns.Select(x => x.Score), Is.EqualTo(new[] { 790, 500, 400, 370 }));
        Assert.That(alns[0].RowA, Is.EqualTo("TTTTCGAACT"));
        Assert.That(alns[0].RowB, Is.EqualTo("TTTTC-AACT"));
        var used = new HashSet<(int, int)>();
        foreach (var x in alns)
        {
            int i = x.StartA, j = x.StartB;
            for (int c = 0; c < x.RowA.Length; c++)
            {
                if (x.RowA[c] != '-' && x.RowB[c] != '-') Assert.That(used.Add((i, j)), Is.True, "an aligned pair is reported once");
                if (x.RowA[c] != '-') i++;
                if (x.RowB[c] != '-') j++;
            }
        }
    }

    // ----- stage c: the guide tree and the alignment ---------------------------------------------------------------------------

    /// <summary>The guide tree (k-tuple distances, NJ, trifurcating root turned into a binary tree) equals the program's, on every oracle set of 3 or more sequences.</summary>
    [Test]
    public void GuideTree_EqualsTheProgram_OnEveryOracleSet()
    {
        var failures = new List<string>();
        var cases = MsaOracle().Where(c => c.Seqs.Length >= 3).ToList();
        foreach (var c in cases)
        {
            string got = TreeText(SequenceAligner.TCoffeeGuideTree(c.Seqs));
            if (got != c.Tree) failures.Add($"{c.Type} {string.Join('|', c.Seqs)}: got {got} expected {c.Tree}");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {cases.Count} trees differ from t_coffee; first: {failures.FirstOrDefault()}");
    }

    /// <summary>
    /// Every row of <c>t_coffee -seq s.fa [-method slow_pair,lalign_id_pair|...] -output fasta_aln</c>: DNA and protein, the plain default (proba_pair), slow_pair + lalign_id_pair
    /// and the other combinations, 2-8 sequences.
    /// </summary>
    [Test]
    public void MultipleAlignTCoffeeProgram_EveryRowEqualsTheProgram()
    {
        var failures = new List<string>();
        var cases = MsaOracle();
        foreach (var c in cases)
        {
            var got = SequenceAligner.MultipleAlignTCoffeeProgram(c.Seqs, new TCoffeeProgramOptions { Methods = Methods(c.Methods) }).AlignedSequences;
            if (!got.SequenceEqual(c.Rows))
                failures.Add($"{c.Type} {c.Methods} {string.Join('|', c.Seqs)}: got {string.Join('|', got)} expected {string.Join('|', c.Rows)}");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {cases.Count} sets differ from t_coffee; first: {failures.FirstOrDefault()}");
    }

    /// <summary>From 100 sequences on the program uses its fast neighbour joining (<c>fast_nj_tree</c>).</summary>
    [Test]
    public void MultipleAlignTCoffeeProgram_HundredSequences_FastNeighborJoining_EqualsTheProgram()
    {
        var cases = MsaOracle("tcoffee_program_big_oracle.tsv");
        Assert.That(cases.Select(c => c.Seqs.Length), Is.EqualTo(new[] { 100, 101 }));
        foreach (var c in cases)
        {
            Assert.That(TreeText(SequenceAligner.TCoffeeGuideTree(c.Seqs)), Is.EqualTo(c.Tree));
            var got = SequenceAligner.MultipleAlignTCoffeeProgram(c.Seqs, new TCoffeeProgramOptions { Methods = Methods(c.Methods) }).AlignedSequences;
            Assert.That(got, Is.EqualTo(c.Rows));
        }
    }

    // ----- API contract ----------------------------------------------------------------------------------------------------------

    [Test]
    public void MultipleAlignTCoffeeProgram_RowsSpellTheInputs_AndHaveEqualLength()
    {
        string[] seqs = { "ACGTTGCATGCCGATAGC", "ACGTGCATGCGATAGGC", "ACTTGCATCCGATTGC", "CGTTGCATGGCCGATAGC" };
        foreach (var methods in new[] { "plain", "slow_pair,lalign_id_pair" })
        {
            var r = SequenceAligner.MultipleAlignTCoffeeProgram(seqs, new TCoffeeProgramOptions { Methods = Methods(methods) });
            Assert.That(r.AlignedSequences.Select(x => x.Length).Distinct().Count(), Is.EqualTo(1));
            Assert.That(r.AlignedSequences.Select(x => x.Replace("-", "")), Is.EqualTo(seqs));
            Assert.That(r.Consensus.Length, Is.EqualTo(r.AlignedSequences[0].Length));
        }
    }

    [Test]
    public void MultipleAlignTCoffeeProgram_Example_EqualsTheProgram()
    {
        // t_coffee -seq s.fa -method slow_pair,lalign_id_pair -output fasta_aln  (program output, 18-nt example of the audit)
        string[] seqs = { "ACGTTGCATGCCGATAGC", "ACGTGCATGCGATAGGC", "ACTTGCATCCGATTGC", "CGTTGCATGGCCGATAGC" };
        var r = SequenceAligner.MultipleAlignTCoffeeProgram(seqs, new TCoffeeProgramOptions { Methods = new[] { TCoffeePairMethod.SlowPair, TCoffeePairMethod.LalignIdPair } });
        Assert.That(r.AlignedSequences, Is.EqualTo(new[] { "ACGTTGCAT-GCCGATAG-C", "ACGT-GCAT-GC-GATAGGC", "AC-TTGCAT--CCGATTG-C", "-CGTTGCATGGCCGATAG-C" }));
    }

    [Test]
    public void MultipleAlignTCoffeeProgram_DegenerateInputs()
    {
        Assert.That(SequenceAligner.MultipleAlignTCoffeeProgram(Array.Empty<string>()).AlignedSequences, Is.Empty);
        Assert.That(SequenceAligner.MultipleAlignTCoffeeProgram(new[] { "ACGT" }).AlignedSequences, Is.EqualTo(new[] { "ACGT" }));
        var two = SequenceAligner.MultipleAlignTCoffeeProgram(new[] { "ACGT", "ACGT" }).AlignedSequences;
        Assert.That(two, Is.EqualTo(new[] { "ACGT", "ACGT" }));
    }

    [Test]
    public void MultipleAlignTCoffeeProgram_Contract_Errors()
    {
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignTCoffeeProgram(null!));
        Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignTCoffeeProgram(new[] { "ACGT", "" }));
        Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignTCoffeeProgram(new[] { "ACGT", "AC-GT" }));
        Assert.Throws<ArgumentException>(() => SequenceAligner.TCoffeePrimaryLibraryEntries(new[] { "ACGT", "" }, new[] { TCoffeePairMethod.SlowPair }));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.TCoffeeLalignTopAlignments("ACGT", "ACGT", 0));
        // RNA (a U) is supported since A1-13c: the RNA pair-HMM
        Assert.That(SequenceAligner.MultipleAlignTCoffeeProgram(new[] { "ACGUACGU", "ACGUUCGU" }).AlignedSequences, Has.Length.EqualTo(2));
    }

    [Test]
    public void TCoffeeKtupSimilarity_IdenticalSequencesAre100_DisjointTuplesAre70()
    {
        // pscore = (int)((s0 - common) / s0 * 3 * 10 + 0.5): 0 for identical sequences, 30 without a common 6-tuple; similarity = 100 - min(pscore(i, j), pscore(j, i))
        var s = SequenceAligner.TCoffeeKtupSimilarity(new[] { "AAAAAAA", "AAAAAAA", "CCCCCCC" });
        Assert.That(s[0, 1], Is.EqualTo(100));
        Assert.That(s[0, 2], Is.EqualTo(70));
        Assert.That(s[2, 0], Is.EqualTo(70));
        Assert.That(s[1, 1], Is.EqualTo(100));
    }
}
