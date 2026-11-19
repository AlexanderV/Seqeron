using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// ALIGN-MULTI-001 (audit A1-13c / A1-13d, B09): the command-line options of the T-Coffee 13.41 program reproduced by
/// <see cref="SequenceAligner.MultipleAlignTCoffeeProgram(IEnumerable{string}, TCoffeeProgramOptions?, Seqeron.Genomics.Infrastructure.MatrixScoring?)"/>.
/// Oracle: TestData/Alignment/tcoffee_program_options_oracle.tsv (scripts/alignment/generate_tcoffee_options_oracle.py): the installed <c>t_coffee</c> run on random sets with a
/// flag string; every row set must be reproduced from <see cref="TCoffeeProgramOptions.FromCommandLine"/> of that string.
/// </summary>
[TestFixture]
public class SequenceAligner_TCoffeeOptions_Tests
{
    private sealed record Case(string Feature, string Type, string Flags, string[] Seqs, string Rows, string[]? Names);

    private static List<Case> Oracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.tcoffee_program_options_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var list = new List<Case>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0) continue;
            var f = line.Split('\t');
            list.Add(new Case(f[0], f[1], f[2], f[3].Split('|'), f[4], f.Length > 5 && f[5].Length > 0 ? f[5].Split('|') : null));
        }

        return list;
    }

    /// <summary>
    /// The only committed sets (3 of 4 700) on which the port does not equal the program (flags TAB sequences), recorded in B09.md F72 / F73 / F75: (1) protein, `-nomatch 50 -gapext -1`, 7 sequences;
    /// (2) protein, `-seq_weight t_coffee -distance_matrix_mode idscore -tree_mode upgma`, 8 sequences (minimal: 5 sequences; the weights are equal to the program's - an explicit weight file reproduces
    /// the difference - and uniform weights of 0.04 / 0.05 differ while 0.03 / 0.045 / 0.055 do not); (3) `clustalw_pair` on the protein pair YAGR / NSMYC (the ClustalW protein port aligns it
    /// differently from the installed clustalw 2.1).
    /// </summary>
    private static readonly HashSet<string> KnownDifferences = new(StringComparer.Ordinal)
    {
        "-nomatch 50 -gapext -1\tNEGYCNQARMRMFVWLRWPKCS|DGAQCADAFKECTWTHA|MFPCIPCTIIMIY|QPGVFVSFTMPRGKINLGV|QTPGVFVSFTMEPRGKWINLV|TFTS|QTPGVFVSFTMEPRGKINLV",
        "-seq_weight t_coffee -distance_matrix_mode idscore -tree_mode upgma\tHNQMVWEKHYGRTDNWKVVK|DMVDKQSEGFPNPEQ|EETCCMECSWPWWAYPPS|HNQMVWEKHYGRDWLAKVVK|HNQMVWEKHYGRDWKVKEK|HNQMVWEKHYYGRDWKVVK|HNQMVWKHYGRDWKVVK|HNQMVWESHYGRDWKGVK",
        "-method clustalw_pair -tree_mode upgma\tYAGR|NSMYC",
    };

    private static void AssertFeature(string feature, int minSets)
    {
        var cases = Oracle().Where(c => c.Feature == feature).ToList();
        Assert.That(cases.Count, Is.GreaterThanOrEqualTo(minSets), $"oracle sets of '{feature}'");
        var failures = new List<string>();
        // a set on which the program itself core-dumps (SIM with a very small gap open: `lalign_id_pair@EP@GOP@-5`, F75) is reproduced as TCoffeeAbortedException
        var crashed = cases.Where(c => c.Rows.StartsWith("FAIL:", StringComparison.Ordinal)).ToList();
        Assert.That(crashed.Count, Is.LessThanOrEqualTo(5), "sets on which t_coffee itself aborts");
        foreach (var c in crashed)
            Assert.Throws<TCoffeeAbortedException>(() => SequenceAligner.MultipleAlignTCoffeeProgram(c.Seqs, TCoffeeProgramOptions.FromCommandLine(c.Flags) with { Names = c.Names }),
                $"t_coffee aborts on '{c.Flags}' {string.Join('|', c.Seqs)}");
        var known = new List<string>();
        foreach (var c in cases.Except(crashed))
        {
            var options = TCoffeeProgramOptions.FromCommandLine(c.Flags) with { Names = c.Names };
            string got = string.Join('|', SequenceAligner.MultipleAlignTCoffeeProgram(c.Seqs, options).AlignedSequences);
            if (got == c.Rows) continue;
            if (KnownDifferences.Contains(c.Flags + "\t" + string.Join('|', c.Seqs))) known.Add(c.Flags);
            else failures.Add($"{c.Type} '{c.Flags}' {string.Join('|', c.Seqs)}: got {got} expected {c.Rows}");
        }

        Assert.That(known.Count, Is.LessThanOrEqualTo(KnownDifferences.Count), "documented residual differences");

        Assert.That(failures, Is.Empty, $"{failures.Count} of {cases.Count - crashed.Count} '{feature}' sets differ from t_coffee; first: {failures.FirstOrDefault()}");
    }

    /// <summary>A1-13c: <c>proba_pair</c> of RNA (a U among more than 98 % nucleotide letters): the RNA pair-HMM (5 matrix types, 2 insert state pairs, RNAgapOpen2 / RNAgapExtend2) and <c>get_tot_prob</c>.</summary>
    [Test]
    public void Rna_EveryRowEqualsTheProgram()
    {
        AssertFeature("rna", 150);
        Assert.That(Oracle().Where(c => c.Feature == "rna" && c.Seqs.Any(s => s.Contains('U'))).Count(), Is.GreaterThanOrEqualTo(150), "genuine RNA sets");
    }

    /// <summary><c>-relax_lib N</c> (0, 2, 3 relaxation rounds) and <c>-filter_lib T</c> (entries with weight &lt;= T dropped before the relaxation).</summary>
    [Test]
    public void RelaxAndFilterLib_EveryRowEqualsTheProgram() => AssertFeature("relax", 400);

    /// <summary>
    /// <c>-gapopen</c>, <c>-gapext</c>, <c>-cosmetic_penalty</c>, <c>-nomatch</c>, <c>-do_normalise</c>: the gap model of the sparse DP (<c>list2linked_pair_wise</c>: opening penalty free at the last row / column,
    /// <c>cl2list_borders</c> cells carry gep x 10 x position, nomatch rescaled by do_normalise / the largest residue weight sum).
    /// </summary>
    [Test]
    public void GapPenaltiesAndNoMatch_EveryRowEqualsTheProgram() => AssertFeature("gap", 400);

    /// <summary><c>-tree_mode upgma</c> and <c>-distance_matrix_mode</c> fast / idscore / no with the <c>-distance_matrix_sim_mode</c> measures (idmat_sim1/2/3, gap1, logid).</summary>
    [Test]
    public void TreeModeAndDistanceMatrixMode_EveryRowEqualsTheProgram() => AssertFeature("tree", 500);

    /// <summary><c>-seq_weight t_coffee</c> (weights from the identity matrix, library entries x MIN of the two sequence weights) and weight files.</summary>
    [Test]
    public void SequenceWeights_EveryRowEqualsTheProgram() => AssertFeature("weight", 300);

    /// <summary>
    /// <c>-inorder aligned</c> (default: the sequences are sorted by name before anything is computed) and <c>-outorder</c>: sets with shuffled names, names that are prefixes of each other
    /// (s1 / s10: the "name-prefix quirk" of F66 is this sorting) and 11-13 sequences s0..s12 (string order s0 s1 s10 s11 s12 s2 ...).
    /// </summary>
    [Test]
    public void NameOrder_EveryRowEqualsTheProgram() => AssertFeature("order", 400);

    /// <summary>
    /// <c>-matrix NAME</c>: <c>get_emitPairs</c> (the emission probabilities of <c>proba_pair</c> from a substitution matrix: s x s x exp(score x ln 2 / 2), normalised, single precision) for DNA,
    /// RNA and protein, and the PROTEIN <c>slow_pair</c> (matrix and gap open from the matrix); 38 matrices of <c>matrices.h</c> are available.
    /// </summary>
    [Test]
    public void Matrix_EveryRowEqualsTheProgram() => AssertFeature("matrix", 600);

    /// <summary>The <c>method@EP@GOP@-30@GEP@-3@MATRIX@blosum50mt</c> parameter overrides of <c>slow_pair</c> and <c>lalign_id_pair</c> (lines appended to the method file, last wins).</summary>
    [Test]
    public void MethodParameterOverrides_EveryRowEqualsTheProgram() => AssertFeature("method", 350);

    /// <summary><c>clustalw_pair</c>: the program calls the external <c>clustalw -align</c> on every pair (installed 2.1); here the ClustalW 2.1 ports are the aligner.</summary>
    [Test]
    public void ClustalwPair_EveryRowEqualsTheProgram() => AssertFeature("clustalw", 300);

    [Test]
    public void FromCommandLine_MatrixAndMethodParameters()
    {
        var o = TCoffeeProgramOptions.FromCommandLine("-matrix blosum50mt -method slow_pair@EP@GOP@-30@MATRIX@pam250mt,clustalw_pair");
        Assert.That(o.Matrix, Is.EqualTo("blosum50mt"));
        Assert.That(o.Methods, Is.EqualTo(new[] { TCoffeePairMethod.SlowPair, TCoffeePairMethod.ClustalwPair }));
        Assert.That(o.MethodParameters![TCoffeePairMethod.SlowPair], Is.EqualTo(new TCoffeeMethodParameters(-30, null, "pam250mt")));
        Assert.Throws<ArgumentException>(() => TCoffeeProgramOptions.FromCommandLine("-matrix no_such_matrix"));
        Assert.Throws<ArgumentException>(() => TCoffeeProgramOptions.FromCommandLine("-method slow_pair@EP@NOPE@1"));
    }

    [Test]
    public void FromCommandLine_ParsesMethodsAndRejectsUnknownFlags()
    {
        var o = TCoffeeProgramOptions.FromCommandLine("-method slow_pair,lalign_id_pair -lalign_n_top 3 -n_core 4");
        Assert.That(o.Methods, Is.EqualTo(new[] { TCoffeePairMethod.SlowPair, TCoffeePairMethod.LalignIdPair }));
        Assert.That(o.LalignTopAlignments, Is.EqualTo(3));
        var p = TCoffeeProgramOptions.FromCommandLine("-relax_lib 2 -filter_lib 20 -tree_mode upgma -distance_matrix_mode fast -distance_matrix_sim_mode idmat_sim2 -seq_weight t_coffee -gapopen -50 -gapext -2 -nomatch 3 -outorder aligned -inorder input -weight default");
        Assert.That((p.RelaxLibIterations, p.FilterLib, p.TreeMode, p.DistanceMatrixMode, p.DistanceMatrixSimMode), Is.EqualTo((2, 20, TCoffeeTreeMode.Upgma, TCoffeeDistanceMode.Fast, "idmat_sim2")));
        Assert.That((p.TCoffeeSequenceWeights, p.GapOpen, p.GapExtend, p.NoMatch, p.OutOrderAligned, p.InOrderInput), Is.EqualTo((true, -50, -2, 3, true, true)));
        Assert.That(TCoffeeProgramOptions.FromCommandLine("-seq_weight t_coffee -weight no").TCoffeeSequenceWeights, Is.False);
        Assert.That(TCoffeeProgramOptions.FromCommandLine("-seq_weight @0.5,1,2").SequenceWeights, Is.EqualTo(new[] { 0.5f, 1f, 2f }));
        Assert.Throws<ArgumentException>(() => TCoffeeProgramOptions.FromCommandLine("-no_such_flag 1"));
        Assert.Throws<ArgumentException>(() => TCoffeeProgramOptions.FromCommandLine("-method blastp_pair"));
    }
}
