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

    private static bool IsIdscoreFast(string flags) => flags.Contains("-distance_matrix_mode fast", StringComparison.Ordinal) && flags.Contains("idscore", StringComparison.Ordinal);

    private static bool HasSecondRowGap(Case c, TCoffeeProgramOptions options)
    {
        bool protein = options.Protein ?? SequenceAligner.TcIsProtein(c.Seqs);
        var order = options.InOrderInput || c.Names == null ? Enumerable.Range(0, c.Seqs.Length) : Enumerable.Range(0, c.Seqs.Length).OrderBy(i => c.Names[i], StringComparer.Ordinal);
        var seqs = order.Select(i => c.Seqs[i]).ToArray();
        for (int a = 0; a < seqs.Length; a++)
            for (int b = a + 1; b < seqs.Length; b++)
                if (SequenceAligner.TCoffeeSlowPairAlign(seqs[a], seqs[b], protein).RowB.Contains('-')) return true;
        return false;
    }

    private static void AssertFeature(string feature, int minSets)
    {
        var cases = Oracle().Where(c => c.Feature == feature).ToList();
        Assert.That(cases.Count, Is.GreaterThanOrEqualTo(minSets), $"oracle sets of '{feature}'");
        var failures = new List<string>();
        int heapDependent = 0;
        // a set on which the program itself core-dumps (SIM's edit script overflows its heap chunk, `lalign_id_pair@EP@GOP@-5`, F75) is reproduced as TCoffeeAbortedException:
        // the rule is deterministic (script length > the usable bytes of the chunk), so "the program aborts" and "the port throws" must coincide on EVERY committed set
        foreach (var c in cases)
        {
            var options = TCoffeeProgramOptions.FromCommandLine(c.Flags) with { Names = c.Names };
            if (c.Rows.StartsWith("FAIL:", StringComparison.Ordinal))
            {
                Assert.Throws<TCoffeeAbortedException>(() => SequenceAligner.MultipleAlignTCoffeeProgram(c.Seqs, options), $"t_coffee aborts on '{c.Flags}' {string.Join('|', c.Seqs)}");
                continue;
            }

            try
            {
                string got = string.Join('|', SequenceAligner.MultipleAlignTCoffeeProgram(c.Seqs, options).AlignedSequences);
                if (got == c.Rows) continue;
                // idscore_pairseq reads m[s1[i] - 'a'][s2[j] - 'a']: with a gap in the SECOND aligned row that is the heap 208 bytes before a matrix row - 0 in 23 of the 24 committed sets that have such a gap,
                // some other value in one (undefined behaviour of the program, not reproducible)
                if (IsIdscoreFast(c.Flags) && HasSecondRowGap(c, options)) heapDependent++;
                else failures.Add($"{c.Type} '{c.Flags}' {string.Join('|', c.Seqs)}: got {got} expected {c.Rows}");
            }
            catch (TCoffeeAbortedException e) when (e.HeapDependent)
            {
                // undefined behaviour of the program (a wild pointer read whose outcome depends on the heap): it may also run through; only such sets may differ this way
                heapDependent++;
            }
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {cases.Count} '{feature}' sets differ from t_coffee; first: {failures.FirstOrDefault()}");
        Assert.That(heapDependent, Is.LessThanOrEqualTo(6), "sets on which the heap-dependent undefined behaviour of the idscore sim mode gave another result (5 aborts that the program survived + 1 differing set of 540)");
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

    /// <summary>
    /// A2-7: <c>clustalw_pair</c> on sets whose sequences come from different alphabets: the external clustalw guesses the residue type of every pair file itself (the SECOND sequence, 85 % of A C G T U N = DNA),
    /// which differs from T-Coffee's type of the whole set (the pair YAGR / NSMYC is nucleotide letters for T-Coffee but protein for clustalw); in DNA mode a protein letter scores 0 against everything.
    /// </summary>
    [Test]
    public void ClustalwPair_MixedAlphabets_EveryRowEqualsTheProgram() => AssertFeature("clustalw_mixed", 150);

    /// <summary>
    /// A2-7: the pair list of the progressive stage sums single-precision terms in the order of the group's sequences, which <c>recode_tree</c> fills in ASCENDING index order (not in the order of the
    /// tree walk): on 6-10 sequences with a nomatch penalty a cell whose exact value is an integer (22.0) is 21.999998 or 22.000002 depending on that order and moves the traceback.
    /// </summary>
    [Test]
    public void NoMatchOnLargeSets_FloatSumOrder_EveryRowEqualsTheProgram() => AssertFeature("floatorder", 250);

    /// <summary>A2-3: <c>biphasic_pair</c> (<c>biphasic_pair_wise</c>: five states, two affine gap models, the program's traceback that starts from the predecessor of the best final state), DNA / protein / RNA.</summary>
    [Test]
    public void BiphasicPair_EveryRowEqualsTheProgram() => AssertFeature("biphasic", 250);

    /// <summary>A2-4: <c>lalign_rs_s_pair</c> (<c>ALN_MODE s_pairwise</c>): <c>lalign_id_pair</c> on every pair including each sequence against itself (SIM sees two sequences, the main diagonal is a candidate).</summary>
    [Test]
    public void LalignRsSPair_EveryRowEqualsTheProgram() => AssertFeature("lalign_rs_s", 250);

    /// <summary>
    /// A2-4: <c>promo_pair</c> / <c>blastr_pair</c> (<c>slow_pair</c> with <c>EXTEND_SEQ 1</c>, <c>WEIGHT extend_sim</c>, matrices promoter_tf1 / blosumR): the dinucleotide alphabet of <c>extend_seq</c>, the lossy
    /// <c>unextend_seq</c> on the program's shared sequence array (protein letters become j for good) and the weight of the extended rows.
    /// </summary>
    [Test]
    public void PromoAndBlastrPair_EveryRowEqualsTheProgram() => AssertFeature("promo_blastr", 250);

    /// <summary>A2-3 / A2-4: the methods the program refuses or dies on (clean_slow_pair, hash_pair, subop1_pair, subop2_pair, cfast_pair): every oracle set is a failure of the program and raises <see cref="TCoffeeAbortedException"/>.</summary>
    [Test]
    public void RefusedMethods_RaiseTCoffeeAbortedException_AsTheProgramFails()
    {
        var cases = Oracle().Where(c => c.Feature == "refused").ToList();
        Assert.That(cases.Count, Is.GreaterThanOrEqualTo(150));
        Assert.That(cases.All(c => c.Rows.StartsWith("FAIL:", StringComparison.Ordinal)), Is.True, "the program fails on every set");
        AssertFeature("refused", 150);
        // one sequence: no pair job and no refusal
        Assert.That(SequenceAligner.MultipleAlignTCoffeeProgram(new[] { "ACGT" }, TCoffeeProgramOptions.FromCommandLine("-method hash_pair")).AlignedSequences, Is.EqualTo(new[] { "ACGT" }));
        var e = Assert.Throws<TCoffeeAbortedException>(() => SequenceAligner.MultipleAlignTCoffeeProgram(new[] { "ACGT", "ACGA" }, TCoffeeProgramOptions.FromCommandLine("-method clean_slow_pair")))!;
        Assert.That((e.ExitStatus, e.Reason), Is.EqualTo((1, "ERROR: Weight Mode clean is unknown [FATAL:T-COFFEE] (seqpair2weight, util_constraints_list.c; COREDUMP in the pair job)")));
    }

    /// <summary>
    /// A2-6: <c>-distance_matrix_sim_mode</c> containing <c>idscore</c> with <c>-distance_matrix_mode fast</c>: <c>get_seq_sim</c> hands the ALIGNED rows to <c>idscore_pairseq</c>, which indexes the matrix with a gap character:
    /// a gap in the first row of a pair is a wild pointer (exit status 11), a gap in the second row reads 0. 260 sets with indels (the program dies in 177) and 280 gap-free sets (rows equal).
    /// </summary>
    [Test]
    public void DistanceSimModeIdscore_GapsInTheFirstRowSegfault_OtherwiseEveryRowEqualsTheProgram()
    {
        AssertFeature("simidscore", 150);
        AssertFeature("simidscore_gapless", 150);
        Assert.That(Oracle().Count(c => c.Feature.StartsWith("simidscore", StringComparison.Ordinal) && c.Rows.StartsWith("FAIL:", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(150), "sets on which the program dies");
        var e = Assert.Throws<TCoffeeAbortedException>(() => SequenceAligner.MultipleAlignTCoffeeProgram(new[] { "MKVLAAGIVGLLLAQP", "MKVLAAGIVGLLLAQPSTF" },
            TCoffeeProgramOptions.FromCommandLine("-distance_matrix_mode fast -distance_matrix_sim_mode idscore")))!;
        Assert.That((e.ExitStatus, e.HeapDependent), Is.EqualTo((11, false)));
        // a gap only in the SECOND row reads zero: no crash
        Assert.DoesNotThrow(() => SequenceAligner.MultipleAlignTCoffeeProgram(new[] { "MKVLAAGIVGLLLAQPSTF", "MKVLAAGIVGLLLAQP" }, TCoffeeProgramOptions.FromCommandLine("-distance_matrix_mode fast -distance_matrix_sim_mode idscore")));
    }

    /// <summary>
    /// A2-6: the residue groups of <c>make_group_aa</c> selected by the text before the first '_' of <c>-distance_matrix_sim_mode</c> (simple, mafft, clustalw, polarity, vasiliky, make_all, an unknown name = no group at all)
    /// under <c>-distance_matrix_mode fast</c>; the groups are built by the first call of the main process, which is the distance matrix (also with a library method: sets of T / U letters, A1-13d had assumed idmat groups there).
    /// </summary>
    [Test]
    public void DistanceSimModeGroups_EveryRowEqualsTheProgram() => AssertFeature("simgroups", 250);

    /// <summary>
    /// A2-6: the tree and distance modes the program does not run: <c>-tree_mode cwph</c> / an unknown name (core dump, exit status 11, for three or more sequences; two sequences use a dummy tree), <c>kmeans</c> (the
    /// program ends with status 0 and writes no alignment), <c>-distance_matrix_mode aln</c> (segmentation fault), <c>ktup2</c> and unknown names (exit status 1): every failure row raises <see cref="TCoffeeAbortedException"/>,
    /// the two-sequence rows of the tree modes equal the program.
    /// </summary>
    [Test]
    public void TreeAndDistanceModesTheProgramRefuses_RaiseTCoffeeAbortedException()
    {
        AssertFeature("treemodes", 150);
        AssertFeature("distrefused", 150);
        var e = Assert.Throws<TCoffeeAbortedException>(() => SequenceAligner.MultipleAlignTCoffeeProgram(new[] { "ACGT", "ACGA", "ACGG" }, TCoffeeProgramOptions.FromCommandLine("-tree_mode kmeans")))!;
        Assert.That(e.ExitStatus, Is.Zero, "the program ends normally without an alignment");
        Assert.Throws<NotSupportedException>(() => SequenceAligner.MultipleAlignTCoffeeProgram(new[] { "ACGT", "ACGA", "ACGG" }, TCoffeeProgramOptions.FromCommandLine("-distance_matrix_mode random")));
    }

    /// <summary><c>extend_seq</c> / <c>unextend_seq</c> of reformat.c on hand-derived values: dinucleotide codes, n / x, u, gaps, the last residue, lossy protein letters.</summary>
    [Test]
    public void ExtendAndUnextendSeq_FollowReformatC()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceAligner.TcExtendSeq("ACGT"), Is.EqualTo("ELRT")); // ac -> e, cg -> l, gt -> r, the last residue keeps its letter
            Assert.That(SequenceAligner.TcExtendSeq("acgtu"), Is.EqualTo("elryu")); // lower case stays lower case; tu is tt (u -> t) -> y; the last u is not changed
            Assert.That(SequenceAligner.TcExtendSeq("AN"), Is.EqualTo("AN")); // a second residue n keeps the first letter
            Assert.That(SequenceAligner.TcExtendSeq("NA"), Is.EqualTo("XA")); // a first residue n gives x
            Assert.That(SequenceAligner.TcExtendSeq("A-C-G"), Is.EqualTo("E-L-G")); // the neighbour is the next UNGAPPED residue
            Assert.That(SequenceAligner.TcExtendSeq("MKV"), Is.EqualTo("JJV")); // letters outside a c g t u n x give j
            Assert.That(SequenceAligner.TcUnextendSeq("ELRT"), Is.EqualTo("ACGT"));
            Assert.That(SequenceAligner.TcUnextendSeq("XA"), Is.EqualTo("NA")); // x -> n
            Assert.That(SequenceAligner.TcUnextendSeq("JJV"), Is.EqualTo("JJV")); // lossy: j stays j
        });
    }

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
