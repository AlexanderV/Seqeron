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

    /// <summary>The parsed oracle, ONE copy per fixture (audit A3-15): the method used to re-read and re-parse the whole embedded TSV (6.0 MB / 3.2 MB as UTF-16 strings, 20 321 / several thousand rows) on every call, i.e. in every test.</summary>
    private static List<Case> Oracle() => OracleRows.Value;

    private static readonly Lazy<List<Case>> OracleRows = new(ReadOracleRows);

    private static List<Case> ReadOracleRows()
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

    private static bool IsGotohDp(string flags) => System.Text.RegularExpressions.Regex.IsMatch(flags, "-dp_mode (gotoh_pair_wise|gotoh_pair_wise_test|gotoh_pair_wise_lgp|ssec_pair_wise|biphasic_pair_wise|glocal_pair_wise|viterbi_pair_wise)( |$)");

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

    /// <summary>
    /// <c>test_pair_wise</c>: (the port aborts, the rows overflow the row buffers). The program pads every row to l0 + l1 - n columns with n = l0 / 2 (l0 &lt; 5) or 5 and reads <c>generate_null (l1 - n)</c>: a second group with fewer than n columns
    /// is a deterministic segmentation fault (the port raises <see cref="TCoffeeAbortedException"/>; 98 of 98 committed sets that have such a merge die). Otherwise the rows are written into buffers of max_len + 1 bytes (<c>seq2aln</c>): up to
    /// one column more than the longest sequence the program runs through and writes exactly the port's rows (50 of 50 committed sets); two or more columns beyond them overflow the heap (valgrind: "Invalid write of size 1 ... in strcpy" of the
    /// <c>sprintf (seq, "%s", buf)</c>): 231 of 272 such sets die, 41 differ (32 lose the last column of every row, 9 have another row clobbered). Returns the largest overflow (columns written beyond the buffer) of the merges.
    /// </summary>
    private static (bool Aborted, bool Overflow, int Excess) TestPairDp(Case c, TCoffeeProgramOptions options)
    {
        SequenceAligner.TcMergeTrace = new();
        try
        {
            SequenceAligner.MultipleAlignTCoffeeProgram(c.Seqs, options);
            int excess = SequenceAligner.TcMergeTrace.Count == 0 ? 0 : SequenceAligner.TcMergeTrace.Max(m => m.Written - m.MaxLen);
            return (false, SequenceAligner.TcMergeTrace.Any(m => m.Written - m.MaxLen >= 2), excess);
        }
        catch (TCoffeeAbortedException)
        {
            return (true, false, 0);
        }
        finally
        {
            SequenceAligner.TcMergeTrace = null;
        }
    }

    private static bool TestPairOverflow(Case c, TCoffeeProgramOptions options, out int excess)
    {
        var (_, overflow, e) = TestPairDp(c, options);
        excess = e;
        return overflow;
    }

    /// <summary>
    /// The rows of <c>tcoffee_program_options_oracle.tsv</c> on which the PROGRAM is in undefined behaviour, enumerated (audit A3-14a, B09 F143): <c>tcoffee_options_ub_cases.tsv</c> = feature, 1-based data row of the oracle, category,
    /// detail. Every other row is compared exactly (rows equal, or the abort is raised). Causes, each proven on the installed t_coffee 13.41 (Debian 13.45.0.4846264+really13.41.0.28bdc39+dfsg-1) with valgrind (Evidence F113 / F115):
    /// <c>gotoh_row_overflow_*</c> (row writer of gotoh / biphasic / glocal_pair_wise writes the terminator one byte past a heap block when the alignment has 8 (mod 16) columns: abort, or the last column lost), <c>test_pair_*</c> (null
    /// allocation / <c>sprintf</c> overflow of <c>seq2aln</c>), <c>sw_uninitialised_write_*</c> (gotoh_pair_wise_sw: "Invalid write of size 4 ... 4 bytes after a block of size 124"), <c>relative_triplet_ub</c> (a read past the end of a library
    /// list, SIGFPE / wild read), <c>idscore_heap_*</c> (idscore_pairseq reads <c>m[s1[i] - 'a'][s2[j] - 'a']</c> with a gap in the second row: the heap 208 bytes before a matrix row). The port reproduces the common outcome (or refuses); the
    /// detail column records WHAT it does, so that a changed behaviour of the port fails the test. Regenerate with the environment variable SEQERON_TCOFFEE_UB_CASES_OUT=&lt;file&gt; (the observed lines are appended there).
    /// </summary>
    private static readonly Lazy<Dictionary<string, List<string>>> UbCases = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.tcoffee_options_ub_cases.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded data '{name}' not found.");
        using var reader = new StreamReader(stream);
        var map = new Dictionary<string, List<string>>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0 || line[0] == '#' || line.StartsWith("feature\t", StringComparison.Ordinal))
                continue;
            int tab = line.IndexOf('\t');
            string feature = line[..tab];
            if (!map.TryGetValue(feature, out var list))
                map[feature] = list = [];
            list.Add(line[(tab + 1)..]);
        }

        return map;
    });

    private static void AssertFeature(string feature, int minSets)
    {
        var cases = Oracle().Select((c, i) => (Case: c, Row: i + 1)).Where(x => x.Case.Feature == feature).ToList();
        Assert.That(cases.Count, Is.GreaterThanOrEqualTo(minSets), $"oracle sets of '{feature}'");
        var failures = new List<string>();
        var observed = new List<string>();
        foreach (var (c, row) in cases)
        {
            var options = TCoffeeProgramOptions.FromCommandLine(c.Flags) with { Names = c.Names };
            if (c.Rows.StartsWith("FAIL:", StringComparison.Ordinal))
            {
                if (c.Flags.Contains("relative_triplet", StringComparison.Ordinal))
                {
                    // the program divides by the total weight of an empty / short list (SIGFPE, exit status 8) or dies on the read past the list: the port stops with its UB refusal
                    string outcome = "returned";
                    try
                    {
                        SequenceAligner.MultipleAlignTCoffeeProgram(c.Seqs, options);
                    }
                    catch (NotSupportedException)
                    {
                        outcome = nameof(NotSupportedException);
                    }
                    catch (TCoffeeAbortedException)
                    {
                        outcome = nameof(TCoffeeAbortedException);
                    }

                    Assert.That(outcome, Is.Not.EqualTo("returned"), $"relative_triplet UB row {row}");
                    observed.Add($"{row}\trelative_triplet_ub\t{outcome}");
                    continue;
                }

                // gotoh_pair_wise(_test) writes the terminator of an alignment row one byte past a heap block when the alignment is 8 (mod 16) columns long (24, 40, 56: the committed aborts of the program, exit status 6, all have
                // such a length; a 1-byte overflow into the slack of a larger block is harmless, so only some sets of that length die): the port returns the alignment, the abort is checked against that necessary condition
                if (c.Rows is "FAIL:6" or "FAIL:-6" && IsGotohDp(c.Flags))
                {
                    string[] rows = SequenceAligner.MultipleAlignTCoffeeProgram(c.Seqs, options).AlignedSequences.ToArray();
                    int columns = rows[0].Length;
                    // the row writer overflows on the 8 (mod 16)-column MERGES of the progressive alignment (an intermediate merge can lie below the final length), so the necessary condition on the final alignment is only
                    // its minimum (8); the exact length of the port's alignment is recorded per row in tcoffee_options_ub_cases.tsv
                    Assert.That(columns, Is.GreaterThanOrEqualTo(8), $"heap-overflow abort of the row writer on '{c.Flags}' {string.Join('|', c.Seqs)}");
                    observed.Add($"{row}\tgotoh_row_overflow_abort\tcolumns={columns}");
                    continue;
                }

                if (c.Flags.Contains("test_pair_wise", StringComparison.Ordinal))
                {
                    var (aborted, overflow, excess) = TestPairDp(c, options);
                    Assert.That(aborted || overflow, Is.True, $"test_pair_wise dies without a null allocation or a row overflow: {string.Join('|', c.Seqs)}");
                    observed.Add(aborted ? $"{row}\ttest_pair_null_allocation_abort\tport aborts" : $"{row}\ttest_pair_row_overflow_abort\texcess={excess}");
                    continue;
                }

                // gotoh_pair_wise_sw: 49 of 500 committed sets die (exit status 6, -6, -9 or -11) although the port's rows equal the program's on the 450 others: undefined behaviour of the program (valgrind: "Invalid write of size 4 ... 4 bytes
                // after a block of size 124" inside the local alignment routine on the set CVWVCITFIAFQYTADPVMHMVE|DKDGSPLSVL, exit status 6 on every run); the port returns the alignment
                if (c.Flags.Contains("gotoh_pair_wise_sw", StringComparison.Ordinal))
                {
                    var swRows = SequenceAligner.MultipleAlignTCoffeeProgram(c.Seqs, options).AlignedSequences.ToArray();
                    Assert.That(swRows, Is.Not.Empty);
                    observed.Add($"{row}\tsw_uninitialised_write_abort\tport returns {swRows.Length} rows of {swRows[0].Length} columns");
                    continue;
                }

                Assert.Throws<TCoffeeAbortedException>(() => SequenceAligner.MultipleAlignTCoffeeProgram(c.Seqs, options), $"t_coffee aborts on '{c.Flags}' {string.Join('|', c.Seqs)}");
                continue;
            }

            try
            {
                string got = string.Join('|', SequenceAligner.MultipleAlignTCoffeeProgram(c.Seqs, options).AlignedSequences);
                if (got == c.Rows) continue;
                // idscore_pairseq reads m[s1[i] - 'a'][s2[j] - 'a']: with a gap in the SECOND aligned row that is the heap 208 bytes before a matrix row - 0 in 23 of the 24 committed sets that have such a gap,
                // some other value in one (undefined behaviour of the program, not reproducible)
                if (IsIdscoreFast(c.Flags) && HasSecondRowGap(c, options)) observed.Add($"{row}\tidscore_heap_value_diff\tport returns other rows");
                // the same one-byte row overflow can also leave the program's rows without their last column (1 committed biphasic_pair_wise set)
                else if (IsGotohDp(c.Flags) && c.Rows.Split('|').Zip(got.Split('|'), (e, g) => g.Length == e.Length + 1 && g.StartsWith(e, StringComparison.Ordinal)).All(x => x)) observed.Add($"{row}\tgotoh_row_overflow_last_column\tthe program's rows lack the last column");
                else if (c.Flags.Contains("test_pair_wise", StringComparison.Ordinal) && TestPairOverflow(c, options, out int excessColumns)) observed.Add($"{row}\ttest_pair_row_overflow_diff\texcess={excessColumns}"); // 32 sets lose the last column of every row, 9 have another row clobbered by the overflow
                else if (c.Flags.Contains("gotoh_pair_wise_sw", StringComparison.Ordinal)) observed.Add($"{row}\tsw_uninitialised_write_diff\tport rows differ");
                else failures.Add($"{c.Type} '{c.Flags}' {string.Join('|', c.Seqs)}: got {got} expected {c.Rows}");
            }
            catch (TCoffeeAbortedException e) when (e.HeapDependent)
            {
                // undefined behaviour of the program (a wild pointer read whose outcome depends on the heap): it may also run through; only such sets may differ this way
                observed.Add($"{row}\tidscore_heap_abort\tport aborts (heap dependent)");
            }
            catch (NotSupportedException) when (c.Flags.Contains("relative_triplet", StringComparison.Ordinal))
            {
                // residue_pair_relative_extended_list reads past the end of a library list (undefined behaviour): not reproducible, the port says so
                observed.Add($"{row}\trelative_triplet_ub\tNotSupportedException (the program ran through)");
            }
        }

        if (Environment.GetEnvironmentVariable("SEQERON_TCOFFEE_UB_CASES_OUT") is { Length: > 0 } outPath)
            File.AppendAllLines(outPath, observed.Select(o => feature + "\t" + o));
        Assert.That(failures, Is.Empty, $"{failures.Count} of {cases.Count} '{feature}' sets differ from t_coffee; first: {failures.FirstOrDefault()}");
        Assert.That(observed, Is.EqualTo(UbCases.Value.TryGetValue(feature, out var expected) ? expected : []), $"the undefined-behaviour rows of '{feature}' are exactly the enumerated ones (tcoffee_options_ub_cases.tsv)");
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
    /// A2-11 / A2-19: <c>-method mafft_pair</c> runs the external MAFFT 7.505 (<c>mafft --anysymbol in &gt; out</c>) on every pair (the installed program works once <c>$HOME/.t_coffee/plugins/linux</c> holds mafft's helper binaries:
    /// the oracle script links them); the pair stage of MAFFT is ported (<see cref="SequenceAligner.MultipleAlignMafft"/>, F119): DNA, protein and mixed alphabets, with a second method, UPGMA and no relaxation.
    /// </summary>
    [Test]
    public void MafftPair_EveryRowEqualsTheProgram()
    {
        var cases = Oracle().Where(c => c.Feature is "mafft_pair" or "mafft_mixed").ToList();
        Assert.That(cases.Count, Is.GreaterThanOrEqualTo(440), "oracle sets of mafft_pair");
        Assert.That(cases.Count(c => c.Rows.StartsWith("FAIL:", StringComparison.Ordinal)), Is.EqualTo(0), "the program runs mafft_pair on every set");
        AssertFeature("mafft_pair", 350);
        AssertFeature("mafft_mixed", 100);
    }

    /// <summary>
    /// A2-11: <c>-method muscle_pair</c> = the external MUSCLE 3.8.1551 on every pair (the installed t_coffee finds it only through <c>$HOME/.t_coffee/plugins/linux/muscle</c>, the oracle script links muscle3 there; the original
    /// "rc 1" of round 1 was the missing helper binaries of mafft / the missing <c>muscle</c> name, not an impossibility); DNA, protein and mixed alphabets, with a second method, UPGMA and no relaxation.
    /// </summary>
    [Test]
    public void MusclePair_EveryRowEqualsTheProgram()
    {
        AssertFeature("muscle_pair", 350);
        AssertFeature("muscle_mixed", 100);
    }

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
        // -distance_matrix_mode random: implemented since A3-6 (glibc rand, seed 10 = the program's DEBUG_SRAND run), see RandomDistance_EveryRowEqualsTheProgram_WithDebugSrand
        Assert.DoesNotThrow(() => SequenceAligner.MultipleAlignTCoffeeProgram(new[] { "ACGT", "ACGA", "ACGG" }, TCoffeeProgramOptions.FromCommandLine("-distance_matrix_mode random")));
    }

    /// <summary>
    /// A2-2: <c>fast_pair</c> = <c>fasta_gotoh_pair_wise</c> (util_dp_fasta_nw.c): k-tuple diagonals (2 protein / 5 nucleotide), the z-score threshold with T = max z / 25, the banded Gotoh DP whose arrays hold the program's initial
    /// gap scores outside the band, terminal gaps paying the extension only. <c>-nomatch</c> and <c>-ndiag</c> are copied into the program's CL AFTER the library methods were configured and therefore do not reach <c>fast_pair</c>
    /// (the oracle contains both flags: 0 / 510 differences).
    /// </summary>
    [Test]
    public void FastPair_EveryRowEqualsTheProgram() => AssertFeature("fast_pair", 500);

    /// <summary>
    /// A2-6: the distance modes that need the RELAXED library, computed between the relaxation and the re-weighting as the program does: <c>slow</c> (linked pair-wise DP score per column, rescaled by do_normalise, float),
    /// <c>cscore</c> (the same score in float without the rescaling), <c>clustalw</c> (the DP score is the identity) and <c>very_fast</c> / <c>-quicktree</c> (<c>very_fast_gotoh_pair_wise</c>: diagonals scored by library entries,
    /// the 20 best with a shadow of 5, <c>-tg_mode</c>), with t_coffee sequence weights, UPGMA, filtered libraries and several methods.
    /// </summary>
    [Test]
    public void DistanceModesOfTheRelaxedLibrary_EveryRowEqualsTheProgram()
    {
        AssertFeature("dist_slow", 400);
        AssertFeature("dist_cscore", 350);
        AssertFeature("dist_clustalw", 300);
        AssertFeature("dist_veryfast", 400);
    }

    /// <summary>A2-6: <c>-distance_matrix_mode geometric</c>: the identity vectors of the <c>fasta_pair_wise</c> alignments (pam250mt, gap open -10 / extend -1, k-tuple 3, no terminal gap penalty) against the first three sequences.</summary>
    [Test]
    public void GeometricDistance_EveryRowEqualsTheProgram() => AssertFeature("dist_geometric", 300);

    /// <summary>
    /// <c>-distance_matrix_mode random</c> (audit A3-6, B09 F135): <c>id = rand() % 100</c> per pair. T-Coffee seeds with <c>vsrand(0)</c> (/dev/urandom: two runs differ) unless <c>DEBUG_SRAND</c> is set (<c>srand(10)</c>);
    /// the 310 committed sets were recorded with <c>DEBUG_SRAND=1</c> and equal the port's default <see cref="TCoffeeProgramOptions.RandomSeed"/> 10 (glibc <c>rand</c> emulation).
    /// </summary>
    [Test]
    public void RandomDistance_EveryRowEqualsTheProgram_WithDebugSrand() => AssertFeature("dist_random", 300);

    [Test]
    public void RandomDistance_SeedIsAnOption_AndTheDefaultIsTheProgramsDebugSrand()
    {
        string[] seqs = ["ACGTACGTTAGCACGTAGCT", "ACGTACGTAGCACGTTGCT", "ACGTTCGTTAGCAAGTAGCT", "TTGTACGTTAGCACGTAGCA", "ACGAACGTTGCACGTAGCT", "ACGTACCTTAGCACGAAGCT"];
        var random = TCoffeeProgramOptions.FromCommandLine("-distance_matrix_mode random -seq_weight t_coffee");
        Assert.Multiple(() =>
        {
            Assert.That(random.RandomSeed, Is.EqualTo(10u));
            var a = SequenceAligner.MultipleAlignTCoffeeProgram(seqs, random).AlignedSequences;
            Assert.That(SequenceAligner.MultipleAlignTCoffeeProgram(seqs, random).AlignedSequences, Is.EqualTo(a), "deterministic for a seed");
            Assert.That(SequenceAligner.MultipleAlignTCoffeeProgram(seqs, random with { RandomSeed = 10 }).AlignedSequences, Is.EqualTo(a));
            bool anyDifferent = Enumerable.Range(1, 40).Any(seed => !SequenceAligner.MultipleAlignTCoffeeProgram(seqs, random with { RandomSeed = (uint)seed }).AlignedSequences.SequenceEqual(a));
            Assert.That(anyDifferent, Is.True, "the seed changes the random distances (and so the tree / weights)");
        });
    }

    /// <summary>
    /// A2-4: <c>exon_pair</c> / <c>exon2_pair</c> / <c>exon3_pair</c> (<c>SEQ_TYPE G</c>): the Myers-Miller alignment of the genomic (<c>_G_</c>) templates with the matrix exon2mt (o / b / j = 100, x rows zeroed), expressed in the residues of
    /// the sequences (<c>fix_seq_seq</c>), only pairs whose two sequences have a template; exon2 / exon3 keep the columns whose letters are both in the quirky set "subset_objOBJ-"; a reconciling alignment inside <c>fix_seq_seq</c> calls
    /// <c>aln2sim(B, "idmat")</c> first and so fixes the process-wide residue groups of the weights to 26 single letters (T and U differ, J / O equal themselves) instead of the groups of the matrix idmat. Rows without any pair raise
    /// <see cref="TCoffeeAbortedException"/> as the program stops ("You have not provided any method or enough Sequences").
    /// </summary>
    [Test]
    public void ExonPairs_EveryRowEqualsTheProgram()
    {
        AssertFeature("exon", 400);
        AssertFeature("exon_sparse", 200);
        Assert.That(Oracle().Count(c => c.Feature.StartsWith("exon", StringComparison.Ordinal) && c.Rows.StartsWith("FAIL:", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(30), "sets on which the program stops");
        Assert.That(Oracle().Count(c => c.Feature.StartsWith("exon", StringComparison.Ordinal) && c.Flags.Contains("exon3_pair", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(150));
        Assert.That(SequenceAligner.TCoffeePrimaryLibraryEntries(new[] { "TATTTGTTATT", "TAGCAGACGCGACCTGTAGGTACCACTTCTTGC" }, new[] { TCoffeePairMethod.ExonPair }, genomicTemplates: new[] { "TAtTtgTtatt", "TtaugAjcgcgCucggtaogtAbctoajtTtg" })
            .Select(e => (e.PosA + 1, e.PosB + 1, e.Weight)), Is.EqualTo(new[] { (2, 22, 450), (3, 23, 450), (4, 24, 450), (7, 25, 450), (8, 26, 450), (9, 27, 450), (10, 28, 450), (11, 31, 450) }), "library of t_coffee -lib_only -out_lib");
    }

    /// <summary>
    /// A2-5: <c>-dp_mode gotoh_pair_wise</c> (+ <c>ssec_pair_wise</c>, <c>gotoh_pair_wise_test</c>), <c>myers_miller_pair_wise</c> (+ <c>default</c>), <c>fasta_pair_wise</c>, <c>cfasta_pair_wise</c>, <c>very_fast_pair_wise</c> on the library column scores
    /// of <c>-extend_mode</c> (very_fast_triplet: the <c>-max_n_pair</c> most similar pairs of sequences by the identity of the distance matrix; slow_triplet; triplet: nomatch twice; fast_triplet), with gap penalties, <c>-tg_mode</c>,
    /// <c>-ndiag</c>, <c>-diag_threshold</c>, sequence weights, UPGMA and several methods. 2 070 sets; the two sets on which the program dies with exit status 6 are the heap overflow of gotoh_pair_wise described in
    /// <see cref="AssertFeature"/>.
    /// </summary>
    [Test]
    public void GenericDpModes_EveryRowEqualsTheProgram()
    {
        AssertFeature("dp_gotoh", 700);
        AssertFeature("dp_mm", 450);
        AssertFeature("dp_fasta", 700);
    }

    /// <summary>A2-5: <c>-dp_mode biphasic_pair_wise</c> (two affine gap models) and <c>glocal_pair_wise</c> (match / deletion / insertion / unmatched pair) on the column scores of <c>-extend_mode</c>; 550 sets, 6 of them hit the row-writer heap overflow.</summary>
    [Test]
    public void BiphasicAndGlocalDp_EveryRowEqualsTheProgram() => AssertFeature("dp_other", 450);

    /// <summary>
    /// A2-5 (WP30): <c>-dp_mode gotoh_pair_wise_sw</c> (the Smith-Waterman form: every merge keeps only the columns of the best local alignment of the two profiles, so the rows lose residues and have different lengths; 500 sets, 450 equal,
    /// 50 undefined behaviour of the program) and <c>test_pair_wise</c> (a debugging routine without an alignment: rows padded with gaps; 420 sets, see <see cref="TestPairDp"/>).
    /// </summary>
    [Test]
    public void LocalAndTestDp_EveryRowEqualsTheProgram()
    {
        AssertFeature("dp_sw", 450);
        AssertFeature("dp_test", 400);
    }

    /// <summary>A2-5 (WP30): <c>-dp_mode gotoh_pair_wise_lgp</c> (the Gotoh DP with the position specific gap opening / closing penalties of <c>aln2local_penalties</c>) on the column scores of <c>-extend_mode</c>; 620 sets.</summary>
    [Test]
    public void GotohLgpDp_EveryRowEqualsTheProgram() => AssertFeature("dp_lgp", 600);

    /// <summary>
    /// A2-5 (WP30): <c>-dp_mode ktup_pair_wise</c> computes a k-tuple score only (the rows stay as they are, unpadded) - 250 sets; <c>precomputed_pair_wise</c> stops with "improper use of the function precomputed pairwise"
    /// (exit status 1) whenever the first rows of two groups differ in length (120 of 120 sets of unequal sequences) and leaves the rows alone otherwise (220 sets of equal length).
    /// </summary>
    [Test]
    public void KtupAndPrecomputedDp_EveryRowEqualsTheProgram()
    {
        AssertFeature("dp_ktup", 240);
        AssertFeature("dp_precomputed_equal", 200);
        AssertFeature("dp_precomputed", 100);
        Assert.That(Oracle().Where(c => c.Feature == "dp_precomputed").All(c => c.Rows == "FAIL:1"), Is.True);
    }

    /// <summary>
    /// A2-5 (WP30): the column scores <c>-extend_mode pc | mixt | g_coffee | g_coffee_quadruplets</c> (<c>slow_get_dp_cost</c> over <c>residue_pair_extended_list_pc</c> / <c>_mixt</c> / <c>_g_coffee</c> / <c>_g_coffee_quadruplet</c>) on
    /// <c>gotoh_pair_wise</c>, <c>myers_miller_pair_wise</c> and <c>fasta_pair_wise</c>, with gap penalties, nomatch, sequence weights and UPGMA; 760 sets (6 = row-writer overflow of the program).
    /// </summary>
    [Test]
    public void PcMixtAndGCoffeeExtendModes_EveryRowEqualsTheProgram() => AssertFeature("dp_extend2", 700);

    /// <summary>A2-5 (WP30): <c>-dp_mode procoffee_pair_wise</c> = the sparse linked DP over the scored column pairs with the older <c>cl2diag_cap390</c> (cl2list_borders first, borders inside the cap) and the two groups swapped; 460 sets.</summary>
    [Test]
    public void ProcoffeeDp_EveryRowEqualsTheProgram() => AssertFeature("dp_procoffee", 450);

    /// <summary>A2-5 (WP30): <c>-dp_mode viterbi_pair_wise</c> = the Viterbi path of the five-state pair-HMM of proba_pair (protein model, whatever the sequence type) between the residue frequency profiles of the two groups; 420 sets.</summary>
    [Test]
    public void ViterbiDp_EveryRowEqualsTheProgram() => AssertFeature("dp_viterbi", 400);

    /// <summary>A2-5 (WP30): <c>-extend_mode quadruplet</c> (<c>get_dp_cost_quadruplet</c>: the quadruplet score of every residue pair when a group has one sequence, else the hashes of the partners of the partners of the two columns; nomatch subtracted twice); 360 sets (1 = row-writer overflow).</summary>
    [Test]
    public void QuadrupletExtendMode_EveryRowEqualsTheProgram() => AssertFeature("dp_quad", 350);

    /// <summary>
    /// A2-5 (WP30): <c>-extend_mode test</c> (<c>slow_get_dp_cost_test</c>: the raw sum of the smaller weights x 10 / (ns1 x ns2), no nomatch) and <c>rna</c> / <c>rna1</c> / <c>rna2</c> / <c>rna3</c> / <c>rna4</c> without a contact library (the column score of rna0 = slow_triplet) on
    /// <c>gotoh_pair_wise</c> / <c>myers_miller_pair_wise</c> / <c>fasta_pair_wise</c>; 470 sets (3 = row-writer overflow of the program).
    /// </summary>
    [Test]
    public void TestAndRnaExtendModes_EveryRowEqualsTheProgram() => AssertFeature("dp_extend3", 450);

    /// <summary>
    /// A2-5: <c>-extend_mode rna0</c> (= slow_triplet), the extend modes without effect on the sparse linked DP, and <c>relative_triplet</c>, whose <c>residue_pair_relative_extended_list</c> reads the weights of the first residue's library list
    /// at the positions of the second residue's (longer) list, past the array: undefined behaviour (the port raises <see cref="NotSupportedException"/> on such a read and equals the program on every set where none happens).
    /// </summary>
    [Test]
    public void ExtendModes_EveryRowEqualsTheProgram_RelativeTripletIsUndefinedBehaviour() => AssertFeature("dp_extend", 150);

    /// <summary>
    /// A2-5: the dynamic programs and extend modes the program refuses or dies on, measured on 2-5 DNA / protein sequences: <c>-dp_mode</c> idscore / subop1 / subop2 / proba / pavie / domain / fasta_cdna / cfasta_cdna (exit status 11),
    /// gotoh_pair_wise_lalign (6), seq_viterbi / hh / cfasta_sw / an unknown name (1); <c>-extend_mode</c> unknown / matrix (1, for any number of sequences), ssp / tm (11) and curvature (1) as soon as there are two sequences.
    /// </summary>
    [Test]
    public void DpAndExtendModesTheProgramRefuses_RaiseTCoffeeAbortedException()
    {
        AssertFeature("dp_refused", 600);
        Assert.That(Oracle().Where(c => c.Feature == "dp_refused").All(c => c.Rows.StartsWith("FAIL:", StringComparison.Ordinal)), Is.True);
        // one sequence: no pair job - but an unknown extend mode stops choose_extension_mode anyway, ssp / tm / curvature need a pair
        Assert.That(SequenceAligner.MultipleAlignTCoffeeProgram(new[] { "ACGT" }, TCoffeeProgramOptions.FromCommandLine("-dp_mode idscore_pair_wise -extend_mode ssp")).AlignedSequences, Is.EqualTo(new[] { "ACGT" }));
        Assert.Throws<TCoffeeAbortedException>(() => SequenceAligner.MultipleAlignTCoffeeProgram(new[] { "ACGT" }, TCoffeeProgramOptions.FromCommandLine("-extend_mode bogus")));
        Assert.Throws<NotSupportedException>(() => SequenceAligner.MultipleAlignTCoffeeProgram(new[] { "ACGT", "ACGA" }, TCoffeeProgramOptions.FromCommandLine("-dp_mode viterbiL_pair_wise")));
    }

    [Test]
    public void FromCommandLine_FastaOptionsQuicktreeAndTemplates()
    {
        // -quicktree replaces the distance mode wherever it stands on the line
        Assert.That(TCoffeeProgramOptions.FromCommandLine("-quicktree -distance_matrix_mode ktup").DistanceMatrixMode, Is.EqualTo(TCoffeeDistanceMode.VeryFast));
        var o = TCoffeeProgramOptions.FromCommandLine("-method exon3_pair,fast_pair -template_file @AC,,GT -ndiag 7 -ktuple 3 -tg_mode 0");
        Assert.That(o.Methods, Is.EqualTo(new[] { TCoffeePairMethod.Exon3Pair, TCoffeePairMethod.FastPair }));
        Assert.That(o.GenomicTemplates, Is.EqualTo(new string?[] { "AC", null, "GT" }));
        Assert.That((o.DiagonalCount, o.Ktuple, o.TgMode), Is.EqualTo((7, 3, 0)));
        Assert.Throws<ArgumentException>(() => TCoffeeProgramOptions.FromCommandLine("-template_file somefile.fa"));
        Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignTCoffeeProgram(new[] { "ACGT", "ACGA" }, TCoffeeProgramOptions.FromCommandLine("-method exon_pair -template_file @AC")));
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
