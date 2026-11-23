using System.Globalization;

namespace Seqeron.Genomics.Alignment;

public sealed partial record TCoffeeProgramOptions
{
    /// <summary>
    /// Reads the flags of a <c>t_coffee</c> command line (everything after the program name and without the input / output file flags
    /// <c>-seq</c>, <c>-outfile</c>, <c>-output</c>, <c>-n_core</c>, <c>-no_warning</c>, which have no influence on the rows) into options.
    /// A token starting with '-' and a letter is a flag; the following tokens are its values (negative numbers are values). Several <c>-method</c> values may be
    /// given as separate tokens or comma separated. A flag that has no implementation here raises <see cref="ArgumentException"/> naming it.
    /// </summary>
    /// <param name="commandLine">For example <c>-method slow_pair,lalign_id_pair -relax_lib 0</c>.</param>
    public static TCoffeeProgramOptions FromCommandLine(string commandLine) => new TCoffeeProgramOptions().WithCommandLine(commandLine);

    /// <summary>Applies the flags of a <c>t_coffee</c> command line on top of these options (see <see cref="FromCommandLine"/>).</summary>
    /// <param name="commandLine">The flags.</param>
    public TCoffeeProgramOptions WithCommandLine(string commandLine)
    {
        ArgumentNullException.ThrowIfNull(commandLine);
        var tokens = commandLine.Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
        var options = this;
        bool weightNo = false, quickTree = false;
        int i = 0;
        while (i < tokens.Length)
        {
            string flag = tokens[i++];
            if (flag.Length < 2 || flag[0] != '-' || !char.IsLetter(flag[1]))
                throw new ArgumentException($"Unexpected token '{flag}' in the command line.", nameof(commandLine));
            var values = new List<string>();
            while (i < tokens.Length && !(tokens[i].Length > 1 && tokens[i][0] == '-' && char.IsLetter(tokens[i][1]))) values.Add(tokens[i++]);
            if (flag == "-quicktree")
            {
                quickTree = true; // `-quicktree` replaces the distance mode by very_fast wherever it stands on the command line (t_coffee.c L2103)
                continue;
            }

            if (flag == "-weight")
            {
                string w = StringValue(flag, values);
                if (w == "no") weightNo = true;
                else if (w != "default") throw new ArgumentException($"-weight {w} is not supported (library weights of the methods).");
                continue;
            }

            options = options.Apply(flag, values);
        }

        if (quickTree) options = options with { DistanceMatrixMode = TCoffeeDistanceMode.VeryFast };

        // `-weight no` skips the sequence weighting altogether
        return weightNo ? options with { TCoffeeSequenceWeights = false, SequenceWeights = null, WeightNo = true } : options;
    }

    private static int IntValue(string flag, List<string> values)
    {
        if (values.Count != 1 || !int.TryParse(values[0], NumberStyles.Integer, CultureInfo.InvariantCulture, out int v))
            throw new ArgumentException($"The flag {flag} needs one integer value.");
        return v;
    }

    private static int MethodInt(string spec, string key, string text)
    {
        if (!int.TryParse(text, NumberStyles.Integer, CultureInfo.InvariantCulture, out int v))
            throw new ArgumentException($"-method {spec}: {key} needs an integer value, got '{text}'.");
        return v;
    }

    private static string StringValue(string flag, List<string> values)
    {
        if (values.Count != 1) throw new ArgumentException($"The flag {flag} needs one value.");
        return values[0];
    }

    private TCoffeeProgramOptions Apply(string flag, List<string> values)
    {
        switch (flag)
        {
            case "-seq":
            case "-outfile":
            case "-output":
            case "-n_core":
            case "-no_warning":
            case "-quiet":
            case "-multi_core":
                return this; // no influence on the rows (checked: -n_core 1, 2, 3, 4, 8)
            case "-mode":
            case "-special_mode":
                // 3dcoffee = `-in Msap_pair -template_file SELF_P_ -profile_template_file SELF_P_`: the method list is sap_pair ALONE (the default proba_pair is not added); the templates are looked up by the sequence NAMES (network): the structures of
                // the sequences are given as TCoffeeProgramOptions.StructureTemplates, and when the names are no PDB entries (the oracle: q000, ...) every pair is aligned by proba_pair. expresso fetches the structures from the RCSB (network)
                return StringValue(flag, values) switch
                {
                    "3dcoffee" => this with { Methods = [TCoffeePairMethod.SapPair] },
                    "expresso" => throw new NotSupportedException("-mode expresso downloads the structures of the sequences from the RCSB PDB (blast against the PDB, network): not reproducible offline."),
                    var m => throw new ArgumentException($"The t_coffee mode {m} is not supported (3dcoffee)."),
                };
            case "-method":
                {
                    var methods = new List<TCoffeePairMethod>();
                    var parameters = new Dictionary<TCoffeePairMethod, TCoffeeMethodParameters>();
                    foreach (string spec in values.SelectMany(v => v.Split(',', StringSplitOptions.RemoveEmptyEntries)))
                    {
                        // name@EP@KEY@value@KEY@value: every pair appends "KEY value" to the method file
                        var parts = spec.Split('@');
                        var method = ParseMethod(parts[0]);
                        methods.Add(method);
                        if (parts.Length > 1)
                        {
                            if (parts[1] != "EP" || parts.Length % 2 != 0) throw new ArgumentException($"-method {spec}: expected name@EP@KEY@value[@KEY@value...].");
                            var mp = new TCoffeeMethodParameters();
                            for (int a = 2; a < parts.Length; a += 2)
                                mp = parts[a] switch
                                {
                                    "GOP" => mp with { GapOpen = MethodInt(spec, "GOP", parts[a + 1]) },
                                    "GEP" => mp with { GapExtend = MethodInt(spec, "GEP", parts[a + 1]) },
                                    "MATRIX" => mp with { Matrix = SequenceAligner.TcResolveMatrix(parts[a + 1]) is null ? null : parts[a + 1] },
                                    _ => throw new ArgumentException($"-method {spec}: the method parameter {parts[a]} is not supported (GOP, GEP, MATRIX)."),
                                };
                            parameters.TryAdd(method, mp);
                        }
                    }

                    if (methods.Count == 0) throw new ArgumentException("-method needs at least one method.");
                    return this with { Methods = methods, MethodParameters = parameters.Count > 0 ? parameters : null };
                }
            case "-extend_mode": // no effect on the default dp_mode linked_pair_wise (relaxation + sparse DP), 0 / 120 sets differ; it scores the profile columns of the generic dp modes
                return StringValue(flag, values) switch
                {
                    "very_fast_triplet" => this with { ExtendMode = TCoffeeExtendMode.VeryFastTriplet },
                    "fast_triplet" => this with { ExtendMode = TCoffeeExtendMode.FastTriplet },
                    "slow_triplet" => this with { ExtendMode = TCoffeeExtendMode.SlowTriplet },
                    "triplet" => this with { ExtendMode = TCoffeeExtendMode.Triplet },
                    "relative_triplet" => this with { ExtendMode = TCoffeeExtendMode.RelativeTriplet },
                    "pc" => this with { ExtendMode = TCoffeeExtendMode.Pc },
                    "test" => this with { ExtendMode = TCoffeeExtendMode.Test },
                    "quadruplet" => this with { ExtendMode = TCoffeeExtendMode.Quadruplet },
                    "mixt" => this with { ExtendMode = TCoffeeExtendMode.Mixt },
                    "g_coffee" => this with { ExtendMode = TCoffeeExtendMode.GCoffee },
                    "g_coffee_quadruplets" => this with { ExtendMode = TCoffeeExtendMode.GCoffeeQuadruplets },
                    var m => this with { ExtendMode = TCoffeeExtendMode.Other, ExtendModeText = m },
                };
            case "-dp_mode":
                return StringValue(flag, values) switch
                {
                    "linked_pair_wise" => this with { DpMode = TCoffeeDpMode.Linked },
                    "gotoh_pair_wise" or "ssec_pair_wise" => this with { DpMode = TCoffeeDpMode.Gotoh },
                    "gotoh_pair_wise_test" => this with { DpMode = TCoffeeDpMode.GotohTest },
                    "gotoh_pair_wise_lgp" => this with { DpMode = TCoffeeDpMode.GotohLgp },
                    "ktup_pair_wise" => this with { DpMode = TCoffeeDpMode.Ktup },
                    "precomputed_pair_wise" => this with { DpMode = TCoffeeDpMode.Precomputed },
                    "gotoh_pair_wise_sw" => this with { DpMode = TCoffeeDpMode.GotohSw },
                    "test_pair_wise" => this with { DpMode = TCoffeeDpMode.Test },
                    "procoffee_pair_wise" => this with { DpMode = TCoffeeDpMode.Procoffee },
                    "viterbi_pair_wise" => this with { DpMode = TCoffeeDpMode.Viterbi },
                    "myers_miller_pair_wise" or "default" => this with { DpMode = TCoffeeDpMode.MyersMiller },
                    "fasta_pair_wise" => this with { DpMode = TCoffeeDpMode.Fasta },
                    "cfasta_pair_wise" => this with { DpMode = TCoffeeDpMode.CFasta },
                    "very_fast_pair_wise" => this with { DpMode = TCoffeeDpMode.VeryFast },
                    "biphasic_pair_wise" => this with { DpMode = TCoffeeDpMode.Biphasic },
                    "glocal_pair_wise" => this with { DpMode = TCoffeeDpMode.Glocal },
                    var m => this with { DpMode = TCoffeeDpMode.Other, DpModeText = m },
                };
            case "-max_n_pair":
                return this with { MaxNPair = IntValue(flag, values) };
            case "-diag_threshold":
                return this with { DiagonalThreshold = IntValue(flag, values) };
            case "-extend":
            case "-evaluate_mode":
            case "-method_evaluate_mode":
            case "-shrink_lib": // only recounts the entries (shrink_constraint_list sets CL->ne), 0 / 120 sets differ
            case "-tg_mode": // the linked DP has its own gap rule (terminal gaps free of the opening penalty): no effect there (0 / 120 sets differ); it rules the terminal gaps of the heuristic DP modes
                return this with { TgMode = IntValue(flag, values) };
            case "-ndiag":
                return this with { DiagonalCount = IntValue(flag, values) };
            case "-ktuple":
                return this with { Ktuple = IntValue(flag, values) };
            case "-tree":
            case "-newtree":
            case "-outseqweight":
            case "-out_lib":
            case "-out_lib_mode":
            case "-clean_seq_name":
            case "-display":
                return this;
            case "-matrix":
                {
                    string m = StringValue(flag, values);
                    SequenceAligner.TcResolveMatrix(m); // rejects a name the program does not know
                    return this with { Matrix = m == "default" ? null : m };
                }
            case "-template_file":
                {
                    // the in-memory form of a template list: @t0,t1,... = the genomic (_G_) template sequence of every input sequence in input order (empty = none)
                    string v = StringValue(flag, values);
                    if (!v.StartsWith('@')) throw new ArgumentException("-template_file: only the in-memory form @template0,template1,... (one genomic template per sequence, empty = none) is supported.");
                    return this with { GenomicTemplates = v[1..].Split(',').Select(t => t.Length == 0 ? null : t).ToArray() };
                }
            case "-lalign_n_top":
                return this with { LalignTopAlignments = IntValue(flag, values) };
            case "-relax_lib":
                {
                    int n = IntValue(flag, values);
                    if (n < 0) throw new ArgumentException("-relax_lib must not be negative.");
                    return this with { RelaxLibIterations = n };
                }
            case "-filter_lib":
                return this with { FilterLib = IntValue(flag, values) };
            case "-tree_mode":
                return StringValue(flag, values) switch
                {
                    "nj" => this with { TreeMode = TCoffeeTreeMode.Nj },
                    "upgma" => this with { TreeMode = TCoffeeTreeMode.Upgma },
                    "cwph" => this with { TreeMode = TCoffeeTreeMode.Cwph },
                    "kmeans" => this with { TreeMode = TCoffeeTreeMode.Kmeans },
                    var m => this with { TreeMode = TCoffeeTreeMode.Other, TreeModeText = m },
                };
            case "-distance_matrix_mode":
                return StringValue(flag, values) switch
                {
                    "ktup" => this with { DistanceMatrixMode = TCoffeeDistanceMode.Ktup },
                    "fast" => this with { DistanceMatrixMode = TCoffeeDistanceMode.Fast },
                    "idscore" => this with { DistanceMatrixMode = TCoffeeDistanceMode.Idscore },
                    "no" => this with { DistanceMatrixMode = TCoffeeDistanceMode.No },
                    "aln" => this with { DistanceMatrixMode = TCoffeeDistanceMode.Aln },
                    "ktup2" => this with { DistanceMatrixMode = TCoffeeDistanceMode.Ktup2 },
                    "slow" => this with { DistanceMatrixMode = TCoffeeDistanceMode.Slow },
                    "cscore" => this with { DistanceMatrixMode = TCoffeeDistanceMode.Cscore },
                    "geometric" => this with { DistanceMatrixMode = TCoffeeDistanceMode.Geometric },
                    "random" => this with { DistanceMatrixMode = TCoffeeDistanceMode.Random },
                    "clustalw" => this with { DistanceMatrixMode = TCoffeeDistanceMode.Clustalw },
                    "very_fast" => this with { DistanceMatrixMode = TCoffeeDistanceMode.VeryFast },
                    var m => this with { DistanceMatrixMode = TCoffeeDistanceMode.Other, DistanceModeText = m },
                };
            case "-distance_matrix_sim_mode":
                return this with { DistanceMatrixSimMode = StringValue(flag, values) };
            case "-seq_weight":
                {
                    string v = values.Count == 0 ? "t_coffee" : StringValue(flag, values);
                    if (v == "t_coffee") return this with { TCoffeeSequenceWeights = true };
                    if (v.StartsWith('@'))
                    {
                        // the in-memory form of a weight file: @w1,w2,... in input order
                        var w = v[1..].Split(',').Select(x => float.TryParse(x, NumberStyles.Float, CultureInfo.InvariantCulture, out float f)
                            ? f
                            : throw new ArgumentException($"-seq_weight: '{x}' is not a number.")).ToArray();
                        return this with { SequenceWeights = w };
                    }

                    return this with { TCoffeeSequenceWeights = false, SequenceWeights = null }; // any other word: no sequence weights
                }
            case "-gapopen":
                return this with { GapOpen = IntValue(flag, values) };
            case "-gapext":
                return this with { GapExtend = IntValue(flag, values) };
            case "-cosmetic_penalty":
                return this with { CosmeticPenalty = IntValue(flag, values) };
            case "-do_normalise":
                return this with { DoNormalise = IntValue(flag, values) };
            case "-nomatch":
                return this with { NoMatch = IntValue(flag, values) };
            case "-outorder":
                return StringValue(flag, values) switch
                {
                    "aligned" => this with { OutOrderAligned = true },
                    "input" => this with { OutOrderAligned = false },
                    var m => throw new ArgumentException($"-outorder {m} is not supported (input, aligned)."),
                };
            case "-inorder":
                return StringValue(flag, values) switch
                {
                    "aligned" => this with { InOrderInput = false },
                    "input" => this with { InOrderInput = true },
                    var m => throw new ArgumentException($"-inorder {m} is not supported (aligned, input)."),
                };
            default:
                throw new ArgumentException($"The t_coffee flag {flag} is not supported.");
        }
    }

    private static TCoffeePairMethod ParseMethod(string name) => name switch
    {
        "slow_pair" => TCoffeePairMethod.SlowPair,
        "lalign_id_pair" => TCoffeePairMethod.LalignIdPair,
        "proba_pair" => TCoffeePairMethod.ProbaPair,
        "clustalw_pair" => TCoffeePairMethod.ClustalwPair,
        "biphasic_pair" => TCoffeePairMethod.BiphasicPair,
        "lalign_rs_s_pair" => TCoffeePairMethod.LalignRsSPair,
        "promo_pair" => TCoffeePairMethod.PromoPair,
        "blastr_pair" => TCoffeePairMethod.BlastrPair,
        "clean_slow_pair" => TCoffeePairMethod.CleanSlowPair,
        "hash_pair" => TCoffeePairMethod.HashPair,
        "subop1_pair" => TCoffeePairMethod.Subop1Pair,
        "subop2_pair" => TCoffeePairMethod.Subop2Pair,
        "cfast_pair" => TCoffeePairMethod.CfastPair,
        "fast_pair" => TCoffeePairMethod.FastPair,
        "exon_pair" => TCoffeePairMethod.ExonPair,
        "exon2_pair" => TCoffeePairMethod.Exon2Pair,
        "exon3_pair" => TCoffeePairMethod.Exon3Pair,
        "muscle_pair" => TCoffeePairMethod.MusclePair,
        "mafft_pair" => TCoffeePairMethod.MafftPair,
        "TMalign_pair" or "tmalign_pair" => TCoffeePairMethod.TMalignPair,
        "mustang_pair" => TCoffeePairMethod.MustangPair,
        "blastp_msa" => TCoffeePairMethod.BlastpMsa,
        "tblastx_msa" => TCoffeePairMethod.TblastxMsa,
        "sap_pair" => TCoffeePairMethod.SapPair,
        "dali_pair" => TCoffeePairMethod.DaliPair,
        "fugue_pair" => TCoffeePairMethod.FuguePair,
        _ => throw new ArgumentException($"The t_coffee method {name} is not supported."),
    };
}
