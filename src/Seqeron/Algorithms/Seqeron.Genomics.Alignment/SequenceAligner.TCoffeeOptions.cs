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
        bool weightNo = false;
        int i = 0;
        while (i < tokens.Length)
        {
            string flag = tokens[i++];
            if (flag.Length < 2 || flag[0] != '-' || !char.IsLetter(flag[1]))
                throw new ArgumentException($"Unexpected token '{flag}' in the command line.", nameof(commandLine));
            var values = new List<string>();
            while (i < tokens.Length && !(tokens[i].Length > 1 && tokens[i][0] == '-' && char.IsLetter(tokens[i][1]))) values.Add(tokens[i++]);
            if (flag == "-weight")
            {
                string w = StringValue(flag, values);
                if (w == "no") weightNo = true;
                else if (w != "default") throw new ArgumentException($"-weight {w} is not supported (library weights of the methods).");
                continue;
            }

            options = options.Apply(flag, values);
        }

        // `-weight no` skips the sequence weighting altogether
        return weightNo ? options with { TCoffeeSequenceWeights = false, SequenceWeights = null } : options;
    }

    private static int IntValue(string flag, List<string> values)
    {
        if (values.Count != 1 || !int.TryParse(values[0], NumberStyles.Integer, CultureInfo.InvariantCulture, out int v))
            throw new ArgumentException($"The flag {flag} needs one integer value.");
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
                                    "GOP" => mp with { GapOpen = int.Parse(parts[a + 1], CultureInfo.InvariantCulture) },
                                    "GEP" => mp with { GapExtend = int.Parse(parts[a + 1], CultureInfo.InvariantCulture) },
                                    "MATRIX" => mp with { Matrix = SequenceAligner.TcResolveMatrix(parts[a + 1]) is null ? null : parts[a + 1] },
                                    _ => throw new ArgumentException($"-method {spec}: the method parameter {parts[a]} is not supported (GOP, GEP, MATRIX)."),
                                };
                            parameters.TryAdd(method, mp);
                        }
                    }

                    if (methods.Count == 0) throw new ArgumentException("-method needs at least one method.");
                    return this with { Methods = methods, MethodParameters = parameters.Count > 0 ? parameters : null };
                }
            case "-extend_mode": // checked: no effect on the default dp_mode linked_pair_wise (relaxation + sparse DP), 0 / 120 sets differ
            case "-extend":
            case "-evaluate_mode":
            case "-method_evaluate_mode":
            case "-shrink_lib": // only recounts the entries (shrink_constraint_list sets CL->ne), 0 / 120 sets differ
            case "-tg_mode": // the linked DP has its own gap rule (terminal gaps free of the opening penalty), 0 / 120 sets differ
            case "-ktuple":
            case "-ndiag":
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
                    var m => throw new ArgumentException($"-tree_mode {m} is not supported (nj, upgma)."),
                };
            case "-distance_matrix_mode":
                return StringValue(flag, values) switch
                {
                    "ktup" => this with { DistanceMatrixMode = TCoffeeDistanceMode.Ktup },
                    "fast" => this with { DistanceMatrixMode = TCoffeeDistanceMode.Fast },
                    "idscore" => this with { DistanceMatrixMode = TCoffeeDistanceMode.Idscore },
                    "no" => this with { DistanceMatrixMode = TCoffeeDistanceMode.No },
                    var m => throw new ArgumentException($"-distance_matrix_mode {m} is not supported."),
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
                        var w = v[1..].Split(',').Select(x => float.Parse(x, NumberStyles.Float, CultureInfo.InvariantCulture)).ToArray();
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
        _ => throw new ArgumentException($"The t_coffee method {name} is not supported."),
    };
}
