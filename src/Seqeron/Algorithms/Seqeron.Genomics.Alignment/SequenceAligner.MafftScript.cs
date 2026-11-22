using System.Globalization;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: the script mafft (option -> binary pipeline)

    internal static class MafftScript
    {
        /// <summary>awk's number to text conversion of <c>print</c> (OFMT = %.6g, integral values as integers).</summary>
        internal static string AwkPrint(double v)
        {
            if (double.IsNaN(v))
                return "nan";
            if (double.IsInfinity(v))
                return v > 0 ? "inf" : "-inf";
            if (v == Math.Floor(v) && Math.Abs(v) < 1e16)
            {
                if (v == 0.0 && double.IsNegative(v))
                    return "-0";
                return ((long)v).ToString(CultureInfo.InvariantCulture);
            }
            return FormatG6(v);
        }

        /// <summary>printf("%.6g").</summary>
        internal static string FormatG6(double v)
        {
            if (v == 0.0)
                return double.IsNegative(v) ? "-0" : "0";
            // the decimal exponent after rounding to 6 significant digits
            string e = v.ToString("E5", CultureInfo.InvariantCulture); // d.ddddde+xxx
            int ePos = e.IndexOf('E');
            int exp = int.Parse(e[(ePos + 1)..], CultureInfo.InvariantCulture);
            if (exp < -4 || exp >= 6)
            {
                string mant = e[..ePos];
                if (mant.Contains('.'))
                    mant = mant.TrimEnd('0').TrimEnd('.');
                string sign = exp < 0 ? "-" : "+";
                return mant + "e" + sign + Math.Abs(exp).ToString("00", CultureInfo.InvariantCulture);
            }
            string f = v.ToString("F" + Math.Max(0, 5 - exp).ToString(CultureInfo.InvariantCulture), CultureInfo.InvariantCulture);
            if (f.Contains('.'))
                f = f.TrimEnd('0').TrimEnd('.');
            return f;
        }

        /// <summary><c>$model</c> of the script: the scoring matrix (<c>-b</c> / <c>-j</c> / <c>-m</c>), <c>-k</c>, <c>-a</c>, <c>-:</c>.</summary>
        internal static List<string> ModelArguments(MafftOptions o)
        {
            var args = new List<string>();
            switch (o.ScoringModel)
            {
                case MafftScoringModel.Blosum:
                    args.Add("-b");
                    args.Add(o.Blosum.ToString(CultureInfo.InvariantCulture));
                    break;
                case MafftScoringModel.Jtt:
                    args.Add("-j");
                    args.Add(o.Pam.ToString(CultureInfo.InvariantCulture));
                    break;
                case MafftScoringModel.Tm:
                    args.Add("-m");
                    args.Add(o.Pam.ToString(CultureInfo.InvariantCulture));
                    break;
            }
            if (o.Kappa is { } kappa)
            {
                args.Add("-k");
                args.Add(kappa.ToString(CultureInfo.InvariantCulture));
            }
            if (o.FrequencyModel)
                args.Add("-a");
            if (o.WildcardN)
                args.Add("-:");
            return args;
        }

        /// <summary>The number of refinement iterations that the script passes to <c>dvtditr</c> (<c>--maxiterate</c>, at most 16).</summary>
        internal static int Iterations(MafftOptions o) => Math.Min(o.MaxIterate, 16);

        /// <summary><c>dndpre $seqtype $model -M 2 -C $numthreads</c>.</summary>
        internal static List<string> DndpreArguments(MafftOptions o)
        {
            var args = new List<string>();
            if (o.SequenceType == MafftSequenceType.Nucleotide)
                args.Add("-D");
            else if (o.SequenceType == MafftSequenceType.Amino)
                args.Add("-P");
            args.AddRange(ModelArguments(o));
            args.Add("-M");
            args.Add("2");
            args.Add("-C");
            args.Add((o.Threads < 0 ? Environment.ProcessorCount : o.Threads).ToString(CultureInfo.InvariantCulture));
            return args;
        }

        /// <summary><c>dvtditr -W $minimumweight -E $fixthreshold -s $unalignlevel $legacygapopt -C $numthreadsit -t $randomseed -F -z 50 $seqtype $model -f -$gop -Q $spfactor -h $aof -I $iterate $weightopt -p BAATARI2 -K $nadd</c>.</summary>
        internal static List<string> DvtditrArguments(MafftOptions o, int nadd = 0)
        {
            string gop = Text(o.GapOpen);
            string aof = AwkPrint(0.0 + 0.0 + MafftGlobals.CAtof(AwkPrint(-1.0 * o.Offset)));
            if (o.Offset == 0.0)
                aof = "0";
            if (o.PairLibrary == MafftPairLibrary.Genaf)
                aof = "0.000"; // E-INS: the offset is 0 whatever --ep says
            int threads = o.Threads < 0 ? Environment.ProcessorCount : o.Threads;
            int threadsIt = o.ThreadsIterate ?? (threads < 8 ? threads : 8);
            var args = new List<string> { "-W", "0.00001", "-E", "0.0", "-s", "0.0" };
            if (o.LegacyGapPenalty)
                args.Add("-L");
            args.Add("-C");
            args.Add(threadsIt.ToString(CultureInfo.InvariantCulture));
            args.Add("-t");
            args.Add("0");
            args.Add("-F");
            if (o.PairLibrary != MafftPairLibrary.None)
            {
                args.Add("-l");
                args.Add(Text(o.WeightI));
            }
            args.Add("-z");
            args.Add("50");
            if (o.SequenceType == MafftSequenceType.Nucleotide)
                args.Add("-D");
            else if (o.SequenceType == MafftSequenceType.Amino)
                args.Add("-P");
            args.AddRange(ModelArguments(o));
            args.Add("-f");
            args.Add("-" + gop);
            args.Add("-Q");
            args.Add("100.0");
            args.Add("-h");
            args.Add(aof);
            args.Add("-I");
            args.Add(Iterations(o).ToString(CultureInfo.InvariantCulture));
            if (o.Unweight)
                args.Add("-u");
            args.Add("-X");
            args.Add(o.TreeMixing is { } mix ? Text(mix) : "0.1");
            args.Add("-p");
            args.Add("BAATARI2");
            args.Add("-K");
            args.Add(nadd.ToString(CultureInfo.InvariantCulture));
            return args;
        }

        /// <summary>The name of the strategy and its description that the script prints at the end of the run (<c>strategy</c>, <c>explanation</c>, <c>performance</c> of <c>mafft.tmpl</c>).</summary>
        internal static (string Strategy, string Explanation, string Performance) Strategy(MafftOptions o, int nseq)
        {
            int cycle = o.Retree;
            if (nseq == 2)
                cycle = 1;
            if (cycle > 3)
                cycle = 3;
            int iterate = Iterations(o);
            string prefix = o.PairLibrary switch { MafftPairLibrary.Global => "G-", MafftPairLibrary.Local => "L-", MafftPairLibrary.Genaf => "E-", _ => o.Fft ? "FFT-" : "NW-" };
            if (o.PairLibrary != MafftPairLibrary.None)
            {
                cycle = 1; // the script sets the number of tree builds to 1 when a pair library is used
                if (o.WeightI > 0.0)
                    prefix += "I";
            }
            string strategy = prefix + "NS-" + (iterate > 0 ? "i" : cycle.ToString(CultureInfo.InvariantCulture));
            string explanation = "?";
            string performance = "Not tested.";
            if (strategy == "G-INS-i")
            {
                explanation = "Iterative refinement method (<" + iterate.ToString(CultureInfo.InvariantCulture) + ") with GLOBAL pairwise alignment information";
                performance = "Suitable for sequences of similar lengths, very slow";
            }
            else if (strategy == "L-INS-i")
            {
                explanation = "Iterative refinement method (<" + iterate.ToString(CultureInfo.InvariantCulture) + ") with LOCAL pairwise alignment information";
                performance = "Probably most accurate, very slow";
            }
            else if (strategy == "E-INS-i")
            {
                explanation = "Iterative refinement method (<" + iterate.ToString(CultureInfo.InvariantCulture) + ") with LOCAL pairwise alignment with generalized affine gap costs (Altschul 1998)";
                performance = "Suitable for sequences with long unalignable regions, very slow";
            }
            else if (strategy == "L-INS-1")
            {
                explanation = "Progressive method incorporating LOCAL pairwise alignment information";
            }
            else if (strategy == "G-INS-1")
            {
                explanation = "Progressive method incorporating GLOBAL pairwise alignment information";
            }
            else if (strategy is "FFT-NS-i" or "NW-NS-i")
            {
                explanation = "Iterative refinement method (max. " + iterate.ToString(CultureInfo.InvariantCulture) + " iterations)";
                performance = iterate > 2 ? "Accurate but slow" : "Standard";
            }
            else if (strategy is "FFT-NS-2" or "NW-NS-2")
            {
                explanation = "Progressive method (guide trees were built " + cycle.ToString(CultureInfo.InvariantCulture) + " times.)";
                performance = "Fast but rough";
            }
            else if (strategy is "FFT-NS-1" or "NW-NS-1")
            {
                explanation = "Progressive method (rough guide tree was used.)";
                performance = "Very fast but very rough";
            }
            return (strategy, explanation, performance);
        }

        /// <summary>The text that the script <c>mafft</c> writes to its progress stream after the binaries have run (the strategy and the notes).</summary>
        internal static string Epilogue(MafftOptions o, int nseq)
        {
            var (strategy, explanation, performance) = Strategy(o, nseq);
            return "\nStrategy:\n " + strategy + " (" + performance + ")\n " + explanation + "\n\n"
                + "If unsure which option to use, try 'mafft --auto input > output'.\n"
                + "For more information, see 'mafft --help', 'mafft --man' and the mafft page.\n\n"
                + "The default gap scoring scheme has been changed in version 7.110 (2013 Oct).\n"
                + "It tends to insert more gaps into gap-rich regions than previous versions.\n"
                + "To disable this change, add the --leavegappyregion option.\n"
                + GenafNote(o)
                + "\n";
        }

        private static string GenafNote(MafftOptions o)
        {
            if (o.PairLibrary != MafftPairLibrary.Genaf)
                return "";
            return o.OldGenafParameters
                ? "\nObsolete parameters used for this calculation.\nAlso try the new parameters for E-INS-i, by not specifying --oldgenafpair.\n"
                : "\nParameters for the E-INS-i option have been changed in version 7.243 (2015 Jun).\nTo switch to the old parameters, use --oldgenafpair, instead of --genafpair.\n";
        }

        /// <summary>
        /// The command <c>tbfast _ (pair library arguments) _ -+ n (arguments of the progressive stage)</c> that the script builds for <c>--globalpair</c> / <c>--localpair</c> / <c>--genafpair</c>: the two argument lists.
        /// </summary>
        internal static (List<string> Pair, List<string> Tb) TbfastArguments(MafftOptions o, int nseq)
        {
            _ = nseq;
            int threads = o.Threads < 0 ? Environment.ProcessorCount : o.Threads;
            int threadsTb = threads < 16 ? threads : 16;
            string gop = Text(o.GapOpen);
            string gopdist = o.GapOpenDistance is { } gd ? Text(gd) : gop;
            string aof = AwkPrint(0.0 + 0.0 + MafftGlobals.CAtof(AwkPrint(-1.0 * o.Offset)));
            if (o.Offset == 0.0)
                aof = "0";
            string weighti = Text(o.WeightI);
            bool genaf = o.PairLibrary == MafftPairLibrary.Genaf;
            if (genaf)
                aof = "0.000"; // E-INS: the offset of the progressive stage is 0 whatever --ep says
            var pair = new List<string> { "-u", "0.0", "-l", weighti, "-C", threads.ToString(CultureInfo.InvariantCulture) };
            if (o.SequenceType == MafftSequenceType.Nucleotide)
                pair.Add("-D");
            else if (o.SequenceType == MafftSequenceType.Amino)
                pair.Add("-P");
            pair.AddRange(ModelArguments(o));
            switch (o.PairLibrary)
            {
                case MafftPairLibrary.Global:
                    pair.AddRange(new[] { "-g", Text(o.GlobalPairExtension), "-f", Text(o.GlobalPairOpen), "-Q", "100.0", "-h", AwkPrint(0.0 + 0.0 + o.GlobalPairOffset), "-A" });
                    break;
                case MafftPairLibrary.Local:
                    pair.AddRange(new[] { "-g", Text(o.LocalPairExtension), "-f", Text(o.LocalPairOpen), "-Q", "100.0", "-h", AwkPrint(0.0 + 0.0 + o.LocalPairOffset), "-L" });
                    break;
                default:
                    // the new parameters of E-INS-i (7.243): no offset, no extension in the pair stage and the score of the pair is the plain sum of the matrix scores of its aligned columns (-Z)
                    pair.AddRange(o.OldGenafParameters
                        ? new[] { "-g", Text(o.LocalPairExtension), "-f", Text(o.LocalPairOpen), "-Q", "100.0", "-h", AwkPrint(0.0 + 0.0 + o.LocalPairOffset), "-O", Text(o.GenafOpen), "-E", Text(o.GenafExtension), "-N" }
                        : new[] { "-g", "0.0", "-f", Text(o.LocalPairOpen), "-Q", "100.0", "-h", AwkPrint(0.0 + 0.0 + 0.0), "-O", Text(o.GenafOpen), "-E", Text(o.GenafExtension), "-N", "-Z" });
                    break;
            }
            var tb = new List<string> { "-+", Iterations(o).ToString(CultureInfo.InvariantCulture), "-W", "0.00001", "-V", "-" + gopdist, "-s", "0.0" };
            if (o.LegacyGapPenalty)
                tb.Add("-L");
            if (o.PairLibrary != MafftPairLibrary.Global && !o.TerminalGapPenalty)
                tb.Add("-O");
            tb.Add("-C");
            tb.Add(threadsTb.ToString(CultureInfo.InvariantCulture));
            if (o.Unweight)
                tb.Add("-u");
            if (o.SequenceType == MafftSequenceType.Nucleotide)
                tb.Add("-D");
            else if (o.SequenceType == MafftSequenceType.Amino)
                tb.Add("-P");
            tb.AddRange(ModelArguments(o));
            tb.AddRange(new[] { "-f", "-" + gop, "-Q", "100.0", "-h", aof });
            if (o.ForceFft)
                tb.Add("-G");
            else if (o.Fft)
                tb.Add("-F");
            tb.AddRange(new[] { "-l", weighti, "-X", o.TreeMixing is { } mix ? Text(mix) : "0.1" });
            return (pair, tb);
        }

        /// <summary>The text of a number as the user wrote it on the command line (the script passes it on unchanged).</summary>
        internal static string Text(double v) => v.ToString("R", CultureInfo.InvariantCulture);

        /// <summary>
        /// The arguments that <c>mafft</c> (7.505, <c>core/mafft.tmpl</c>) passes to <c>disttbfast</c> for the plain progressive methods (the default FFT-NS-2, <c>--retree 1</c>, <c>--nofft</c>, <c>--6merpair</c>, ...).
        /// <paramref name="nseq"/> is the number of sequences (the script sets the number of tree builds to 1 for two sequences).
        /// </summary>
        internal static List<string> DisttbfastArguments(MafftOptions o, int nseq)
        {
            int cycle = o.Retree;
            if (nseq == 2)
                cycle = 1;
            if (cycle > 3)
                cycle = 3;
            int cycledisttbfast = cycle == 0 ? 1 : cycle;
            string gop = Text(o.GapOpen);
            string gopdist = o.GapOpenDistance is { } gd ? Text(gd) : gop;
            string aof = AwkPrint(0.0 + 0.0 + MafftGlobals.CAtof(AwkPrint(-1.0 * o.Offset)));
            if (o.Offset == 0.0)
                aof = "0";
            string gexp = AwkPrint(-1.0 * o.GapExtension);
            var args = new List<string>
            {
                "-q", "0",
                "-E", cycledisttbfast.ToString(CultureInfo.InvariantCulture),
                "-V", "-" + gopdist,
                "-s", "0.0",
            };
            if (o.LegacyGapPenalty)
                args.Add("-L");
            args.Add("-W");
            args.Add(o.TupleSize.ToString(CultureInfo.InvariantCulture));
            if (!o.TerminalGapPenalty)
                args.Add("-O");
            int threads = o.Threads < 0 ? Environment.ProcessorCount : o.Threads;
            int threadsTb = threads < 16 ? threads : 16;
            args.Add("-C");
            args.Add(threads.ToString(CultureInfo.InvariantCulture) + "-" + threadsTb.ToString(CultureInfo.InvariantCulture));
            if (o.Unweight)
                args.Add("-u");
            if (o.SequenceType == MafftSequenceType.Nucleotide)
                args.Add("-D");
            else if (o.SequenceType == MafftSequenceType.Amino)
                args.Add("-P");
            args.AddRange(ModelArguments(o));
            args.Add("-g");
            args.Add(gexp);
            args.Add("-f");
            args.Add("-" + gop);
            args.Add("-Q");
            args.Add("100.0");
            args.Add("-h");
            args.Add(aof);
            if (o.ForceFft)
                args.Add("-G");
            else if (o.Fft)
                args.Add("-F");
            args.Add("-X");
            args.Add(o.TreeMixing is { } mix ? Text(mix) : "0.1");
            args.Add("-x");
            args.Add("1000");
            return args;
        }
    }

    #endregion
}
