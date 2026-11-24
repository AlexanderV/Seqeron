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
        internal static int Iterations(MafftOptions o)
        {
            if (o.Retree == 0)
                return 0;
            int limit = o.ParallelizationStrategy is "BESTFIRST" or "BAATARI0" ? 254 : 16;
            return Math.Min(o.MaxIterate, limit);
        }

        /// <summary>The command line of <c>setcore</c> (<c>--core</c>): <c>-w corewin -i corethr [-c]</c>.</summary>
        internal static List<string> SetcoreArguments(MafftOptions o)
        {
            var args = new List<string> { "-w", o.CoreWindow, "-i", o.CoreThreshold };
            if (o.CoreExtend)
                args.Add("-c");
            return args;
        }

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
            nadd = Math.Max(nadd, o.AddCount); // --add / --addprofile with --maxiterate: dvtditr -K n
            string gop = Text(o.GapOpen);
            string aof = AwkPrint(0.0 + 0.0 + MafftGlobals.CAtof(AwkPrint(-1.0 * o.Offset)));
            if (o.Offset == 0.0)
                aof = "0";
            if (o.PairLibrary == MafftPairLibrary.Genaf)
                aof = "0.000"; // E-INS: the offset is 0 whatever --ep says
            int threads = o.Threads < 0 ? Environment.ProcessorCount : o.Threads;
            int threadsIt = o.ThreadsIterate ?? (threads < 8 ? threads : 8);
            var args = new List<string> { "-W", Text(o.MinimumWeight) };
            if (o.NoSegmentation)
                args.Add("-B");
            args.AddRange(new[] { "-E", Text(o.FixThreshold), "-s", Text(o.UnalignLevel) });
            if (o.LegacyGapPenalty)
                args.Add("-L");
            args.AddRange(MergeArguments(o));
            if (o.OutNumber)
                args.Add("-n");
            args.Add("-C");
            args.Add(threadsIt.ToString(CultureInfo.InvariantCulture));
            args.Add("-t");
            args.Add(o.RandomSeed.ToString(CultureInfo.InvariantCulture));
            args.Add("-F");
            args.AddRange(MemoryArguments(o));
            if (o.NoScore)
                args.Add("-Z");
            if (o.PairLibrary != MafftPairLibrary.None || o.SeedTable != null)
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
            args.Add(Text(o.ShiftPenalty));
            args.Add("-h");
            args.Add(aof);
            args.Add("-I");
            args.Add(Iterations(o).ToString(CultureInfo.InvariantCulture));
            if (o.Unweight)
                args.Add("-u");
            if (o.GuideTree != MafftGuideTree.Default)
                args.Add("-U");
            args.Add("-X");
            args.Add(o.TreeMixing is { } mix ? Text(mix) : "0.1");
            args.Add("-p");
            args.Add(o.ParallelizationStrategy);
            if (o.ScoreOut)
                args.AddRange(new[] { "-S", "-B" });
            args.Add("-K");
            args.Add(nadd.ToString(CultureInfo.InvariantCulture));
            return args;
        }

        /// <summary>The command of the later passes of a seed run (<c>tbfast -W $minimumweight -V -$gopdist -s $unalignlevel $legacygapopt $mergearg $termgapopt $outnum -C $numthreadstb $weightopt $treeoutopt $distoutopt $memopt $seqtype $model -f -$gop -Q 100.0 -h $aof $param_fft -l $weighti -J $treealg $scoreoutarg</c>).</summary>
        internal static List<string> TbfastJArguments(MafftOptions o)
        {
            int threads = o.Threads < 0 ? Environment.ProcessorCount : o.Threads;
            int threadsTb = o.ThreadsProgressive ?? (threads < 16 ? threads : 16);
            string gop = Text(o.GapOpen);
            string gopdist = o.GapOpenDistance is { } gd ? Text(gd) : gop;
            string aof = AwkPrint(0.0 + 0.0 + MafftGlobals.CAtof(AwkPrint(-1.0 * o.Offset)));
            if (o.Offset == 0.0)
                aof = "0";
            var args = new List<string> { "-W", Text(o.MinimumWeight), "-V", "-" + gopdist, "-s", Text(o.UnalignLevel) };
            if (o.LegacyGapPenalty)
                args.Add("-L");
            args.AddRange(MergeArguments(o));
            if (!o.TerminalGapPenalty)
                args.Add("-O");
            if (o.OutNumber)
                args.Add("-n");
            args.Add("-C");
            args.Add(threadsTb.ToString(CultureInfo.InvariantCulture));
            if (o.Unweight)
                args.Add("-u");
            args.AddRange(TreeOutArguments(o, 1));
            args.AddRange(MemoryArguments(o));
            if (o.SequenceType == MafftSequenceType.Nucleotide)
                args.Add("-D");
            else if (o.SequenceType == MafftSequenceType.Amino)
                args.Add("-P");
            args.AddRange(ModelArguments(o));
            args.AddRange(new[] { "-f", "-" + gop, "-Q", Text(o.ShiftPenalty), "-h", aof });
            if (o.ForceFft)
                args.Add("-G");
            else if (o.Fft)
                args.Add("-F");
            args.AddRange(new[] { "-l", Text(o.WeightI), "-J", "-X", o.TreeMixing is { } mix ? Text(mix) : "0.1" });
            if (o.ScoreOut)
                args.AddRange(new[] { "-S", "-B" });
            return args;
        }

        /// <summary><c>$memopt</c>: <c>--memsave</c> = <c>-M -B</c>, <c>--nomemsave</c> = <c>-N</c>.</summary>
        internal static IEnumerable<string> MemoryArguments(MafftOptions o) => o.MemoryMode switch
        {
            MafftMemoryMode.MemSave => new[] { "-M", "-B" },
            MafftMemoryMode.NoMemSave => new[] { "-N" },
            _ => Array.Empty<string>(),
        };

        /// <summary>The number of passes of <c>splittbfast</c> (<c>cycle</c>) and the group size the script gives it, after the rules of <c>mafft.tmpl</c> for <c>--retree 0</c> (with <c>--treeout</c>: no alignment, <c>-s 1</c>), two sequences, at most 3 passes, <c>--dpparttree</c> / <c>--fastaparttree</c> / <c>--groupsize</c> (one pass).</summary>
        internal static (int Cycle, int GroupSize) PartTreeCycle(MafftOptions o, int nseq)
        {
            int cycle = o.Retree;
            int groupsize = o.GroupSize;
            if (cycle == 0 && (o.TreeOut || o.NodeOut))
                groupsize = 1;
            if (nseq == 2)
                cycle = 1;
            if (cycle > 3)
                cycle = 3;
            if (o.PartTree is MafftPartTree.DynamicProgramming or MafftPartTree.Fasta)
                cycle = 1;
            if (groupsize > -1)
                cycle = 1;
            return (cycle, groupsize);
        }

        /// <summary><c>splittbfast $legacygapopt [-Z] $splitopt $partorderopt $parttreeoutopt $seqtype $model -f -gop -Q 100.0 -h aof -p partsize -s groupsize $treealg -i infile</c>.</summary>
        internal static List<string> SplittbfastArguments(MafftOptions o, int nseq, bool fromAlignment, int cycle, int groupsize)
        {
            _ = nseq;
            _ = cycle;
            string aof = AwkPrint(0.0 + 0.0 + MafftGlobals.CAtof(AwkPrint(-1.0 * o.Offset)));
            if (o.Offset == 0.0)
                aof = "0";
            var args = new List<string>();
            if (o.LegacyGapPenalty)
                args.Add("-L");
            if (fromAlignment)
                args.Add("-Z");
            if (o.PartTree == MafftPartTree.DynamicProgramming)
                args.Add("-U");
            else if (o.PartTree == MafftPartTree.Fasta)
                args.Add("-S");
            if (!o.Reorder)
                args.Add("-x");
            if (o.TreeOut || o.NodeOut)
                args.Add("-t");
            if (o.SequenceType == MafftSequenceType.Nucleotide)
                args.Add("-D");
            else if (o.SequenceType == MafftSequenceType.Amino)
                args.Add("-P");
            args.AddRange(ModelArguments(o));
            args.AddRange(new[] { "-f", "-" + Text(o.GapOpen), "-Q", Text(o.ShiftPenalty), "-h", aof, "-p", o.PartSize.ToString(CultureInfo.InvariantCulture), "-s", groupsize.ToString(CultureInfo.InvariantCulture),
                "-X", o.TreeMixing is { } mix ? Text(mix) : "0.1", "-i", "infile" });
            return args;
        }

        /// <summary>
        /// The name of the strategy and its description that the script prints at the end of the run (<c>strategy</c>, <c>explanation</c>, <c>performance</c> of <c>mafft.tmpl</c>), from the variables of the script.
        /// </summary>
        internal static string StrategyName(string treeext, string distance, bool sw, bool fft, string weighti, bool memsavetree, int iterate, string partdist, int fragment, int cycle, out string explanation, out string performance)
        {
            string strategy;
            if (treeext == "pileup")
                strategy = "Pileup-";
            else if (treeext == "randomchain")
                strategy = "Randomchain-";
            else if (distance == "fasta" && !sw)
                strategy = "F-";
            else if (distance == "fasta" && sw)
                strategy = "H-";
            else if (distance == "blast")
                strategy = "B-";
            else if (distance is "global" or "distonly")
                strategy = "G-";
            else if (distance == "local")
                strategy = "L-";
            else if (distance == "last")
                strategy = "Last-";
            else if (distance == "hybrid")
                strategy = "Hybrid-";
            else if (distance == "multi")
                strategy = "Multi-";
            else if (distance == "lastmulti")
                strategy = "LastMulti-";
            else if (distance == "localgenaf")
                strategy = "E-";
            else if (distance == "globalgenaf")
                strategy = "K-";
            else if (fft)
                strategy = "FFT-";
            else
                strategy = "NW-";
            if (memsavetree)
                strategy += "large-";
            if (MafftGlobals.CAtof(weighti) > 0.0)
                strategy += "I";
            strategy += "NS-";
            string c = cycle.ToString(CultureInfo.InvariantCulture);
            if (iterate > 0)
                strategy += "i";
            else if (distance == "parttree")
                strategy += PartTreeName(partdist) + c;
            else if (fragment == 1)
                strategy += "fragment";
            else if (fragment == -1)
                strategy += "full";
            else if (fragment == -2)
                strategy += "long";
            else if (fragment == -3)
                strategy += "top";
            else
                strategy += c;
            string it = iterate.ToString(CultureInfo.InvariantCulture);
            explanation = "?";
            performance = "Not tested.";
            switch (strategy)
            {
                case "F-INS-i":
                    explanation = "Iterative refinement method (<" + it + ") with LOCAL pairwise alignment information";
                    performance = "Most accurate, but very slow";
                    break;
                case "L-INS-i":
                    explanation = "Iterative refinement method (<" + it + ") with LOCAL pairwise alignment information";
                    performance = "Probably most accurate, very slow";
                    break;
                case "E-INS-i":
                    explanation = "Iterative refinement method (<" + it + ") with LOCAL pairwise alignment with generalized affine gap costs (Altschul 1998)";
                    performance = "Suitable for sequences with long unalignable regions, very slow";
                    break;
                case "G-INS-i":
                    explanation = "Iterative refinement method (<" + it + ") with GLOBAL pairwise alignment information";
                    performance = "Suitable for sequences of similar lengths, very slow";
                    break;
                case "X-INS-i":
                    explanation = "RNA secondary structure information is taken into account.";
                    performance = "For short RNA sequences only, extremely slow";
                    break;
                case "F-INS-1":
                case "L-INS-1":
                    explanation = "Progressive method incorporating LOCAL pairwise alignment information";
                    break;
                case "G-INS-1":
                    explanation = "Progressive method incorporating GLOBAL pairwise alignment information";
                    break;
                case "FFT-NS-i":
                case "NW-NS-i":
                    explanation = "Iterative refinement method (max. " + it + " iterations)";
                    performance = iterate > 2 ? "Accurate but slow" : "Standard";
                    break;
                case "FFT-NS-2":
                case "NW-NS-2":
                    explanation = "Progressive method (guide trees were built " + c + " times.)";
                    performance = "Fast but rough";
                    break;
                case "FFT-NS-1":
                case "NW-NS-1":
                    explanation = "Progressive method (rough guide tree was used.)";
                    performance = "Very fast but very rough";
                    break;
            }
            return strategy;
        }

        private static string PartTreeName(string partdist)
        {
            if (partdist == "fasta")
                return "FastaPartTree-";
            return partdist == "localalign" ? "DPPartTree-" : "PartTree-";
        }

        /// <summary>The name of the strategy and its description for the library options (<see cref="MafftOptions"/>).</summary>
        internal static (string Strategy, string Explanation, string Performance) Strategy(MafftOptions o, int nseq)
        {
            int cycle = o.Retree;
            if (nseq == 2)
                cycle = 1;
            if (cycle > 3)
                cycle = 3;
            string distance = o.PartTree != MafftPartTree.None ? "parttree" : o.PairLibrary switch { MafftPairLibrary.Global => "global", MafftPairLibrary.Local => "local", MafftPairLibrary.Genaf => "localgenaf", _ => "ktuples" };
            string partdist = o.PartTree switch { MafftPartTree.Fasta => "fasta", MafftPartTree.DynamicProgramming => "localalign", _ => "ktuples" };
            if (o.PartTree != MafftPartTree.None)
                cycle = PartTreeCycle(o, nseq).Cycle;
            else if (o.PairLibrary != MafftPairLibrary.None)
                cycle = 1; // the script sets the number of tree builds to 1 when a pair library is used
            string weighti = o.PairLibrary != MafftPairLibrary.None ? Text(o.WeightI) : "0.0";
            string strategy = StrategyName("none", distance, false, o.Fft, weighti, false, Iterations(o), partdist, 0, cycle, out string explanation, out string performance);
            return (strategy, explanation, performance);
        }

        /// <summary>The text that the script <c>mafft</c> writes to its progress stream after the binaries have run (the strategy and the notes).</summary>
        internal static string Epilogue(MafftOptions o, int nseq)
        {
            var (strategy, explanation, performance) = Strategy(o, nseq);
            string before = o.PartTree == MafftPartTree.Fasta ? FastaNote : "";
            return EpilogueText(strategy, explanation, performance, before, GenafNote(o));
        }

        internal const string FastaNote = "Pairwise alignments were computed by FASTA\n(Pearson & Lipman, 1988, PNAS 85:2444-2448)\n";

        /// <summary>The closing text of the script: <c>echo ''</c>, the lines about the pairwise programs (<paramref name="before"/>), the strategy, the notes and the E-INS-i note (<paramref name="genafNote"/>, with its leading blank line).</summary>
        internal static string EpilogueText(string strategy, string explanation, string performance, string before, string genafNote) =>
            "\n" + before + "Strategy:\n " + strategy + " (" + performance + ")\n " + explanation + "\n\n"
            + "If unsure which option to use, try 'mafft --auto input > output'.\n"
            + "For more information, see 'mafft --help', 'mafft --man' and the mafft page.\n\n"
            + "The default gap scoring scheme has been changed in version 7.110 (2013 Oct).\n"
            + "It tends to insert more gaps into gap-rich regions than previous versions.\n"
            + "To disable this change, add the --leavegappyregion option.\n"
            + genafNote
            + "\n";

        private static string GenafNote(MafftOptions o) => GenafNoteText(o.PairLibrary == MafftPairLibrary.Genaf, o.OldGenafParameters);

        internal static string GenafNoteText(bool genaf, bool old)
        {
            if (!genaf)
                return "";
            return old
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
            var pair = new List<string> { "-u", Text(o.UnalignLevel), "-l", weighti, "-C", threads.ToString(CultureInfo.InvariantCulture) };
            if (o.SequenceType == MafftSequenceType.Nucleotide)
                pair.Add("-D");
            else if (o.SequenceType == MafftSequenceType.Amino)
                pair.Add("-P");
            pair.AddRange(ModelArguments(o));
            switch (o.PairLibrary)
            {
                case MafftPairLibrary.Global:
                    pair.AddRange(new[] { "-g", Text(o.GlobalPairExtension), "-f", Text(o.GlobalPairOpen), "-Q", Text(o.ShiftPenalty), "-h", AwkPrint(0.0 + 0.0 + o.GlobalPairOffset), "-A" });
                    break;
                case MafftPairLibrary.Local:
                    pair.AddRange(new[] { "-g", Text(o.LocalPairExtension), "-f", Text(o.LocalPairOpen), "-Q", Text(o.ShiftPenalty), "-h", AwkPrint(0.0 + 0.0 + o.LocalPairOffset), "-L" });
                    break;
                default:
                    // the new parameters of E-INS-i (7.243): no offset, no extension in the pair stage and the score of the pair is the plain sum of the matrix scores of its aligned columns (-Z)
                    pair.AddRange(o.OldGenafParameters
                        ? new[] { "-g", Text(o.LocalPairExtension), "-f", Text(o.LocalPairOpen), "-Q", Text(o.ShiftPenalty), "-h", AwkPrint(0.0 + 0.0 + o.LocalPairOffset), "-O", Text(o.GenafOpen), "-E", Text(o.GenafExtension), "-N" }
                        : new[] { "-g", "0.0", "-f", Text(o.LocalPairOpen), "-Q", Text(o.ShiftPenalty), "-h", AwkPrint(0.0 + 0.0 + 0.0), "-O", Text(o.GenafOpen), "-E", Text(o.GenafExtension), "-N", "-Z" });
                    break;
            }
            var tb = new List<string>();
            if (o.GuideTree is not (MafftGuideTree.MemSaveTree or MafftGuideTree.StepAdd))
                tb.AddRange(new[] { "-+", Iterations(o).ToString(CultureInfo.InvariantCulture) });
            tb.AddRange(new[] { "-W", Text(o.MinimumWeight), "-V", "-" + gopdist, "-s", Text(o.UnalignLevel) });
            if (o.LegacyGapPenalty)
                tb.Add("-L");
            if (o.PairLibrary != MafftPairLibrary.Global && !o.TerminalGapPenalty)
                tb.Add("-O");
            if (o.OutNumber)
                tb.Add("-n");
            tb.AddRange(AddArguments(o));
            tb.AddRange(MergeArguments(o));
            tb.Add("-C");
            tb.Add((o.ThreadsProgressive ?? threadsTb).ToString(CultureInfo.InvariantCulture));
            if (o.Unweight)
                tb.Add("-u");
            if (o.GuideTree != MafftGuideTree.Default)
                tb.Add("-U");
            tb.AddRange(TreeOutArguments(o, o.Retree));
            if (o.SequenceType == MafftSequenceType.Nucleotide)
                tb.Add("-D");
            else if (o.SequenceType == MafftSequenceType.Amino)
                tb.Add("-P");
            tb.AddRange(ModelArguments(o));
            tb.AddRange(new[] { "-f", "-" + gop, "-Q", Text(o.ShiftPenalty), "-h", aof });
            if (o.ForceFft)
                tb.Add("-G");
            else if (o.Fft)
                tb.Add("-F");
            tb.AddRange(new[] { "-l", weighti, "-X", o.TreeMixing is { } mix ? Text(mix) : "0.1" });
            if (o.ScoreOut)
                tb.AddRange(new[] { "-S", "-B" });
            return (pair, tb);
        }

        /// <summary><c>$lhlimit -u ... $localparam -C n $seqtype $model -g -f -Q -h -A|-L|-N $usenaivepairscore $focusarg $treeinopt $treeoutopt -i infile</c>: the command of <c>nodepair</c> (the arguments of the pair stage of <c>tbfast</c> with the tree options).</summary>
        internal static List<string> NodepairArguments(MafftOptions o)
        {
            var (pair, _) = TbfastArguments(o, 0);
            var args = new List<string>();
            if (o.LhLimit != int.MaxValue)
            {
                args.Add("-q");
                args.Add(o.LhLimit.ToString(CultureInfo.InvariantCulture));
            }
            args.AddRange(pair);
            if (o.GuideTreeText != null)
                args.Add("-U");
            args.AddRange(TreeOutArguments(o, o.Retree));
            args.Add("-i");
            args.Add("infile");
            return args;
        }

        /// <summary>The file <c>_guidetree</c> that the script writes for the method of <see cref="MafftOptions.GuideTree"/> (<c>--memsavetreex</c> with the script's initial RAM usage of 20GB).</summary>
        internal static string GuideTreeFile(MafftOptions o) => o.GuideTree switch
        {
            MafftGuideTree.Pileup => "pileup\n",
            MafftGuideTree.RandomChain => "shuffle " + o.RandomSeed.ToString(CultureInfo.InvariantCulture) + "\n",
            MafftGuideTree.MemSaveTree => "very compact\n",
            MafftGuideTree.MemSaveTreeX => "compact  20GB\n",
            MafftGuideTree.StepAdd => "stepadd\n",
            MafftGuideTree.YoungestLinkage => "youngestlinkage\n",
            MafftGuideTree.UserTree => NormalizeScriptFile(o.GuideTreeText ?? ""),
            _ => "",
        };

        /// <summary><c>cat file | tr "\r" "\n" | grep -v "^$"</c>: carriage returns are line ends and empty lines are dropped.</summary>
        internal static string NormalizeScriptFile(string text)
        {
            var sb = new System.Text.StringBuilder();
            foreach (string line in text.Replace('\r', '\n').Split('\n'))
            {
                if (line.Length == 0)
                    continue;
                sb.Append(line).Append('\n');
            }
            return sb.ToString();
        }

        /// <summary>The text of a number as the user wrote it on the command line (the script passes it on unchanged).</summary>
        internal static string Text(double v) => v.ToString("R", CultureInfo.InvariantCulture);

        /// <summary><c>$treeoutopt $distoutopt</c> of the script: <c>--retree 0</c> computes only the guide tree (<c>-t -T</c> / <c>-^ -T</c>) or the distances (<c>-y -T</c>), otherwise the tree and the distances are written along with the alignment.</summary>
        internal static List<string> TreeOutArguments(MafftOptions o, int cycle)
        {
            var args = new List<string>();
            bool tree = o.TreeOut || o.NodeOut;
            if (cycle == 0)
            {
                if (tree)
                    args.AddRange(o.NodeOut ? new[] { "-^", "-T" } : new[] { "-t", "-T" });
                if (o.DistOut)
                    args.AddRange(new[] { "-y", "-T" });
                else if (!tree)
                    args.AddRange(new[] { "-t", "-T" });
            }
            else
            {
                if (o.NodeOut)
                    args.Add("-^");
                else if (o.TreeOut)
                    args.Add("-t");
                if (o.DistOut)
                    args.Add("-y");
            }
            return args;
        }

        /// <summary><c>$mergearg</c> of the script: <c>-H 0</c> with <c>--merge</c> (the script sets it before the number of seed sequences is known).</summary>
        private static IEnumerable<string> MergeArguments(MafftOptions o)
        {
            if (o.MergeTable == null)
                yield break;
            yield return "-H";
            yield return "0";
        }

        /// <summary><c>$addarg $add2ndhalfarg</c> of the script: <c>-K -I n</c> (<c>--add</c>) or <c>-I n</c> (<c>--addprofile</c>), <c>-p</c> (<c>--smoothing</c>), <c>-Y</c> (<c>--keeplength</c>), <c>-Z -Y</c> (<c>--mapout</c>), <c>-z -Y</c> (<c>--compactmapout</c>).</summary>
        private static IEnumerable<string> AddArguments(MafftOptions o)
        {
            if (o.AddCount <= 0)
                yield break;
            if (!o.AddProfile)
                yield return "-K";
            yield return "-I";
            yield return o.AddCount.ToString(CultureInfo.InvariantCulture);
            if (o.Smoothing)
                yield return "-p";
            if (o.KeepLength && o.MapOut == 0)
                yield return "-Y";
            if (o.MapOut == 1)
            {
                yield return "-Z";
                yield return "-Y";
            }
            else if (o.MapOut == 2)
            {
                yield return "-z";
                yield return "-Y";
            }
        }

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
            if (o.SeedTable != null)
                cycledisttbfast = 1; // a seed: the other passes are made by tbfast -J
            string gop = Text(o.GapOpen);
            string gopdist = o.GapOpenDistance is { } gd ? Text(gd) : gop;
            string aof = AwkPrint(0.0 + 0.0 + MafftGlobals.CAtof(AwkPrint(-1.0 * o.Offset)));
            if (o.Offset == 0.0)
                aof = "0";
            string gexp = AwkPrint(-1.0 * o.GapExtension);
            var args = new List<string>
            {
                "-q", o.SparsePickup.ToString(CultureInfo.InvariantCulture),
                "-E", cycledisttbfast.ToString(CultureInfo.InvariantCulture),
                "-V", "-" + gopdist,
                "-s", Text(o.UnalignLevel),
            };
            if (o.LegacyGapPenalty)
                args.Add("-L");
            args.Add("-W");
            args.Add(o.TupleSize.ToString(CultureInfo.InvariantCulture));
            if (!o.TerminalGapPenalty)
                args.Add("-O");
            int threads = o.Threads < 0 ? Environment.ProcessorCount : o.Threads;
            int threadsTb = o.ThreadsProgressive ?? (threads < 16 ? threads : 16);
            if (o.OutNumber)
                args.Add("-n");
            args.AddRange(AddArguments(o));
            args.AddRange(MergeArguments(o));
            args.Add("-C");
            args.Add(threads.ToString(CultureInfo.InvariantCulture) + "-" + threadsTb.ToString(CultureInfo.InvariantCulture));
            args.AddRange(MemoryArguments(o));
            if (o.Unweight)
                args.Add("-u");
            if (o.GuideTree != MafftGuideTree.Default)
                args.Add("-U");
            args.AddRange(TreeOutArguments(o, o.Retree));
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
            args.Add(Text(o.ShiftPenalty));
            args.Add("-h");
            args.Add(aof);
            if (o.ForceFft)
                args.Add("-G");
            else if (o.Fft)
                args.Add("-F");
            args.Add("-X");
            args.Add(o.TreeMixing is { } mix ? Text(mix) : "0.1");
            if (o.ScoreOut)
                args.AddRange(new[] { "-S", "-B" });
            args.Add("-x");
            args.Add("1000");
            if (o.OneIteration)
                args.Add("-r");
            return args;
        }
    }

    #endregion
}
