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

        /// <summary>The name of the strategy and its description that the script prints at the end of the run (<c>strategy</c>, <c>explanation</c>, <c>performance</c> of <c>mafft.tmpl</c>).</summary>
        internal static (string Strategy, string Explanation, string Performance) Strategy(MafftOptions o, int nseq)
        {
            int cycle = o.Retree;
            if (nseq == 2)
                cycle = 1;
            if (cycle > 3)
                cycle = 3;
            string strategy = (o.Fft ? "FFT-" : "NW-") + "NS-" + cycle.ToString(CultureInfo.InvariantCulture);
            string explanation = "?";
            string performance = "Not tested.";
            if (strategy is "FFT-NS-2" or "NW-NS-2")
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
                + "To disable this change, add the --leavegappyregion option.\n\n";
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
