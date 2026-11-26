#pragma warning disable S3626, S3871, S1643, S2325, IDE0052, IDE0051, IDE0059, S127, S3358, S1066, S4487, S1144, CA1823, CS0414, CA1870, S1125, S3776, S107, S138, S1541, S1854
using System.Globalization;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    internal sealed partial class MafftScriptRun
    {
        /// <summary>
        /// The RNA options in an installation that has the extensions (<see cref="MafftInstallation.Mxscarna"/>), as far as they stop before the base-pairing probabilities are needed (A2-19c-rna, run 2 of the Docker package, MAFFT 7.505 built with its extensions):
        /// <list type="bullet">
        /// <item><c>--mccaskill</c> with protein input: <c>mccaskillwrap -s -C n -d prefix -i infile &gt; hat4</c> prints its arguments, <c>getnumlen</c> finds no nucleotides and it says "nuc only" (exit 1; <c>mccaskillwrap.c</c>, <c>main</c> after <c>arguments</c>).</item>
        /// <item><c>--scarnapair</c> without <c>--mccaskill</c> (<c>--dafs</c>, <c>--contrafold</c>): the script calls <c>tbfast _ ... -s -d prefix _ -+ n ...</c> directly; the progressive program prints its arguments, <c>treein</c> / <c>compacttree</c>, the stack size, <c>whereispairalign</c>,
        /// the scoring matrix, "All-to-all alignment.", <c>pairlocalalign</c> (alg 's') prepares the base-pairing probabilities with "Preparing bpp" and finds no file <c>hat4</c> (<c>preparebpp</c>, <c>pairlocalalign.c</c>): "Cannot open hat4", exit 1 - the file is written only by <c>mccaskillwrap</c> /
        /// <c>contrafoldwrap</c> of the three options that set it.</item>
        /// </list>
        /// Everything else (the base-pairing probabilities of <c>mxscarnamod -m -writebpp</c>, the profile scores of <c>RNAalifold</c>) is not ported: <see cref="NeedPort"/>.
        /// </summary>
        /// <returns>The exit status of the stop, or null when the run does not stop here.</returns>
        private int? RnaExtensionStop(List<MafftProgramIo.Record> recs, int nseq)
        {
            if (!_installation.Mxscarna || _nadd > 0 || _seed != "x" || _seedtable != "x" || _fragment != 0)
                return null;
            bool scarnaDirect = _distance == "scarna" && !_mccaskill && !_dafs && !_contrafold;
            bool dafspairDirect = _distance == "dafs" && !_mccaskill && !_dafs && !_contrafold;
            if (!_mccaskill && !_dafs && !scarnaDirect && !dafspairDirect)
                return null;
            var originals = recs.Select(r => MafftEngine.PlainOrAnyResidues(r.Raw, _anysymbol)).ToList();
            bool nuc = MafftEngine.GuessNucleotide(originals, new MafftOptions { SequenceType = _seqtype == "-D" ? MafftSequenceType.Nucleotide : _seqtype == "-P" ? MafftSequenceType.Amino : MafftSequenceType.Auto });
            if (_mccaskill || _dafs)
            {
                // mccaskillwrap -s|-G -C n -d prefix -i infile > hat4: the arguments are printed in the order they are read
                string args = "nthread = " + _numthreads.ToString(CultureInfo.InvariantCulture) + "\nwhereismccaskillmea = " + _installation.LibraryDirectory + "\ninputfile = infile\n";
                if (!nuc)
                {
                    if (_anysymbol)
                        Prog(MafftEngine.ReplaceuLog(originals, nuc));
                    Err(args + "nuc only\n");
                    return 1;
                }
                if (!_dafs || _installation.DafsMafftInterface)
                    return null;
                // --dafs: every sequence is cut from its gaps and written to _mccaskillin; `dafs --mafft-out _mccaskillout _mccaskillin` is rejected by the DAFS that is installed (see MafftInstallation.DafsMafftInterface), so the first sequence
                // that has a residue ends the run ("1 / n" for the program without threads; the threads of mccaskillwrap race for the jobs, the order of their messages is not defined: not ported)
                if (_numthreads > 1)
                    return null;
                var g0 = MafftGlobals.Parse(Array.Empty<string>());
                g0.Dorp = 'd';
                g0.Constants(Array.Empty<string>());
                var sb = new System.Text.StringBuilder(args).Append(g0.Log).Append("Running DAFS (Sato et al. 2012; http://www.ncrna.org/).\n");
                int njob = originals.Count;
                for (int i = 0; i < njob; i++)
                {
                    bool hasResidue = originals[i].Any(c => c != '-');
                    if (_numthreads > 0)
                    {
                        if (!hasResidue)
                            continue; // athread: an empty sequence is skipped before its message
                        sb.Append(CultureInfo.InvariantCulture, $"{i + 1} / {njob} (by thread {0,4})\n");
                    }
                    else
                    {
                        sb.Append(CultureInfo.InvariantCulture, $"{i + 1} / {njob}\n");
                    }
                    if (hasResidue)
                    {
                        if (_anysymbol)
                            Prog(MafftEngine.ReplaceuLog(originals, nuc));
                        Err(sb.Append("ERROR IN mccaskill_mea\n").ToString());
                        return 1;
                    }
                }
                return null;
            }
            if (_anysymbol)
                Prog(MafftEngine.ReplaceuLog(originals, nuc));
            var o = BuildOptions(recs.Select(r => r.Name).ToList(), nuc) with { PairLibrary = MafftPairLibrary.Local };
            var (_, tb) = MafftScript.TbfastArguments(o, nseq);
            var g = MafftGlobals.Parse(MafftScript.ModelArguments(o));
            g.Dorp = nuc ? 'd' : 'p';
            g.Constants(Array.Empty<string>());
            if (dafspairDirect)
            {
                // pairlocalalign ... -G -d prefix: the same missing hat4, without the arguments of tbfast
                Err("whereispairalign = " + _installation.LibraryDirectory + "\n" + g.Log + "All-to-all alignment.\nPreparing bpp\nCannot open hat4\n");
                return 1;
            }
            Err(MafftGlobals.ParseTbfast(tb).Log + "treein = " + (o.GuideTree != MafftGuideTree.Default ? "1" : "0") + "\ncompacttree = 0\n"
                + MafftRun.StackSizeMessage(8192, nseq) + "whereispairalign = " + _installation.LibraryDirectory + "\n" + g.Log + "All-to-all alignment.\nPreparing bpp\nCannot open hat4\n");
            return 1;
        }

        /// <summary>The conditions under which the port runs the CONTRAfold branch: the installation has <c>contrafoldwrap</c> and <c>contrafold</c>, a pair library (global, local or generalized affine) and the progressive method (no iterative refinement, no <c>--add</c>, <c>--seed</c>, fragments).</summary>
        private bool ContrafoldSupported() =>
            _installation.Contrafold && !_mccaskill && !_dafs && !_rnaalifold && _nadd == 0 && _seed == "x" && _seedtable == "x" && _fragment == 0 && _iterate == 0
            && _distance is "global" or "local" or "localgenaf" && _mergetable == "/dev/null";

        /// <summary>
        /// <c>contrafoldwrap -d prefix -i infile &gt; hat4</c> (<c>contrafoldwrap.c</c>): for every sequence (letters other than A, C, G, U / T are N, gaps are removed) it runs <c>contrafold predict _contrafoldin --posteriors 0.01 _contrafoldout</c>
        /// (<see cref="ContrafoldProgram"/>), reads the file as <c>readcontrafold</c> does (the last line is read twice at the end of the file) and writes the pairs <c>i j %f</c> (i &lt; j, zero-based positions) to <c>hat4</c>.
        /// </summary>
        /// <returns>The pairing (hat4 and the messages), or <c>null</c> when the program stopped (the messages are written to standard error).</returns>
        private (MafftRnaPairing? Pairing, int Rc) ContrafoldWrap(List<MafftProgramIo.Record> recs)
        {
            var originals = recs.Select(r => MafftEngine.PlainOrAnyResidues(r.Raw, _anysymbol)).ToList();
            bool nuc = MafftEngine.GuessNucleotide(originals, new MafftOptions { SequenceType = _seqtype == "-D" ? MafftSequenceType.Nucleotide : _seqtype == "-P" ? MafftSequenceType.Amino : MafftSequenceType.Auto });
            var sb = new System.Text.StringBuilder("whereiscontrafold = " + _installation.LibraryDirectory + "\ninputfile = infile\n");
            if (_anysymbol)
                Prog(MafftEngine.ReplaceuLog(originals, nuc));
            if (!nuc)
            {
                Err(sb.Append("nuc only\n").ToString());
                return (null, 1);
            }
            var g0 = MafftGlobals.Parse(Array.Empty<string>());
            g0.Dorp = 'd';
            g0.Constants(Array.Empty<string>());
            sb.Append(g0.Log).Append("running contrafold\n");
            var predictor = _installation.ContrafoldPredictor ?? ContrafoldProgram.PredictPosteriors;
            var hat4 = new System.Text.StringBuilder();
            int njob = originals.Count;
            for (int i = 0; i < njob; i++)
            {
                sb.Append(CultureInfo.InvariantCulture, $"{i + 1} / {njob}\n");
                var nogap = new System.Text.StringBuilder();
                foreach (char c in originals[i])
                {
                    if (c == '-')
                        continue;
                    nogap.Append(char.ToLowerInvariant(c) switch { 'a' => 'A', 't' or 'u' => 'U', 'g' => 'G', 'c' => 'C', _ => 'N' });
                }
                string posteriors;
                try
                {
                    posteriors = predictor(nogap.ToString(), 0.01);
                }
                catch (ArgumentException)
                {
                    Err(sb.Append("ERROR: Zero-length sequence read.\nerror in contrafold\n=================================================================\n=================================================================\n==\n== This version of MAFFT supports CONTRAfold v2.02.\n== If you have a lower version of CONTRAfold installed in the\n== " + _installation.LibraryDirectory + " directory,\n== please update it!\n==\n=================================================================\n=================================================================\n").ToString());
                    return (null, 1);
                }
                var pairprob = new List<RnaPair>[nogap.Length];
                for (int k = 0; k < pairprob.Length; k++)
                    pairprob[k] = new List<RnaPair>();
                // readcontrafold: the loop tests feof before it reads, so the last line is processed twice
                var lines = posteriors.Split('\n', StringSplitOptions.RemoveEmptyEntries);
                string last = "";
                for (int n = 0; n <= lines.Length; n++)
                {
                    string line = n < lines.Length ? lines[n] : last;
                    last = line;
                    int sp = line.IndexOf(' ');
                    if (sp < 0)
                        continue;
                    int left = MafftGlobals.CAtoi(line[..sp]) - 1;
                    int p2 = line.IndexOf(' ', sp + 1);
                    while (p2 >= 0)
                    {
                        int next = line.IndexOf(' ', p2 + 1);
                        string token = next >= 0 ? line[(p2 + 1)..next] : line[(p2 + 1)..];
                        int colon = token.IndexOf(':');
                        if (colon > 0)
                        {
                            int right = MafftGlobals.CAtoi(token[..colon]) - 1;
                            double prob = MafftGlobals.CAtof(token[(colon + 1)..]);
                            pairprob[left].Add(new RnaPair(right, prob));
                            pairprob[right].Add(new RnaPair(left, prob));
                        }
                        p2 = next;
                    }
                }
                hat4.Append(CultureInfo.InvariantCulture, $">{i}\n");
                for (int k = 0; k < pairprob.Length; k++)
                    foreach (var pt in pairprob[k])
                        if (pt.Pos > k)
                            hat4.Append(CultureInfo.InvariantCulture, $"{k} {pt.Pos} {pt.Score:F6}\n");
            }
            return (new MafftRnaPairing(hat4.ToString(), Messages: sb.ToString()), 0);
        }

        /// <summary>The conditions under which the port runs the McCaskill branch (protein input stopped before in <c>RnaExtensionStop</c>): MXSCARNA is installed and given (<see cref="MafftInstallation.MxscarnaPredictor"/>), a pair library (global, local or generalized affine), the progressive method, at most one thread (the messages of several threads are a race).</summary>
        private bool McCaskillSupported() =>
            _installation.Mxscarna && _installation.MxscarnaPredictor != null && !_dafs && !_rnaalifold && !_contrafold && _nadd == 0 && _seed == "x" && _seedtable == "x" && _fragment == 0 && _iterate == 0
            && _distance is "global" or "local" or "localgenaf" && _mergetable == "/dev/null" && _numthreads <= 1;


        /// <summary>
        /// <c>mccaskillwrap -s -C n -d prefix -i infile &gt; hat4</c> (<c>mccaskillwrap.c</c>): for every sequence without its gaps the program <c>mxscarnamod -m -writebpp</c> (<see cref="MafftInstallation.MxscarnaPredictor"/>), <c>readrawmccaskill</c> (the threshold 0.01,
        /// both directions; a last line without a line end is dropped) and <c>outmccaskill</c> (<c>"%d %d %50.40f"</c>). The thread of <c>-C 1</c> numbers its messages ("i / n (by thread %4d)") and skips an empty sequence before the message.
        /// </summary>
        private (MafftRnaPairing? Pairing, int Rc) McCaskillWrap(List<MafftProgramIo.Record> recs)
        {
            var originals = recs.Select(r => MafftEngine.PlainOrAnyResidues(r.Raw, _anysymbol)).ToList();
            var sb = new System.Text.StringBuilder("nthread = " + _numthreads.ToString(CultureInfo.InvariantCulture) + "\nwhereismccaskillmea = " + _installation.LibraryDirectory + "\ninputfile = infile\n");
            if (_anysymbol)
                Prog(MafftEngine.ReplaceuLog(originals, true));
            var g0 = MafftGlobals.Parse(Array.Empty<string>());
            g0.Dorp = 'd';
            g0.Constants(Array.Empty<string>());
            sb.Append(g0.Log).Append("Running mxscarna with the mccaskill_mea mode.\n");
            var hat4 = new System.Text.StringBuilder();
            int njob = originals.Count;
            int nlenmax = originals.Max(o => o.Length);
            for (int i = 0; i < njob; i++)
            {
                string nogap = originals[i].Replace("-", "", StringComparison.Ordinal);
                if (_numthreads > 0)
                {
                    if (nogap.Length == 0)
                    {
                        hat4.Append(CultureInfo.InvariantCulture, $">{i}\n");
                        continue;
                    }
                    sb.Append(CultureInfo.InvariantCulture, $"{i + 1} / {njob} (by thread {0,4})\n");
                }
                else
                {
                    sb.Append(CultureInfo.InvariantCulture, $"{i + 1} / {njob}\n");
                    if (nogap.Length == 0)
                    {
                        hat4.Append(CultureInfo.InvariantCulture, $">{i}\n");
                        continue;
                    }
                }
                string output;
                try
                {
                    output = _installation.MxscarnaPredictor!(nogap);
                }
                catch (ArgumentException)
                {
                    Err(sb.Append("ERROR IN mccaskill_mea\n").ToString());
                    return (null, 1);
                }
                var pairprob = new List<RnaPair>[nlenmax];
                for (int k = 0; k < nlenmax; k++)
                    pairprob[k] = new List<RnaPair>();
                var lines = output.Split('\n');
                for (int n = 0; n < lines.Length - 1; n++) // a last line without a line end is read together with the end of the file and dropped
                {
                    string line = lines[n];
                    if (line.StartsWith('>'))
                        continue;
                    var f = line.Split(' ', StringSplitOptions.RemoveEmptyEntries);
                    if (f.Length < 3)
                        continue;
                    int left = MafftGlobals.CAtoi(f[0]);
                    int right = MafftGlobals.CAtoi(f[1]);
                    double prob = MafftGlobals.CAtof(f[2]);
                    if (prob < 0.01)
                        continue;
                    if (left != right && prob > 0.0)
                    {
                        if (left < 0 || right < 0 || left >= nlenmax || right >= nlenmax)
                            throw new NotSupportedException("mxscarnamod wrote a position outside the sequences: the program reads outside its table there.");
                        pairprob[left].Add(new RnaPair(right, prob));
                        pairprob[right].Add(new RnaPair(left, prob));
                    }
                }
                hat4.Append(CultureInfo.InvariantCulture, $">{i}\n");
                for (int k = 0; k < nlenmax; k++)
                    foreach (var pt in pairprob[k])
                        if (pt.Pos > k)
                            hat4.Append(CultureInfo.InvariantCulture, $"{k} {pt.Pos} {pt.Score.ToString("F40", CultureInfo.InvariantCulture),50}\n");
            }
            sb.Append(CultureInfo.InvariantCulture, $"{_numthreads} thread(s)\n");
            return (new MafftRnaPairing(hat4.ToString(), Messages: sb.ToString()), 0);
        }
    }
}
