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
            if (!_mccaskill && !scarnaDirect)
                return null;
            var originals = recs.Select(r => MafftEngine.PlainOrAnyResidues(r.Raw, _anysymbol)).ToList();
            bool nuc = MafftEngine.GuessNucleotide(originals, new MafftOptions { SequenceType = _seqtype == "-D" ? MafftSequenceType.Nucleotide : _seqtype == "-P" ? MafftSequenceType.Amino : MafftSequenceType.Auto });
            if (_mccaskill && !nuc)
            {
                if (_anysymbol)
                    Prog(MafftEngine.ReplaceuLog(originals, nuc));
                Err("nthread = " + _numthreads.ToString(CultureInfo.InvariantCulture) + "\nwhereismccaskillmea = " + _installation.LibraryDirectory + "\ninputfile = infile\nnuc only\n");
                return 1;
            }
            if (!scarnaDirect)
                return null;
            if (_anysymbol)
                Prog(MafftEngine.ReplaceuLog(originals, nuc));
            var o = BuildOptions(recs.Select(r => r.Name).ToList(), nuc) with { PairLibrary = MafftPairLibrary.Local };
            var (_, tb) = MafftScript.TbfastArguments(o, nseq);
            var g = MafftGlobals.Parse(MafftScript.ModelArguments(o));
            g.Dorp = nuc ? 'd' : 'p';
            g.Constants(Array.Empty<string>());
            Err(MafftGlobals.ParseTbfast(tb).Log + "treein = " + (o.GuideTree != MafftGuideTree.Default ? "1" : "0") + "\ncompacttree = 0\n"
                + MafftRun.StackSizeMessage(8192, nseq) + "whereispairalign = " + _installation.LibraryDirectory + "\n" + g.Log + "All-to-all alignment.\nPreparing bpp\nCannot open hat4\n");
            return 1;
        }
    }
}
