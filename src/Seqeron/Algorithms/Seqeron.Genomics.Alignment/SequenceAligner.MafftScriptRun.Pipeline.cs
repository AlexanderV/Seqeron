#pragma warning disable S3626, S3871, S1643, S2325, IDE0052, IDE0051, IDE0059, S127, S3358, S1066, S4487, S1144, CA1823, CS0414, CA1870, S1125, S3776, S107, S138, S1541, S1854
using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    internal sealed partial class MafftScriptRun
    {
        /// <summary>The library options that the variables of the script stand for (the arguments the script passes to the binaries are built from them by <see cref="MafftScript"/>).</summary>
        private MafftOptions BuildOptions(IReadOnlyList<string> names, bool nuc)
        {
            _ = nuc;
            var o = new MafftOptions
            {
                SequenceType = _seqtype == "-D" ? MafftSequenceType.Nucleotide : _seqtype == "-P" ? MafftSequenceType.Amino : MafftSequenceType.Auto,
                GapOpen = Num(_gop),
                GapOpenDistance = _opdistspecified ? Num(_gopdist) : null,
                Offset = -Num(_aof),
                GapExtension = -Num(_gexp),
                AnySymbol = _anysymbol,
                Retree = _retreeUser,
                Fft = _fft,
                ForceFft = _forcefft,
                TerminalGapPenalty = _termgapopt.Trim().Length == 0,
                Unweight = _weightopt.Trim() == "-u",
                TupleSize = _tuplesize,
                TreeMixing = TreeMixing(),
                LegacyGapPenalty = _legacygapopt == "-L",
                Threads = _numthreads,
                ThreadsIterate = _numthreadsit,
                MaxIterate = _iterate,
                GuideTree = GuideTreeMethod(),
                GuideTreeText = (GuideTreeMethod() == MafftGuideTree.UserTree || (GuideTreeMethod() is MafftGuideTree.MemSaveTree or MafftGuideTree.StepAdd && _treein)) && _files.TryGetValue(_treeinfile, out var treeText) ? treeText : null,
                LhLimit = _lhlimit.Trim().Length == 0 ? int.MaxValue : MafftGlobals.CAtoi(_lhlimit.Trim().Split(' ')[1]),
                Reorder = _distance == "parttree" ? _partorderopt.Trim().Length == 0 : _outorder == "aligned",
                PartSize = _partsize,
                GroupSize = _groupsize,
                TreeOut = _treeout && !_nodeout,
                NodeOut = _nodeout,
                DistOut = _distout,
                Names = names,
                OldGenafParameters = _oldgenafparam,
                AddCount = _nadd,
                MergeTable = _mergetable != "/dev/null" && _files.TryGetValue(_mergetable, out var mergeText) ? MafftScript.NormalizeScriptFile(mergeText) : null,
                AddProfile = _addarg0.Trim() == "-I",
                KeepLength = _add2ndhalfarg.Contains("-Y", StringComparison.Ordinal),
                MapOut = MapOutMode(),
                Smoothing = _add2ndhalfarg.Contains("-p", StringComparison.Ordinal),
                NamesAreProgramTitles = true,
                MinimumWeight = Num(_minimumweight),
                FixThreshold = Num(_fixthreshold),
                AaMatrix = _sbstmodel.Contains("-b -1", StringComparison.Ordinal) ? (_files.TryGetValue(_scorematrix, out var aaText) ? MafftProgramIo.DropEmptyLines(MafftProgramIo.CarriageReturnsToNewlines(aaText)) : "")
                    : _sbstmodel.Contains("-b -2", StringComparison.Ordinal) ? (_files.TryGetValue(_scorematrix, out var textMatrix) ? MafftProgramIo.DropEmptyLines(MafftProgramIo.CarriageReturnsToNewlines(textMatrix)) : "") : null,
                TextMode = _f2clext == "-E",
                CoreOut = _coreout,
                CoreExtend = _coreext.Trim() == "-c",
                CoreThreshold = _corethr,
                CoreWindow = _corewin,
                RandomSeed = _randomseed,
                NoSegmentation = _bunkatsuopt.Trim() == "-B",
                OneIteration = _oneiterationopt.Trim() == "-r",
                NoScore = _scorecalcopt.Trim() == "-Z",
                ParallelizationStrategy = _parallelizationstrategy,
                ScoreOut = _scoreoutarg.Trim().Length > 0,
                OutNumber = _outnum.Trim() == "-n",
                SparsePickup = _npickup,
                ThreadsProgressive = _numthreadstb,
                MemoryMode = _memopt.Trim() == "-M -B" ? MafftMemoryMode.MemSave : _memopt.Trim() == "-N" ? MafftMemoryMode.NoMemSave : MafftMemoryMode.Default,
                PartTree = _distance != "parttree" ? MafftPartTree.None : _partdist == "fasta" ? MafftPartTree.Fasta : _partdist == "localalign" ? MafftPartTree.DynamicProgramming : MafftPartTree.KmerDistance,
                PairLibrary = _distance switch { "global" => MafftPairLibrary.Global, "local" => MafftPairLibrary.Local, "localgenaf" => MafftPairLibrary.Genaf, _ => MafftPairLibrary.None },
                WeightI = Num(_weighti),
                GlobalPairOpen = Num(_pggop),
                GlobalPairOffset = Num(_pgaof),
                GlobalPairExtension = Num(_pgexp),
                LocalPairOpen = Num(_lgop),
                LocalPairOffset = Num(_laof),
                LocalPairExtension = Num(_lexp),
                GenafOpen = Num(_LGOP),
                GenafExtension = Num(_LEXP),
            };
            return ApplyModel(o);
        }

        /// <summary>The method that the script lets <c>-U</c> read: <c>--pileup</c> / <c>--randomchain</c> / <c>--large</c> ... write their name to <c>_guidetree</c> (only with the 6-mer distances), <c>--treein</c> gives the tree of the user.</summary>
        private MafftGuideTree GuideTreeMethod()
        {
            if (_treeext is "memsavetree" or "stepadd" && _distance is "global" or "local" or "localgenaf")
                return _treeext == "stepadd" ? MafftGuideTree.StepAdd : MafftGuideTree.MemSaveTree;
            if (_treeext != "none" && _distance == "ktuples")
            {
                return _treeext switch
                {
                    "randomchain" => MafftGuideTree.RandomChain,
                    "pileup" => MafftGuideTree.Pileup,
                    "memsavetree" => MafftGuideTree.MemSaveTree,
                    "memsavetreex" => MafftGuideTree.MemSaveTreeX,
                    "stepadd" => MafftGuideTree.StepAdd,
                    "youngestlinkage" => MafftGuideTree.YoungestLinkage,
                    _ => MafftGuideTree.Default,
                };
            }
            return _treein ? MafftGuideTree.UserTree : MafftGuideTree.Default;
        }

        /// <summary>The last of <c>-Z</c> (1) / <c>-z</c> (2) in the options that the script passes on for the added sequences.</summary>
        private int MapOutMode()
        {
            int mode = 0;
            foreach (var t in _add2ndhalfarg.Split(' ', StringSplitOptions.RemoveEmptyEntries))
            {
                if (t == "-Z")
                    mode = 1;
                else if (t == "-z")
                    mode = 2;
            }
            return mode;
        }

        private double? TreeMixing()
        {
            int i = _treealg.IndexOf("-X", StringComparison.Ordinal);
            string v = _treealg[(i + 2)..].Trim();
            return v == "0.1" ? null : Num(v);
        }

        /// <summary><c>$model</c>: <c>-b</c> / <c>-j</c> / <c>-m</c>, <c>-k</c>, <c>-a</c>, <c>-:</c>.</summary>
        private MafftOptions ApplyModel(MafftOptions o)
        {
            var w = _sbstmodel.Split(' ', StringSplitOptions.RemoveEmptyEntries);
            if (w.Length >= 2)
            {
                int v = MafftGlobals.CAtoi(w[1]);
                o = w[0] switch
                {
                    "-b" => o with { ScoringModel = MafftScoringModel.Blosum, Blosum = v },
                    "-j" => o with { ScoringModel = MafftScoringModel.Jtt, Pam = v },
                    "-m" => o with { ScoringModel = MafftScoringModel.Tm, Pam = v },
                    _ => o,
                };
            }
            var k = _kappa.Split(' ', StringSplitOptions.RemoveEmptyEntries);
            if (k.Length >= 2)
                o = o with { Kappa = MafftGlobals.CAtoi(k[1]) };
            return o with { FrequencyModel = _fmodel.Contains("-a", StringComparison.Ordinal) || w.Contains("-a"), WildcardN = _nmodel.Contains("-:", StringComparison.Ordinal) };
        }

        /// <summary>The title of sequence <paramref name="i"/> in the files of the programs.</summary>
        private string Title(int i, string name)
        {
            string t = _anysymbol ? Truncate(Prefix(i) + name, MafftProgramIo.NameLimit) : name;
            if (_outnum.Trim() == "-n")
            {
                int c = t.IndexOf("_numo_e_", StringComparison.Ordinal);
                string tail = c >= 0 ? t[(c + 8)..] : t;
                t = Truncate("_numo_s_" + (i + 1).ToString("D8", CultureInfo.InvariantCulture) + "_numo_e_" + tail, MafftProgramIo.NameLimit);
            }
            return t;
        }

        /// <summary>
        /// <c>--blastpair</c> and <c>--dash</c> need programs that the reference installation does not have (<c>blastall</c> for <c>dndblast</c>, the <c>dash_client</c>): the messages of the script and of <c>dndblast</c> after the run of <c>replaceu</c>.
        /// </summary>
        private int PairProgramMissing(List<MafftProgramIo.Record> recs)
        {
            var originals = recs.Select(r => MafftEngine.PlainOrAnyResidues(r.Raw, _anysymbol)).ToList();
            bool nuc1 = MafftEngine.GuessNucleotide(originals, new MafftOptions { SequenceType = _seqtype == "-D" ? MafftSequenceType.Nucleotide : _seqtype == "-P" ? MafftSequenceType.Amino : MafftSequenceType.Auto });
            if (_anysymbol)
                Prog(MafftEngine.ReplaceuLog(originals, nuc1));
            if (_newdash)
            {
                var plain = recs.Select(r => MafftEngine.PlainResidues(r.Raw)).ToList();
                if (MafftEngine.GuessNucleotide(plain, new MafftOptions()) && _seqtype != "-P")
                {
                    Prog("\nError: This data seems to be nucleotide sequences.\nAdd the --amino flag if this is surely protein.\n\n");
                    return 1;
                }
                Err("\n== Install DASH client =====================================================\nTo use this feature, uncomment the following line in Makefile\nDASH_CLIENT = dash_client\nand re-compile the source.\nNote that it requires the 'Go' compiler.\n============================================================================\n\n");
                return 1;
            }
            // dndblast: constants(), then "formatting .. " and the missing data file of blast
            var g = MafftGlobals.Parse(MafftScript.ModelArguments(BuildOptions(Array.Empty<string>(), nuc1)));
            g.Dorp = nuc1 ? 'd' : 'p';
            g.Constants(Array.Empty<string>());
            Prog(g.Log.ToString());
            Prog("formatting .. Cannot open datafile.\n");
            return 1;
        }

        private static string Prefix(int i) => "_os_" + i.ToString(CultureInfo.InvariantCulture) + "_oe_";

        private static string Truncate(string s, int n) => s.Length > n ? s[..n] : s;

        /// <summary><c>restoreu</c>: the text from "_os_" to the first "_oe_" is removed from the title.</summary>
        private static string RemovePrefix(string name)
        {
            int os = name.IndexOf("_os_", StringComparison.Ordinal);
            int oe = name.IndexOf("_oe_", StringComparison.Ordinal);
            return os < 0 || oe < 0 ? name : name[..os] + name[(oe + 4)..];
        }

        private int Pipeline(string infile, int nseq, int nadd, string strategy, string explanation, string performance, bool memsavetree, int cycletbfast, string localparam)
        {
            _ = nseq;
            _ = nadd;
            _ = localparam;


            if (_mergetable != "/dev/null" && _nadd > 0)
            {
                // mafft.tmpl: --merge with --add / --addprofile, after the run of replaceu
                var orig = MafftProgramIo.Records(infile).Select(r => MafftEngine.PlainOrAnyResidues(r.Raw, _anysymbol)).ToList();
                if (_anysymbol)
                    Prog(MafftEngine.ReplaceuLog(orig, MafftEngine.GuessNucleotide(orig, new MafftOptions { SequenceType = _seqtype == "-D" ? MafftSequenceType.Nucleotide : _seqtype == "-P" ? MafftSequenceType.Amino : MafftSequenceType.Auto })));
                Err("Impossible\n");
                return 1;
            }
            if (_distance == "blast" || _newdash)
                return PairProgramMissing(MafftProgramIo.Records(infile));
            if (_distance == "ktuplesmulti" && _fragment == 0)
            {
                // the chain of the script has no branch for it: "Impossible" after replaceu
                var orig = MafftProgramIo.Records(infile).Select(r => MafftEngine.PlainOrAnyResidues(r.Raw, _anysymbol)).ToList();
                if (_anysymbol)
                    Prog(MafftEngine.ReplaceuLog(orig, MafftEngine.GuessNucleotide(orig, new MafftOptions { SequenceType = _seqtype == "-D" ? MafftSequenceType.Nucleotide : _seqtype == "-P" ? MafftSequenceType.Amino : MafftSequenceType.Auto })));
                Err("Impossible\n");
                return 1;
            }
            bool engine = _fragment == 0 && _seed == "x" && _seedtable == "x" && !_newdash && _adjustdirection == 0 && _anchoropt.Length == 0
                && _distance is "ktuples" or "global" or "local" or "localgenaf" or "parttree";
            if (_nadd > 0 && _distance == "parttree")
                throw new NotSupportedException("--add with --parttree: the script refuses it (Impossible).");
            if (Num(_unalignlevel) > 0.0 || Num(_spfactor) < 100.0)
                throw new NotSupportedException("--allowshift / --unalignlevel / --shiftpenalty (the warp and the unaligned level of the alignment programs) are not ported yet.");
            if (!engine || cycletbfast > 1 && _distance != "parttree")
                throw new NotSupportedException("The strategy '" + strategy + "' of the script (distance " + _distance + ") is not ported yet.");
            var recs = MafftProgramIo.Records(infile);
            string pre;
            var names = new List<string>();
            foreach (var r in recs)
                names.Add(r.Name);
            // the title as the programs hold it (replaceu: prefix "_os_<i>_oe_"; --outnum: prefix "_numo_s_<i>_numo_e_"; both inside the 254 characters)
            var shown = new List<string>();
            for (int i = 0; i < recs.Count; i++)
                shown.Add(Title(i, names[i]));
            var o = BuildOptions(shown, false);
            MafftRunResult res;
            var seqs = recs.Select(r => r.Raw).ToList();
            if (_anysymbol && seqs.Any(s => s.IndexOfAny(new[] { '=', '<', '>' }) >= 0))
            {
                Prog("inputfile = orig\n\nCharacters '= < >' can be used only in the title lines in the --anysymbol or --text mode.\n\n");
                return 1;
            }
            string progress;
            if (recs.Count == 0)
            {
                res = new MafftRunResult(MultipleAlignmentResult.Empty, "");
                progress = "nthread = " + _numthreads.ToString(CultureInfo.InvariantCulture) + "\nnthreadpair = " + _numthreads.ToString(CultureInfo.InvariantCulture) + "\nnthreadtb = " + _numthreadstb.ToString(CultureInfo.InvariantCulture) + "\nppenalty_ex = 0\nWarning: Only 0 sequence found.\n";
            }
            else if (recs.Count == 1)
            {
                res = new MafftRunResult(MafftEngine.SingleSequence(seqs[0], o), "") { InputIndex = new[] { 0 } };
                progress = "nthread = " + _numthreads.ToString(CultureInfo.InvariantCulture) + "\nnthreadpair = " + _numthreads.ToString(CultureInfo.InvariantCulture) + "\nnthreadtb = " + _numthreadstb.ToString(CultureInfo.InvariantCulture) + "\nppenalty_ex = 0\nWarning: Only 1 sequence found.\n";
            }
            else
            {
                try
                {
                    res = MafftEngine.AlignMany(seqs, o, epilogue: false);
                }
                catch (MafftProgramFailureException f)
                {
                    Prog(f.Text);
                    return f.ExitCode;
                }
                catch (ArgumentException e)
                {
                    throw new NotSupportedException("The program refuses this input in a way that is not ported yet: " + e.Message, e);
                }
                progress = res.StandardError;
            }
            Prog(progress);
            if (res.TerminalText.Length > 0)
                Err(res.TerminalText); // the standard error of setcore is not redirected to the progress file
            if (_anysymbol && !_coreout && (res.NotUsedFile != null || res.InputIndex.Length == 0) && res.InputIndex.Length != recs.Count)
            {
                // restoreu: the alignment file holds the picked sequences only
                Err("ERROR!!\nIn input file,\nnjob = " + recs.Count.ToString(CultureInfo.InvariantCulture) + "\nbut in alignment file,\nnjob = " + res.InputIndex.Length.ToString(CultureInfo.InvariantCulture) + "\n");
                return 1;
            }
            string before = _partdist == "fasta" && _distance == "parttree" ? "Pairwise alignments were computed by FASTA\n(Pearson & Lipman, 1988, PNAS 85:2444-2448)\n" : "";
            Prog(MafftScript.EpilogueText(strategy, explanation, performance, before, MafftScript.GenafNoteText(_distance is "localgenaf" or "globalgenaf", _oldgenafparam)));

            // the output (the file pre, then f2cl for the other formats)
            var rows = res.Alignment.AlignedSequences.ToList();
            var index = res.InputIndex.Length == rows.Count ? res.InputIndex : Enumerable.Range(0, rows.Count).ToArray();
            string TitleOf(int k) => res.ChangedNames?[k] is { } changed ? changed[1..] : shown[k];
            var outNames = index.Select(k => _anysymbol ? RemovePrefix(TitleOf(k)) : TitleOf(k)).ToList();
            bool coreRaw = res.CoreNames != null && _outorder != "aligned" && _outputformat != "clustal" && _outputformat != "phylip";
            if (res.CoreNames != null)
                outNames = res.CoreNames.Select(nm => coreRaw ? nm : Truncate(nm, MafftProgramIo.NameLimit)).ToList(); // setcore writes the names as it read them: with the '=' of the programs and the prefix of replaceu; f2cl cuts them at 254 characters with the '='.
            if (res.GroupLabels is { } labels)
            {
                // splittbfast: sprintf( tmpname, "Group-%d %s", ... ), tmpname[B-1] = 0, then restoreu removes the prefix of replaceu
                for (int k = 0; k < outNames.Count; k++)
                {
                    outNames[k] = RemovePrefix(Truncate(labels[k] + " " + shown[index[k]], MafftProgramIo.NameLimit + 1));
                }
            }

            bool nuc = MafftEngine.GuessNucleotide(rows.Select(r => r.Replace("-", "", StringComparison.Ordinal)).ToList(), new MafftOptions());
            if (_outputformat == "clustal")
                pre = MafftProgramIo.Clustal(strategy, outNames, rows, _namelength, nuc, _f2clext == "-E");
            else if (_outputformat == "phylip")
                pre = MafftProgramIo.Phylip(outNames, rows, _namelength);
            else if (coreRaw)
                pre = string.Concat(outNames.Select((nm, k) => ">" + nm + "\n" + rows[k] + "\n")); // cat of the file written by setcore: one line per sequence
            else
                pre = (_outorder == "aligned" ? "" : res.OutputPrefix) + MafftProgramIo.Fasta(outNames, rows);
            if (_outputfile == "")
                _stdout.Append(pre);
            else
                _outFiles[_outputfile] = pre;
            if (!string.IsNullOrEmpty(res.MapFile))
                _outFiles[_mapoutfile == "/dev/null" ? _addfile + ".map" : _mapoutfile] = res.MapFile;
            if (_npickup != 0 && res.NotUsedFile != null)
                _outFiles[_infilename + ".notused"] = res.NotUsedFile;
            if (!string.IsNullOrEmpty(res.GuideTreeCopy))
                _outFiles["GuideTree"] = res.GuideTreeCopy;
            if (_treeout && res.TreeFile != null)
                _outFiles[_infilename + ".tree"] = res.TreeFile;
            if (_distout && res.DistanceFile != null)
                _outFiles[_infilename + ".hat2"] = res.DistanceFile;
            return 0;
        }
    }
}
