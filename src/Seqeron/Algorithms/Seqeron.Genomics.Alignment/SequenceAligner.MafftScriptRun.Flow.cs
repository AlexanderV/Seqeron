#pragma warning disable S3626, S3871, S1643, S2325, IDE0052, IDE0051, IDE0059, S127, S3358, S1066, S4487, S1144, CA1823, CS0414, CA1870, S1125, S3776, S107, S138, S1541, S1854
using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    internal sealed partial class MafftScriptRun
    {
        private int Rest(string infile, string addfileText)
        {
            string progname = "mafft";
            _ = progname;
            // ---- the number of threads (mafft.tmpl 1196 - 1243)
            if (_numthreads < 0)
            {
                int cores = Environment.ProcessorCount;
                Prog("OS = linux\n");
                Prog("The number of physical cores =  " + cores.ToString(CultureInfo.InvariantCulture) + "\n");
                _numthreads = cores;
            }
            if (_numthreadstb < 0)
                _numthreadstb = _numthreads < 16 || _mpi ? _numthreads : 16;
            if (_numthreadsit < 0)
                _numthreadsit = _numthreads < 8 ? _numthreads : 8;
            if (_numthreadsit == 0 && _parallelizationstrategy == "BESTFIRST")
            {
                Err("Impossible\n");
                Exit(1);
            }

            // ---- --add*, --auto (mafft.tmpl 1250 - 1385)
            int nadd;
            if (_addarg0 != " ")
            {
                nadd = MafftProgramIo.Records(addfileText).Count;
                if (nadd == 0)
                {
                    Err("Check " + _addfile + "\n");
                    Exit(1);
                }
                if (_seed != "x" || _seedtable != "x")
                {
                    Err("Impossible\nUse either ONE of --seed, --seedtable, --addprofile and --add.\n");
                    Exit(1);
                }
            }
            else
            {
                nadd = 0;
            }
            if (_auto)
            {
                if (infile.Length == 0)
                    throw new NotSupportedException("--auto with an unreadable input file: the script reads a file of its random temporary directory.");
                var recs = MafftProgramIo.Records(infile);
                int nseqCl = recs.Count;
                int nlenCl = recs.Count == 0 ? 0 : recs.Max(r => MafftEngine.PlainResidues(r.Raw).Replace("-", "", StringComparison.Ordinal).Length);
                if (nlenCl < 3000 && nseqCl < 100)
                {
                    _distance = "local"; _iterate = 1000; _cycle = 1; _treeext = "none";
                }
                else if (nlenCl < 1000 && nseqCl < 200)
                {
                    _distance = "local"; _iterate = 2; _cycle = 1; _treeext = "none";
                }
                else if (nlenCl < 10000 && nseqCl < 500)
                {
                    _distance = "ktuples"; _iterate = 2; _cycle = 2; _treeext = "none";
                }
                else if (nseqCl < 20000)
                {
                    _distance = "ktuples"; _iterate = 0; _cycle = 2; _treeext = "none";
                }
                else if (nseqCl < 100000)
                {
                    _distance = "ktuples"; _iterate = 0; _cycle = 2;
                    if (_fragment == 0 && _mergetable == "/dev/null")
                    {
                        _treeext = "memsavetree";
                    }
                }
                else if (nseqCl < 200000)
                {
                    _distance = "ktuples"; _iterate = 0; _cycle = 1;
                    if (_fragment == 0 && _mergetable == "/dev/null")
                    {
                        _treeext = "memsavetree";
                    }
                }
                else if (nlenCl < 3000)
                {
                    _distance = "parttree"; _partdist = "localalign"; _algopt = " "; _algoptit = " "; _cycle = 1;
                }
                else
                {
                    _distance = "parttree"; _partdist = "ktuples"; _algopt = " "; _algoptit = " "; _cycle = 1;
                }
                if (_fragment != 0)
                {
                    int norg = nseqCl - nadd;
                    long npair = (long)norg * nadd;
                    Prog("nadd = " + nadd.ToString(CultureInfo.InvariantCulture) + "\n");
                    Prog("npair = " + npair.ToString(CultureInfo.InvariantCulture) + "\n");
                    Prog("nseq = " + nseqCl.ToString(CultureInfo.InvariantCulture) + "\n");
                    Prog("nlen = " + nlenCl.ToString(CultureInfo.InvariantCulture) + "\n");
                    if (norg == 0)
                    {
                        Prog("\nThe reference sequence was removed because of ambiguous letters?\n\n");
                        Exit(1);
                    }
                    if (npair > 10000000 || nlenCl > 5000 || nadd > 500000)
                    {
                        _distance = "ktuples";
                        Prog("use ktuples, size=" + _tuplesize.ToString(CultureInfo.InvariantCulture) + "!\n");
                    }
                    else if (npair > 3000000 || nlenCl > 5000)
                    {
                        _distance = "multi";
                        _weighti = "0.0";
                        Prog("use multipair, weighti=0.0!\n");
                    }
                    else
                    {
                        _distance = "multi";
                        Prog("use multipair, weighti=" + _weighti + "!\n");
                    }
                    _pairspecified = true;
                }
            }

            // ---- the checks of the numbers (mafft.tmpl 1386 - 1488)
            if (Num(_sueff) < 0.0 || Num(_sueff) > 1.0)
            {
                Prog("\nThe argument of --mixedlinkage must be between 0.0 and 1.0\n\n");
                Exit(1);
            }
            if (Num(_maxambiguous) < 0.0 || Num(_maxambiguous) > 1.0)
            {
                Prog("\nThe argument of --maxambiguous must be between 0.0 and 1.0\n\n");
                Exit(1);
            }
            if (_allowshift)
            {
                if (!_unalignspecified)
                    _unalignlevel = "0.8";
                if (!_shiftpenaltyspecified)
                    _spfactor = "2.00";
            }
            if (!_opdistspecified)
                _gopdist = _gop;
            if (_unalignlevel != "0.0" || Num(_spfactor) < 100.0)
            {
                _nmodel = " -: ";
                _termgapopt = " ";
                if (_distance == "localgenaf")
                {
                    Prog("\nThe combination of --allowshift and --genafpair (E-INS-i/-1) is not supported.\nInstead, please try --allowshift --globalpair (G-INS-i/-1 in the web version),\nwhich covers the situation for --genafpair (E-INS-i/-1), too.\n\n");
                    Exit(1);
                }
                if (_distance != "global" || Num(_weighti) < 1.0)
                {
                    Prog("\nAt present, --unalignlevel # or --allowshift is supported only with the --globalpair option.\n\n");
                    Exit(1);
                }
                if (_fragment != 0)
                {
                    Prog("\nAt present, --unalignlevel # or --allowshift is not supported with the --addfragments option.\n\n");
                    Exit(1);
                }
            }
            if (Num(_spfactor) < 1.0)
            {
                Prog("\nshiftpenalty must be >1.\n");
                Exit(1);
            }
            if (Num(_fixthreshold) < 0.0)
            {
                Prog("\nThe 'fix' parameter must be >= 0.0\n\n");
                Exit(1);
            }
            if (Num(_unalignlevel) < 0.0 || Num(_unalignlevel) > 1.0)
            {
                Prog("\nThe 'unalignlevel' parameter must be between 0.0 and 1.0\n\n");
                Exit(1);
            }
            if (Num(_unalignlevel) > 0.0)
            {
                _laof = "0"; _lexp = "0"; _pgaof = "0"; _pgexp = "0"; _LEXP = "0"; _GEXP = "0";
                _termgapopt = " ";
                if (_fragment != 0)
                {
                    Prog("\nAt present, the 'unalignlevel > 0' mode is not supported with the --addfragments option.\n\n");
                    Exit(1);
                }
                if (_distance == "parttree")
                {
                    Prog("\nAt present, the 'unalignlevel > 0' mode is not supported in the (dp)parttree option.\n\n");
                    Exit(1);
                }
                if (_distance == "localgenaf")
                {
                    Prog("\nThe --genafpair is not supported in the 'unalignlevel > 0' mode.\nInstead, please try --unalignlevel xx --globalpair,\nwhich covers the situation for --genafpair (E-INS-i), too.\n\n");
                    Exit(1);
                }
            }
            if (Num(_similarityoffset) != 0.0 && Num(_unalignlevel) != 0.0)
            {
                Prog("\nDo not simultaneously specify --similaritylevel and --unalignlevel\n\n");
                Exit(1);
            }
            if (Num(_similarityoffset) < -1.0 || Num(_similarityoffset) > 1.0)
            {
                Prog("\nSimilarity must be between -1.0 and +1.0\n\n");
                Exit(1);
            }
            _aof = MafftScript.AwkPrint(0.0 + Num(_similarityoffset) + Num(_aof));
            _laof = MafftScript.AwkPrint(0.0 + Num(_similarityoffset) + Num(_laof));
            _pgaof = MafftScript.AwkPrint(0.0 + Num(_similarityoffset) + Num(_pgaof));
            int iteratelimit = _parallelizationstrategy is "BESTFIRST" or "BAATARI0" ? 254 : 16;
            if (_iterate > iteratelimit)
                _iterate = iteratelimit;
            if (_sprigorous)
            {
                _algopt = " -@ ";
                if (_iterate > 0)
                {
                    if (_numthreadsit == 0)
                    {
                        _algoptit = " -@ -B -Z -z 1000 ";
                    }
                    else
                    {
                        Prog("\nAt present, the combination of --sp and iterative refinement is supported only in a single thread.\nPlease try \"--thread -1 --threadit 0\", which runs the iterative refinment calculation on a single thread.\n\n");
                        Exit(1);
                    }
                }
                _termgapopt = " ";
                _fft = false;
                _memopt = " -N ";
            }
            if (_er == 1)
            {
                Err(Usage);
                Exit(1);
            }
            _ = _sw;

            // ---- the programs that must be installed (mafft.tmpl 1640 - 1740): fasta34 is part of the port (the installation of the oracle has it)
            if (_distance is "last" or "lastmulti")
            {
                Err("\n== Install LAST ============================================================\nLAST (Kielbasa, Wan, Sato, Horton, Frith 2011 Genome Res. 21:487) is required.\nhttp://last.cbrc.jp/\nhttp://mafft.cbrc.jp/alignment/software/xxxxxxx.html \n============================================================================\n\n");
                Exit(1);
            }
            if (_distance is "lara" or "slara")
            {
                Err("\n== Install LaRA =========================================================\nThis option requires LaRA (Bauer et al. http://www.planet-lisa.net/).\nThe executable have to be renamed to 'mafft_lara' and installed into \nthe " + _installation.LibraryDirectory + " directory. \nA configuration file of LaRA also have to be given\nmafft-xinsi --larapair --laraparams parameter_file\nmafft-xinsi --slarapair --laraparams parameter_file\n=========================================================================\n\n");
                Exit(1);
            }
            if (_distance is "foldalignlocal" or "foldalignglobal")
            {
                Err("\n== Install FOLDALIGN ====================================================\nThis option requires FOLDALIGN (Havgaard et al. http://foldalign.ku.dk/).\nThe executable have to be renamed to 'foldalign210' and installed into \nthe " + _installation.LibraryDirectory + " directory. \n=========================================================================\n\n");
                Exit(1);
            }
            if ((_distance == "scarna" || _mccaskill) && !_installation.Mxscarna)
            {
                Err("\n== Install MXSCARNA ======================================================\nMXSCARNA (Tabei et al. BMC Bioinformatics 2008 9:33) is required.\nPlease 'make' at the 'extensions' directory of the MAFFT source package,\nwhich contains the modified version of MXSCARNA.\nhttp://mafft.cbrc.jp/alignment/software/source.html \n==========================================================================\n\n");
                Exit(1);
            }
            if ((_distance == "dafs" || _dafs) && !_installation.Dafs)
            {
                Err("\n== Install DAFS===========================================================\nDAFS (Sato et al. Journal 2012 issue:page) is required.\nhttp://www.ncrna.org/ \n==========================================================================\n\n");
                Exit(1);
            }
            if (_contrafold && !_installation.Contrafold)
            {
                Err("\n== Install CONTRAfold ===================================================\nThis option requires CONTRAfold\n(Do et al. http://contra.stanford.edu/contrafold/).\nThe executable 'contrafold' have to be installed into \nthe " + _installation.LibraryDirectory + " directory. \n=========================================================================\n\n");
                Exit(1);
            }

            // ---- --retree 0, the tree and distance files (mafft.tmpl 1745 - 1790)
            _retreeUser = _cycle;
            if (_cycle == 0)
            {
                _iterate = 0;
                _weighti = "0.0";
                if (_treeout)
                    _groupsize = 1;
            }
            else if (_nodeout)
            {
                if (_iterate > 0)
                {
                    Err("The --nodeout option supports only progressive method (--maxiterate 0) for now.\n");
                    Exit(1);
                }
            }

            // ---- the input (mafft.tmpl 1806 - 1830)
            var lines = MafftProgramIo.SplitLines(infile);
            var bad = lines.Select((l, i) => (l, i)).Where(t => t.l.Length > 0 && IsBlankLed(t.l)).ToList();
            if (bad.Count > 0)
            {
                Err("The first character of a description line must be \nthe greater-than (>) symbol, not a blank.\nPlease check the format around the following line(s):\n");
                foreach (var (l, i) in bad)
                    Err((i + 1).ToString(CultureInfo.InvariantCulture) + ":" + l + "\n");
                Exit(1);
            }
            int nseq = MafftProgramIo.CountRecordLines(infile);
            if (nseq == 2)
                _cycle = 1;
            if (_cycle > 3)
                _cycle = 3;
            if (nseq > 60000 && _iterate > 1)
            {
                Err("Too many sequences to perform iterative refinement!\nPlease use a progressive method.\n");
                Exit(1);
            }
            if (_distance is "lastmulti" or "multi")
            {
                if (_fragment == 0)
                {
                    Err("Specify --addfragments too\n");
                    Exit(1);
                }
            }
            if (_fragment != 0)
            {
                if (!_pairspecified)
                    _distance = "multi";
                if (_distance is not ("multi" or "hybrid" or "lastmulti" or "local" or "last" or "ktuples" or "ktuplesmulti"))
                {
                    Err("Specify --multipair, --lastmultipair, --lastpair, --localpair, --6merpair, --6mermultipair or --hybridpair\n");
                    Exit(1);
                }
            }
            if (_memopt == " -M -B " && _distance != "ktuples")
            {
                Err("Impossible\n");
                Exit(1);
            }
            string splitopt = "";
            if (_distance == "parttree")
            {
                if (_mergetable != "/dev/null")
                {
                    Err("The combination of (dp)parttree and merge is Impossible.  \n");
                    Exit(1);
                }
                if (_addfile != "/dev/null")
                {
                    Err("The combination of (dp)parttree and add(fragments) is Impossible.  \n");
                    Exit(1);
                }
                if (_seed != "x" || _seedtable != "x")
                {
                    Err("Impossible\n");
                    Exit(1);
                }
                if (_iterate > 1)
                {
                    Err("Impossible\n");
                    Exit(1);
                }
                if (_outorder == "aligned")
                    _outorder = "input";
                _outorder = "input";
                if (_partdist == "localalign")
                {
                    splitopt = " -U ";
                    _cycle = 1;
                }
                else if (_partdist == "fasta")
                {
                    splitopt = " -S ";
                    _cycle = 1;
                }
            }
            _ = splitopt;
            bool seedLike = _seed != "x" || _seedtable != "x" || _ownlist != "/dev/null" || _pdblist != "/dev/null" || _newdash;
            string localparam;
            if (_distance is "ktuples" or "ktuplesmulti" && !seedLike)
            {
                localparam = "";
                _weighti = "0.0";
            }
            else if (_distance is "ktuples" or "ktuplesmulti")
            {
                if (_cycle < 2)
                    _cycle = 2;
                if (_iterate < 2)
                {
                    Err("############################################################################\n# Warning:\n#   Progressive alignment method is incompatible with the --seed option.\n#   Automatically switched to the iterative refinement method.\n#   \n# Also consider using the '--add' option, which is compatible with\n#   the progressive method and FASTER than the '--seed' option.\n#   Usage is:\n#   % mafft --add newSequences existingAlignment > output\n############################################################################\n");
                    _iterate = 2;
                }
                localparam = "-l " + _weighti;
            }
            else if (_distance == "parttree")
            {
                localparam = "";
                _weighti = "0.0";
                if (_groupsize > -1)
                    _cycle = 1;
            }
            else
            {
                localparam = "-l " + _weighti;
                if (_cycle > 1)
                    _cycle = 1;
            }
            if (_distance is "localgenaf" or "globalgenaf")
            {
                _aof = "0.000";
                if (!_oldgenafparam)
                {
                    _laof = "0.0";
                    _lexp = "0.0";
                    _usenaivepairscore = "-Z";
                }
            }
            if (_forcefft)
                _fft = true;
            if (_seed != "x" && _seedtable != "x")
            {
                Err("Use either one of seedtable and seed.  Not both.\n");
                Exit(1);
            }
            if (_f2clext == "-E" && _anysymbol)
            {
                Err("\nThe combination of --text and ( --anysymbol or --preservecase ) is impossible.\n\n");
                Exit(1);
            }
            bool memsavetree = false;
            if (_treeext != "none")
            {
                if (_distance == "ktuples")
                {
                    _treein = true;
                    _treeinopt = " -U ";
                    if (_treeext is "randomchain" or "pileup")
                        _cycle = 1;
                    else
                        memsavetree = true;
                }
                else
                {
                    if (_treeext is "memsavetree" or "stepadd")
                    {
                        memsavetree = true;
                    }
                    else
                    {
                        Prog("With globalpair, localpair or genafpair,\nUse --large, --minimumlinkage, --averagelinkage or --mixedlinkage.\n--" + _treeext + " is supported only with --6merpair.\n\n");
                        Exit(0);
                    }
                }
                if (_iterate > 0)
                {
                    Prog("Iterative refinment is not supported for --large or --" + _treeext + "\n\n");
                    Exit(1);
                }
                if (_fragment != 0)
                {
                    Prog("--addfragments, --addfull or --addlong is not yet supported for --large or --" + _treeext + "\nUse --add newsequences --" + _treeext + "\nOr, --addfragments (long, full) newsequences, without --" + _treeext + "\n\n");
                    Exit(1);
                }
                if (_mergetable != "/dev/null")
                {
                    Prog("--merge is not yet supported for --large or --" + _treeext + "\nUse --merge without --" + _treeext + "\n\n");
                    Exit(1);
                }
            }
            if (nadd > 0)
            {
                if (_fragment == 1)
                {
                    _addarg = _addarg0 + " " + nadd.ToString(CultureInfo.InvariantCulture) + " -g -0.01";
                    _addsinglearg = "";
                    _cycle = 1;
                    _iterate = 0;
                }
                else if (_fragment == -1)
                {
                    _addarg = _addarg0 + " " + nadd.ToString(CultureInfo.InvariantCulture);
                    _addsinglearg = "-V";
                    _cycle = 1;
                    _iterate = 0;
                }
                else if (_fragment == -2)
                {
                    _addarg = _addarg0 + " " + nadd.ToString(CultureInfo.InvariantCulture);
                    _addsinglearg = "-V";
                    _add2ndhalfarg += " -p ";
                    _cycle = 1;
                    _usenaivepairscore = "-Z";
                    _laof = "0.0";
                    _lexp = "0.0";
                    _iterate = 0;
                }
                else if (_fragment == -3)
                {
                    _addarg = _addarg0 + " " + nadd.ToString(CultureInfo.InvariantCulture);
                    _addsinglearg = "-x";
                    _cycle = 1;
                    _iterate = 0;
                }
                else
                {
                    _addarg = _addarg0 + " " + nadd.ToString(CultureInfo.InvariantCulture);
                    _addsinglearg = "";
                    _bunkatsuopt = " -B ";
                    if (_add2ndhalfarg != " ")
                    {
                        if (_auto || _iterate > 0)
                            _iterate = 0;
                    }
                }
                _nadd = nadd;
                if (_fragment != 0 && _distance is not ("ktuples" or "ktuplesmulti" or "multi" or "hybrid" or "local"))
                    NeedPort("--addfragments / --addfull / --addlong / --addtotop with --lastpair / --lastmultipair (needs the program lastal)");
            }
            else if (_add2ndhalfarg != " ")
            {
                Prog("\nThe --keeplength and --mapout options are supported\nonly with --add, --addfragments or --addlong.\n\n");
                Exit(1);
            }
            if (_codonposfile != "/dev/null" || _codonscorefile != "/dev/null")
            {
                if (nadd == 0 || _fragment == 0)
                {
                    Prog("\n'--codonpos' and '--codonscore' options are supported only with the '--6merpair --addfragments' option.\n\n");
                    Exit(1);
                }
                if (_distance != "ktuples")
                {
                    Prog("\n'--codonpos' and '--codonscore' options are supported only with the '--6merpair --addfragments' option, at this point.\n\n");
                    Exit(1);
                }
            }
            int cycletbfast, cycledisttbfast;
            if (localparam.Length == 0 && _fragment == 0 && _distance != "parttree")
            {
                cycletbfast = 1;
                cycledisttbfast = _cycle;
                if (cycledisttbfast == 0)
                    cycledisttbfast = 1;
            }
            else
            {
                cycletbfast = _cycle;
                cycledisttbfast = 1;
            }
            _ = cycledisttbfast;
            if (_adjustdirection > 0 && _seed != "x")
            {
                Err("\nThe combination of --adjustdirection(accurately) and --seed is not supported.\n\n");
                Exit(1);
            }
            if (_mccaskill || _dafs || _rnaalifold || _contrafold)
            {
                if (_distance == "ktuples")
                {
                    Err("Not supported.\nPlease add --globalpair, --localpair, --scarnapair, --dafspair\n--larapair, --slarapair, --foldalignlocalpair or --foldalignglobalpair\n");
                    Exit(1);
                }
                if (_f2clext == "-E")
                {
                    Err("\nFor RNA alignment, the --text mode is impossible.\n\n");
                    Exit(1);
                }
            }

            // ---- the strategy (mafft.tmpl 2080 - 2200)
            string strategy = MafftScript.StrategyName(_treeext, _distance, _sw, _fft, _weighti, memsavetree, _iterate, _partdist, _fragment, _cycle, out string explanation, out string performance);
            if (_outorder == "aligned" && _outputformat == "pir")
                _ = 0;

            if (_portNeeded.Count > 0)
                throw new NotSupportedException("The option " + _portNeeded[0] + " of the script mafft is not ported yet.");
            if (nseq < 0)
                return 1;
            return Pipeline(infile, nseq, nadd, strategy, explanation, performance, memsavetree, cycletbfast, localparam);
        }

        private static bool IsBlankLed(string line)
        {
            int i = 0;
            while (i < line.Length && (line[i] == ' ' || line[i] == '\t'))
                i++;
            return i > 0 && i < line.Length && line[i] == '>';
        }

        private const string Usage =
            "------------------------------------------------------------------------------\n" +
            "  MAFFT " + Version + "\n" +
            "  https://mafft.cbrc.jp/alignment/software/\n" +
            "  MBE 30:772-780 (2013), NAR 30:3059-3066 (2002)\n" +
            "------------------------------------------------------------------------------\n" +
            "High speed:\n" +
            "  % mafft in > out\n" +
            "  % mafft --retree 1 in > out (fast)\n" +
            "\n" +
            "High accuracy (for <~200 sequences x <~2,000 aa/nt):\n" +
            "  % mafft --maxiterate 1000 --localpair  in > out (% linsi in > out is also ok)\n" +
            "  % mafft --maxiterate 1000 --genafpair  in > out (% einsi in > out)\n" +
            "  % mafft --maxiterate 1000 --globalpair in > out (% ginsi in > out)\n" +
            "\n" +
            "If unsure which option to use:\n" +
            "  % mafft --auto in > out\n" +
            "\n" +
            "--op # :         Gap opening penalty, default: 1.53\n" +
            "--ep # :         Offset (works like gap extension penalty), default: 0.0\n" +
            "--maxiterate # : Maximum number of iterative refinement, default: 0\n" +
            "--clustalout :   Output: clustal format, default: fasta\n" +
            "--reorder :      Outorder: aligned, default: input order\n" +
            "--quiet :        Do not report progress\n" +
            "--thread # :     Number of threads (if unsure, --thread -1)\n" +
            "--dash :         Add structural information (Rozewicki et al, submitted)\n";
    }
}
