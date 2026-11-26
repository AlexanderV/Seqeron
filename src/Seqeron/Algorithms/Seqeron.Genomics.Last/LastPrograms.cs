// SPDX-License-Identifier: GPL-3.0-or-later
// Derived from LAST 1256 (https://gitlab.com/mcfrith/last: lastdb.cc, lastal.cc, LastdbArguments.cc), Copyright 2008 - 2019 Martin C. Frith; GPL-3.0-or-later.
// Part of Seqeron.Genomics.Last, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Globalization;
using System.Text;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Last;

/// <summary>
/// The programs <c>lastdb</c> and <c>lastal</c> of LAST 1256 (DNA databases with the default <c>YASS</c> seed, adaptive seeds, gapless and gapped x-drop extension, TAB and MAF output),
/// run on a virtual working directory (file name -> text). They implement <see cref="ILastPrograms"/> of the MIT assembly, which is how <c>mafft --lastpair</c> reaches them.
/// </summary>
public sealed class LastPrograms : ILastPrograms
{
    /// <summary>The version that the oracle of the port ran (<c>lastal 1256</c>).</summary>
    public const string Version = "1256";

    /// <inheritdoc/>
    public LastProcessResult Run(string programPath, IReadOnlyList<string> arguments, IDictionary<string, string> workingDirectory, string standardInput = "")
    {
        ArgumentNullException.ThrowIfNull(programPath);
        ArgumentNullException.ThrowIfNull(arguments);
        ArgumentNullException.ThrowIfNull(workingDirectory);
        string name = programPath[(programPath.LastIndexOf('/') + 1)..];
        var argv = new List<string> { programPath };
        argv.AddRange(arguments);
        var stdout = new StringBuilder();
        try
        {
            switch (name)
            {
                case "lastdb":
                    return Lastdb(argv, workingDirectory, stdout);
                case "lastal":
                    return Lastal(argv, workingDirectory, standardInput, stdout);
                default:
                    throw new ArgumentException("LAST program not ported: " + name, nameof(programPath));
            }
        }
        catch (LastProgramExitException e)
        {
            return new LastProcessResult(e.Output, "", 0);
        }
        catch (LastProgramException e)
        {
            return new LastProcessResult(stdout.ToString(), programPath + ": " + e.Message + "\n", 1);
        }
    }

    // ------------------------------------------------------------------------------------------------------ lastdb

    private static LastProcessResult Lastdb(List<string> argv, IDictionary<string, string> files, StringBuilder stdout)
    {
        string programName = argv[0];
        var (options, rest, error) = LastGetopt.Scan(argv.Skip(1).ToList(), "hVpqR:cm:d:S:s:w:W:P:u:a:i:b:B:C:xDvQ:", programName);
        bool keepLowercase = false;
        int tantan = -1;
        string stderr = error ?? "";
        foreach (var (c, arg) in options)
        {
            switch (c)
            {
                case 'h':
                    throw new LastProgramExitException("Usage: " + programName + " [options] output-name fasta-sequence-file(s)\nPrepare sequences for subsequent alignment with lastal.\n");
                case 'V':
                    throw new LastProgramExitException("lastdb " + Version + "\n");
                case 'v':
                    break;
                case 'R':
                    if (arg is null || arg.Length != 2 || arg[0] < '0' || arg[0] > '1' || arg[1] < '0' || arg[1] > '3')
                        throw new LastProgramException("bad option value: -R " + arg);
                    keepLowercase = arg[0] == '1';
                    tantan = arg[1] - '0';
                    break;
                case '?':
                    return new LastProcessResult("", stderr + programName + ": bad option\n", 1);
                default:
                    throw new LastProgramException("option -" + c + " is not ported (protein and other alphabets, seeds, volumes, masking, strands)");
            }
        }

        if (rest.Count == 0)
            throw new LastProgramException("please give me an output name and sequence file(s)\n\nUsage: " + programName + " [options] output-name fasta-sequence-file(s)\nPrepare sequences for subsequent alignment with lastal.\n");
        string dbName = rest[0];
        if (tantan < 0)
            tantan = 1; // setDefaults: DNA, not -c
        var alph = new LastAlphabet(LastAlphabet.Dna);
        var multi = new LastMultiSequence();
        var inputs = rest.Skip(1).ToList();
        if (inputs.Count == 0)
            inputs.Add("-");
        bool first = true;
        var errOut = new StringBuilder(stderr);
        foreach (string input in inputs)
        {
            if (!files.TryGetValue(input, out string? text))
                throw new LastProgramException("can't open file: " + input);
            foreach (var (name, residues) in LastFasta.Read(text))
            {
                var bytes = Encoding.Latin1.GetBytes(residues);
                alph.Translate(bytes, 0, bytes.Length, keepLowercase);
                if (first)
                {
                    first = false;
                    if (IsDubiousDna(alph, bytes))
                        errOut.Append(programName).Append(": that's some funny-lookin DNA\n");
                }

                multi.BeginSequence(name);
                multi.Seq.AddRange(bytes);
                multi.FinishSequence();
            }
        }

        // the files of the database (the formats are private to this port: the program reads them back; the .prj file has the fields of lastdb)
        var letterCounts = new long[alph.Size];
        long letterTotal = 0;
        long maxSeqLen = 0;
        for (int i = 0; i < multi.FinishedSequences; i++)
        {
            long total = 0;
            for (int p = multi.SeqBeg(i); p < multi.SeqEnd(i); p++)
            {
                int u = alph.NumbersToUppercase[multi.Seq[p]];
                if (u < alph.Size)
                {
                    letterCounts[u]++;
                    total++;
                }
            }

            maxSeqLen = Math.Max(maxSeqLen, total);
            letterTotal += total;
        }

        var prj = new StringBuilder();
        prj.Append("version=").Append(Version).Append('\n');
        prj.Append("alphabet=").Append(alph.Letters).Append('\n');
        prj.Append("numofsequences=").Append(multi.FinishedSequences).Append('\n');
        prj.Append("numofletters=").Append(letterTotal).Append('\n');
        prj.Append("maxsequenceletters=").Append(maxSeqLen).Append('\n');
        prj.Append("letterfreqs=").Append(string.Join(' ', letterCounts)).Append('\n');
        prj.Append("maxunsortedinterval=0\n");
        prj.Append("keeplowercase=").Append(keepLowercase ? 1 : 0).Append('\n');
        if (tantan != 0)
            prj.Append("tantansetting=").Append(tantan).Append('\n');
        prj.Append("masklowercase=0\n");
        prj.Append("numofindexes=1\n");
        prj.Append("integersize=32\n");
        files[dbName + ".prj"] = prj.ToString();
        files[dbName + ".tis"] = Encoding.Latin1.GetString(multi.Seq.ToArray());
        files[dbName + ".ssp"] = string.Join('\n', multi.Ends);
        files[dbName + ".des"] = string.Join('\n', multi.Names.Select(n => n)) + (multi.Names.Count > 0 ? "\n" : "");
        return new LastProcessResult("", errOut.ToString(), 0);
    }

    /// <summary><c>isDubiousDna</c>: the first 100 letters of the first sequence have more than 10% unexpected letters.</summary>
    private static bool IsDubiousDna(LastAlphabet alph, byte[] seq)
    {
        int dnaCount = 0;
        for (int i = 0; i < 100; i++)
        {
            int c = i < seq.Length ? alph.NumbersToUppercase[seq[i]] : alph.Size; // the delimiter after the sequence
            if (c == alph.Size)
                return false;
            if (c < alph.Size || c == alph.Encode['N'])
                dnaCount++;
        }

        return dnaCount < 90;
    }

    // ------------------------------------------------------------------------------------------------------ lastal

    private static LastProcessResult Lastal(List<string> argv, IDictionary<string, string> files, string standardInput, StringBuilder stdout)
    {
        var args = new LastalArguments();
        args.FromArgs(argv);
        if (args.Message is not null)
        {
            // glibc getopt wrote its message to the standard error; the program then stops with "bad option"
            return new LastProcessResult("", args.Message + argv[0] + ": bad option\n", 1);
        }

        var run = new LastalRun(args, files, standardInput);
        run.Execute(stdout);
        return new LastProcessResult(stdout.ToString(), run.StandardError.ToString(), 0);
    }

    private sealed class LastalRun
    {
        private readonly LastalArguments _args;
        private readonly IDictionary<string, string> _files;
        private readonly string _stdin;
        private readonly LastAlphabet _alph = new(LastAlphabet.Dna);
        private readonly LastScoreMatrix _matrix = new();
        private readonly LastMatrixStats _stats = new();
        private readonly LastGappedAligner _aligner = new();
        private LastGapCosts _gap = null!;
        private Database _db = null!;
        private LastMultiSequence _query = new();
        private int _minScoreGapless;

        private sealed class Database
        {
            public byte[] Text = [];
            public LastMultiSequence Multi = new();
            public long Sequences, Letters;
            public bool KeepLowercase;
            public int Tantan;
            public bool CaseSensitiveSeeds;
            public LastSuffixIndex Index = null!;
        }

        private static Database ReadDatabase(IDictionary<string, string> files, string name)
        {
            if (!files.TryGetValue(name + ".prj", out string? prj))
                throw new LastProgramException("can't open file: " + name + ".prj");
            var db = new Database { Tantan = 0 };
            string? letters = null;
            bool seen = false;
            foreach (string line in prj.Split('\n'))
            {
                int eq = line.IndexOf('=');
                if (eq < 0)
                    continue;
                string key = line[..eq], value = line[(eq + 1)..].Trim();
                switch (key)
                {
                    case "alphabet":
                        letters = value;
                        break;
                    case "numofsequences":
                        db.Sequences = long.Parse(value, CultureInfo.InvariantCulture);
                        seen = true;
                        break;
                    case "numofletters":
                        db.Letters = long.Parse(value, CultureInfo.InvariantCulture);
                        break;
                    case "keeplowercase":
                        db.KeepLowercase = value == "1";
                        break;
                    case "tantansetting":
                        db.Tantan = int.Parse(value, CultureInfo.InvariantCulture);
                        break;
                    case "masklowercase":
                        db.CaseSensitiveSeeds = value == "1";
                        break;
                }
            }

            if (letters != LastAlphabet.Dna || !seen)
                throw new LastProgramException("can't read file: " + name + ".prj");
            if (db.CaseSensitiveSeeds)
                throw new LastProgramException("databases made with lastdb -c are not ported");
            if (!files.TryGetValue(name + ".tis", out string? tis) || !files.TryGetValue(name + ".ssp", out string? ssp) || !files.TryGetValue(name + ".des", out string? des))
                throw new LastProgramException("can't open file: " + name + ".tis");
            db.Text = Encoding.Latin1.GetBytes(tis);
            db.Multi.Seq.Clear();
            db.Multi.Seq.AddRange(db.Text);
            db.Multi.Ends.Clear();
            foreach (string e in ssp.Split('\n'))
                db.Multi.Ends.Add(int.Parse(e, CultureInfo.InvariantCulture));
            var names = des.Length == 0 ? [] : des.TrimEnd('\n').Split('\n');
            foreach (string n in names)
                db.Multi.BeginSequence(n);
            return db;
        }


        public StringBuilder StandardError { get; } = new();

        public LastalRun(LastalArguments args, IDictionary<string, string> files, string standardInput)
        {
            _args = args;
            _files = files;
            _stdin = standardInput;
        }

        private void Log(string message)
        {
            if (_args.Verbosity > 0)
                StandardError.Append(_args.ProgramName).Append(": ").Append(message).Append('\n');
        }

        public void Execute(StringBuilder stdout)
        {
            var a = _args;
            _db = ReadDatabase(_files, a.LastdbName);
            // the command line overrides the options stored with the database (the port stores none): isKeepLowercase / tantan settings
            if (a.OutputType == 0)
                throw new LastProgramException("option -j 0 (match counts) is not ported");
            if (a.OutputType > 3)
                throw new LastProgramException("option -j " + a.OutputType + " (probabilistic alignments) is not ported");
            if (a.Globality != 0)
                throw new LastProgramException("option -T 1 (overlap alignments) is not ported");
            if (a.OutputFormat is 'b' or 'B')
                throw new LastProgramException("output formats BlastTab and BlastTab+ are not ported");
            if (a.CullingLimitForGaplessAlignments != 0 || a.CullingLimitForFinalAlignments + 1 != 0)
                throw new LastProgramException("options -C and -K are not ported");
            if (a.MinimizerWindow > 1)
                throw new LastProgramException("option -W is not ported");
            if (a.MaxHitDepth + 1 != 0 || a.OneHitMultiplicity <= 0)
                throw new LastProgramException("options -L and -m 0 are not ported");
            int threads = a.NumOfThreads;
            if (threads < 1)
                throw new LastProgramException("option -P 0 is not ported");
            a.SetDefaultsFromAlphabet(_db.KeepLowercase, _db.Tantan, _db.CaseSensitiveSeeds, threads);
            if (a.MaskLowercase > 0)
                throw new LastProgramException("option -u is not ported");
            if (a.OutputType < 3 && a.ScoreType != 0)
                throw new LastProgramException("option -J is not ported");

            // makeScoreMatrix
            if (a.MatrixFile.Length > 0)
            {
                if (!_files.TryGetValue(a.MatrixFile, out string? matrixText))
                    throw new LastProgramException("can't open file: " + a.MatrixFile);
                _matrix.FromString(matrixText);
            }
            else
            {
                _matrix.SetMatchMismatch(a.MatchScore, a.MismatchCost, _alph.Letters);
            }

            _matrix.Init(_alph.Encode);
            if (a.OutputType > 0)
            {
                CalculateSubstitutionScoreMatrixStatistics();
                if (!_stats.IsBad)
                {
                    _matrix.AddAmbiguousScores(a.AmbiguousLetterOpt % 2 == 1, a.AmbiguousLetterOpt / 2 == 1, _alph.Encode, _stats.Lambda, _stats.LetterProbs1, _stats.LetterProbs2);
                    _matrix.Init(_alph.Encode);
                }
            }

            _gap = new LastGapCosts(a.DelOpenCosts[0], a.DelGrowCosts[0], a.InsOpenCosts[0], a.InsGrowCosts[0], a.GapPairCost);

            // calculateScoreStatistics: the gapless Karlin-Altschul parameters (for the default minimum score of gapless alignments)
            LastGaplessEvaluer? gapless = null;
            if (a.OutputType > 0 && !_stats.IsBad)
                gapless = LastGaplessEvaluer.TryCreate(_matrix.CaseSensitive, _alph.Size, _stats.LetterProbs1, _stats.LetterProbs2);
            double minScore = -1, eg2 = -1;
            // The E-value parameters of gapped alignments (the Sls / ALP library of LAST) are not ported: no E-value is reported and no threshold can be derived from one.
            a.SetDefaultsFromMatrix(_stats.Lambda, minScore, eg2);
            _minScoreGapless = CalcMinScoreGapless(gapless, _db.Letters);
            a.MinScoreGapless = _minScoreGapless;

            _alph.Translate(_query.Seq.ToArray(), 0, 0, a.IsKeepLowercase);
            _db.Index = new LastSuffixIndex(_db.Text, _db.Text.Length, LastSeed.Yass(_alph));

            WriteHeader(stdout);
            ReadAndAlignQueries(stdout, gapless);
        }

        private void CalculateSubstitutionScoreMatrixStatistics()
        {
            var scoreMat = new int[_alph.Size][];
            for (int i = 0; i < _alph.Size; i++)
                scoreMat[i] = _matrix.CaseSensitive[i].Take(_alph.Size).ToArray();
            if (_args.Temperature >= 0)
                throw new LastProgramException("option -t is not ported");
            _stats.CalcUnbiased(scoreMat, _alph.Size);
            if (_stats.IsBad)
                Log("can't calculate probabilities: maybe the mismatch costs are too weak");
        }

        private int CalcMinScoreGapless(LastGaplessEvaluer? gapless, double numLettersInReference)
        {
            var a = _args;
            if (a.MinScoreGapless >= 0)
                return a.MinScoreGapless;
            if (gapless is null)
                return (int)Math.Ceiling(a.MinScoreGapped);
            double n = a.MaxGaplessAlignmentsPerQueryPosition;
            if (a.MaxGaplessAlignmentsPerQueryPosition + 1 == 0)
                n = 10;
            double propConst = 0.0004;
            double e = n * 1 * propConst;
            double s = gapless.MinScore(e, numLettersInReference);
            s = Math.Max(1.0, s);
            return (int)Math.Ceiling(Math.Min(s, a.MinScoreGapped));
        }

        private void WriteHeader(StringBuilder o)
        {
            var a = _args;
            o.Append("# LAST version ").Append(Version).Append('\n');
            o.Append("#\n");
            o.Append(a.WriteCommented());
            o.Append("# Reference sequences=").Append(_db.Sequences.ToString(CultureInfo.InvariantCulture)).Append(" normal letters=").Append(_db.Letters.ToString(CultureInfo.InvariantCulture)).Append('\n');
            o.Append("#\n");
            o.Append(_matrix.WriteCommented());
            o.Append("#\n");
            o.Append("# Coordinates are 0-based.  For - strand matches, coordinates\n");
            o.Append("# in the reverse complement of the 2nd sequence are used.\n");
            o.Append("#\n");
            if (a.OutputFormat == 't')
                o.Append("# score\tname1\tstart1\talnSize1\tstrand1\tseqSize1\tname2\tstart2\talnSize2\tstrand2\tseqSize2\tblocks\n");
            if (a.OutputFormat == 'm')
                o.Append("# name start alnSize strand seqSize alignment\n#\n");
        }

        private void ReadAndAlignQueries(StringBuilder o, LastGaplessEvaluer? gapless)
        {
            var a = _args;
            long queryBatchCount = 0;
            long sequenceCount = 0;
            var inputs = a.Inputs.Count > 0 ? a.Inputs : ["-"];
            // all the records of all the inputs, read as the program reads them: sequence by sequence, a batch ends when the next sequence does not fit into the buffer
            var records = new List<(string Name, byte[] Residues)>();
            foreach (string input in inputs)
            {
                string text;
                if (input == "-")
                    text = _stdin;
                else if (!_files.TryGetValue(input, out string? t))
                    throw new LastProgramException("can't open file: " + input);
                else
                    text = t;
                foreach (var (name, residues) in LastFasta.Read(text))
                    records.Add((name, Encoding.Latin1.GetBytes(residues)));
            }

            _query = new LastMultiSequence();
            long maxSeqLen = a.BatchSize;
            foreach (var (name, residues) in records)
            {
                // appendFromFasta: the letters are read while the buffer is smaller than the limit (the first sequence of a batch has no limit);
                // the sequence is finished only if it was read completely and the pad still fits
                long limit = _query.FinishedSequences == 0 ? uint.MaxValue : maxSeqLen;
                bool finished = _query.Seq.Count + (long)residues.Length <= limit - LastMultiSequence.PadSize;
                if (!finished)
                {
                    // the batch of the finished sequences is aligned, the sequence that did not fit starts the next batch
                    o.Append("# batch ").Append(queryBatchCount++.ToString(CultureInfo.InvariantCulture)).Append('\n');
                    ScanAllVolumes(o, gapless);
                    _query.ReinitForAppending();
                }

                _query.BeginSequence(name);
                var buf = (byte[])residues.Clone();
                _alph.Translate(buf, 0, buf.Length, a.IsKeepLowercase);
                _query.Seq.AddRange(buf);
                _query.FinishSequence();
                sequenceCount++;
            }

            if (_query.FinishedSequences > 0)
            {
                o.Append("# batch ").Append(queryBatchCount.ToString(CultureInfo.InvariantCulture)).Append('\n');
                ScanAllVolumes(o, gapless);
            }

            o.Append("# Query sequences=").Append(sequenceCount.ToString(CultureInfo.InvariantCulture)).Append('\n');
        }

        private void ScanAllVolumes(StringBuilder o, LastGaplessEvaluer? gapless)
        {
            for (int i = 0; i < _query.FinishedSequences; i++)
                AlignOneQuery(o, i);
        }

        private void AlignOneQuery(StringBuilder o, int queryNum)
        {
            var a = _args;
            if (a.Strand != 0)
                TranslateAndScan(o, queryNum, false);
            if (a.Strand == 2 || a.Strand == 0)
                _query.ReverseComplementOneSequence(queryNum, _alph.Complement);
            if (a.Strand != 1)
                TranslateAndScan(o, queryNum, true);
        }

        private void TranslateAndScan(StringBuilder o, int queryNum, bool reverse)
        {
            Scan(o, queryNum);
        }

        private void Scan(StringBuilder o, int queryNum)
        {
            var a = _args;
            byte[] queryAll = _query.Seq.ToArray();
            int origin = _query.PadBeg(queryNum);
            byte[] query = queryAll[origin..];
            // the query array of LAST starts at the pad of the sequence
            var gaplessAlns = new LastSegmentPairPot();
            AlignGapless(gaplessAlns, queryNum, query, origin);
            if (a.OutputType == 1)
            {
                foreach (var sp in gaplessAlns.Items)
                    WriteSegmentPair(o, sp, queryNum);
                return;
            }

            if (gaplessAlns.Size == 0)
                return;
            var gappedAlns = new List<LastAlignment>();
            if (a.MaxDropFinal != a.MaxDropGapped)
            {
                AlignGapped(gappedAlns, gaplessAlns, query, a.MaxDropGapped, pregapped: true);
                gaplessAlns.Items = gaplessAlns.Items.Where(x => !LastSegmentPairPot.IsNotMarkedAsGood(x)).ToList();
            }

            AlignGapped(gappedAlns, gaplessAlns, query, a.MaxDropFinal, pregapped: false);
            if (gappedAlns.Count == 0)
                return;
            if (a.OutputType > 2)
                gappedAlns = EraseSuboptimal(gappedAlns);
            gappedAlns.Sort(MoreScore);
            foreach (var aln in gappedAlns)
                WriteAlignment(o, aln, queryNum);
        }

        private static int Cmp(uint x, uint y) => x.CompareTo(y);

        private static int LessBeg(LastAlignment x, LastAlignment y)
        {
            int c = Cmp(x.Beg1, y.Beg1);
            if (c != 0)
                return c;
            c = Cmp(x.Beg2, y.Beg2);
            if (c != 0)
                return c;
            if (x.Score != y.Score)
                return y.Score.CompareTo(x.Score);
            c = Cmp(x.Seed.Beg1, y.Seed.Beg1);
            return c != 0 ? c : Cmp(x.Seed.Beg2, y.Seed.Beg2);
        }

        private static int LessEnd(LastAlignment x, LastAlignment y)
        {
            int c = Cmp(x.End1, y.End1);
            if (c != 0)
                return c;
            c = Cmp(x.End2, y.End2);
            if (c != 0)
                return c;
            if (x.Score != y.Score)
                return y.Score.CompareTo(x.Score);
            c = Cmp(x.Seed.Beg1, y.Seed.Beg1);
            return c != 0 ? c : Cmp(x.Seed.Beg2, y.Seed.Beg2);
        }

        private static int MoreScore(LastAlignment x, LastAlignment y) => x.Score != y.Score ? y.Score.CompareTo(x.Score) : LessBeg(x, y);

        /// <summary><c>AlignmentPot::eraseSuboptimal</c>: of alignments with the same start (or the same end) only the first of the sort order stays.</summary>
        private static List<LastAlignment> EraseSuboptimal(List<LastAlignment> items)
        {
            if (items.Count == 0)
                return items;
            var l = items.OrderBy(x => x, Comparer<LastAlignment>.Create(LessBeg)).ToList();
            for (int i = 1; i < l.Count; i++)
                if (l[i].Beg1 == l[i - 1].Beg1 && l[i].Beg2 == l[i - 1].Beg2)
                    l[i].Marked = true;
            l = l.OrderBy(x => x, Comparer<LastAlignment>.Create(LessEnd)).ToList();
            for (int i = 1; i < l.Count; i++)
                if (l[i].End1 == l[i - 1].End1 && l[i].End2 == l[i - 1].End2)
                    l[i].Marked = true;
            return l.Where(x => !x.Marked).ToList();
        }

        private void AlignGapless(LastSegmentPairPot gaplessAlns, int queryNum, byte[] query, int origin)
        {
            var a = _args;
            var dt = new LastDiagonalTable();
            int loopBeg = _query.SeqBeg(queryNum) - origin;
            int loopEnd = _query.SeqEnd(queryNum) - origin;
            long minDepth = 1;
            if (a.MinHitDepth > minDepth)
                loopEnd -= (int)Math.Min(a.MinHitDepth - minDepth, loopEnd);
            long maxSignificantAlignments = a.MaxAlignmentsPerQueryStrand != 0 ? a.MaxAlignmentsPerQueryStrand : 1;
            if (a.MaxAlignmentsPerQueryStrand + 1 == 0)
                maxSignificantAlignments = long.MaxValue;
            int d = a.MaxDropGapless;
            int[][] m = _matrix.CaseInsensitive;
            byte[] text = _db.Text;
            for (int qryPos = loopBeg; qryPos < loopEnd; qryPos += (int)a.QueryStep)
            {
                var (beg, end) = _db.Index.Match(query, qryPos, (int)a.OneHitMultiplicity, (int)a.MinHitDepth, int.MaxValue);
                long maxAlignments = a.MaxGaplessAlignmentsPerQueryPosition + 1 == 0 ? long.MaxValue : a.MaxGaplessAlignmentsPerQueryPosition;
                for (int k = beg; k < end; k++)
                {
                    if (maxAlignments == 0)
                        break;
                    uint refPos = (uint)_db.Index.Positions[k];
                    if (dt.IsCovered((uint)qryPos, refPos))
                        continue;
                    int fs = LastGapless.ForwardScore(text, (int)refPos, query, qryPos, m, d);
                    int rs = LastGapless.ReverseScore(text, (int)refPos, query, qryPos, m, d);
                    int score = fs + rs;
                    if (score < _minScoreGapless)
                        continue;
                    int tEnd = LastGapless.ForwardEnd(text, (int)refPos, query, qryPos, m, fs);
                    int tBeg = LastGapless.ReverseEnd(text, (int)refPos, query, qryPos, m, rs);
                    int qBeg = qryPos - ((int)refPos - tBeg);
                    if (!LastGapless.IsOptimal(text, tBeg, tEnd, query, qBeg, m, d))
                        continue;
                    var sp = new LastSegmentPair((uint)tBeg, (uint)qBeg, (uint)(tEnd - tBeg), score);
                    dt.AddEndpoint(sp.End2, sp.End1);
                    gaplessAlns.Add(sp);
                    maxAlignments--;
                    if (score >= a.MinScoreGapped && --maxSignificantAlignments == 0)
                        return;
                }
            }
        }

        private void AlignGapped(List<LastAlignment> gappedAlns, LastSegmentPairPot gaplessAlns, byte[] query, int maxDrop, bool pregapped)
        {
            var a = _args;
            int[][] m = _matrix.CaseInsensitive;
            byte[] text = _db.Text;
            for (int i = 0; i < gaplessAlns.Items.Count; i++)
            {
                var sp = gaplessAlns.Items[i];
                int fs = LastGapless.ForwardScore(text, (int)sp.Beg1, query, (int)sp.Beg2, m, maxDrop);
                int rs = LastGapless.ReverseScore(text, (int)sp.Beg1, query, (int)sp.Beg2, m, maxDrop);
                int tEnd = LastGapless.ForwardEnd(text, (int)sp.Beg1, query, (int)sp.Beg2, m, fs);
                int tBeg = LastGapless.ReverseEnd(text, (int)sp.Beg1, query, (int)sp.Beg2, m, rs);
                int qBeg = (int)sp.Beg2 - ((int)sp.Beg1 - tBeg);
                sp.Start1 = (uint)tBeg;
                sp.Start2 = (uint)qBeg;
                sp.Size = (uint)(tEnd - tBeg);
                sp.Score = fs + rs;
                if (!LastGapless.IsOptimal(text, tBeg, tEnd, query, qBeg, m, maxDrop))
                    LastSegmentPairPot.Mark(sp);
            }

            gaplessAlns.Items = gaplessAlns.Items.Where(x => !LastSegmentPairPot.IsMarked(x)).ToList();
            gaplessAlns.Sort();
            int gappedAlignmentCount = 0;
            for (int i = 0; i < gaplessAlns.Size; i++)
            {
                var sp = gaplessAlns.Get(i);
                if (LastSegmentPairPot.IsMarked(sp))
                    continue;
                var aln = new LastAlignment { Seed = sp.Clone() };
                aln.Seed.MaxIdenticalRun(text, query, _alph.NumbersToUppercase, _alph.NumbersToUppercase);
                aln.Seed.Score = LastGapless.AlignmentScore(text, (int)aln.Seed.Beg1, (int)aln.Seed.End1, query, (int)aln.Seed.Beg2, m);
                aln.MakeXdrop(_aligner, text, query, m, _matrix.MaxScore, _gap, maxDrop);
                if (aln.Score < a.MinScoreGapped)
                    continue;
                if (!aln.IsOptimal(text, query, m, maxDrop, _gap))
                    continue;
                gaplessAlns.MarkAllOverlaps(aln.Blocks);
                gaplessAlns.MarkTandemRepeats(aln.Seed, (uint)a.MaxRepeatDistance);
                if (!pregapped)
                    gappedAlns.Add(aln);
                else
                    LastSegmentPairPot.MarkAsGood(sp);
                gappedAlignmentCount++;
                if (a.MaxAlignmentsPerQueryStrand + 1 != 0 && gappedAlignmentCount >= a.MaxAlignmentsPerQueryStrand)
                    break;
            }
        }

        private void WriteSegmentPair(StringBuilder o, LastSegmentPair sp, int queryNum)
        {
            var aln = new LastAlignment { Seed = sp, Score = sp.Score };
            aln.Blocks.Add(sp.Clone());
            WriteAlignment(o, aln, queryNum);
        }

        private void WriteAlignment(StringBuilder o, LastAlignment aln, int queryNum)
        {
            if (_args.OutputFormat == 'm')
                WriteMaf(o, aln, queryNum);
            else
                WriteTab(o, aln, queryNum);
        }

        private static string ScoreText(double score) => score.ToString("F0", CultureInfo.InvariantCulture);

        private void WriteTab(StringBuilder o, LastAlignment aln, int queryNum)
        {
            int seqNum1 = _db.Multi.WhichSequence((int)aln.Beg1);
            long seqStart1 = _db.Multi.SeqBeg(seqNum1);
            long seqLen1 = _db.Multi.SeqLen(seqNum1);
            long seqStart2 = _query.SeqBeg(queryNum) - _query.PadBeg(queryNum);
            long seqLen2 = _query.SeqLen(queryNum);
            o.Append(ScoreText(aln.Score)).Append('\t');
            o.Append(_db.Multi.SeqName(seqNum1)).Append('\t').Append((aln.Beg1 - seqStart1).ToString(CultureInfo.InvariantCulture)).Append('\t').Append((aln.End1 - aln.Beg1).ToString(CultureInfo.InvariantCulture)).Append('\t').Append('+').Append('\t').Append(seqLen1.ToString(CultureInfo.InvariantCulture)).Append('\t');
            o.Append(_query.SeqName(queryNum)).Append('\t').Append((aln.Beg2 - seqStart2).ToString(CultureInfo.InvariantCulture)).Append('\t').Append((aln.End2 - aln.Beg2).ToString(CultureInfo.InvariantCulture)).Append('\t').Append(_query.Strand(queryNum)).Append('\t').Append(seqLen2.ToString(CultureInfo.InvariantCulture)).Append('\t');
            o.Append(aln.BlockText()).Append('\n');
        }

        private void WriteMaf(StringBuilder o, LastAlignment aln, int queryNum)
        {
            int seqNum1 = _db.Multi.WhichSequence((int)aln.Beg1);
            long seqStart1 = _db.Multi.SeqBeg(seqNum1);
            long seqLen1 = _db.Multi.SeqLen(seqNum1);
            long seqStart2 = _query.SeqBeg(queryNum) - _query.PadBeg(queryNum);
            long seqLen2 = _query.SeqLen(queryNum);
            string n1 = _db.Multi.SeqName(seqNum1), n2 = _query.SeqName(queryNum);
            string b1 = (aln.Beg1 - seqStart1).ToString(CultureInfo.InvariantCulture), b2 = (aln.Beg2 - seqStart2).ToString(CultureInfo.InvariantCulture);
            string r1 = (aln.End1 - aln.Beg1).ToString(CultureInfo.InvariantCulture), r2 = (aln.End2 - aln.Beg2).ToString(CultureInfo.InvariantCulture);
            string s1 = seqLen1.ToString(CultureInfo.InvariantCulture), s2 = seqLen2.ToString(CultureInfo.InvariantCulture);
            int nw = Math.Max(n1.Length, n2.Length), bw = Math.Max(b1.Length, b2.Length), rw = Math.Max(r1.Length, r2.Length), sw = Math.Max(s1.Length, s2.Length);
            o.Append("a score=").Append(ScoreText(aln.Score)).Append('\n');
            byte[] queryAll = _query.Seq.ToArray();
            int origin = _query.PadBeg(queryNum);
            var top = new StringBuilder();
            var bot = new StringBuilder();
            for (int i = 0; i < aln.Blocks.Count; i++)
            {
                var y = aln.Blocks[i];
                if (i > 0)
                {
                    var x = aln.Blocks[i - 1];
                    for (uint p = x.End1; p < y.Beg1; p++)
                        top.Append((char)_alph.Decode[_db.Text[p]]);
                    top.Append('-', (int)(y.Beg2 > x.End2 ? y.Beg2 - x.End2 : 0));
                    bot.Append('-', (int)(y.Beg1 - x.End1));
                    for (uint p = x.End2; p < y.Beg2; p++)
                        bot.Append((char)_alph.Decode[queryAll[origin + (int)p]]);
                }

                for (uint p = y.Beg1; p < y.End1; p++)
                    top.Append((char)_alph.Decode[_db.Text[p]]);
                for (uint p = y.Beg2; p < y.End2; p++)
                    bot.Append((char)_alph.Decode[queryAll[origin + (int)p]]);
            }

            o.Append("s ").Append(n1.PadRight(nw)).Append(' ').Append(b1.PadLeft(bw)).Append(' ').Append(r1.PadLeft(rw)).Append(' ').Append('+').Append(' ').Append(s1.PadLeft(sw)).Append(' ').Append(top).Append('\n');
            o.Append("s ").Append(n2.PadRight(nw)).Append(' ').Append(b2.PadLeft(bw)).Append(' ').Append(r2.PadLeft(rw)).Append(' ').Append(_query.Strand(queryNum)).Append(' ').Append(s2.PadLeft(sw)).Append(' ').Append(bot).Append('\n');
            o.Append('\n');
        }
    }
}

/// <summary>The reader of FASTA records that <c>MultiSequence::appendFromFasta</c> implements: a record starts at '&gt;', the name is the first word of the line, the sequence is every character above the blank outside the header lines.</summary>
internal static class LastFasta
{
    public static IEnumerable<(string Name, string Residues)> Read(string text)
    {
        int i = 0;
        while (true)
        {
            while (i < text.Length && text[i] <= ' ')
                i++;
            if (i >= text.Length)
                yield break;
            if (text[i] != '>')
                throw new LastProgramException("bad FASTA sequence data: missing '>'");
            i++;
            int eol = text.IndexOf('\n', i);
            if (eol < 0)
                eol = text.Length;
            string header = text[i..eol];
            i = Math.Min(eol + 1, text.Length);
            string[] words = header.Split([' ', '\t', '\r', '\f', '\v'], StringSplitOptions.RemoveEmptyEntries);
            string name = words.Length > 0 ? words[0] : "";
            var sb = new StringBuilder();
            while (i < text.Length && text[i] != '>')
            {
                if (text[i] > ' ')
                    sb.Append(text[i]);
                i++;
            }

            yield return (name, sb.ToString());
        }
    }
}
