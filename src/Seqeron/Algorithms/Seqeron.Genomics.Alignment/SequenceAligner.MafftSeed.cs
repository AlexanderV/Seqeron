#pragma warning disable S3626, S3871, S1643, S2325, IDE0052, IDE0051, IDE0059, S127, S3358, S1066, S4487, S1144, CA1823, CS0414, CA1870, S1125, S3776, S107, S138, S1541, S1854
using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: --seed / --seedtable (multi2hat3s, the hat3 table with the blocks of kind 'k', the structure weights of tbfast -J and dvtditr)

    /// <summary>The result of <c>multi2hat3s</c>: its standard error, the sequences it writes (<c>&gt;_seed_name</c>, gaps kept) and the file <c>hat3</c>.</summary>
    internal sealed record Multi2Hat3Result(string Log, string Rows, string Hat3, MafftProgramFailureException? Failure);

    internal sealed partial class MafftRun
    {
        /// <summary>
        /// <c>multi2hat3s -t nhomologs -o seedoffset -i file</c> (<c>multi2hat3s.c</c>): the seed alignment is read (the alphabet guessed from the letters of the file itself), written as <c>&gt;_seed_name</c> with its gaps, and the ungapped blocks of every pair
        /// of its rows are written to <c>hat3</c> as blocks of kind 'k' with the weight <c>nhomologs^2 * 100</c> and the sequence numbers shifted by <paramref name="seedoffset"/>.
        /// </summary>
        internal static Multi2Hat3Result Multi2Hat3s(string seedText, int nhomologs, int seedoffset, string fileName)
        {
            var g = new MafftGlobals { Nthread = 1, Ppenalty = MafftGlobals.NotSpecified, PpenaltyEx = MafftGlobals.NotSpecified, Poffset = MafftGlobals.NotSpecified, PpenaltyDist = MafftGlobals.NotSpecified, KimuraR = MafftGlobals.NotSpecified, PamN = MafftGlobals.NotSpecified, Alg = 'A', Outgap = 1 };
            g.Report("nhomologs = " + nhomologs.ToString(CultureInfo.InvariantCulture) + "\nseedoffset = " + seedoffset.ToString(CultureInfo.InvariantCulture) + "\nseed = " + fileName + "\n");
            var run = new MafftRun(g);
            var recs = MafftProgramIo.Records(seedText);
            int njob = recs.Count;
            if (njob < 2)
                return new Multi2Hat3Result(g.Log.ToString(), "", "", new MafftProgramFailureException("At least 2 sequences should be input!\nOnly " + njob.ToString(CultureInfo.InvariantCulture) + " sequence found.\n"));
            var plain = recs.Select(r => MafftEngine.PlainResidues(r.Raw)).ToList();
            bool nuc = MafftEngine.GuessNucleotide(plain, new MafftOptions());
            g.Dorp = nuc ? 'd' : 'p';
            var seq = plain.Select(x => nuc ? x.ToLowerInvariant() : x.ToUpperInvariant()).ToArray();
            run._njob = njob;
            run._nlenmax = seq.Max(x => x.Length);
            g.Constants(seq);
            g.Report("Gap Penalty = " + MafftGlobals.F((double)g.Ppenalty / 1000, 5, 2, true) + ", " + MafftGlobals.F((double)g.PpenaltyEx / 1000, 5, 2, true) + ", " + MafftGlobals.F((double)g.Poffset / 1000, 5, 2, true) + "\n");
            if (run.CheckAlphabet(seq, "Illeagal") is { } unknown)
                return new Multi2Hat3Result(g.Log.ToString(), "", "", new MafftProgramFailureException(unknown.Text));
            var rows = new StringBuilder();
            for (int i = 0; i < njob; i++)
                rows.Append(">_seed_").Append(recs[i].Name).Append('\n').Append(seq[i]).Append('\n');
            if (nhomologs < 1)
                nhomologs = 1;
            double tsuyosa = (double)nhomologs * nhomologs * 100;
            g.Report("tsuyosa = " + MafftGlobals.F(tsuyosa, 0, 6) + "\n");
            var hat3 = new StringBuilder();
            g.Report("\n");
            for (int i = 0; i < njob - 1; i++)
            {
                for (int j = i + 1; j < njob; j++)
                {
                    if (seq[i].Length != seq[j].Length)
                        return new Multi2Hat3Result(g.Log.ToString(), rows.ToString(), hat3.ToString(), new MafftProgramFailureException("## ERROR  ###\nNot aligned,  =" + recs[i].Name + " - =" + recs[j].Name + "\n## ERROR  ###\n"));
                    g.Report("adding " + i.ToString(CultureInfo.InvariantCulture) + "-" + j.ToString(CultureInfo.InvariantCulture) + "\r");
                    var head = new MafftLocalHom();
                    run.PutLocalHom2(seq[i], seq[j], head, 0, 0, 'k');
                    for (var t = head; t != null; t = t.Next)
                    {
                        if (t.Opt == -1.0 || t.Start1 == -1)
                            continue;
                        hat3.Append((i + seedoffset).ToString(CultureInfo.InvariantCulture)).Append(' ').Append((j + seedoffset).ToString(CultureInfo.InvariantCulture)).Append(' ').Append(t.Overlapaa.ToString(CultureInfo.InvariantCulture)).Append(' ')
                            .Append(MafftGlobals.F(t.Opt * tsuyosa, 6, 3)).Append(' ').Append(t.Start1.ToString(CultureInfo.InvariantCulture)).Append(' ').Append(t.End1.ToString(CultureInfo.InvariantCulture)).Append(' ')
                            .Append(t.Start2.ToString(CultureInfo.InvariantCulture)).Append(' ').Append(t.End2.ToString(CultureInfo.InvariantCulture)).Append(" k\n");
                    }
                }
            }
            g.Report("\n");
            g.Report("multi2hat3s (" + g.KindName + ") Version 7.505\nalg=" + g.Alg + ", model=" + g.ModelName + ", amax=" + g.SpecificityConsideration.ToString("F1", CultureInfo.InvariantCulture) + "\n" + g.Nthread.ToString(CultureInfo.InvariantCulture) + " thread(s)\n\n");
            return new Multi2Hat3Result(g.Log.ToString(), rows.ToString(), hat3.ToString(), null);
        }

        /// <summary>
        /// <c>readlocalhomtable2_half</c> (<c>io.c</c>): the lines "i j overlapaa opt start1 end1 start2 end2 kind" of the file hat3 (i &lt; j, 0-based) become the chain of blocks of the pair; <paramref name="kozo"/> marks the sequences that a block of
        /// kind 'k' (a structure) names. The scores are rescaled as <c>opt / 5.8 * 600</c>. A line with fewer fields (the program then reuses the values of the line before) and a number outside the table (the program writes outside its table) are refused.
        /// </summary>
        internal static MafftLocalHom[][] ReadLocalHomTable2Half(string hat3Text, int njob, bool[] kozo, MafftLocalHom[][]? existing = null, bool afterFirstCharacter = false)
        {
            // the blocks are appended to the chains that the table holds already (the pair stage of tbfast; the seed table is read after it), a fresh table otherwise
            var table = existing ?? new MafftLocalHom[njob][];
            var nokori = new int[njob][];
            var last = new MafftLocalHom[njob][];
            for (int i = 0; i < njob; i++)
            {
                if (existing == null)
                {
                    table[i] = new MafftLocalHom[njob - i];
                    for (int j = 0; j < table[i].Length; j++)
                        table[i][j] = new MafftLocalHom();
                }
                nokori[i] = new int[njob - i];
                last[i] = new MafftLocalHom[njob - i];
                for (int j = 0; j < last[i].Length; j++)
                {
                    var head = table[i][j];
                    var tail = head;
                    int count = head.Opt == -1.0 && head.Next == null ? 0 : 1;
                    for (var t = head.Next; t != null; t = t.Next)
                    {
                        tail = t;
                        count++;
                    }
                    nokori[i][j] = count;
                    last[i][j] = tail;
                }
            }
            foreach (string buff in SplitFgetsLines(hat3Text, 256))
            {
                var sc = new CScanner(buff);
                if (!sc.Int(out int i) || !sc.Int(out int j) || !sc.Int(out int overlapaa) || !sc.Double(out double opt) || !sc.Int(out int start1) || !sc.Int(out int end1) || !sc.Int(out int start2) || !sc.Int(out int end2) || !sc.Word(out string infor))
                    throw new NotSupportedException("A line of the file hat3 / of the seed table has fewer than nine fields: the program reads the variables of the line before (undefined on the first line).");
                if (i < 0 || i >= njob || j >= njob && i < njob)
                {
                    if (i >= njob || i < 0)
                        throw new MafftProgramFailureException("Check hat3.  The first sequence must be younger than the second one.\n");
                    if (!afterFirstCharacter)
                        throw new NotSupportedException("The seed table names a sequence number that is not in the alignment: the program writes behind the table of the pair (undefined behaviour).");
                }
                if (infor[0] == 'k')
                {
                    // the flags of the sequences (a vector of the length of the whole input; a number behind it is written outside the vector, which the program does not read again)
                    if (i < kozo.Length)
                        kozo[i] = true;
                    if (j < kozo.Length)
                        kozo[j] = true;
                }
                if (j <= i)
                    throw new MafftProgramFailureException("Check hat3.  The first sequence must be younger than the second one.\n");
                if (j >= njob)
                    continue; // nodepair branch: the first character of the file was read before (the first line of the file is read from its second character: "0 1 20 ..." gives i = 1, j = 20), the block goes behind the row of its pair
                MafftLocalHom node;
                if (nokori[i][j - i]++ > 0)
                {
                    node = new MafftLocalHom { Importance = 0.0 };
                    last[i][j - i].Next = node;
                    last[i][j - i] = node;
                }
                else
                {
                    node = table[i][j - i];
                }
                node.Start1 = start1;
                node.Start2 = start2;
                node.End1 = end1;
                node.End2 = end2;
                node.Opt = (opt + 0.00) / 5.8 * 600;
                node.Overlapaa = overlapaa;
                node.Korh = infor[0];
            }
            return table;
        }
    }

    #endregion
}
