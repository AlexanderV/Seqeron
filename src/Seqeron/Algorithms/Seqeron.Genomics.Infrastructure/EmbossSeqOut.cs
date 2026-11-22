namespace Seqeron.Genomics.Infrastructure;

/// <summary>
/// A sequence as <c>seqClone</c> (<c>ajseqwrite.c</c>) leaves it in the sequence output object (<c>AjPSeqout</c>) that the writers print: the fields of the sequence object cut to its range, with the
/// accuracy values and the feature table. A5-8c3h, B09 F265.
/// </summary>
internal sealed class EmbossSeqOut
{
    public string Name = "", Acc = "", Sv = "", Gi = "", Desc = "", Db = "", Setdb = "", Tax = "", Taxcommon = "", Taxid = "", Organelle = "", Type = "", Seq = "";
    public List<string> Acclist = [], Keylist = [], Taxlist = [];

    /// <summary><c>Fulldesc-&gt;Multi</c>; null when the reader made no full description.</summary>
    public List<string>? FullDescMulti;

    /// <summary><c>outseq-&gt;Accuracy</c> (null when the output object never had any): NaN is memory that was never written.</summary>
    public float[]? Accuracy;

    /// <summary>The feature table (<c>outseq-&gt;Fttable</c>) the writers that keep the features in the sequence file print; the range, the reversal and the circular flag are applied by <c>seqClone</c>.</summary>
    public EmbossSeqFeatureTable? Fttable;
    public bool Circular, Rev;

    /// <summary><c>outseq-&gt;Offset</c>: the start of the range minus one.</summary>
    public int Offset;

    public bool IsProtein => Type.Length > 0 && Type[0] == 'P';
}

/// <summary>The state of the output object of a program run (<c>AjPSeqout</c>) and of the writers (the <c>static</c> variables of <c>ajseqwrite.c</c>) that lasts from one sequence to the next.</summary>
internal sealed class EmbossSeqWriteRun
{
    /// <summary>The counters of the names <c>gnl|BL_ORD_ID|N</c> of the ncbi and gifasta writers.</summary>
    public int NcbiOrdinal, GifastaOrdinal;

    /// <summary>The number of sequences written to the file that is open (<c>outseq-&gt;Count</c>), and to all the files of the run: the DAS writers write their header when it is 0.</summary>
    public int Count, HeaderCount;

    /// <summary>The date of the run (the header of mase).</summary>
    public DateTime RunDate = DateTime.Now;

    /// <summary>The name of the program and its arguments (<c>ajUtilGetProgram</c> / <c>ajUtilGetCmdline</c>: the <c>@PG</c> line of SAM files).</summary>
    public string Program = "seqret", CommandLine = "";

    /// <summary>The accuracy values of the output object (<c>outseq-&gt;Accuracy</c>, <c>Qualsize</c>): it keeps its buffer from one sequence to the next, <c>seqClone</c> and <c>seqDeclone</c> zero it.</summary>
    public float[]? OutAccuracy;
    public int OutQualsize;

    /// <summary>The features of the sequences are written to the sequence file (<c>seqoutUfoLocal</c>: the feature output is wanted and no feature file is named).</summary>
    public bool UfoLocal;

    /// <summary>The sequences saved for the formats that write a set (<c>outseq-&gt;Savelist</c>), and the name of the output file (<c>%F</c>).</summary>
    public List<EmbossSeqSaved> Saved = [];

    public string FileName = "";

    /// <summary><c>outseq-&gt;Type</c>: the type of the last sequence cloned (<c>seqClone</c> sets it when the sequence has one).</summary>
    public string OutType = "";

    /// <summary>The counter of <c>ajSeqSetNameMulti</c> (a <c>static</c>).</summary>
    public int NameMultiCount;

    /// <summary>The name of the database for the names of the sequences (<c>-osdbname</c>, <c>outseq-&gt;Setoutdb</c>).</summary>
    public string Setoutdb = "";
}

/// <summary>A sequence in the save list of the output object (<c>outseq-&gt;Savelist</c>, <c>seqWriteListAppend</c>): the formats that write a set of sequences keep the sequences until the file is closed.</summary>
internal sealed class EmbossSeqSaved
{
    public string Name = "", Desc = "", Type = "", Seq = "";
    public float Weight = 1.0f;

    /// <summary>The name of the first gene (<c>Genelist</c>); null when the sequence has none.</summary>
    public string? Gene;
}
