#pragma warning disable S1121, S3358, S2325, S1066, S3776, S107, S1854, S1481, S2589, S3626 // the readers keep the control flow of ajseqread.c

namespace Seqeron.Genomics.Infrastructure;

// The accuracy values of the sequence object of the stream (thys->Accuracy, thys->Qualsize of AjPSeq): the object (and so its buffer) lives as long as the stream, ajSeqClear zeroes the first Qualsize values, AJCRESIZE
// (realloc) grows or shrinks the buffer to the size asked for (the new values are zero); a value read beyond the buffer is a NaN here: a writer that reads one has no defined output.
internal sealed partial class EmbossSeqin
{
    private float[]? _acc;
    private int _qualsize;

    // the tables of the classes of ajseq.c are made by the first call of a class function (ajSeqclsSetGb), and the sequence object keeps its date object once a reader made one (ajSeqClear clears it, it is not freed)
    private bool _clsTables, _dateObject;

    // a sequence that comes after one that had a date has an empty date object
    private void DateObject(EmbossSeq t)
    {
        if (_dateObject)
            t.Date = new EmbossSeqDate { VersionsExist = true };
    }

    // ajSeqClear
    private void AccuracyClear()
    {
        if (_acc is not null && _qualsize > 0)
            Array.Clear(_acc, 0, Math.Min(_qualsize, _acc.Length));
    }

    // AJCRESIZE(thys->Accuracy, n)
    private void AccuracyResize(int n)
    {
        if (_acc is null)
            _acc = new float[n];
        else if (n != _acc.Length)
            Array.Resize(ref _acc, n);
        // the library's AJCRESIZE leaves the new values zero here (valgrind reports no use of uninitialised values on these runs)
    }

    // the Accuracy of the sequence that is returned (a copy: the buffer of the stream changes with the next sequence)
    private void AccuracySnapshot(EmbossSeq t)
    {
        t.Accuracy = _acc is null ? null : (float[])_acc.Clone();
        t.Qualsize = _qualsize;
    }

    // the value of a quality character: FASTQ-sanger (phred, qmin 33), FASTQ-illumina (qmin 64), FASTQ-solexa (qmin 59) and SAM (phred)
    private static float QualityValue(int c, string kind) => kind switch
    {
        "SOLEXA" => (float)EmbossQuality.Solexa[c],
        "ILLUMINA" => c >= 64 ? c - 64 : 0,
        _ => c >= 33 ? c - 33 : 0,
    };
}
