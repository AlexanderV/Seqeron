using System.Text.RegularExpressions;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// The <c>seqall</c> input of an EMBOSS program (<c>acdSetSeqall</c> / <c>ajSeqAllRead</c> / <c>ajSeqallNext</c> of the AJAX library): the USA (<c>[format::]file</c>, <c>asis::SEQUENCE</c>), the reading of every sequence of the file by the
/// sequence readers of <see cref="EmbossSeqin"/> and the messages of the library. A5-8c3g, B09 F248.
/// </summary>
internal sealed class EmbossSeqall
{
    private static readonly Regex RegAsis = new("^[Aa][Ss][Ii][Ss]:+([^:].*)$", RegexOptions.CultureInvariant | RegexOptions.Singleline);
    private static readonly Regex RegFormat = new("^([A-Za-z0-9-]*)::(.*)$", RegexOptions.CultureInvariant | RegexOptions.Singleline);

    private readonly EmbossSeqin _input;
    private EmbossSeq? _first;
    private bool _firstReturned;

    private EmbossSeqall(EmbossSeqin input, EmbossSeq first)
    {
        _input = input;
        _first = first;
    }

    /// <summary>
    /// Opens the input and reads its first sequence (<c>ajSeqAllRead</c> of the ACD). Null: the program dies ("Unable to read sequence", the messages of the library are in <paramref name="context"/>; the caller adds the <c>Died:</c> line).
    /// </summary>
    /// <exception cref="InvalidOperationException">The input is one on which the reference program never stops or dies by a fault.</exception>
    public static EmbossSeqall? Open(string program, string spec, IReadOnlyDictionary<string, string> files, EmbossSequenceQualifiers q, bool features, EmbossFeatureReadContext context, string sequenceType = "gapany")
    {
        context.Program = program;
        string? prefix = null;
        string file = spec.Trim(' ', '\t', '\n');
        string? text;
        var am = RegAsis.Match(file);
        var pm = RegFormat.Match(file);
        if (am.Success)
        {
            // asis::SEQUENCE: the text itself is the sequence, read with the format asis
            prefix = "asis";
            text = am.Groups[1].Value;
            file = "asis";
        }
        else
        {
            if (pm.Success)
            {
                prefix = pm.Groups[1].Value;
                file = pm.Groups[2].Value;
            }

            if (!files.TryGetValue(file, out text))
            {
                context.Error($"Failed to open filename '{file}'");
                context.Error($"Unable to read sequence '{spec}'");
                return null;
            }
        }

        var options = new EmbossSeqinOptions { Format = q.Format, PrefixFormat = prefix, Filename = file, Features = features, UfoGiven = !string.IsNullOrEmpty(q.FeaturesUfo), Lower = q.Lower, Upper = q.Upper, IsNuc = q.Nucleotide, IsProt = q.Protein, SequenceType = sequenceType };
        var input = new EmbossSeqin(text, options, context, name => files.TryGetValue(name, out string? t) ? t : null);
        var first = input.Next();
        if (first is null)
        {
            if (prefix == "asis")
            {
                // an asis:: text with no residue: the library drops the warning about the empty sequence and fails to open the "file"
                context.Messages.RemoveAll(m => m.StartsWith("Warning: Sequence 'asis::", StringComparison.Ordinal) && m.EndsWith("has zero length, ignored", StringComparison.Ordinal));
                context.Error("Failed to open filename 'asis'");
            }

            context.Error($"Unable to read sequence '{spec}'");
            return null;
        }

        return new EmbossSeqall(input, first);
    }

    /// <summary>The first sequence (the one the ACD read when it opened the input) without taking it.</summary>
    public EmbossSeq Peek() => _first ?? throw new InvalidOperationException("The first sequence was taken.");

    /// <summary><c>ajSeqallNext</c>: the next sequence, or null (the library reads it, and prints its messages, when it is asked).</summary>
    public EmbossSeq? Next()
    {
        if (!_firstReturned)
        {
            _firstReturned = true;
            var f = _first;
            _first = null;
            return f;
        }

        return _input.Next();
    }
}
