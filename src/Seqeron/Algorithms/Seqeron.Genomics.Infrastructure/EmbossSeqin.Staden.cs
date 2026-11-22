#pragma warning disable S1121, S3358, S2325, S1066, S3776, S107, S1854, S1481, S2589, S3626 // the readers keep the control flow of ajseqread.c

using System.Text;

namespace Seqeron.Genomics.Infrastructure;

// seqReadExperiment (the Staden experiment file) and seqReadAce (the consed ACE format) of ajseqread.c.
internal sealed partial class EmbossSeqin
{
    // seqReadExperiment: the lines of an entry that starts with an ID line (EN / EX experiment, AV accuracy values, AC accession, SV version, DE description, SQ and the sequence up to //); the accuracy values are checked
    // against the sequence (the library warns about missing and bad values)
    private bool ReadExperiment(EmbossSeq t)
    {
        if (!ReadLine())
            return false;
        Sw lineprefix = PrefixSwiss(_line);
        if (lineprefix != Sw.Id)
        {
            StoreReset();
            return false;
        }

        AssignText();
        var h = new StrTok(_line, " \n\r\t");
        h.Next(out _);
        h.Next(out string token);
        SetName(t, token);
        bool ok = ReadLine();
        lineprefix = PrefixSwiss(_line);
        string accval = "";
        string desc = "";
        while (ok && lineprefix != Sw.Sq)
        {
            if (lineprefix == Sw.Ex)
            {
                var a = new StrTok(_line, " ");
                a.Next(out _);
                a.Next(out token, "\n\r");
                desc = desc.Length > 0 ? new StringBuilder(desc).Append(' ').Append(token).ToString() : token;
            }

            if (lineprefix == Sw.Av)
            {
                var a = new StrTok(_line, " ");
                a.Next(out _);
                a.Next(out token, "\n\r");
                accval = accval.Length > 0 ? new StringBuilder(accval).Append(' ').Append(token).ToString() : token;
            }

            if (lineprefix == Sw.Ac)
            {
                var a = new StrTok(_line, " ;\n\r");
                a.Next(out _);
                while (a.Next(out token))
                    AccSave(t, token);
            }

            if (lineprefix == Sw.Sv)
            {
                var a = new StrTok(_line, " \n\r");
                a.Next(out _);
                a.Next(out token);
                SvSave(t, token);
            }

            if (lineprefix == Sw.De)
            {
                var a = new StrTok(_line, " ");
                a.Next(out _);
                a.Next(out token, "\n\r");
                desc = desc.Length > 0 ? new StringBuilder(desc).Append(' ').Append(token).ToString() : token;
            }

            ok = ReadLine();
            lineprefix = PrefixSwiss(_line);
        }

        t.Desc = desc;
        ok = ReadLine();
        lineprefix = PrefixSwiss(_line);
        while (ok && lineprefix != Sw.End)
        {
            Append(t, _line);
            ok = ReadLine();
            lineprefix = PrefixSwiss(_line);
        }

        while (ok && lineprefix != Sw.Id)
        {
            ok = ReadLine();
            lineprefix = PrefixSwiss(_line);
        }

        if (ok)
            StoreClear(1);
        else
            _buff.Clear(0);
        if (accval.Length > 0)
        {
            int ilen = t.Seq.Length;
            var tokens = new StrTok(accval, " ");
            for (int i = 0; i < ilen; i++)
            {
                if (!tokens.Next(out string value))
                {
                    _context.Warn($"Missing accuracy for base {i + 1} in experiment format\n");
                    break;
                }

                var inner = new StrTok(value, ",");
                while (inner.Next(out string part))
                {
                    if (!CStr.ToInt(part, out _))
                    {
                        _context.Warn($"Bad accuracy '{value}' for base {i + 1} in experiment format\n");
                        break;
                    }
                }
            }
        }

        return true;
    }

    // seqReadAce: the library reads the line "AS contigs reads" and then takes the names of the contigs from a table it has not filled (Names is a null pointer: the process dies by a segmentation fault)
    private bool ReadAce(EmbossSeq t)
    {
        _ = t;
        if (!ReadLine())
            return false;
        var h = new StrTok(_line, " \n\r");
        h.Next(out string token, " \t");
        if (!string.Equals(token, "AS", StringComparison.OrdinalIgnoreCase))
        {
            StoreReset();
            return false;
        }

        throw new InvalidOperationException("The ACE reader of the library dereferences a null table of names (the reference program dies by a segmentation fault on every ACE file).");
    }
}
