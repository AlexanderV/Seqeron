#pragma warning disable S1121, S3358, S2325, S2234 // keeps the control flow of ajfeatread.c
using System.Xml;
using System.Xml.Linq;

namespace Seqeron.Genomics.Infrastructure;

// featReadDasgff of ajfeatread.c: a DAS GFF (version 1.5 / 1.6) XML document read with the DOM of ajax. A5-8c3e, B09 F241.
internal sealed partial class EmbossSeqFeatureTable
{
    private static string? Attr(XElement? e, string name) => e?.Attribute(name)?.Value;

    // ajDomElementGetText: the text of the FIRST child when it is a text node; null otherwise
    private static string? ElementText(XElement? e)
    {
        if (e is null)
            return null;
        return e.FirstNode is XText t ? t.Value : null;
    }

    private static XElement? FirstChild(XElement e, string name) => e.Elements(name).FirstOrDefault();

    /// <summary><c>featReadDasgff</c>: false when the rest of the file is no XML document, has no SEGMENT with an id, or a FEATURE is incomplete.</summary>
    internal bool ReadDasgff(EmbossFeatureLines file)
    {
        var sb = new System.Text.StringBuilder();
        string? line;
        while ((line = file.Next()) is not null)
            sb.Append(line).Append('\n');
        XDocument doc;
        try
        {
            using var reader = XmlReader.Create(new StringReader(sb.ToString()), new XmlReaderSettings { DtdProcessing = DtdProcessing.Ignore });
            doc = XDocument.Load(reader);
        }
        catch (XmlException)
        {
            return false;
        }

        var segments = doc.Descendants("SEGMENT").ToList();
        if (segments.Count == 0)
            return false;
        string id = Attr(segments[0], "id") ?? "";
        if (id.Length == 0)
            return false;
        SeqId = id;
        CStr.ToUint(Attr(segments[0], "start"), out Start);
        Start = unchecked(Start - 1);
        CStr.ToUint(Attr(segments[0], "stop"), out End);
        End = unchecked(End - 1);
        Length = unchecked((int)(End - Start));
        bool ret = true;
        foreach (var das in doc.Descendants("FEATURE"))
        {
            CStr.ToUint(ElementText(FirstChild(das, "START")), out uint start);
            CStr.ToUint(ElementText(FirstChild(das, "END")), out uint end);
            var feature = NewII(unchecked((int)start), unchecked((int)end));
            string fid = Attr(das, "id") ?? "";
            if (fid.Length == 0)
            {
                ret = false;
                break;
            }

            GfftagAdd(feature, "ID", fid);
            var type = FirstChild(das, "TYPE");
            if (type is null)
            {
                ret = false;
                break;
            }

            feature.Type = Attr(type, "id") ?? "";
            string cv = Attr(type, "cvId") ?? "";
            if (cv.Length > 0)
                GfftagAdd(feature, "Ontology_term", cv);
            string cat = Attr(type, "category") ?? "";
            if (cat.Length > 0)
                GfftagAdd(feature, "category", cat);
            var method = FirstChild(das, "METHOD");
            if (method is null)
            {
                ret = false;
                break;
            }

            feature.Source = Attr(method, "id") ?? "";
            string mcv = Attr(method, "cvId") ?? "";
            if (mcv.Length > 0)
                feature.Source = mcv;
            var score = FirstChild(das, "SCORE");
            if (score is not null)
            {
                CStr.ToFloat(ElementText(score), out float sc);
                feature.Score = sc;
            }

            var orient = FirstChild(das, "ORIENTATION");
            if (orient is not null)
            {
                string t = ElementText(orient) ?? "";
                feature.Strand = t.Length > 0 ? t[0] : '\0';
            }

            var phase = FirstChild(das, "PHASE");
            if (phase is not null)
            {
                string t = ElementText(phase) ?? "";
                if (t.Length > 0 && t[0] == '-')
                    feature.Frame = 0;
                CStr.ToInt(t, out int fr);
                feature.Frame = fr;
            }

            foreach (var note in das.Descendants("NOTE"))
            {
                string? text = ElementText(note);
                if (text is not null)
                    GfftagAdd(feature, "Note", text);
            }

            foreach (var link in das.Descendants("LINK"))
                GfftagAdd(feature, "dasgff_link", Attr(link, "href") ?? "");
            foreach (var target in das.Descendants("TARGET"))
            {
                CStr.ToUint(Attr(target, "start"), out uint ts);
                CStr.ToUint(Attr(target, "stop"), out uint te);
                string tid = Attr(target, "id") ?? "";
                if (te != 0)
                    tid += " " + ts + " " + te;
                GfftagAdd(feature, "Target", tid);
            }

            bool bad = false;
            foreach (var group in das.Descendants("GROUP"))
            {
                string gid = Attr(group, "id") ?? "";
                if (gid.Length == 0)
                {
                    bad = true;
                    break;
                }

                string gtype = Attr(group, "type") ?? "";
                if (gtype.Length > 0)
                    gid += " type=" + gtype;
                string glabel = Attr(group, "label") ?? "";
                if (glabel.Length > 0)
                    gid += " label=" + glabel;
                GfftagAdd(feature, "group", gid);
            }

            if (bad)
                ret = false;
            foreach (var parent in das.Descendants("PARENT"))
            {
                string pid = Attr(parent, "id") ?? "";
                if (pid.Length == 0)
                {
                    ret = false;
                    break;
                }

                GfftagAdd(feature, "Parent", pid);
            }

            foreach (var part in das.Descendants("PART"))
            {
                string pid = Attr(part, "id") ?? "";
                if (pid.Length == 0)
                {
                    ret = false;
                    break;
                }

                TagAdd(feature, "part", pid);
            }
        }

        return ret;
    }

    // ajFeatNewII: a feature of the type misc_feature from Start to End (a start beyond the end is the reverse strand)
    private EmbossSeqFeature NewII(int start, int end)
    {
        if (start > end)
            return NewNuc("", "misc_feature", end, start, 0f, '-', 0, 0, 0, 0, null, null, 0, null)!;
        return NewNuc("", "misc_feature", start, end, 0f, '\0', 0, 0, 0, 0, null, null, 0, null)!;
    }
}
