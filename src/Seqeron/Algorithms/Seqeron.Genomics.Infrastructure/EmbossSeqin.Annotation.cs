#pragma warning disable S1121, S3358, S2325, S1066, S3776, S107, S1854, S1481, S2589, S3626, S3241 // the readers keep the control flow of ajseqread.c

namespace Seqeron.Genomics.Infrastructure;

// The annotation lines of the EMBL and Swiss-Prot entries (KW, OS, OC, OG, CC, DR, RN, RG, RX, RP, RA, RT, RL, RC, DT, XX) that seqReadEmbl and seqReadSwiss keep in the sequence object: keywords, taxonomy, comments,
// cross-references, citations and dates, which the entry writers print. A5-8c3h, B09 F268.
internal sealed partial class EmbossSeqin
{
    private sealed class AnnotState
    {
        public EmbossSeqRef? Ref;
        public string Cmt = "";
        public string Datestr = "", Relstr = "";
    }

    // ajSeqAddKey / ajSeqAddRef ...; the line is one of the annotation prefixes: true when it was taken
    private bool AnnotationLine(Sw lineprefix, EmbossSeq t, AnnotState st)
    {
        string token;
        switch (lineprefix)
        {
            case Sw.Kw:
            {
                var a = new StrTok(_line, " \n\r");
                a.Next(out _);
                while (a.Next(out token, ".;\n\r"))
                    t.Keylist.Add(token.Trim(' ', '\t', '\n', '\r', '\f', '\v'));
                return true;
            }

            case Sw.Os:
            {
                var a = new StrTok(_line, " \n\r");
                a.Next(out _);
                while (a.Next(out token, ".;\n\r"))
                {
                    token = token.Trim(' ', '\t', '\n', '\r', '\f', '\v');
                    var b = new StrTok(token, "()");
                    int itaxtype = 1;
                    while (b.Next(out string token2))
                    {
                        TaxSave(t, token2.Trim(' ', '\t', '\n', '\r', '\f', '\v'), itaxtype);
                        itaxtype = 3;
                    }
                }

                return true;
            }

            case Sw.Oc:
            case Sw.Og:
            {
                var a = new StrTok(_line, " \n\r");
                a.Next(out _);
                while (a.Next(out token, ".;\n\r"))
                    TaxSave(t, token.Trim(' ', '\t', '\n', '\r', '\f', '\v'), lineprefix == Sw.Oc ? 0 : 2);
                return true;
            }

            case Sw.Cc:
            {
                var a = new StrTok(_line, " ");
                a.Next(out _);
                a.Next(out token, "\n\r");
                if (st.Cmt.Length > 0)
                    st.Cmt += "\n";
                st.Cmt += token;
                return true;
            }

            case Sw.Dr:
            {
                var x = new EmbossSeqXref();
                var a = new StrTok(_line, " ;\n\r");
                a.Next(out _);
                a.Next(out token, ";\n\r");
                x.Db = token.Trim(' ', '\t', '\n', '\r', '\f', '\v');
                a.Next(out token);
                x.Id = token.Trim(' ', '\t', '\n', '\r', '\f', '\v');
                a.Next(out token);
                if (token.Length == 0)
                {
                    if (x.Id.EndsWith('.'))
                        x.Id = x.Id[..^1];
                }
                else
                {
                    if (token.EndsWith('.'))
                        token = token[..^1];
                    x.Secid = token.Trim(' ', '\t', '\n', '\r', '\f', '\v');
                    a.Next(out token);
                    if (token.Length == 0)
                    {
                        if (x.Secid.EndsWith('.'))
                            x.Secid = x.Secid[..^1];
                    }
                    else
                    {
                        if (token.EndsWith('.'))
                            token = token[..^1];
                        x.Terid = token.Trim(' ', '\t', '\n', '\r', '\f', '\v');
                        a.Next(out token);
                        if (token.Length == 0)
                        {
                            if (x.Terid.EndsWith('.'))
                                x.Terid = x.Terid[..^1];
                        }
                        else
                        {
                            if (token.EndsWith('.'))
                                token = token[..^1];
                            x.Quatid = token.Trim(' ', '\t', '\n', '\r', '\f', '\v');
                        }
                    }
                }

                x.Type = 1;
                t.Xreflist.Add(x);
                return true;
            }

            case Sw.Rn:
            {
                if (st.Ref is not null)
                {
                    st.Ref.Standard();
                    t.Reflist.Add(st.Ref);
                }

                st.Ref = new EmbossSeqRef();
                var a = new StrTok(_line, " ");
                a.Next(out _);
                a.Next(out token, "\n\r");
                string num = token.Length >= 2 ? token[1..^1] : "";
                CStr.ToUint(num, out uint refnum);
                st.Ref.Number = refnum;
                return true;
            }

            case Sw.Rg:
            case Sw.Rx:
            case Sw.Rp:
            case Sw.Ra:
            case Sw.Rt:
            case Sw.Rl:
            case Sw.Rc:
            {
                st.Ref ??= new EmbossSeqRef();
                var a = new StrTok(_line, " ");
                a.Next(out _);
                a.Next(out token, lineprefix == Sw.Ra ? "\n\r;" : "\n\r");
                switch (lineprefix)
                {
                    case Sw.Rg: st.Ref.AppendGroupname(token); break;
                    case Sw.Rx: st.Ref.AppendXref(token); break;
                    case Sw.Rp: st.Ref.AppendPosition(token); break;
                    case Sw.Ra: st.Ref.AppendAuthors(token); break;
                    case Sw.Rt:
                        if (token != ";")
                            st.Ref.AppendTitle(token);
                        break;
                    case Sw.Rl: st.Ref.AppendLocation(token); break;
                    default: st.Ref.AppendComment(token); break;
                }

                return true;
            }

            case Sw.Dt:
            {
                t.Date ??= new EmbossSeqDate();
                var a = new StrTok(_line, " (),");
                int icount = 0;
                while (a.Next(out token))
                {
                    icount++;
                    if (icount == 2)
                        st.Datestr = token;
                    else if (icount == 4)
                        st.Relstr = token;
                    else if (icount == 5)
                    {
                        if (token == "Created")
                        {
                            t.Date.SetCreate(st.Datestr);
                            t.Date.CreRel = st.Relstr;
                        }
                    }
                    else if (icount == 8)
                    {
                        t.Date.SetModify(st.Datestr);
                        t.Date.ModRel = st.Relstr;
                        t.Date.ModVer = token;
                    }
                }

                return true;
            }

            case Sw.Xx:
                if (st.Ref is not null)
                {
                    st.Ref.Standard();
                    t.Reflist.Add(st.Ref);
                    st.Ref = null;
                }

                if (st.Cmt.Length > 0)
                {
                    t.Cmtlist.Add(st.Cmt);
                    st.Cmt = "";
                }

                return true;
            default:
                return false;
        }
    }

    // ajSeqreflistGetXrefs: the cross-references of the RX lines (database and identifier) are added to the list of the cross-references (type RX); the list is longer than the DR lines it holds
    private static void RefXrefs(EmbossSeq t)
    {
        foreach (var r in t.Reflist)
        {
            if (r.Xref.Length == 0)
                continue;
            var handle = new StrTok(r.Xref, " =;\r\n");
            while (handle.Next(out string dbtoken, " =;\r\n"))
            {
                handle.Next(out string token, " ;\r\n");
                if (token.EndsWith('.'))
                    token = token[..^1];
                if (dbtoken.Length == 0 || token.Length == 0)
                    continue;
                string db = dbtoken;
                foreach (string known in new[] { "MEDLINE", "PubMed", "DOI", "AGRICOLA" })
                {
                    if (string.Equals(dbtoken, known, StringComparison.OrdinalIgnoreCase))
                        db = known;
                }

                t.Xreflist.Add(new EmbossSeqXref { Type = 6, Db = db, Id = token });
            }
        }
    }
}
