// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Infrastructure;

/// <summary>What <c>ajSeqParseNcbi</c> stores of a FASTA header with an NCBI-style identifier (<c>sp|P12345|NAME_HUMAN desc</c>).</summary>
/// <param name="Id">The identifier (the sequence name before <c>seqSetName</c>); empty when the program leaves the sequence unnamed.</param>
/// <param name="Accession">The accession, or empty.</param>
/// <param name="SequenceVersion">The accession.version, or empty.</param>
/// <param name="Gi">The GI number, or empty.</param>
/// <param name="Db">The database prefix (<c>sp</c>, <c>gb</c>, <c>lcl</c>, ...), or empty.</param>
/// <param name="Description">The description (the rest of the header line).</param>
public sealed record EmbossNcbiHeader(string Id, string Accession, string SequenceVersion, string Gi, string Db, string Description)
{
    /// <summary>
    /// <c>ajSeqParseNcbi</c> (EMBOSS 6.6.0 <c>ajax/core/ajseqread.c</c>), the parser of the <c>fasta</c> / <c>ncbi</c> input format. Returns null when the header is not an NCBI-style header (no <c>|</c> in the
    /// first word, or a PIR header with <c>;</c> at position 3: the plain <see cref="EmbossFastaReader"/> branch applies). The program tokenises with <c>ajStrToken*</c>: leading delimiters are skipped when a string is
    /// assigned and runs of delimiters are skipped after every token.
    /// </summary>
    internal static EmbossNcbiHeader? Parse(string header)
    {
        if (header.Length > 3 && header[3] == ';')
            return null;
        if (header.Length == 0 || header[0] != '>')
            return null;
        var idTok = new Tokenizer(header, "> \t\r\n");
        string idstr = idTok.Next(out _) ? idTok.Last : "";
        idTok.Delimiters = "\r\n";
        string reststr = idTok.Next(out _) ? idTok.Last : "";
        if (idstr.Length == 0 || !idstr.Contains('|'))
            return null;
        return Core(idstr, reststr);
    }

    private static EmbossNcbiHeader Core(string idstr, string reststr)
    {
        string db = "", gi = "", acc = "", sv = "", id = "";
        var t = new Tokenizer(idstr, "|");
        t.Next(out _);
        string prefix = t.Last;
        if (prefix.StartsWith("gi", StringComparison.Ordinal))
        {
            t.Next(out _);
            gi = t.Last;
            if (!t.Next(out _))
                return new EmbossNcbiHeader(gi, "", "", gi, "", reststr);
            prefix = t.Last;
        }

        if (idstr.StartsWith("gnl|BL_ORD_ID|", StringComparison.Ordinal))
        {
            // formatdb: the rest of the line is the real header (parsed again); numtoken is never set by the program, so a failed parse leaves the id empty
            string inner = ">" + reststr;
            var rec = ParseAny(inner);
            return rec with { Db = rec.Db.Length == 0 ? "BL_ORD_ID" : rec.Db };
        }

        switch (prefix)
        {
            case "bbs":
            case "lcl":
                if (prefix == "lcl")
                    db = prefix;
                t.Next(out _);
                return new EmbossNcbiHeader(t.Last, "", "", gi, db, reststr);
            case "gnl":
            case "pat":
                if (prefix == "gnl")
                {
                    t.Next(out _);
                    db = t.Last;
                }
                else
                {
                    t.Step();
                }

                t.Next(out _);
                return new EmbossNcbiHeader(t.Last, "", "", gi, db, reststr);
            case "pdb":
            {
                db = prefix;
                t.Next(out _);
                id = t.Last;
                if (t.Next(out _))
                    id += t.Last;
                return new EmbossNcbiHeader(id, "", "", gi, db, reststr);
            }

            case "gb":
            case "emb":
            case "dbj":
            case "tpd":
            case "tpg":
            case "sp":
            case "ref":
            {
                db = prefix;
                t.Next(out _);
                string token = t.Last;
                string? vacc = EmbossFastaReader.SeqVersionAccession(token);
                if (vacc is not null)
                {
                    sv = token;
                    acc = vacc;
                }
                else if (EmbossFastaReader.IsAccession(token))
                {
                    acc = token;
                }

                id = t.Next(out _) ? t.Last : token;
                return new EmbossNcbiHeader(id, acc, sv, gi, db, reststr);
            }

            case "pir":
            case "prf":
                db = prefix;
                t.Next(out _);
                return new EmbossNcbiHeader(t.Last, "", "", gi, db, reststr);
        }

        // else: the last two or three barred tokens are [db|]acc|id
        var counter = new Tokenizer(idstr, "|");
        uint nt = 0;
        while (counter.Step())
            nt++;
        if (idstr[^1] == '|')
            nt++;
        if (nt < 2)
            return new EmbossNcbiHeader("", "", "", gi, "", reststr);   // not reachable: the id holds a '|'

        var u = new Tokenizer(idstr, "|");
        // for (i = 0; i < nt - 3; ++i) with unsigned arithmetic: nt = 2 wraps to 4294967295 steps, which consume every token (the program needs ~28 s for it)
        uint skips = unchecked(nt - 3);
        for (uint i = 0; i < skips; i++)
        {
            if (!u.Step())
                break;   // the token is reset (empty); the remaining steps of the program's loop change nothing
        }

        u.Next(out _);
        db = u.Last;
        u.Next(out _);
        string tok = u.Last;
        string? v = EmbossFastaReader.SeqVersionAccession(tok);
        if (v is not null)
        {
            sv = tok;
            acc = v;
            id = v;
        }
        else if (EmbossFastaReader.IsAccession(tok))
        {
            acc = tok;
            id = tok;
        }
        else
        {
            id = tok;
        }

        u.Delimiters = " \n\t\r";
        if (u.Next(out _) && u.Last.Length > 0)
            id = u.Last;
        return new EmbossNcbiHeader(id, acc, sv, gi, db, reststr);
    }

    // the recursion of gnl|BL_ORD_ID: ajSeqParseNcbi on the rest of the line; a header without '|' goes through ajSeqParseFasta
    private static EmbossNcbiHeader ParseAny(string header)
    {
        var ncbi = Parse(header);
        if (ncbi is not null)
            return ncbi;
        var (name, accession, description) = EmbossFastaReader.ParsePlainRaw(header);
        return new EmbossNcbiHeader(name, accession, "", "", "", description);
    }

    // ajStrToken: Pos starts after the leading delimiters; Next returns the token at Pos and skips the delimiters behind it; at the end the token is reset (empty string)
    private sealed class Tokenizer
    {
        private string _s;
        private int _pos;

        public Tokenizer(string s, string delimiters)
        {
            _s = s;
            Delimiters = delimiters;
            while (_pos < _s.Length && Delimiters.Contains(_s[_pos]))
                _pos++;
        }

        public string Delimiters { get; set; }

        public string Last { get; private set; } = "";

        private void Reset()
        {
            _s = "";
            _pos = 0;
        }

        public bool Next(out string token)
        {
            token = Last = "";
            if (_pos >= _s.Length)
            {
                Reset();
                return false;
            }

            int start = _pos;
            while (_pos < _s.Length && !Delimiters.Contains(_s[_pos]))
                _pos++;
            token = Last = _s[start.._pos];
            while (_pos < _s.Length && Delimiters.Contains(_s[_pos]))
                _pos++;
            return true;
        }

        public bool Step()
        {
            if (_pos >= _s.Length)
            {
                Reset();
                return false;
            }

            while (_pos < _s.Length && !Delimiters.Contains(_s[_pos]))
                _pos++;
            while (_pos < _s.Length && Delimiters.Contains(_s[_pos]))
                _pos++;
            return true;
        }
    }
}
