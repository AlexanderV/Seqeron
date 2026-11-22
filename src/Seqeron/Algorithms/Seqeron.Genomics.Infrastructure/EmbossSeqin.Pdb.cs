#pragma warning disable S1121, S3358, S2325, S1066, S3776, S107, S1854, S1481, S2589, S3626 // the readers keep the control flow of ajseqread.c

using System.Text;

namespace Seqeron.Genomics.Infrastructure;

// seqReadPdb / seqReadPdbnuc (the residues of the ATOM records: the sequence observed in the structure) and seqReadPdbseq / seqReadPdbnucseq (the SEQRES records: the original sequence) of ajseqread.c.
internal sealed partial class EmbossSeqin
{
    private static readonly string[] BaseAaTable =
    [
        "ALA", "ASX", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "---", "LYS", "LEU", "MET", "ASN", "---", "PRO", "GLN", "ARG", "SER", "THR", "---", "VAL", "TRP", "XAA", "TYR", "GLX"
    ];

    private static readonly string[] BaseNucTable =
    [
        "DA", "DB", "DC", "DD", "--", "--", "DG", "DH", "--", "--", "DK", "--", "DM", "DN", "--", "--", "--", "DR", "DS", "DT", "DU", "DV", "DW", "DX", "DY", "--"
    ];

    // ajResidueFromTriplet / ajBaseFromDoublet: true when the code is recognised
    private static bool ResidueFromCode(string code, bool nuc, out char c)
    {
        var table = nuc ? BaseNucTable : BaseAaTable;
        for (int i = 0; i < 26; i++)
        {
            if (code == table[i])
            {
                c = (char)(i + 'A');
                return true;
            }
        }

        if (code == "UNK")
        {
            c = nuc ? 'N' : 'X';
            return true;
        }

        c = nuc ? 'N' : 'X';
        return false;
    }

    // ajStrAssignSubS (the positions are clamped to the string)
    private static string SubS(string s, int pos1, int pos2)
    {
        int ibegin = PosStart(s.Length, 0, pos1);
        int iend = PosStart(s.Length, ibegin, pos2);
        if (iend == s.Length)
            iend--;
        return s.Substring(ibegin, iend - ibegin + 1);
    }

    // ajStrKeepRange(&str, 0, 71)
    private static string KeepRange71(string s) => s.Length > 72 ? s[..72] : s;

    private static string TrimWhite(string s) => s.Trim(' ', '\t', '\n', '\r', '\f', '\v');

    private sealed class PdbData
    {
        public readonly Dictionary<string, AlnItem> Table = [];
        public string?[] Names = [];
        public int Nseq, Count;
    }

    private PdbData? _pdb;

    private bool PdbReturn(EmbossSeq t)
    {
        var d = _pdb!;
        string? name = d.Count < d.Names.Length ? d.Names[d.Count] : null;
        if (name is null || !d.Table.TryGetValue(name, out var item))
            throw new InvalidOperationException("The reference program dereferences a null sequence (the name of the entry is not set).");
        t.Name = name;
        t.Seq.Clear().Append(item.Seq.ToString());
        d.Count++;
        return true;
    }

    // seqReadPdb (nuc false) / seqReadPdbnuc (nuc true)
    private bool ReadPdbAtoms(EmbossSeq t, bool nuc)
    {
        if (_pdb is null)
        {
            _ = ReadLine();
            if (!PrefixC(_line, "HEADER    "))
            {
                StoreReset();
                return false;
            }

            string name = TrimWhite(SubS(_line, 62, 71));
            var data = new PdbData();
            _pdb = data;
            var list = new List<string>();
            bool ok = ReadLine();
            AlnItem? alnitem = null;
            string model = "", chain = "";
            long lastaa = 0;
            int nseq = 0;
            while (ok && _line != "END")
            {
                string line = _line;
                if (PrefixC(line, "MODEL"))
                    model = TrimWhite(SubS(line, 7, 14));
                else if (PrefixC(line, "ATOM"))
                {
                    alnitem ??= new AlnItem();
                    line = KeepRange71(line);
                    string aa3 = nuc ? SubS(line, 18, 19) : SubS(line, 17, 19);
                    chain = SubS(line, 21, 21);
                    CStr.ToUint(SubS(line, 22, 25), out uint iaa);
                    if (iaa > lastaa)
                    {
                        if (ResidueFromCode(aa3, nuc, out char aa))
                            alnitem.Seq.Append(aa);
                        lastaa = iaa;
                    }
                }
                else if (PrefixC(line, "TER"))
                {
                    if (nuc && alnitem is null)
                        throw new InvalidOperationException("The reference program dereferences a null sequence (a TER record with no residue before it).");
                    if (alnitem is not null && alnitem.Seq.Length == 0)
                        alnitem = null;
                    else
                    {
                        if (alnitem is null)
                            throw new InvalidOperationException("The reference program dereferences a null sequence (a TER record with no residue before it).");
                        nseq++;
                        string token = name + "_" + chain;
                        if (model.Length > 0)
                            token += model;
                        alnitem.Name = ItemName(token);
                        alnitem.Weight = 1.0f;
                        data.Table[alnitem.Name] = alnitem;
                        list.Add(alnitem.Name);
                        alnitem = null;
                    }

                    lastaa = 0;
                }

                ok = ReadLine();
            }

            if (nuc && alnitem is not null && alnitem.Seq.Length == 0 && list.Count > 0)
                list.RemoveAt(list.Count - 1);
            if (nseq == 0)
            {
                _pdb = null;
                StoreReset();
                return false;
            }

            data.Names = new string?[nseq];
            for (int i = 0; i < nseq && i < list.Count; i++)
                data.Names[i] = list[i];
            data.Nseq = nseq;
            data.Count = 0;
        }

        var d = _pdb;
        if (d.Count >= d.Nseq)
        {
            _buff.Clear(0);
            _pdb = null;
            return false;
        }

        return PdbReturn(t);
    }

    // seqReadPdbseq (nuc false) / seqReadPdbnucseq (nuc true)
    private bool ReadPdbSeqres(EmbossSeq t, bool nuc)
    {
        if (_pdb is not null && _pdb.Nseq > 0 && _pdb.Count >= _pdb.Nseq)
        {
            _buff.Clear(0);
            _pdb = null;
        }

        if (_pdb is null)
        {
            bool ok = ReadLine();
            if (!ok || !PrefixC(_line, "HEADER    "))
            {
                StoreReset();
                return false;
            }

            string name = TrimWhite(SubS(_line, 62, 71));
            var data = new PdbData();
            _pdb = data;
            var list = new List<string>();
            ok = ReadLine();
            AlnItem? alnitem = null;
            int nseq = 0;
            string chain = "";
            while (ok && _line != "END")
            {
                if (PrefixC(_line, "SEQRES"))
                {
                    string line = KeepRange71(_line);
                    var tok = new StrTok(line, " \n\r");
                    tok.Next(out _);
                    tok.Next(out string number);
                    CStr.ToUint(number, out uint iseq);
                    if (tok.Next(out string ch))
                        chain = ch;
                    if (iseq == 1)
                    {
                        if (alnitem is not null && alnitem.Seq.Length == 0)
                        {
                            nseq--;
                            RemovePdbItem(data, list, alnitem);
                            alnitem = null;
                        }

                        nseq++;
                        alnitem = new AlnItem { Name = ItemName(name + "_" + chain), Weight = 1.0f };
                        data.Table[alnitem.Name] = alnitem;
                        list.Add(alnitem.Name);
                    }

                    while (tok.Next(out string code))
                    {
                        if (ResidueFromCode(code, nuc, out char aa))
                        {
                            if (alnitem is null)
                                throw new InvalidOperationException("The reference program dereferences a null sequence (a SEQRES record of a chain that does not start at 1).");
                            alnitem.Seq.Append(aa);
                        }
                    }
                }

                ok = ReadLine();
            }

            if (alnitem is not null && alnitem.Seq.Length == 0)
            {
                nseq--;
                RemovePdbItem(data, list, alnitem);
            }

            if (nseq == 0)
            {
                _pdb = null;
                StoreReset();
                return false;
            }

            data.Names = new string?[nseq];
            for (int i = 0; i < nseq && i < list.Count; i++)
                data.Names[i] = list[i];
            data.Nseq = nseq;
            data.Count = 0;
        }

        return PdbReturn(t);
    }

    // the item of a chain with no residue is dropped: the last name of the list is popped and the key of the item removed from the table
    private static void RemovePdbItem(PdbData data, List<string> list, AlnItem item)
    {
        if (list.Count > 0)
            list.RemoveAt(list.Count - 1);
        data.Table.Remove(item.Name);
    }
}
