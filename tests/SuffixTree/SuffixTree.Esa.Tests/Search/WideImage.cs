using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.Fuzz;
using SuffixTree.Esa.Tests.Views;

namespace SuffixTree.Esa.Tests.Search;

/// <summary>
/// The sections of a text at one layout (work package S2.4): the stored width (4 or 5 bytes), the kind of the LCP array (U8, U16, U32, U64) and the CHILD codec (the tagged one exists over 4 bytes only), each in memory that ends against an
/// inaccessible page (<see cref="GuardedImage"/>, <see cref="ViewImage"/>). <see cref="Sections"/> are what the long descent reads; <see cref="Index"/> is an <see cref="EsaIndex"/> over the same sections with no container (the header
/// reader refuses width 5 and LCP codec 4 until S2.7), so that every operation of the index runs over a 5-byte suffix array and a 64-bit LCP array at a size where the in-memory index is the oracle.
/// </summary>
internal sealed unsafe class WideImage : IDisposable
{
    /// <summary>The layouts of the matrix: width, LCP kind, CHILD codec (the tagged one over 4 bytes only). 12 rows.</summary>
    internal static readonly (int Width, LcpKind Kind, bool Tagged)[] Layouts =
    [
        (4, LcpKind.U8, false), (4, LcpKind.U8, true), (4, LcpKind.U16, false), (4, LcpKind.U16, true),
        (4, LcpKind.U32, false), (4, LcpKind.U32, true), (4, LcpKind.U64, false), (4, LcpKind.U64, true),
        (5, LcpKind.U8, false), (5, LcpKind.U16, false), (5, LcpKind.U32, false), (5, LcpKind.U64, false),
    ];

    private readonly ViewImage _views;
    private readonly GuardedImage _text;
    private readonly byte[] _packed;
    private readonly LcpChildCase _case;
    private readonly int _width;
    private readonly LcpKind _kind;
    private readonly bool _tagged;

    internal WideImage(LcpChildCase c, int width, LcpKind kind, bool tagged)
    {
        _case = c;
        _width = width;
        _kind = kind;
        _tagged = tagged;
        _views = ViewImage.Of(c, width, kind, tagged);
        _packed = new byte[PackedTextFormat.TotalBytes(c.N, c.Bits)];
        PackedTextFormat.PackInto(c.Symbols, c.Bits, _packed);
        _text = new GuardedImage(_packed, flush: true);
    }

    /// <summary>The sections for the long descent.</summary>
    internal WideSections Sections(bool isChecked)
    {
        CompactLcp lcp = Lcp(isChecked);
        return new WideSections(isChecked, _width, _kind, _tagged, _views.SaPointer, _views.ChildPointer, _views.LcpBasePointer, lcp, _text.Data, _case.Bits, _case.N);
    }

    private CompactLcp Lcp(bool isChecked) =>
        CompactLcp.FromMappedWide(_kind, _case.N + 1, _views.LcpBasePointer, _views.PairsPointer, _views.Pairs.Count, _width, isChecked);

    /// <summary>
    /// The index over these sections: <see cref="EsaIndex.FromViews"/> with the statistics of the in-memory index of the text (they are what the header of a container would say). <paramref name="psi"/> and <paramref name="rmq"/> give it the PSI and
    /// RMQ sections at the width of the layout (the navigation pack then reads them in place); without them the pack builds the inverse suffix array and the pyramid.
    /// </summary>
    internal EsaIndex Index(bool isChecked, char[] alphabet, bool psi = false, bool rmq = false)
    {
        CompactLcp lcp = Lcp(isChecked);
        var reader = MappedTextReader.FromMemory(_text.Data, _packed.Length, _case.N, _case.Bits);
        var text = EsaText.FromPacked(reader, alphabet, _case.N);
        SaView sa = SaView.FromMapped(_views.SaPointer, _case.N + 1, _width, isChecked);
        ChildTable child = ChildTable.FromMapped(MappedChild.FromMappedWide(_views.ChildPointer, _case.N + 1, _width, lcp, isChecked, _tagged));
        var search = new MappedSearch(WideSearcher.Create(new WideSections(isChecked, _width, _kind, _tagged, _views.SaPointer, _views.ChildPointer, _views.LcpBasePointer, lcp, _text.Data, _case.Bits, _case.N)), lcp, _case.Bits);
        SaView? psiView = psi ? SaView.FromMapped(_views.PsiPointer, _case.N + 1, _width, isChecked) : null;
        SaView? rmqView = rmq && _views.RmqCount > 0 ? SaView.FromMapped(_views.RmqPointer, _views.RmqCount, _width, isChecked) : null;
        NavSections? navigation = psi || rmq ? new NavSections(psiView, rmqView, SesaFormat.RmqFanout, isChecked, 0) : null;
        return EsaIndex.FromViews(text, sa, lcp, child, _case.Index.Statistics, search, navigation);
    }

    public void Dispose()
    {
        _views.Dispose();
        _text.Dispose();
    }
}
