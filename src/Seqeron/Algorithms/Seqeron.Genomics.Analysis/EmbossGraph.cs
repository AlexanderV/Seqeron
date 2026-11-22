#pragma warning disable S3776, S107, S1854, S1066, S3358, S1125, S127, S1905, S1075, S1643
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>An object of a graph (<c>AjPGraphobj</c>): a line, a rectangle or a text.</summary>
internal sealed record EmbossGraphObject(string Kind, float X1, float Y1, float X2, float Y2, int Colour, string Text = "");

/// <summary>One set of points of an XY graph (<c>AjPGraphdata</c> of <c>ajgraph.c</c>).</summary>
internal sealed class EmbossGraphData
{
    public string Type = "", Dataname = "", Title = "", Subtitle = "", Xaxis = "", Yaxis = "";
    public float[] X = [], Y = [];
    public float TminX, TmaxX, TminY, TmaxY, MinX, MaxX, MinY, MaxY;
    public readonly List<EmbossGraphObject> Objects = [];

    /// <summary><c>ajGraphdataSetTruescale</c>: the true scale; a range of 0 width is widened.</summary>
    public void SetTruescale(float xmin, float xmax, float ymin, float ymax)
    {
        TminX = xmin;
        TminY = ymin;
        TmaxX = xmax;
        TmaxY = ymax;
        if (TminX == TmaxX)
        {
            if (TminX == 0)
                TmaxX = 0.5f;
            else if (TminX > 0)
            {
                TmaxX *= 1.1f;
                TminX /= 1.1f;
            }
            else
            {
                TmaxX /= 1.1f;
                TminX *= 1.1f;
            }
        }

        if (TminY == TmaxY)
        {
            if (TminY == 0)
                TmaxY = 0.5f;
            else if (TminY > 0)
            {
                TmaxY *= 1.1f;
                TminY /= 1.1f;
            }
            else
            {
                TmaxY /= 1.1f;
                TminY *= 1.1f;
            }
        }
    }

    /// <summary><c>ajGraphdataAddposLine</c>.</summary>
    public void AddLine(float x, float y, float x2, float y2, int colour) => Objects.Add(new EmbossGraphObject("Line", x, y, x2, y2, colour));
}

/// <summary>
/// The XY graph of an ACD <c>xygraph</c> item (<c>AjPGraph</c> of <c>ajgraph.c</c>, B09 WP63, F282): the device is checked as <c>ajGraphxySetDevicetype</c> does (an exact name or a unique prefix of one, in the order of the table of the library),
/// and the display of the devices that write text is ported: <c>data</c> (<c>GraphxyDisplayToData</c>), <c>bedgraph</c> (<c>GraphxyDisplayToBedgraph</c>) and the devices that draw nothing (<c>none</c>, <c>null</c>, <c>text</c>). The devices that draw with plplot
/// (<c>ps</c>, <c>png</c>, <c>svg</c>, ...) need the drivers of that library and are refused.
/// </summary>
internal sealed class EmbossGraphXy
{
    // the table graphType[] of ajgraph.c of this build (X11 and the png / pdf drivers present): name, device, shows in the list of the devices (it has an open function)
    private static readonly (string Name, string Device, bool Alias, bool Listed)[] Types =
    [
        ("ps", "ps", false, true), ("postscript", "ps", true, true), ("hpgl", "lj_hpgl", false, true), ("hp7470", "hp7470", false, true), ("hp7580", "hp7580", false, true), ("meta", "plmeta", false, true),
        ("cps", "psc", false, true), ("colourps", "psc", true, true), ("x11", "xwin", false, true), ("xwindows", "xwin", true, true), ("tek", "tek4107t", false, true), ("tek4107t", "tek4107t", true, true),
        ("tekt", "tekt", false, true), ("tektronics", "tekt", true, true), ("none", "null", false, true), ("null", "null", true, true), ("text", "null", true, true), ("data", "data", false, true),
        ("das", "das", false, false), ("wig", "wig", false, false), ("bedgraph", "bedgraph", false, false), ("xterm", "xterm", false, true), ("png", "png", false, true), ("gif", "gif", false, true),
        ("pdf", "pdf", false, true), ("svg", "svg", false, true),
    ];

    private readonly string _program;
    private readonly int _multi;
    private string _kind = "";

    /// <summary>The title, the description, the subtitle, the axis labels, the name of the data set and the base of the output file names (<c>ajGraphxyNewI</c>, <c>GraphxyNewPlplot</c>: the program name).</summary>
    public string Title = "", Desc = "", Subtitle = "", Dataname = "", OutputFile;

    /// <summary>The title of the graph as the ACD set it (<c>-gtitle</c>): every plot starts from it (<c>ajGraphSetTitleS</c> with the title the program took at the start).</summary>
    public string InitialTitle = "";

    /// <summary>The labels of the axes; null = never set (the data file prints <c>&lt;null&gt;</c>).</summary>
    public string? Xaxis, Yaxis;

    /// <summary>The graphs of the object (<c>ajGraphDataReplace</c> / <c>ajGraphDataReplaceI</c>).</summary>
    public readonly List<EmbossGraphData> Graphs = [];

    private readonly List<string> _created = [];
    private readonly Dictionary<string, string> _files = [];

    private EmbossGraphXy(string program, int multi)
    {
        _program = program;
        _multi = multi;
        OutputFile = program;
    }

    /// <summary>The files the graph wrote, by name.</summary>
    public IReadOnlyDictionary<string, string> Files => _files;

    /// <summary>The device of the graph: it draws with plplot (the program cannot be run).</summary>
    public bool NeedsPlplot => _kind is not ("data" or "bedgraph" or "none");

    /// <summary>
    /// <c>acdSetGraphxy</c> of a value: the device name must be one of the table (an exact name, else a unique prefix: <c>ajGraphxySetDevicetype</c>); an ambiguous prefix is an error of the library, an unknown or ambiguous name dumps the devices.
    /// </summary>
    /// <returns>The graph; null when the value is refused (the messages are in <paramref name="context"/>).</returns>
    public static EmbossGraphXy? Create(string program, string value, int multi, EmbossFeatureReadContext context)
    {
        int j = -1;
        for (int i = 0; i < Types.Length; i++)
        {
            string name = Types[i].Name;
            if (string.Equals(value, name, StringComparison.OrdinalIgnoreCase))
            {
                j = i;
                break;
            }

            if (value.Length <= name.Length && name.StartsWith(value, StringComparison.OrdinalIgnoreCase))
            {
                if (j < 0)
                {
                    j = i;
                }
                else
                {
                    var aliases = Types.Where(t => value.Length <= t.Name.Length && t.Name.StartsWith(value, StringComparison.OrdinalIgnoreCase)).Select(t => t.Name);
                    context.Error($"Ambiguous graph device name '{value}' ({string.Join(", ", aliases)})");
                    j = -2;
                    break;
                }
            }
        }

        if (j < 0)
        {
            DumpDevices(context);
            context.Error($"Invalid XY graph value '{value}'");
            return null;
        }

        return new EmbossGraphXy(program, Math.Max(1, multi)) { _kind = Types[j].Device is "null" ? "none" : Types[j].Name is "data" or "bedgraph" ? Types[j].Name : Types[j].Device };
    }

    // ajGraphicsDumpDevices: ajUserDumpC / ajUser print on the error stream without a prefix
    private static void DumpDevices(EmbossFeatureReadContext context)
    {
        context.Messages.Add("Devices allowed (with alternative names) are:-");
        foreach (var t in Types)
        {
            if (!t.Listed || t.Alias)
                continue;
            var aliases = Types.Where(o => o.Alias && o.Device == t.Device).Select(o => o.Name).ToList();
            context.Messages.Add(aliases.Count > 0 ? $"{t.Name} ({string.Join(", ", aliases)})" : t.Name);
        }
    }

    /// <summary><c>ajGraphSetOutdirS</c>: the directory is put in front of the output file name.</summary>
    public void SetOutdir(string dir)
    {
        if (dir.Length == 0)
            return;
        OutputFile = (dir.EndsWith('/') ? dir : dir + "/") + OutputFile;
    }

    /// <summary><c>ajGraphSetOutfileS</c>.</summary>
    public void SetOutfile(string name)
    {
        if (name.Length > 0)
            OutputFile = name;
    }

    /// <summary><c>ajGraphAppendTitleS</c>: a graph without a title is titled "description of what".</summary>
    public void AppendTitle(string what)
    {
        if (Title.Length > 0)
            return;
        Title = Desc.Length > 0 ? $"{Desc} of {what}" : $"{_program} of {what}";
    }

    /// <summary><c>ajGraphDataReplace</c>: the first graph.</summary>
    public void Replace(EmbossGraphData g)
    {
        Graphs.Clear();
        Graphs.Add(g);
    }

    /// <summary><c>ajGraphDataReplaceI</c>: the graph number <paramref name="num"/> (the next one when it is the number of the graphs).</summary>
    public bool Replace(EmbossGraphData g, int num)
    {
        if (Graphs.Count > num)
            Graphs[num] = g;
        else if (Graphs.Count < _multi)
            Graphs.Add(g);
        else
            return false;
        return true;
    }

    /// <summary>
    /// <c>ajGraphxyDisplay</c> (<paramref name="closeit"/>: the call closes the output, which prints the files created: a display that could not open its file returns before that). The standard output of the closing is <see cref="Stdout"/>.
    /// </summary>
    /// <exception cref="InvalidOperationException">The device draws with plplot, or the display is one the reference program does not define.</exception>
    public void Display(EmbossFeatureReadContext context, bool closeit = false)
    {
        bool ok;
        switch (_kind)
        {
            case "data":
                ok = ToData(context);
                break;
            case "bedgraph":
                ok = ToBedgraph(context);
                if (ok && closeit && Graphs.Count > 1)
                    throw new InvalidOperationException("GraphxyDisplayToBedgraph pushes the one file name of the display once for every graph on the list of the files created and GraphClose deletes it every time: a string is used after it was freed (valgrind: invalid read): the reference program has no defined result for several graphs.");
                break;
            case "none":
                // GraphxyDisplayXwin on the null device draws the objects: the first filled rectangle makes plplot warn that the driver has no hardware fills
                if (Graphs.Any(g => g.Objects.Any(o => o.Kind == "Filled Rectangle")) && !_warnedFill)
                {
                    _warnedFill = true;
                    context.Messages.Add("\n*** PLPLOT WARNING ***\nDriver does not support hardware solid fills, switching to software fill.\n");
                }

                ok = true;
                break;
            default:
                throw new InvalidOperationException($"The graph device '{_kind}' draws with plplot (the graph is not ported: use -graph data, bedgraph or none).");
        }

        if (ok && closeit)
            Stdout += Close();
    }

    private bool _warnedFill;

    /// <summary>What closing the graph printed on the standard output (<c>Created name</c>).</summary>
    public string Stdout { get; private set; } = "";

    // ajFmtPrintF("%f"): a float as a double with 6 decimals
    private static string F(float v) => CFormat.Fixed(v, 6);

    private bool Open(string name, EmbossFeatureReadContext context, out StringBuilder? sb)
    {
        sb = new StringBuilder();
        string? slash = name.Contains('/', StringComparison.Ordinal) ? name[..name.LastIndexOf('/')] : null;
        if (slash is not null && !(ExistingDirectories?.Contains(slash) ?? true))
        {
            context.Warn($"Failed to open output file '{name}' error:2 'No such file or directory'");
            context.Error($"Could not open graph file {name}\n");
            sb = null;
            return false;
        }

        return true;
    }

    /// <summary>The directories that exist where the program runs; null = all.</summary>
    public IReadOnlyCollection<string>? ExistingDirectories { get; set; }

    // GraphxyDisplayToData
    private bool ToData(EmbossFeatureReadContext context)
    {
        float minxa = float.MaxValue, minya = float.MaxValue, maxxa = -float.MaxValue, maxya = -float.MaxValue;
        foreach (var g in Graphs)
        {
            minxa = minxa < g.MinX ? minxa : g.MinX;
            minya = minya < g.MinY ? minya : g.MinY;
            maxxa = maxxa > g.MaxX ? maxxa : g.MaxX;
            maxya = maxya > g.MaxY ? maxya : g.MaxY;
        }

        for (int i = 0; i < Graphs.Count; i++)
        {
            var g = Graphs[i];
            string temp = OutputFile + (i + 1).ToString(System.Globalization.CultureInfo.InvariantCulture) + ".dat";
            _created.Add(temp);
            if (!Open(temp, context, out var o))
                return false;
            if (Title.Length <= 1)
                Title = string.Concat(Title, _program);
            o!.Append("##").Append(g.Type).Append('\n');
            o.Append("##Title ").Append(Title).Append('\n');
            o.Append("##Graphs ").Append(Graphs.Count).Append('\n');
            o.Append("##Number ").Append(i + 1).Append('\n');
            o.Append("##Points ").Append(g.X.Length).Append('\n');
            o.Append("##XminA ").Append(F(minxa)).Append(" XmaxA ").Append(F(maxxa)).Append(" YminA ").Append(F(minya)).Append(" YmaxA ").Append(F(maxya)).Append('\n');
            o.Append("##Xmin ").Append(F(g.TminX)).Append(" Xmax ").Append(F(g.TmaxX)).Append(" Ymin ").Append(F(g.TminY)).Append(" Ymax ").Append(F(g.TmaxY)).Append('\n');
            o.Append("##ScaleXmin ").Append(F(g.MinX)).Append(" ScaleXmax ").Append(F(g.MaxX)).Append(" ScaleYmin ").Append(F(g.MinY)).Append(" ScaleYmax ").Append(F(g.MaxY)).Append('\n');
            o.Append("##Maintitle ").Append(g.Title).Append('\n');
            if (Graphs.Count == 1)
            {
                o.Append("##Xtitle ").Append(Xaxis ?? "<null>").Append('\n');
                o.Append("##Ytitle ").Append(Yaxis ?? "<null>").Append('\n');
            }
            else
            {
                o.Append("##Subtitle ").Append(g.Subtitle).Append('\n');
                o.Append("##Xtitle ").Append(g.Xaxis).Append('\n');
                o.Append("##Ytitle ").Append(g.Yaxis).Append('\n');
            }

            for (int j = 0; j < g.X.Length; j++)
                o.Append(F(g.X[j])).Append('\t').Append(F(g.Y[j])).Append('\n');
            o.Append("##DataObjects\n##Number ").Append(g.Objects.Count).Append('\n');
            foreach (var obj in g.Objects)
                Obj(o, obj);
            o.Append("##GraphObjects\n##Number 0\n");
            _files[temp] = o.ToString();
        }

        return true;
    }

    private static void Obj(StringBuilder o, EmbossGraphObject p)
    {
        if (p.Kind == "Text")
        {
            o.Append("Textline x1 ").Append(F(p.X1)).Append(" y1 ").Append(F(p.Y1)).Append(" x2 ").Append(F(p.X2)).Append(" y2 ").Append(F(p.Y2)).Append(" colour ").Append(p.Colour).Append(" size 1.0 ").Append(p.Text).Append('\n');
            return;
        }

        o.Append(p.Kind).Append(" x1 ").Append(F(p.X1)).Append(" y1 ").Append(F(p.Y1)).Append(" x2 ").Append(F(p.X2)).Append(" y2 ").Append(F(p.Y2)).Append(" colour ").Append(p.Colour).Append('\n');
    }

    // GraphxyDisplayToBedgraph
    private bool ToBedgraph(EmbossFeatureReadContext context)
    {
        string temp = OutputFile + ".bedgraph";
        if (!Open(temp, context, out var o))
        {
            _created.Add(temp);
            return false;
        }

        o!.Append("track type=bedGraph graphType=points\n");
        foreach (var g in Graphs)
        {
            _created.Add(temp);
            for (int j = 0; j < g.X.Length; j++)
            {
                int ipos = (int)g.X[j] - 1;
                string name = g.Dataname.Length > 0 ? g.Dataname : Dataname.Length > 0 ? Dataname : "unknown";
                o.Append(name).Append(' ').Append(ipos).Append(' ').Append(ipos).Append(' ').Append(F(g.Y[j])).Append('\n');
            }
        }

        _files[temp] = o.ToString();
        return true;
    }

    /// <summary><c>ajGraphicsClose</c>: "Created name" for every file the displays pushed on the list (<c>ajFmtPrint</c>: the standard output).</summary>
    public string Close()
    {
        var sb = new StringBuilder();
        foreach (string f in _created)
            sb.Append("Created ").Append(f).Append('\n');
        _created.Clear();
        return sb.ToString();
    }
}
