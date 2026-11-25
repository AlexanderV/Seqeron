// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.EmbossSequenceFormats, licensed under the terms stated in LICENSE (this directory).
using System.ComponentModel;
using System.Text;
using ModelContextProtocol.Server;
using Seqeron.Mcp.Analysis.Tools;
using static Seqeron.Mcp.Analysis.Tools.AnalysisTools;

namespace Seqeron.Mcp.EmbossSequenceFormats.Tools;

public partial class EmbossSequenceFormatsTools
{
    [McpServerTool(Name = "drfinddata", Title = "Data — EMBOSS drfinddata (find data resources by the EDAM data type of their output)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 drfinddata: the data resources of the catalogue DRCAT whose EDAM data terms match the keywords of query (comma separated; the identifier, the accession or a word of the name or a synonym of a term of the ontology edam; * is a wildcard): every term found and, with subclasses (default), all its descendants (is_a); every catalogue entry that lists such a term as EDAMdat is written once, in the order the indexes of the installation return them, in the outputFormat drcat (the entries as in DRCAT.dat), basic (identifier, name and URL) or list (drcat:IDENTIFIER). The EDAM ontology and the DRCAT catalogue are those of the EMBOSS 6.6.0 installation, embedded in the library. Equals the installed program on 233 runs (standard output, error stream, exit status and written files byte for byte); an output directory, or an output file in a directory that does not exist (the reference program dies by a segmentation fault), is refused. B09 WP67, F309.")]
    public static EmbossProgramToolResult Drfinddata(
        [Description("The keywords, comma separated: an EDAM identifier (1300 or EDAM_data:1300), a word of a term name (transcript) or a pattern with * (Al*); upper or lower case.")] string query,
        [Description("Also find the descendants of the terms (-subclasses, default true; false = -nosubclasses).")] bool subclasses = true,
        [Description("The output format (-oformat): drcat (default), basic or list.")] string? outputFormat = null,
        CancellationToken cancellationToken = default)
    {
        if (string.IsNullOrWhiteSpace(query))
            throw new ArgumentException("The query cannot be empty", nameof(query));
        cancellationToken.ThrowIfCancellationRequested();
        var args = new List<string> { "-query", query };
        if (!subclasses)
            args.Add("-nosubclasses");
        if (!string.IsNullOrEmpty(outputFormat))
            args.AddRange(["-oformat", outputFormat]);
        args.AddRange(["-outfile", "stdout"]);
        return RunEmbossMultiProgram("drfinddata", args, new Dictionary<string, string>(), nameof(query), "stdout");
    }
    [McpServerTool(Name = "dbtell", Title = "Data — EMBOSS dbtell (the definition of a database of the EMBOSS configuration)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 dbtell: the definition of the database as the standard configuration of the installation holds it (the ontologies chebi, eco, edam, edam_data, edam_format, edam_identifier, edam_operation, edam_topic, go, go_component, go_function, go_process, pw, ro, so, swo, the data resource catalogue drcat and the taxonomy taxon): the file that defines it, the access levels (id, query, all), the attributes; an alias (any case other than the lower case name) is reported as such; a name that only the DRCAT catalogue defines (embl, uniprot, pdb ...) gives the line 'defined in DRCAT'; any other name gives no output. The databases of a user configuration file (.embossrc) and of a server (-server) are not read. Equals the installed program on 302 runs (standard output, error stream, exit status and written files byte for byte). B09 WP67, F309.")]
    public static EmbossProgramToolResult Dbtell(
        [Description("The name of the database (any case), e.g. go_process, drcat, EMBL.")] string database,
        [Description("Verbose output (-full): the aliases of the database follow its definition.")] bool full = false,
        CancellationToken cancellationToken = default)
    {
        cancellationToken.ThrowIfCancellationRequested();
        var args = new List<string> { "-database", database ?? "" };
        if (full)
            args.Add("-full");
        args.AddRange(["-outfile", "stdout"]);
        return RunEmbossMultiProgram("dbtell", args, new Dictionary<string, string>(), nameof(database), "stdout");
    }
    [McpServerTool(Name = "showdb", Title = "Data — EMBOSS showdb (the databases of the EMBOSS configuration)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 showdb: the table of the databases of the standard configuration of the installation (taxon, drcat and the ontologies), grouped by type and sorted by definition file and name, with the columns selected by qualifiers (-full all columns; -only with -type -id -query -all -comment -methods -taxscope -fields -numfields -aliases -numaliases -examples -numexamples -defined -release -heading -access to choose the columns; -noprotein -nonucleic -nosequence -nofeature -notext -notaxonomy -noresource -noobo -noxml and -nostandard -nouser -noinclude to leave databases out; -html for an HTML table); with database only that database (a name the configuration does not define stops the program). The databases of a user configuration file (.embossrc) are not read. Equals the installed program on 301 runs (standard output, error stream, exit status and written files byte for byte). B09 WP67, F309.")]
    public static EmbossProgramToolResult Showdb(
        [Description("The name of one database to give information on (-database); default: all.")] string? database = null,
        [Description("Qualifiers as typed on the command line, e.g. \"-full -nodefined\" or \"-only -type -comment -noheading\".")] string? qualifiers = null,
        CancellationToken cancellationToken = default)
    {
        cancellationToken.ThrowIfCancellationRequested();
        var args = new List<string>();
        if (!string.IsNullOrEmpty(database))
            args.AddRange(["-database", database]);
        AddQualifierTokens(args, qualifiers);
        args.AddRange(["-outfile", "stdout"]);
        return RunEmbossMultiProgram("showdb", args, new Dictionary<string, string>(), nameof(database), "stdout");
    }
}
