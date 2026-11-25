// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossSequenceFormats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.EmbossSequenceFormats.Tests;

/// <summary>B09 WP67 (F309): the MCP tools of the data resource programs drfinddata, dbtell and showdb; the expected text is what the installed EMBOSS 6.6.0 program writes for the example of the tool document.</summary>
[TestFixture]
public class EmbossDataProgramToolsTests
{
    [Test]
    public void Drfinddata_FindsTheResourcesOfAnEdamTerm()
    {
        var r = EmbossSequenceFormatsTools.Drfinddata("microarray", subclasses: false, outputFormat: "list");
        Assert.That(r.Text, Is.EqualTo("drcat:CleanEx\ndrcat:ArrayExpress\n"));
        Assert.That(EmbossSequenceFormatsTools.Drfinddata("microarray", subclasses: false, outputFormat: "basic").Text, Does.StartWith("ID      CleanEx\nName    "));
    }

    [Test]
    public void Dbtell_WritesTheDefinitionOfAStandardDatabaseAndTheLineOfADrcatEntry()
    {
        Assert.That(EmbossSequenceFormatsTools.Dbtell("go_process").Text, Is.EqualTo("# go_process is defined in standard\n# access levels id: Yes query: Yes all: Yes\n\nDBNAME go_process [\n   format:         \"obo\"\n   method:         \"emboss\"\n   type:           \"obo\"\n   comment:        \"Gene Ontology (biological processes)\"\n   directory:      \"/usr/share/EMBOSS/data/OBO\"\n   field:          \"id ! identifier without the prefix\"\n   field:          \"acc ! full name and any alternate identifier(s)\"\n   field:          \"nam ! words in the name\"\n   field:          \"isa ! parent identifier from is_a relation(s)\"\n   field:          \"des ! words in the description\"\n   field:          \"ns ! namespace\"\n   indexdirectory: \"/usr/share/EMBOSS/index\"\n   namespace:      \"biological_process\"\n   release:        \"1.1.3300\"\n   edamdat:        \"0966 ! Ontology term\"\n   edamfmt:        \"2196 ! OBO format\"\n   edamid:         \"1176 ! GO term ID\"\n   edamtpc:        \"0089 ! Ontologies\"\n   taxon:          \"1 all\"\n]\n"));
        Assert.That(EmbossSequenceFormatsTools.Dbtell("EMBL").Text, Is.EqualTo("DBNAME EMBL [\n   comment: \"defined in DRCAT\"\n]\n"));
        Assert.That(EmbossSequenceFormatsTools.Dbtell("nosuchdatabase").Text, Is.Empty);
    }

    [Test]
    public void Showdb_WritesTheTableOfTheDatabases()
    {
        Assert.That(EmbossSequenceFormatsTools.Showdb("taxon").Text, Is.EqualTo("# Name         Type     Comment\n# ============ ======== =======\ntaxon          Taxonomy NCBI taxonomy\n"));
        Assert.That(EmbossSequenceFormatsTools.Showdb(qualifiers: "-only -type -noheading -noobo -noresource").Text, Is.EqualTo("taxon           Taxonomy \n"));
    }

    [Test]
    public void Tools_RefuseWhatTheLibraryRefuses()
    {
        Assert.Throws<ArgumentException>(() => EmbossSequenceFormatsTools.Drfinddata(""));
        var e = Assert.Throws<ArgumentException>(() => EmbossSequenceFormatsTools.Drfinddata("microarray", outputFormat: "html"));
        Assert.That(e!.Message, Does.Contain("Unknown output format 'html'"));
        var f = Assert.Throws<ArgumentException>(() => EmbossSequenceFormatsTools.Showdb("nosuchdatabase"));
        Assert.That(f!.Message, Does.Contain("The database 'nosuchdatabase' does not exist"));
    }
}
