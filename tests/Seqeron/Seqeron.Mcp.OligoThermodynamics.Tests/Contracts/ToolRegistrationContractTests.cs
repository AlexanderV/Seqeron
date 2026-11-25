// SPDX-License-Identifier: GPL-2.0-or-later
// Part of Seqeron.Mcp.OligoThermodynamics (Primer3-derived primer / probe design); licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the server project).
using Seqeron.Mcp.OligoThermodynamics.Tools;
using Seqeron.Mcp.Tests.Contracts;
using NUnit.Framework;

namespace Seqeron.Mcp.OligoThermodynamics.Tests;

/// <summary>Pins the Seqeron.Mcp.OligoThermodynamics server's public MCP tool surface (assertions in <see cref="ToolRegistrationContractTestsBase"/>).</summary>
[TestFixture]
public sealed class ToolRegistrationContractTests : ToolRegistrationContractTestsBase
{
    protected override Type ToolsType => typeof(OligoThermodynamicsTools);

    protected override string[] ExpectedToolNames { get; } =
    {
        "analyze_oligo",
        "design_antisense_probes",
        "design_molecular_beacon",
        "design_primers",
        "design_probes",
        "design_probes_primer3",
        "design_tiling_probes",
        "evaluate_primer",
        "generate_primer_candidates",
        "primer_dimer",
        "validate_probe",
    };
}
