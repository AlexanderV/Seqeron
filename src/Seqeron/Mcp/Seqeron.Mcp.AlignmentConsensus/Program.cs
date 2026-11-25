// SPDX-License-Identifier: GPL-3.0-only
// This file is part of Seqeron.Mcp.AlignmentConsensus, licensed under the terms stated in LICENSE (this directory).
using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Hosting;
using Seqeron.Mcp.AlignmentConsensus.Tools;

var builder = Host.CreateApplicationBuilder(args);

builder.Services
    .AddMcpServer(options =>
    {
        options.ServerInfo = new()
        {
            Name = "Seqeron.Mcp.AlignmentConsensus",
            Version = "1.0.0"
        };
    })
    .WithStdioServerTransport()
    .WithTools<AlignmentConsensusTools>();

await builder.Build().RunAsync();
