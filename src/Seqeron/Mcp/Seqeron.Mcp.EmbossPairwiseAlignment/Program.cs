// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.EmbossPairwiseAlignment, licensed under the terms stated in LICENSE (this directory).
using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Hosting;
using Seqeron.Mcp.EmbossPairwiseAlignment.Tools;

Seqeron.Mcp.EmbossSequenceFormats.Tools.EmbossSequenceFormatsTools.InstallAlignmentFormats();

var builder = Host.CreateApplicationBuilder(args);

builder.Services
    .AddMcpServer(options =>
    {
        options.ServerInfo = new()
        {
            Name = "Seqeron.Mcp.EmbossPairwiseAlignment",
            Version = "1.0.0"
        };
    })
    .WithStdioServerTransport()
    .WithTools<EmbossPairwiseAlignmentTools>();

await builder.Build().RunAsync();
