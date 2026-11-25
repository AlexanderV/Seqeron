// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.EmbossSequenceFormats, licensed under the terms stated in LICENSE (this directory).
using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Hosting;
using Seqeron.Mcp.EmbossSequenceFormats.Tools;

Seqeron.Mcp.EmbossSequenceFormats.Tools.EmbossSequenceFormatsTools.InstallAlignmentFormats();

var builder = Host.CreateApplicationBuilder(args);

builder.Services
    .AddMcpServer(options =>
    {
        options.ServerInfo = new()
        {
            Name = "Seqeron.Mcp.EmbossSequenceFormats",
            Version = "1.0.0"
        };
    })
    .WithStdioServerTransport()
    .WithTools<EmbossSequenceFormatsTools>();

await builder.Build().RunAsync();
