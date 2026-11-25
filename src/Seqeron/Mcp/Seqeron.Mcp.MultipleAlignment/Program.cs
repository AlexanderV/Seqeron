// SPDX-License-Identifier: GPL-3.0-or-later
// This file is part of Seqeron.Mcp.MultipleAlignment, licensed under the terms stated in LICENSE (this directory).
using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Hosting;
using Seqeron.Mcp.MultipleAlignment.Tools;

Seqeron.Mcp.EmbossSequenceFormats.Tools.EmbossSequenceFormatsTools.InstallAlignmentFormats();

var builder = Host.CreateApplicationBuilder(args);

builder.Services
    .AddMcpServer(options =>
    {
        options.ServerInfo = new()
        {
            Name = "Seqeron.Mcp.MultipleAlignment",
            Version = "1.0.0"
        };
    })
    .WithStdioServerTransport()
    .WithTools<MultipleAlignmentTools>();

await builder.Build().RunAsync();
