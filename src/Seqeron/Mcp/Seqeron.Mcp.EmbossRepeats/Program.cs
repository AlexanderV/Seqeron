// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.EmbossRepeats, licensed under the terms stated in LICENSE (this directory).
using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Hosting;
using Seqeron.Mcp.EmbossRepeats.Tools;

var builder = Host.CreateApplicationBuilder(args);

builder.Services
    .AddMcpServer(options =>
    {
        options.ServerInfo = new()
        {
            Name = "Seqeron.Mcp.EmbossRepeats",
            Version = "1.0.0"
        };
    })
    .WithStdioServerTransport()
    .WithTools<EmbossRepeatsTools>();

await builder.Build().RunAsync();
