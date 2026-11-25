// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.OligoThermodynamics, licensed under the terms stated in LICENSE (this directory).
using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Hosting;
using Seqeron.Mcp.OligoThermodynamics.Tools;

var builder = Host.CreateApplicationBuilder(args);

builder.Services
    .AddMcpServer(options =>
    {
        options.ServerInfo = new()
        {
            Name = "Seqeron.Mcp.OligoThermodynamics",
            Version = "1.0.0"
        };
    })
    .WithStdioServerTransport()
    .WithTools<OligoThermodynamicsTools>();

await builder.Build().RunAsync();
