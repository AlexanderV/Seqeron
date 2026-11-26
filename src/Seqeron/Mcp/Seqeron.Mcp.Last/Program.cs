// SPDX-License-Identifier: GPL-3.0-or-later
// This file is part of Seqeron.Mcp.Last, licensed under the terms stated in LICENSE (this directory).
using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Hosting;
using Seqeron.Mcp.Last.Tools;

var builder = Host.CreateApplicationBuilder(args);

builder.Services
    .AddMcpServer(options =>
    {
        options.ServerInfo = new()
        {
            Name = "Seqeron.Mcp.Last",
            Version = "1.0.0"
        };
    })
    .WithStdioServerTransport()
    .WithTools<LastTools>();

await builder.Build().RunAsync();
