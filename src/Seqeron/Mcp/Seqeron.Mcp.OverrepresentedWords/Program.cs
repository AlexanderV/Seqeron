// SPDX-License-Identifier: AGPL-3.0-only
// This file is part of Seqeron.Mcp.OverrepresentedWords, licensed under the terms stated in LICENSE (this directory).
using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Hosting;
using Seqeron.Mcp.OverrepresentedWords.Tools;

var builder = Host.CreateApplicationBuilder(args);

builder.Services
    .AddMcpServer(options =>
    {
        options.ServerInfo = new()
        {
            Name = "Seqeron.Mcp.OverrepresentedWords",
            Version = "1.0.0"
        };
    })
    .WithStdioServerTransport()
    .WithTools<OverrepresentedWordsTools>();

await builder.Build().RunAsync();
