using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Hosting;
using ModelContextProtocol.Server;
using Seqeron.Mcp.MolTools.Tools;

var builder = Host.CreateApplicationBuilder(args);

builder.Services
    .AddMcpServer(options =>
    {
        options.ServerInfo = new()
        {
            Name = "Seqeron.Mcp.MolTools",
            Version = "1.0.0"
        };
    })
    .WithStdioServerTransport()
    .WithTools<MolToolsTools>();

await builder.Build().RunAsync();
