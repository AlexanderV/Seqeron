using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Hosting;
using ModelContextProtocol.Server;
using SuffixTree.Mcp.Core.Tools;

var builder = Host.CreateApplicationBuilder(args);

builder.Services
    .AddMcpServer(options =>
    {
        options.ServerInfo = new()
        {
            Name = "SuffixTree.Mcp.Core",
            Version = "1.0.0"
        };
    })
    .WithStdioServerTransport()
    .WithTools<SuffixTreeCoreTools>()
    .WithTools<SuffixTreeGenomicsTools>();

await builder.Build().RunAsync();
