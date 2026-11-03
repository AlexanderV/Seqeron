using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Hosting;
using SuffixTree.Mcp.Core;
using SuffixTree.Mcp.Core.Tools;

var builder = Host.CreateApplicationBuilder(args);

// Per-call size limit of the suffix-tree tools (default 5 000 000 characters; see SuffixTreeMcpLimits).
SuffixTreeMcpLimits.Configure(builder.Configuration);

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
