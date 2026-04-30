using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Hosting;
using ModelContextProtocol.Server;
using Seqeron.Mcp.Phylogenetics.Tools;

var builder = Host.CreateApplicationBuilder(args);

builder.Services
    .AddMcpServer(options =>
    {
        options.ServerInfo = new()
        {
            Name = "Seqeron.Mcp.Phylogenetics",
            Version = "1.0.0"
        };
    })
    .WithStdioServerTransport()
    .WithTools<PhylogeneticsTools>();

await builder.Build().RunAsync();
