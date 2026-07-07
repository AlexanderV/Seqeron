using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Hosting;
using Seqeron.Mcp.Metagenomics.Tools;

var builder = Host.CreateApplicationBuilder(args);

builder.Services
    .AddMcpServer(options =>
    {
        options.ServerInfo = new()
        {
            Name = "Seqeron.Mcp.Metagenomics",
            Version = "1.0.0"
        };
    })
    .WithStdioServerTransport()
    .WithTools<MetagenomicsTools>();

await builder.Build().RunAsync();
