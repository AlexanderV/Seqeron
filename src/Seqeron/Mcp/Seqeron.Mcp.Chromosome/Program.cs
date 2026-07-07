using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Hosting;
using Seqeron.Mcp.Chromosome.Tools;

var builder = Host.CreateApplicationBuilder(args);

builder.Services
    .AddMcpServer(options =>
    {
        options.ServerInfo = new()
        {
            Name = "Seqeron.Mcp.Chromosome",
            Version = "1.0.0"
        };
    })
    .WithStdioServerTransport()
    .WithTools<ChromosomeTools>();

await builder.Build().RunAsync();
