using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Hosting;
using ModelContextProtocol.Server;
using Seqeron.Mcp.Population.Tools;

var builder = Host.CreateApplicationBuilder(args);

builder.Services
    .AddMcpServer(options =>
    {
        options.ServerInfo = new()
        {
            Name = "Seqeron.Mcp.Population",
            Version = "1.0.0"
        };
    })
    .WithStdioServerTransport()
    .WithTools<PopulationTools>();

await builder.Build().RunAsync();
