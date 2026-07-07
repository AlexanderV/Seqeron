using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Hosting;
using Seqeron.Mcp.Analysis.Tools;

var builder = Host.CreateApplicationBuilder(args);

builder.Services
    .AddMcpServer(options =>
    {
        options.ServerInfo = new()
        {
            Name = "Seqeron.Mcp.Analysis",
            Version = "1.0.0"
        };
    })
    .WithStdioServerTransport()
    .WithTools<AnalysisTools>();

await builder.Build().RunAsync();
