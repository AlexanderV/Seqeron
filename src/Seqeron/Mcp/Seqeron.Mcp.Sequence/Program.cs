using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Hosting;
using ModelContextProtocol.Server;
using Seqeron.Mcp.Sequence.Tools;

var builder = Host.CreateApplicationBuilder(args);

builder.Services
    .AddMcpServer(options =>
    {
        options.ServerInfo = new()
        {
            Name = "Seqeron.Mcp.Sequence",
            Version = "1.0.0"
        };
    })
    .WithStdioServerTransport()
    .WithTools<SequenceTools>();

await builder.Build().RunAsync();
