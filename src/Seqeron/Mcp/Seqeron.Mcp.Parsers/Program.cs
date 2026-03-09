using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Hosting;
using ModelContextProtocol.Server;
using Seqeron.Mcp.Parsers.Tools;

var builder = Host.CreateApplicationBuilder(args);

builder.Services
    .AddMcpServer(options =>
    {
        options.ServerInfo = new()
        {
            Name = "Seqeron.Mcp.Parsers",
            Version = "1.0.0"
        };
    })
    .WithStdioServerTransport()
    .WithTools<ParsersTools>();

await builder.Build().RunAsync();
