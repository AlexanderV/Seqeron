using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Hosting;
using Seqeron.Mcp.Alignment.Tools;

var builder = Host.CreateApplicationBuilder(args);

builder.Services
    .AddMcpServer(options =>
    {
        options.ServerInfo = new()
        {
            Name = "Seqeron.Mcp.Alignment",
            Version = "1.0.0"
        };
    })
    .WithStdioServerTransport()
    .WithTools<AlignmentTools>();

await builder.Build().RunAsync();
