using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Hosting;
using ModelContextProtocol.Server;
using Seqeron.Mcp.Annotation.Tools;

var builder = Host.CreateApplicationBuilder(args);

builder.Services
    .AddMcpServer(options =>
    {
        options.ServerInfo = new()
        {
            Name = "Seqeron.Mcp.Annotation",
            Version = "1.0.0"
        };
    })
    .WithStdioServerTransport()
    .WithTools<AnnotationTools>();

await builder.Build().RunAsync();
