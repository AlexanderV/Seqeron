using System.Runtime.CompilerServices;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// The library default is <see cref="LimitationMode.Strict"/> (limitation branches throw). The
/// canonical / discipline fixtures exercise those best-effort branches directly, so this assembly
/// runs under <see cref="LimitationMode.Permissive"/>. The default Strict behaviour is covered
/// explicitly by <c>LimitationPolicy_Strict_Tests</c> via scoped <c>UseStrict()</c> overrides.
/// </summary>
internal static class LimitationPolicyTestBootstrap
{
    [ModuleInitializer]
    internal static void Init() => LimitationPolicy.DefaultMode = LimitationMode.Permissive;
}
