using System.Diagnostics;

namespace Sdcb.Arithmetic.Examples;

/// <summary>Warns at startup when timings will not be representative.</summary>
internal static class PerfWarning
{
    public static void Print()
    {
#if DEBUG
        const bool debugBuild = true;
#else
        const bool debugBuild = false;
#endif
        if (!debugBuild && !Debugger.IsAttached)
        {
            return;
        }
        string reason = (debugBuild, Debugger.IsAttached) switch
        {
            (true, true) => "Debug build with a debugger attached",
            (true, false) => "Debug build",
            _ => "debugger attached",
        };
        ConsoleColor old = Console.ForegroundColor;
        Console.ForegroundColor = ConsoleColor.Yellow;
        Console.WriteLine($"WARNING: {reason}; timings will be much slower than normal.");
        Console.WriteLine("         For real performance, run with -c Release and without a debugger.");
        Console.ForegroundColor = old;
        Console.WriteLine();
    }
}