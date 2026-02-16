using UnrealBuildTool;

public class SimplexNoise : ModuleRules
{
	public SimplexNoise(ReadOnlyTargetRules Target) : base(Target)
	{
		PCHUsage = PCHUsageMode.UseExplicitOrSharedPCHs;

        PublicDependencyModuleNames.AddRange([
	        "Core",
			"CoreUObject",
			"Engine"
        ]);
	}
}
