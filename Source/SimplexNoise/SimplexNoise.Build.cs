/**
* SimplexNoise 1.3.0
* DevDad - Afan Olovcic @ www.art-and-code.com - 2015
* Solessfir - 2024
*/

using UnrealBuildTool;

public class SimplexNoise : ModuleRules
{
	public SimplexNoise(ReadOnlyTargetRules Target) : base(Target)
	{
		PCHUsage = PCHUsageMode.UseExplicitOrSharedPCHs;
		
        PublicDependencyModuleNames.AddRange(new[]
        { 
			"Core", 
			"CoreUObject", 
			"Engine"
        });
	}
}
