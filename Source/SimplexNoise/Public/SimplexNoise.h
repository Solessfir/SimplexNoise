/*
 SimplexNoise 1.3.0
 -----
 DevDad @ art-and-code.com - 2024
*/

#pragma once

#include "Modules/ModuleManager.h"

class FSimplexNoiseModule final : public IModuleInterface
{
public:
	virtual void StartupModule() override;
	
	virtual void ShutdownModule() override;
};
