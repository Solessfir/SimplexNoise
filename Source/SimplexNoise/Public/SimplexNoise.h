/**
* SimplexNoise 2.0.0
* DevDad - Afan Olovcic @ www.art-and-code.com - 2015
* Solessfir - 2026
*/

#pragma once

#include "Modules/ModuleManager.h"

class FSimplexNoiseModule final : public IModuleInterface
{
public:
	virtual void StartupModule() override;

	virtual void ShutdownModule() override;
};
