/**
* SimplexNoise 1.3.0
* DevDad - Afan Olovcic @ www.art-and-code.com - 2015
* Solessfir - 2025
*/

#pragma once

#include "Kismet/BlueprintFunctionLibrary.h"
#include "SimplexNoiseBlueprintFunctionLibrary.generated.h"

UCLASS()
class SIMPLEXNOISE_API USimplexNoiseBlueprintFunctionLibrary : public UBlueprintFunctionLibrary
{
	GENERATED_BODY()

public:
	/**
	* Simplex Noise 1D
	* @param X Input coordinate
	* @param Scale Frequency multiplier (higher = bigger patterns)
	* @param MinRange Minimum Range Value
	* @param MaxRange Maximum Range Value
	* @return Noise value remapped from [-1, 1] to [MinRange, MaxRange]
	*/
	UFUNCTION(BlueprintPure, Meta = (AdvancedDisplay = 2), Category = "SimplexNoise")
	static double SimplexNoise1D(const double X, const double Scale = 0.1f, const double MinRange = -1.f, const double MaxRange = 1.f);

	/**
	* Simplex Noise 2D
	* @param X Input coordinate
	* @param Y Input coordinate
	* @param Scale Frequency multiplier (higher = bigger patterns)
	* @param MinRange Minimum Range Value
	* @param MaxRange Maximum Range Value
	* @return Noise value remapped from [-1, 1] to [MinRange, MaxRange]
	*/
	UFUNCTION(BlueprintPure, Meta = (AdvancedDisplay = 3), Category = "SimplexNoise")
	static double SimplexNoise2D(const double X, const double Y, const double Scale = 0.1f, const double MinRange = -1.f, const double MaxRange = 1.f);

	/**
	* Simplex Noise 3D
	* @param X Input coordinate
	* @param Y Input coordinate
	* @param Z Input coordinate
	* @param Scale Frequency multiplier (higher = bigger patterns)
	* @param MinRange Minimum Range Value
	* @param MaxRange Maximum Range Value
	* @return Noise value remapped from [-1, 1] to [MinRange, MaxRange]
	*/
	UFUNCTION(BlueprintPure, Meta = (AdvancedDisplay = 4), Category = "SimplexNoise")
	static double SimplexNoise3D(const double X, const double Y, const double Z, const double Scale = 0.1f, const double MinRange = -1.f, const double MaxRange = 1.f);

	/**
	* Simplex Noise 4D
	* @param X Input coordinate
	* @param Y Input coordinate
	* @param Z Input coordinate
	* @param W Input coordinate
	* @param Scale Frequency multiplier (higher = bigger patterns)
	* @param MinRange Minimum Range Value
	* @param MaxRange Maximum Range Value
	* @return Noise value remapped from [-1, 1] to [MinRange, MaxRange]
	*/
	UFUNCTION(BlueprintPure, Meta = (AdvancedDisplay = 5), Category = "SimplexNoise")
	static double SimplexNoise4D(const double X, const double Y, const double Z, const double W, double Scale = 0.1f, const double MinRange = -1.f, const double MaxRange = 1.f);
	
	/**
	* Fractional Brownian Motion (FBM) summation of 1D Simplex Noise
	* @param X Input coordinate
	* @param Scale Frequency multiplier (higher = bigger patterns)
	* @param Octaves Number of noise layers to combine (complexity/Detail level)
	* @param Lacunarity Specifies the frequency multiplier between successive octaves
	* @param Persistence Loss of amplitude between successive octaves. Typical usage uses 1/Lacunarity but can be adjusted for different effects
	* @param MinRange Minimum Range Value
	* @param MaxRange Maximum Range Value
	* @return Noise value remapped from [-1, 1] to [MinRange, MaxRange]
	*/
	UFUNCTION(BlueprintPure, Meta = (AdvancedDisplay = 2), DisplayName = "Simplex Noise 1D (FBM)", Category = "Simplex Noise|FBM")
	static double SimplexNoise1D_FBM(const double X, const double Scale = 0.1f, const int32 Octaves = 4, const double Lacunarity = 2.3f, const double Persistence = 0.44f, const double MinRange = -1.f, const double MaxRange = 1.f);

	/**
	* Fractional Brownian Motion (FBM) summation of 2D Simplex Noise
	* @param X Input coordinate
	* @param Y Input coordinate
	* @param Scale Frequency multiplier (higher = bigger patterns)
	* @param Octaves Number of noise layers to combine (complexity/Detail level)
	* @param Lacunarity Specifies the frequency multiplier between successive octaves
	* @param Persistence Loss of amplitude between successive octaves. Typical usage uses 1/Lacunarity but can be adjusted for different effects
	* @param MinRange Minimum Range Value
	* @param MaxRange Maximum Range Value
	* @return Noise value remapped from [-1, 1] to [MinRange, MaxRange]
	*/
	UFUNCTION(BlueprintPure, Meta = (AdvancedDisplay = 3), DisplayName = "Simplex Noise 2D (FBM)", Category = "Simplex Noise|FBM")
	static double SimplexNoise2D_FBM(const double X, const double Y, const double Scale = 0.1f, const int32 Octaves = 4, const double Lacunarity = 2.3f, const double Persistence = 0.44f, const double MinRange = -1.f, const double MaxRange = 1.f);

	/**
	* Fractional Brownian Motion (FBM) summation of 3D Simplex Noise
	* @param X Input coordinate
	* @param Y Input coordinate
	* @param Z Input coordinate
	* @param Scale Frequency multiplier (higher = bigger patterns)
	* @param Octaves Number of noise layers to combine (complexity/Detail level)
	* @param Lacunarity Specifies the frequency multiplier between successive octaves
	* @param Persistence Loss of amplitude between successive octaves. Typical usage uses 1/Lacunarity but can be adjusted for different effects
	* @param MinRange Minimum Range Value
	* @param MaxRange Maximum Range Value
	* @return Noise value remapped from [-1, 1] to [MinRange, MaxRange]
	*/
	UFUNCTION(BlueprintPure, Meta = (AdvancedDisplay = 4), DisplayName = "Simplex Noise 3D (FBM)", Category = "Simplex Noise|FBM")
	static double SimplexNoise3D_FBM(const double X, const double Y, const double Z, const double Scale = 0.1f, const int32 Octaves = 4, const double Lacunarity = 2.3f, const double Persistence = 0.44f, const double MinRange = -1.f, const double MaxRange = 1.f);

	/**
	* Fractional Brownian Motion (FBM) summation of 4D Simplex Noise
	* @param X Input coordinate
	* @param Y Input coordinate
	* @param Z Input coordinate
	* @param W Input coordinate
	* @param Scale Frequency multiplier (higher = bigger patterns)
	* @param Octaves Number of noise layers to combine (complexity/Detail level)
	* @param Lacunarity Specifies the frequency multiplier between successive octaves
	* @param Persistence Loss of amplitude between successive octaves. Typical usage uses 1/Lacunarity but can be adjusted for different effects
	* @param MinRange Minimum Range Value
	* @param MaxRange Maximum Range Value
	* @return Noise value remapped from [-1, 1] to [MinRange, MaxRange]
	*/
	UFUNCTION(BlueprintPure, Meta = (AdvancedDisplay = 5), DisplayName = "Simplex Noise 4D (FBM)", Category = "Simplex Noise|FBM")
	static double SimplexNoise4D_FBM(const double X, const double Y, const double Z, const double W, const double Scale = 0.1f, const int32 Octaves = 4, const double Lacunarity = 2.3f, const double Persistence = 0.44f, const double MinRange = -1.f, const double MaxRange = 1.f);

	/**
	* Reseeds permutation table using Fisher-Yates shuffle
	* Affects all Noise functions
	*/
	UFUNCTION(BlueprintCallable, Category = "SimplexNoise|Seed")
	static void SetSimplexNoiseSeed(const int32 NewSeed);

private: /* Raw noise implementations returning values in [-1,1] range */
	static constexpr double SimplexNoise1D_Internal(const double X);
	
	static constexpr double SimplexNoise2D_Internal(const double X, const double Y);
	
	static constexpr double SimplexNoise3D_Internal(const double X, const double Y, const double Z);
	
	static constexpr double SimplexNoise4D_Internal(const double X, const double Y, const double Z, const double W);
};
