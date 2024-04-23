/**
* SimplexNoise 1.3.0
* DevDad - Afan Olovcic @ www.art-and-code.com - 2015
* Solessfir - 2024
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
	* @param X Coord of the distance to the corner
	* @param Scale Overall Noise Scale
	* @param MaxRange Maximum Range Value
	* @param MinRange Minimum Range Value
	* @return Noise value in the range of MinRange to MaxRange
	*/
	UFUNCTION(BlueprintPure, Meta = (AdvancedDisplay = 2), Category = "SimplexNoise")
	static double SimplexNoise1D(const double X, const double Scale = 0.1f, const double MaxRange = 1.f, const double MinRange = -1.f);

	/**
	* Simplex Noise 2D
	* @param X Coord of the distance to the corner
	* @param Y Coord of the distance to the corner
	* @param Scale Overall Noise Scale
	* @param MaxRange Maximum Range Value
	* @param MinRange Minimum Range Value
	* @return Noise value in the range of MinRange to MaxRange
	*/
	UFUNCTION(BlueprintPure, Meta = (AdvancedDisplay = 3), Category = "SimplexNoise")
	static double SimplexNoise2D(const double X, const double Y, const double Scale = 0.1f, const double MaxRange = 1.f, const double MinRange = -1.f);

	/**
	* Simplex Noise 3D
	* @param X Coord of the distance to the corner
	* @param Y Coord of the distance to the corner
	* @param Z Coord of the distance to the corner
	* @param Scale Overall Noise Scale
	* @param MaxRange Maximum Range Value
	* @param MinRange Minimum Range Value
	* @return Noise value in the range of MinRange to MaxRange
	*/
	UFUNCTION(BlueprintPure, Meta = (AdvancedDisplay = 4), Category = "SimplexNoise")
	static double SimplexNoise3D(const double X, const double Y, const double Z, const double Scale = 0.1f, const double MaxRange = 1.f, const double MinRange = -1.f);

	/**
	* Simplex Noise 4D
	* @param X Coord of the distance to the corner
	* @param Y Coord of the distance to the corner
	* @param Z Coord of the distance to the corner
	* @param W Coord of the distance to the corner
	* @param Scale Overall Noise Scale
	* @param MaxRange Maximum Range Value
	* @param MinRange Minimum Range Value
	* @return Noise value in the range of MinRange to MaxRange
	*/
	UFUNCTION(BlueprintPure, Meta = (AdvancedDisplay = 5), Category = "SimplexNoise")
	static double SimplexNoise4D(const double X, const double Y, const double Z, const double W, double Scale = 0.1f, const double MaxRange = 1.f, const double MinRange = -1.f);
	
	/**
	* Fractional Brownian Motion (FBM) summation of 1D Simplex Noise
	* @param X Coord of the distance to the corner
	* @param Scale Overall Noise Scale
	* @param Octaves Number of fraction of noise to sum
	* @param Lacunarity Specifies the frequency multiplier between successive octaves
	* @param Persistence Loss of amplitude between successive octaves. Usually: 1 / Lacunarity
	* @param MaxRange Maximum Range Value
	* @param MinRange Minimum Range Value
	* @return Noise value in the range of MinRange to MaxRange
	*/
	UFUNCTION(BlueprintPure, Meta = (AdvancedDisplay = 2), DisplayName = "Simplex Noise 1D (FBM)", Category = "SimplexNoise|FBM")
	static double SimplexNoise1D_FBM(const double X, const double Scale = 0.1f, const int32 Octaves = 4, const double Lacunarity = 2.3f, const double Persistence = 0.44f, const double MaxRange = 1.f, const double MinRange = -1.f);

	/**
	* Fractional Brownian Motion (FBM) summation of 2D Simplex Noise
	* @param X Coord of the distance to the corner
	* @param Y Coord of the distance to the corner
	* @param Scale Overall Noise Scale
	* @param Octaves Number of fraction of noise to sum
	* @param Lacunarity Specifies the frequency multiplier between successive octaves
	* @param Persistence Loss of amplitude between successive octaves. Usually: 1 / Lacunarity
	* @param MaxRange Maximum Range Value
	* @param MinRange Minimum Range Value
	* @return Noise value in the range of MinRange to MaxRange
	*/
	UFUNCTION(BlueprintPure, Meta = (AdvancedDisplay = 3), DisplayName = "Simplex Noise 2D (FBM)", Category = "SimplexNoise|FBM")
	static double SimplexNoise2D_FBM(const double X, const double Y, const double Scale = 0.1f, const int32 Octaves = 4, const double Lacunarity = 2.3f, const double Persistence = 0.44f, const double MaxRange = 1.f, const double MinRange = -1.f);

	/**
	* Fractional Brownian Motion (FBM) summation of 3D Simplex Noise
	* @param X Coord of the distance to the corner
	* @param Y Coord of the distance to the corner
	* @param Z Coord of the distance to the corner
	* @param Scale Overall Noise Scale
	* @param Octaves Number of fraction of noise to sum
	* @param Lacunarity Specifies the frequency multiplier between successive octaves
	* @param Persistence Loss of amplitude between successive octaves. Usually: 1 / Lacunarity
	* @param MaxRange Maximum Range Value
	* @param MinRange Minimum Range Value
	* @return Noise value in the range of MinRange to MaxRange
	*/
	UFUNCTION(BlueprintPure, Meta = (AdvancedDisplay = 4), DisplayName = "Simplex Noise 3D (FBM)", Category = "SimplexNoise|FBM")
	static double SimplexNoise3D_FBM(const double X, const double Y, const double Z, const double Scale = 0.1f, const int32 Octaves = 4, const double Lacunarity = 2.3f, const double Persistence = 0.44f, const double MaxRange = 1.f, const double MinRange = -1.f);

	/**
	* Fractional Brownian Motion (FBM) summation of 4D Simplex Noise
	* @param X Coord of the distance to the corner
	* @param Y Coord of the distance to the corner
	* @param Z Coord of the distance to the corner
	* @param W Coord of the distance to the corner
	* @param Scale Overall Noise Scale
	* @param Octaves Number of fraction of noise to sum
	* @param Lacunarity Specifies the frequency multiplier between successive octaves
	* @param Persistence Loss of amplitude between successive octaves. Usually: 1 / Lacunarity
	* @param MaxRange Maximum Range Value
	* @param MinRange Minimum Range Value
	* @return Noise value in the range of MinRange to MaxRange
	*/
	UFUNCTION(BlueprintPure, Meta = (AdvancedDisplay = 5), DisplayName = "Simplex Noise 4D (FBM)", Category = "SimplexNoise|FBM")
	static double SimplexNoise4D_FBM(const double X, const double Y, const double Z, const double W, const double Scale = 0.1f, const int32 Octaves = 4, const double Lacunarity = 2.3f, const double Persistence = 0.44f, const double MaxRange = 1.f, const double MinRange = -1.f);

	UFUNCTION(BlueprintCallable, Category = "SimplexNoise|Seed")
	static void SetSimplexNoiseSeed(const int32 NewSeed);
	
private:
	static constexpr double SimplexNoise1D_Internal(const double X);
	
	static constexpr double SimplexNoise2D_Internal(const double X, const double Y);
	
	static constexpr double SimplexNoise3D_Internal(const double X, const double Y, const double Z);
	
	static constexpr double SimplexNoise4D_Internal(const double X, const double Y, const double Z, const double W);
};
