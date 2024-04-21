/**
* SimplexNoise 1.3.0
* DevDad - Afan Olovcic @ www.art-and-code.com - 2024
*
* This algorithm was originally designed by Ken Perlin.
* This code has been adapted and extended from the implementation written by Stefan Gustavson (stegu@itn.liu.se) and modified to fit to Unreal Engine.
*
* This is a clean, fast, modern and free Perlin Simplex noise function.
* If we change float to double it could be even faster but there is no double type in Blueprint
* All Public Functions are BlueprintCallable, so they can be used in every blueprint
*
* From DevDad and Dedicated to you and Unreal Community. Use it free for what ever you want.
* I only request that you mention me in the credits for your game in the way that feels most appropriate to you.
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
	static float SimplexNoise1D(const float X, const float Scale = 0.1f, const float MaxRange = 1.f, const float MinRange = -1.f);

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
	static float SimplexNoise2D(const float X, const float Y, const float Scale = 0.1f, const float MaxRange = 1.f, const float MinRange = -1.f);

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
	static float SimplexNoise3D(const float X, const float Y, const float Z, const float Scale = 0.1f, const float MaxRange = 1.f, const float MinRange = -1.f);

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
	static float SimplexNoise4D(const float X, const float Y, const float Z, const float W, float Scale = 0.1f, const float MaxRange = 1.f, const float MinRange = -1.f);
	
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
	static float SimplexNoise1D_FBM(const float X, const float Scale = 0.1f, const int32 Octaves = 4, const float Lacunarity = 2.3f, const float Persistence = 0.44f, const float MaxRange = 1.f, const float MinRange = -1.f);

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
	static float SimplexNoise2D_FBM(const float X, const float Y, const float Scale = 0.1f, const int32 Octaves = 4, const float Lacunarity = 2.3f, const float Persistence = 0.44f, const float MaxRange = 1.f, const float MinRange = -1.f);

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
	static float SimplexNoise3D_FBM(const float X, const float Y, const float Z, const float Scale = 0.1f, const int32 Octaves = 4, const float Lacunarity = 2.3f, const float Persistence = 0.44f, const float MaxRange = 1.f, const float MinRange = -1.f);

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
	static float SimplexNoise4D_FBM(const float X, const float Y, const float Z, const float W, const float Scale = 0.1f, const int32 Octaves = 4, const float Lacunarity = 2.3f, const float Persistence = 0.44f, const float MaxRange = 1.f, const float MinRange = -1.f);

	UFUNCTION(BlueprintCallable, Category = "SimplexNoise|Seed")
	static void SetSimplexNoiseSeed(const int32 NewSeed);
	
private:
	static constexpr float SimplexNoise1D_Internal(const float X);
	
	static constexpr float SimplexNoise2D_Internal(const float X, const float Y);
	
	static constexpr float SimplexNoise3D_Internal(const float X, const float Y, const float Z);
	
	static constexpr float SimplexNoise4D_Internal(const float X, const float Y, const float Z, const float W);
};
