/**
* SimplexNoise 2.0.0
* DevDad - Afan Olovcic @ www.art-and-code.com - 2015
* Solessfir - 2026
*/

#pragma once

#include "Kismet/BlueprintFunctionLibrary.h"
#include "SimplexNoiseBlueprintFunctionLibrary.generated.h"

/**
* Provides Simplex Noise sampling for Blueprints and C++.
*
* Supports 1D through 4D noise and Fractional Brownian Motion (FBM) variants.
* All Blueprint functions remap the raw [-1, 1] noise output to a custom range.
*
* For large workloads (terrain generation, texture baking), use the C++ batch functions.
* They process many points at once and use SIMD instructions on supported hardware to run significantly faster than looping over the scalar functions individually.
*
* Each thread keeps its own private noise seed, so seeding on one thread does not change results on other threads.
* Call SetSimplexNoiseSeed() separately on every thread that needs a non-default seed.
*/
UCLASS()
class SIMPLEXNOISE_API USimplexNoiseBlueprintFunctionLibrary : public UBlueprintFunctionLibrary
{
	GENERATED_BODY()

public:

	/**
	* Evaluates 1D Simplex Noise at position X.
	* @param X          Input position along a line.
	* @param Scale      Controls the zoom level. Higher values produce larger, smoother features.
	* @param MinRange   Minimum value of the output range.
	* @param MaxRange   Maximum value of the output range.
	* @return           Noise value remapped to [MinRange, MaxRange].
	*/
	UFUNCTION(BlueprintPure, Meta = (AdvancedDisplay = 2), Category = "SimplexNoise")
	static double SimplexNoise1D(const double X, const double Scale = 0.1, const double MinRange = -1.0, const double MaxRange = 1.0);

	/**
	* Evaluates 2D Simplex Noise at position (X, Y).
	* Useful for heightmaps, flat terrain, and surface textures.
	* @param X          Input X position.
	* @param Y          Input Y position.
	* @param Scale      Controls the zoom level. Higher values produce larger, smoother features.
	* @param MinRange   Minimum value of the output range.
	* @param MaxRange   Maximum value of the output range.
	* @return           Noise value remapped to [MinRange, MaxRange].
	*/
	UFUNCTION(BlueprintPure, Meta = (AdvancedDisplay = 3), Category = "SimplexNoise")
	static double SimplexNoise2D(const double X, const double Y, const double Scale = 0.1, const double MinRange = -1.0, const double MaxRange = 1.0);

	/**
	* Evaluates 3D Simplex Noise at position (X, Y, Z).
	* Useful for volumetric effects such as clouds, fog, or 3D terrain carving.
	* @param X          Input X position.
	* @param Y          Input Y position.
	* @param Z          Input Z position.
	* @param Scale      Controls the zoom level. Higher values produce larger, smoother features.
	* @param MinRange   Minimum value of the output range.
	* @param MaxRange   Maximum value of the output range.
	* @return           Noise value remapped to [MinRange, MaxRange].
	*/
	UFUNCTION(BlueprintPure, Meta = (AdvancedDisplay = 4), Category = "SimplexNoise")
	static double SimplexNoise3D(const double X, const double Y, const double Z, const double Scale = 0.1, const double MinRange = -1.0, const double MaxRange = 1.0);

	/**
	* Evaluates 4D Simplex Noise at position (X, Y, Z, W).
	* The W axis is commonly used as a time offset to animate 3D noise smoothly without visible repetition or sliding.
	* @param X          Input X position.
	* @param Y          Input Y position.
	* @param Z          Input Z position.
	* @param W          Input W position (commonly used as time for animation).
	* @param Scale      Controls the zoom level. Higher values produce larger, smoother features.
	* @param MinRange   Minimum value of the output range.
	* @param MaxRange   Maximum value of the output range.
	* @return           Noise value remapped to [MinRange, MaxRange].
	*/
	UFUNCTION(BlueprintPure, Meta = (AdvancedDisplay = 5), Category = "SimplexNoise")
	static double SimplexNoise4D(const double X, const double Y, const double Z, const double W, const double Scale = 0.1, const double MinRange = -1.0, const double MaxRange = 1.0);

	/**
	* Evaluates 1D Simplex Noise using Fractional Brownian Motion (FBM).
	*
	* FBM layers multiple passes of noise on top of each other.
	* Each layer (octave) is finer and quieter than the previous one, building up natural-looking complexity similar to real terrain or clouds.
	*
	* @param X            Input position along a line.
	* @param Scale        Controls the zoom level of the base layer.
	* @param Octaves      Number of noise layers to stack. More octaves = more fine detail.
	* @param Lacunarity   How much finer (higher frequency) each next layer becomes. Typically ~2.0.
	* @param Persistence  How much quieter (lower amplitude) each next layer becomes. Typically ~0.5.
	* @param MinRange     Minimum value of the output range.
	* @param MaxRange     Maximum value of the output range.
	* @return             Noise value remapped to [MinRange, MaxRange].
	*/
	UFUNCTION(BlueprintPure, Meta = (AdvancedDisplay = 2), DisplayName = "Simplex Noise 1D (FBM)", Category = "Simplex Noise|FBM")
	static double SimplexNoise1D_FBM(const double X, const double Scale = 0.1, const int32 Octaves = 4, const double Lacunarity = 2.3, const double Persistence = 0.44, const double MinRange = -1.0, const double MaxRange = 1.0);

	/**
	* Evaluates 2D Simplex Noise using Fractional Brownian Motion (FBM).
	*
	* FBM layers multiple passes of noise on top of each other.
	* Each layer (octave) is finer and quieter than the previous one, building up natural-looking complexity similar to real terrain or clouds.
	*
	* @param X            Input X position.
	* @param Y            Input Y position.
	* @param Scale        Controls the zoom level of the base layer.
	* @param Octaves      Number of noise layers to stack. More octaves = more fine detail.
	* @param Lacunarity   How much finer (higher frequency) each next layer becomes. Typically ~2.0.
	* @param Persistence  How much quieter (lower amplitude) each next layer becomes. Typically ~0.5.
	* @param MinRange     Minimum value of the output range.
	* @param MaxRange     Maximum value of the output range.
	* @return             Noise value remapped to [MinRange, MaxRange].
	*/
	UFUNCTION(BlueprintPure, Meta = (AdvancedDisplay = 3), DisplayName = "Simplex Noise 2D (FBM)", Category = "Simplex Noise|FBM")
	static double SimplexNoise2D_FBM(const double X, const double Y, const double Scale = 0.1, const int32 Octaves = 4, const double Lacunarity = 2.3, const double Persistence = 0.44, const double MinRange = -1.0, const double MaxRange = 1.0);

	/**
	* Evaluates 3D Simplex Noise using Fractional Brownian Motion (FBM).
	*
	* FBM layers multiple passes of noise on top of each other.
	* Each layer (octave) is finer and quieter than the previous one, building up natural-looking complexity similar to real terrain or clouds.
	*
	* @param X            Input X position.
	* @param Y            Input Y position.
	* @param Z            Input Z position.
	* @param Scale        Controls the zoom level of the base layer.
	* @param Octaves      Number of noise layers to stack. More octaves = more fine detail.
	* @param Lacunarity   How much finer (higher frequency) each next layer becomes. Typically ~2.0.
	* @param Persistence  How much quieter (lower amplitude) each next layer becomes. Typically ~0.5.
	* @param MinRange     Minimum value of the output range.
	* @param MaxRange     Maximum value of the output range.
	* @return             Noise value remapped to [MinRange, MaxRange].
	*/
	UFUNCTION(BlueprintPure, Meta = (AdvancedDisplay = 4), DisplayName = "Simplex Noise 3D (FBM)", Category = "Simplex Noise|FBM")
	static double SimplexNoise3D_FBM(const double X, const double Y, const double Z, const double Scale = 0.1, const int32 Octaves = 4, const double Lacunarity = 2.3, const double Persistence = 0.44, const double MinRange = -1.0, const double MaxRange = 1.0);

	/**
	* Evaluates 4D Simplex Noise using Fractional Brownian Motion (FBM).
	*
	* FBM layers multiple passes of noise on top of each other.
	* Each layer (octave) is finer and quieter than the previous one, building up natural-looking complexity similar to real terrain or clouds.
	*
	* @param X            Input X position.
	* @param Y            Input Y position.
	* @param Z            Input Z position.
	* @param W            Input W position (commonly used as time for animation).
	* @param Scale        Controls the zoom level of the base layer.
	* @param Octaves      Number of noise layers to stack. More octaves = more fine detail.
	* @param Lacunarity   How much finer (higher frequency) each next layer becomes. Typically ~2.0.
	* @param Persistence  How much quieter (lower amplitude) each next layer becomes. Typically ~0.5.
	* @param MinRange     Minimum value of the output range.
	* @param MaxRange     Maximum value of the output range.
	* @return             Noise value remapped to [MinRange, MaxRange].
	*/
	UFUNCTION(BlueprintPure, Meta = (AdvancedDisplay = 5), DisplayName = "Simplex Noise 4D (FBM)", Category = "Simplex Noise|FBM")
	static double SimplexNoise4D_FBM(const double X, const double Y, const double Z, const double W, const double Scale = 0.1, const int32 Octaves = 4, const double Lacunarity = 2.3, const double Persistence = 0.44, const double MinRange = -1.0, const double MaxRange = 1.0);

	/**
	* Re-seeds the noise generator for the calling thread using a Fisher-Yates shuffle.
	* Affects all noise functions called on this thread.
	* Must be called separately on each thread that needs a unique seed.
	*/
	UFUNCTION(BlueprintCallable, Category = "SimplexNoise|Seed")
	static void SetSimplexNoiseSeed(const int32 NewSeed);

	/**
	* C++ batch functions evaluate many noise points at once.
	* They are significantly faster than calling the scalar functions in a loop because they amortize overhead and exploit CPU caches better.
	*
	* The 2D functions additionally use AVX + FMA SIMD on Clang x86-64 hardware, processing four samples per clock cycle in the floating-point stage.
	*
	* All input array views must have the same length as OutValues.
	*/

	static void SimplexNoise1D_Batch(TConstArrayView<double> InX, double Scale, TArrayView<double> OutValues, double MinRange = 0.0, double MaxRange = 1.0);

	static void SimplexNoise2D_Batch(TConstArrayView<double> InX, TConstArrayView<double> InY, double Scale, TArrayView<double> OutValues, double MinRange = 0.0, double MaxRange = 1.0);

	static void SimplexNoise3D_Batch(TConstArrayView<double> InX, TConstArrayView<double> InY, TConstArrayView<double> InZ, double Scale, TArrayView<double> OutValues, double MinRange = 0.0, double MaxRange = 1.0);

	static void SimplexNoise4D_Batch(TConstArrayView<double> InX, TConstArrayView<double> InY, TConstArrayView<double> InZ, TConstArrayView<double> InW, double Scale, TArrayView<double> OutValues, double MinRange = 0.0, double MaxRange = 1.0);

	static void SimplexNoise2D_FBM_Batch(TConstArrayView<double> InX, TConstArrayView<double> InY, double Scale, int32 Octaves, double Lacunarity, double Persistence, TArrayView<double> OutValues, double MinRange = 0.0, double MaxRange = 1.0);

private:
	/**
	* Raw noise evaluators.
	* Output is in approximately [-1, 1] before any range remapping.
	*/
	static double SimplexNoise1D_Internal(double X);
	static double SimplexNoise2D_Internal(double X, double Y);
	static double SimplexNoise3D_Internal(double X, double Y, double Z);
	static double SimplexNoise4D_Internal(double X, double Y, double Z, double W);

	/**
	* Shared compute function used by all 2D batch entry points.
	* Writes raw [-1, 1] results for Count samples into OutValues.
	* Dispatches to the SIMD kernel on supported hardware, otherwise runs a scalar loop.
	*/
	static void EvaluateBatch2D(const double* RESTRICT InX, const double* RESTRICT InY, int32 Count, double* RESTRICT OutValues);

	/**
	* AVX + FMA SIMD kernel. Processes four 2D samples per iteration.
	* Only compiled when targeting Clang on x86-64.
	* Always guarded by a runtime CPU feature check before being called.
	*/
	static void SimplexNoise2D_Batch_SIMD(const double* RESTRICT InX, const double* RESTRICT InY, int32 Count, double* RESTRICT OutValues);
};
