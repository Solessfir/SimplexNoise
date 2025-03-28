/**
* SimplexNoise 1.3.0
* DevDad - Afan Olovcic @ www.art-and-code.com - 2015
* Solessfir - 2025
*/

#include "SimplexNoiseBlueprintFunctionLibrary.h"

namespace
{
	/**
	* 4D simplex traversal lookup table. Contains 24 valid entries (non-zero vectors) and 40 padding zeros.
	* Indexed via c which ranges 0-63 (6-bit ordering code).
	*/
	constexpr uint8 Simplex[64][4] =
	{
		{0, 1, 2, 3}, {0, 1, 3, 2}, {0, 0, 0, 0}, {0, 2, 3, 1}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {1, 2, 3, 0},
		{0, 2, 1, 3}, {0, 0, 0, 0}, {0, 3, 1, 2}, {0, 3, 2, 1}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {1, 3, 2, 0},
		{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0},
		{1, 2, 0, 3}, {0, 0, 0, 0}, {1, 3, 0, 2}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {2, 3, 0, 1}, {2, 3, 1, 0},
		{1, 0, 2, 3}, {1, 0, 3, 2}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {2, 0, 3, 1}, {0, 0, 0, 0}, {2, 1, 3, 0},
		{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0},
		{2, 0, 1, 3}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {3, 0, 1, 2}, {3, 0, 2, 1}, {0, 0, 0, 0}, {3, 1, 2, 0},
		{2, 1, 0, 3}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {3, 1, 0, 2}, {0, 0, 0, 0}, {3, 2, 0, 1}, {3, 2, 1, 0}
	};
	
	/**
	* Thread-local permutation table allows concurrent noise generation with different seeds across threads.
	* Initialized via SetSimplexNoiseSeed().
	* Default contains Ken Perlin's original permutation (1983).
	*/
	thread_local uint8 PermutationTable[256] =
	{
		151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 
		140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23, 190, 6, 148, 
		247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 
		57, 177, 33, 88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 
		74, 165, 71, 134, 139, 48, 27, 166, 77, 146, 158, 231, 83, 111, 229, 122, 
		60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244, 102, 143, 54, 
		65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169, 
		200, 196, 135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 
		52, 217, 226, 250, 124, 123, 5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 
		207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42, 223, 183, 170, 213, 
		119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9, 
		129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 
		218, 246, 97, 228, 251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241, 
		81, 51, 145, 235, 249, 14, 239, 107, 49, 192, 214, 31, 181, 199, 106, 157, 
		184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254, 138, 236, 205, 93, 
		222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180
	};
	
	/**
	* Helper function to hash an integer using the above permutation table
	* This function costs around 1ns, and is called N+1 times for a noise of N dimension.
	*
	* Using a real hash function would be better to improve the "repeatability of 256" of the above permutation table,
	* but fast integer Hash functions uses more time and have bad random properties.
	*/
	constexpr uint8 Hash(const int32 Value)
	{
		return PermutationTable[static_cast<uint8>(Value)];
	}
	
	/**
	* Computes the largest integer value not greater than the double one
	* This method is faster than using (int32_t)std::floor(fp).
	*
	* Measured it to be approximately twice as fast:
	* double:  ~18.4ns instead of ~39.6ns on an AMD APU),
	* double: ~20.6ns instead of ~36.6ns on an AMD APU),
	* Reference: http://www.codeproject.com/Tips/700780/Fast-floor-ceiling-functions
	*/
	constexpr int32 FastFloor(const double Value)
	{
		const int32 IntValue = static_cast<int32>(Value);
		return Value < IntValue ? IntValue - 1 : IntValue;
	}
	
	/**
	* Helper functions to compute gradients-dot-residual vectors
	*
	* Note that these generate gradients of more than unit length.
	* To make a close match with the value range of classic Perlin noise,
	* the final noise values need to be rescaled to fit nicely within [-1, 1].
	* The simplex noise functions as such also have different scaling.
	* Note also that these noise functions are the most practical and useful signed version of Perlin noise.
	*
	* X, Y, Z, W are coords of the distance to the corner
	*/
	
	constexpr double ComputeGradient(const int32 Hash, const double X)
	{
		double Gradient = 1.0 + (Hash & 0x7); // Gradient value 1.0, 2.0, ..., 8.0
		if (Hash & 8) Gradient = -Gradient; // Set a random sign for the gradient
		return Gradient * X;
	}
	
	constexpr double ComputeGradient(const int32 Hash, const double X, const double Y)
	{
		const double U = Hash & 1 ? X : Y;
		const double V = Hash & 1 ? Y : X;
		return ((Hash & 2) ? -U : U) + ((Hash & 4) ? -2.0 * V : 2.0 * V); // and compute the dot product with (x,y)
	}

	constexpr double ComputeGradient(const int32 Hash, const double X, const double Y, const double Z)
	{
		// Convert low 4 bits of hash code into 12 gradient directions
		const int32 h = Hash & 15;
		const double u = h < 8 ? X : Y;
		const double v = h < 4 ? Y : (h == 12 || h == 14) ? X : Z;
		return ((h & 1) ? -u : u) + ((h & 2) ? -v : v);
	}
	
	constexpr double ComputeGradient(const int32 Hash, const double X, const double Y, const double Z, const double W)
	{
		// Selects 3 coordinates from XYZW based on hash bits:
		// - h >> 3 determines coordinate set (4 possible combinations)
		// - h bits 0-2 determine sign for each component
		// Results in 32 unique gradient directions

		const int32 h = Hash & 31;
		double a = Y, b = Z, c = W; // Default to YZW

		switch (h >> 3)
		{
			case 1:
				a = W; b = X; c = Y;
				break;
			case 2:
				a = Z; b = W; c = X;
				break;
			case 3:
				a = Y; b = W; c = Z;
				break;
			default:  // case 0
				break; // Keep default YZW values
		}

		return ((h & 4) ? -a : a) + ((h & 2) ? -b : b) + ((h & 1) ? -c : c);
	}

	/**
	* Remaps noise from [-1,1] to [MinRange, MaxRange]
	* Assumes input value is already normalized by noise functions
	*/
	constexpr double RecomputeRange(const double Value, const double MinRange, const double MaxRange)
	{
		return (Value * 0.5 + 0.5) * (MaxRange - MinRange) + MinRange;
	}
}

double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise1D(const double X, const double Scale, const double MinRange, const double MaxRange)
{
	return RecomputeRange(SimplexNoise1D_Internal(X * Scale), MinRange, MaxRange);
}

double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise2D(const double X, const double Y, const double Scale, const double MinRange, const double MaxRange)
{
	return RecomputeRange(SimplexNoise2D_Internal(X * Scale, Y * Scale), MinRange, MaxRange);
}

double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise3D(const double X, const double Y, const double Z, const double Scale, const double MinRange, const double MaxRange)
{
	return RecomputeRange(SimplexNoise3D_Internal(X * Scale, Y * Scale, Z * Scale), MinRange, MaxRange);
}

double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise4D(const double X, const double Y, const double Z, const double W, const double Scale, const double MinRange, const double MaxRange)
{
	return RecomputeRange(SimplexNoise4D_Internal(X * Scale, Y * Scale, Z * Scale, W * Scale), MinRange, MaxRange);
}

double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise1D_FBM(const double X, const double Scale, const int32 Octaves, const double Lacunarity, const double Persistence, const double MinRange, const double MaxRange)
{
	double Result = 0.0;
	double Denominator = 0.0;
	double Frequency = 1.0;
	double Amplitude = 1.0;

	const double ClampedScale = FMath::Max<double>(Scale, UE_DOUBLE_SMALL_NUMBER);
	const int32 ClampedOctaves = FMath::Max<int32>(Octaves, 1);
	
	for (int32 Index = 0; Index < ClampedOctaves; ++Index)
	{
		Result += Amplitude * SimplexNoise1D_Internal(X * ClampedScale * Frequency);
		Denominator += Amplitude;
		Frequency *= Lacunarity;
		Amplitude *= Persistence;
	}
	
	return RecomputeRange(Result / Denominator, MinRange, MaxRange);
}

double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise2D_FBM(const double X, const double Y, const double Scale, const int32 Octaves, const double Lacunarity, const double Persistence, const double MinRange, const double MaxRange)
{
	double Result = 0.0;
	double Denominator = 0.0;
	double Frequency = 1.0;
	double Amplitude = 1.0;

	const double ClampedScale = FMath::Max<double>(Scale, UE_DOUBLE_SMALL_NUMBER);
	const int32 ClampedOctaves = FMath::Max<int32>(Octaves, 1);
	
	for (int32 Index = 0; Index < ClampedOctaves; ++Index)
	{
		Result += Amplitude * SimplexNoise2D_Internal(X * ClampedScale * Frequency, Y * ClampedScale * Frequency);
		Denominator += Amplitude;
		Frequency *= Lacunarity;
		Amplitude *= Persistence;
	}
	
	return RecomputeRange(Result / Denominator, MinRange, MaxRange);
}

double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise3D_FBM(const double X, const double Y, const double Z, const double Scale, const int32 Octaves, const double Lacunarity, const double Persistence, const double MinRange, const double MaxRange)
{
	double Result = 0.0;
	double Denominator = 0.0;
	double Frequency = 1.0;
	double Amplitude = 1.0;

	const double ClampedScale = FMath::Max<double>(Scale, UE_DOUBLE_SMALL_NUMBER);
	const int32 ClampedOctaves = FMath::Max<int32>(Octaves, 1);
	
	for (int32 Index = 0; Index < ClampedOctaves; ++Index)
	{
		Result += Amplitude * SimplexNoise3D_Internal(X * ClampedScale * Frequency, Y * ClampedScale * Frequency, Z * ClampedScale * Frequency);
		Denominator += Amplitude;
		Frequency *= Lacunarity;
		Amplitude *= Persistence;
	}
	
	return RecomputeRange(Result / Denominator, MinRange, MaxRange);
}

double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise4D_FBM(const double X, const double Y, const double Z, const double W, const double Scale, const int32 Octaves, const double Lacunarity, const double Persistence, const double MinRange, const double MaxRange)
{
	double Result = 0.0;
	double Denominator = 0.0;
	double Frequency = 1.0;
	double Amplitude = 1.0;

	const double ClampedScale = FMath::Max<double>(Scale, UE_DOUBLE_SMALL_NUMBER);
	const int32 ClampedOctaves = FMath::Max<int32>(Octaves, 1);
	
	for (int32 Index = 0; Index < ClampedOctaves; ++Index)
	{
		Result += Amplitude * SimplexNoise4D_Internal(X * ClampedScale * Frequency, Y * ClampedScale * Frequency, Z * ClampedScale * Frequency, W * ClampedScale * Frequency);
		Denominator += Amplitude;
		Frequency *= Lacunarity;
		Amplitude *= Persistence;
	}
	
	return RecomputeRange(Result / Denominator, MinRange, MaxRange);
}

void USimplexNoiseBlueprintFunctionLibrary::SetSimplexNoiseSeed(const int32 NewSeed)
{
	const FRandomStream RandStream(NewSeed);

	for (int32 Index = 0; Index < 256; ++Index)
	{
		PermutationTable[Index] = Index;
	}

	for (int32 Index = 255; Index > 0; --Index)
	{
		const int32 SwapIndex = RandStream.RandRange(0, Index);
		Swap(PermutationTable[Index], PermutationTable[SwapIndex]);
	}
}

constexpr double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise1D_Internal(const double X)
{
	// Corners coordinates (nearest integer values):
	const int32 i0 = FastFloor(X);
	const int32 i1 = i0 + 1;
	
	// Distances to corners (between 0 and 1):
	const double x0 = X - i0;
	const double x1 = x0 - 1.0f;

	// Calculate the contribution from the first corner
	double t0 = 1.0f - x0 * x0;
	t0 *= t0;

	// Noise contribution from the two "corners"
	const double n0 = t0 * t0 * ComputeGradient(Hash(i0), x0);

	// Calculate the contribution from the second corner
	double t1 = 1.0f - x1 * x1;
	t1 *= t1;

	// Noise contribution from the two "corners"
	const double n1 = t1 * t1 * ComputeGradient(Hash(i1), x1);

	// 1D scaling factor (empirically determined)
	constexpr double Scale1D = 0.395f;
	// Normalizes 1D noise to [-1,1]

	return Scale1D * (n0 + n1);
}

constexpr double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise2D_Internal(const double X, const double Y)
{
	// 2D Skew/Unskew factors (sqrt(3)-1)/2 and (3-sqrt(3))/6
	constexpr double F2 = 0.366025403784438;
	constexpr double G2 = 0.211324865405187;

	// Noise contributions from the three corners
	double n0, n1, n2;
	
	// Skew the input space to determine which simplex cell we're in
	const double s = (X + Y) * F2; // Hairy factor for 2D
	const double xs = X + s;
	const double ys = Y + s;
	const int32 i = FastFloor(xs);
	const int32 j = FastFloor(ys);

	// Unskew the cell origin back to (x,y) space
	const double t = static_cast<double>(i + j) * G2;
	const double X0 = i - t;
	const double Y0 = j - t;
	const double x0 = X - X0; // The x,y distances from the cell origin
	const double y0 = Y - Y0;
	
	// For the 2D case, the simplex shape is an equilateral triangle. 

	// Offsets for second (middle) corner of simplex in (i,j) coords
	int32 i1, j1;

	// Determine which simplex we are in.
	if (x0 > y0)
	{
		// lower triangle, XY order: (0,0)->(1,0)->(1,1)
		i1 = 1;
		j1 = 0;
	} 
	else
	{
		// upper triangle, YX order: (0,0)->(0,1)->(1,1)
		i1 = 0;
		j1 = 1;
	}

	// A step of (1,0) in (i,j) means a step of (1-c,-c) in (x,y), and
	// a step of (0,1) in (i,j) means a step of (-c,1-c) in (x,y), where
	// c = (3-sqrt(3))/6

	const double x1 = x0 - i1 + G2; // Offsets for middle corner in (x,y) unskewed coords
	const double y1 = y0 - j1 + G2;
	const double x2 = x0 - 1.0 + 2.0 * G2; // Offsets for last corner in (x,y) unskewed coords
	const double y2 = y0 - 1.0 + 2.0 * G2;
	
	// Work out the hashed gradient indices of the three simplex corners
	const int32 gi0 = Hash(i + Hash(j));
	const int32 gi1 = Hash(i + i1 + Hash(j + j1));
	const int32 gi2 = Hash(i + 1 + Hash(j + 1));
	
	// Calculate the contribution from the first corner
	if (double t0 = 0.5f - x0 * x0 - y0 * y0; t0 < 0.0)
	{
		n0 = 0.0;
	}
	else
	{
		t0 *= t0;
		n0 = t0 * t0 * ComputeGradient(gi0, x0, y0);
	}

	// Calculate the contribution from the second corner
	if (double t1 = 0.5f - x1 * x1 - y1 * y1; t1 < 0.0)
	{
		n1 = 0.0;
	}
	else
	{
		t1 *= t1;
		n1 = t1 * t1 * ComputeGradient(gi1, x1, y1);
	}

	// Calculate the contribution from the third corner
	if (double t2 = 0.5f - x2 * x2 - y2 * y2; t2 < 0.0)
	{
		n2 = 0.0;
	}
	else
	{
		t2 *= t2;
		n2 = t2 * t2 * ComputeGradient(gi2, x2, y2);
	}

	// Add contributions from each corner to get the final noise value.
	// The result is scaled to return values in the interval [-1, 1].
	return 45.23065f * (n0 + n1 + n2);
}

constexpr double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise3D_Internal(const double X, const double Y, const double Z)
{
	// 3D Skew/Unskew factors from simplex mathematics
	constexpr double F3 = 0.333333333333333; // (1.0/3.0) - Skews input to 3D grid
	constexpr double G3 = 0.166666666666666; // (1.0/6.0) - Unskews back to space

	double n0, n1, n2, n3; // Noise contributions from four simplex corners

	// Skew input space to simplified 3D grid
	const double Skew = (X + Y + Z) * F3;
	const int32 i = FastFloor(X + Skew);
	const int32 j = FastFloor(Y + Skew);
	const int32 k = FastFloor(Z + Skew);

	// Calculate unskew factor for cell origin
	const double t = (i + j + k) * G3;

	// Cell origin coordinates in unskewed space
	const double X0 = i - t;
	const double Y0 = j - t;
	const double Z0 = k - t;

	// Distance vectors from cell origin to input point
	const double x0 = X - X0;
	const double y0 = Y - Y0;
	const double z0 = Z - Z0;

	// 3D Simplex Cell Identification
	// ------------------------------
	int32 i1, j1, k1; // Offsets for second corner
	int32 i2, j2, k2; // Offsets for third corner

	if (x0 >= y0)
	{
		if (y0 >= z0)
		{	// X >= Y >= Z
			i1 = 1; j1 = 0; k1 = 0;  // X Y Z order
			i2 = 1; j2 = 1; k2 = 0;
		}
		else if (x0 >= z0)
		{	// X >= Z > Y
			i1 = 1; j1 = 0; k1 = 0;  // X Z Y order
			i2 = 1; j2 = 0; k2 = 1;
		}
		else
		{	// Z > X >= Y
			i1 = 0; j1 = 0; k1 = 1;  // Z X Y order
			i2 = 1; j2 = 0; k2 = 1;
		}
	}
	else
	{	// Primary comparison: Y > X
		if (y0 < z0)
		{	// Y > X < Z
			i1 = 0; j1 = 0; k1 = 1;  // Z Y X order
			i2 = 0; j2 = 1; k2 = 1;
		}
		else if (x0 < z0)
		{	// Y >= Z > X
			i1 = 0; j1 = 1; k1 = 0;  // Y Z X order
			i2 = 0; j2 = 1; k2 = 1;
		}
		else
		{	// Y > X >= Z
			i1 = 0; j1 = 1; k1 = 0;  // Y X Z order
			i2 = 1; j2 = 1; k2 = 0;
		}
	}

	// Calculate corner offsets using simplex pattern
	const double x1 = x0 - i1 + G3; // Second corner (1x G3)
	const double y1 = y0 - j1 + G3;
	const double z1 = z0 - k1 + G3;

	const double x2 = x0 - i2 + 2.0 * G3; // Third corner (2x G3)
	const double y2 = y0 - j2 + 2.0 * G3;
	const double z2 = z0 - k2 + 2.0 * G3;

	const double x3 = x0 - 1.0 + 3.0 * G3; // Fourth corner (3x G3)
	const double y3 = y0 - 1.0 + 3.0 * G3;
	const double z3 = z0 - 1.0 + 3.0 * G3;

	// Generate gradient indices through permutation hashing
	const int32 gi0 = Hash(i + Hash(j + Hash(k))); // Origin corner
	const int32 gi1 = Hash(i + i1 + Hash(j + j1 + Hash(k + k1))); // Second
	const int32 gi2 = Hash(i + i2 + Hash(j + j2 + Hash(k + k2))); // Third
	const int32 gi3 = Hash(i + 1 + Hash(j + 1 + Hash(k + 1))); // Last

	// Calculate noise contributions
	if (double t0 = 0.6 - x0 * x0 - y0 * y0 - z0 * z0; t0 < 0.0)
	{
		n0 = 0.0;
	}
	else
	{
		t0 *= t0;
		n0 = t0 * t0 * ComputeGradient(gi0, x0, y0, z0);
	}

	if (double t1 = 0.6 - x1 * x1 - y1 * y1 - z1 * z1; t1 < 0.0)
	{
		n1 = 0.0;
	}
	else
	{
		t1 *= t1;
		n1 = t1 * t1 * ComputeGradient(gi1, x1, y1, z1);
	}

	if (double t2 = 0.6 - x2 * x2 - y2 * y2 - z2 * z2; t2 < 0.0)
	{
		n2 = 0.0;
	}
	else
	{
		t2 *= t2;
		n2 = t2 * t2 * ComputeGradient(gi2, x2, y2, z2);
	}

	if (double t3 = 0.6 - x3 * x3 - y3 * y3 - z3 * z3; t3 < 0.0)
	{
		n3 = 0.0;
	}
	else
	{
		t3 *= t3;
		n3 = t3 * t3 * ComputeGradient(gi3, x3, y3, z3);
	}

	// Combine and scale contributions
	return 32.0 * (n0 + n1 + n2 + n3);
}

constexpr double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise4D_Internal(const double X, const double Y, const double Z, const double W)
{
	// 4D Skew/Unskew factors from simplex mathematics
	constexpr double F4 = 0.309016994374947; // (sqrt(5)-1)/4 - Skews input to grid
	constexpr double G4 = 0.138196601125010; // (5-sqrt(5))/20 - Unskews back to space

	double n0, n1, n2, n3, n4; // Noise contributions from five simplex corners

	// Skew input space to 4D hypergrid
	const double Skew = (X + Y + Z + W) * F4;
	const int32 i = FastFloor(X + Skew);
	const int32 j = FastFloor(Y + Skew);
	const int32 k = FastFloor(Z + Skew);
	const int32 l = FastFloor(W + Skew);

	// Unskewing factor for cell origin
	const double t = (i + j + k + l) * G4;

	// Cell origin coordinates in unskewed space
	const double X0 = i - t;
	const double Y0 = j - t;
	const double Z0 = k - t;
	const double W0 = l - t;

	// Distance vectors from cell origin to input point
	const double x0 = X - X0;
	const double y0 = Y - Y0;
	const double z0 = Z - Z0;
	const double w0 = W - W0;

	// 4D Simplex Cell Identification
	// ------------------------------
	// 1. Determine magnitude ordering of coordinates through pairwise comparisons
	// 2. Create 6-bit signature (0-63) where each bit represents a > comparison
	//    Bit weights: x>y(32), x>z(16), y>z(8), x>w(4), y>w(2), z>w(1)
	// 3. Use signature to lookup simplex vertex permutation pattern

	const int32 c1 = (x0 > y0) ? 32 : 0;  // x vs y
	const int32 c2 = (x0 > z0) ? 16 : 0;  // x vs z
	const int32 c3 = (y0 > z0) ? 8 : 0;   // y vs z
	const int32 c4 = (x0 > w0) ? 4 : 0;   // x vs w
	const int32 c5 = (y0 > w0) ? 2 : 0;   // y vs w
	const int32 c6 = (z0 > w0) ? 1 : 0;   // z vs w
	const int32 c = c1 + c2 + c3 + c4 + c5 + c6;

	// Validate simplex lookup index (should never fail with proper inputs)
	checkf(c >= 0 && c < 64, TEXT("Invalid simplex index c=%d in 4D noise"), c);
	checkf(Simplex[c][0] != 0 || Simplex[c][1] != 0 || Simplex[c][2] != 0 || Simplex[c][3] != 0, TEXT("Corrupted simplex table entry at c=%d"), c);

	// Determine vertex offsets using simplex permutation pattern
	// ----------------------------------------------------------
	// Simplex[c] contains priority levels (0-3) for each coordinate:
	// 3 = largest component, 0 = smallest
	// Threshold checks convert priority levels to offset multipliers

	// First vertex offsets (largest components)
	const int32 i1 = Simplex[c][0] >= 3 ? 1 : 0;
	const int32 j1 = Simplex[c][1] >= 3 ? 1 : 0;
	const int32 k1 = Simplex[c][2] >= 3 ? 1 : 0;
	const int32 l1 = Simplex[c][3] >= 3 ? 1 : 0;

	// Second vertex offsets (second largest components)
	const int32 i2 = Simplex[c][0] >= 2 ? 1 : 0;
	const int32 j2 = Simplex[c][1] >= 2 ? 1 : 0;
	const int32 k2 = Simplex[c][2] >= 2 ? 1 : 0;
	const int32 l2 = Simplex[c][3] >= 2 ? 1 : 0;

	// Third vertex offsets (second smallest components)
	const int32 i3 = Simplex[c][0] >= 1 ? 1 : 0;
	const int32 j3 = Simplex[c][1] >= 1 ? 1 : 0;
	const int32 k3 = Simplex[c][2] >= 1 ? 1 : 0;
	const int32 l3 = Simplex[c][3] >= 1 ? 1 : 0;

	// Calculate corner offsets
	// ------------------------
	// Each corner position is calculated by:
	// 1. Subtracting the vertex offset (0 or 1)
	// 2. Adding multiple of unskew factor G4
	// 3. Progressively increasing G4 multiple for farther corners

	const double x1 = x0 - i1 + G4;    // Second corner (1x G4)
	const double y1 = y0 - j1 + G4;
	const double z1 = z0 - k1 + G4;
	const double w1 = w0 - l1 + G4;

	const double x2 = x0 - i2 + 2 * G4;  // Third corner (2x G4)
	const double y2 = y0 - j2 + 2 * G4;
	const double z2 = z0 - k2 + 2 * G4;
	const double w2 = w0 - l2 + 2 * G4;

	const double x3 = x0 - i3 + 3 * G4;  // Fourth corner (3x G4)
	const double y3 = y0 - j3 + 3 * G4;
	const double z3 = z0 - k3 + 3 * G4;
	const double w3 = w0 - l3 + 3 * G4;

	const double x4 = x0 - 1.0 + 4 * G4; // Fifth corner (4x G4, full offset)
	const double y4 = y0 - 1.0 + 4 * G4;
	const double z4 = z0 - 1.0 + 4 * G4;
	const double w4 = w0 - 1.0 + 4 * G4;

	// Wrap indices to permutation table size (256 entries)
	const int32 ii = i & 0xff;  // Equivalent to i % 256
	const int32 jj = j & 0xff;
	const int32 kk = k & 0xff;
	const int32 ll = l & 0xff;

	// Calculate noise contributions
	// ----------------------------
	// Each corner's contribution is calculated as:
	// 1. Distance falloff: (0.6 - dist²)^4
	// 2. Gradient value from permutation table lookup
	// 3. Dot product between gradient and distance vector

	// Calculate the contribution from the five corners
	if (double t0 = 0.6f - x0 * x0 - y0 * y0 - z0 * z0 - w0 * w0; t0 < 0.0)
	{
		n0 = 0.0;
	}
	else
	{
		t0 *= t0;
		n0 = t0 * t0 * ComputeGradient(PermutationTable[ii + PermutationTable[jj + PermutationTable[kk + PermutationTable[ll]]]], x0, y0, z0, w0);
	}

	if (double t1 = 0.6f - x1 * x1 - y1 * y1 - z1 * z1 - w1 * w1; t1 < 0.0)
	{
		n1 = 0.0;
	}
	else
	{
		t1 *= t1;
		n1 = t1 * t1 * ComputeGradient(PermutationTable[ii + i1 + PermutationTable[jj + j1 + PermutationTable[kk + k1 + PermutationTable[ll + l1]]]], x1, y1, z1, w1);
	}

	if (double t2 = 0.6f - x2 * x2 - y2 * y2 - z2 * z2 - w2 * w2; t2 < 0.0)
	{
		n2 = 0.0;
	}
	else
	{
		t2 *= t2;
		n2 = t2 * t2 * ComputeGradient(PermutationTable[ii + i2 + PermutationTable[jj + j2 + PermutationTable[kk + k2 + PermutationTable[ll + l2]]]], x2, y2, z2, w2);
	}

	if (double t3 = 0.6f - x3 * x3 - y3 * y3 - z3 * z3 - w3 * w3; t3 < 0.0)
	{
		n3 = 0.0;
	}
	else
	{
		t3 *= t3;
		n3 = t3 * t3 * ComputeGradient(PermutationTable[ii + i3 + PermutationTable[jj + j3 + PermutationTable[kk + k3 + PermutationTable[ll + l3]]]], x3, y3, z3, w3);
	}

	if (double t4 = 0.6f - x4 * x4 - y4 * y4 - z4 * z4 - w4 * w4; t4 < 0.0)
	{
		n4 = 0.0;
	}
	else
	{
		t4 *= t4;
		n4 = t4 * t4 * ComputeGradient(PermutationTable[ii + 1 + PermutationTable[jj + 1 + PermutationTable[kk + 1 + PermutationTable[ll + 1]]]], x4, y4, z4, w4);
	}

	// Final scaling to normalize output to [-1,1]
	// Empirical factor determined through noise range analysis
	return 27.0 * (n0 + n1 + n2 + n3 + n4);
}
