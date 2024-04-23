/**
* SimplexNoise 1.3.0
* DevDad - Afan Olovcic @ www.art-and-code.com - 2015
* Solessfir - 2024
*/

#include "SimplexNoiseBlueprintFunctionLibrary.h"

namespace
{
	/**
	* Used only in SimplexNoise4D
	* A lookup table to traverse the simplex around a given point in 4D
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
	* This is just a random jumble of all numbers 0-255
	*
	* This produce a repeatable pattern of 256, but Ken Perlin stated
	* that it is not a problem for graphic texture as the noise features disappear
	* at a distance far enough to be able to see a repeatable pattern of 256.
	*
	* This needs to be exactly the same for all instances on all platforms,
	* so it's easiest to just keep it as static explicit data.
	* This also removes the need for any initialisation of this class.
	*
	* Note that making this an uint32[] instead of an uint8[] might make the
	* code run faster on platforms with a high penalty for unaligned single
	* byte addressing. Intel x86 is generally single-byte-friendly, but
	* some other CPUs are faster with 4-aligned reads.
	* However, a char[] is smaller, which avoids cache trashing, and that
	* is probably the most important aspect on most architectures.
	* This array is accessed a *lot* by the noise functions.
	* A vector-valued noise over 3D accesses it 96 times, and a
	* double-valued 4D noise 64 times. We want this to fit in the cache!
	*/
	
	uint8 PermutationTable[256] = 
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
		return ((Hash & 2) ? -U : U) + ((Hash & 4) ? -2.0 * V : 2.0 * V); // and compute the dot product with (x,y).
	}

	constexpr double ComputeGradient(const int32 Hash, const double X, const double Y, const double Z)
	{
		const uint8 ConvertedHash = Hash & 31;
		double U;
		double V;
		
		if (ConvertedHash < 11)
		{
			U = X;
			V = Y;
		} 
		else if (ConvertedHash < 22)
		{
			U = X;
			V = Z;
		} 
		else
		{
			U = Y;
			V = Z;
		}
		
		return ((ConvertedHash & 32) ? -U : U) + ((ConvertedHash & 64) ? -V : V);
	}
	
	constexpr double ComputeGradient(const int32 Hash, const double X, const double Y, const double Z, const double W)
	{
		const int32 ConvertedHash = Hash & 31;		// Convert low 5 bits of hash code into 32 simple
		const double U = ConvertedHash < 24 ? X : Y; // gradient directions, and compute dot product.
		const double V = ConvertedHash < 16 ? Y : Z;
		const double T = ConvertedHash < 8 ? Z : W;
		return ((ConvertedHash & 1) ? -U : U) + ((ConvertedHash & 2) ? -V : V) + ((ConvertedHash & 4) ? -T : T);
	}
	
	constexpr double RecomputeRange(const double Value, const double MaxRange, const double MinRange)
	{
		return (Value + 1.0) * 0.5f * (FMath::Max(MaxRange, MinRange + 1.0) - MinRange) + MinRange;
	}
}

double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise1D(const double X, const double Scale, const double MaxRange, const double MinRange)
{
	return RecomputeRange(SimplexNoise1D_Internal(X * Scale), MaxRange, MinRange);
}

double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise2D(const double X, const double Y, const double Scale, const double MaxRange, const double MinRange)
{
	return RecomputeRange(SimplexNoise2D_Internal(X * Scale, Y * Scale), MaxRange, MinRange);
}

double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise3D(const double X, const double Y, const double Z, const double Scale, const double MaxRange, const double MinRange)
{
	return RecomputeRange(SimplexNoise3D_Internal(X * Scale, Y * Scale, Z * Scale), MaxRange, MinRange);
}

double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise4D(const double X, const double Y, const double Z, const double W, const double Scale, const double MaxRange, const double MinRange)
{
	return RecomputeRange(SimplexNoise4D_Internal(X * Scale, Y * Scale, Z * Scale, W * Scale), MaxRange, MinRange);
}

double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise1D_FBM(const double X, const double Scale, const int32 Octaves, const double Lacunarity, const double Persistence, const double MaxRange, const double MinRange)
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
	
	return RecomputeRange(Result / Denominator, MaxRange, MinRange);
}

double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise2D_FBM(const double X, const double Y, const double Scale, const int32 Octaves, const double Lacunarity, const double Persistence, const double MaxRange, const double MinRange)
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
	
	return RecomputeRange(Result / Denominator, MaxRange, MinRange);
}

double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise3D_FBM(const double X, const double Y, const double Z, const double Scale, const int32 Octaves, const double Lacunarity, const double Persistence, const double MaxRange, const double MinRange)
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
	
	return RecomputeRange(Result / Denominator, MaxRange, MinRange);
}

double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise4D_FBM(const double X, const double Y, const double Z, const double W, const double Scale, const int32 Octaves, const double Lacunarity, const double Persistence, const double MaxRange, const double MinRange)
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
	
	return RecomputeRange(Result / Denominator, MaxRange, MinRange);
}

void USimplexNoiseBlueprintFunctionLibrary::SetSimplexNoiseSeed(const int32 NewSeed)
{
	const FRandomStream RandomStream = FRandomStream(NewSeed);
	for (auto& TableValue : PermutationTable)
	{
		TableValue = RandomStream.RandRange(0, 255);
	}
}

constexpr double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise1D_Internal(const double X)
{
	// Noise contributions from the two "corners"
	double n0, n1;
	
	// Corners coordinates (nearest integer values):
	int32_t i0 = FastFloor(X);
	int32_t i1 = i0 + 1;
	
	// Distances to corners (between 0 and 1):
	double x0 = X - i0;
	double x1 = x0 - 1.0f;

	// Calculate the contribution from the first corner
	double t0 = 1.0f - x0 * x0;
	t0 *= t0;
	n0 = t0 * t0 * ComputeGradient(Hash(i0), x0);

	// Calculate the contribution from the second corner
	double t1 = 1.0f - x1 * x1;
	t1 *= t1;
	n1 = t1 * t1 * ComputeGradient(Hash(i1), x1);

	// The maximum value of this noise is 8 * (3 / 4) ^ 4 = 2.53125
	// A factor of 0.395 scales to fit exactly within [-1,1]
	return 0.395f * (n0 + n1);
}

constexpr double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise2D_Internal(const double X, const double Y)
{
	// Skewing/Unskewing factors for 2D
	constexpr double F2 = 0.366025403784438; // (FMath::Sqrt(3.0) - 1.0) / 2.0
	constexpr double G2 = 0.211324865405187; // (3.0 - FMath::Sqrt(3.0)) / 6.0

	// Noise contributions from the three corners
	double n0, n1, n2;
	
	// Skew the input space to determine which simplex cell we're in
	double s = (X + Y) * F2; // Hairy factor for 2D
	double xs = X + s;
	double ys = Y + s;
	int32 i = FastFloor(xs);
	int32 j = FastFloor(ys);

	// Unskew the cell origin back to (x,y) space
	double t = static_cast<double>(i + j) * G2;
	double X0 = i - t;
	double Y0 = j - t;
	double x0 = X - X0; // The x,y distances from the cell origin
	double y0 = Y - Y0;
	
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
    double t0 = 0.5f - x0 * x0 - y0 * y0;
	
	if (t0 < 0.0)
	{
		n0 = 0.0;
	}
	else
	{
		t0 *= t0;
		n0 = t0 * t0 * ComputeGradient(gi0, x0, y0);
	}

	// Calculate the contribution from the second corner
	double t1 = 0.5f - x1 * x1 - y1 * y1;
	
	if (t1 < 0.0)
	{
		n1 = 0.0;
	}
	else
	{
		t1 *= t1;
		n1 = t1 * t1 * ComputeGradient(gi1, x1, y1);
	}

	// Calculate the contribution from the third corner
	double t2 = 0.5f - x2 * x2 - y2 * y2;
	
	if (t2 < 0.0)
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
	// Simple skewing factors for the 3D case
	constexpr double F3 = 0.333333333333333; // 1.0 / 3.0
	constexpr double G3 = 0.166666666666666; // 1.0 / 6.0

	double n0, n1, n2, n3; // Noise contributions from the four corners

	// Skew the input space to determine which simplex cell we're in
	const double Skew = (X + Y + Z) * F3;
	const int32 i = FastFloor(X + Skew);
	const int32 j = FastFloor(Y + Skew);
	const int32 k = FastFloor(Z + Skew);

	// Factor for 3D unskewing
	const double t = (i + j + k) * G3;

	// Unskew the cell origin back to (x,y,z) space
	const double X0 = i - t;
	const double Y0 = j - t;
	const double Z0 = k - t;

	// The x,y,z distances from the cell origin
	const double x0 = X - X0;
	const double y0 = Y - Y0;
	const double z0 = Z - Z0;

	// For the 3D case, the simplex shape is a slightly irregular tetrahedron.
	// Determine which simplex we are in.
	int32 i1, j1, k1; // Offsets for second corner of simplex in (i,j,k) coords
	int32 i2, j2, k2; // Offsets for third corner of simplex in (i,j,k) coords

	/* This code would benefit from a backport from the GLSL version! */
	if (x0 >= y0)
	{
		if (y0 >= z0)
		{
			// X Y Z order
			i1 = 1;
			j1 = 0;
			k1 = 0;
			i2 = 1;
			j2 = 1;
			k2 = 0;
		} 
		else if (x0 >= z0)
		{
			// X Z Y order
			i1 = 1;
			j1 = 0;
			k1 = 0;
			i2 = 1;
			j2 = 0;
			k2 = 1;
		}
		else
		{
			// Z X Y order
			i1 = 0;
			j1 = 0;
			k1 = 1;
			i2 = 1;
			j2 = 0;
			k2 = 1;
		}
	}
	else // x0 < y0
	{
		if (y0 < z0)
		{
			// Z Y X order
			i1 = 0;
			j1 = 0;
			k1 = 1;
			i2 = 0;
			j2 = 1;
			k2 = 1;
		}
		else if (x0 < z0)
		{
			// Y Z X order
			i1 = 0;
			j1 = 1;
			k1 = 0;
			i2 = 0;
			j2 = 1;
			k2 = 1;
		}
		else
		{
			// Y X Z order
			i1 = 0;
			j1 = 1;
			k1 = 0;
			i2 = 1;
			j2 = 1;
			k2 = 0;
		}
	}

	// A step of (1,0,0) in (i,j,k) means a step of (1-c,-c,-c) in (x,y,z),
	// a step of (0,1,0) in (i,j,k) means a step of (-c,1-c,-c) in (x,y,z), and
	// a step of (0,0,1) in (i,j,k) means a step of (-c,-c,1-c) in (x,y,z), where
	// c = 1/6.

	const double x1 = x0 - i1 + G3; // Offsets for second corner in (x,y,z) coords
	const double y1 = y0 - j1 + G3;
	const double z1 = z0 - k1 + G3;
	const double x2 = x0 - i2 + 2.0 * G3; // Offsets for third corner in (x,y,z) coords
	const double y2 = y0 - j2 + 2.0 * G3;
	const double z2 = z0 - k2 + 2.0 * G3;
	const double x3 = x0 - 1.0 + 3.0 * G3; // Offsets for last corner in (x,y,z) coords
	const double y3 = y0 - 1.0 + 3.0 * G3;
	const double z3 = z0 - 1.0 + 3.0 * G3;

	// Work out the hashed gradient indices of the four simplex corners
	const int32 gi0 = Hash(i + Hash(j + Hash(k)));
	const int32 gi1 = Hash(i + i1 + Hash(j + j1 + Hash(k + k1)));
	const int32 gi2 = Hash(i + i2 + Hash(j + j2 + Hash(k + k2)));
	const int32 gi3 = Hash(i + 1 + Hash(j + 1 + Hash(k + 1)));

	// Calculate the contribution from the four corners
	double t0 = 0.6f - x0 * x0 - y0 * y0 - z0 * z0;
	
	if (t0 < 0.0)
	{
		n0 = 0.0;
	}
	else
	{
		t0 *= t0;
		n0 = t0 * t0 * ComputeGradient(gi0, x0, y0, z0);
	}

	double t1 = 0.6f - x1 * x1 - y1 * y1 - z1 * z1;
	
	if (t1 < 0.0)
	{
		n1 = 0.0;
	}
	else
	{
		t1 *= t1;
		n1 = t1 * t1 * ComputeGradient(gi1, x1, y1, z1);
	}

	double t2 = 0.6f - x2 * x2 - y2 * y2 - z2 * z2;
	
	if (t2 < 0.0)
	{
		n2 = 0.0;
	}
	else
	{
		t2 *= t2;
		n2 = t2 * t2 * ComputeGradient(gi2, x2, y2, z2);
	}

	double t3 = 0.6f - x3 * x3 - y3 * y3 - z3 * z3;
	
	if (t3 < 0.0)
	{
		n3 = 0.0;
	}
	else
	{
		t3 *= t3;
		n3 = t3 * t3 * ComputeGradient(gi3, x3, y3, z3);
	}

	// Add contributions from each corner to get the final noise value.
	// The result is scaled to stay just inside [-1,1]
	return 32.0f * (n0 + n1 + n2 + n3);
}

constexpr double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise4D_Internal(const double X, const double Y, const double Z, const double W)
{
	constexpr double F4 = 0.309016994374947; // (FMath::Sqrt(5.0) - 1.0) / 4.0
	constexpr double G4 = 0.138196601125010; // (5.0 - FMath::Sqrt(5.0)) / 20.0
	
	double n0, n1, n2, n3, n4; // Noise contributions from the five corners

	// Skew the (x,y,z,w) space to determine which cell of 24 simplices we're in
	double Skew = (X + Y + Z + W) * F4; // Factor for 4D skewing
	const int32 i = FastFloor(X + Skew);
	const int32 j = FastFloor(Y + Skew);
	const int32 k = FastFloor(Z + Skew);
	const int32 l = FastFloor(W + Skew);

	// Factor for 4D unskewing
	const double t = (i + j + k + l) * G4;
	
	// Unskew the cell origin back to (x,y,z,w) space
	const double X0 = i - t;
	const double Y0 = j - t;
	const double Z0 = k - t;
	const double W0 = l - t;
	
	// The x,y,z,w distances from the cell origin
	const double x0 = X - X0;
	const double y0 = Y - Y0;
	const double z0 = Z - Z0;
	const double w0 = W - W0;

	// For the 4D case, the simplex is a 4D shape I won't even try to describe.
	// To find out which of the 24 possible simplices we're in, we need to
	// determine the magnitude ordering of x0, y0, z0 and w0.
	// The method below is a good way of finding the ordering of x,y,z,w and
	// then find the correct traversal order for the simplex were in.
	// First, six pair-wise comparisons are performed between each possible pair
	// of the four coordinates, and the results are used to add up binary bits
	// for an integer index.
	int32 c1 = (x0 > y0) ? 32 : 0;
	int32 c2 = (x0 > z0) ? 16 : 0;
	int32 c3 = (y0 > z0) ? 8 : 0;
	int32 c4 = (x0 > w0) ? 4 : 0;
	int32 c5 = (y0 > w0) ? 2 : 0;
	int32 c6 = (z0 > w0) ? 1 : 0;
	int32 c = c1 + c2 + c3 + c4 + c5 + c6;

	int32 i1, j1, k1, l1; // The integer offsets for the second simplex corner
	int32 i2, j2, k2, l2; // The integer offsets for the third simplex corner
	int32 i3, j3, k3, l3; // The integer offsets for the fourth simplex corner

	// simplex[c] is a 4-vector with the numbers 0, 1, 2 and 3 in some order.
	// Many values of c will never occur, since e.g. x>y>z>w makes x<z, y<w and x<w impossible.
	// Only the 24 indices which have non-zero entries make any sense.
	// We use a thresholding to set the coordinates in turn from the largest magnitude.
	// The number 3 in the "simplex" array is at the position of the largest coordinate.
	i1 = Simplex[c][0] >= 3 ? 1 : 0;
	j1 = Simplex[c][1] >= 3 ? 1 : 0;
	k1 = Simplex[c][2] >= 3 ? 1 : 0;
	l1 = Simplex[c][3] >= 3 ? 1 : 0;
	// The number 2 in the "simplex" array is at the second - largest coordinate.
	i2 = Simplex[c][0] >= 2 ? 1 : 0;
	j2 = Simplex[c][1] >= 2 ? 1 : 0;
	k2 = Simplex[c][2] >= 2 ? 1 : 0;
	l2 = Simplex[c][3] >= 2 ? 1 : 0;
	// The number 1 in the "simplex" array is at the second - smallest coordinate.
	i3 = Simplex[c][0] >= 1 ? 1 : 0;
	j3 = Simplex[c][1] >= 1 ? 1 : 0;
	k3 = Simplex[c][2] >= 1 ? 1 : 0;
	l3 = Simplex[c][3] >= 1 ? 1 : 0;
	// The fifth corner has all coordinate offsets = 1, so no need to look that up.

	double x1 = x0 - i1 + G4; // Offsets for second corner in (x,y,z,w) coords
	double y1 = y0 - j1 + G4;
	double z1 = z0 - k1 + G4;
	double w1 = w0 - l1 + G4;
	double x2 = x0 - i2 + 2.0 * G4; // Offsets for third corner in (x,y,z,w) coords
	double y2 = y0 - j2 + 2.0 * G4;
	double z2 = z0 - k2 + 2.0 * G4;
	double w2 = w0 - l2 + 2.0 * G4;
	double x3 = x0 - i3 + 3.0 * G4; // Offsets for fourth corner in (x,y,z,w) coords
	double y3 = y0 - j3 + 3.0 * G4;
	double z3 = z0 - k3 + 3.0 * G4;
	double w3 = w0 - l3 + 3.0 * G4;
	double x4 = x0 - 1.0 + 4.0 * G4; // Offsets for last corner in (x,y,z,w) coords
	double y4 = y0 - 1.0 + 4.0 * G4;
	double z4 = z0 - 1.0 + 4.0 * G4;
	double w4 = w0 - 1.0 + 4.0 * G4;
	
	// Wrap the integer indices at 256, to avoid indexing PermutationTable[] out of bounds
	int32 ii = i & 0xff;
	int32 jj = j & 0xff;
	int32 kk = k & 0xff;
	int32 ll = l & 0xff;

	// Calculate the contribution from the five corners
	double t0 = 0.6f - x0 * x0 - y0 * y0 - z0 * z0 - w0 * w0;
	
	if (t0 < 0.0)
	{
		n0 = 0.0;
	}
	else
	{
		t0 *= t0;
		n0 = t0 * t0 * ComputeGradient(PermutationTable[ii + PermutationTable[jj + PermutationTable[kk + PermutationTable[ll]]]], x0, y0, z0, w0);
	}

	double t1 = 0.6f - x1 * x1 - y1 * y1 - z1 * z1 - w1 * w1;
	
	if (t1 < 0.0)
	{
		n1 = 0.0;
	}
	else
	{
		t1 *= t1;
		n1 = t1 * t1 * ComputeGradient(PermutationTable[ii + i1 + PermutationTable[jj + j1 + PermutationTable[kk + k1 + PermutationTable[ll + l1]]]], x1, y1, z1, w1);
	}

	double t2 = 0.6f - x2 * x2 - y2 * y2 - z2 * z2 - w2 * w2;
	
	if (t2 < 0.0)
	{
		n2 = 0.0;
	}
	else
	{
		t2 *= t2;
		n2 = t2 * t2 * ComputeGradient(PermutationTable[ii + i2 + PermutationTable[jj + j2 + PermutationTable[kk + k2 + PermutationTable[ll + l2]]]], x2, y2, z2, w2);
	}

	double t3 = 0.6f - x3 * x3 - y3 * y3 - z3 * z3 - w3 * w3;
	
	if (t3 < 0.0)
	{
		n3 = 0.0;
	}
	else
	{
		t3 *= t3;
		n3 = t3 * t3 * ComputeGradient(PermutationTable[ii + i3 + PermutationTable[jj + j3 + PermutationTable[kk + k3 + PermutationTable[ll + l3]]]], x3, y3, z3, w3);
	}

	double t4 = 0.6f - x4 * x4 - y4 * y4 - z4 * z4 - w4 * w4;
	
	if (t4 < 0.0)
	{
		n4 = 0.0;
	}
	else
	{
		t4 *= t4;
		n4 = t4 * t4 * ComputeGradient(PermutationTable[ii + 1 + PermutationTable[jj + 1 + PermutationTable[kk + 1 + PermutationTable[ll + 1]]]], x4, y4, z4, w4);
	}

	// Sum up and scale the result to cover the range [-1,1]
	return 27.0 * (n0 + n1 + n2 + n3 + n4);
}
