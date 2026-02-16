/**
* SimplexNoise 2.0.0
* DevDad - Afan Olovcic @ www.art-and-code.com - 2015
* Solessfir - 2026
*/

#include "SimplexNoiseBlueprintFunctionLibrary.h"

#if PLATFORM_COMPILER_CLANG && PLATFORM_CPU_X86_FAMILY
#include <immintrin.h>
#endif

namespace
{
	/**
	* Lookup table for 4D simplex noise. Tells us the order in which to visit the five corners of a 4D simplex shape for any given input position.
	*
	* There are 24 valid entries.
	* The remaining 40 entries are zero-padded and correspond to combinations of comparisons that can never occur in practice.
	* The index into this table is a 6-bit code built from comparing the four input coordinates pairwise (see SimplexNoise4D_Internal).
	*/
	constexpr uint8 SimplexTraversalTable[64][4] =
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
	* Per-thread permutation table used to look up gradient directions.
	*
	* Each thread has its own copy so that different threads can run with different seeds simultaneously without interfering with each other.
	* This table is initialized with Ken Perlin's original 1983 permutation.
	* Call SetSimplexNoiseSeed() to replace it with a shuffled version.
	*/
	thread_local uint8 PermutationTable[256] =
	{
		151, 160, 137,  91,  90,  15, 131,  13, 201,  95,  96,  53, 194, 233,   7, 225,
		140,  36, 103,  30,  69, 142,   8,  99,  37, 240,  21,  10,  23, 190,   6, 148,
		247, 120, 234,  75,   0,  26, 197,  62,  94, 252, 219, 203, 117,  35,  11,  32,
		 57, 177,  33,  88, 237, 149,  56,  87, 174,  20, 125, 136, 171, 168,  68, 175,
		 74, 165,  71, 134, 139,  48,  27, 166,  77, 146, 158, 231,  83, 111, 229, 122,
		 60, 211, 133, 230, 220, 105,  92,  41,  55,  46, 245,  40, 244, 102, 143,  54,
		 65,  25,  63, 161,   1, 216,  80,  73, 209,  76, 132, 187, 208,  89,  18, 169,
		200, 196, 135, 130, 116, 188, 159,  86, 164, 100, 109, 198, 173, 186,   3,  64,
		 52, 217, 226, 250, 124, 123,   5, 202,  38, 147, 118, 126, 255,  82,  85, 212,
		207, 206,  59, 227,  47,  16,  58,  17, 182, 189,  28,  42, 223, 183, 170, 213,
		119, 248, 152,   2,  44, 154, 163,  70, 221, 153, 101, 155, 167,  43, 172,   9,
		129,  22,  39, 253,  19,  98, 108, 110,  79, 113, 224, 232, 178, 185, 112, 104,
		218, 246,  97, 228, 251,  34, 242, 193, 238, 210, 144,  12, 191, 179, 162, 241,
		 81,  51, 145, 235, 249,  14, 239, 107,  49, 192, 214,  31, 181, 199, 106, 157,
		184,  84, 204, 176, 115, 121,  50,  45, 127,   4, 150, 254, 138, 236, 205,  93,
		222, 114,  67,  29,  24,  72, 243, 141, 128, 195,  78,  66, 215,  61, 156, 180
	};

	/**
	* Looks up a value from the permutation table.
	* Casting to uint8 wraps the index into [0, 255], preventing out-of-bounds access when the input value equals 256 (e.g. when a grid corner index is 255 + 1).
	*/
	uint8 Hash(const int32 Value)
	{
		return PermutationTable[static_cast<uint8>(Value)];
	}

	/**
	* Rounds a floating-point number down to the nearest integer.
	* Faster than std::floor() because it avoids a function call and branch overhead.
	* Reference: https://en.wikipedia.org/wiki/Floor_and_ceiling_functions
	*/
	int32 FastFloor(const double Value)
	{
		const int32 Truncated = static_cast<int32>(Value);
		return (Value < static_cast<double>(Truncated)) ? Truncated - 1 : Truncated;
	}

	/**
	* Computes the directional contribution of a 1D lattice corner to the noise value.
	* The gradient direction is derived from the hash: magnitude 1-8, random sign.
	* X is the distance from the corner to the sample point.
	*/
	double ComputeGradient(const int32 Hash, const double X)
	{
		double Gradient = 1.0 + static_cast<double>(Hash & 0x7);
		if (Hash & 8)
		{
			Gradient = -Gradient;
		}
		return Gradient * X;
	}

	/**
	* Computes the directional contribution of a 2D lattice corner to the noise value.
	* Selects one of 12 gradient directions based on the low 3 bits of the hash.
	* X and Y are the distances from the corner to the sample point.
	*/
	double ComputeGradient(const int32 Hash, const double X, const double Y)
	{
		const double U = (Hash & 1) ? X : Y;
		const double V = (Hash & 1) ? Y : X;
		return ((Hash & 2) ? -U : U) + ((Hash & 4) ? -2.0 * V : 2.0 * V);
	}

	/**
	* Computes the directional contribution of a 3D lattice corner to the noise value.
	* Selects one of 12 gradient directions based on the low 4 bits of the hash.
	* X, Y, Z are the distances from the corner to the sample point.
	*/
	double ComputeGradient(const int32 Hash, const double X, const double Y, const double Z)
	{
		const int32 H = Hash & 15;
		const double U = (H < 8) ? X : Y;
		const double V = (H < 4) ? Y : ((H == 12 || H == 14) ? X : Z);
		return ((H & 1) ? -U : U) + ((H & 2) ? -V : V);
	}

	/**
	* Computes the directional contribution of a 4D lattice corner to the noise value.
	* Selects one of 32 gradient directions based on the low 5 bits of the hash.
	* X, Y, Z, W are the distances from the corner to the sample point.
	*/
	double ComputeGradient(const int32 Hash, const double X, const double Y, const double Z, const double W)
	{
		const int32 H = Hash & 31;
		double A = Y;
		double B = Z;
		double C = W;

		switch (H >> 3)
		{
			case 1: A = W; B = X; C = Y; break;
			case 2: A = Z; B = W; C = X; break;
			case 3: A = Y; B = W; C = Z; break;
			default: break;
		}

		return ((H & 4) ? -A : A) + ((H & 2) ? -B : B) + ((H & 1) ? -C : C);
	}

	/**
	* Remaps a noise value from [-1, 1] to [MinRange, MaxRange].
	* Used by all public-facing functions to convert raw noise into a useful range.
	*/
	double RemapToRange(const double Value, const double MinRange, const double MaxRange)
	{
		return (Value * 0.5 + 0.5) * (MaxRange - MinRange) + MinRange;
	}

#if PLATFORM_COMPILER_CLANG && PLATFORM_CPU_X86_FAMILY

	// Evaluates the 2D gradient dot product for four samples simultaneously.
	// Scalar version: u = (h&1) ? x : y; v = (h&1) ? y : x; result = ((h&2) ? -u : u) + ((h&4) ? -2v : 2v)
	// SIMD version: uses _mm256_blendv_pd to select x or y per lane, and sign coefficient vectors to apply the conditional negation.
	__attribute__((target("avx,fma")))
	static __m256d ComputeVecGradient2D(const __m256d VecX, const __m256d VecY, const int32 H0, const int32 H1, const int32 H2, const int32 H3)
	{
		// Build a mask where each lane is all-ones (-0.0 bit pattern) when bit 1 of H is set.
		// _mm256_blendv_pd picks from VecX when the mask MSB is 1, VecY otherwise.
		const __m256d VecMaskBit1 = _mm256_set_pd(
			(H3 & 1) ? -0.0 : 0.0,
			(H2 & 1) ? -0.0 : 0.0,
			(H1 & 1) ? -0.0 : 0.0,
			(H0 & 1) ? -0.0 : 0.0);

		const __m256d VecU = _mm256_blendv_pd(VecY, VecX, VecMaskBit1);
		const __m256d VecV = _mm256_blendv_pd(VecX, VecY, VecMaskBit1);

		// Build sign coefficients: +1 or -1 for U (bit 2), +2 or -2 for V (bit 4).
		const __m256d VecSignU = _mm256_set_pd(
			(H3 & 2) ? -1.0 : 1.0,
			(H2 & 2) ? -1.0 : 1.0,
			(H1 & 2) ? -1.0 : 1.0,
			(H0 & 2) ? -1.0 : 1.0);

		const __m256d VecSignV = _mm256_set_pd(
			(H3 & 4) ? -2.0 : 2.0,
			(H2 & 4) ? -2.0 : 2.0,
			(H1 & 4) ? -2.0 : 2.0,
			(H0 & 4) ? -2.0 : 2.0);

		// Final: SignU * U + SignV * V, computed with a fused multiply-add.
		return _mm256_fmadd_pd(VecSignV, VecV, _mm256_mul_pd(VecSignU, VecU));
	}

	// Computes the smooth falloff weight for four corners simultaneously.
	// T = 0.5 - x^2 - y^2, then T^4.
	// Negative values are zeroed out (no contribution).
	// The negative check is branchless: we AND the result with a comparison mask.
	__attribute__((target("avx,fma")))
	static __m256d ComputeVecFalloff2D(const __m256d VecX, const __m256d VecY, const __m256d VecZero, const __m256d Vec0_5)
	{
		// T = 0.5 - x^2 - y^2, using FMA: -(x*x) + 0.5, then -(y*y) + that result.
		__m256d VecT = _mm256_fnmadd_pd(VecX, VecX, Vec0_5);
		VecT = _mm256_fnmadd_pd(VecY, VecY, VecT);

		// Zero out any lane where T < 0 (the corner is too far away to contribute).
		const __m256d VecPositiveMask = _mm256_cmp_pd(VecT, VecZero, _CMP_GE_OQ);
		VecT = _mm256_and_pd(VecT, VecPositiveMask);

		// T^2, then T^4.
		VecT = _mm256_mul_pd(VecT, VecT);
		return _mm256_mul_pd(VecT, VecT);
	}

#endif // PLATFORM_COMPILER_CLANG && PLATFORM_CPU_X86_FAMILY

} // namespace

double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise1D(const double X, const double Scale, const double MinRange, const double MaxRange)
{
	return RemapToRange(SimplexNoise1D_Internal(X * Scale), MinRange, MaxRange);
}

double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise2D(const double X, const double Y, const double Scale, const double MinRange, const double MaxRange)
{
	return RemapToRange(SimplexNoise2D_Internal(X * Scale, Y * Scale), MinRange, MaxRange);
}

double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise3D(const double X, const double Y, const double Z, const double Scale, const double MinRange, const double MaxRange)
{
	return RemapToRange(SimplexNoise3D_Internal(X * Scale, Y * Scale, Z * Scale), MinRange, MaxRange);
}

double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise4D(const double X, const double Y, const double Z, const double W, const double Scale, const double MinRange, const double MaxRange)
{
	return RemapToRange(SimplexNoise4D_Internal(X * Scale, Y * Scale, Z * Scale, W * Scale), MinRange, MaxRange);
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

	return RemapToRange(Result / Denominator, MinRange, MaxRange);
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

	return RemapToRange(Result / Denominator, MinRange, MaxRange);
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

	return RemapToRange(Result / Denominator, MinRange, MaxRange);
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

	return RemapToRange(Result / Denominator, MinRange, MaxRange);
}

void USimplexNoiseBlueprintFunctionLibrary::SetSimplexNoiseSeed(const int32 NewSeed)
{
	const FRandomStream RandStream(NewSeed);

	for (int32 Index = 0; Index < 256; ++Index)
	{
		PermutationTable[Index] = static_cast<uint8>(Index);
	}

	// Fisher-Yates shuffle: walks backward through the table, swapping each entry with a randomly chosen earlier entry.
	// Produces a uniformly random permutation.
	for (int32 Index = 255; Index > 0; --Index)
	{
		const int32 SwapIndex = RandStream.RandRange(0, Index);
		Swap(PermutationTable[Index], PermutationTable[SwapIndex]);
	}
}

double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise1D_Internal(const double X)
{
	// Find the two integer lattice points on either side of X.
	const int32 I0 = FastFloor(X);
	const int32 I1 = I0 + 1;

	// Distance from X to each lattice point.
	const double X0 = X - static_cast<double>(I0);
	const double X1 = X0 - 1.0;

	// Compute the falloff weight for each corner.
	// The weight is (1 - dist^2)^4, which smoothly fades to zero as we move away from the corner.
	double T0 = 1.0 - X0 * X0;
	T0 *= T0;
	const double N0 = T0 * T0 * ComputeGradient(Hash(I0), X0);

	double T1 = 1.0 - X1 * X1;
	T1 *= T1;
	const double N1 = T1 * T1 * ComputeGradient(Hash(I1), X1);

	// Empirically determined scale factor that maps the output to approximately [-1, 1].
	constexpr double Scale1D = 0.395;
	return Scale1D * (N0 + N1);
}

double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise2D_Internal(const double X, const double Y)
{
	// Skew and unskew factors for 2D simplex noise. Skewing maps the square grid to a simplex (triangle) grid and back.
	// F2 = (sqrt(3) - 1) / 2, G2 = (3 - sqrt(3)) / 6
	constexpr double F2 = 0.366025403784438;
	constexpr double G2 = 0.211324865405187;

	// Skew the input point into the simplex grid to find which cell it falls in.
	const double Skew = (X + Y) * F2;
	const int32 CellI = FastFloor(X + Skew);
	const int32 CellJ = FastFloor(Y + Skew);

	// Unskew the cell origin back to regular (X, Y) space.
	const double Unskew = static_cast<double>(CellI + CellJ) * G2;
	const double X0 = X - (static_cast<double>(CellI) - Unskew);
	const double Y0 = Y - (static_cast<double>(CellJ) - Unskew);

	// The 2D simplex is a triangle.
	// Determine which triangle we are in by comparing X0 and Y0.
	// This decides which of the two possible triangles within the skewed cell contains our sample point.
	int32 I1;
	int32 J1;

	if (X0 > Y0)
	{
		// Lower triangle: traverse corners in order (0,0) -> (1,0) -> (1,1).
		I1 = 1;
		J1 = 0;
	}
	else
	{
		// Upper triangle: traverse corners in order (0,0) -> (0,1) -> (1,1).
		I1 = 0;
		J1 = 1;
	}

	// Compute the positions of the second and third corners relative to the first.
	const double X1 = X0 - static_cast<double>(I1) + G2;
	const double Y1 = Y0 - static_cast<double>(J1) + G2;
	const double X2 = X0 - 1.0 + 2.0 * G2;
	const double Y2 = Y0 - 1.0 + 2.0 * G2;

	// Look up a gradient direction for each corner using the permutation table.
	const int32 GradIndex0 = Hash(CellI + Hash(CellJ));
	const int32 GradIndex1 = Hash(CellI + I1 + Hash(CellJ + J1));
	const int32 GradIndex2 = Hash(CellI + 1 + Hash(CellJ + 1));

	// Compute the noise contribution from each corner.
	// The contribution is the falloff weight (t^4) multiplied by the gradient dot product.
	// If the falloff would be negative, the corner is too far away and contributes nothing.
	double N0;
	double N1;
	double N2;

	if (double T0 = 0.5 - X0 * X0 - Y0 * Y0; T0 < 0.0)
	{
		N0 = 0.0;
	}
	else
	{
		T0 *= T0;
		N0 = T0 * T0 * ComputeGradient(GradIndex0, X0, Y0);
	}

	if (double T1 = 0.5 - X1 * X1 - Y1 * Y1; T1 < 0.0)
	{
		N1 = 0.0;
	}
	else
	{
		T1 *= T1;
		N1 = T1 * T1 * ComputeGradient(GradIndex1, X1, Y1);
	}

	if (double T2 = 0.5 - X2 * X2 - Y2 * Y2; T2 < 0.0)
	{
		N2 = 0.0;
	}
	else
	{
		T2 *= T2;
		N2 = T2 * T2 * ComputeGradient(GradIndex2, X2, Y2);
	}

	// Empirically determined scale factor that maps the output to approximately [-1, 1].
	constexpr double Scale2D = 45.23065;
	return Scale2D * (N0 + N1 + N2);
}

double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise3D_Internal(const double X, const double Y, const double Z)
{
	// Skew and unskew factors for 3D simplex noise. F3 = 1/3, G3 = 1/6
	constexpr double F3 = 0.333333333333333;
	constexpr double G3 = 0.166666666666666;

	// Skew the input point into the simplex grid to find which cell it falls in.
	const double Skew = (X + Y + Z) * F3;
	const int32 CellI = FastFloor(X + Skew);
	const int32 CellJ = FastFloor(Y + Skew);
	const int32 CellK = FastFloor(Z + Skew);

	// Unskew the cell origin back to regular (X, Y, Z) space.
	const double Unskew = static_cast<double>(CellI + CellJ + CellK) * G3;
	const double X0 = X - (static_cast<double>(CellI) - Unskew);
	const double Y0 = Y - (static_cast<double>(CellJ) - Unskew);
	const double Z0 = Z - (static_cast<double>(CellK) - Unskew);

	// The 3D simplex is a tetrahedron (4 corners).
	// Determine the traversal order by ranking the magnitudes of X0, Y0, Z0.
	// This tells us which edges of the unit cube our sample point is closest to.
	int32 I1, J1, K1;
	int32 I2, J2, K2;

	if (X0 >= Y0)
	{
		if (Y0 >= Z0)
		{
			I1 = 1; J1 = 0; K1 = 0;
			I2 = 1; J2 = 1; K2 = 0;
		}
		else if (X0 >= Z0)
		{
			I1 = 1; J1 = 0; K1 = 0;
			I2 = 1; J2 = 0; K2 = 1;
		}
		else
		{
			I1 = 0; J1 = 0; K1 = 1;
			I2 = 1; J2 = 0; K2 = 1;
		}
	}
	else
	{
		if (Y0 < Z0)
		{
			I1 = 0; J1 = 0; K1 = 1;
			I2 = 0; J2 = 1; K2 = 1;
		}
		else if (X0 < Z0)
		{
			I1 = 0; J1 = 1; K1 = 0;
			I2 = 0; J2 = 1; K2 = 1;
		}
		else
		{
			I1 = 0; J1 = 1; K1 = 0;
			I2 = 1; J2 = 1; K2 = 0;
		}
	}

	// Positions of the second, third, and fourth corners relative to the first.
	const double X1 = X0 - static_cast<double>(I1) + G3;
	const double Y1 = Y0 - static_cast<double>(J1) + G3;
	const double Z1 = Z0 - static_cast<double>(K1) + G3;
	const double X2 = X0 - static_cast<double>(I2) + 2.0 * G3;
	const double Y2 = Y0 - static_cast<double>(J2) + 2.0 * G3;
	const double Z2 = Z0 - static_cast<double>(K2) + 2.0 * G3;
	const double X3 = X0 - 1.0 + 3.0 * G3;
	const double Y3 = Y0 - 1.0 + 3.0 * G3;
	const double Z3 = Z0 - 1.0 + 3.0 * G3;

	// Look up gradient directions for all four corners.
	const int32 GradIndex0 = Hash(CellI + Hash(CellJ + Hash(CellK)));
	const int32 GradIndex1 = Hash(CellI + I1 + Hash(CellJ + J1 + Hash(CellK + K1)));
	const int32 GradIndex2 = Hash(CellI + I2 + Hash(CellJ + J2 + Hash(CellK + K2)));
	const int32 GradIndex3 = Hash(CellI + 1 + Hash(CellJ + 1 + Hash(CellK + 1)));

	// Accumulate contributions from each corner. Corners too far away (T < 0) contribute nothing.
	double N0;
	double N1;
	double N2;
	double N3;

	if (double T0 = 0.6 - X0 * X0 - Y0 * Y0 - Z0 * Z0; T0 < 0.0)
	{
		N0 = 0.0;
	}
	else
	{
		T0 *= T0;
		N0 = T0 * T0 * ComputeGradient(GradIndex0, X0, Y0, Z0);
	}

	if (double T1 = 0.6 - X1 * X1 - Y1 * Y1 - Z1 * Z1; T1 < 0.0)
	{
		N1 = 0.0;
	}
	else
	{
		T1 *= T1;
		N1 = T1 * T1 * ComputeGradient(GradIndex1, X1, Y1, Z1);
	}

	if (double T2 = 0.6 - X2 * X2 - Y2 * Y2 - Z2 * Z2; T2 < 0.0)
	{
		N2 = 0.0;
	}
	else
	{
		T2 *= T2;
		N2 = T2 * T2 * ComputeGradient(GradIndex2, X2, Y2, Z2);
	}

	if (double T3 = 0.6 - X3 * X3 - Y3 * Y3 - Z3 * Z3; T3 < 0.0)
	{
		N3 = 0.0;
	}
	else
	{
		T3 *= T3;
		N3 = T3 * T3 * ComputeGradient(GradIndex3, X3, Y3, Z3);
	}

	// Empirically determined scale factor that maps the output to approximately [-1, 1].
	constexpr double Scale3D = 32.0;
	return Scale3D * (N0 + N1 + N2 + N3);
}

double USimplexNoiseBlueprintFunctionLibrary::SimplexNoise4D_Internal(const double X, const double Y, const double Z, const double W)
{
	// Skew and unskew factors for 4D simplex noise. F4 = (sqrt(5) - 1) / 4, G4 = (5 - sqrt(5)) / 20
	constexpr double F4 = 0.309016994374947;
	constexpr double G4 = 0.138196601125010;

	// Skew the input point into the simplex grid to find which cell it falls in.
	const double Skew = (X + Y + Z + W) * F4;
	const int32 CellI = FastFloor(X + Skew);
	const int32 CellJ = FastFloor(Y + Skew);
	const int32 CellK = FastFloor(Z + Skew);
	const int32 CellL = FastFloor(W + Skew);

	// Unskew the cell origin back to regular space.
	const double Unskew = static_cast<double>(CellI + CellJ + CellK + CellL) * G4;
	const double X0 = X - (static_cast<double>(CellI) - Unskew);
	const double Y0 = Y - (static_cast<double>(CellJ) - Unskew);
	const double Z0 = Z - (static_cast<double>(CellK) - Unskew);
	const double W0 = W - (static_cast<double>(CellL) - Unskew);

	// Build a 6-bit index by comparing the four coordinates pairwise. Each bit records which coordinate is larger in a given pair.
	// The resulting index selects a row from SimplexTraversalTable, which tells us the order in which to visit the five simplex corners.
	// Bit weights: X>Y = 32, X>Z = 16, Y>Z = 8, X>W = 4, Y>W = 2, Z>W = 1
	const int32 OrderingCode =
		((X0 > Y0) ? 32 : 0) |
		((X0 > Z0) ? 16 : 0) |
		((Y0 > Z0) ?  8 : 0) |
		((X0 > W0) ?  4 : 0) |
		((Y0 > W0) ?  2 : 0) |
		((Z0 > W0) ?  1 : 0);

	// Ensure the ordering code maps to one of the 24 valid simplex traversal entries.
	checkf(OrderingCode >= 0 && OrderingCode < 64, TEXT("SimplexNoise4D: ordering code out of range (%d)"), OrderingCode);
	checkf(SimplexTraversalTable[OrderingCode][0] != 0 || SimplexTraversalTable[OrderingCode][1] != 0 || SimplexTraversalTable[OrderingCode][2] != 0 || SimplexTraversalTable[OrderingCode][3] != 0, TEXT("SimplexNoise4D: invalid traversal table entry at index %d"), OrderingCode);

	// The traversal table stores a priority rank (0-3) for each coordinate. Rank 3 = this coordinate gets offset to 1 first (it changed most).
	// We convert each rank to a binary offset (0 or 1) for each of the four intermediate simplex corners we visit between the origin and the far corner.
	const int32 I1 = (SimplexTraversalTable[OrderingCode][0] >= 3) ? 1 : 0;
	const int32 J1 = (SimplexTraversalTable[OrderingCode][1] >= 3) ? 1 : 0;
	const int32 K1 = (SimplexTraversalTable[OrderingCode][2] >= 3) ? 1 : 0;
	const int32 L1 = (SimplexTraversalTable[OrderingCode][3] >= 3) ? 1 : 0;

	const int32 I2 = (SimplexTraversalTable[OrderingCode][0] >= 2) ? 1 : 0;
	const int32 J2 = (SimplexTraversalTable[OrderingCode][1] >= 2) ? 1 : 0;
	const int32 K2 = (SimplexTraversalTable[OrderingCode][2] >= 2) ? 1 : 0;
	const int32 L2 = (SimplexTraversalTable[OrderingCode][3] >= 2) ? 1 : 0;

	const int32 I3 = (SimplexTraversalTable[OrderingCode][0] >= 1) ? 1 : 0;
	const int32 J3 = (SimplexTraversalTable[OrderingCode][1] >= 1) ? 1 : 0;
	const int32 K3 = (SimplexTraversalTable[OrderingCode][2] >= 1) ? 1 : 0;
	const int32 L3 = (SimplexTraversalTable[OrderingCode][3] >= 1) ? 1 : 0;

	// Positions of corners 2 through 5 relative to corner 1 (the origin).
	const double X1 = X0 - static_cast<double>(I1) + G4;
	const double Y1 = Y0 - static_cast<double>(J1) + G4;
	const double Z1 = Z0 - static_cast<double>(K1) + G4;
	const double W1 = W0 - static_cast<double>(L1) + G4;

	const double X2 = X0 - static_cast<double>(I2) + 2.0 * G4;
	const double Y2 = Y0 - static_cast<double>(J2) + 2.0 * G4;
	const double Z2 = Z0 - static_cast<double>(K2) + 2.0 * G4;
	const double W2 = W0 - static_cast<double>(L2) + 2.0 * G4;

	const double X3 = X0 - static_cast<double>(I3) + 3.0 * G4;
	const double Y3 = Y0 - static_cast<double>(J3) + 3.0 * G4;
	const double Z3 = Z0 - static_cast<double>(K3) + 3.0 * G4;
	const double W3 = W0 - static_cast<double>(L3) + 3.0 * G4;

	const double X4 = X0 - 1.0 + 4.0 * G4;
	const double Y4 = Y0 - 1.0 + 4.0 * G4;
	const double Z4 = Z0 - 1.0 + 4.0 * G4;
	const double W4 = W0 - 1.0 + 4.0 * G4;

	// Look up gradient directions for all five corners.
	// All lookups now go through Hash(), which wraps the index to [0, 255] via a uint8 cast.
	const int32 GradIndex0 = Hash(CellI + Hash(CellJ + Hash(CellK + Hash(CellL))));
	const int32 GradIndex1 = Hash(CellI + I1 + Hash(CellJ + J1 + Hash(CellK + K1 + Hash(CellL + L1))));
	const int32 GradIndex2 = Hash(CellI + I2 + Hash(CellJ + J2 + Hash(CellK + K2 + Hash(CellL + L2))));
	const int32 GradIndex3 = Hash(CellI + I3 + Hash(CellJ + J3 + Hash(CellK + K3 + Hash(CellL + L3))));
	const int32 GradIndex4 = Hash(CellI + 1 + Hash(CellJ + 1 + Hash(CellK + 1 + Hash(CellL + 1))));

	// Accumulate contributions from each corner. Corners too far away (T < 0) contribute nothing.
	double N0;
	double N1;
	double N2;
	double N3;
	double N4;

	if (double T0 = 0.6 - X0 * X0 - Y0 * Y0 - Z0 * Z0 - W0 * W0; T0 < 0.0)
	{
		N0 = 0.0;
	}
	else
	{
		T0 *= T0;
		N0 = T0 * T0 * ComputeGradient(GradIndex0, X0, Y0, Z0, W0);
	}

	if (double T1 = 0.6 - X1 * X1 - Y1 * Y1 - Z1 * Z1 - W1 * W1; T1 < 0.0)
	{
		N1 = 0.0;
	}
	else
	{
		T1 *= T1;
		N1 = T1 * T1 * ComputeGradient(GradIndex1, X1, Y1, Z1, W1);
	}

	if (double T2 = 0.6 - X2 * X2 - Y2 * Y2 - Z2 * Z2 - W2 * W2; T2 < 0.0)
	{
		N2 = 0.0;
	}
	else
	{
		T2 *= T2;
		N2 = T2 * T2 * ComputeGradient(GradIndex2, X2, Y2, Z2, W2);
	}

	if (double T3 = 0.6 - X3 * X3 - Y3 * Y3 - Z3 * Z3 - W3 * W3; T3 < 0.0)
	{
		N3 = 0.0;
	}
	else
	{
		T3 *= T3;
		N3 = T3 * T3 * ComputeGradient(GradIndex3, X3, Y3, Z3, W3);
	}

	if (double T4 = 0.6 - X4 * X4 - Y4 * Y4 - Z4 * Z4 - W4 * W4; T4 < 0.0)
	{
		N4 = 0.0;
	}
	else
	{
		T4 *= T4;
		N4 = T4 * T4 * ComputeGradient(GradIndex4, X4, Y4, Z4, W4);
	}

	// Empirically determined scale factor that maps the output to approximately [-1, 1].
	constexpr double Scale4D = 27.0;
	return Scale4D * (N0 + N1 + N2 + N3 + N4);
}

#if PLATFORM_COMPILER_CLANG && PLATFORM_CPU_X86_FAMILY

/**
* AVX + FMA SIMD kernel for 2D simplex noise.
*
* Processes four samples at a time using 256-bit AVX registers.
* The __attribute__((target)) tells the compiler to emit AVX + FMA instructions for this function only, without requiring the whole project to be compiled with AVX flags.
*
* Structure of each loop iteration:
*   1. Skew and floor to find the simplex cell (scalar, due to table lookups).
*   2. Unskew and compute corner offsets (SIMD).
*   3. Compute the smooth falloff weight: T = (0.5 - x^2 - y^2)^4 (SIMD).
*   4. Evaluate the gradient dot product for each corner (SIMD).
*   5. Sum and scale the contributions (SIMD).
*   6. Store four results.
* Any leftover samples that do not fill a group of four are handled by a scalar tail loop.
*/
__attribute__((target("avx,fma")))
void USimplexNoiseBlueprintFunctionLibrary::SimplexNoise2D_Batch_SIMD(const double* RESTRICT InX, const double* RESTRICT InY, const int32 Count, double* RESTRICT OutValues)
{
	constexpr double F2 = 0.366025403784438;
	constexpr double G2 = 0.211324865405187;
	constexpr double Scale2D = 45.23065;

	const __m256d VecF2 = _mm256_set1_pd(F2);
	const __m256d VecG2 = _mm256_set1_pd(G2);
	const __m256d Vec2G2 = _mm256_set1_pd(2.0 * G2);
	const __m256d Vec0_5 = _mm256_set1_pd(0.5);
	const __m256d Vec1_0 = _mm256_set1_pd(1.0);
	const __m256d VecScale = _mm256_set1_pd(Scale2D);
	const __m256d VecZero = _mm256_setzero_pd();

	int32 Index = 0;

	for (; Index + 4 <= Count; Index += 4)
	{
		const __m256d VecX = _mm256_loadu_pd(InX + Index);
		const __m256d VecY = _mm256_loadu_pd(InY + Index);

		// Skew the four sample points into the simplex grid.
		const __m256d VecSkew = _mm256_mul_pd(_mm256_add_pd(VecX, VecY), VecF2);
		const __m256d VecXSkewed = _mm256_add_pd(VecX, VecSkew);
		const __m256d VecYSkewed = _mm256_add_pd(VecY, VecSkew);

		// Floor the skewed coordinates to get the cell indices.
		const __m256d VecCellXF = _mm256_floor_pd(VecXSkewed);
		const __m256d VecCellYF = _mm256_floor_pd(VecYSkewed);

		// Convert floored doubles to int32 for permutation table lookups.
		// _mm256_cvttpd_epi32 packs four int32s into a 128-bit register.
		int32 CellI[4];
		int32 CellJ[4];
		_mm_storeu_si128(reinterpret_cast<__m128i*>(CellI), _mm256_cvttpd_epi32(VecCellXF));
		_mm_storeu_si128(reinterpret_cast<__m128i*>(CellJ), _mm256_cvttpd_epi32(VecCellYF));

		// Unskew the cell origin and compute distance from it to the sample point.
		const __m256d VecUnskew = _mm256_mul_pd(_mm256_add_pd(VecCellXF, VecCellYF), VecG2);
		const __m256d VecX0 = _mm256_sub_pd(VecX, _mm256_sub_pd(VecCellXF, VecUnskew));
		const __m256d VecY0 = _mm256_sub_pd(VecY, _mm256_sub_pd(VecCellYF, VecUnskew));

		// Determine which triangle within the cell contains each sample point.
		// We need this as scalars to feed the permutation table hash lookups.
		double X0s[4];
		double Y0s[4];
		_mm256_storeu_pd(X0s, VecX0);
		_mm256_storeu_pd(Y0s, VecY0);

		int32 I1[4];
		int32 J1[4];
		int32 GradIndex0[4];
		int32 GradIndex1[4];
		int32 GradIndex2[4];

		for (int32 Lane = 0; Lane < 4; ++Lane)
		{
			const bool bXGreaterThanY = X0s[Lane] > Y0s[Lane];
			I1[Lane] = bXGreaterThanY ? 1 : 0;
			J1[Lane] = bXGreaterThanY ? 0 : 1;

			GradIndex0[Lane] = Hash(CellI[Lane] + Hash(CellJ[Lane]));
			GradIndex1[Lane] = Hash(CellI[Lane] + I1[Lane] + Hash(CellJ[Lane] + J1[Lane]));
			GradIndex2[Lane] = Hash(CellI[Lane] + 1 + Hash(CellJ[Lane] + 1));
		}

		// Reload the per-lane triangle offsets as SIMD vectors for the floating-point stage.
		const __m256d VecI1 = _mm256_set_pd(
			static_cast<double>(I1[3]), static_cast<double>(I1[2]),
			static_cast<double>(I1[1]), static_cast<double>(I1[0]));
		const __m256d VecJ1 = _mm256_set_pd(
			static_cast<double>(J1[3]), static_cast<double>(J1[2]),
			static_cast<double>(J1[1]), static_cast<double>(J1[0]));

		// Compute corner positions relative to the cell origin (SIMD).
		const __m256d VecX1 = _mm256_add_pd(_mm256_sub_pd(VecX0, VecI1), VecG2);
		const __m256d VecY1 = _mm256_add_pd(_mm256_sub_pd(VecY0, VecJ1), VecG2);
		const __m256d VecX2 = _mm256_add_pd(_mm256_sub_pd(VecX0, Vec1_0), Vec2G2);
		const __m256d VecY2 = _mm256_add_pd(_mm256_sub_pd(VecY0, Vec1_0), Vec2G2);

		// Evaluate noise contributions from all three corners (SIMD).
		const __m256d VecN0 = _mm256_mul_pd(
			ComputeVecFalloff2D(VecX0, VecY0, VecZero, Vec0_5),
			ComputeVecGradient2D(VecX0, VecY0, GradIndex0[0], GradIndex0[1], GradIndex0[2], GradIndex0[3]));

		const __m256d VecN1 = _mm256_mul_pd(
			ComputeVecFalloff2D(VecX1, VecY1, VecZero, Vec0_5),
			ComputeVecGradient2D(VecX1, VecY1, GradIndex1[0], GradIndex1[1], GradIndex1[2], GradIndex1[3]));

		const __m256d VecN2 = _mm256_mul_pd(
			ComputeVecFalloff2D(VecX2, VecY2, VecZero, Vec0_5),
			ComputeVecGradient2D(VecX2, VecY2, GradIndex2[0], GradIndex2[1], GradIndex2[2], GradIndex2[3]));

		// Sum all three contributions and apply the normalization scale.
		const __m256d VecResult = _mm256_mul_pd(VecScale, _mm256_add_pd(VecN0, _mm256_add_pd(VecN1, VecN2)));

		_mm256_storeu_pd(OutValues + Index, VecResult);
	}

	// Handle any remaining samples that did not fill a full group of four.
	for (; Index < Count; ++Index)
	{
		OutValues[Index] = SimplexNoise2D_Internal(InX[Index], InY[Index]);
	}
}

#else

// Non-Clang or non-x86 platform: SIMD function is never called, but must be defined to satisfy the linker. The dispatch in EvaluateBatch2D skips it entirely.
void USimplexNoiseBlueprintFunctionLibrary::SimplexNoise2D_Batch_SIMD(const double* RESTRICT InX, const double* RESTRICT InY, const int32 Count, double* RESTRICT OutValues)
{
	for (int32 Index = 0; Index < Count; ++Index)
	{
		OutValues[Index] = SimplexNoise2D_Internal(InX[Index], InY[Index]);
	}
}

#endif // PLATFORM_COMPILER_CLANG && PLATFORM_CPU_X86_FAMILY

void USimplexNoiseBlueprintFunctionLibrary::EvaluateBatch2D(const double* RESTRICT InX, const double* RESTRICT InY, const int32 Count, double* RESTRICT OutValues)
{
#if PLATFORM_COMPILER_CLANG && PLATFORM_CPU_X86_FAMILY
	// Check at runtime whether the CPU supports AVX and FMA.
	// This result is cached after the first call so there is no repeated overhead.
	static const bool bSIMDSupported = __builtin_cpu_supports("avx") && __builtin_cpu_supports("fma");

	if (bSIMDSupported)
	{
		SimplexNoise2D_Batch_SIMD(InX, InY, Count, OutValues);
		return;
	}
#endif

	for (int32 Index = 0; Index < Count; ++Index)
	{
		OutValues[Index] = SimplexNoise2D_Internal(InX[Index], InY[Index]);
	}
}

void USimplexNoiseBlueprintFunctionLibrary::SimplexNoise1D_Batch(TConstArrayView<double> InX, const double Scale, TArrayView<double> OutValues, const double MinRange, const double MaxRange)
{
	check(InX.Num() == OutValues.Num());

	const double ClampedScale = FMath::Max<double>(Scale, UE_DOUBLE_SMALL_NUMBER);

	for (int32 Index = 0; Index < InX.Num(); ++Index)
	{
		OutValues[Index] = RemapToRange(SimplexNoise1D_Internal(InX[Index] * ClampedScale), MinRange, MaxRange);
	}
}

void USimplexNoiseBlueprintFunctionLibrary::SimplexNoise2D_Batch(TConstArrayView<double> InX, TConstArrayView<double> InY, const double Scale, TArrayView<double> OutValues, const double MinRange, const double MaxRange)
{
	check(InX.Num() == InY.Num() && InX.Num() == OutValues.Num());

	const int32 Count = InX.Num();
	const double ClampedScale = FMath::Max<double>(Scale, UE_DOUBLE_SMALL_NUMBER);

	// Scale the input coordinates before handing them to the compute function.
	TArray<double> ScaledX;
	TArray<double> ScaledY;
	ScaledX.SetNumUninitialized(Count);
	ScaledY.SetNumUninitialized(Count);

	for (int32 Index = 0; Index < Count; ++Index)
	{
		ScaledX[Index] = InX[Index] * ClampedScale;
		ScaledY[Index] = InY[Index] * ClampedScale;
	}

	// Evaluate raw [-1, 1] noise values (dispatches to SIMD where available).
	EvaluateBatch2D(ScaledX.GetData(), ScaledY.GetData(), Count, OutValues.GetData());

	// Remap all results from [-1, 1] to [MinRange, MaxRange].
	for (int32 Index = 0; Index < Count; ++Index)
	{
		OutValues[Index] = RemapToRange(OutValues[Index], MinRange, MaxRange);
	}
}

void USimplexNoiseBlueprintFunctionLibrary::SimplexNoise3D_Batch(TConstArrayView<double> InX, TConstArrayView<double> InY, TConstArrayView<double> InZ, const double Scale, TArrayView<double> OutValues, const double MinRange, const double MaxRange)
{
	check(InX.Num() == InY.Num() && InX.Num() == InZ.Num() && InX.Num() == OutValues.Num());

	const double ClampedScale = FMath::Max<double>(Scale, UE_DOUBLE_SMALL_NUMBER);

	for (int32 Index = 0; Index < InX.Num(); ++Index)
	{
		OutValues[Index] = RemapToRange(SimplexNoise3D_Internal(InX[Index] * ClampedScale, InY[Index] * ClampedScale, InZ[Index] * ClampedScale), MinRange, MaxRange);
	}
}

void USimplexNoiseBlueprintFunctionLibrary::SimplexNoise4D_Batch(TConstArrayView<double> InX, TConstArrayView<double> InY, TConstArrayView<double> InZ, TConstArrayView<double> InW, const double Scale, TArrayView<double> OutValues, const double MinRange, const double MaxRange)
{
	check(InX.Num() == InY.Num() && InX.Num() == InZ.Num() && InX.Num() == InW.Num() && InX.Num() == OutValues.Num());

	const double ClampedScale = FMath::Max<double>(Scale, UE_DOUBLE_SMALL_NUMBER);

	for (int32 Index = 0; Index < InX.Num(); ++Index)
	{
		OutValues[Index] = RemapToRange(SimplexNoise4D_Internal(InX[Index] * ClampedScale, InY[Index] * ClampedScale, InZ[Index] * ClampedScale, InW[Index] * ClampedScale), MinRange, MaxRange);
	}
}

void USimplexNoiseBlueprintFunctionLibrary::SimplexNoise2D_FBM_Batch(TConstArrayView<double> InX, TConstArrayView<double> InY, const double Scale, const int32 Octaves, const double Lacunarity, const double Persistence, TArrayView<double> OutValues, const double MinRange, const double MaxRange)
{
	check(InX.Num() == InY.Num() && InX.Num() == OutValues.Num());

	const int32 Count = InX.Num();
	const double ClampedScale = FMath::Max<double>(Scale, UE_DOUBLE_SMALL_NUMBER);
	const int32 ClampedOctaves = FMath::Max<int32>(Octaves, 1);

	// Accumulates the weighted sum of all octave contributions per sample.
	TArray<double> AccumulatedNoise;
	AccumulatedNoise.SetNumZeroed(Count);

	// Accumulates the sum of all amplitudes so we can normalize at the end.
	TArray<double> Denominators;
	Denominators.SetNumZeroed(Count);

	// Temporary arrays for the scaled input coordinates of each octave.
	TArray<double> OctaveX;
	TArray<double> OctaveY;
	TArray<double> OctaveNoise;
	OctaveX.SetNumUninitialized(Count);
	OctaveY.SetNumUninitialized(Count);
	OctaveNoise.SetNumUninitialized(Count);

	double Frequency = 1.0;
	double Amplitude = 1.0;

	for (int32 OctaveIndex = 0; OctaveIndex < ClampedOctaves; ++OctaveIndex)
	{
		const double FrequencyScale = ClampedScale * Frequency;

		for (int32 Index = 0; Index < Count; ++Index)
		{
			OctaveX[Index] = InX[Index] * FrequencyScale;
			OctaveY[Index] = InY[Index] * FrequencyScale;
		}

		// Evaluate one octave worth of raw [-1, 1] noise (SIMD where available).
		EvaluateBatch2D(OctaveX.GetData(), OctaveY.GetData(), Count, OctaveNoise.GetData());

		for (int32 Index = 0; Index < Count; ++Index)
		{
			AccumulatedNoise[Index] += Amplitude * OctaveNoise[Index];
			Denominators[Index] += Amplitude;
		}

		Frequency *= Lacunarity;
		Amplitude *= Persistence;
	}

	// Normalize each sample by the sum of its amplitudes, then remap to the output range.
	for (int32 Index = 0; Index < Count; ++Index)
	{
		OutValues[Index] = RemapToRange(AccumulatedNoise[Index] / Denominators[Index], MinRange, MaxRange);
	}
}
