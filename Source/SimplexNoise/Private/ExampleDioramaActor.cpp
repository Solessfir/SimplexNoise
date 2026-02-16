/**
* BiomeDioramaActor.cpp
*
* HOW THE GENERATION PIPELINE WORKS
* --------------------------------------------------------
* Step 1 - Seed
*   SetSimplexNoiseSeed() is called once. This scrambles the internal
*   permutation table so every noise sample produces a unique result for
*   this seed. Change NoiseSeed in the Details panel to get a new landscape.
*
* Step 2 - Heightmap  (2D FBM, one sample per column)
*   For each (X, Y) position we call SimplexNoise2D_FBM. The output is
*   remapped directly to [TerrainHeightMin, TerrainHeightMax] in blocks.
*   This single value tells us how tall that column of terrain is.
*
*   Why FBM? Plain simplex noise produces gentle, single-frequency bumps.
*   FBM stacks several octaves (layers) of noise at increasing frequencies
*   and decreasing amplitudes, which adds small hills on top of large hills,
*   producing terrain that looks like it belongs in nature.
*
* Step 3 - Biome blend  (2D simplex, one sample per column, different scale)
*   A second 2D noise lookup at a much lower frequency produces a smooth map
*   that varies between -1 and +1 across the diorama. Positive values shift
*   the surface block toward one look (e.g. sand), negative toward another
*   (e.g. grass). Combined with altitude rules this gives beach and snow zones.
*
* Step 4 - Cave carving  (3D FBM, one sample per voxel)
*   For every solid block we check a 3D noise value. If it rises above
*   CaveThreshold the block is removed. Because this is 3D noise, the carved
*   spaces form connected tunnels and caverns rather than flat holes.
*   Cave carving is suppressed above CaveCeilingZ to protect the surface.
*
* Step 5 - Layer assignment
*   Solid, non-cave blocks are assigned a Role based on:
*     - Is this the topmost block in the column? -> Surface / Sand / Snow
*     - Is this within SubsurfaceDepth of the top?  -> Subsurface
*     - Everything else?                            -> Deep
*     - Is this block exposed to a cave above it?   -> CaveWall (optional)
*
* Step 6 - Instancing
*   Each block role maps to an entry in the BlockTypes array which has a
*   mesh and optional material. All blocks of the same type share one
*   UInstancedStaticMeshComponent, keeping draw calls minimal.
*/

#include "ExampleDioramaActor.h"
#include "SimplexNoiseBlueprintFunctionLibrary.h"

AExampleDioramaActor::AExampleDioramaActor()
{
	PrimaryActorTick.bCanEverTick = false;

	USceneComponent* SceneRoot = CreateDefaultSubobject<USceneComponent>(TEXT("SceneRoot"));
	SetRootComponent(SceneRoot);
}

void AExampleDioramaActor::BeginPlay()
{
	Super::BeginPlay();

	// Optionally auto-generate on play
	// GenerateDiorama();
}

void AExampleDioramaActor::ClearDiorama()
{
	// Destroy all ISMC components that were created during the last generation pass
	for (UInstancedStaticMeshComponent* Component : BlockComponents)
	{
		if (IsValid(Component))
		{
			Component->DestroyComponent();
		}
	}

	BlockComponents.Empty();
}

void AExampleDioramaActor::GenerateDiorama()
{
	ClearDiorama();

	if (BlockTypes.IsEmpty())
	{
		UE_LOG(LogTemp, Warning, TEXT("BiomeDioramaActor: BlockTypes array is empty. Add at least one block type before generating."));
		return;
	}

	// Clamp world dimensions to be safe even if the designer entered odd values
	const int32 SizeX = FMath::Max(DioramaSize.X, 4);
	const int32 SizeY = FMath::Max(DioramaSize.Y, 4);
	const int32 SizeZ = FMath::Max(DioramaSize.Z, 4);

	/**
	* STEP 1: SEED
	*
	* Seed the noise generator. This affects all subsequent noise calls on this thread.
	* The seed is per-thread, so editor and game threads are independent.
	*/
	USimplexNoiseBlueprintFunctionLibrary::SetSimplexNoiseSeed(NoiseSeed);

	// Pre-calculate valid height range
	// We clamp the max to one below SizeZ so there is always at least one air block above the tallest column, which keeps the diorama open at the top
	const int32 ClampedHeightMin = FMath::Clamp(TerrainHeightMin, 1, SizeZ - 1);
	const int32 ClampedHeightMax = FMath::Clamp(TerrainHeightMax, ClampedHeightMin + 1, SizeZ - 1);

	// Create one ISMC per block type entry.
	// Blocks with no mesh assigned are skipped during instance placement (role still resolves, mesh just does not appear)
	// Keeping empty slots allows designers to comment out a block type without reshuffling the whole array
	BlockComponents.SetNumZeroed(BlockTypes.Num());

	for (int32 TypeIndex = 0; TypeIndex < BlockTypes.Num(); ++TypeIndex)
	{
		const FBiomeBlockType& BlockType = BlockTypes[TypeIndex];

		if (!IsValid(BlockType.Mesh))
		{
			continue;
		}

		UInstancedStaticMeshComponent* ISMC = NewObject<UInstancedStaticMeshComponent>(this);
		ISMC->SetStaticMesh(BlockType.Mesh);

		if (IsValid(BlockType.Material))
		{
			ISMC->SetMaterial(0, BlockType.Material);
		}

		// Register and attach the new component so it is part of the actor
		ISMC->RegisterComponent();
		ISMC->AttachToComponent(GetRootComponent(), FAttachmentTransformRules::SnapToTargetNotIncludingScale);

		BlockComponents[TypeIndex] = ISMC;
	}

	/**
	* GENERATION LOOP
	*
	* We iterate over every column (X, Y) in the diorama. For each column we:
	*   1. Sample 2D FBM noise to get the terrain height for that column.
	*   2. Sample the biome blend noise to decide surface block appearance.
	*   3. Walk vertically from Z = 0 up to the column height.
	*   4. For each solid voxel, sample 3D cave noise to decide if it is hollow.
	*   5. Assign a block role and add an instance to the matching ISMC.
	*/
	for (int32 GridX = 0; GridX < SizeX; ++GridX)
	{
		for (int32 GridY = 0; GridY < SizeY; ++GridY)
		{
			/**
			* STEP 2: HEIGHTMAP
			*
			* SimplexNoise2D_FBM returns a value in [MinRange, MaxRange].
			* We pass our height limits directly as the range so the output
			* is already in block units - no manual remapping needed.
			*
			* TerrainScale controls the spatial frequency. Smaller values
			* zoom the noise out, making features broader. The X and Y
			* coordinates are passed in block units; the scale handles the rest.
			*/
			const double HeightNoise = USimplexNoiseBlueprintFunctionLibrary::SimplexNoise2D_FBM(
				static_cast<double>(GridX),
				static_cast<double>(GridY),
				TerrainScale,
				TerrainOctaves,
				TerrainLacunarity,
				TerrainPersistence,
				static_cast<double>(ClampedHeightMin),
				static_cast<double>(ClampedHeightMax));

			const int32 ColumnHeight = FMath::RoundToInt(HeightNoise);

			/**
			* STEP 3: BIOME BLEND
			*
			* A second, lower-frequency 2D noise sample produces a value in [-1, 1].
			* We use it together with altitude rules to decide what the surface
			* block looks like. This is sampled at a very different scale from the
			* height map so biome patches do not simply mirror the terrain.
			*
			* Note: We offset the input by (1000, 1000) to decorelate this noise
			* from the height noise even though they share the same seed.
			*/
			const double BiomeNoise = USimplexNoiseBlueprintFunctionLibrary::SimplexNoise2D(
				static_cast<double>(GridX) + 1000.0,
				static_cast<double>(GridY) + 1000.0,
				BiomeBlendScale,
				-1.0,
				1.0);

			// Determine surface appearance rules for this column
			//   Beach: column height is close to sea level
			//   Snow:  column height is at or above SnowAltitude, or the biome noise pushes a high column over the threshold
			const bool bIsBeachColumn = (ColumnHeight <= SeaLevel + BeachWidth);
			const bool bIsSnowColumn = (ColumnHeight >= SnowAltitude) || (ColumnHeight >= SnowAltitude - 2 && BiomeNoise > 0.4);

			// Walk vertically from the bottom of the diorama up to (and including) the column's terrain height
			// Everything above ColumnHeight is air
			for (int32 GridZ = 0; GridZ < ColumnHeight; ++GridZ)
			{
				/**
				* STEP 4: CAVE CARVING
				*
				* SimplexNoise3D_FBM evaluates noise at the 3D block position.
				* The output is in [0, 1] (MinRange = 0, MaxRange = 1).
				* If the value exceeds CaveThreshold this voxel becomes hollow.
				*
				* Cave carving is disabled above CaveCeilingZ to preserve the surface.
				* Without this, caves would punch through the ground.
				*
				* CaveScale controls how wide the cave tunnels are.
				* Smaller scale = more open cave systems.
				* Larger scale = tighter, craggier tunnels.
				*/
				bool bIsCave = false;

				if (GridZ < CaveCeilingZ)
				{
					const double CaveNoise = USimplexNoiseBlueprintFunctionLibrary::SimplexNoise3D_FBM(
						static_cast<double>(GridX),
						static_cast<double>(GridY),
						static_cast<double>(GridZ),
						CaveScale,
						CaveOctaves,
						2.0,
						0.5,
						0.0,
						1.0);

					bIsCave = (CaveNoise > CaveThreshold);
				}

				/**
				* STEP 5: LAYER ASSIGNMENT
				*
				* ResolveBlockIndex figures out which block type from BlockTypes should appear here.
				* Cave voxels return INDEX_NONE (empty).
				*/
				const int32 BlockIndex = ResolveBlockIndex(ColumnHeight, GridZ, bIsCave, bIsBeachColumn, bIsSnowColumn);

				if (BlockIndex == INDEX_NONE)
				{
					continue;
				}

				if (!BlockComponents.IsValidIndex(BlockIndex) || !IsValid(BlockComponents[BlockIndex]))
				{
					continue;
				}

				/**
				* STEP 6: INSTANCING
				*
				* Add one instance at the local-space position of this block.
				* All instances in the same ISMC share a single draw call, so even thousands of blocks are rendered efficiently.
				*/
				BlockComponents[BlockIndex]->AddInstance(GetBlockTransform(GridX, GridY, GridZ), false);
			}
		}
	}

	UE_LOG(LogTemp, Log, TEXT("BiomeDioramaActor: Generation complete. Seed=%d  Size=%d x %d x %d"), NoiseSeed, SizeX, SizeY, SizeZ);
}

int32 AExampleDioramaActor::ResolveBlockIndex(const int32 ColumnHeight, const int32 BlockZ, const bool bIsCave, const bool bIsBeach, const bool bIsSnow) const
{
	// Cave voxels are always empty air
	if (bIsCave)
	{
		return INDEX_NONE;
	}

	const int32 DepthFromSurface = ColumnHeight - 1 - BlockZ;

	/**
	* Surface block: the very top of the terrain column.
	* The appearance depends on altitude and the biome blend value for this column.
	*
	*   Snow column  -> Snow role   (mountain peak)
	*   Beach column -> Sand role   (near sea level)
	*   Otherwise    -> Surface role (grass / default)
	*
	* If a preferred role is not found in BlockTypes, we fall back to the next
	* best option so the terrain is never invisible due to a missing block type.
	*/
	if (DepthFromSurface == 0)
	{
		if (bIsSnow)
		{
			const int32 SnowIndex = FindBlockIndexByRole(EBiomeBlockRole::Snow);
			if (SnowIndex != INDEX_NONE)
			{
				return SnowIndex;
			}
		}

		if (bIsBeach)
		{
			const int32 SandIndex = FindBlockIndexByRole(EBiomeBlockRole::Sand);
			if (SandIndex != INDEX_NONE)
			{
				return SandIndex;
			}
		}

		// Default surface (grass). Falls back to Subsurface then Deep if not found
		const int32 SurfaceIndex = FindBlockIndexByRole(EBiomeBlockRole::Surface);
		if (SurfaceIndex != INDEX_NONE)
		{
			return SurfaceIndex;
		}
	}

	// Subsurface band: a few blocks directly below the surface (dirt layer)
	// SubsurfaceDepth controls how thick this band is
	if (DepthFromSurface > 0 && DepthFromSurface <= SubsurfaceDepth)
	{
		const int32 SubIndex = FindBlockIndexByRole(EBiomeBlockRole::Subsurface);
		if (SubIndex != INDEX_NONE)
		{
			return SubIndex;
		}
	}

	/**
	* Deep block: stone or bedrock below the dirt layer.
	* This is the most common block type and fills the bulk of the terrain.
	*
	* We also try CaveWall here as a decorative variant for blocks just below the surface inside cave zones but since we already returned for confirmed caves above,
	* this only applies to narrow edge cases where the designer might want a transition material.
	* Kept for extensibility.
	*/
	const int32 DeepIndex = FindBlockIndexByRole(EBiomeBlockRole::Deep);
	if (DeepIndex != INDEX_NONE)
	{
		return DeepIndex;
	}

	// Last resort: return the first valid block type in the array so something always appears rather than leaving invisible holes in the terrain
	for (int32 FallbackIndex = 0; FallbackIndex < BlockTypes.Num(); ++FallbackIndex)
	{
		if (IsValid(BlockTypes[FallbackIndex].Mesh))
		{
			return FallbackIndex;
		}
	}

	return INDEX_NONE;
}

int32 AExampleDioramaActor::FindBlockIndexByRole(const EBiomeBlockRole InRole) const
{
	for (int32 Index = 0; Index < BlockTypes.Num(); ++Index)
	{
		if (BlockTypes[Index].Role == InRole && IsValid(BlockTypes[Index].Mesh))
		{
			return Index;
		}
	}

	return INDEX_NONE;
}

FTransform AExampleDioramaActor::GetBlockTransform(const int32 GridX, const int32 GridY, const int32 GridZ) const
{
	// Convert integer grid coordinates to world-space centimetres
	// BlockSizeCm should match the actual rendered size of your mesh (typically 100 cm)

	// bCenterBlocks shifts the position by half a block so the mesh pivot (if it is at the mesh center) lands exactly on the grid point
	// Leave false if your mesh pivot is at a corner or bottom face

	const double HalfBlock = bCenterBlocks ? BlockSizeCm * 0.5 : 0.0;

	const double CenterOffsetX = bCenterDiorama ? static_cast<double>(DioramaSize.X) * 0.5 : 0.0;
	const double CenterOffsetY = bCenterDiorama ? static_cast<double>(DioramaSize.Y) * 0.5 : 0.0;

	const FVector LocalPosition
	(
		(static_cast<double>(GridX) - CenterOffsetX) * BlockSizeCm + HalfBlock,
		(static_cast<double>(GridY) - CenterOffsetY) * BlockSizeCm + HalfBlock,
		static_cast<double>(GridZ) * BlockSizeCm + HalfBlock
	);

	return FTransform(FRotator::ZeroRotator, LocalPosition, FVector::OneVector);
}
