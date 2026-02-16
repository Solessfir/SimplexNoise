/**
* OVERVIEW
* --------
* This actor procedurally builds a small voxel diorama (a bounded 3D scene made of 1 x 1 x 1 metre cubes) using the SimplexNoise plugin.
* It is intended as a tutorial and showcase example demonstrating how different noise functions work together to produce natural-looking terrain.
*
* NOISE USAGE SUMMARY
* -------------------
*  - SimplexNoise2D_FBM  : Generates the terrain heightmap. Each (X, Y) column is sampled once to decide how tall that column of blocks should be.
*                          Multiple octaves stack detail at different scales, producing rolling hills, valleys, and peaks.
*
*  - SimplexNoise3D_FBM  : Carves caves through the terrain. Every block position (X, Y, Z) is sampled.
*                          If the returned value exceeds CaveThreshold the block is removed, hollowing out tunnels and caverns.
*
*  - SimplexNoise2D      : Generates a separate "biome blend" map that shifts surface blocks between grass,
*                          sand, and snow depending on moisture/temperature without requiring a full biome system.
*
* BLOCK TYPE ARRAY SETUP
* ----------------------
* Assign elements to the BlockTypes array in the Details panel. Each element has a Role that tells the generator when to use it:
*
*   Index | Role          | Suggested mesh / material
*   ------|---------------|---------------------------
*     0   | Surface       | Grass cube
*     1   | Subsurface    | Dirt cube
*     2   | Deep          | Stone cube
*     3   | Sand          | Sand cube (replaces Surface near sea level)
*     4   | Snow          | Snow cube (replaces Surface at high altitude)
*     5   | CaveWall      | Mossy stone / cave rock cube (optional accent)
*
* Only the roles present in the array are used.
* Roles that are missing fall back to the closest defined role (e.g. no Sand -> Surface is used on beaches instead).
*
* MESH REQUIREMENTS
* -----------------
* Use a Static Mesh that is exactly 100 x 100 x 100 cm (1 m^3) with its pivot at the bottom-left-front corner, or at the center - just set bCenterBlocks accordingly.
* Unreal's default Engine/BasicShapes/Cube.uasset works well when scaled to 1 m.
*
* QUICK START
* -----------
*  1. Place a BiomeDioramaActor in your level.
*  2. Populate the BlockTypes array (see table above).
*  3. Adjust DioramaSize, NoiseSeed, and any noise parameters.
*  4. Click "Generate Diorama" in the Details panel (or call GenerateDiorama() from Blueprint).
*  5. Iterate: change seed or parameters, click Generate again. ClearDiorama() removes all blocks.
*
* PERFORMANCE NOTE
* ----------------
* All blocks of the same type share one UInstancedStaticMeshComponent so the entire diorama renders in as few draw calls as there are distinct block types.
* A 32x32x24 diorama typically produces ~8 000 visible blocks; this is comfortable for editor use and lightweight enough for a real-time preview.
*/

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "Components/InstancedStaticMeshComponent.h"
#include "ExampleDioramaActor.generated.h"

/**
* Decides which layer of the terrain a block type belongs to.
* The generator picks the most appropriate role for each voxel automatically.
*/
UENUM(BlueprintType)
enum class EBiomeBlockRole : uint8
{
	Surface     UMETA(DisplayName = "Surface (e.g. Grass)"),
	Subsurface  UMETA(DisplayName = "Subsurface (e.g. Dirt)"),
	Deep        UMETA(DisplayName = "Deep (e.g. Stone)"),
	Sand        UMETA(DisplayName = "Sand (Beach near sea level)"),
	Snow        UMETA(DisplayName = "Snow (Mountain peaks)"),
	CaveWall    UMETA(DisplayName = "Cave Wall (inside cave ceilings)"),
};

/**
* Defines one type of block used in the diorama.
* Add one of these per material/appearance you want to see in the generated scene.
*/
USTRUCT(BlueprintType)
struct FBiomeBlockType
{
	GENERATED_BODY()

	// Human-readable label shown in the Details panel
	// Does not affect generation
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Block")
	FName BlockName = NAME_None;

	// Where in the terrain this block appears
	// The generator automatically places blocks based on depth, altitude, and proximity to sea level
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Block")
	EBiomeBlockRole Role = EBiomeBlockRole::Deep;

	// The static mesh to spawn for this block. Should be 100 x 100 x 100 cm
	// All blocks sharing this mesh are batched into one draw call via ISMC
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Block")
	TObjectPtr<UStaticMesh> Mesh;

	// Optional material override applied to the mesh
	// Leave empty to use the mesh's default material
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Block")
	TObjectPtr<UMaterialInterface> Material;
};

/**
* Procedural voxel diorama built with 1 m cubes and SimplexNoise.
*
* See the file header comment for a full description, noise breakdown,
* and quick-start guide.
*/
UCLASS(Blueprintable, BlueprintType)
class SIMPLEXNOISE_API AExampleDioramaActor : public AActor
{
	GENERATED_BODY()

public:
	AExampleDioramaActor();

protected:
	virtual void BeginPlay() override;

public:
	// -------------------------------------------------------------------------
	//  Editor Actions
	// -------------------------------------------------------------------------

	/** Clears any existing blocks and rebuilds the entire diorama from scratch using the current settings and seed. Safe to call multiple times. */
	UFUNCTION(BlueprintCallable, CallInEditor, Category = "Diorama|Generation")
	void GenerateDiorama();

	/** Destroys all spawned block instances and ISMC components. */
	UFUNCTION(BlueprintCallable, CallInEditor, Category = "Diorama|Generation")
	void ClearDiorama();

	// -------------------------------------------------------------------------
	//  Block Types
	// -------------------------------------------------------------------------

	/**
	* List of block types available to the generator.
	*
	* Recommended setup (add in this order for best results):
	*   [0] Surface    - Grass cube
	*   [1] Subsurface - Dirt cube
	*   [2] Deep       - Stone cube
	*   [3] Sand       - Sand cube
	*   [4] Snow       - Snow cube
	*   [5] CaveWall   - Cave rock cube (optional)
	*
	* The generator searches this array by Role, so order does not matter, but having all six roles filled produces the richest scene.
	*/
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Meta = (TitleProperty = "BlockName"), Category = "Diorama|Blocks")
	TArray<FBiomeBlockType> BlockTypes;

	// -------------------------------------------------------------------------
	//  World Shape
	// -------------------------------------------------------------------------

	/**
	* The number of blocks along each axis.
	* X = width (East), Y = depth (North), Z = height.
	* Default 32 x 32 x 24 gives a diorama roughly the size of a small room.
	* Increase with care: generation time and instance count scale as X * Y * Z.
	*/
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Diorama|Shape", Meta = (ClampMin = "4"))
	FIntVector DioramaSize = FIntVector(32, 32, 24);

	/**
	* Z level (in blocks) where water/sea begins. Blocks at or below this height
	* on low-lying terrain will be considered beach/sea-level and may receive
	* sand blocks if the Sand role is defined.
	*/
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Diorama|Shape", Meta = (ClampMin = "0"))
	int32 SeaLevel = 8;

	/**
	* Z blocks at and above this altitude receive snow surface blocks instead
	* of grass, simulating mountain peaks. Must be above SeaLevel.
	*/
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Diorama|Shape", Meta = (ClampMin = "0"))
	int32 SnowAltitude = 17;

	/**
	* How many blocks wide the beach band is around sea level.
	* Columns whose terrain height falls within SeaLevel + BeachWidth
	* get sand surface blocks instead of grass.
	*/
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Diorama|Shape", Meta = (ClampMin = "0"))
	int32 BeachWidth = 2;

	/**
	* How many dirt (Subsurface) blocks appear directly below the surface layer
	* before transitioning to stone (Deep).
	*/
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Diorama|Shape", Meta = (ClampMin = "1"))
	int32 SubsurfaceDepth = 3;

	// -------------------------------------------------------------------------
	//  Terrain Noise (2D FBM)
	// -------------------------------------------------------------------------

	/**
	* Seed used for all noise functions. Change this to get a completely
	* different landscape while keeping the same shape parameters.
	*/
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Diorama|Terrain Noise")
	int32 NoiseSeed = 42;

	/**
	* Controls how zoomed in the terrain features are. Lower values = broader,
	* smoother hills. Higher values = tighter, noisier terrain. Try 0.05 – 0.15.
	*/
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Diorama|Terrain Noise", Meta = (ClampMin = "0.001"))
	double TerrainScale = 0.07;

	/** Minimum terrain height in blocks. Should be above 0. */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Diorama|Terrain Noise", Meta = (ClampMin = "1"))
	int32 TerrainHeightMin = 4;

	/** Maximum terrain height in blocks. Should be below DioramaSize.Z. */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Diorama|Terrain Noise", Meta = (ClampMin = "2"))
	int32 TerrainHeightMax = 20;

	/**
	* Number of noise layers stacked for terrain. More octaves = more fine
	* surface detail (small bumps on top of large hills). 4 is a good default.
	*/
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Diorama|Terrain Noise", Meta = (ClampMin = "1", ClampMax = "8"))
	int32 TerrainOctaves = 4;

	/**
	* How much finer each terrain octave is compared to the previous one.
	* 2.0 = each layer is twice as detailed. Typical range: 1.8 – 2.5.
	*/
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Diorama|Terrain Noise", Meta = (ClampMin = "1.0"))
	double TerrainLacunarity = 2.0;

	/**
	* How much quieter each terrain octave is compared to the previous one.
	* 0.5 = each layer has half the influence. Typical range: 0.4 – 0.6.
	* Higher values make the finer layers more prominent (rougher terrain).
	*/
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Diorama|Terrain Noise", Meta = (ClampMin = "0.01", ClampMax = "1.0"))
	double TerrainPersistence = 0.5;

	// -------------------------------------------------------------------------
	//  Cave Noise (3D FBM)
	// -------------------------------------------------------------------------

	/**
	* Scale for the 3D cave noise. Lower values = larger, more open caverns.
	* Higher values = smaller, tighter tunnels. Try 0.08 – 0.2.
	*/
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Diorama|Cave Noise", Meta = (ClampMin = "0.001"))
	double CaveScale = 0.12;

	/**
	* Number of FBM octaves used for cave carving. Fewer octaves = simpler,
	* more tube-like caves. More octaves = rough, rocky cavern walls.
	*/
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Diorama|Cave Noise", Meta = (ClampMin = "1", ClampMax = "6"))
	int32 CaveOctaves = 3;

	/**
	* A block is carved out (becomes a cave) when its 3D noise value exceeds
	* this threshold. The noise output is in [0, 1] after remapping.
	* 0.65 = roughly 35% of deep blocks are hollow. Try 0.55 – 0.75.
	* Lower values make caves larger and more frequent.
	*/
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Diorama|Cave Noise", Meta = (ClampMin = "0.0", ClampMax = "1.0"))
	double CaveThreshold = 0.65;

	/**
	* Caves will not appear above this Z level (in blocks).
	* This prevents caves from carving through the surface layer and creating
	* ugly holes in the top of the terrain. Should be a few blocks below SeaLevel.
	*/
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Diorama|Cave Noise", Meta = (ClampMin = "0"))
	int32 CaveCeilingZ = 6;

	// -------------------------------------------------------------------------
	//  Biome Blend Noise (2D)
	// -------------------------------------------------------------------------

	/**
	* Scale for the 2D biome blend noise that shifts surface blocks between
	* grass, sand, and snow across the diorama. Lower = broader biome patches.
	* This is sampled independently of the terrain height noise.
	*/
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Diorama|Biome Noise", Meta = (ClampMin = "0.001"))
	double BiomeBlendScale = 0.04;

	// -------------------------------------------------------------------------
	//  Placement
	// -------------------------------------------------------------------------

	/**
	* Size of each cube in Unreal units (centimetres). Set to 100 if your mesh
	* is a 1 m cube scaled to default. Change to match your mesh size.
	*/
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Diorama|Placement", Meta = (ClampMin = "1.0"))
	double BlockSizeCm = 100.0;

	/**
	* When true, block instances are placed so their center aligns to the grid
	* (offset by half a block). When false, the bottom face of each block sits
	* on the grid. Match this to your mesh's pivot point.
	*/
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Diorama|Placement")
	bool bCenterBlocks = true;

	/**
	* When true, the diorama is centered on the actor's origin along X and Y.
	* When false, the diorama extends in the +X +Y direction from the origin.
	*/
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Diorama|Placement")
	bool bCenterDiorama = true;;

private:

	/**
	* One ISMC is created per block type during generation. All blocks sharing
	* the same mesh and material are batched into a single component for efficiency.
	* These components are owned by this actor and destroyed in ClearDiorama().
	*/
	UPROPERTY()
	TArray<TObjectPtr<UInstancedStaticMeshComponent>> BlockComponents;

	/**
	* Finds the first block type in BlockTypes whose Role matches the given role.
	* Returns INDEX_NONE if no matching entry exists.
	*/
	int32 FindBlockIndexByRole(EBiomeBlockRole InRole) const;

	/**
	* Given a column height, the current Z level, and whether the point is
	* inside a cave, returns the block role that should be placed here - or
	* INDEX_NONE if the voxel should remain empty (air or cave).
	*
	* This is the core of the terrain layer logic. Reading this function is a
	* good place to start when customizing block placement rules.
	*/
	int32 ResolveBlockIndex(int32 ColumnHeight, int32 BlockZ, bool bIsCave, bool bIsBeach, bool bIsSnow) const;

	/**
	* Returns the world-space transform for a block at integer grid position
	* (GridX, GridY, GridZ), accounting for BlockSizeCm and bCenterBlocks.
	*/
	FTransform GetBlockTransform(int32 GridX, int32 GridY, int32 GridZ) const;
};
