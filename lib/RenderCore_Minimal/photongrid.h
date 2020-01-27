#pragma once

class Cell
{
  public:
	vector<long> photonCount;
	vector<float3> dir;
	long total = 0;
	bool p_dist = false;
	vector<float> probs;
	void CalculateProbs();
	long dirCount = 0;
};

//A 3d hashed grid that resembles a photon map
class PhotonGrid
{
  public:
	vector<vector<vector<Cell>>> grid; //3D
	uint lightCount;
	uint k; //dimension size
	aabb boundingBox;
	float x_length, y_length, z_length;
	void InitializeGrid( uint dim, uint lc, aabb bounds );
	void StorePhoton( float3 position, uint lightID );
	void StoreDirection( float3 position, float3 dir );
	void LookupCell( uint &x, uint &y, uint &z, const float3 &position );
	float3 pickDirection( float3 position, float3 norm );
	float3 CosineWeightedDiffuseReflection( float3 N );
};