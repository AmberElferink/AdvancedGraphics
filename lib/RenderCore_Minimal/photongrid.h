#pragma once

class Cell
{
  public: 
	vector<long> photonCount;
	long total = 0;
	bool p_dist = false;
	vector<float> probs;
	void CalculateProbs();
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
	void InitializeGrid(uint dim, uint lc, aabb bounds);
	void StorePhoton( float3 position, uint lightID );
	void LookupCell( uint &x, uint &y, uint &z, const float3 &position );
};