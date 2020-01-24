#include "core_settings.h"

/* This class contains methods that help building a photon map */

/* Method that computes a probability distribution for the different area lights */
void Cell::CalculateProbs()
{
	uint nlights = photonCount.size(); 
	probs.resize( nlights );

	if (total > 10)
	{
		long div = total + nlights;
		for (int i = 0; i < nlights; i++)
			probs[i] = (float)(1 + photonCount[i]) / div; //every light has to have a small probability to be sampled
	}
	else //not enough photons to make a good estimate
	{
		for ( int i = 0; i < nlights; i++ )
			probs[i] = (float) 1 / nlights;
	}
	p_dist = true;
}

/* Method that initializes a grid used for the storage of a photon map. */
void PhotonGrid::InitializeGrid(uint dim, uint lc, aabb bounds)
{
	k = dim;
	dim++;
	lightCount = lc;
	boundingBox = bounds;
	x_length = bounds.Maximum( 0 ) - bounds.Minimum( 0 );
	y_length = bounds.Maximum( 1 ) - bounds.Minimum( 1 );
	z_length = bounds.Maximum( 2 ) - bounds.Minimum( 2 );
    Cell emptyCell;
	emptyCell.photonCount.resize( lc );
	for (int i = 0; i < lc; i++)
	{
		emptyCell.photonCount[i] = 0;
		emptyCell.total = 0;
	}
	grid.resize( dim );
	for (int x = 0; x < dim; x++)
	{
		grid[x].resize( dim );
		for (int y = 0; y < dim; y++)
		{
			grid[x][y].resize( dim );
			for ( int z = 0; z < dim; z++ )
				grid[x][y][z] = emptyCell;
		}
	}
}

/* Method that adds a new photon to the photon map */
void PhotonGrid::StorePhoton(float3 position, uint lightID)
{
	uint x, y, z;
	LookupCell( x, y, z, position );
	grid[x][y][z].total++;
	grid[x][y][z].photonCount[lightID]++;
}

/* Look up the cell in which a point lies according to its coordinates */
void PhotonGrid::LookupCell(uint &x, uint &y, uint &z, const float3 &position)
{
	x = ( ( position.x - boundingBox.bmin[0] ) / x_length ) * k;
	y = ( ( position.y - boundingBox.bmin[1] ) / y_length ) * k;
	z = ( ( position.z - boundingBox.bmin[2] ) / z_length ) * k;
}
