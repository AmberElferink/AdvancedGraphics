#include "core_settings.h"

/* This class contains methods that help building a photon map */

/* Method that computes a probability distribution for the different area lights */
void Cell::CalculateProbs()
{
	uint nlights = photonCount.size();
	probs.resize( nlights );

	if ( total > 10 )
	{
		long div = total + nlights;
		for ( int i = 0; i < nlights; i++ )
			probs[i] = (float)( 1 + photonCount[i] ) / div; //every light has to have a small probability to be sampled
	}
	else //not enough photons to make a good estimate
	{
		for ( int i = 0; i < nlights; i++ )
			probs[i] = (float)1 / nlights;
	}
	p_dist = true;
}

/* Method that initializes a grid used for the storage of a photon map. */
void PhotonGrid::InitializeGrid( uint dim, uint lc, aabb bounds )
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
	for ( int i = 0; i < lc; i++ )
	{
		emptyCell.photonCount[i] = 0;
		emptyCell.total = 0;
	}
	grid.resize( dim );
	for ( int x = 0; x < dim; x++ )
	{
		grid[x].resize( dim );
		for ( int y = 0; y < dim; y++ )
		{
			grid[x][y].resize( dim );
			for ( int z = 0; z < dim; z++ )
				grid[x][y][z] = emptyCell;
		}
	}
}

/* Method that adds a new photon to the photon map */
void PhotonGrid::StorePhoton( float3 position, uint lightID )
{
	uint x, y, z;
	LookupCell( x, y, z, position );
	grid[x][y][z].total++;
	grid[x][y][z].photonCount[lightID]++;
}

void PhotonGrid::StoreDirection( float3 position, float3 direction )
{
	uint x, y, z;
	LookupCell( x, y, z, position );
	grid[x][y][z].dir.push_back( direction );
	grid[x][y][z].dirCount++;
}

/* Look up the cell in which a point lies according to its coordinates */
void PhotonGrid::LookupCell( uint &x, uint &y, uint &z, const float3 &position )
{
	x = ( ( position.x - boundingBox.bmin[0] ) / x_length ) * k;
	y = ( ( position.y - boundingBox.bmin[1] ) / y_length ) * k;
	z = ( ( position.z - boundingBox.bmin[2] ) / z_length ) * k;
}

/* Look up a random direction that is stored in the desired grid cell */
float3 PhotonGrid::pickDirection( float3 position, float3 norm )
{
	uint x, y, z;
	LookupCell( x, y, z, position );
	Cell cell = grid[x][y][z];
	if ( cell.dirCount < 1 )
		return CosineWeightedDiffuseReflection( norm );
	uint randId = cell.dirCount * ((float)rand() / (float)RAND_MAX);
	if ( dot(cell.dir[randId],norm) < 0) 
		return CosineWeightedDiffuseReflection( norm );
	float3 dir = CosineWeightedDiffuseReflection( cell.dir[randId] );
	if ( dot( dir, norm ) < 0 || dot( cell.dir[randId],norm ) < 0 )
		return CosineWeightedDiffuseReflection( norm );
	else
		return dir;
}

/* Random vector on the hemisphere following a cosine weighted distribution */
float3 PhotonGrid::CosineWeightedDiffuseReflection( float3 N )
{
	float r0 = (float)rand() / (float)RAND_MAX;
	float r1 = (float)rand() / (float)RAND_MAX;
	float r = sqrt( r0 );
	float theta = 2 * PI * r1;
	float x = r * cosf( theta );
	float y = r * sinf( theta );
	float3 P = make_float3( x, y, sqrt( 1 - r0 ) );

	//Vector P needs to be transformed to tangent space
	float3 W;
	if ( abs( N.x ) > 0.99 )
		W = make_float3( 0, 1, 0 );
	else
		W = make_float3( 1, 0, 0 );
	float3 T = cross( N, W );
	T = T / length( T );
	float3 B = cross( T, N );
	//return make_float3( dot( P, T ), dot( P, B ), dot( P, N ) );
	return P.x * T + P.y * B + P.z * N;
}
