/* rendercore.cpp - Copyright 2019 Utrecht University

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

	   http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/
#include "core_settings.h"



//#include "common_classes.h"

//  +-----------------------------------------------------------------------------+
//  |  RenderCore::Init                                                           |
//  |  Initialization.                                                      LH2'19|
//  +-----------------------------------------------------------------------------+
void RenderCore::Init()
{

}

void RenderCore::SetMaterials(CoreMaterial* mat, const CoreMaterialEx* matEx, const int materialCount) // textures must be in sync when calling this
{
	// copy the supplied array of materials
	for (int i = 0; i < materialCount; i++)
	{
		Material* m;
		if (i < matList.size()) m = matList[i];
		else matList.push_back(m = new Material());
		m->texture = 0;
		int texID = matEx[i].texture[TEXTURE0];
		if (texID == -1)
		{
			float r = mat[i].diffuse_r, g = mat[i].diffuse_g, b = mat[i].diffuse_b;
			m->diffuse = ((int)(b * 255.0f) << 16) + ((int)(g * 255.0f) << 8) + (int)(r * 255.0f);
		}
		else
		{
			m->texture = texList[texID];
			m->texture->width = mat[i].texwidth0; // we know this only now, so set it properly
			m->texture->height = mat[i].texheight0;
		}
	}
}

//  +-----------------------------------------------------------------------------+
//  |  RenderCore::SetTarget                                                      |
//  |  Set the OpenGL texture that serves as the render target.             LH2'19|
//  +-----------------------------------------------------------------------------+
void RenderCore::SetTarget( GLTexture *target )
{
	//accepts opengl texture
	// synchronize OpenGL viewport
	targetTextureID = target->ID;
	if ( screen != 0 && target->width == screen->width && target->height == screen->height ) return; // nothing changed
	delete screen;
	screen = new Bitmap( target->width, target->height );
}

//  +-----------------------------------------------------------------------------+
//  |  RenderCore::SetGeometry                                                    |
//  |  Set the geometry data for a model.                                   LH2'19|
//  +-----------------------------------------------------------------------------+
//meshIdx can be ignored. vertexData zitten x y en z in, de 4e is voor optimization, negeer die.
//triangle count is vertex count/3.
//core triangles: 3 vertices kun je uit de rest halen, daarmee kun je al een intersection doen, maar de normal ed zit in de coreTriangles
//copy data to new mesh, want je kan niet garanderen dat de data niet verandert terwijl je tekent. Daarom moet je een kopie maken.
void RenderCore::SetGeometry( const int meshIdx, const float4 *vertexData, const int vertexCount, const int triangleCount, const CoreTri *triangleData, const uint *alphaFlags )
{
	Mesh newMesh;
	// copy the supplied vertices; we cannot assume that the render system does not modify
	// the original data after we leave this function.
	newMesh.vertices = new float4[vertexCount];
	newMesh.vcount = vertexCount;
	memcpy( newMesh.vertices, vertexData, vertexCount * sizeof( float4 ) );
	// copy the supplied 'fat triangles'
	newMesh.triangles = new CoreTri[vertexCount / 3];
	memcpy( newMesh.triangles, triangleData, ( vertexCount / 3 ) * sizeof( CoreTri ) );
	meshes.push_back( newMesh );
}

//  +-----------------------------------------------------------------------------+
//  |  RenderCore::Render                                                         |
//  |  Produce one image.                                                   LH2'19|
//  +-----------------------------------------------------------------------------+
// its not a raytracer, maar gaat over de vertices in de meshes en tekent een 2d versie of the level met dots voor je vertices. Zodra hij klaar is, tekent hij het naar de texture die wordt gerendert.
void RenderCore::Render( const ViewPyramid &view, const Convergence converge, const float brightness, const float contrast )
{
	// render
	screen->Clear();

	for ( int i = 0; i < screen->width; i++ ) //TODO niet in loop berekenen
	{
		for ( int j = 0; j < screen->height; j++ )
		{
			float u = (float)i / (float)screen->width;
			float v = (float)j / (float)screen->height;
			float3 P = view.p1 + u * ( view.p2 - view.p1 ) + v * ( view.p3 - view.p1 );

			float3 dir = P - view.pos;
			float3 D = dir / length( dir );

			Ray ray = Ray( view.pos, D );

			Intersection closest;
			closest.t = 10000000;
			closest.material = 0xffffff;

			for (Mesh &mesh : meshes)
			{
			
				int verticeCount = mesh.vcount / 3;
				for (int i = 0; i < verticeCount; i++) //mesh.vcount / 3; i++ ) // TODO
				{
					Intersection intersection;
					if (Intersect(ray, mesh.triangles[i], intersection))
					{
						if (intersection.t < closest.t)
							closest = intersection;
					}
				}
			}


			screen->pixels[i + j * screen->width] = closest.material;
		}
	}

	//	for( Mesh& mesh : meshes ) for( int i = 0; i < mesh.vcount; i++ )
	//	{
	//
	//		// convert a vertex position to a screen coordinate
	//		int screenx = mesh.vertices[i].x / 80 * (float)screen->width + screen->width / 2;
	//		int screeny = mesh.vertices[i].z / 80 * (float)screen->height + screen->height / 2;
	//		screen->Plot( screenx, screeny, 0xffffff /* white */ );
	//}
	// copy pixel buffer to OpenGL render target texture
	glBindTexture( GL_TEXTURE_2D, targetTextureID );
	glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, screen->width, screen->height, 0, GL_RGBA, GL_UNSIGNED_BYTE, screen->pixels );
}

//  +-----------------------------------------------------------------------------+
//  |  RenderCore::Shutdown                                                       |
//  |  Free all resources.                                                  LH2'19|
//  +-----------------------------------------------------------------------------+
void RenderCore::Shutdown()
{
	delete screen;
}

// EOF