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

//  +-----------------------------------------------------------------------------+
//  |  RenderCore::SetTextures                                                    |
//  |  Set the texture data.                                                LH2'19|
//  +-----------------------------------------------------------------------------+
void RenderCore::SetTextures(const CoreTexDesc* tex, const int textures)
{
	Timer timer;
	timer.reset();
	
	// copy the supplied array of texture descriptors
	for (int i = 0; i < textures; i++)
	{
		Texture* t;
		if (i < raytracer.scene.texList.size()) t = raytracer.scene.texList[i];
		else  raytracer.scene.texList.push_back(t = new Texture());
		t->pixels = (uint*)MALLOC64(tex[i].pixelCount * sizeof(uint));
		if (tex[i].idata) memcpy(t->pixels, tex[i].idata, tex[i].pixelCount * sizeof(uint));
		else memcpy(t->pixels, 0, tex[i].pixelCount * sizeof(uint) /* assume integer textures */);
		// Note: texture width and height are not known yet, will be set when we get the materials.
	}
	printf("loaded textures in %5.3fs\n", timer.elapsed());
}

void RenderCore::SetMaterials(CoreMaterial* mat, const CoreMaterialEx* matEx, const int materialCount) // textures must be in sync when calling this
{
	Timer timer;
	timer.reset();
	// copy the supplied array of materials
	for (int i = 0; i < materialCount; i++)
	{
		Material* m;
		if (i < raytracer.scene.matList.size()) m = raytracer.scene.matList[i];
		else  raytracer.scene.matList.push_back(m = new Material());
		m->texture = 0;
		int texID = matEx[i].texture[TEXTURE0];
		if (texID == -1)
		{
			m->diffuse = make_float3(mat[i].diffuse_r, mat[i].diffuse_g, mat[i].diffuse_b);
		}
		else
		{
			m->texture = raytracer.scene.texList[texID];
			m->texture->width = mat[i].texwidth0; // we know this only now, so set it properly
			m->texture->height = mat[i].texheight0;
		}
	}
	printf("loaded materials in %5.3fs\n", timer.elapsed());
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

bool firsttime = true;
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

	if (firsttime) //meshes were loaded two times. For now only load it the first time the program is started.
	{
		Timer timer;
		timer.reset();
		Mesh newMesh;
		// copy the supplied vertices; we cannot assume that the render system does not modify
		// the original data after we leave this function.
		newMesh.vertices = new float4[vertexCount];
		newMesh.vcount = vertexCount;
		memcpy(newMesh.vertices, vertexData, vertexCount * sizeof(float4));
		// copy the supplied 'fat triangles'
		newMesh.triangles = new CoreTri[vertexCount / 3];
		memcpy(newMesh.triangles, triangleData, (vertexCount / 3) * sizeof(CoreTri));
		raytracer.scene.meshList.push_back(newMesh);
		printf("loaded geometry in %5.3fs\n", timer.elapsed());
	}
	firsttime = false;

}

//  +-----------------------------------------------------------------------------+
//  |  RenderCore::SetLights                                                      |
//  |  Set the light data.                                                  LH2'19|
//  +-----------------------------------------------------------------------------+
void RenderCore::SetLights( const CoreLightTri *areaLights, const int areaLightCount,
							const CorePointLight *pointLights, const int pointLightCount,
							const CoreSpotLight *spotLights, const int spotLightCount,
							const CoreDirectionalLight *directionalLights, const int directionalLightCount )
{
	for ( int i = 0; i < pointLightCount; i++ )
	{
		Light l;
		l.position = pointLights[i].position;
		l.radiance = pointLights[i].radiance;
		l.energy = pointLights[i].energy;
		raytracer.scene.lightList.push_back( l );
	}
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
	Timer t;
	t.reset();
	raytracer.rayTrace( screen, view, targetTextureID );
	printf("raytraced in %5.3fs\n", t.elapsed());

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