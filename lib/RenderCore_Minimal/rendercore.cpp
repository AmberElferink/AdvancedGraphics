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
#define THREADING

//#include "common_classes.h"

//  +-----------------------------------------------------------------------------+
//  |  RenderCore::Init                                                           |
//  |  Initialization.                                                      LH2'19|
//  +-----------------------------------------------------------------------------+
Timer t;
void RenderCore::Init()
{
	t.reset();

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

int glassIndex = 0;
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
			//parameters encapsulated in host_material.cpp
			uint metallic = (mat[i].parameters.x & 255);
			if (metallic >= 1)
				m->metallic = true;
			else
				m->metallic = false;

			m->color = make_float3(float(mat[i].diffuse_b), float(mat[i].diffuse_g), float(mat[i].diffuse_r));

			//if (glassIndex == 1)
			//{
			//	m->dielectric = true;
			//	m->metallic = false;
			//	m->indexOfRefraction = 1.05; //glass
			//	m->absorption = make_float3(0.8, 0.8, 0);
			//	m->color = make_float3(1, 1, 1);
			//}
			//glassIndex++;


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
void RenderCore::SetTarget(GLTexture *target)
{
	//accepts opengl texture
	// synchronize OpenGL viewport
	targetTextureID = target->ID;
	if (screen != 0 && target->width == screen->width && target->height == screen->height) return; // nothing changed
	delete screen;
	screen = new Bitmap(target->width, target->height);
	delete raytracer.buffer;
	raytracer.buffer = new Bitmap(screen->width, screen->height);
	printf("screenWidth: %u, screenHeight: %u \n", screen->width, screen->height);
}


//  +-----------------------------------------------------------------------------+
//  |  RenderCore::SetGeometry                                                    |
//  |  Set the geometry data for a model.                                   LH2'19|
//  +-----------------------------------------------------------------------------+
//meshIdx can be ignored. vertexData zitten x y en z in, de 4e is voor optimization, negeer die.
//triangle count is vertex count/3.
//core triangles: 3 vertices kun je uit de rest halen, daarmee kun je al een intersection doen, maar de normal ed zit in de coreTriangles
//copy data to new mesh, want je kan niet garanderen dat de data niet verandert terwijl je tekent. Daarom moet je een kopie maken.
void RenderCore::SetGeometry(const int meshIdx, const float4 *vertexData, const int vertexCount, const int triangleCount, const CoreTri *triangleData, const uint *alphaFlags)
{
	Timer timer;
	timer.reset();
	Mesh newMesh = Mesh();
	// copy the supplied vertices; we cannot assume that the render system does not modify
	// the original data after we leave this function.
	newMesh.vertices.resize(vertexCount);
	newMesh.vcount = vertexCount;
	for (int i = 0; i < vertexCount; i++)
	{
		newMesh.vertices[i] = vertexData[i];
	}
	// copy the supplied 'fat triangles'
	newMesh.triangles.resize( triangleCount );
	for ( int i = 0; i < triangleCount; i++ )
		newMesh.triangles[i] = triangleData[i];
	raytracer.scene.meshList.push_back(newMesh);
	raytracer.storeBVH();
	printf("loaded mesh with %i vertices in %5.3fs\n", vertexCount, timer.elapsed());
	//if (meshIdx == 3)
}


//  +-----------------------------------------------------------------------------+
//  |  RenderCore::SetLights                                                      |
//  |  Set the light data.                                                  LH2'19|
//  +-----------------------------------------------------------------------------+
void RenderCore::SetLights(const CoreLightTri *areaLights, const int areaLightCount,
	const CorePointLight *pointLights, const int pointLightCount,
	const CoreSpotLight *spotLights, const int spotLightCount,
	const CoreDirectionalLight *directionalLights, const int directionalLightCount)
{
	printf("setlights\n");
	for (int i = 0; i < pointLightCount; i++)
	{
		Light l;
		l.position = pointLights[i].position;
		l.radiance = pointLights[i].radiance;
		l.pointLight = true;
		raytracer.scene.lightList.push_back(l);
	}

	for (int i = 0; i < directionalLightCount; i++)
	{
		Light l;
		l.direction = directionalLights[i].direction;
		l.radiance = directionalLights[i].radiance;
		l.directionalLight = true;
		raytracer.scene.lightList.push_back(l);
	}

	for (int i = 0; i < spotLightCount; i++)
	{
		Light l;
		l.direction = spotLights[i].direction;
		l.position = spotLights[i].position;
		l.radiance = spotLights[i].radiance;
		l.cosInner = spotLights[i].cosInner;
		l.cosOuter = spotLights[i].cosOuter;
		l.spotLight = true;
		raytracer.scene.lightList.push_back(l);
	}

	for (int i = 0; i < areaLightCount; i++)
	{
		Light l;
		l.areaLight = true;
		l.triangle = areaLights[i];
		raytracer.scene.lightList.push_back(l);
	}
}

//  +-----------------------------------------------------------------------------+
//  |  RenderCore::Render                                                         |
//  |  Produce one image.                                                   LH2'19|
//  +-----------------------------------------------------------------------------+
int lineNr = 0;
int frameCounter = 0;

#ifdef THREADING
vector<thread> threads;
#endif
void RenderCore::Render(const ViewPyramid& view, const Convergence converge)
{
#ifdef THREADING
	threads.clear();
	t.reset();
	for (int i = 0; i < 4; i++)
	{
		threads.push_back(thread([=]() {
			raytracer.rayTraceBlock(view, screen, 0, i * (screen->height / 4), (i + 1) * (screen->height / 4));
		}));
	}

	for (int i = 0; i < 4; i++)
	{
		if (threads[i].joinable())
			threads[i].join();
	}
		printf("raytracer traced in %f\n", t.elapsed());

#else
	//---------per line raytracing --------

if (lineNr < screen->height)
{
	raytracer.rayTraceLine(screen, view, targetTextureID, lineNr);
	lineNr++;
	//printf("raytraced line in %f\n", t.elapsed());
}
else
{
	lineNr = 0;
	printf("raytraced in %f\n", t.elapsed());
	t.reset();
}

#endif

	//for (int i = 0; i < 4; i++)
	//{
	//	raytracer.rayTraceBlock(view, screen, 0, i * (screen->height / 4), (i + 1) * (screen->height / 4));
	//}


	//raytracer.rayTrace(screen, view, targetTextureID);



	// -----------per block raytracing ---------------
	//raytracer.rayTraceRandom(view, targetTextureID, frameCounter);
	//int screenSize = screen->width * screen->height;
	//for (int j = 0; j < screenSize; j++)
	//{
	//	screen->pixels[j] = raytracer.buffer->pixels[j] / frameCounter;
	//}



	glBindTexture(GL_TEXTURE_2D, targetTextureID);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, screen->width, screen->height, 0, GL_RGBA, GL_UNSIGNED_BYTE, screen->pixels);
}

//  +-----------------------------------------------------------------------------+
//  |  RenderCore::Shutdown                                                       |
//  |  Free all resources.                                                  LH2'19|
//  +-----------------------------------------------------------------------------+
void RenderCore::Shutdown()
{
#ifdef THREADING
	for (int i = 0; i < 4; i++)
	{
		if (threads[i].joinable())
			threads[i].join();
	}
#endif

	delete screen;
}

// EOF