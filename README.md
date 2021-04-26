# How to run?

This project consists from three assignments. You can switch between them via the branches.

If you just want to see the visuals, you can open up the imguiapp.exe in apps > imguiapp > imguiapp.exe
Sadly only Assignment 1 and 3 work like this currently on my computer a year after implementation, probably due to diverse updates. 

To run it from the code:
Due to an update on OpenGL and some Microsoft updates, Assigment 2 and 3 may give issues in compilation here.

If you clone the project, and open the .sln in Visual Studio. It may ask to Retarget the project. It is the safest to not Retarget it, and to install the version it wants in the Visual Studio Installer. It may work on your installed version, but it depends per version. 

In the Solution Explorer fold open render cores > rendercore_minimal. You can find the implemented code there. All files here are implemented by us. The other rendercores should be unloaded, except for rendercore_softrasterizer and rendercore_minimal, to prevent possible errors (right click and select "Unload"). Fold open applications, and set imguiapp as the startup project, and it should run. It may take a while for the image of assignment 1 to render, during which you will see the Lighthouse Technology demo screen, please give it a while depending on your CPU. Do not run it in Debug, since you will be waiting for a very long 



# Assignment 1 - Whitted style raytracer
![Raytraced - Whitted style -CPU 25s - A1](https://user-images.githubusercontent.com/32518317/115625724-c7a6fc80-a2fc-11eb-92ae-140999ff1a0d.png)
Report project 1, Advanced Graphics (2019)
Amber Elferink(5491525) & Ymke Wegereef (3971457)

We have developed a render core using lighthouse 2 with the following functionality:
* correct handling of the ‘ViewPyramid’ handed to you by the RenderSystem
* support for triangles as handed to you via SetGeometry
* a basic but extendible material class
* a Whitted-style ray tracing renderer, running on the CPU, supporting shadows and reflections, to demonstrate and test the core
* A basic but extendible texture class
* read in .png pictures as Texture and display the according texture on triangles
* Loading the metallic setting from the CoreMaterial. Since we could not figure out a
way to get a scene with correct translations with glass, we used another way to load
the glass.
* Handling spot- and directional lights
* Area lights
* Dielectrics: Snell, Fresnel and Beer
* Efficient and generic multi-threading, yielding a 205% speedup on dual core with
hyper threading.

## Division of the work
We have worked together on the first four (mandatory) points above. Once we had
implemented a basic rendercore, Amber started working on the dielectrics and the
multi-threading, while Ymke focussed on the different lights and the textures.
## Sources
* Our intersection method is adapted from the wikipedia page containing the
Möller–Trumbore intersection algorithm:
https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
* For the dielectric and reflection algorithms, the slides were used as a guide.

# Assignment 2 - BVH (Bounding Volume Hierarchy) raytracer
![Raytraced - Whitted style - CPU 180ms - A2](https://user-images.githubusercontent.com/32518317/115626838-56684900-a2fe-11eb-8dc3-60a6055605c3.png)

We implemented a BVH using lighthouse 2 on the Whitted Raytracer implemented in
assignment one with the following functionality:
* Binned building with the surface area heuristic (SAH)
The method “Binning” computes and compares results for each of the three axis. It
then chooses the best split plane based on the surface area heuristic. (There is also
a method called “Partition that greedily compares the centroids of the triangles,
based on their SAH, that returns as soon as an improvement is found.)
* Building a tree for the DeLorean (1082382 vertices in lighthouse = 360794 triangles)
is completed within 0.556s (single core, 3.6 GHz). This is measured with the
“Binning” method.
* Default BVH traversal, which renders the given DeLorean (360794 triangles,
1280x720px) in 180 ms quad core hyperthreaded (no threadpool) and in 2.45
seconds single core.

## Division of work
For this assignment we’ve mostly worked together, as it was harder to split up the work on
the BVH and as it was important for us both to fully understand the process. Hence we did a
lot of pair programming.

## Sources
* Ray-Box intersections:
https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-renderi
ng-simple-shapes/ray-box-intersection
* Our Binning method is based on the article: On fast Construction of SAH-based
Bounding Volume Hierarchies, Wald, 2007

# Report project 3 - Advanced Graphics (2020)
Amber Elferink(5491525) & Ymke Wegereef (3971457)

![Raytraced - Path tracer - packet - A3](https://user-images.githubusercontent.com/32518317/115626679-202ac980-a2fe-11eb-9efe-b0a4f28d5a6d.png)

We converted our existing ray tracer into a path tracer (using lighthouse 2 ) with the following
new functionality:
* Multiple importance sampling
* A photon map (with 1000000 photons) is build during initialization. Photons are
stored in a hashed grid that is build around the scene.
* Photon based Next Event Estimation as described in the master thesis of Andreas
Mikolajewski (see sources).
* Packet Traversal as described in Overbeck et al. in their paper “Large Ray Packets
for Real-time Whitted Ray Tracing”. This was implemented for the whitted raytracer
primary rays with frustra, for the pathtracer for the primary rays and the metallic rays.
The intersection methods for both the nodes and the primitives are fully AVX.
One can switch between ray and pathtracer within the core_settings.h, although the light is a
bit too strong in current scene for the raytracer.

## Division of work
For this assignment we first implemented a basic path tracer together. Once this was
finished, Ymke worked further on the importance sampling and PNEE and Amber worked on
packet traversal.

## Sources
* Efficient data structures and sampling of many light sources for Next Event
Estimation, by Andreas Mikolajewski
* Large Ray Packets for Real-time Whitted Ray Tracing, by Overbeck et al.


# About the framework: lighthouse2 (by Jacco Bikker, our teacher)
https://github.com/jbikker/lighthouse2
Lighthouse 2 framework for real-time ray tracing

This is the public repo for Lighthouse 2, a rendering framework for real-time ray tracing / path tracing experiments. 
Lighthouse 2 uses a state-of-the-art wavefront / streaming ray tracing implementation to reach high ray througput on RTX hardware 
(using Optix 7) and pre-RTX hardware (using Optix 5 Prime) and soon on AMD hardware (using RadeonRays / OpenCL) and CPUs (using Embree).
A software rasterizer is also included, mostly as an example of a minimal API implementation.

![Screenshot](https://user-images.githubusercontent.com/32518317/115627568-89f7a300-a2ff-11eb-9c1c-9c260dcb5ed7.png)

Quick pointers / Important advice:

* Building Lighthouse 2: CUDA 10 currently does *not* properly support vs2019; use vs2017 for now.
* Lighthouse 2 wiki: https://github.com/jbikker/lighthouse2/wiki (early stages)
* Trouble shooting page on the wiki: https://github.com/jbikker/lighthouse2/wiki/TroubleShooting
* Lighthouse 2 forum: https://ompf2.com/viewforum.php?f=18
* Follow the project on Twitter: @j_bikker

Lighthouse 2 uses a highly modular approach to ease the development of renderers.

The main layers are:

1. The application layer, which implements application logic and handles user input;
2. The RenderSystem, which handles scene I/O and host-side scene storage;
3. The render cores, which implement low-level rendering functionality.

Render cores have a common interface and are supplied to the RenderSystem as dlls. The RenderSystem supplies the cores with scene data 
(meshes, instances, triangles, textures, materials, lights) and sparse updates to this data.

The Lighthouse 2 project has the following target audience:

*Researchers*

Lighthouse 2 is designed to be a high-performance starting point for novel algorithms involving real-time ray tracing. This may include
new work on filtering, sampling, materials and lights. The provided ray tracers easily reach hundreds of millions of rays per second 
on NVidia and AMD GPUs. Combined with a generic GPGPU implementation, this enables a high level of freedom in the implementation of 
new code.

*Educators*

The Lighthouse 2 system implements all the boring things such as scene I/O, window management, user interfaces and access to ray tracing
APIs such as Optix, RadeonRays and Embree; your students can get straight to the interesting bits. The architecture of Lighthouse 2 is
carefully designed to be easily accessible. Very fast scene loading and carefully tuned project files ensure quick development cycles.

*Industry*

Lighthouse 2 is an R&D platform. It is however distributed with the Apache 2.0 license, which allows you to use the code in your
own products. Experimental cores can be shared with the community in binary / closed form, and application development is separated
from core development.

<b>What it is not</b>

The ray tracing infrastructure (with related scene management acceleration structure maintenance) should be close to optimal. The implemented estimators however (unidirectional path tracers without filtering and blue noise) are not, and neither is the shading
model (Lambert + speculars). This may or may not change depending on the use cases encountered. This video shows what can be
achieved with the platform: https://youtu.be/uEDTtu2ky3o .

Lighthouse 2 should compile out-of-the-box on Windows using Visual Studio 2017 / 2019. For the CUDA/Optix based cores CUDA 10.2 is required:

https://developer.nvidia.com/cuda-downloads

Optix 5.x, 6.0 and 7.0 libraries are included in the Lighthouse 2 download.

For more information on Lighthouse 2 please visit: http://jacco.ompf2.com.

<b>Credits</b>

Lighthouse 2 was developed at the Utrecht University, The Netherlands.

Lighthouse 2 uses the following libraries:<br>
Dear ImGui https://github.com/ocornut/imgui<br>
FreeImage http://freeimage.sourceforge.net<br>
Glad https://glad.dav1d.de<br>
GLFW https://www.glfw.org<br>
half 1.12 http://half.sourceforge.net<br>
tinygltf https://github.com/syoyo/tinygltf<br>
tinyobj https://github.com/syoyo/tinyobjloader<br>
tinyxml2 https://github.com/leethomason/tinyxml2<br>
zlib https://www.zlib.net

<b>Contributions</b>

* The Lighthouse2 Vulkan core (and sharedBSDF) was developed by Mèir Noordermeer (https://github.com/MeirBon).
* A Linux port by Marijn Suijten (https://github.com/MarijnS95) is being incorporated in the main repo.
* Animation code uses low-level optimizations by Alysha Bogaers and Naraenda Prasetya.
* OptixPrime_BDPT core by Guowei (Peter) Lu (https://github.com/pasu).

<b>Previous Work</b>

Lighthouse 2 implements research by (very incomplete):

* Marsaglia: random numbers
* Van Antwerpen, Laine, Karras, Aila: streaming path tracing
* Aila, Laine: persistent kernels
* Schied et al.: Spatiotemporal Variance-Guided Filtering (SVGF)
* Victor Voorhuis: improved SVGF for specular and glossy reprojection
* Eric Heitz: Blue noise distributions
