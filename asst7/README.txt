CS175 Assignment 7 Submission

Winnie Wu - weiwu@college.harvard.edu
Jesse Chen - jesseyuanchen@college.harvard.edu

Files included:
  animation.txt
  arcball.h
  asst7
  asst7.cpp
  asstcommon.h
  AUTHORS
  cube.mesh
  cvec.h
  drawer.h
  Fieldstone.ppm
  FieldstoneNormal.ppm
  geometry.cpp
  geometry.h
  geometrymaker.h
  glsupport.cpp
  glsupport.h
  LICENSE
  Makefile
  material.cpp
  material.h
  matrix4.h
  mesh.h
  mesh.interface
  picker.cpp
  picker.h
  ppm.cpp
  ppm.h
  quat.h
  README.txt
  renderstates.cpp
  renderstates.h
  rigtform.h
  scenegraph.cpp
  scenegraph.h
  sgutils.h
  texture.cpp
  texture.h
  uniforms.h
  shaders
    basic-gl2.vshader
    basic-gl3.vshader
    diffuse-gl2.fshader 
    diffuse-gl3.fshader
    normal-gl2.fshader
    normal-gl2.vshader
    normal-gl3.fshader
    normal-gl3.vshader
    pick-gl2.fshader
    pick-gl3.fshader
    solid-gl2.fshader
    solid-gl3.fshader 
    specular-gl2.fshader
    specular-gl3.fshader

Modified files from starter code:
  asst7.cpp

Platform used:
  Mac OS X 10.10

Compilation instructions:
  Run 'make' from the command line while in the main directory.

Problem set requirements:
  We followed the spec. We completed Assignment 7.
  In Assignment 7, we implemented Catmull-Clark subdivision on a mesh data structure. All tasks
  are completed and match the spec. 

Design:
  == Assignment 7 == 
  
  == Task 1 == 
  We create a shared_ptr to a SimpleGeometryPN, called g_subdivSurface; this gets the vertices
  that get drawn on the screen with the upload function. 

  == Task 2 == 
  We use a global flag g_smoothShadeOn to tell our getMeshVertices function to set the normals 
  in each face. 

  == Task 3 == 
  In this task we use glutTimerFunc and register our function animateSubdivSurfaceCallback, scaling
  off of the global variable g_subdivSurfaceAnimateSpeed, to allow for animation of the mesh vertices.

  == Task 4 == 
  We use Catmull-Clark rules to get the new subdivided face-vertices, edge-vertices, and vertex-vertices.
 
Usage:
  - Use 'h' to get help directions. 
  - All keyboard controls and mouse actions are as per spec guidelines. Also described in 'h'.
  
