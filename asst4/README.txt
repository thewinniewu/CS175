CS175 Assignment 2 Submission

Winnie Wu - weiwu@college.harvard.edu
Jesse Chen - jesseyuanchen@college.harvard.edu

Files included:
  AUTHORS
  asst2.cpp
  asst2.vcxproj.filters
  glsupport.cpp
  ppm.cpp
  LICENSE
  asst2.sln
  cvec.h
  glsupport.h
  ppm.h
  Makefile
  asst2.vcxproj
  geometrymaker.h
  matrix4.h
  shaders/
    basic-gl2.vshader
    basic-gl3.vshader
    diffuse-gl2.fshader 
    diffuse-gl3.fshader
    solid-gl2.fshader
    solid-gl3.fshader 

Modified files from starter code:
  asst2.cpp
  matrix4.h 

Platform used:
  Mac OS X 10.9.5

Compilation instructions:
  Run 'make' from the command line while in the main directory.

Problem set requirements:
  We followed the spec. There is a weird behavior if you are acting as one of the cubes and rotating yourself, since it
  uses its own axes for rotation, it does not behave like a camera would be expected to. But, this is as per the spec,
  so we didn't change anything (we reversed the direction of the rotation as the spec directed).

Design:
  - Added a second blue cube to the objectRbt array (with the blue color being in the parallel color array) 
  - Added functions to control different keyboard behaviors (for 'o', 'v', and 'm').
  - Kept track of the current view (controlled by 'v') with a global
  - Kept track of the current object being manipulated (controlled by 'o') with a global
  - Added cases to motion function to control the different cases of rotation/translation directions
  - Used linFact and transFact to compose hybrid auxFrames for transformations.

Usage:
  - Use 'h' to get help directions. 
  - All keyboard controls and mouse actions are as per spec guidelines. Also described in 'h'.
  - The default when the program is opened is skycamera view with object being manipulated being
    the red cube. The default when the skycamera is the object being manipulated is the sky sky aFrame,
    which can be changed with 'm' to world sky. 
  - The sequence of control rotation is red cube -> blue cube -> sky camera with both 'v' and 'o'. 
  - Play around. Tested methodically rotation/translation of red/blue cube from skycamera view, 
    rotating/translating skycamera itself, switching to worldsky auxFrame and then rotating/translating
    skycamera, changing view to red cube and rotating/translating blue cube, rotating/translating red
    cube itself, and then repeat for blue cube.
 


