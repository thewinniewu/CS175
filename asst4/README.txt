CS175 Assignment 4 Submission

Winnie Wu - weiwu@college.harvard.edu
Jesse Chen - jesseyuanchen@college.harvard.edu

Files included:
  arcball.h
  asst3.sdf
  asst3.sln
  asst3.v12.suo
  asst3.vcxproj
  asst3.vcxproj.filters
  asst4-GL2.exe
  asst4-GL3.exe
  asst4-mac
  asst4-snippets.cpp
  asst4.cpp
  asstcommon.h
  AUTHORS
  cvec.h
  Debug
  drawer.h
  geometrymaker.h
  glsupport.cpp
  glsupport.h
  LICENSE
  Makefile
  matrix4.h
  output.txt
  picker.cpp
  picker.h
  ppm.cpp
  ppm.h
  quat.h
  README.txt
  rigtform.h
  scenegraph.cpp
  scenegraph.h
  shaders/
    basic-gl2.vshader
    basic-gl3.vshader
    diffuse-gl2.fshader 
    diffuse-gl3.fshader
    pick-gl2.fshader
    pick-gl3.fshader
    solid-gl2.fshader
    solid-gl3.fshader 

Modified files from starter code:
  asst4.cpp
  scenegraph.cpp
  picker.cpp 

Platform used:
  Mac OS X 10.9.5

Compilation instructions:
  Run 'make' from the command line while in the main directory.

Problem set requirements:
  We followed the spec. We completed tasks 1 through 4; in task 1 we added the code from asst4-snippets.cpp 
  to add the torso and right arms of the robot. In task 2 we wrote the necessary code for picker.cpp to allow
  the user to select different parts of the robot. In task 3 we wrote the necessary code for scenegraph and 
  modified our asst4.cpp code such that the user could actually move the robots. In task 4, we added the parts
  of the robot as directed. We started with the solution set from asst3. 

Design:
  - Created a function getNodeForEye that takes an ObjId as input and returns the corresponding Node
  - Used global struct ManipMode to keep track of the location of the arcball
  - Chose to call certain specialized cases of functions (getting nodes) inline instead of creating 
  auxiliary functions because they are part of the class

Usage:
  - Use 'h' to get help directions. 
  - All keyboard controls and mouse actions are as per spec guidelines. Also described in 'h'.
  - The default when the program is opened is skycamera view with the worldview being manipulated.
  - The sequence of viewing is red robot -> blue robot -> sky camera with 'v'. 
  - Play around. Allows the user to press 'p' and left click on part of a robot, which moves the arcball
  there and allows movement. Switching views brings you to the torso of the robot, from which you 
  are allowed to view yourself and pick the other robot as well. 
 


