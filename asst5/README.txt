CS175 Assignment 5 Submission

Winnie Wu - weiwu@college.harvard.edu
Jesse Chen - jesseyuanchen@college.harvard.edu

Files included:
  AUTHORS
  Debug
  LICENSE
  Makefile
  README.txt
  arcball.h
  asst4-snippets.cpp
  asst4.cpp
  asst5
  asst5-GL2.exe
  asst5-GL3.exe
  asstcommon.h
  cvec.h
  drawer.h
  geometrymaker.h
  glsupport.cpp
  glsupport.h
  matrix4.h
  picker.cpp
  picker.h
  ppm.cpp
  ppm.h
  quat.h
  rigtform.h
  scenegraph.cpp
  scenegraph.h
  sgutils.h
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
  quat.h
  rigtform.h

Platform used:
  Mac OS X 10.10

Compilation instructions:
  Run 'make' from the command line while in the main directory.

Problem set requirements:
  We followed the spec. We completed tasks 1 through 3 and built off from our solution to asst4. 
  In task1, we implemented basic keyframe infrastructure for the keys as described in the spec.
  In task2, we added the mathematical equations for lerp and slerp (and helpers) in the quat and rigtform .h files
  In task3, we added play functionality and interpolation with speed up and slow down as described in the spec.

Design:
  == TASK 1 == 
  - We stored keyframes as lists of RigTForm vectors and made a new type RigTFormVector for this purpose
  - RigTFormVectors are parallel to the SgRbtNode vector we get back from using the helper dump function from sgutils.h
  - The list of keyframes is stored as a global 'keyframeList'
  - We keep track of the current keyframe with a global g_currentKeyframe
  - To handle '<' and '>' we made a new enum constants ADVANCE and RETREAT for quick comparison
  - For read and write, we use a txt file. 
    The first line of the file is the number of keyframes (k).
    The second line of the file is the number of nodes per keyframe (n).
    The following k sets of n lines is the numerical representation of the RigTForm for each node in that keyframe.
      The first 3 numbers are the translational portion, the last 3 number describe rotation.
    There is no break between keyframe sets in the file, our read function takes care of parsing based on k and n.
   
  == TASK 2 ==
  - An equality (==) operator was added to quat.h to be used for conditional negation
  - A power (pow) function was added to quat.h to be used for slerping
  - Lerp and slerp were added as functions to rigtform.h

  == TASK 3 ==
  - The animation callback code from the spec was used
  - 2 global iterators are used to denote the frames we are interpolating between
  - Speed of animation and fps are both stored as globals
  - There is a maximum speed of animation (100ms break between keyframes) 
  - There is a minimum speed of animation (1000ms break between keyframes)
  - Changing speed of animation goes by 100ms between keyframes at a time
 
Usage:
  - Use 'h' to get help directions. 
  - All keyboard controls and mouse actions are as per spec guidelines. Also described in 'h'.
  
