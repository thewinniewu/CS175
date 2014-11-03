CS175 Assignment 6_5 Submission

Winnie Wu - weiwu@college.harvard.edu
Jesse Chen - jesseyuanchen@college.harvard.edu

Files included:
  animation.txt
  arcball.h
  asst5.cpp
  asstcommon.h
  AUTHORS
  cvec.h
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
  sgutils.h
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

Modified files from starter code:
  asst4.cpp
  quat.h
  normal-gl3.fshader

Platform used:
  Mac OS X 10.10

Compilation instructions:
  Run 'make' from the command line while in the main directory.

Problem set requirements:
  We followed the spec. We completed Assignments 6 and 6.5. 
  In Assignment 6, we replaced our interpolation from Assignment 5 with Catmull-Rom interpolation for smooth animation.
  In Assignment 6.5, we completed tasks 1 and 2 (Material Infrastructure and Bump Mapping).
  In task 1, we implement the code infrastructure fo different materials. 
  In task 2, we made the light movable and read in the map bump info

Design:
  == Assignment 6 == 
  - We replace the setRbt function from the staff solution with our own Catmull-Rom evaluation function
  - Our Catmull Rom evaluation function takes in four key frames and a time, computes the d and e values for the four
  key frames, and passes them to the Bezier Evaluation function
  - Our Bezier Evaluation function takes the two keyframes we want to interpolate between and the d and e values for them;
  it then uses lerp and slerp to compute the proper RigTForm interpolation and returns it
  - We implement a negation function, cn, in quat.h that negates a quat. We use it to negate what we pass to the power 
  function in d if the first argument of that term is negative (as per lecture notes)
  - In summary, we created Bezier controls to interpolate between two keyframes using the one following and one previous 
  keyframes to the two keyframes; we do this to evaluate the Bezier curves at the times between two keyframes for a smoother
  animation 
   
  == Assigment 6.5 Part 1 ==
  - We migrated the code and followed the instructions in the snippets file to adapt to the Materials system. 
  - We redrew the robots with diffuse shading 
  - We redrew the arcball with wireframe
  - We redrew the ground with the texture

  == Assignment 6.5 Part 2 ==
  - We implemented a bump mapping to improve the quality of our graphics. 
  - We stored the lights as Transform nodes to let us move them 
  - We changed the normal fragment shader to read in the bump information and use them as the normals 
 
Usage:
  - Use 'h' to get help directions. 
  - All keyboard controls and mouse actions are as per spec guidelines. Also described in 'h'.
  
