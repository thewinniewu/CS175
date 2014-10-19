////////////////////////////////////////////////////////////////////////
//
//   Harvard University
//   CS175 : Computer Graphics
//   Professor Steven Gortler
//
////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <memory>
#include <stdexcept>
#if __GNUG__
#   include <tr1/memory>
#endif

#ifdef __MAC__
#   include <OpenGL/gl3.h>
#   include <GLUT/glut.h>
#else
#   include <GL/glew.h>
#   include <GL/glut.h>
#endif

#include "arcball.h"
#include "cvec.h"
#include "matrix4.h"
#include "geometrymaker.h"
#include "ppm.h"
#include "glsupport.h"
#include "rigtform.h"

// from asst4-snippets
#include "asstcommon.h"
#include "scenegraph.h"
#include "drawer.h"
#include "picker.h"

#define CUBES 2
#define SKYCAMERA 0
#define REDCUBE 1
#define BLUECUBE 2

using namespace std;      // for string, vector, iostream, and other standard C++ stuff
using namespace tr1; // for shared_ptr

// G L O B A L S ///////////////////////////////////////////////////

// --------- IMPORTANT --------------------------------------------------------
// Before you start working on this assignment, set the following variable
// properly to indicate whether you want to use OpenGL 2.x with GLSL 1.0 or
// OpenGL 3.x+ with GLSL 1.5.
//
// Set g_Gl2Compatible = true to use GLSL 1.0 and g_Gl2Compatible = false to
// use GLSL 1.5. Use GLSL 1.5 unless your system does not support it.
//
// If g_Gl2Compatible=true, shaders with -gl2 suffix will be loaded.
// If g_Gl2Compatible=false, shaders with -gl3 suffix will be loaded.
// To complete the assignment you only need to edit the shader files that get
// loaded
// ----------------------------------------------------------------------------

const bool g_Gl2Compatible = false;

static const float g_frustMinFov = 60.0;  // A minimal of 60 degree field of view
static float g_frustFovY = g_frustMinFov; // FOV in y direction (updated by updateFrustFovY)

static const float g_frustNear = -0.1;    // near plane
static const float g_frustFar = -50.0;    // far plane
static const float g_groundY = -2.0;      // y coordinate of the ground
static const float g_groundSize = 10.0;   // half the ground length

static int g_windowWidth = 512;
static int g_windowHeight = 512;
static bool g_mouseClickDown = false;    // is the mouse button pressed
static bool g_mouseLClickButton, g_mouseRClickButton, g_mouseMClickButton;
static bool g_spaceDown = false;         // space state, for middle mouse emulation
static int g_mouseClickX, g_mouseClickY; // coordinates for mouse click event
static int g_activeShader = 0;

// from asst4-snippets
static const int PICKING_SHADER = 2; // index of the picking shader is g_shaderFiles
static const int g_numShaders = 3; // 3shaders instead of 2
static const char * const g_shaderFiles[g_numShaders][2] = {
  {"./shaders/basic-gl3.vshader", "./shaders/diffuse-gl3.fshader"},
  {"./shaders/basic-gl3.vshader", "./shaders/solid-gl3.fshader"},
  {"./shaders/basic-gl3.vshader", "./shaders/pick-gl3.fshader"},
};
static const char * const g_shaderFilesGl2[g_numShaders][2] = {
  {"./shaders/basic-gl2.vshader", "./shaders/diffuse-gl2.fshader"},
  {"./shaders/basic-gl2.vshader", "./shaders/solid-gl2.fshader"},
  {"./shaders/basic-gl2.vshader", "./shaders/pick-gl2.fshader"},
};
static vector<shared_ptr<ShaderState> > g_shaderStates; // our global shader states




// --------- Geometry

// Macro used to obtain relative offset of a field within a struct
#define FIELD_OFFSET(StructType, field) &(((StructType *)0)->field)

// A vertex with floating point position and normal
struct VertexPN {
  Cvec3f p, n;

  VertexPN() {}
  VertexPN(float x, float y, float z,
           float nx, float ny, float nz)
    : p(x,y,z), n(nx, ny, nz)
  {}

  // Define copy constructor and assignment operator from GenericVertex so we can
  // use make* functions from geometrymaker.h
  VertexPN(const GenericVertex& v) {
    *this = v;
  }

  VertexPN& operator = (const GenericVertex& v) {
    p = v.pos;
    n = v.normal;
    return *this;
  }
};

struct Geometry {
  GlBufferObject vbo, ibo;
  GlArrayObject vao;
  int vboLen, iboLen;

  Geometry(VertexPN *vtx, unsigned short *idx, int vboLen, int iboLen) {
    this->vboLen = vboLen;
    this->iboLen = iboLen;
	
    // Now create the VBO and IBO
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(VertexPN) * vboLen, vtx, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned short) * iboLen, idx, GL_STATIC_DRAW);
  }

  void draw(const ShaderState& curSS) {
	// bind the object's VAO
	glBindVertexArray(vao);
	
    // Enable the attributes used by our shader
    safe_glEnableVertexAttribArray(curSS.h_aPosition);
    safe_glEnableVertexAttribArray(curSS.h_aNormal);

    // bind vbo
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    safe_glVertexAttribPointer(curSS.h_aPosition, 3, GL_FLOAT, GL_FALSE, sizeof(VertexPN), FIELD_OFFSET(VertexPN, p));
    safe_glVertexAttribPointer(curSS.h_aNormal, 3, GL_FLOAT, GL_FALSE, sizeof(VertexPN), FIELD_OFFSET(VertexPN, n));

    // bind ibo
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);

    // draw!
    glDrawElements(GL_TRIANGLES, iboLen, GL_UNSIGNED_SHORT, 0);

    // Disable the attributes used by our shader
    safe_glDisableVertexAttribArray(curSS.h_aPosition);
    safe_glDisableVertexAttribArray(curSS.h_aNormal);
	
	// disable VAO
	glBindVertexArray(NULL);
  }
};

typedef SgGeometryShapeNode<Geometry> MyShapeNode;

// Vertex buffer and index buffer associated with the ground and cube geometry
static shared_ptr<Geometry> g_ground, g_cube, g_sphere;

// --------- Scene

// from asst4-snippets
static shared_ptr<SgRootNode> g_world;
static shared_ptr<SgRbtNode> g_skyNode, g_groundNode, g_robot1Node, g_robot2Node;
static shared_ptr<SgRbtNode> g_currentPickedRbtNode; // used later when you do picking



static const Cvec3 g_light1(2.0, 3.0, 14.0), g_light2(-2, -3.0, -5.0);  // define two lights positions in world space
static Cvec3f g_objectColors[CUBES] = {Cvec3f(1, 0, 0), Cvec3f(0, 0, 1)};

static RigTForm g_skyRbt = RigTForm(Cvec3(0.0, 0.25, 4.0));
static RigTForm g_worldSkyRbt = RigTForm(linFact(g_skyRbt));

/* Keeps track of current object we are using
 * as the camera. 'v' changes the view.
 * 0 = sky camera
 * 1 = red cube
 * 2 = blue cube
 */
static int g_currentCamera = SKYCAMERA;
static RigTForm g_objectRbt[CUBES] = {RigTForm(Cvec3(-1,0,0)), RigTForm(Cvec3(1,0,0))};

/* Keeps track of current object we are manipulating.
 * 'o' changes the current object.
 * 0 = sky camera
 * 1 = red cube
 * 2 = blue cube
 * */
static int g_currentObj = REDCUBE; 
static RigTForm g_auxFrame = transFact(g_objectRbt[0]) * linFact(g_skyRbt);

static const Cvec3 g_arcballColor = Cvec3(0,0,0);
static double g_arcballScale = 1.0;
static double g_arcballScreenRadius = 1.0;
static RigTForm g_arcballOrigin = g_auxFrame;

static bool g_isWorldSky = false;

///////////////// END OF G L O B A L S //////////////////////////////////////////////////

// some helpful functions
static RigTForm getCurrentView() {
  return (g_currentCamera == SKYCAMERA) ? 
          g_skyRbt 
          : g_objectRbt[g_currentCamera - 1];
}

static bool selfCubeManip() {
  return (g_currentObj > 0 && g_currentCamera == g_currentObj);
}

static bool useArcball() {
  return ((g_currentCamera == SKYCAMERA && g_currentObj > 0) 
          || (g_currentCamera == SKYCAMERA && g_isWorldSky)
          || (g_currentCamera != SKYCAMERA && g_currentCamera != g_currentObj)); 
  
}

// end helpful functions

static void initGround() {
  // A x-z plane at y = g_groundY of dimension [-g_groundSize, g_groundSize]^2
  VertexPN vtx[4] = {
    VertexPN(-g_groundSize, g_groundY, -g_groundSize, 0, 1, 0),
    VertexPN(-g_groundSize, g_groundY,  g_groundSize, 0, 1, 0),
    VertexPN( g_groundSize, g_groundY,  g_groundSize, 0, 1, 0),
    VertexPN( g_groundSize, g_groundY, -g_groundSize, 0, 1, 0),
  };
  unsigned short idx[] = {0, 1, 2, 0, 2, 3};
  g_ground.reset(new Geometry(&vtx[0], &idx[0], 4, 6));
}

static void initCubes() {
  int ibLen, vbLen;
  getCubeVbIbLen(vbLen, ibLen);

  // Temporary storage for cube geometry
  vector<VertexPN> vtx(vbLen);
  vector<unsigned short> idx(ibLen);
  

  makeCube(1, vtx.begin(), idx.begin());
  g_cube.reset(new Geometry(&vtx[0], &idx[0], vbLen, ibLen));
}

static void initSphere() {
  int ibLen, vbLen, slices, stacks;
  //signature: makeSphere(float radius, int slices, int stacks, VtxOutIter vtxIter, IdxOutIter idxIter) { 
  slices = 10;
  stacks = 10; 
  
  getSphereVbIbLen(slices, stacks, vbLen, ibLen);

  vector<VertexPN> vtx(vbLen);
  vector<unsigned short> idx(ibLen);

  makeSphere(1, slices, stacks, vtx.begin(), idx.begin()); 

  g_sphere.reset(new Geometry(&vtx[0], &idx[0], vbLen, ibLen)); 


}

// takes a projection matrix and send to the the shaders
static void sendProjectionMatrix(const ShaderState& curSS, const Matrix4& projMatrix) {
  GLfloat glmatrix[16];
  projMatrix.writeToColumnMajorMatrix(glmatrix); // send projection matrix
  safe_glUniformMatrix4fv(curSS.h_uProjMatrix, glmatrix);
}

// update g_frustFovY from g_frustMinFov, g_windowWidth, and g_windowHeight
static void updateFrustFovY() {
  if (g_windowWidth >= g_windowHeight)
    g_frustFovY = g_frustMinFov;
  else {
    const double RAD_PER_DEG = 0.5 * CS175_PI/180;
    g_frustFovY = atan2(sin(g_frustMinFov * RAD_PER_DEG) * g_windowHeight / g_windowWidth, cos(g_frustMinFov * RAD_PER_DEG)) / RAD_PER_DEG;
  }
}

static Matrix4 makeProjectionMatrix() {
  return Matrix4::makeProjection(
           g_frustFovY, g_windowWidth / static_cast <double> (g_windowHeight),
           g_frustNear, g_frustFar);
}

static void setWrtFrame() {
  if (g_currentObj == SKYCAMERA) {
    if (g_currentCamera == SKYCAMERA) {
      if (g_isWorldSky) {
        g_auxFrame = linFact(g_skyRbt);
      } else {
        g_auxFrame = g_skyRbt;
      }
    }
  } else {
    if (g_currentCamera == SKYCAMERA) {
     g_auxFrame = transFact(g_objectRbt[g_currentObj - 1]) * linFact(g_skyRbt);
    } else {
     g_auxFrame = transFact(g_objectRbt[g_currentObj - 1]) * linFact(getCurrentView());
    }
  }
}

static void drawStuff(const ShaderState& curSS, bool picking) {
  setWrtFrame();       
  
  // short hand for current shader state
  //const ShaderState& curSS = *g_shaderStates[g_activeShader];

  // build & send proj. matrix to vshader
  const Matrix4 projmat = makeProjectionMatrix();
  sendProjectionMatrix(curSS, projmat);

  // use the skyRbt as the eyeRbt
  const RigTForm eyeRbt = getCurrentView();
  const RigTForm invEyeRbt = inv(eyeRbt);

  const Cvec3 eyeLight1 = Cvec3(invEyeRbt * Cvec4(g_light1, 1)); // g_light1 position in eye coordinates
  const Cvec3 eyeLight2 = Cvec3(invEyeRbt * Cvec4(g_light2, 1)); // g_light2 position in eye coordinates
  safe_glUniform3f(curSS.h_uLight, eyeLight1[0], eyeLight1[1], eyeLight1[2]);
  safe_glUniform3f(curSS.h_uLight2, eyeLight2[0], eyeLight2[1], eyeLight2[2]);

  // draw ground
  // ===========
  //
  const Matrix4 groundRbt = Matrix4();  // identity
  Matrix4 MVM = rigTFormToMatrix(invEyeRbt) * groundRbt;
  Matrix4 NMVM = normalMatrix(MVM);
  sendModelViewNormalMatrix(curSS, MVM, NMVM);
  safe_glUniform3f(curSS.h_uColor, 0.1, 0.95, 0.1); // set color
  g_ground->draw(curSS);

  // draw cubes
  // ==========
  if (!picking) {
    Drawer drawer(invEyeRbt, curSS);
    g_world->accept(drawer);
 
    // draw sphere
    // ===============
           
    if (g_currentCamera == SKYCAMERA) {
      if (g_isWorldSky) { 
        g_arcballOrigin = inv(RigTForm());
      } else {
        g_arcballOrigin = g_objectRbt[g_currentObj - 1];
      }
    } else {
      if (g_currentObj != SKYCAMERA) {
        g_arcballOrigin = g_objectRbt[g_currentObj - 1];
      } else {
        g_arcballOrigin = g_objectRbt[g_currentCamera - 1];
      }
    }
         
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    // compute MVM & NMVM
    Matrix4 scale = Matrix4::makeScale(g_arcballScale * g_arcballScreenRadius * 0.008);
    MVM = rigTFormToMatrix(invEyeRbt * g_arcballOrigin) * scale;
    NMVM = normalMatrix(MVM);
    // send in MVM and NMVM
    sendModelViewNormalMatrix(curSS, MVM, NMVM);
    // send uColor
    safe_glUniform3f(curSS.h_uColor, g_arcballColor[0], g_arcballColor[1], g_arcballColor[3]); 
    //draw 
    g_sphere->draw(curSS);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); 
  } else {
    Picker picker(invEyeRbt, curSS);
    g_world->accept(picker);
    glFlush();
    g_currentPickedRbtNode = picker.getRbtNodeAtXY(g_mouseClickX,g_mouseClickY);
    if (g_currentPickedRbtNode == g_groundNode)
      g_currentPickedRbtNode = shared_ptr<SgRbtNode>(); //set to NULL
  }
}

static void display() {
  glUseProgram(g_shaderStates[g_activeShader]->program);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);                   // clear framebuffer color&depth

  drawStuff(*g_shaderStates[g_activeShader], false);

  glutSwapBuffers();                                    // show the back buffer (where we rendered stuff)

  checkGlErrors();
}

static void reshape(const int w, const int h) {
  g_windowWidth = w;
  g_windowHeight = h;
  g_arcballScreenRadius = 0.25 * min(g_windowWidth, g_windowHeight);  
  glViewport(0, 0, w, h);
  cerr << "Size of window is now " << w << "x" << h << endl;
  updateFrustFovY();
  glutPostRedisplay();
}

static RigTForm getCurrentObj() {
  if (g_currentObj == 0) {
      return g_skyRbt;
  } else {
     return g_objectRbt[g_currentObj - 1];
  }
}

static double getZCoord(double x, double y, double pixRad) {
  return sqrt(max(0.0, (pow(pixRad,2) - pow(y , 2) - pow(x ,2))));
}


static RigTForm getArcballTransform(const int x, const int y) {
  RigTForm obj = getCurrentObj(); 
 
  Cvec2 arcballScreenPos;
  
  if (g_isWorldSky) {
    arcballScreenPos = Cvec2((g_windowWidth - 1)/2.0, (g_windowHeight - 1)/2.0);
  } else {
    arcballScreenPos = getScreenSpaceCoord(
      (inv(getCurrentView()) * obj).getTranslation(), 
      makeProjectionMatrix(),
      g_frustNear,
      g_frustFovY, 
      g_windowWidth,
      g_windowHeight
    );
  }
 
  Cvec3 arcballCenter = Cvec3(arcballScreenPos, 0);

  Cvec3 before = Cvec3(g_mouseClickX, g_mouseClickY, 0) - arcballCenter;
  Cvec3 after = Cvec3(x, y, 0) - arcballCenter;

  Cvec3 v1 = normalize(Cvec3(
    before[0], 
    before[1], 
    getZCoord(before[0], before[1], g_arcballScreenRadius)
  ));
  Cvec3 v2 = normalize(Cvec3(
    after[0],
    after[1],
    getZCoord(after[0], after[1], g_arcballScreenRadius)
  ));

  if (g_isWorldSky) {
    return RigTForm(Quat(0, (v1 * -1.0)) * Quat(0, v2));
  } else {
    return RigTForm(Quat(0, v2) * Quat(0, (v1 * -1.0)));
  }

}

static void motion(const int x, const int y) {

  if (g_currentCamera != SKYCAMERA && g_currentObj == SKYCAMERA) {
    return;
  }

  double dx = x - g_mouseClickX;
  double dy = g_windowHeight - y - 1 - g_mouseClickY;
  
  RigTForm m = RigTForm();

  // self-movement flag.
  // 1 = red cube moving itself
  // 2 = blue cube moving itself
  // 3 = skycamera moving itself 
  // 0 = everything else
   

  int flag = 0;
  if (g_currentObj == g_currentCamera) {
    if (g_currentObj == REDCUBE) {
      flag = 1;
    } else if (g_currentObj == BLUECUBE) { 
      flag = 2;
    } else if (g_currentObj == SKYCAMERA) {
      flag = 3;
    }
  }
  double tFactor = 0.01; // translation factor 
  
  if (useArcball()) {
    tFactor = g_arcballScale * 0.008; 
  }
  
  if (g_mouseLClickButton && !g_mouseRClickButton && !g_spaceDown) { // left button down?
    if (useArcball()) {
      m = getArcballTransform(x, g_windowHeight - 1 - y);
    } else if (g_currentObj == 0 || (flag == 1) || (flag == 2)) { 
      m.setRotation(m.getRotation().makeXRotation(dy) * m.getRotation().makeYRotation(-dx));
    } else { 
      m.setRotation(m.getRotation().makeXRotation(-dy) * m.getRotation().makeYRotation(dx));
    }
  }
  else if (g_mouseRClickButton && !g_mouseLClickButton) { // right button down?
    if (g_currentObj == 0 && mequals(rigTFormToMatrix(g_auxFrame), rigTFormToMatrix(g_worldSkyRbt))) {
      m.setTranslation(m.getTranslation() + Cvec3(-dx,-dy,0) * tFactor);
    } else { 
      m.setTranslation(m.getTranslation() + Cvec3(dx, dy, 0) * tFactor);
    } 
  }
  
  else if (g_mouseMClickButton || (g_mouseLClickButton && g_mouseRClickButton) || (g_mouseLClickButton && !g_mouseRClickButton && g_spaceDown)) {  // middle or (left and right, or left + space) button down?
    if (g_currentObj == 0 && mequals(rigTFormToMatrix(g_auxFrame), rigTFormToMatrix(g_worldSkyRbt))) {
      m.setTranslation(m.getTranslation() + Cvec3(0,0,dy) * tFactor);
    } else { 
      m.setTranslation(m.getTranslation() + Cvec3(0, 0, -dy) * tFactor);
    }
    
    g_arcballScale = getScreenToEyeScale(
      (inv(getCurrentView()) * g_arcballOrigin).getTranslation()[2], 
      g_frustFovY,
      g_windowHeight
    ) * 108;
  }
  
  if (g_mouseClickDown) {
    m = g_auxFrame * m * inv(g_auxFrame);

    if (g_currentObj == 0) {
      g_skyRbt = m * g_skyRbt;
    } else {
      g_objectRbt[g_currentObj - 1] = m * g_objectRbt[g_currentObj - 1];
    }

  }
          
          
     if (flag == 3) { // skycamera movement
        g_currentCamera = SKYCAMERA;
      } 
    
    if (flag == 1) { 
      g_currentCamera = REDCUBE;
    } else if (flag == 2) {
      g_currentCamera = BLUECUBE;
    }
  
  glutPostRedisplay(); // we always redraw if we changed the scene
 
  g_mouseClickX = x;
  g_mouseClickY = g_windowHeight - y - 1;
}


static void mouse(const int button, const int state, const int x, const int y) {
  g_mouseClickX = x;
  g_mouseClickY = g_windowHeight - y - 1;  // conversion from GLUT window-coordinate-system to OpenGL window-coordinate-system

  g_mouseLClickButton |= (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN);
  g_mouseRClickButton |= (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN);
  g_mouseMClickButton |= (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN);

  g_mouseLClickButton &= !(button == GLUT_LEFT_BUTTON && state == GLUT_UP);
  g_mouseRClickButton &= !(button == GLUT_RIGHT_BUTTON && state == GLUT_UP);
  g_mouseMClickButton &= !(button == GLUT_MIDDLE_BUTTON && state == GLUT_UP);

  g_mouseClickDown = g_mouseLClickButton || g_mouseRClickButton || g_mouseMClickButton;

  glutPostRedisplay();
}

static void keyboardUp(const unsigned char key, const int x, const int y) {
  switch (key) {
  case ' ':
    g_spaceDown = false;
    break;
  }
  glutPostRedisplay();
}

static void changeView() {
  g_isWorldSky = false;
  cout << "The current view is:";
  if (g_currentCamera == SKYCAMERA) {
    g_currentCamera = REDCUBE;
    cout << "g_objectRbt[0]. You are the red cube.\n"; 
  }
  else if (g_currentCamera == REDCUBE) {
    g_currentCamera = BLUECUBE;
    cout << "g_objectRbt[1]. You are the blue cube.\n"; 
  } else {
    g_currentCamera = SKYCAMERA;
    cout << "g_skyRbt. You are omniscient.\n"; 
  }
  glutPostRedisplay();
}

static void changeCurrentObject() {
  g_currentObj = (g_currentObj + 1) % (CUBES + 1);
  setWrtFrame(); 
}

static void changeSkyCameraAux() {
	if (g_isWorldSky) { 
    g_isWorldSky = false;
  } else {
    g_isWorldSky = true; 
  }
  setWrtFrame();
}

static void keyboard(const unsigned char key, const int x, const int y) {
  switch (key) {
  case 27:
    exit(0);                                  // ESC
  case 'v':
    changeView(); 
    break;
  case 'o':
    changeCurrentObject();
    break;
  case 'm':
    if (g_currentObj == 0) {
      changeSkyCameraAux();
    } 
    break;
  case 'h':
    cout << " ============== H E L P ==============\n\n"
    << "h\t\thelp menu\n"
    << "s\t\tsave screenshot\n"
    << "f\t\tToggle flat shading on/off.\n"
    << "o\t\tCycle object to edit\n"
    << "v\t\tCycle view\n"
    << "drag left mouse to rotate\n" 
    << "drag right mouse to translate left-right or up-down\n"
    << "drag middle mouse to translate in-out\n" << endl;
    break;
  case 's':
    glFlush();
    writePpmScreenshot(g_windowWidth, g_windowHeight, "out.ppm");
    break;
  case 'f':
    g_activeShader ^= 1;
    break;
  case ' ':
    g_spaceDown = true;
    break;
  }
  glutPostRedisplay();
}

static void initGlutState(int argc, char * argv[]) {
  glutInit(&argc, argv);                                  // initialize Glut based on cmd-line args
#ifdef __MAC__
  glutInitDisplayMode(GLUT_3_2_CORE_PROFILE|GLUT_RGBA|GLUT_DOUBLE|GLUT_DEPTH); // core profile flag is required for GL 3.2 on Mac
#else
  glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE|GLUT_DEPTH);  //  RGBA pixel channels and double buffering
#endif
  glutInitWindowSize(g_windowWidth, g_windowHeight);      // create a window
  glutCreateWindow("Assignment 2");                       // title the window

  glutIgnoreKeyRepeat(true);                              // avoids repeated keyboard calls when holding space to emulate middle mouse

  glutDisplayFunc(display);                               // display rendering callback
  glutReshapeFunc(reshape);                               // window reshape callback
  glutMotionFunc(motion);                                 // mouse movement callback
  glutMouseFunc(mouse);                                   // mouse click callback
  glutKeyboardFunc(keyboard);
  glutKeyboardUpFunc(keyboardUp);
}

static void initGLState() {
  glClearColor(128./255., 200./255., 255./255., 0.);
  glClearDepth(0.);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glCullFace(GL_BACK);
  glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_GREATER);
  glReadBuffer(GL_BACK);
  if (!g_Gl2Compatible)
    glEnable(GL_FRAMEBUFFER_SRGB);
}

static void initShaders() {
  g_shaderStates.resize(g_numShaders);
  for (int i = 0; i < g_numShaders; ++i) {
    if (g_Gl2Compatible)
      g_shaderStates[i].reset(new ShaderState(g_shaderFilesGl2[i][0], g_shaderFilesGl2[i][1]));
    else
      g_shaderStates[i].reset(new ShaderState(g_shaderFiles[i][0], g_shaderFiles[i][1]));
  }
}

static void initGeometry() {
  initGround();
  initCubes();
  initSphere();
}

static void constructRobot(shared_ptr<SgTransformNode> base, const Cvec3& color) {
  const double ARM_LEN = 0.7,
               ARM_THICK = 0.25,
               TORSO_LEN = 1.5,
               TORSO_THICK = 0.25,
               TORSO_WIDTH = 1;
  const int NUM_JOINTS = 3,
            NUM_SHAPES = 3;

  struct JointDesc { 
    int parent;
    float x, y, z;
  };

  JointDesc jointDesc[NUM_JOINTS] = {
    {-1}, //torso
    {0, TORSO_WIDTH/2, TORSO_LEN/2, 0}, //upper right arm
    {1, ARM_LEN, 0 ,0}, //lower right arm
  };

  struct ShapeDesc {
    int parentJointId;
    float x, y, z, sx, sy, sz;
    shared_ptr<Geometry> geometry;
  };
  
  ShapeDesc shapeDesc[NUM_SHAPES] = {
    {0, 0,         0, 0, TORSO_WIDTH, TORSO_LEN, TORSO_THICK, g_cube}, //torso
    {1, ARM_LEN/2, 0, 0, ARM_LEN, ARM_THICK, ARM_THICK, g_cube}, // upper right arm
    {2, ARM_LEN/2, 0, 0, ARM_LEN, ARM_THICK, ARM_THICK, g_cube}, //lower right arm
  };

  shared_ptr<SgTransformNode> jointNodes[NUM_JOINTS];

  for (int i = 0; i < NUM_JOINTS; ++i) {
    if (jointDesc[i].parent == -1)
      jointNodes[i] = base;
    else {
      jointNodes[i].reset(new SgRbtNode(RigTForm(Cvec3(jointDesc[i].x, jointDesc[i].y, jointDesc[i].z))));
      jointNodes[jointDesc[i].parent]->addChild(jointNodes[i]);
    }
  }
  for (int i = 0; i < NUM_SHAPES; ++i) {
    shared_ptr<MyShapeNode> shape(
      new MyShapeNode(shapeDesc[i].geometry,
                      color,
                      Cvec3(shapeDesc[i].x, shapeDesc[i].y, shapeDesc[i].z),
                      Cvec3(0, 0, 0),
                      Cvec3(shapeDesc[i].sx, shapeDesc[i].sy, shapeDesc[i].sz)
       )
     );
    jointNodes[shapeDesc[i].parentJointId]->addChild(shape);
  }
}

static void initScene() {
  g_world.reset(new SgRootNode());
  
  g_skyNode.reset(new SgRbtNode(RigTForm(Cvec3(0.0, 0.25, 4.0))));
  
  g_groundNode.reset(new SgRbtNode());
  g_groundNode->addChild(shared_ptr<MyShapeNode>(
                          new MyShapeNode(g_ground, Cvec3(0.1, 0.95, 0.1))));

  g_robot1Node.reset(new SgRbtNode(RigTForm(Cvec3(-2, 1, 0))));
  g_robot2Node.reset(new SgRbtNode(RigTForm(Cvec3(2, 1, 0))));

  constructRobot(g_robot1Node, Cvec3(1, 0, 0)); // a red robot
  constructRobot(g_robot2Node, Cvec3(0, 0, 1)); //a blue robot

  g_world->addChild(g_skyNode);
  g_world->addChild(g_groundNode);
  g_world->addChild(g_robot1Node);
  g_world->addChild(g_robot2Node);

}


int main(int argc, char * argv[]) {
  try {
    initGlutState(argc,argv);

	// on Mac, we shouldn't use GLEW.

#ifndef __MAC__
    glewInit(); // load the OpenGL extensions
#endif

    cout << (g_Gl2Compatible ? "Will use OpenGL 2.x / GLSL 1.0" : "Will use OpenGL 3.x / GLSL 1.5") << endl;

#ifndef __MAC__
    if ((!g_Gl2Compatible) && !GLEW_VERSION_3_0)
      throw runtime_error("Error: card/driver does not support OpenGL Shading Language v1.3");
    else if (g_Gl2Compatible && !GLEW_VERSION_2_0)
      throw runtime_error("Error: card/driver does not support OpenGL Shading Language v1.0");
#endif
	
    initGLState();
    initShaders();
    initGeometry();
    initScene();

    glutMainLoop();
    return 0;
  }
  catch (const runtime_error& e) {
    cout << "Exception caught: " << e.what() << endl;
    return -1;
  }
}
