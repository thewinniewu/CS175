////////////////////////////////////////////////////////////////////////
//
//   Harvard University
//   CS175 : Computer Graphics
//   Professor Steven Gortler
//
////////////////////////////////////////////////////////////////////////

#include <cstddef>
#include <vector>
#include <string>
#include <memory>
#include <stdexcept>
#include <algorithm>
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

#include "ppm.h"
#include "cvec.h"
#include "matrix4.h"
#include "rigtform.h"
#include "glsupport.h"
#include "geometrymaker.h"
#include "arcball.h"

// from asst4-snippets
#include "asstcommon.h"
#include "scenegraph.h"
#include "drawer.h"
#include "picker.h"

// for keyframe animation
#include <list>
#include "sgutils.h"

using namespace std;
using namespace tr1;

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

enum ObjId {SKY=0, OBJECT0=1, OBJECT1=2, LIMB = 3};
enum SkyMode {WORLD_SKY=0, SKY_SKY=1};

static const char * const g_objNames[] = {"Sky", "Object 0", "Object 1"};

static int g_windowWidth = 512;
static int g_windowHeight = 512;
static bool g_mouseClickDown = false;    // is the mouse button pressed
static bool g_mouseLClickButton, g_mouseRClickButton, g_mouseMClickButton;
static bool g_spaceDown = false;         // space state, for middle mouse emulation
static int g_mouseClickX, g_mouseClickY; // coordinates for mouse click event
static int g_activeShader = 0;

static ObjId g_activeEye = SKY;
static SkyMode g_activeCameraFrame = WORLD_SKY;

static bool g_displayArcball = true;
static double g_arcballScreenRadius = 100; // number of pixels
static double g_arcballScale = 1;

// asst-snippets
static const int PICKING_SHADER = 2; // index of the picking shader is g_shaerFiles
static const int g_numShaders = 3; // 3 shaders instead of 2
static const char * const g_shaderFiles[g_numShaders][2] = {
		{ "./shaders/basic-gl3.vshader", "./shaders/diffuse-gl3.fshader" },
		{ "./shaders/basic-gl3.vshader", "./shaders/solid-gl3.fshader" },
		{ "./shaders/basic-gl3.vshader", "./shaders/pick-gl3.fshader" }
};
static const char * const g_shaderFilesGl2[g_numShaders][2] = {
		{ "./shaders/basic-gl2.vshader", "./shaders/diffuse-gl2.fshader" },
		{ "./shaders/basic-gl2.vshader", "./shaders/solid-gl2.fshader" },
		{ "./shaders/basic-gl2.vshader", "./shaders/pick-gl2.fshader" }
};

static vector<shared_ptr<ShaderState> > g_shaderStates; // our global shader states

// --------- Geometry

// Macro used to obtain relative offset of a field within a struct
#define FIELD_OFFSET(StructType, field) ((GLvoid*)offsetof(StructType, field))

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
static shared_ptr<SgRbtNode> g_currentPickedRbtNode; // used later when you do pickingj

// for keyframe animation
typedef std::vector<std::tr1::shared_ptr<SgRbtNode> > SgRbtNodes;
static shared_ptr<SgRbtNodes> g_currentKeyframe; // pointer to vector of SgRbtNodes that represent the current frame 
list<SgRbtNodes> keyframeList;  // list of SgRbtNodes 
// list<SgRbtNodes>::iterator iter = keyframeList.begin(); // iterator through the list
 
static const Cvec3 g_light1(2.0, 3.0, 14.0), g_light2(-2, -3.0, -5.0);  // define two lights positions in world space
static RigTForm g_skyRbt = RigTForm(Cvec3(0.0, 0.25, 4.0));
static RigTForm g_objectRbt[2] = {RigTForm(Cvec3(-1,0,0)), RigTForm(Cvec3(1,0,0))};
static Cvec3f g_objectColors[2] = {Cvec3f(1, 0, 0), Cvec3f(0, 0, 1)};

static bool g_pickerMode = false;

///////////////// END OF G L O B A L S //////////////////////////////////////////////////


static shared_ptr<SgRbtNode> getNodeForEye(ObjId i) {
   switch (i) { 
    case SKY:
      return g_skyNode;
    case OBJECT0:
      return g_robot1Node;
    case OBJECT1:
      return g_robot2Node;
    case LIMB:
      cout << "Shouldn't have this eye ever, something has gone wrong";
      exit(1); 
   }
}

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
  int ibLen, vbLen;
  getSphereVbIbLen(20, 10, vbLen, ibLen);

  // Temporary storage for sphere geometry
  vector<VertexPN> vtx(vbLen);
  vector<unsigned short> idx(ibLen);
  makeSphere(1, 20, 10, vtx.begin(), idx.begin());
  g_sphere.reset(new Geometry(&vtx[0], &idx[0], vtx.size(), idx.size()));
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

enum ManipMode {
  ARCBALL_ON_PICKED,
  ARCBALL_ON_SKY,
  EGO_MOTION
};

static ObjId getCurrentObjId() {
  if (g_currentPickedRbtNode == g_robot1Node)
    return OBJECT0;
  else if (g_currentPickedRbtNode == g_robot2Node)
    return OBJECT1;
  else if (g_currentPickedRbtNode)
    return LIMB;
  else
    return SKY;
}

static ManipMode getManipMode() {
  if (g_currentPickedRbtNode && (g_activeEye != getCurrentObjId())) {
    return ARCBALL_ON_PICKED;
  } else if (g_activeEye == SKY && g_activeCameraFrame == WORLD_SKY)
    return ARCBALL_ON_SKY;
  else {
    return EGO_MOTION;
  }
}

static bool shouldUseArcball() {
  return getManipMode() != EGO_MOTION && (!(g_activeEye != SKY && getCurrentObjId() == SKY));
}

// The translation part of the aux frame either comes from the current
// active object, or is the identity matrix when
static RigTForm getArcballRbt() {
  if (g_currentPickedRbtNode) 
     return getPathAccumRbt(g_world, g_currentPickedRbtNode);
  else
     return RigTForm();
}

static void updateArcballScale() {
  RigTForm arcballEye = inv(getNodeForEye(g_activeEye)->getRbt()) * getArcballRbt();
  double depth = arcballEye.getTranslation()[2];
  if (depth > -CS175_EPS)
    g_arcballScale = 0.02;
  else
    g_arcballScale = getScreenToEyeScale(depth, g_frustFovY, g_windowHeight);
}

static void drawArcBall(const ShaderState& curSS) {
  // switch to wire frame mode
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  RigTForm arcballEye = inv(getNodeForEye(g_activeEye)->getRbt()) * getArcballRbt();
  Matrix4 MVM = rigTFormToMatrix(arcballEye) * Matrix4::makeScale(Cvec3(1,1,1)* g_arcballScale * g_arcballScreenRadius);
  sendModelViewNormalMatrix(curSS, MVM, normalMatrix(MVM));

  safe_glUniform3f(curSS.h_uColor, 0.27, 0.82, 0.35); // set color

  g_sphere->draw(curSS);

  // switch back to solid mode
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

static void drawStuff(const ShaderState& curSS, bool picking) {
  // if we are not translating, update arcball scale
  if (!(g_mouseMClickButton || (g_mouseLClickButton && g_mouseRClickButton) || (g_mouseLClickButton && !g_mouseRClickButton && g_spaceDown)))
    updateArcballScale();

  // build & send proj. matrix to vshader
  const Matrix4 projmat = makeProjectionMatrix();
  sendProjectionMatrix(curSS, projmat);

  const RigTForm eyeRbt = getPathAccumRbt(g_world, getNodeForEye(g_activeEye));
  const RigTForm invEyeRbt = inv(eyeRbt);

  const Cvec3 eyeLight1 = Cvec3(invEyeRbt * Cvec4(g_light1, 1));
  const Cvec3 eyeLight2 = Cvec3(invEyeRbt * Cvec4(g_light2, 1));
  safe_glUniform3f(curSS.h_uLight, eyeLight1[0], eyeLight1[1], eyeLight1[2]);
  safe_glUniform3f(curSS.h_uLight2, eyeLight2[0], eyeLight2[1], eyeLight2[2]);

  if (!picking) {
	  Drawer drawer(invEyeRbt, curSS);
	  g_world->accept(drawer);

	  if (g_displayArcball && shouldUseArcball()) {
		  drawArcBall(curSS);
	  }
  } else {
	  Picker picker(invEyeRbt, curSS);
	  g_world->accept(picker);
	  glFlush();
	  g_currentPickedRbtNode = picker.getRbtNodeAtXY(g_mouseClickX, g_mouseClickY);
	  
	  if (g_currentPickedRbtNode == g_groundNode) {
		  g_currentPickedRbtNode = shared_ptr<SgRbtNode>();   // set to NULL
		  cout << "clicked on ground\n";
	  }
  }
}

static void display() {
  glUseProgram(g_shaderStates[g_activeShader]->program);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);                   // clear framebuffer color&depth

  drawStuff(*g_shaderStates[g_activeShader], g_pickerMode);

  if (!g_pickerMode) {
	  glutSwapBuffers();                                    // show the back buffer (where we rendered stuff)
  }

  checkGlErrors();
}

static void reshape(const int w, const int h) {
  g_windowWidth = w;
  g_windowHeight = h;
  glViewport(0, 0, w, h);
  cerr << "Size of window is now " << w << "x" << h << endl;
  g_arcballScreenRadius = max(1.0, min(h,w) * 0.25);
  updateFrustFovY();
  glutPostRedisplay();
}


static Cvec3 getArcballDirection(const Cvec2& p, const double r) {
  double n2 = norm2(p);
  if (n2 >= r*r)
    return normalize(Cvec3(p, 0));
  else
    return normalize(Cvec3(p, sqrt(r*r - n2)));
}


static RigTForm moveArcball(const Cvec2& p0, const Cvec2& p1) {
  const Matrix4 projMatrix = makeProjectionMatrix();
  const RigTForm eyeInverse = inv(getNodeForEye(g_activeEye)->getRbt());
  const Cvec3 arcballCenter = getArcballRbt().getTranslation();
  const Cvec3 arcballCenter_ec = Cvec3(eyeInverse * Cvec4(arcballCenter, 1));

  if (arcballCenter_ec[2] > -CS175_EPS)
    return RigTForm();

  Cvec2 ballScreenCenter = getScreenSpaceCoord(arcballCenter_ec,
                                               projMatrix, g_frustNear, g_frustFovY, g_windowWidth, g_windowHeight);
  const Cvec3 v0 = getArcballDirection(p0 - ballScreenCenter, g_arcballScreenRadius);
  const Cvec3 v1 = getArcballDirection(p1 - ballScreenCenter, g_arcballScreenRadius);

  return RigTForm(Quat(0.0, v1[0], v1[1], v1[2]) * Quat(0.0, -v0[0], -v0[1], -v0[2]));
}

static RigTForm doMtoOwrtA(const RigTForm& M, const RigTForm& O, const RigTForm& A) {
  return A * M * inv(A) * O;
}

static RigTForm getMRbt(const double dx, const double dy) {
  RigTForm M;

  if (g_mouseLClickButton && !g_mouseRClickButton && !g_spaceDown) {
    if (shouldUseArcball())
      M = moveArcball(Cvec2(g_mouseClickX, g_mouseClickY), Cvec2(g_mouseClickX + dx, g_mouseClickY + dy));
    else
      M = RigTForm(Quat::makeXRotation(-dy) * Quat::makeYRotation(dx));
  }
  else {
    double movementScale = getManipMode() == EGO_MOTION ? 0.02 : g_arcballScale;
    if (g_mouseRClickButton && !g_mouseLClickButton) {
       M = RigTForm(Cvec3(dx, dy, 0) * movementScale);
    }
    else if (g_mouseMClickButton || (g_mouseLClickButton && g_mouseRClickButton) || (g_mouseLClickButton && g_spaceDown)) {
      M = RigTForm(Cvec3(0, 0, -dy) * movementScale);
    }
  }

  switch (getManipMode()) {
  case ARCBALL_ON_PICKED:
    break;
  case ARCBALL_ON_SKY:
    M = inv(M);
    break;
  case EGO_MOTION:
    if (g_mouseLClickButton && !g_mouseRClickButton && !g_spaceDown) // only invert rotation
      M = inv(M);
    break;
  }
  return M;
}

static RigTForm makeMixedFrame(const RigTForm& objRbt, const RigTForm& eyeRbt) {
  return transFact(objRbt) * linFact(eyeRbt);
}

static void pick() {
	// We need to set the clear color to black, for pick rendering.
	// so let's save the clear color
	GLdouble clearColor[4];
	glGetDoublev(GL_COLOR_CLEAR_VALUE, clearColor);

	glClearColor(0, 0, 0, 0);

	// using PICKING_SHADER as the shader
	glUseProgram(g_shaderStates[PICKING_SHADER]->program);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	drawStuff(*g_shaderStates[PICKING_SHADER], true);

	// Uncomment below and comment out the glutPostRedisplay in mouse(...) call back
	// to see result of the pick rendering pass
  
  //glutSwapBuffers();

	//Now set back the clear color
	glClearColor(clearColor[0], clearColor[1], clearColor[2], clearColor[3]);

	checkGlErrors();
}

static void motion(const int x, const int y) {
  if (!g_mouseClickDown)
    return;
  if (getCurrentObjId() == SKY && g_activeEye != SKY)
    return;                  // we do not edit the eye when viewed from the objects
  const double dx = x - g_mouseClickX;
  const double dy = g_windowHeight - y - 1 - g_mouseClickY;

  const RigTForm M = getMRbt(dx, dy);   // the "action" matrix

  RigTForm A;
  // the matrix for the auxiliary frame (the w.r.t.)
  if (g_activeEye == SKY) {
    if (g_currentPickedRbtNode)
      A = transFact(g_currentPickedRbtNode->getRbt()) * linFact(g_skyNode->getRbt()); 
    else
      A = linFact(g_skyNode->getRbt());
  } else {
     RigTForm obj = getPathAccumRbt(g_world, g_currentPickedRbtNode, 0);
     RigTForm parent = getPathAccumRbt(g_world, g_currentPickedRbtNode, 1);
     RigTForm frame = getPathAccumRbt(g_world, getNodeForEye(g_activeEye), 0);   
     A = inv(parent) * transFact(obj) * linFact(frame);
  }
 
  if (g_currentPickedRbtNode) {
    RigTForm O = doMtoOwrtA(M, g_currentPickedRbtNode->getRbt(), A);
    g_currentPickedRbtNode->setRbt(O);   
  } else {
    shared_ptr<SgRbtNode> eyeNode = getNodeForEye(g_activeEye);
    RigTForm O = doMtoOwrtA(M, eyeNode->getRbt(), A); 
    eyeNode->setRbt(O); 
  }

  g_mouseClickX += dx;
  g_mouseClickY += dy;

  glutPostRedisplay();                    // we always redraw if we changed the scene
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

  if (g_pickerMode && g_mouseLClickButton && !g_mouseRClickButton) {
	  pick();
	  g_pickerMode = false;
  }
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

static void initializeNewKeyframe() {
	SgRbtNodes new_keyframe;
	if (keyframeList.empty()) {
		printf("The current frame has not been initialized, creating a new keyframe\n");
		// create a new keyframe from the beginning of the stack
		keyframeList.push_back(new_keyframe);
	}
	else {
		printf("creating new frame after the current keyframe\n");
		for (list<SgRbtNodes>::iterator iter = keyframeList.begin(), end = keyframeList.end(); iter != end; ++iter) {
			if (*iter == *g_currentKeyframe) {
				keyframeList.insert(++iter, new_keyframe);
				break;
			}
		}
	}
	// copy scene graph rbt data to new key frame
	dumpSgRbtNodes(g_world, new_keyframe);
	// set current key frame to this new frame
	*g_currentKeyframe = new_keyframe;
}

static void copyCurrentKeyframe() {
	// copy current key frame RBT data to the scene graph if current frame is defined
	if (keyframeList.empty()) {
		// no frames, do nothing
		printf("no key frames\n");
	}
	else {
		// TODO: how to do? maybe declare a new class? 
	}
}

static void deleteCurrentFrame() {
	// if our stack is empty, 
	if ((*g_currentKeyframe).empty()) {
		printf("there's no frame, nothing to delete\n");
	}
	else {
		printf("current frame was set, deleting now\n");
		SgRbtNodes temp = *g_currentKeyframe;
		for (list<SgRbtNodes>::iterator iter = keyframeList.begin(), end = keyframeList.end(); iter != end; ++iter) {
			if (*iter == *g_currentKeyframe) {
				// current frame is the first one, set current to right after it
				if (iter == keyframeList.begin()) {
					++iter;
					*g_currentKeyframe = *iter;
					// if there was actually only one frame and by incrementing we got to the end symbol, set to undefined instead
					if (iter == keyframeList.end()) {
						g_currentKeyframe = shared_ptr<SgRbtNodes>();
					}
					break; 
				}
				else {
					// current frame was not the first one, set current right before it
					--iter;
					*g_currentKeyframe = *iter;
					break;
				}
			}
		}
		// remove the old keyframe
		keyframeList.remove(temp);

		// copy RBT data from new current frame to scene graph
		copyCurrentKeyframe();
	}
}

static void shiftKeyframe(char* direction) {
	for (list<SgRbtNodes>::iterator iter = keyframeList.begin(), end = keyframeList.end(); iter != end; ++iter) {
		if (*iter == *g_currentKeyframe) {
			if (direction == "advance") {
				++iter;
				if (iter == keyframeList.end()) {
					printf("Cannot move forward, reached end of stack\n");
					return;
				}
			}
			else if (direction == "retreat") {
				if (iter == keyframeList.begin()) {
					printf("Cannot move backwards, at beginning of stack\n");
					return;
				}
				--iter;
			} 
			copyCurrentKeyframe();
			break;
		}
	}
}

static void updateScene() {
	// if empty stack, make a new keyframe
	if (keyframeList.empty()) {
		initializeNewKeyframe();
	}
	// otherwise, copy the scene graph rbt to current frame 
	else {
		dumpSgRbtNodes(g_world, *g_currentKeyframe);
	}
}

static void keyboard(const unsigned char key, const int x, const int y) {
  switch (key) {
  case 27:
    exit(0);                                  // ESC
  case 'h':
    cout << " ============== H E L P ==============\n\n"
    << "h\t\thelp menu\n"
    << "s\t\tsave screenshot\n"
    << "f\t\tToggle flat shading on/off.\n"
    << "p\t\tEnter picker mode\n"
    << "v\t\tCycle view\n"
    << "drag left mouse to rotate\n"
    << "drag right or middle mouse to translate\n" << endl;
    break;
  case 's':
    glFlush();
    writePpmScreenshot(g_windowWidth, g_windowHeight, "out.ppm");
    break;
  case 'f':
    g_activeShader ^= 1;
    break;
  case 'v':
    g_activeEye = ObjId((g_activeEye+1) % 3);
    cerr << "Active eye is " << g_objNames[g_activeEye] << endl;
    break;
  case 'm':
    g_activeCameraFrame = SkyMode((g_activeCameraFrame+1) % 2);
    cerr << "Editing sky eye w.r.t. " << (g_activeCameraFrame == WORLD_SKY ? "world-sky frame\n" : "sky-sky frame\n") << endl;
    break;
  case 'p':
	g_pickerMode = !g_pickerMode;
	cout << "Picking mode: " << g_pickerMode << "\n";
	break;
  case ' ':
	g_spaceDown = true;
	break;

// keyframe animation
  case 'c':
	printf("c was pressed\n");
	copyCurrentKeyframe();
	break;
  case 'u':
	printf("u was pressed\n");
	updateScene();
	break;
  case '>':
	printf("> was pressed\n");
	shiftKeyframe("advance");
	break;
  case '<':
	printf("< was pressed\n");
	shiftKeyframe("retreat");
	break;
  case 'd':
	printf("d was pressed\n");
	deleteCurrentFrame();
	break;
  case 'n':
	printf("n was pressed\n");
	initializeNewKeyframe();
	break;
  case 'i':
	printf("i was pressed\n");
	break;
  case 'w':
	printf("w was pressed\n");
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
  glutCreateWindow("Assignment 4");                       // title the window
 
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
		TORSO_WIDTH = 1,
		LEG_LEN = 0.7,
		LEG_THICK = 0.25,
		HEAD_WIDTH = 0.5,
		HEAD_HEIGHT = 0.5;
	const int NUM_JOINTS = 10,
		NUM_SHAPES = 10;

	struct JointDesc {
		int parent;
		float x, y, z;
	};

	JointDesc jointDesc[NUM_JOINTS] = {
			{ -1 }, // torso
			{ 0, TORSO_WIDTH / 2, TORSO_LEN / 2, 0 }, // upper right arm
			{ 1, ARM_LEN, 0, 0 }, // lower right arm
			{ 0, TORSO_WIDTH / 2, -TORSO_LEN / 2, 0 }, // upper right leg
			{ 3, 0, -LEG_LEN, 0 }, // lower right leg
			{ 0, -TORSO_WIDTH / 2, TORSO_LEN / 2, 0 }, // upper left arm
			{ 5, -ARM_LEN, 0, 0 }, // lower left arm
			{ 0, -TORSO_WIDTH / 2, -TORSO_LEN / 2, 0 }, // upper left leg
			{ 7, 0, -LEG_LEN, 0 }, // lower left leg
			{ 0, 0, TORSO_LEN / 2, 0 }, // head
	};

	struct ShapeDesc {
		int parentJointId;
		float x, y, z, sx, sy, sz;
		shared_ptr<Geometry> geometry;
	};

	ShapeDesc shapeDesc[NUM_SHAPES] = {
			{ 0, 0, 0, 0, TORSO_WIDTH, TORSO_LEN, TORSO_THICK, g_cube }, // torso
			{ 1, ARM_LEN / 2, 0, 0, ARM_LEN, ARM_THICK, ARM_THICK, g_cube }, // upper right arm
			{ 2, ARM_LEN / 2, 0, 0, ARM_LEN, ARM_THICK, ARM_THICK, g_cube }, // lower right arm
			{ 3, 0, -LEG_LEN / 2, 0, LEG_THICK, LEG_LEN, LEG_THICK, g_cube }, // upper right leg
			{ 4, 0, -LEG_LEN / 2, 0, LEG_THICK, LEG_LEN, LEG_THICK, g_cube }, // lower right leg
			{ 5, - ARM_LEN / 2, 0, 0, ARM_LEN, ARM_THICK, ARM_THICK, g_cube }, // upper left arm
			{ 6, - ARM_LEN / 2, 0, 0, ARM_LEN, ARM_THICK, ARM_THICK, g_cube }, // lower left arm
			{ 7, 0, -LEG_LEN / 2, 0, LEG_THICK, LEG_LEN, LEG_THICK, g_cube }, // upper left leg
			{ 8, 0, -LEG_LEN / 2, 0, LEG_THICK, LEG_LEN, LEG_THICK, g_cube }, // lower left leg
			{ 9, 0, HEAD_HEIGHT, 0, HEAD_WIDTH, HEAD_HEIGHT, HEAD_WIDTH, g_sphere }, // head
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
			Cvec3(shapeDesc[i].sx, shapeDesc[i].sy, shapeDesc[i].sz)));
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

	constructRobot(g_robot1Node, Cvec3(1, 0, 0)); // a Red robot
	constructRobot(g_robot2Node, Cvec3(0, 0, 1)); // a Blue robot

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