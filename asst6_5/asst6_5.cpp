////////////////////////////////////////////////////////////////////////
//
//   Harvard University
//   CS175 : Computer Graphics
//   Professor Steven Gortler
//
////////////////////////////////////////////////////////////////////////

#include <iostream>
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
#include "geometry.h"

// from asst4-snippets
#include "asstcommon.h"
#include "scenegraph.h"
#include "drawer.h"
#include "picker.h"

// for keyframe animation
#include <list>
#include "sgutils.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

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
static shared_ptr<Material> g_redDiffuseMat,
                            g_blueDiffuseMat,
                            g_bumpFloorMat,
                            g_arcballMat,
                            g_pickingMat,
                            g_lightMat;

shared_ptr<Material> g_overridingMaterial;

// --------- Geometry
typedef SgGeometryShapeNode MyShapeNode;
/// Vertex buffer and index buffer associated with the ground and cube geometry
static shared_ptr<Geometry> g_ground, g_cube, g_sphere;

// --------- Scene

static shared_ptr<SgRootNode> g_world;
static shared_ptr<SgRbtNode> g_skyNode, g_groundNode, g_robot1Node, g_robot2Node, g_light1Node, g_light2Node;
static shared_ptr<SgRbtNode> g_currentPickedRbtNode; // used later when you do pickingj

// for keyframe animation
static int g_msBetweenKeyFrames = 2000; // 2 seconds between keyframes
static int g_animateFramesPerSecond = 60; // frames to render per second
static bool g_isPlayingAnimation = false; // whether or not animation is currently playing

enum KeyframeId { ADVANCE = 0, RETREAT = 1 };
typedef std::vector<RigTForm> RigTFormVector;
list<RigTFormVector> keyframeList;  // list of RigTFormVector 
list<RigTFormVector>::iterator g_currentKeyframe = keyframeList.begin(); // pointer to vector of RigTFormVector that represent the current frame 

list<RigTFormVector>::iterator g_previousPlayingFromKeyframe; // the frame before the frame we are interpolating from
list<RigTFormVector>::iterator g_currentPlayingFromKeyframe; // the frame we are interpolating from 
list<RigTFormVector>::iterator g_currentPlayingToKeyframe; // the frame we are interpolating toward
list<RigTFormVector>::iterator g_followingPlayingToKeyframe; // the frame after the frame we are interpolating from 
int g_mostRecentPlayedKeyframe = 0; // keep track of when to switch the above two variables in interpolation 

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
 int ibLen, vbLen;
 getPlaneVbIbLen(vbLen, ibLen);
 
 //Temporary storage for cube Geometry
 vector<VertexPNTBX> vtx(vbLen);
 vector<unsigned short> idx(ibLen);

 makePlane(g_groundSize*2, vtx.begin(), idx.begin());
 g_ground.reset(new SimpleIndexedGeometryPNTBX(&vtx[0], &idx[0], vbLen, ibLen));

}

static void initCubes() {
  int ibLen, vbLen;
  getCubeVbIbLen(vbLen, ibLen);

  //Temporary storage for cube Geometry
  vector<VertexPNTBX> vtx(vbLen);
  vector<unsigned short> idx(ibLen);

  makeCube(1, vtx.begin(), idx.begin());
  g_cube.reset(new SimpleIndexedGeometryPNTBX(&vtx[0], &idx[0], vbLen, ibLen));
  
}

static void initSphere() {
  int ibLen, vbLen;
  getSphereVbIbLen(20, 10, vbLen, ibLen);

  // Temporary storage for sphere geometry
  vector<VertexPNTBX> vtx(vbLen);
  vector<unsigned short> idx(ibLen);
  makeSphere(1, 20, 10, vtx.begin(), idx.begin());
  g_sphere.reset(new SimpleIndexedGeometryPNTBX(&vtx[0], &idx[0], vtx.size(), idx.size()));
}

// takes a projection matrix and send to the the shaders
static void sendProjectionMatrix(Uniforms& uniforms, const Matrix4& projMatrix) {
  uniforms.put("uProjMatrix", projMatrix);
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

static void drawArcBall(Uniforms& uniforms) {
  // switch to wire frame mode
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  RigTForm arcballEye = inv(getNodeForEye(g_activeEye)->getRbt()) * getArcballRbt();
  Matrix4 MVM = rigTFormToMatrix(arcballEye) * Matrix4::makeScale(Cvec3(1,1,1)* g_arcballScale * g_arcballScreenRadius);
  sendModelViewNormalMatrix(uniforms, MVM, normalMatrix(MVM));

  g_arcballMat->draw(*g_sphere, uniforms);

}

static void drawStuff(bool picking) {
 
  // Declare an empty uniforms
  Uniforms uniforms;
  
  // if we are not translating, update arcball scale
  if (!(g_mouseMClickButton || (g_mouseLClickButton && g_mouseRClickButton) || (g_mouseLClickButton && !g_mouseRClickButton && g_spaceDown)))
    updateArcballScale();

  // build & send proj. matrix to vshader
  const Matrix4 projmat = makeProjectionMatrix();
  sendProjectionMatrix(uniforms, projmat);

  const RigTForm eyeRbt = getPathAccumRbt(g_world, getNodeForEye(g_activeEye));
  const RigTForm invEyeRbt = inv(eyeRbt);

  const Cvec3 eyeLight1 = getPathAccumRbt(g_world, g_light1Node).getTranslation();
  const Cvec3 eyeLight2 = getPathAccumRbt(g_world, g_light2Node).getTranslation();
  
  uniforms.put("uLight", eyeLight1);
  uniforms.put("uLight2", eyeLight2);
  
  if (!picking) {
	  Drawer drawer(invEyeRbt, uniforms);
	  g_world->accept(drawer);

	  if (g_displayArcball && shouldUseArcball()) {
		  drawArcBall(uniforms);
	  }
  } else {
	  Picker picker(invEyeRbt, uniforms);
	  
    // set overriding material to our picking material
    g_overridingMaterial = g_pickingMat;
    
    g_world->accept(picker);
    
    // unset overriding material 
    g_overridingMaterial.reset();
	  
    glFlush();
	  g_currentPickedRbtNode = picker.getRbtNodeAtXY(g_mouseClickX, g_mouseClickY);
	  
	  if (g_currentPickedRbtNode == g_groundNode) {
		  g_currentPickedRbtNode = shared_ptr<SgRbtNode>();   // set to NULL
		  cout << "clicked on ground\n";
	  }
  }
}

static void display() {
  //glUseProgram(g_shaderStates[g_activeShader]->program);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);                   // clear framebuffer color&depth

  drawStuff(g_pickerMode);

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
	//glUseProgram(g_shaderStates[PICKING_SHADER]->program);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	drawStuff(true);

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
	RigTFormVector new_keyframe;
	
	if (keyframeList.empty()) {
		// nothing initialized yet, do so now
		printf("initializing first frame\n");
		keyframeList.push_front(new_keyframe); 
		g_currentKeyframe = keyframeList.begin();
	}
	else if (g_currentKeyframe == keyframeList.end()) {
		printf("g_currentKeyframe is undefined, it returned the end character\n");
		// create a new keyframe at the end of the stack
		keyframeList.push_back(new_keyframe);
		// point to last element
		--g_currentKeyframe;
	}
	else {
		printf("current key frame initialized, creating a new keyframe after it\n");

		// iter must always be where g_currentKeyframe is
		++g_currentKeyframe;
		keyframeList.insert(g_currentKeyframe, new_keyframe);
		--g_currentKeyframe;
	}
	// copy scene graph rbt data to new key frame		
	
	std::vector<std::tr1::shared_ptr<SgRbtNode> > RbtNodes;
	dumpSgRbtNodes(g_world, RbtNodes);
	(*g_currentKeyframe).clear();
	for (int i = 0; i < RbtNodes.size(); ++i) {
		(*g_currentKeyframe).push_back(RbtNodes[i]->getRbt());
	}
}

static void copyCurrentKeyframe() {
	// copy current key frame RBT data to the scene graph if current frame is defined
	if (keyframeList.empty()) {
		// no frames, do nothing
		printf("no key frames\n");
	}
	else {
		// grab pointers to current scene data with the dump function, then copy our current key frame rbts there
		std::vector<std::tr1::shared_ptr<SgRbtNode> > RbtNodes;
		dumpSgRbtNodes(g_world, RbtNodes);
		
		std::vector<shared_ptr<SgRbtNode> >::iterator iter1 = RbtNodes.begin();
		std::vector<RigTForm>::iterator iter2 = (*g_currentKeyframe).begin();

		while (true) {
			if (iter1 != RbtNodes.end() && iter2 != (*g_currentKeyframe).end()) {
				(**iter1).setRbt(*iter2);
				++iter1;
				++iter2;
			}
			else {
				break;
			}
		}
	}
}

static void deleteCurrentFrame() {
	// if our stack is empty, 
	if (keyframeList.empty()) {
		printf("there's no frame, nothing to delete\n");
	}
	else {
		printf("current frame was set, deleting now\n");
		
		list<RigTFormVector>::iterator temp = g_currentKeyframe;
		if (g_currentKeyframe == keyframeList.begin()) {
			if (++g_currentKeyframe == keyframeList.end()) {
         g_currentKeyframe = keyframeList.end();
			}
		}
		else {
			--g_currentKeyframe;
		}
		
		keyframeList.erase(temp);
		
		// update scene
		copyCurrentKeyframe();
	}
}

static RigTForm evaluateBezier(RigTForm from, RigTForm to, RigTForm d, RigTForm e, float alpha) {
	Cvec3 f_t = lerp(from.getTranslation(), d.getTranslation(), alpha);
	Cvec3 g_t = lerp(d.getTranslation(), e.getTranslation(), alpha);
	Cvec3 h_t = lerp(e.getTranslation(), to.getTranslation(), alpha);
	Cvec3 m_t = lerp(f_t, g_t, alpha);
	Cvec3 n_t = lerp(g_t, h_t, alpha);
	Cvec3 c_t = lerp(m_t, n_t, alpha);

	Quat f_q = slerp(from.getRotation(), d.getRotation(), alpha);
	Quat g_q = slerp(d.getRotation(), e.getRotation(), alpha);
	Quat h_q = slerp(e.getRotation(), to.getRotation(), alpha);
	Quat m_q = slerp(f_q, g_q, alpha);
	Quat n_q = slerp(g_q, h_q, alpha);
	Quat c_q = slerp(m_q, n_q, alpha);

	return RigTForm(c_t, c_q);
}

static RigTForm evaluateCatmull_Rom(RigTForm prev, RigTForm from, RigTForm to, RigTForm post, float alpha) {
	// TODO: sanity check for quat negation if first coordinate of the part of D is negative
	Cvec3 BC_CvecD = (to.getTranslation() - prev.getTranslation()) * (1 / 6) + from.getTranslation();
	Cvec3 BC_CvecE = (post.getTranslation() - from.getTranslation()) * (-1 / 6) + to.getTranslation();
	
	//Quat d_pow = to.getRotation() * inv(prev.getRotation());
	//Quat e_pow = post.getRotation() * inv(from.getRotation());
	/*
	if (d_pow[0] < 0) {
		d_pow = d_pow * (-1);
		printf("%i", d_pow[0]);
	}*/

	//Quat BC_QuatD = (to.getRotation() * inv(prev.getRotation())).pow(1 / 6) * from.getRotation();
	//Quat BC_QuatE = (post.getRotation() * inv(from.getRotation())).pow(-1 / 6) * to.getRotation();
	//Quat BC_QuatD = d_pow.pow(1 / 6) * from.getRotation();
	//Quat BC_QuatE = e_pow.pow(-1 / 6) * to.getRotation();
	
	Quat BC_QuatD = Quat(1, 0, 0, 0);
	Quat BC_QuatE = Quat(1, 0, 0, 0);
	RigTForm BC_D = RigTForm(BC_CvecD, BC_QuatD);
	RigTForm BC_E = RigTForm(BC_CvecE, BC_QuatE);

	return evaluateBezier(from, to, BC_D, BC_E, alpha);
}

bool interpolateAndDisplay(float t) {

  // get alpha for slerping and lerping 
  float alpha = t - floor(t);

  int keyframe = floor(t);
  
  if (keyframe == 0) {
	  g_previousPlayingFromKeyframe = keyframeList.begin();

	  // for Catmull Rom we begin interpolating from the second keyframe
	  g_currentPlayingFromKeyframe = keyframeList.begin();
	  ++g_currentPlayingFromKeyframe;

	  g_currentPlayingToKeyframe = g_currentPlayingFromKeyframe;
	  ++g_currentPlayingToKeyframe;

	  // for Catmull Rom we keep track of the keyframe after the frame we're playing to
	  g_followingPlayingToKeyframe = g_currentPlayingToKeyframe;
	  ++g_followingPlayingToKeyframe;
  }

  // interpolate between the next pair of keyframes if we have passed them
  if (keyframe != g_mostRecentPlayedKeyframe) {
	  g_mostRecentPlayedKeyframe = keyframe;
	  ++g_previousPlayingFromKeyframe;
	  ++g_currentPlayingFromKeyframe;
	  ++g_currentPlayingToKeyframe;
	  ++g_followingPlayingToKeyframe;
  }

  // get the second to last keyframe
  list<RigTFormVector>::iterator secondToLastKeyframe = keyframeList.end();
  --secondToLastKeyframe;
  --secondToLastKeyframe;

  // stop animation when we reach the second to last keyframe
  if (g_currentPlayingToKeyframe == secondToLastKeyframe) {
	  g_isPlayingAnimation = false;
	  g_mostRecentPlayedKeyframe = 0;
	  return true;
  }

  // make a list for interpolations
  vector<RigTForm> interpolations;

  //interpolate between the two
  int size = (*g_currentPlayingFromKeyframe).size();
  for (int i = 0; i < size; i++) {
	  // Catmull Rom Interpolation
	  RigTForm interpolation = evaluateCatmull_Rom((*g_previousPlayingFromKeyframe)[i], 
		  (*g_currentPlayingFromKeyframe)[i], (*g_currentPlayingToKeyframe)[i], (*g_followingPlayingToKeyframe)[i], alpha);
	  interpolations.push_back(interpolation);
  }
  /*
  int keyframe = floor(t);
  if (keyframe == 0) {
    g_currentPlayingFromKeyframe = keyframeList.begin();
    g_currentPlayingToKeyframe = g_currentPlayingFromKeyframe;
    ++g_currentPlayingToKeyframe;
  }
 
  // interpolate between the next pair of keyframes if we have passed them
  if (keyframe != g_mostRecentPlayedKeyframe) {
    g_mostRecentPlayedKeyframe = keyframe;
    ++g_currentPlayingFromKeyframe; 
    ++g_currentPlayingToKeyframe; 
  }

  // stop animation when we are done
  if (g_currentPlayingToKeyframe == keyframeList.end()) {
    g_isPlayingAnimation = false;
    g_mostRecentPlayedKeyframe = 0;
    return true;
  }
  
  // make a list for interpolations
  vector<RigTForm> interpolations;

  //interpolate between the two
  int size = (*g_currentPlayingFromKeyframe).size(); 
  for (int i = 0; i < size; i++) { 
    // lerp translation and slerp rotation
    RigTForm interpolation = 
      RigTForm(
        lerp(
          (*g_currentPlayingFromKeyframe)[i].getTranslation(),
          (*g_currentPlayingToKeyframe)[i].getTranslation(),
          alpha
        ),
        slerp(
          (*g_currentPlayingFromKeyframe)[i].getRotation(),
          (*g_currentPlayingToKeyframe)[i].getRotation(),
          alpha
        )
      );
    interpolations.push_back(interpolation); 
  }*/
 
  // get current nodes in scene
	std::vector<std::tr1::shared_ptr<SgRbtNode> > current_nodes;
	dumpSgRbtNodes(g_world, current_nodes);

  // apply the transformation
  for (int i = 0; i < current_nodes.size(); i++) {
    current_nodes[i]->setRbt(interpolations[i]);  
  }
 
  // redraw scene 
  glutPostRedisplay();

  // the animation is not done, so return false
  return false;
} 

static void animateTimerCallback(int ms) {
  float t = (float) ms / (float) g_msBetweenKeyFrames;

  bool endReached = interpolateAndDisplay(t);
  if (!endReached) {
          glutTimerFunc(1000/g_animateFramesPerSecond,
                        animateTimerCallback,
                        ms + 1000/g_animateFramesPerSecond
          );
  } else {
    g_isPlayingAnimation = false; 
  }
}

static void controlAnimation() {
  if (!g_isPlayingAnimation) {
    // check that there are 4 keyframes
    if (keyframeList.size() < 4) {
      cout << "Unable to start animation -- you only have "
           << keyframeList.size()
           << " keyframes, but we need at least 4!\n";
      return;
    }
    // toggle global flag for animation play
    g_isPlayingAnimation = true;

    // call animateTimerCallback(0)
    animateTimerCallback(0);

  } else {
    // animation is alrady playing, so stop 
    g_isPlayingAnimation = false; 
  }
}

static void shiftKeyframe(KeyframeId i) {
	if (i == ADVANCE && g_currentKeyframe != keyframeList.end()) {
		++g_currentKeyframe;
		if (g_currentKeyframe == keyframeList.end()) {
			printf("Cannot move forward, reached end of stack\n");
			--g_currentKeyframe;
			return;
		}
	}
	else if (i == RETREAT) {
		if (g_currentKeyframe == keyframeList.begin()) {
			printf("Cannot move backwards, at beginning of stack\n");
			return;
		}
		--g_currentKeyframe;
	}
	copyCurrentKeyframe();
}

static void updateScene() {
	// if empty stack, make a new keyframe
	if (keyframeList.empty()) {
		initializeNewKeyframe();
	}
	// otherwise, copy the scene graph rbt to current frame 
	else {
		std::vector<std::tr1::shared_ptr<SgRbtNode> > RbtNodes;
		dumpSgRbtNodes(g_world, RbtNodes);
		(*g_currentKeyframe).clear();
		for (int i = 0; i < RbtNodes.size(); ++i) {
			(*g_currentKeyframe).push_back(RbtNodes[i]->getRbt());
		}
	}
}

static void readFrameDataFromFile() {
	ifstream myfile; 
	string line; 
	myfile.open("animation.txt");
	getline(myfile, line);
	int num_frames = atoi(line.c_str());
	getline(myfile, line);
	int num_nodes = atoi(line.c_str());
	printf("%d %d\n", num_frames, num_nodes);
	
	keyframeList.clear();
	for (int i = 0; i < num_frames; i++) {
		RigTFormVector frame; 
		for (int j = 0; j < num_nodes; j++) {
			RigTForm node; 
			Cvec3 translation;
			Quat rotation;
			std::string line;
			std::getline(myfile, line);
			std::stringstream stream(line);
			for (int k = 0; k < 3; k++) {
				double cvec_double;
				stream >> cvec_double;
				translation[k] = cvec_double;
			}
			for (int l = 0; l < 4; l++) {
				double quat_double;
				stream >> quat_double;
				rotation[l] = quat_double;
			}
			printf("\n");
			node.setTranslation(translation); 
			node.setRotation(rotation);
			frame.push_back(node);
		}
		keyframeList.push_back(frame);
	}
	g_currentKeyframe = keyframeList.begin();

	copyCurrentKeyframe();

	myfile.close();
}

static void writeFrameDataToFile() {
	if (!keyframeList.empty()) {
		ofstream myfile;
		myfile.open("animation.txt");
		// myfile << "Writing this to a file.\n";
		// myfile << "testing, testing \n";
		string output = "";
		int num_frames = keyframeList.size();
		list<RigTFormVector>::iterator iter = keyframeList.begin();
		int num_nodes = (*iter).size();
    
    std::ostringstream s;
    s << num_frames
      << "\n"
      << num_nodes
      << "\n";
    output.append(s.str());
		for (iter; iter != keyframeList.end(); ++iter) {
			RigTFormVector scene = (*iter);
			for (int i = 0; i < scene.size(); ++i) {
				RigTForm rigTForm = scene[i];
				Cvec3 translation = rigTForm.getTranslation();
				Quat rotation = rigTForm.getRotation();
				for (int j = 0; j < 3; j++) {
				  s.str(""); 
          s.clear(); 
          s << translation[j] << " ";
          output.append(s.str());
				}
				for (int k = 0; k < 4; k++) {
				  s.str(""); 
          s.clear(); 
          s << rotation[k] << " ";
          output.append(s.str());
				}
				output.append("\n");
			}
		}
		myfile << output;

		myfile.close();
		printf("Writing to animation.txt\n");
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
    << "drag right or middle mouse to translate\n\n"
    << "\t---ANIMATION CONTROLS---\n"
    << "n\t\tMake new keyframe\n"
    << "c\t\tCopy current keyframe\n"
    << "u\t\tUpdate keyframe\n"
    << "d\t\tDelete current keyframe\n" 
    << ">\t\tGo to next keyframe\n"
    << "<\t\tGo to previous keyframe\n"
    << "w\t\tWrite keyframes to file\n"
    << "i\t\tRead keyframes from file\n"
    << "y\t\tPlay/Stop playing keyframes\n"
    << "+\t\tSpeed up animation\n"
    << "-\t\tSlow down animation\n"
    << endl;
    break;
  case 's':
    glFlush();
    writePpmScreenshot(g_windowWidth, g_windowHeight, "out.ppm");
    break;
  /*case 'f':
    g_activeShader ^= 1;
    break;*/
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
	  shiftKeyframe(ADVANCE);
	  break;
  case '<':
	  printf("< was pressed\n");
	  shiftKeyframe(RETREAT);
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
	readFrameDataFromFile();
	break;
  case 'w':
	printf("w was pressed\n");
	writeFrameDataToFile();
	break;
  case 'y':
    printf("y was pressed\n");
    controlAnimation();
    break;
  case '+':
    printf("+ was pressed\n");
    g_msBetweenKeyFrames = max(100, g_msBetweenKeyFrames - 100);
    cout << "The new speed is: " << g_msBetweenKeyFrames << "ms between keyframes\n"; 
   break;
  case '-':
    printf("- was pressed\n"); 
    g_msBetweenKeyFrames = min(1000, g_msBetweenKeyFrames + 100);
    cout << "The new speed is: " << g_msBetweenKeyFrames << "ms between keyframes\n"; 
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

static void initMaterials() {
  // Create some prototype materials
  Material diffuse("./shaders/basic-gl3.vshader", "./shaders/diffuse-gl3.fshader");
  Material solid("./shaders/basic-gl3.vshader", "./shaders/solid-gl3.fshader");

  // copy diffuse prototype and set red color
  g_redDiffuseMat.reset(new Material(diffuse));
  g_redDiffuseMat->getUniforms().put("uColor", Cvec3f(1, 0, 0));

  // copy diffuse prototype and set blue color
  g_blueDiffuseMat.reset(new Material(diffuse));
  g_blueDiffuseMat->getUniforms().put("uColor", Cvec3f(0, 0, 1));

  // normal mapping material
  g_bumpFloorMat.reset(new Material("./shaders/normal-gl3.vshader", "./shaders/normal-gl3.fshader"));
  g_bumpFloorMat->getUniforms().put("uTexColor", shared_ptr<ImageTexture>(new ImageTexture("Fieldstone.ppm", true)));
  g_bumpFloorMat->getUniforms().put("uTexNormal", shared_ptr<ImageTexture>(new ImageTexture("FieldstoneNormal.ppm", false)));

  // copy solid prototype, and set to wireframed rendering
  g_arcballMat.reset(new Material(solid));
  g_arcballMat->getUniforms().put("uColor", Cvec3f(0.27f, 0.82f, 0.35f));
  g_arcballMat->getRenderStates().polygonMode(GL_FRONT_AND_BACK, GL_LINE);

  // copy solid prototype, and set to color white
  g_lightMat.reset(new Material(solid));
  g_lightMat->getUniforms().put("uColor", Cvec3f(1, 1, 1));

  // pick shader
  g_pickingMat.reset(new Material("./shaders/basic-gl3.vshader", "./shaders/pick-gl3.fshader"));
}

static void initGeometry() {
  initGround();
  initCubes();
  initSphere();
}

static void constructRobot(shared_ptr<SgTransformNode> base, shared_ptr<Material> material) {

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
			material,
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
		new MyShapeNode(g_ground, g_bumpFloorMat, Cvec3(0, g_groundY, 0))));

	g_robot1Node.reset(new SgRbtNode(RigTForm(Cvec3(-2, 1, 0))));
	g_robot2Node.reset(new SgRbtNode(RigTForm(Cvec3(2, 1, 0))));

	constructRobot(g_robot1Node, g_redDiffuseMat); // a Red robot
	constructRobot(g_robot2Node, g_blueDiffuseMat); // a Blue robot

  /* lights */
  Cvec3 light1_pos = Cvec3(4.5, 3.0, 5.0);
  Cvec3 light2_pos = Cvec3(-4.5, 0, -5.0);

  g_light1Node.reset(new SgRbtNode(RigTForm(light1_pos)));
  g_light1Node->addChild(shared_ptr<MyShapeNode>(
    new MyShapeNode(g_sphere, g_lightMat, Cvec3(0, 0, 0))));

  g_light2Node.reset(new SgRbtNode(RigTForm(light2_pos)));
  g_light2Node->addChild(shared_ptr<MyShapeNode>(
    new MyShapeNode(g_sphere, g_lightMat, Cvec3(0, 0, 0))));

	g_world->addChild(g_skyNode);
	g_world->addChild(g_groundNode);

  g_world->addChild(g_light1Node);
  g_world->addChild(g_light2Node);

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
//    initShaders();
    initMaterials(); 
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
