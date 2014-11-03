8a9
> #include <iostream>
11d11
< #include <list>
14,15d13
< #include <map>
< #include <fstream>
16a15
> #include <algorithm>
36,37c35
< #include "scenegraph.h"
< #include "sgutils.h"
---
> #include "geometry.h"
38a37
> // from asst4-snippets
39a39
> #include "scenegraph.h"
42a43,50
> // for keyframe animation
> #include <list>
> #include "sgutils.h"
> #include <iostream>
> #include <fstream>
> #include <string>
> #include <sstream>
> 
61d68
< const bool g_Gl2Compatible = false;
62a70
> const bool g_Gl2Compatible = false;
71a80
> enum ObjId {SKY=0, OBJECT0=1, OBJECT1=2, LIMB = 3};
73a83,84
> static const char * const g_objNames[] = {"Sky", "Object 0", "Object 1"};
> 
81a93
> static ObjId g_activeEye = SKY;
88,90c100,106
< static bool g_pickingMode = false;
< 
< static bool g_playingAnimation = false;
---
> // asst-snippets
> static shared_ptr<Material> g_redDiffuseMat,
>                             g_blueDiffuseMat,
>                             g_bumpFloorMat,
>                             g_arcballMat,
>                             g_pickingMat,
>                             g_lightMat;
92,105c108
< // -------- Shaders
< static const int g_numShaders = 3, g_numRegularShaders = 2;
< static const int PICKING_SHADER = 2;
< static const char * const g_shaderFiles[g_numShaders][2] = {
<   {"./shaders/basic-gl3.vshader", "./shaders/diffuse-gl3.fshader"},
<   {"./shaders/basic-gl3.vshader", "./shaders/solid-gl3.fshader"},
<   {"./shaders/basic-gl3.vshader", "./shaders/pick-gl3.fshader"}
< };
< static const char * const g_shaderFilesGl2[g_numShaders][2] = {
<   {"./shaders/basic-gl2.vshader", "./shaders/diffuse-gl2.fshader"},
<   {"./shaders/basic-gl2.vshader", "./shaders/solid-gl2.fshader"},
<   {"./shaders/basic-gl2.vshader", "./shaders/pick-gl2.fshader"}
< };
< static vector<shared_ptr<ShaderState> > g_shaderStates; // our global shader states
---
> shared_ptr<Material> g_overridingMaterial;
108,183c111,112
< 
< 
< // Macro used to obtain relative offset of a field within a struct
< #define FIELD_OFFSET(StructType, field) ((GLvoid*)offsetof(StructType, field))
< 
< // A vertex with floating point position and normal
< struct VertexPN {
<   Cvec3f p, n;
< 
<   VertexPN() {}
<   VertexPN(float x, float y, float z,
<            float nx, float ny, float nz)
<     : p(x,y,z), n(nx, ny, nz)
<   {}
< 
<   // Define copy constructor and assignment operator from GenericVertex so we can
<   // use make* functions from geometrymaker.h
<   VertexPN(const GenericVertex& v) {
<     *this = v;
<   }
< 
<   VertexPN& operator = (const GenericVertex& v) {
<     p = v.pos;
<     n = v.normal;
<     return *this;
<   }
< };
< 
< struct Geometry {
<   GlBufferObject vbo, ibo;
<   GlArrayObject vao;
<   int vboLen, iboLen;
< 
<   Geometry(VertexPN *vtx, unsigned short *idx, int vboLen, int iboLen) {
<     this->vboLen = vboLen;
<     this->iboLen = iboLen;
< 
<     // Now create the VBO and IBO
<     glBindBuffer(GL_ARRAY_BUFFER, vbo);
<     glBufferData(GL_ARRAY_BUFFER, sizeof(VertexPN) * vboLen, vtx, GL_STATIC_DRAW);
< 
<     glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
<     glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned short) * iboLen, idx, GL_STATIC_DRAW);
<   }
< 
<   void draw(const ShaderState& curSS) {
<     // bind the object's VAO
<     glBindVertexArray(vao);
< 
<     // Enable the attributes used by our shader
<     safe_glEnableVertexAttribArray(curSS.h_aPosition);
<     safe_glEnableVertexAttribArray(curSS.h_aNormal);
< 
<     // bind vbo
<     glBindBuffer(GL_ARRAY_BUFFER, vbo);
<     safe_glVertexAttribPointer(curSS.h_aPosition, 3, GL_FLOAT, GL_FALSE, sizeof(VertexPN), FIELD_OFFSET(VertexPN, p));
<     safe_glVertexAttribPointer(curSS.h_aNormal, 3, GL_FLOAT, GL_FALSE, sizeof(VertexPN), FIELD_OFFSET(VertexPN, n));
< 
<     // bind ibo
<     glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
< 
<     // draw!
<     glDrawElements(GL_TRIANGLES, iboLen, GL_UNSIGNED_SHORT, 0);
< 
<     // Disable the attributes used by our shader
<     safe_glDisableVertexAttribArray(curSS.h_aPosition);
<     safe_glDisableVertexAttribArray(curSS.h_aNormal);
< 
<     // disable VAO
<     glBindVertexArray(NULL);
<   }
< };
< 
< typedef SgGeometryShapeNode<Geometry> MyShapeNode;
< 
< // Vertex buffer and index buffer associated with the ground and cube geometry
---
> typedef SgGeometryShapeNode MyShapeNode;
> /// Vertex buffer and index buffer associated with the ground and cube geometry
188,189d116
< static const Cvec3 g_light1(2.0, 3.0, 14.0), g_light2(-2, -3.0, -5.0);  // define two lights positions in world space
< 
191,195c118,119
< static shared_ptr<SgRbtNode> g_skyNode, g_groundNode, g_robot1Node, g_robot2Node;
< 
< static shared_ptr<SgRbtNode> g_currentCameraNode;
< static shared_ptr<SgRbtNode> g_currentPickedRbtNode;
< 
---
> static shared_ptr<SgRbtNode> g_skyNode, g_groundNode, g_robot1Node, g_robot2Node, g_light1Node, g_light2Node;
> static shared_ptr<SgRbtNode> g_currentPickedRbtNode; // used later when you do pickingj
197,213c121,124
< static RigTForm evaluateBezier(RigTForm from, RigTForm to, RigTForm d, RigTForm e, double alpha) {
< 	Cvec3 f_t = lerp(from.getTranslation(), d.getTranslation(), alpha);
< 	Cvec3 g_t = lerp(d.getTranslation(), e.getTranslation(), alpha);
< 	Cvec3 h_t = lerp(e.getTranslation(), to.getTranslation(), alpha);
< 	Cvec3 m_t = lerp(f_t, g_t, alpha);
< 	Cvec3 n_t = lerp(g_t, h_t, alpha);
< 	Cvec3 c_t = lerp(m_t, n_t, alpha);
< 
< 	Quat f_q = slerp(from.getRotation(), d.getRotation(), alpha);
< 	Quat g_q = slerp(d.getRotation(), e.getRotation(), alpha);
< 	Quat h_q = slerp(e.getRotation(), to.getRotation(), alpha);
< 	Quat m_q = slerp(f_q, g_q, alpha);
< 	Quat n_q = slerp(g_q, h_q, alpha);
< 	Quat c_q = slerp(m_q, n_q, alpha);
< 
< 	return RigTForm(c_t, c_q);
< }
---
> // for keyframe animation
> static int g_msBetweenKeyFrames = 2000; // 2 seconds between keyframes
> static int g_animateFramesPerSecond = 60; // frames to render per second
> static bool g_isPlayingAnimation = false; // whether or not animation is currently playing
215,218c126,135
< static RigTForm evaluateCatmull_Rom(RigTForm prev, RigTForm from, RigTForm to, RigTForm post, double alpha) {
< 	
< 	Cvec3 BC_CvecD = (to.getTranslation() - prev.getTranslation()) * (1 / 6) + from.getTranslation();
< 	Cvec3 BC_CvecE = (post.getTranslation() - from.getTranslation()) * (-1 / 6) + to.getTranslation();
---
> enum KeyframeId { ADVANCE = 0, RETREAT = 1 };
> typedef std::vector<RigTForm> RigTFormVector;
> list<RigTFormVector> keyframeList;  // list of RigTFormVector 
> list<RigTFormVector>::iterator g_currentKeyframe = keyframeList.begin(); // pointer to vector of RigTFormVector that represent the current frame 
> 
> list<RigTFormVector>::iterator g_previousPlayingFromKeyframe; // the frame before the frame we are interpolating from
> list<RigTFormVector>::iterator g_currentPlayingFromKeyframe; // the frame we are interpolating from 
> list<RigTFormVector>::iterator g_currentPlayingToKeyframe; // the frame we are interpolating toward
> list<RigTFormVector>::iterator g_followingPlayingToKeyframe; // the frame after the frame we are interpolating from 
> int g_mostRecentPlayedKeyframe = 0; // keep track of when to switch the above two variables in interpolation 
220,221c137,140
< 	Quat d_pow = to.getRotation() * inv(prev.getRotation());
< 	Quat e_pow = post.getRotation() * inv(from.getRotation());
---
> static const Cvec3 g_light1(2.0, 3.0, 14.0), g_light2(-2, -3.0, -5.0);  // define two lights positions in world space
> static RigTForm g_skyRbt = RigTForm(Cvec3(0.0, 0.25, 4.0));
> static RigTForm g_objectRbt[2] = {RigTForm(Cvec3(-1,0,0)), RigTForm(Cvec3(1,0,0))};
> static Cvec3f g_objectColors[2] = {Cvec3f(1, 0, 0), Cvec3f(0, 0, 1)};
223,225c142
< 	if (d_pow[0] < 0) {
< 		d_pow = cn(d_pow);
< 	}
---
> static bool g_pickerMode = false;
227,228c144
< 	Quat BC_QuatD = pow(d_pow, 1.0 / 6.0) * from.getRotation();
< 	Quat BC_QuatE = pow(e_pow, -1.0 / 6.0) * to.getRotation();
---
> ///////////////// END OF G L O B A L S //////////////////////////////////////////////////
230,231d145
< 	RigTForm BC_D = RigTForm(BC_CvecD, BC_QuatD);
< 	RigTForm BC_E = RigTForm(BC_CvecE, BC_QuatE);
233c147,158
< 	return evaluateBezier(from, to, BC_D, BC_E, alpha);
---
> static shared_ptr<SgRbtNode> getNodeForEye(ObjId i) {
>    switch (i) { 
>     case SKY:
>       return g_skyNode;
>     case OBJECT0:
>       return g_robot1Node;
>     case OBJECT1:
>       return g_robot2Node;
>     case LIMB:
>       cout << "Shouldn't have this eye ever, something has gone wrong";
>       exit(1); 
>    }
235,299d159
< class Animator {
< public:
<   typedef vector<shared_ptr<SgRbtNode> > SgRbtNodes;
<   typedef vector<RigTForm> KeyFrame;
<   typedef list<KeyFrame> KeyFrames;
<   typedef KeyFrames::iterator KeyFrameIter;
< 
< private:
<   SgRbtNodes nodes_;
<   KeyFrames keyFrames_;
< 
< public:
<   void attachSceneGraph(shared_ptr<SgNode> root) {
<     nodes_.clear();
<     keyFrames_.clear();
<     dumpSgRbtNodes(root, nodes_);
<   }
< 
<   void loadAnimation(const char *filename) {
<     ifstream f(filename, ios::binary);
<     if (!f)
<       throw runtime_error(string("Cannot load ") + filename);
<     int numFrames, numRbtsPerFrame;
<     f >> numFrames >> numRbtsPerFrame;
<     if (numRbtsPerFrame != nodes_.size()) {
<       cerr << "Number of Rbt per frame in " << filename
<            <<" does not match number of SgRbtNodes in the current scene graph.";
<       return;
<     }
< 
<     Cvec3 t;
<     Quat r;
<     keyFrames_.clear();
<     for (int i = 0; i < numFrames; ++i) {
<       keyFrames_.push_back(KeyFrame());
<       keyFrames_.back().reserve(numRbtsPerFrame);
<       for (int j = 0; j < numRbtsPerFrame; ++j) {
<         f >> t[0] >> t[1] >> t[2] >> r[0] >> r[1] >> r[2] >> r[3];
<         keyFrames_.back().push_back(RigTForm(t, r));
<       }
<     }
<   }
< 
<   void saveAnimation(const char *filename) {
<     ofstream f(filename, ios::binary);
<     int numRbtsPerFrame = nodes_.size();
<     f << getNumKeyFrames() << ' ' << numRbtsPerFrame << '\n';
<     for (KeyFrames::const_iterator frameIter = keyFrames_.begin(), e = keyFrames_.end(); frameIter != e; ++frameIter) {
<       for (int j = 0; j < numRbtsPerFrame; ++j) {
<         const RigTForm& rbt = (*frameIter)[j];
<         const Cvec3& t = rbt.getTranslation();
<         const Quat& r = rbt.getRotation();
<         f << t[0] << ' ' << t[1] << ' ' << t[2] << ' '
<         << r[0] << ' ' << r[1] << ' ' << r[2] << ' ' << r[3] << '\n';
<       }
<     }
<   }
< 
<   int getNumKeyFrames() const {
<     return keyFrames_.size();
<   }
< 
<   int getNumRbtNodes() const {
<     return nodes_.size();
<   }
301,369c161,167
<   // t can be in the range [0, keyFrames_.size()-3]. Fractional amount like 1.5 is allowed.
<   void animate(double t) {
<     if (t < 0 || t > keyFrames_.size() - 3)
<       throw runtime_error("Invalid animation time parameter. Must be in the range [0, numKeyFrames - 3]");
< 
<     t += 1; // interpret the key frames to be at t= -1, 0, 1, 2, ...
<     const int integralT = int(floor(t));
<     const double fraction = t - integralT;
< 
<     KeyFrameIter f0 = getNthKeyFrame(integralT), f1 = f0, f_1 = f0, f2 = f1;
<     ++f1;
< 	--f_1;
< 	++f2;
< 
<     for (int i = 0, n = nodes_.size(); i < n; ++i) {
< 		// Evaluate Catmull_Rom
< 		nodes_[i]->setRbt(evaluateCatmull_Rom((*f_1)[i],
< 			(*f0)[i], (*f1)[i], (*f2)[i], fraction));
<     }
<   }
< 
<   KeyFrameIter keyFramesBegin() {
<     return keyFrames_.begin();
<   }
< 
<   KeyFrameIter keyFramesEnd() {
<     return keyFrames_.end();
<   }
< 
<   KeyFrameIter getNthKeyFrame(int n) {
<     KeyFrameIter frameIter = keyFrames_.begin();
<     advance(frameIter, n);
<     return frameIter;
<   }
< 
<   void deleteKeyFrame(KeyFrameIter keyFrameIter) {
<     keyFrames_.erase(keyFrameIter);
<   }
< 
<   void pullKeyFrameFromSg(KeyFrameIter keyFrameIter) {
<     for (int i = 0, n = nodes_.size(); i < n; ++i) {
<       (*keyFrameIter)[i] = nodes_[i]->getRbt();
<     }
<   }
< 
<   void pushKeyFrameToSg(KeyFrameIter keyFrameIter) {
<     for (int i = 0, n = nodes_.size(); i < n; ++i) {
<       nodes_[i]->setRbt((*keyFrameIter)[i]);
<     }
<   }
< 
<   KeyFrameIter insertEmptyKeyFrameAfter(KeyFrameIter beforeFrame) {
<     if (beforeFrame != keyFrames_.end())
<       ++beforeFrame;
< 
<     KeyFrameIter frameIter = keyFrames_.insert(beforeFrame, KeyFrame());
<     frameIter->resize(nodes_.size());
<     return frameIter;
<   }
< 
< };
< 
< static int g_msBetweenKeyFrames = 2000; // 2 seconds between keyframes
< static int g_animateFramesPerSecond = 60; // frames to render per second during animation playback
< 
< 
< static Animator g_animator;
< static Animator::KeyFrameIter g_curKeyFrame;
< static int g_curKeyFrameNum;
---
> static void initGround() {
>  int ibLen, vbLen;
>  getPlaneVbIbLen(vbLen, ibLen);
>  
>  //Temporary storage for cube Geometry
>  vector<VertexPNTBX> vtx(vbLen);
>  vector<unsigned short> idx(ibLen);
371c169,170
< ///////////////// END OF G L O B A L S //////////////////////////////////////////////////
---
>  makePlane(g_groundSize*2, vtx.begin(), idx.begin());
>  g_ground.reset(new SimpleIndexedGeometryPNTBX(&vtx[0], &idx[0], vbLen, ibLen));
373,382d171
< static void initGround() {
<   // A x-z plane at y = g_groundY of dimension [-g_groundSize, g_groundSize]^2
<   VertexPN vtx[4] = {
<     VertexPN(-g_groundSize, g_groundY, -g_groundSize, 0, 1, 0),
<     VertexPN(-g_groundSize, g_groundY,  g_groundSize, 0, 1, 0),
<     VertexPN( g_groundSize, g_groundY,  g_groundSize, 0, 1, 0),
<     VertexPN( g_groundSize, g_groundY, -g_groundSize, 0, 1, 0)
<   };
<   unsigned short idx[] = {0, 1, 2, 0, 2, 3};
<   g_ground.reset(new Geometry(&vtx[0], &idx[0], 4, 6));
389,391c178,179
< 
<   // Temporary storage for cube geometry
<   vector<VertexPN> vtx(vbLen);
---
>   //Temporary storage for cube Geometry
>   vector<VertexPNTBX> vtx(vbLen);
395c183,184
<   g_cube.reset(new Geometry(&vtx[0], &idx[0], vbLen, ibLen));
---
>   g_cube.reset(new SimpleIndexedGeometryPNTBX(&vtx[0], &idx[0], vbLen, ibLen));
>   
403c192
<   vector<VertexPN> vtx(vbLen);
---
>   vector<VertexPNTBX> vtx(vbLen);
406,410c195
<   g_sphere.reset(new Geometry(&vtx[0], &idx[0], vtx.size(), idx.size()));
< }
< 
< static void initRobots() {
<   // Init whatever geometry needed for the robots
---
>   g_sphere.reset(new SimpleIndexedGeometryPNTBX(&vtx[0], &idx[0], vtx.size(), idx.size()));
414,417c199,200
< inline void sendProjectionMatrix(const ShaderState& curSS, const Matrix4& projMatrix) {
<   GLfloat glmatrix[16];
<   projMatrix.writeToColumnMajorMatrix(glmatrix); // send projection matrix
<   safe_glUniformMatrix4fv(curSS.h_uProjMatrix, glmatrix);
---
> static void sendProjectionMatrix(Uniforms& uniforms, const Matrix4& projMatrix) {
>   uniforms.put("uProjMatrix", projMatrix);
442,449c225,231
< static ManipMode getManipMode() {
<   // if nothing is picked or the picked transform is the transfrom we are viewing from
<   if (g_currentPickedRbtNode == NULL || g_currentPickedRbtNode == g_currentCameraNode) {
<     if (g_currentCameraNode == g_skyNode && g_activeCameraFrame == WORLD_SKY)
<       return ARCBALL_ON_SKY;
<     else
<       return EGO_MOTION;
<   }
---
> static ObjId getCurrentObjId() {
>   if (g_currentPickedRbtNode == g_robot1Node)
>     return OBJECT0;
>   else if (g_currentPickedRbtNode == g_robot2Node)
>     return OBJECT1;
>   else if (g_currentPickedRbtNode)
>     return LIMB;
450a233,237
>     return SKY;
> }
> 
> static ManipMode getManipMode() {
>   if (g_currentPickedRbtNode && (g_activeEye != getCurrentObjId())) {
451a239,243
>   } else if (g_activeEye == SKY && g_activeCameraFrame == WORLD_SKY)
>     return ARCBALL_ON_SKY;
>   else {
>     return EGO_MOTION;
>   }
455c247
<   return getManipMode() != EGO_MOTION;
---
>   return getManipMode() != EGO_MOTION && (!(g_activeEye != SKY && getCurrentObjId() == SKY));
461,470c253,256
<   switch (getManipMode()) {
<   case ARCBALL_ON_PICKED:
<     return getPathAccumRbt(g_world, g_currentPickedRbtNode);
<   case ARCBALL_ON_SKY:
<     return RigTForm();
<   case EGO_MOTION:
<     return getPathAccumRbt(g_world, g_currentCameraNode);
<   default:
<     throw runtime_error("Invalid ManipMode");
<   }
---
>   if (g_currentPickedRbtNode) 
>      return getPathAccumRbt(g_world, g_currentPickedRbtNode);
>   else
>      return RigTForm();
474c260
<   RigTForm arcballEye = inv(getPathAccumRbt(g_world, g_currentCameraNode)) * getArcballRbt();
---
>   RigTForm arcballEye = inv(getNodeForEye(g_activeEye)->getRbt()) * getArcballRbt();
482c268
< static void drawArcBall(const ShaderState& curSS) {
---
> static void drawArcBall(Uniforms& uniforms) {
486,488c272,274
<   RigTForm arcballEye = inv(getPathAccumRbt(g_world, g_currentCameraNode)) * getArcballRbt();
<   Matrix4 MVM = rigTFormToMatrix(arcballEye) * Matrix4::makeScale(Cvec3(1, 1, 1) * g_arcballScale * g_arcballScreenRadius);
<   sendModelViewNormalMatrix(curSS, MVM, normalMatrix(MVM));
---
>   RigTForm arcballEye = inv(getNodeForEye(g_activeEye)->getRbt()) * getArcballRbt();
>   Matrix4 MVM = rigTFormToMatrix(arcballEye) * Matrix4::makeScale(Cvec3(1,1,1)* g_arcballScale * g_arcballScreenRadius);
>   sendModelViewNormalMatrix(uniforms, MVM, normalMatrix(MVM));
490,491c276
<   safe_glUniform3f(curSS.h_uColor, 0.27, 0.82, 0.35); // set color
<   g_sphere->draw(curSS);
---
>   g_arcballMat->draw(*g_sphere, uniforms);
493,494d277
<   // switch back to solid mode
<   glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
497c280,284
< static void drawStuff(const ShaderState& curSS, bool picking) {
---
> static void drawStuff(bool picking) {
>  
>   // Declare an empty uniforms
>   Uniforms uniforms;
>   
504c291
<   sendProjectionMatrix(curSS, projmat);
---
>   sendProjectionMatrix(uniforms, projmat);
506c293
<   const RigTForm eyeRbt = getPathAccumRbt(g_world, g_currentCameraNode);
---
>   const RigTForm eyeRbt = getPathAccumRbt(g_world, getNodeForEye(g_activeEye));
509,513c296,301
<   const Cvec3 eyeLight1 = Cvec3(invEyeRbt * Cvec4(g_light1, 1));
<   const Cvec3 eyeLight2 = Cvec3(invEyeRbt * Cvec4(g_light2, 1));
<   safe_glUniform3f(curSS.h_uLight, eyeLight1[0], eyeLight1[1], eyeLight1[2]);
<   safe_glUniform3f(curSS.h_uLight2, eyeLight2[0], eyeLight2[1], eyeLight2[2]);
< 
---
>   const Cvec3 eyeLight1 = getPathAccumRbt(g_world, g_light1Node).getTranslation();
>   const Cvec3 eyeLight2 = getPathAccumRbt(g_world, g_light2Node).getTranslation();
>   
>   uniforms.put("uLight", eyeLight1);
>   uniforms.put("uLight2", eyeLight2);
>   
515,516c303,304
<     Drawer drawer(invEyeRbt, curSS);
<     g_world->accept(drawer);
---
> 	  Drawer drawer(invEyeRbt, uniforms);
> 	  g_world->accept(drawer);
518,522c306,314
<     if (g_displayArcball && shouldUseArcball())
<       drawArcBall(curSS);
<   }
<   else {
<     Picker picker(invEyeRbt, curSS);
---
> 	  if (g_displayArcball && shouldUseArcball()) {
> 		  drawArcBall(uniforms);
> 	  }
>   } else {
> 	  Picker picker(invEyeRbt, uniforms);
> 	  
>     // set overriding material to our picking material
>     g_overridingMaterial = g_pickingMat;
>     
523a316,319
>     
>     // unset overriding material 
>     g_overridingMaterial.reset();
> 	  
525,529c321,326
<     g_currentPickedRbtNode = picker.getRbtNodeAtXY(g_mouseClickX, g_mouseClickY);
<     if (g_currentPickedRbtNode == g_groundNode)
<       g_currentPickedRbtNode = shared_ptr<SgRbtNode>(); // set to NULL
< 
<     cout << (g_currentPickedRbtNode ? "Part picked" : "No part picked") << endl;
---
> 	  g_currentPickedRbtNode = picker.getRbtNodeAtXY(g_mouseClickX, g_mouseClickY);
> 	  
> 	  if (g_currentPickedRbtNode == g_groundNode) {
> 		  g_currentPickedRbtNode = shared_ptr<SgRbtNode>();   // set to NULL
> 		  cout << "clicked on ground\n";
> 	  }
534,566c331,332
<   glUseProgram(g_shaderStates[g_activeShader]->program);
<   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
< 
<   drawStuff(*g_shaderStates[g_activeShader], false);
< 
<   glutSwapBuffers();
< 
<   checkGlErrors();
< }
< 
< static void pick() {
<   // We need to set the clear color to black, for pick rendering.
<   // so let's save the clear color
<   GLdouble clearColor[4];
<   glGetDoublev(GL_COLOR_CLEAR_VALUE, clearColor);
< 
<   glClearColor(0, 0, 0, 0);
< 
<   // using PICKING_SHADER as the shader
<   glUseProgram(g_shaderStates[PICKING_SHADER]->program);
< 
<   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
<   drawStuff(*g_shaderStates[PICKING_SHADER], true);
< 
<   // Uncomment below and comment out the glutPostRedisplay in mouse(...) call back
<   // to see result of the pick rendering pass
<   // glutSwapBuffers();
< 
<   //Now set back the clear color
<   glClearColor(clearColor[0], clearColor[1], clearColor[2], clearColor[3]);
< 
<   checkGlErrors();
< }
---
>   //glUseProgram(g_shaderStates[g_activeShader]->program);
>   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);                   // clear framebuffer color&depth
568,573c334
< bool interpolateAndDisplay(float t) {
<   if (t > g_animator.getNumKeyFrames() - 3)
<     return true;
<   g_animator.animate(t);
<   return false;
< }
---
>   drawStuff(g_pickerMode);
575,579c336,337
< static void animateTimerCallback(int ms) {
<   double t = (double)ms / g_msBetweenKeyFrames;
<   bool endReached = interpolateAndDisplay(t);
<   if (g_playingAnimation && !endReached) {
<     glutTimerFunc(1000/g_animateFramesPerSecond, animateTimerCallback, ms + 1000/g_animateFramesPerSecond);
---
>   if (!g_pickerMode) {
> 	  glutSwapBuffers();                                    // show the back buffer (where we rendered stuff)
581,586d338
<   else {
<     cerr << "Finished playing animation" << endl;
<     g_curKeyFrame = g_animator.keyFramesEnd();
<     advance(g_curKeyFrame, -2);
<     g_animator.pushKeyFrameToSg(g_curKeyFrame);
<     g_playingAnimation = false;
588,591c340
<     g_curKeyFrameNum = g_animator.getNumKeyFrames() - 2;
<     cerr << "Now at frame [" << g_curKeyFrameNum << "]" << endl;
<   }
<   display();
---
>   checkGlErrors();
599c348
<   g_arcballScreenRadius = max(1.0, min(h, w) * 0.25);
---
>   g_arcballScreenRadius = max(1.0, min(h,w) * 0.25);
603a353
> 
611a362
> 
614c365
<   const RigTForm eyeInverse = inv(getPathAccumRbt(g_world, g_currentCameraNode));
---
>   const RigTForm eyeInverse = inv(getNodeForEye(g_activeEye)->getRbt());
645c396
<       M = RigTForm(Cvec3(dx, dy, 0) * movementScale);
---
>        M = RigTForm(Cvec3(dx, dy, 0) * movementScale);
670,674c421,444
< // l = w X Y Z
< // o = l O
< // a = w A = l (Z Y X)^1 A = l A'
< // o = a (A')^-1 O
< //   => a M (A')^-1 O = l A' M (A')^-1 O
---
> static void pick() {
> 	// We need to set the clear color to black, for pick rendering.
> 	// so let's save the clear color
> 	GLdouble clearColor[4];
> 	glGetDoublev(GL_COLOR_CLEAR_VALUE, clearColor);
> 
> 	glClearColor(0, 0, 0, 0);
> 
> 	// using PICKING_SHADER as the shader
> 	//glUseProgram(g_shaderStates[PICKING_SHADER]->program);
> 
> 	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
> 	drawStuff(true);
> 
> 	// Uncomment below and comment out the glutPostRedisplay in mouse(...) call back
> 	// to see result of the pick rendering pass
>   
>   //glutSwapBuffers();
> 
> 	//Now set back the clear color
> 	glClearColor(clearColor[0], clearColor[1], clearColor[2], clearColor[3]);
> 
> 	checkGlErrors();
> }
679c449,450
< 
---
>   if (getCurrentObjId() == SKY && g_activeEye != SKY)
>     return;                  // we do not edit the eye when viewed from the objects
684a456
>   RigTForm A;
686,698c458,476
<   RigTForm A = makeMixedFrame(getArcballRbt(), getPathAccumRbt(g_world, g_currentCameraNode));
< 
<   shared_ptr<SgRbtNode> target;
<   switch (getManipMode()) {
<   case ARCBALL_ON_PICKED:
<     target = g_currentPickedRbtNode;
<     break;
<   case ARCBALL_ON_SKY:
<     target = g_skyNode;
<     break;
<   case EGO_MOTION:
<     target = g_currentCameraNode;
<     break;
---
>   if (g_activeEye == SKY) {
>     if (g_currentPickedRbtNode)
>       A = transFact(g_currentPickedRbtNode->getRbt()) * linFact(g_skyNode->getRbt()); 
>     else
>       A = linFact(g_skyNode->getRbt());
>   } else {
>      RigTForm obj = getPathAccumRbt(g_world, g_currentPickedRbtNode, 0);
>      RigTForm parent = getPathAccumRbt(g_world, g_currentPickedRbtNode, 1);
>      RigTForm frame = getPathAccumRbt(g_world, getNodeForEye(g_activeEye), 0);   
>      A = inv(parent) * transFact(obj) * linFact(frame);
>   }
>  
>   if (g_currentPickedRbtNode) {
>     RigTForm O = doMtoOwrtA(M, g_currentPickedRbtNode->getRbt(), A);
>     g_currentPickedRbtNode->setRbt(O);   
>   } else {
>     shared_ptr<SgRbtNode> eyeNode = getNodeForEye(g_activeEye);
>     RigTForm O = doMtoOwrtA(M, eyeNode->getRbt(), A); 
>     eyeNode->setRbt(O); 
701,704d478
<   A = inv(getPathAccumRbt(g_world, target, 1)) * A;
< 
<   target->setRbt(doMtoOwrtA(M, target->getRbt(), A));
< 
707c481,482
<   glutPostRedisplay();  // we always redraw if we changed the scene
---
> 
>   glutPostRedisplay();                    // we always redraw if we changed the scene
709a485
> 
724,728c500,502
<   if (g_pickingMode && button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
<     pick();
<     g_pickingMode = false;
<     cerr << "Picking mode is off" << endl;
<     glutPostRedisplay(); // request redisplay since the arcball will have moved
---
>   if (g_pickerMode && g_mouseLClickButton && !g_mouseRClickButton) {
> 	  pick();
> 	  g_pickerMode = false;
741a516,922
> static void initializeNewKeyframe() {
> 	RigTFormVector new_keyframe;
> 	
> 	if (keyframeList.empty()) {
> 		// nothing initialized yet, do so now
> 		printf("initializing first frame\n");
> 		keyframeList.push_front(new_keyframe); 
> 		g_currentKeyframe = keyframeList.begin();
> 	}
> 	else if (g_currentKeyframe == keyframeList.end()) {
> 		printf("g_currentKeyframe is undefined, it returned the end character\n");
> 		// create a new keyframe at the end of the stack
> 		keyframeList.push_back(new_keyframe);
> 		// point to last element
> 		--g_currentKeyframe;
> 	}
> 	else {
> 		printf("current key frame initialized, creating a new keyframe after it\n");
> 
> 		// iter must always be where g_currentKeyframe is
> 		++g_currentKeyframe;
> 		keyframeList.insert(g_currentKeyframe, new_keyframe);
> 		--g_currentKeyframe;
> 	}
> 	// copy scene graph rbt data to new key frame		
> 	
> 	std::vector<std::tr1::shared_ptr<SgRbtNode> > RbtNodes;
> 	dumpSgRbtNodes(g_world, RbtNodes);
> 	(*g_currentKeyframe).clear();
> 	for (int i = 0; i < RbtNodes.size(); ++i) {
> 		(*g_currentKeyframe).push_back(RbtNodes[i]->getRbt());
> 	}
> }
> 
> static void copyCurrentKeyframe() {
> 	// copy current key frame RBT data to the scene graph if current frame is defined
> 	if (keyframeList.empty()) {
> 		// no frames, do nothing
> 		printf("no key frames\n");
> 	}
> 	else {
> 		// grab pointers to current scene data with the dump function, then copy our current key frame rbts there
> 		std::vector<std::tr1::shared_ptr<SgRbtNode> > RbtNodes;
> 		dumpSgRbtNodes(g_world, RbtNodes);
> 		
> 		std::vector<shared_ptr<SgRbtNode> >::iterator iter1 = RbtNodes.begin();
> 		std::vector<RigTForm>::iterator iter2 = (*g_currentKeyframe).begin();
> 
> 		while (true) {
> 			if (iter1 != RbtNodes.end() && iter2 != (*g_currentKeyframe).end()) {
> 				(**iter1).setRbt(*iter2);
> 				++iter1;
> 				++iter2;
> 			}
> 			else {
> 				break;
> 			}
> 		}
> 	}
> }
> 
> static void deleteCurrentFrame() {
> 	// if our stack is empty, 
> 	if (keyframeList.empty()) {
> 		printf("there's no frame, nothing to delete\n");
> 	}
> 	else {
> 		printf("current frame was set, deleting now\n");
> 		
> 		list<RigTFormVector>::iterator temp = g_currentKeyframe;
> 		if (g_currentKeyframe == keyframeList.begin()) {
> 			if (++g_currentKeyframe == keyframeList.end()) {
>          g_currentKeyframe = keyframeList.end();
> 			}
> 		}
> 		else {
> 			--g_currentKeyframe;
> 		}
> 		
> 		keyframeList.erase(temp);
> 		
> 		// update scene
> 		copyCurrentKeyframe();
> 	}
> }
> 
> static RigTForm evaluateBezier(RigTForm from, RigTForm to, RigTForm d, RigTForm e, float alpha) {
> 	Cvec3 f_t = lerp(from.getTranslation(), d.getTranslation(), alpha);
> 	Cvec3 g_t = lerp(d.getTranslation(), e.getTranslation(), alpha);
> 	Cvec3 h_t = lerp(e.getTranslation(), to.getTranslation(), alpha);
> 	Cvec3 m_t = lerp(f_t, g_t, alpha);
> 	Cvec3 n_t = lerp(g_t, h_t, alpha);
> 	Cvec3 c_t = lerp(m_t, n_t, alpha);
> 
> 	Quat f_q = slerp(from.getRotation(), d.getRotation(), alpha);
> 	Quat g_q = slerp(d.getRotation(), e.getRotation(), alpha);
> 	Quat h_q = slerp(e.getRotation(), to.getRotation(), alpha);
> 	Quat m_q = slerp(f_q, g_q, alpha);
> 	Quat n_q = slerp(g_q, h_q, alpha);
> 	Quat c_q = slerp(m_q, n_q, alpha);
> 
> 	return RigTForm(c_t, c_q);
> }
> 
> static RigTForm evaluateCatmull_Rom(RigTForm prev, RigTForm from, RigTForm to, RigTForm post, float alpha) {
> 	// TODO: sanity check for quat negation if first coordinate of the part of D is negative
> 	Cvec3 BC_CvecD = (to.getTranslation() - prev.getTranslation()) * (1.0 / 6.0) + from.getTranslation();
> 	Cvec3 BC_CvecE = (post.getTranslation() - from.getTranslation()) * (-1.0 / 6.0) + to.getTranslation();
>   
>   Quat q = prev.getRotation();
>   cout << "GETPREVROTATION: " << q[0] << q[1] << q[2] << q[3] << endl;
> 
> 
> 	Quat d_pow = cn(to.getRotation() * inv(prev.getRotation()));
> 	Quat e_pow = post.getRotation() * inv(from.getRotation());
>   cout << "=======\n"; 
>   cout << d_pow[0] << d_pow[1] << d_pow[2] << d_pow[3] << endl;
>     cout << e_pow[0] << e_pow[1] << e_pow[2] << e_pow[3] << endl;
>     cout << "------\n";
> 	/*if (d_pow[0] < 0) {
> 		d_pow = d_pow.cn;
> 		printf("%i", d_pow[0]);
> 	}*/
> 
> 	//Quat BC_QuatD = (to.getRotation() * inv(prev.getRotation())).quat_pow(1 / 6) * from.getRotation();
> 	//Quat BC_QuatE = (post.getRotation() * inv(from.getRotation())).quat_pow(-1 / 6) * to.getRotation();
> 	Quat BC_QuatD = d_pow.quat_pow(1.0 / 6.0) * from.getRotation();
> 	Quat BC_QuatE = e_pow.quat_pow(-1.0 / 6.0) * to.getRotation();
> 	
> //	Quat BC_QuatD = Quat(1, 0, 0, 0);
> //	Quat BC_QuatE = Quat(1, 0, 0, 0);
> 
>     RigTForm BC_D = RigTForm(BC_CvecD, BC_QuatD);
> 	RigTForm BC_E = RigTForm(BC_CvecE, BC_QuatE);
> 
> 	return evaluateBezier(from, to, BC_D, BC_E, alpha);
> }
> 
> bool interpolateAndDisplay(float t) {
> 
>   // get alpha for slerping and lerping 
>   float alpha = t - floor(t);
> 
>   int keyframe = floor(t);
>   
>   if (keyframe == 0) {
> 	  g_previousPlayingFromKeyframe = keyframeList.begin();
> 
> 	  // for Catmull Rom we begin interpolating from the second keyframe
> 	  g_currentPlayingFromKeyframe = keyframeList.begin();
> 	  ++g_currentPlayingFromKeyframe;
> 
> 	  g_currentPlayingToKeyframe = g_currentPlayingFromKeyframe;
> 	  ++g_currentPlayingToKeyframe;
> 
> 	  // for Catmull Rom we keep track of the keyframe after the frame we're playing to
> 	  g_followingPlayingToKeyframe = g_currentPlayingToKeyframe;
> 	  ++g_followingPlayingToKeyframe;
>   }
> 
>   // interpolate between the next pair of keyframes if we have passed them
>   if (keyframe != g_mostRecentPlayedKeyframe) {
> 	  g_mostRecentPlayedKeyframe = keyframe;
> 	  ++g_previousPlayingFromKeyframe;
> 	  ++g_currentPlayingFromKeyframe;
> 	  ++g_currentPlayingToKeyframe;
> 	  ++g_followingPlayingToKeyframe;
>   }
> 
>   // stop animation when we reach the second to last keyframe
>   if (g_followingPlayingToKeyframe == keyframeList.end()) { 
> 	  g_isPlayingAnimation = false;
> 	  g_mostRecentPlayedKeyframe = 0;
> 	  return true;
>   }
> 
>   // make a list for interpolations
>   vector<RigTForm> interpolations;
> 
>   //interpolate between the two
>   int size = (*g_currentPlayingFromKeyframe).size();
>   for (int i = 0; i < size; i++) {
>     // Catmull Rom Interpolation
> 	  RigTForm interpolation = evaluateCatmull_Rom((*g_previousPlayingFromKeyframe)[i], 
> 		  (*g_currentPlayingFromKeyframe)[i], (*g_currentPlayingToKeyframe)[i], (*g_followingPlayingToKeyframe)[i], alpha);
> 	  interpolations.push_back(interpolation);
>   }
>   /*
>   int keyframe = floor(t);
>   if (keyframe == 0) {
>     g_currentPlayingFromKeyframe = keyframeList.begin();
>     g_currentPlayingToKeyframe = g_currentPlayingFromKeyframe;
>     ++g_currentPlayingToKeyframe;
>   }
>  
>   // interpolate between the next pair of keyframes if we have passed them
>   if (keyframe != g_mostRecentPlayedKeyframe) {
>     g_mostRecentPlayedKeyframe = keyframe;
>     ++g_currentPlayingFromKeyframe; 
>     ++g_currentPlayingToKeyframe; 
>   }
> 
>   // stop animation when we are done
>   if (g_currentPlayingToKeyframe == keyframeList.end()) {
>     g_isPlayingAnimation = false;
>     g_mostRecentPlayedKeyframe = 0;
>     return true;
>   }
>   
>   // make a list for interpolations
>   vector<RigTForm> interpolations;
> 
>   //interpolate between the two
>   int size = (*g_currentPlayingFromKeyframe).size(); 
>   for (int i = 0; i < size; i++) { 
>     // lerp translation and slerp rotation
>     RigTForm interpolation = 
>       RigTForm(
>         lerp(
>           (*g_currentPlayingFromKeyframe)[i].getTranslation(),
>           (*g_currentPlayingToKeyframe)[i].getTranslation(),
>           alpha
>         ),
>         slerp(
>           (*g_currentPlayingFromKeyframe)[i].getRotation(),
>           (*g_currentPlayingToKeyframe)[i].getRotation(),
>           alpha
>         )
>       );
>     interpolations.push_back(interpolation); 
>   }*/
>  
>   // get current nodes in scene
> 	std::vector<std::tr1::shared_ptr<SgRbtNode> > current_nodes;
> 	dumpSgRbtNodes(g_world, current_nodes);
> 
>   // apply the transformation
>   for (int i = 0; i < current_nodes.size(); i++) {
>     current_nodes[i]->setRbt(interpolations[i]);  
>   }
>  
>   // redraw scene 
>   glutPostRedisplay();
> 
>   // the animation is not done, so return false
>   return false;
> } 
> 
> static void animateTimerCallback(int ms) {
>   float t = (float) ms / (float) g_msBetweenKeyFrames;
> 
>   bool endReached = interpolateAndDisplay(t);
>   if (!endReached) {
>           glutTimerFunc(1000/g_animateFramesPerSecond,
>                         animateTimerCallback,
>                         ms + 1000/g_animateFramesPerSecond
>           );
>   } else {
>     g_isPlayingAnimation = false; 
>   }
> }
> 
> static void controlAnimation() {
>   if (!g_isPlayingAnimation) {
>     // check that there are 4 keyframes
>     if (keyframeList.size() < 4) {
>       cout << "Unable to start animation -- you only have "
>            << keyframeList.size()
>            << " keyframes, but we need at least 4!\n";
>       return;
>     }
>     // toggle global flag for animation play
>     g_isPlayingAnimation = true;
> 
>     // call animateTimerCallback(0)
>     animateTimerCallback(0);
> 
>   } else {
>     // animation is alrady playing, so stop 
>     g_isPlayingAnimation = false; 
>   }
> }
> 
> static void shiftKeyframe(KeyframeId i) {
> 	if (i == ADVANCE && g_currentKeyframe != keyframeList.end()) {
> 		++g_currentKeyframe;
> 		if (g_currentKeyframe == keyframeList.end()) {
> 			printf("Cannot move forward, reached end of stack\n");
> 			--g_currentKeyframe;
> 			return;
> 		}
> 	}
> 	else if (i == RETREAT) {
> 		if (g_currentKeyframe == keyframeList.begin()) {
> 			printf("Cannot move backwards, at beginning of stack\n");
> 			return;
> 		}
> 		--g_currentKeyframe;
> 	}
> 	copyCurrentKeyframe();
> }
> 
> static void updateScene() {
> 	// if empty stack, make a new keyframe
> 	if (keyframeList.empty()) {
> 		initializeNewKeyframe();
> 	}
> 	// otherwise, copy the scene graph rbt to current frame 
> 	else {
> 		std::vector<std::tr1::shared_ptr<SgRbtNode> > RbtNodes;
> 		dumpSgRbtNodes(g_world, RbtNodes);
> 		(*g_currentKeyframe).clear();
> 		for (int i = 0; i < RbtNodes.size(); ++i) {
> 			(*g_currentKeyframe).push_back(RbtNodes[i]->getRbt());
> 		}
> 	}
> }
> 
> static void readFrameDataFromFile() {
> 	ifstream myfile; 
> 	string line; 
> 	myfile.open("animation.txt");
> 	getline(myfile, line);
> 	int num_frames = atoi(line.c_str());
> 	getline(myfile, line);
> 	int num_nodes = atoi(line.c_str());
> 	printf("%d %d\n", num_frames, num_nodes);
> 	
> 	keyframeList.clear();
> 	for (int i = 0; i < num_frames; i++) {
> 		RigTFormVector frame; 
> 		for (int j = 0; j < num_nodes; j++) {
> 			RigTForm node; 
> 			Cvec3 translation;
> 			Quat rotation;
> 			std::string line;
> 			std::getline(myfile, line);
> 			std::stringstream stream(line);
> 			for (int k = 0; k < 3; k++) {
> 				double cvec_double;
> 				stream >> cvec_double;
> 				translation[k] = cvec_double;
> 			}
> 			for (int l = 0; l < 4; l++) {
> 				double quat_double;
> 				stream >> quat_double;
> 				rotation[l] = quat_double;
> 			}
> 			printf("\n");
> 			node.setTranslation(translation); 
> 			node.setRotation(rotation);
> 			frame.push_back(node);
> 		}
> 		keyframeList.push_back(frame);
> 	}
> 	g_currentKeyframe = keyframeList.begin();
> 
> 	copyCurrentKeyframe();
> 
> 	myfile.close();
> }
> 
> static void writeFrameDataToFile() {
> 	if (!keyframeList.empty()) {
> 		ofstream myfile;
> 		myfile.open("animation.txt");
> 		// myfile << "Writing this to a file.\n";
> 		// myfile << "testing, testing \n";
> 		string output = "";
> 		int num_frames = keyframeList.size();
> 		list<RigTFormVector>::iterator iter = keyframeList.begin();
> 		int num_nodes = (*iter).size();
>     
>     std::ostringstream s;
>     s << num_frames
>       << "\n"
>       << num_nodes
>       << "\n";
>     output.append(s.str());
> 		for (iter; iter != keyframeList.end(); ++iter) {
> 			RigTFormVector scene = (*iter);
> 			for (int i = 0; i < scene.size(); ++i) {
> 				RigTForm rigTForm = scene[i];
> 				Cvec3 translation = rigTForm.getTranslation();
> 				Quat rotation = rigTForm.getRotation();
> 				for (int j = 0; j < 3; j++) {
> 				  s.str(""); 
>           s.clear(); 
>           s << translation[j] << " ";
>           output.append(s.str());
> 				}
> 				for (int k = 0; k < 4; k++) {
> 				  s.str(""); 
>           s.clear(); 
>           s << rotation[k] << " ";
>           output.append(s.str());
> 				}
> 				output.append("\n");
> 			}
> 		}
> 		myfile << output;
> 
> 		myfile.close();
> 		printf("Writing to animation.txt\n");
> 	}
> }
> 
744,746d924
<   case ' ':
<     g_spaceDown = true;
<     break;
754c932
<     << "p\t\tUse mouse to pick a part to edit\n"
---
>     << "p\t\tEnter picker mode\n"
757,766c935,947
<     << "a\t\tToggle display arcball\n"
<     << "w\t\tWrite animation to animation.txt\n"
<     << "i\t\tRead animation from animation.txt\n"
<     << "c\t\tCopy frame to scene\n"
<     << "u\t\tCopy sceneto frame\n"
<     << "n\t\tCreate new frame after current frame and copy scene to it\n"
<     << "d\t\tDelete frame\n"
<     << ">\t\tGo to next frame\n"
<     << "<\t\tGo to prev. frame\n"
<     << "y\t\tPlay/Stop animation\n"
---
>     << "drag right or middle mouse to translate\n\n"
>     << "\t---ANIMATION CONTROLS---\n"
>     << "n\t\tMake new keyframe\n"
>     << "c\t\tCopy current keyframe\n"
>     << "u\t\tUpdate keyframe\n"
>     << "d\t\tDelete current keyframe\n" 
>     << ">\t\tGo to next keyframe\n"
>     << "<\t\tGo to previous keyframe\n"
>     << "w\t\tWrite keyframes to file\n"
>     << "i\t\tRead keyframes from file\n"
>     << "y\t\tPlay/Stop playing keyframes\n"
>     << "+\t\tSpeed up animation\n"
>     << "-\t\tSlow down animation\n"
773,775c954,956
<   case 'f':
<     g_activeShader = (g_activeShader + 1) % g_numRegularShaders;
<     break;
---
>   /*case 'f':
>     g_activeShader ^= 1;
>     break;*/
777,789c958,959
<   {
<     shared_ptr<SgRbtNode> viewers[] = {g_skyNode, g_robot1Node, g_robot2Node};
<     for (int i = 0; i < 3; ++i) {
<       if (g_currentCameraNode == viewers[i]) {
<         g_currentCameraNode = viewers[(i+1)%3];
<         break;
<       }
<     }
<   }
<   break;
<   case 'p':
<     g_pickingMode = !g_pickingMode;
<     cerr << "Picking mode is " << (g_pickingMode ? "on" : "off") << endl;
---
>     g_activeEye = ObjId((g_activeEye+1) % 3);
>     cerr << "Active eye is " << g_objNames[g_activeEye] << endl;
795,802c965,971
<   case 'a':
<     g_displayArcball = !g_displayArcball;
<     break;
<   case 'u':
<     if (g_playingAnimation) {
<       cerr << "Cannot operate when playing animation" << endl;
<       break;
<     }
---
>   case 'p':
> 	  g_pickerMode = !g_pickerMode;
> 	  cout << "Picking mode: " << g_pickerMode << "\n";
> 	  break;
>   case ' ':
> 	  g_spaceDown = true;
> 	  break;
804,822c973
<     if (g_curKeyFrame == g_animator.keyFramesEnd()) { // only possible when frame list is empty
<       cerr << "Create new frame [0]."  << endl;
<       g_curKeyFrame = g_animator.insertEmptyKeyFrameAfter(g_animator.keyFramesBegin());
<       g_curKeyFrameNum = 0;
<     }
<     cerr << "Copying scene graph to current frame [" << g_curKeyFrameNum << "]" << endl;
<     g_animator.pullKeyFrameFromSg(g_curKeyFrame);
<     break;
<   case 'n':
<     if (g_playingAnimation) {
<       cerr << "Cannot operate when playing animation" << endl;
<       break;
<     }
<     if (g_animator.getNumKeyFrames() != 0)
<       ++g_curKeyFrameNum;
<     g_curKeyFrame = g_animator.insertEmptyKeyFrameAfter(g_curKeyFrame);
<     g_animator.pullKeyFrameFromSg(g_curKeyFrame);
<     cerr << "Create new frame [" << g_curKeyFrameNum << "]" << endl;
<     break;
---
>   // keyframe animation
824,863c975,981
<     if (g_playingAnimation) {
<       cerr << "Cannot operate when playing animation" << endl;
<       break;
<     }
<     if (g_curKeyFrame != g_animator.keyFramesEnd()) {
<       cerr << "Loading current key frame [" << g_curKeyFrameNum << "] to scene graph" << endl;
<       g_animator.pushKeyFrameToSg(g_curKeyFrame);
<     }
<     else {
<       cerr << "No key frame defined" << endl;
<     }
<     break;
<   case 'd':
<     if (g_playingAnimation) {
<       cerr << "Cannot operate when playing animation" << endl;
<       break;
<     }
<     if (g_curKeyFrame != g_animator.keyFramesEnd()) {
<       Animator::KeyFrameIter newCurKeyFrame = g_curKeyFrame;
<       cerr << "Deleting current frame [" << g_curKeyFrameNum << "]" << endl;;
<       if (g_curKeyFrame == g_animator.keyFramesBegin()) {
<         ++newCurKeyFrame;
<       }
<       else {
<         --newCurKeyFrame;
<         --g_curKeyFrameNum;
<       }
<       g_animator.deleteKeyFrame(g_curKeyFrame);
<       g_curKeyFrame = newCurKeyFrame;
<       if (g_curKeyFrame != g_animator.keyFramesEnd()) {
<         g_animator.pushKeyFrameToSg(g_curKeyFrame);
<         cerr << "Now at frame [" << g_curKeyFrameNum << "]" << endl;
<       }
<       else
<         cerr << "No frames defined" << endl;
<     }
<     else {
<       cerr << "Frame list is now EMPTY" << endl;
<     }
<     break;
---
> 	  printf("c was pressed\n");
> 	  copyCurrentKeyframe();
> 	  break;
>   case 'u':
> 	  printf("u was pressed\n");
> 	  updateScene();
> 	  break;
865,878c983,985
<     if (g_playingAnimation) {
<       cerr << "Cannot operate when playing animation" << endl;
<       break;
<     }
<     if (g_curKeyFrame != g_animator.keyFramesEnd()) {
<       if (++g_curKeyFrame == g_animator.keyFramesEnd())
<         --g_curKeyFrame;
<       else {
<         ++g_curKeyFrameNum;
<         g_animator.pushKeyFrameToSg(g_curKeyFrame);
<         cerr << "Stepped forward to frame [" << g_curKeyFrameNum <<"]" << endl;
<       }
<     }
<     break;
---
> 	  printf("> was pressed\n");
> 	  shiftKeyframe(ADVANCE);
> 	  break;
880,894c987,997
<     if (g_playingAnimation) {
<       cerr << "Cannot operate when playing animation" << endl;
<       break;
<     }
<     if (g_curKeyFrame != g_animator.keyFramesBegin()) {
<       --g_curKeyFrame;
<       --g_curKeyFrameNum;
<       g_animator.pushKeyFrameToSg(g_curKeyFrame);
<       cerr << "Stepped backward to frame [" << g_curKeyFrameNum << "]" << endl;
<     }
<     break;
<   case 'w':
<     cerr << "Writing animation to animation.txt\n";
<     g_animator.saveAnimation("animation.txt");
<     break;
---
> 	  printf("< was pressed\n");
> 	  shiftKeyframe(RETREAT);
> 	break;
>   case 'd':
> 	  printf("d was pressed\n");
> 	  deleteCurrentFrame();
> 	  break;
>   case 'n':
> 	  printf("n was pressed\n");
> 	  initializeNewKeyframe();
> 	  break;
896,912c999,1008
<     if (g_playingAnimation) {
<       cerr << "Cannot operate when playing animation" << endl;
<       break;
<     }
<     cerr << "Reading animation from animation.txt\n";
<     g_animator.loadAnimation("animation.txt");
<     g_curKeyFrame = g_animator.keyFramesBegin();
<     cerr << g_animator.getNumKeyFrames() << " frames read.\n";
<     if (g_curKeyFrame != g_animator.keyFramesEnd()) {
<       g_animator.pushKeyFrameToSg(g_curKeyFrame);
<       cerr << "Now at frame [0]" << endl;
<     }
<     g_curKeyFrameNum = 0;
<     break;
<   case '-':
<     g_msBetweenKeyFrames = min(g_msBetweenKeyFrames + 100, 10000);
<     cerr << g_msBetweenKeyFrames << " ms between keyframes.\n";
---
> 	printf("i was pressed\n");
> 	readFrameDataFromFile();
> 	break;
>   case 'w':
> 	printf("w was pressed\n");
> 	writeFrameDataToFile();
> 	break;
>   case 'y':
>     printf("y was pressed\n");
>     controlAnimation();
915,932c1011,1018
<     g_msBetweenKeyFrames = max(g_msBetweenKeyFrames - 100, 100);
<     cerr << g_msBetweenKeyFrames << " ms between keyframes.\n";
<     break;
<   case 'y':
<     if (!g_playingAnimation) {
<       if (g_animator.getNumKeyFrames() < 4) {
<         cerr << " Cannot play animation with less than 4 keyframes." << endl;
<       }
<       else {
<         g_playingAnimation = true;
<         cerr << "Playing animation... "<< endl;
<         animateTimerCallback(0);
<       }
<     }
<     else {
<       cerr << "Stopping animation... " << endl;
<       g_playingAnimation = false;
<     }
---
>     printf("+ was pressed\n");
>     g_msBetweenKeyFrames = max(100, g_msBetweenKeyFrames - 100);
>     cout << "The new speed is: " << g_msBetweenKeyFrames << "ms between keyframes\n"; 
>    break;
>   case '-':
>     printf("- was pressed\n"); 
>     g_msBetweenKeyFrames = min(1000, g_msBetweenKeyFrames + 100);
>     cout << "The new speed is: " << g_msBetweenKeyFrames << "ms between keyframes\n"; 
935,939d1020
< 
<   // Sanity check that our g_curKeyFrameNum is in sync with the g_curKeyFrame
<   if (g_animator.getNumKeyFrames() > 0)
<     assert(g_animator.getNthKeyFrame(g_curKeyFrameNum) == g_curKeyFrame);
< 
951,952c1032,1033
<   glutCreateWindow("Assignment 5");                       // title the window
< 
---
>   glutCreateWindow("Assignment 4");                       // title the window
>  
961a1043
>  
977,984c1059,1087
< static void initShaders() {
<   g_shaderStates.resize(g_numShaders);
<   for (int i = 0; i < g_numShaders; ++i) {
<     if (g_Gl2Compatible)
<       g_shaderStates[i].reset(new ShaderState(g_shaderFilesGl2[i][0], g_shaderFilesGl2[i][1]));
<     else
<       g_shaderStates[i].reset(new ShaderState(g_shaderFiles[i][0], g_shaderFiles[i][1]));
<   }
---
> static void initMaterials() {
>   // Create some prototype materials
>   Material diffuse("./shaders/basic-gl3.vshader", "./shaders/diffuse-gl3.fshader");
>   Material solid("./shaders/basic-gl3.vshader", "./shaders/solid-gl3.fshader");
> 
>   // copy diffuse prototype and set red color
>   g_redDiffuseMat.reset(new Material(diffuse));
>   g_redDiffuseMat->getUniforms().put("uColor", Cvec3f(1, 0, 0));
> 
>   // copy diffuse prototype and set blue color
>   g_blueDiffuseMat.reset(new Material(diffuse));
>   g_blueDiffuseMat->getUniforms().put("uColor", Cvec3f(0, 0, 1));
> 
>   // normal mapping material
>   g_bumpFloorMat.reset(new Material("./shaders/normal-gl3.vshader", "./shaders/normal-gl3.fshader"));
>   g_bumpFloorMat->getUniforms().put("uTexColor", shared_ptr<ImageTexture>(new ImageTexture("Fieldstone.ppm", true)));
>   g_bumpFloorMat->getUniforms().put("uTexNormal", shared_ptr<ImageTexture>(new ImageTexture("FieldstoneNormal.ppm", false)));
> 
>   // copy solid prototype, and set to wireframed rendering
>   g_arcballMat.reset(new Material(solid));
>   g_arcballMat->getUniforms().put("uColor", Cvec3f(0.27f, 0.82f, 0.35f));
>   g_arcballMat->getRenderStates().polygonMode(GL_FRONT_AND_BACK, GL_LINE);
> 
>   // copy solid prototype, and set to color white
>   g_lightMat.reset(new Material(solid));
>   g_lightMat->getUniforms().put("uColor", Cvec3f(1, 1, 1));
> 
>   // pick shader
>   g_pickingMat.reset(new Material("./shaders/basic-gl3.vshader", "./shaders/pick-gl3.fshader"));
991d1093
<   initRobots();
994,1061c1096,1165
< static void constructRobot(shared_ptr<SgTransformNode> base, const Cvec3& color) {
<   const double ARM_LEN = 0.7,
<                ARM_THICK = 0.25,
<                LEG_LEN = 1,
<                LEG_THICK = 0.25,
<                TORSO_LEN = 1.5,
<                TORSO_THICK = 0.25,
<                TORSO_WIDTH = 1,
<                HEAD_SIZE = 0.7;
<   const int NUM_JOINTS = 10,
<             NUM_SHAPES = 10;
< 
<   struct JointDesc {
<     int parent;
<     float x, y, z;
<   };
< 
<   JointDesc jointDesc[NUM_JOINTS] = {
<     {-1}, // torso
<     {0,  TORSO_WIDTH/2, TORSO_LEN/2, 0}, // upper right arm
<     {0, -TORSO_WIDTH/2, TORSO_LEN/2, 0}, // upper left arm
<     {1,  ARM_LEN, 0, 0}, // lower right arm
<     {2, -ARM_LEN, 0, 0}, // lower left arm
<     {0, TORSO_WIDTH/2-LEG_THICK/2, -TORSO_LEN/2, 0}, // upper right leg
<     {0, -TORSO_WIDTH/2+LEG_THICK/2, -TORSO_LEN/2, 0}, // upper left leg
<     {5, 0, -LEG_LEN, 0}, // lower right leg
<     {6, 0, -LEG_LEN, 0}, // lower left
<     {0, 0, TORSO_LEN/2, 0} // head
<   };
< 
<   struct ShapeDesc {
<     int parentJointId;
<     float x, y, z, sx, sy, sz;
<     shared_ptr<Geometry> geometry;
<   };
< 
<   ShapeDesc shapeDesc[NUM_SHAPES] = {
<     {0, 0,         0, 0, TORSO_WIDTH, TORSO_LEN, TORSO_THICK, g_cube}, // torso
<     {1, ARM_LEN/2, 0, 0, ARM_LEN/2, ARM_THICK/2, ARM_THICK/2, g_sphere}, // upper right arm
<     {2, -ARM_LEN/2, 0, 0, ARM_LEN/2, ARM_THICK/2, ARM_THICK/2, g_sphere}, // upper left arm
<     {3, ARM_LEN/2, 0, 0, ARM_LEN, ARM_THICK, ARM_THICK, g_cube}, // lower right arm
<     {4, -ARM_LEN/2, 0, 0, ARM_LEN, ARM_THICK, ARM_THICK, g_cube}, // lower left arm
<     {5, 0, -LEG_LEN/2, 0, LEG_THICK/2, LEG_LEN/2, LEG_THICK/2, g_sphere}, // upper right leg
<     {6, 0, -LEG_LEN/2, 0, LEG_THICK/2, LEG_LEN/2, LEG_THICK/2, g_sphere}, // upper left leg
<     {7, 0, -LEG_LEN/2, 0, LEG_THICK, LEG_LEN, LEG_THICK, g_cube}, // lower right leg
<     {8, 0, -LEG_LEN/2, 0, LEG_THICK, LEG_LEN, LEG_THICK, g_cube}, // lower left leg
<     {9, 0, HEAD_SIZE/2 * 1.5, 0, HEAD_SIZE/2, HEAD_SIZE/2, HEAD_SIZE/2, g_sphere}, // head
<   };
< 
<   shared_ptr<SgTransformNode> jointNodes[NUM_JOINTS];
< 
<   for (int i = 0; i < NUM_JOINTS; ++i) {
<     if (jointDesc[i].parent == -1)
<       jointNodes[i] = base;
<     else {
<       jointNodes[i].reset(new SgRbtNode(RigTForm(Cvec3(jointDesc[i].x, jointDesc[i].y, jointDesc[i].z))));
<       jointNodes[jointDesc[i].parent]->addChild(jointNodes[i]);
<     }
<   }
<   for (int i = 0; i < NUM_SHAPES; ++i) {
<     shared_ptr<MyShapeNode> shape(
<       new MyShapeNode(shapeDesc[i].geometry,
<                       color,
<                       Cvec3(shapeDesc[i].x, shapeDesc[i].y, shapeDesc[i].z),
<                       Cvec3(0, 0, 0),
<                       Cvec3(shapeDesc[i].sx, shapeDesc[i].sy, shapeDesc[i].sz)));
<     jointNodes[shapeDesc[i].parentJointId]->addChild(shape);
<   }
---
> static void constructRobot(shared_ptr<SgTransformNode> base, shared_ptr<Material> material) {
> 
> 	const double ARM_LEN = 0.7,
> 		ARM_THICK = 0.25,
> 		TORSO_LEN = 1.5,
> 		TORSO_THICK = 0.25,
> 		TORSO_WIDTH = 1,
> 		LEG_LEN = 0.7,
> 		LEG_THICK = 0.25,
> 		HEAD_WIDTH = 0.5,
> 		HEAD_HEIGHT = 0.5;
> 	const int NUM_JOINTS = 10,
> 		NUM_SHAPES = 10;
> 
> 	struct JointDesc {
> 		int parent;
> 		float x, y, z;
> 	};
> 
> 	JointDesc jointDesc[NUM_JOINTS] = {
> 			{ -1 }, // torso
> 			{ 0, TORSO_WIDTH / 2, TORSO_LEN / 2, 0 }, // upper right arm
> 			{ 1, ARM_LEN, 0, 0 }, // lower right arm
> 			{ 0, TORSO_WIDTH / 2, -TORSO_LEN / 2, 0 }, // upper right leg
> 			{ 3, 0, -LEG_LEN, 0 }, // lower right leg
> 			{ 0, -TORSO_WIDTH / 2, TORSO_LEN / 2, 0 }, // upper left arm
> 			{ 5, -ARM_LEN, 0, 0 }, // lower left arm
> 			{ 0, -TORSO_WIDTH / 2, -TORSO_LEN / 2, 0 }, // upper left leg
> 			{ 7, 0, -LEG_LEN, 0 }, // lower left leg
> 			{ 0, 0, TORSO_LEN / 2, 0 }, // head
> 	};
> 
> 	struct ShapeDesc {
> 		int parentJointId;
> 		float x, y, z, sx, sy, sz;
> 		shared_ptr<Geometry> geometry;
> 	};
> 
> 	ShapeDesc shapeDesc[NUM_SHAPES] = {
> 			{ 0, 0, 0, 0, TORSO_WIDTH, TORSO_LEN, TORSO_THICK, g_cube }, // torso
> 			{ 1, ARM_LEN / 2, 0, 0, ARM_LEN, ARM_THICK, ARM_THICK, g_cube }, // upper right arm
> 			{ 2, ARM_LEN / 2, 0, 0, ARM_LEN, ARM_THICK, ARM_THICK, g_cube }, // lower right arm
> 			{ 3, 0, -LEG_LEN / 2, 0, LEG_THICK, LEG_LEN, LEG_THICK, g_cube }, // upper right leg
> 			{ 4, 0, -LEG_LEN / 2, 0, LEG_THICK, LEG_LEN, LEG_THICK, g_cube }, // lower right leg
> 			{ 5, - ARM_LEN / 2, 0, 0, ARM_LEN, ARM_THICK, ARM_THICK, g_cube }, // upper left arm
> 			{ 6, - ARM_LEN / 2, 0, 0, ARM_LEN, ARM_THICK, ARM_THICK, g_cube }, // lower left arm
> 			{ 7, 0, -LEG_LEN / 2, 0, LEG_THICK, LEG_LEN, LEG_THICK, g_cube }, // upper left leg
> 			{ 8, 0, -LEG_LEN / 2, 0, LEG_THICK, LEG_LEN, LEG_THICK, g_cube }, // lower left leg
> 			{ 9, 0, HEAD_HEIGHT, 0, HEAD_WIDTH, HEAD_HEIGHT, HEAD_WIDTH, g_sphere }, // head
> 	};
> 
> 	shared_ptr<SgTransformNode> jointNodes[NUM_JOINTS];
> 
> 	for (int i = 0; i < NUM_JOINTS; ++i) {
> 		if (jointDesc[i].parent == -1)
> 			jointNodes[i] = base;
> 		else {
> 			jointNodes[i].reset(new SgRbtNode(RigTForm(Cvec3(jointDesc[i].x, jointDesc[i].y, jointDesc[i].z))));
> 			jointNodes[jointDesc[i].parent]->addChild(jointNodes[i]);
> 		}
> 	}
> 	for (int i = 0; i < NUM_SHAPES; ++i) {
> 		shared_ptr<MyShapeNode> shape(
> 			new MyShapeNode(shapeDesc[i].geometry,
> 			material,
> 			Cvec3(shapeDesc[i].x, shapeDesc[i].y, shapeDesc[i].z),
> 			Cvec3(0, 0, 0),
> 			Cvec3(shapeDesc[i].sx, shapeDesc[i].sy, shapeDesc[i].sz)));
> 		jointNodes[shapeDesc[i].parentJointId]->addChild(shape);
> 	}
1065c1169
<   g_world.reset(new SgRootNode());
---
> 	g_world.reset(new SgRootNode());
1067c1171
<   g_skyNode.reset(new SgRbtNode(RigTForm(Cvec3(0.0, 0.25, 4.0))));
---
> 	g_skyNode.reset(new SgRbtNode(RigTForm(Cvec3(0.0, 0.25, 4.0))));
1068a1173,1175
> 	g_groundNode.reset(new SgRbtNode());
> 	g_groundNode->addChild(shared_ptr<MyShapeNode>(
> 		new MyShapeNode(g_ground, g_bumpFloorMat, Cvec3(0, g_groundY, 0))));
1070,1072c1177,1178
<   g_groundNode.reset(new SgRbtNode());
<   g_groundNode->addChild(shared_ptr<MyShapeNode>(
<                            new MyShapeNode(g_ground, Cvec3(0.1, 0.95, 0.1))));
---
> 	g_robot1Node.reset(new SgRbtNode(RigTForm(Cvec3(-2, 1, 0))));
> 	g_robot2Node.reset(new SgRbtNode(RigTForm(Cvec3(2, 1, 0))));
1074,1075c1180,1181
<   g_robot1Node.reset(new SgRbtNode(RigTForm(Cvec3(-2, 1, 0))));
<   g_robot2Node.reset(new SgRbtNode(RigTForm(Cvec3(2, 1, 0))));
---
> 	constructRobot(g_robot1Node, g_redDiffuseMat); // a Red robot
> 	constructRobot(g_robot2Node, g_blueDiffuseMat); // a Blue robot
1077,1078c1183,1185
<   constructRobot(g_robot1Node, Cvec3(1, 0, 0)); // a Red robot
<   constructRobot(g_robot2Node, Cvec3(0, 0, 1)); // a Blue robot
---
>   /* lights */
>   Cvec3 light1_pos = Cvec3(4.5, 3.0, 5.0);
>   Cvec3 light2_pos = Cvec3(-4.5, 0, -5.0);
1080,1083c1187,1189
<   g_world->addChild(g_skyNode);
<   g_world->addChild(g_groundNode);
<   g_world->addChild(g_robot1Node);
<   g_world->addChild(g_robot2Node);
---
>   g_light1Node.reset(new SgRbtNode(RigTForm(light1_pos)));
>   g_light1Node->addChild(shared_ptr<MyShapeNode>(
>     new MyShapeNode(g_sphere, g_lightMat, Cvec3(0, 0, 0))));
1085,1086c1191,1196
<   g_currentCameraNode = g_skyNode;
< }
---
>   g_light2Node.reset(new SgRbtNode(RigTForm(light2_pos)));
>   g_light2Node->addChild(shared_ptr<MyShapeNode>(
>     new MyShapeNode(g_sphere, g_lightMat, Cvec3(0, 0, 0))));
> 
> 	g_world->addChild(g_skyNode);
> 	g_world->addChild(g_groundNode);
1088,1090c1198,1203
< static void initAnimation() {
<   g_animator.attachSceneGraph(g_world);
<   g_curKeyFrame = g_animator.keyFramesBegin();
---
>   g_world->addChild(g_light1Node);
>   g_world->addChild(g_light2Node);
> 
>   g_world->addChild(g_robot1Node);
> 	g_world->addChild(g_robot2Node);
>   
1097c1210
<     // on Mac, we shouldn't use GLEW.
---
> 	// on Mac, we shouldn't use GLEW.
1113c1226,1227
<     initShaders();
---
> //    initShaders();
>     initMaterials(); 
1115,1116c1229
<     initScene();
<     initAnimation();
---
> 	  initScene();
