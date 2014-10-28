#ifndef KEYFRAME_H
#define KEYFRAME_H

#include <vector>

#include "sgutils.h"
#include "scenegraph.h"

typedef std::vector<std::tr1::shared_ptr<SgRbtNode> > SgRbtNodes;

class KeyframeList : public SgNodeVisitor {
protected:
	shared_ptr<SgRbtNodes> g_currentKeyframe;
	list<SgRbtNodes> keyframeList_;

public:

}


inline void copyCurrentKeyframe(){

}

/*
struct RbtNodesScanner : public SgNodeVisitor {
	typedef std::vector<std::tr1::shared_ptr<SgRbtNode> > SgRbtNodes;

	SgRbtNodes& nodes_;

	RbtNodesScanner(SgRbtNodes& nodes) : nodes_(nodes) {}

	virtual bool visit(SgTransformNode& node) {
		using namespace std;
		using namespace tr1;
		shared_ptr<SgRbtNode> rbtPtr = dynamic_pointer_cast<SgRbtNode>(node.shared_from_this());
		if (rbtPtr)
			nodes_.push_back(rbtPtr);
		return true;
	}
};

inline void dumpSgRbtNodes(std::tr1::shared_ptr<SgNode> root, std::vector<std::tr1::shared_ptr<SgRbtNode> >& rbtNodes) {
	RbtNodesScanner scanner(rbtNodes);
	root->accept(scanner);
}

*/
#endif