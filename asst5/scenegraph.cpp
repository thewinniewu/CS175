#include <algorithm>
#include <cassert>

#include "scenegraph.h"

using namespace std;
using namespace std::tr1;

bool SgTransformNode::accept(SgNodeVisitor& visitor) {
  if (!visitor.visit(*this))
    return false;
  for (int i = 0, n = children_.size(); i < n; ++i) {
    if (!children_[i]->accept(visitor))
      return false;
  }
  return visitor.postVisit(*this);
}

void SgTransformNode::addChild(shared_ptr<SgNode> child) {
  children_.push_back(child);
}

void SgTransformNode::removeChild(shared_ptr<SgNode> child) {
  children_.erase(find(children_.begin(), children_.end(), child));
}

bool SgShapeNode::accept(SgNodeVisitor& visitor) {
  if (!visitor.visit(*this))
    return false;
  return visitor.postVisit(*this);
}

class RbtAccumVisitor : public SgNodeVisitor {
protected:
  vector<RigTForm> rbtStack_;
  SgTransformNode& target_;
  bool found_;
public:
  RbtAccumVisitor(SgTransformNode& target)
    : target_(target)
    , found_(false) {}

  const RigTForm getAccumulatedRbt(int offsetFromStackTop = 0) {
    // TODO
	  int index = rbtStack_.size() - offsetFromStackTop;
	  //printf("%i", index);
	  RigTForm init = RigTForm();
	  for (int i = 0; i < index; i++) {
		//  printf("%i", i);
		//  printRigTForm(rbtStack_[i]);
		  init = init * rbtStack_[i];
	  }
	  return init;
  }

  virtual bool visit(SgTransformNode& node) {
    // TODO
//	 printRigTForm(node.getRbt());
	rbtStack_.push_back(node.getRbt());

	found_ = (target_ != node);
	return found_;
	/*
	if (target_ == node) {
		return false;
	}
	else {
		return true;
	}*/
  }

  virtual bool postVisit(SgTransformNode& node) {
    // TODO
	  rbtStack_.pop_back();
	  return true;
  }
};

RigTForm getPathAccumRbt(
  shared_ptr<SgTransformNode> source,
  shared_ptr<SgTransformNode> destination,
  int offsetFromDestination) {
 
  // Ensure source and destination aren't null ptrs
  assert(source);
  assert(destination);

  RbtAccumVisitor accum(*destination);
  source->accept(accum);
  return accum.getAccumulatedRbt(offsetFromDestination);
}
