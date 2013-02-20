/*
 * slp.cpp
 *
 *  Created on: Feb 20, 2013
 *      Author: dpantele
 */


#include "slp.h"

namespace crag {
namespace slp {

const Vertex Vertex::Null;
const LongInteger Vertex::LongZero;

Vertex internal::BasicNonterminalVertex::negate() const {
  return NonterminalVertex(::std::make_shared<BasicNonterminalVertex>(
      ::std::shared_ptr<NonterminalVertexNodeData>(node_data_ptr_),
      !negate_node_
  ));
}

const ::std::hash<std::shared_ptr<internal::NonterminalVertexNodeData>> internal::BasicNonterminalVertex::ptr_hash;

const inspector::InspectorTask inspector::InspectorTask::DO_NOTHING = inspector::InspectorTask();

}
}


