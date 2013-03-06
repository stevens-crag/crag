/*
 * slp_mapper.h
 *
 *  Created on: Mar 6, 2013
 *      Author: pmorar
 */

#pragma once
#ifndef CRAG_FREEGROUP_SLP_MAPPER_H_
#define CRAG_FREEGROUP_SLP_MAPPER_H_

#include <unordered_map>
#include "slp_vertex.h"

namespace crag {
namespace slp {



namespace mapper {

typedef typename std::unordered_map<slp::Vertex, slp::Vertex> VerticesMap;

class SkipMappedVerticesAcceptor: public inspector::TaskAcceptor {
  public:
  SkipMappedVerticesAcceptor(const VerticesMap const* vertices_map)
      : vertices_map_(vertices_map) {}

    //overriding accept
    bool accept(const inspector::InspectorTask& task) {
      return inspector::TaskAcceptor::accept(task)
        && vertices_map_->find(task.vertex) == vertices_map_->end();//do not accept if mapped vertex already
    }

    ::std::unique_ptr<inspector::TaskAcceptor> clone() const {
      return ::std::unique_ptr<inspector::TaskAcceptor>(new TaskAcceptor(*this));
    }
  private:
    const VerticesMap const* vertices_map_;
};


}//namespace mapper

//! Maps #root and its descendants using #f.
/*
 * It inspects vertices and if a vertex is already a key in #vertices_map then it is skipped
 * along with its descendants, otherwise after its descendants processed  the mapping {@code vertex => f(vertex)}
 * is added to #vertices_map.
 * @tparam VertexMapping function or functor mapping a vertex to a vertex
 */
template<typename VertexMapping>
void map_vertices(const slp::Vertex& root, std::unordered_map<slp::Vertex, slp::Vertex>* vertices_map, VertexMapping f) {
  if (vertices_map->find(root) != vertices_map->end())
    return;//root is already mapped

  slp::PostorderInspector inspector(root,
      std::unique_ptr<mapper::SkipMappedVerticesAcceptor>(new mapper::SkipMappedVerticesAcceptor(vertices_map)));
  while (!inspector.stopped()) {
    auto new_entry = std::make_pair(inspector.vertex(), f(inspector.vertex()));
    vertices_map->insert(new_entry);
    inspector.next();
  }
}


#endif /* CRAG_FREEGROUP_SLP_MAPPER_H_ */
