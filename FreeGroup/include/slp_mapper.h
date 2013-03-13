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
#include <functional>
#include "slp_vertex.h"

namespace crag {
namespace slp {

namespace mapper {

//! Functor that accepts only vertices that were not mapped already.
/**
 * @tparam MappedVertices collection of vertices such that we can test the membership
 *                        by checking whether obj.find(key) == obj.end(key) as for std::map
 */
template<typename MappedVertices>
class SkipMappedVerticesAcceptor {
  public:
    //default copy constructor/assignment
    SkipMappedVerticesAcceptor()
      : mapped_vertices_(nullptr)
    { }

    SkipMappedVerticesAcceptor(const MappedVertices& mapped_vertices)
      : mapped_vertices_(&mapped_vertices) {}

    bool operator()(const inspector::InspectorTask& task) {
      return mapped_vertices_->find(task.vertex) == mapped_vertices_->end();//do not accept if vertex is mapped already
    }
  private:
    const MappedVertices* mapped_vertices_;
};
}//namespace mapper


//! Maps #root and its descendants using #f.
/*
 * It inspects vertices and if a vertex is already a key in #vertices_map then it is skipped
 * along with its descendants, otherwise after its descendants processed  the mapping {@code vertex => f(vertex)}
 * is added to #vertices_map.
 * @tparam Func function or functor accepting signature ImageType (const slp::Vertex&, const std::unordered_map<slp::Vertex, ImageType>&), i.e. mapping slp::Vertex to ImageType
 * @tparam ImageType type of object which is an image of slp::Vertex under Func
 */
template<typename ImageType, typename Func>
void map_vertices(const slp::Vertex& root, std::unordered_map<slp::Vertex, ImageType>* p_images, Func f) {
  if (p_images->find(root) != p_images->end())
    return;//root is already mapped

  typedef typename mapper::SkipMappedVerticesAcceptor<std::unordered_map<slp::Vertex, ImageType>> Acceptor;
  Acceptor acceptor(*p_images);
  slp::Inspector<inspector::Postorder, Acceptor> inspector(root, acceptor);
  while (!inspector.stopped()) {
    auto img = f(inspector.vertex(), p_images);
    auto new_entry = std::make_pair(inspector.vertex(), img);
    p_images->insert(new_entry);
    inspector.next();
  }
}

}//namespace slp
}//namespace crag

#endif /* CRAG_FREEGROUP_SLP_MAPPER_H_ */
