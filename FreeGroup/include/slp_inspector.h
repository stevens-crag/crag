/**
 * \file slp_inspector.h
 * \brief Definition of the utilities for tree traversal
 *
 * Definition of several classes to perform SLP vertex traversal in abstract way. See \ref slp_description for examples.
 */

#pragma once
#ifndef CRAG_FREEGROUP_SLP_INSPECTOR_H_
#define CRAG_FREEGROUP_SLP_INSPECTOR_H_

#include <vector>
#include <functional>
#include <type_traits>
#include "slp_vertex.h"
#include "boost/pool/pool_alloc.hpp"

namespace crag {
namespace slp {

namespace inspector{

//! Internal structure used in path policy
struct InspectorTask {
    Vertex vertex;

    enum class Command {
      GO_LEFT,
      GO_RIGHT,
      VISIT,
    };
    Command command;

    LongInteger left_siblings_length;

    const static InspectorTask DO_NOTHING;

    InspectorTask()
      : vertex()
      , command(Command::VISIT)
      , left_siblings_length()
    { }

    InspectorTask(Vertex vertex, Command command, LongInteger left_siblings_length)
      : vertex(std::move(vertex))
      , command(command)
      , left_siblings_length(std::move(left_siblings_length))
    { }

    static InspectorTask for_current(const InspectorTask& current, Command command) {
      return InspectorTask(current.vertex, command, current.left_siblings_length);
    }

    static InspectorTask for_left_child(const InspectorTask& current, Command command) {
      return InspectorTask(current.vertex.left_child(), command, current.left_siblings_length);
    }

    static InspectorTask for_right_child(const InspectorTask& current, Command command) {
      return InspectorTask(current.vertex.right_child(), command, current.left_siblings_length + current.vertex.left_child().length());
    }

    explicit operator bool() const {
      return static_cast<bool>(vertex);
    }
};

inline bool operator==(const InspectorTask& first, const InspectorTask& second) {
  return (!first.vertex && !second.vertex) ||
      (first.vertex == second.vertex && first.command == second.command && first.left_siblings_length == second.left_siblings_length);
}

inline bool operator!=(const InspectorTask& first, const InspectorTask& second) {
  return !(first == second);
}

template <typename Function>
struct AcceptPolicy {
  public:
    Function task_accept_functor_;
  protected:
    AcceptPolicy() {}
    AcceptPolicy(Function&& functor_accepts) : task_accept_functor_(::std::move(functor_accepts)) {}

    bool is_task_accepted(const InspectorTask& task) {
      return task && task.vertex && task_accept_functor_(task);
    }
};

template <>
struct AcceptPolicy<std::nullptr_t> {
  protected:
    AcceptPolicy() {}
    AcceptPolicy(std::nullptr_t&&) {}
    bool is_task_accepted(const InspectorTask& task) {
      return task && task.vertex;
    }
};

//! Core of inspector customizations
/**
 * This class is the base class for path strategy (see pattern strategy, or policy).
 * It defines the next task based on the current one. All other routine (calling
 * process()) is implemented in the #crag::slp::Inspector. If you
 * want to redefine this, see crag:slp::inspector::Preorder, crag:slp::inspector::Postorder,
 * crag:slp::inspector::Inorder, crag::slp::PatternMatchesGenerator, crag::slp::TerminalVertexPath
 * for examples.
 */
template <typename AcceptFunctor>
class InspectorPath : public AcceptPolicy<AcceptFunctor> {
  //redefine these methods
  public:
    InspectorTask process(const InspectorTask& current_task);
    InspectorTask initial_task(const Vertex& root) const;

  //here is just common part
  public:
    InspectorPath()
      : AcceptPolicy<AcceptFunctor>()
    { }

    InspectorTask pop_scheduled() {
      InspectorTask result;
      do {
        if (scheduled_tasks_.empty()) {
          return InspectorTask();
        }
        result = std::move(scheduled_tasks_.back());
        scheduled_tasks_.pop_back();
      } while(!is_task_accepted(result));
      return result;
    }

    void clear_scheduled() {
      scheduled_tasks_.clear();
    }

    using AcceptPolicy<AcceptFunctor>::is_task_accepted;
  protected:
    InspectorPath(AcceptFunctor&& accept_functor)
      : AcceptPolicy<AcceptFunctor>(std::move(accept_functor))
    { }

    const InspectorTask& schedule(InspectorTask&& task) {
      scheduled_tasks_.push_back(std::move(task));
      return scheduled_tasks_.back();
    }

    const InspectorTask& last_scheduled() const {
      return scheduled_tasks_.back();
    }

  private:
//    std::vector<InspectorTask> scheduled_tasks_;
    std::vector<InspectorTask, boost::pool_allocator<InspectorTask, boost::default_user_allocator_malloc_free, boost::details::pool::null_mutex, 128, 0>> scheduled_tasks_;
};


//! Define preorder tree traversal
template <typename AcceptFunctor>
class Preorder : public InspectorPath<AcceptFunctor> {
  public:
    Preorder() {}
    Preorder(AcceptFunctor&& accept_functor)
      : InspectorPath<AcceptFunctor>(std::move(accept_functor))
    { }
    using InspectorPath<AcceptFunctor>::schedule;

    InspectorTask process(const InspectorTask& current_task) {
      schedule(InspectorTask::for_right_child(current_task, InspectorTask::Command::VISIT));
      return InspectorTask::for_left_child(current_task, InspectorTask::Command::VISIT);
    }

    InspectorTask initial_task(Vertex&& root) const {
      return InspectorTask(std::move(root), InspectorTask::Command::VISIT, LongInteger());
    }


};

//! Define inorder tree traversal
template <typename AcceptFunctor>
class Inorder : public InspectorPath<AcceptFunctor> {
  public:
    Inorder() {}
    Inorder(AcceptFunctor&& accept_functor)
      : InspectorPath<AcceptFunctor>(std::move(accept_functor))
    { }
    using InspectorPath<AcceptFunctor>::schedule;

    InspectorTask process(const InspectorTask& current_task) {
      if (current_task.command == InspectorTask::Command::VISIT) {
        return InspectorTask::for_right_child(current_task, InspectorTask::Command::GO_LEFT);
      } else {
        schedule(InspectorTask::for_current(current_task, InspectorTask::Command::VISIT));
        return InspectorTask::for_left_child(current_task, InspectorTask::Command::GO_LEFT);
      }
    }
    InspectorTask initial_task(Vertex&& root) const {
      return InspectorTask(std::move(root), InspectorTask::Command::GO_LEFT, LongInteger());
    }
};

//! Define postorder tree traversal
template <typename AcceptFunctor>
class Postorder : public InspectorPath<AcceptFunctor> {
  public:
    Postorder() {}
    Postorder(AcceptFunctor&& accept_functor)
      : InspectorPath<AcceptFunctor>(std::move(accept_functor))
    { }
    using InspectorPath<AcceptFunctor>::schedule;

    InspectorTask process(const InspectorTask& current_task) {
      if (current_task.command == InspectorTask::Command::GO_LEFT) {
        schedule(InspectorTask::for_current(current_task, InspectorTask::Command::GO_RIGHT));
        return InspectorTask::for_left_child(current_task, InspectorTask::Command::GO_LEFT);
      } else if (current_task.command == InspectorTask::Command::GO_RIGHT) {
        schedule(InspectorTask::for_current(current_task, InspectorTask::Command::VISIT));
        return InspectorTask::for_right_child(current_task, InspectorTask::Command::GO_LEFT);
      } else {
        return InspectorTask();
      }
    }

    InspectorTask initial_task(Vertex&& root) const {
      return InspectorTask(::std::move(root), InspectorTask::Command::GO_LEFT, LongInteger());
    }
};

}//namespace inspector

//! Main class to walk through all vertices of SLP
/**
 * This class provides all methods to deal with SLP traversal
 *
 * @tparam PathPolicy The policy defining the function which basically provides the next vertex
 * @tparam AcceptFunctor Addition to PathPolicy which provides the possibility to consider some
 *                       vertices as Null vertices, so omitting them and all their children.
 */
template <template <typename> class PathPolicy, typename AcceptFunctor = std::nullptr_t>
class Inspector {
  public:
    //! Traverse the subtree of \a root using \a inspector_path to choose next vertex
    Inspector(Vertex root, PathPolicy<AcceptFunctor> inspector_path)
      : inspector_path_(::std::move(inspector_path))
      , current_task_(inspector_path_.initial_task(::std::move(root)))
    {
      init();
    }

    //! Traverse the subtree of \a root using \a accept_functor to consider some vertices as Null and default PathPolicy()
    Inspector(Vertex root, AcceptFunctor accept_functor)
      : inspector_path_(std::move(accept_functor))
      , current_task_(inspector_path_.initial_task(::std::move(root)))
    {
      init();
    }

    //! Traverse the subtree of \a root using default PathPolicy()
    Inspector(Vertex root)
      : inspector_path_()
      , current_task_(inspector_path_.initial_task(::std::move(root)))
    {
      init();
    }

    //! Empty inspector (stops immediately). Useful for 'end' iterators, and other 'Null' objects
    Inspector()
      : inspector_path_()
      , current_task_(InspectorTask())
    { }

    //! Check whether the position of this inspector is the same as the other one
    bool operator==(const Inspector& other) const {
      return stopped() ? other.stopped() : !other.stopped() &&
              current_task_.left_siblings_length == other.current_task_.left_siblings_length &&
              vertex() == other.vertex();
    }

    bool operator!=(const Inspector& other) const {
      return !(*this==other);
    }

    //! Moves to the next vertex.
    Inspector& next() {
      if (stopped()) {
        return *this;
      }

      do {
        process_current_task();
      } while(!stopped() && current_task_.command != InspectorTask::Command::VISIT);

      return *this;
    }

    //! Returns current vertex.
    const Vertex& vertex() const {
      return current_task_.vertex;
    }

    //! Returns the number of terminal vertices which the inspector, going from left to right and not skipping any vertices, would visit.
    const LongInteger& vertex_left_siblings_length() const {
      return current_task_.left_siblings_length;
    }

    //! Returns true if the inspection is ended.
    bool stopped() const {
      return !current_task_;
    }

    //! Alias for next()
    Inspector& operator++() {
      return next();
    }

    Inspector operator++(int) {
      Inspector copy(*this);
      ++*this;
      return copy;
    }

    //! Alias for vertex()
    const Vertex& operator*() const {
      return vertex();
    }

  public:
    PathPolicy<AcceptFunctor> inspector_path_;

  protected:
    typedef inspector::InspectorTask InspectorTask;

    void init() {
      if (!inspector_path_.is_task_accepted(current_task_)) {
        current_task_ = InspectorTask();
      }

      if (current_task_.command != InspectorTask::Command::VISIT) {
        next();
      }
    }

    InspectorTask current_task_;
    void process_current_task() {
      current_task_ = inspector_path_.process(current_task_);

      if (!inspector_path_.is_task_accepted(current_task_)) {
        current_task_ = inspector_path_.pop_scheduled();
      }
    }

};

//! Shortcut to inspect SLP in Post order (left child, right child, vertex)
typedef Inspector<inspector::Postorder> PostorderInspector;
//! Shortcut to inspect SLP in Pre order (vertex, left child, right child)
typedef Inspector<inspector::Preorder> PreorderInspector;
//! Shortcut to inspect SLP in In order (left child, vertex, right child)
typedef Inspector<inspector::Inorder> InorderInspector;

} //namespace slp
} //namespace crag

#endif /* CRAG_FREEGROUP_SLP_INSPECTOR_H_ */
