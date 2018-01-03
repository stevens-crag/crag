/*
 * slp_inspector.h
 *
 *  Created on: Feb 11, 2013
 *      Author: dpantele
 */

#pragma once
#ifndef CRAG_FREEGROUP_SLP_INSPECTOR_H_
#define CRAG_FREEGROUP_SLP_INSPECTOR_H_

#include <vector>
#include <functional>
#include <type_traits>
#include "slp_vertex.h"
#include "boost/pool/pool_alloc.hpp"

/**
 * Module defining inspector over SLP representation.
 *
 * An inspector follows given order policy, there are several predefined ones and the corresponding inspectors {@link PostorderInspector}, ...
 *
 * An inspector action (to visit a vertex, to go to the left from a vertex, etc) is described by InspectorTask.
 * When inspector is about to perform an action the corresponding InspectorTask object is passed to the method #accept
 * of the associated with this inspector TaskAcceptor, and if the call returns true the inspector performs the action, and skips it otherwise moving
 * to the next action defined by the order policy.
 * To change the default acceptor, extend class TaskAcceptor and override its method #accept, and then pass it to the corresponding inspector constructor.
 *
 * Example:
 * {@code
 *  PostorderInspector inspector(slp_vertex, ::std::unique_ptr<MyTaskAcceptor>(new MyTaskAcceptor(...));
 *  while (!inspector.stopped()) {
 *    //some code here using inspector.vertex()
 *    ...
 *    inspector.next();
 *  }
 * }
 */

namespace crag {
namespace slp {

namespace inspector{

//!
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


template <template <typename> class PathPolicy, typename AcceptFunctor = std::nullptr_t>
class Inspector {
  public:
    Inspector(Vertex root, PathPolicy<AcceptFunctor> inspector_path)
      : inspector_path_(::std::move(inspector_path))
      , current_task_(inspector_path_.initial_task(::std::move(root)))
    {
      init();
    }

    Inspector(Vertex root, AcceptFunctor accept_functor)
      : inspector_path_(std::move(accept_functor))
      , current_task_(inspector_path_.initial_task(::std::move(root)))
    {
      init();
    }

    Inspector(Vertex root)
      : inspector_path_()
      , current_task_(inspector_path_.initial_task(::std::move(root)))
    {
      init();
    }

    Inspector()
      : inspector_path_()
      , current_task_(InspectorTask())
    { }

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

    //! Returns the current vertex.
    const Vertex& vertex() const {
      return current_task_.vertex;
    }

    //! Returns the number of times an inspector going from left to right not skipping any vertices would visit a terminal vertex.
    const LongInteger& vertex_left_siblings_length() const {
      return current_task_.left_siblings_length;
    }

    //! Returns true if the inspection is ended.
    bool stopped() const {
      return !current_task_;
    }

    Inspector& operator++() {
      return next();
    }

    Inspector operator++(int) {
      Inspector copy(*this);
      ++*this;
      return copy;
    }

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

typedef Inspector<inspector::Postorder> PostorderInspector;
typedef Inspector<inspector::Preorder> PreorderInspector;
typedef Inspector<inspector::Inorder> InorderInspector;

} //namespace slp
} //namespace crag

#endif /* CRAG_FREEGROUP_SLP_INSPECTOR_H_ */
