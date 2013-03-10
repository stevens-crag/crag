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
      : vertex(::std::move(vertex))
      , command(::std::move(command))
      , left_siblings_length(::std::move(left_siblings_length))
    { }

    static InspectorTask for_current(const InspectorTask& current, Command&& command) {
      return InspectorTask(current.vertex, std::move(command), current.left_siblings_length);
    }

    static InspectorTask for_left_child(const InspectorTask& current, Command&& command) {
      return InspectorTask(current.vertex.left_child(), std::move(command), current.left_siblings_length);
    }

    static InspectorTask for_right_child(const InspectorTask& current, Command&& command) {
      return InspectorTask(current.vertex.right_child(), std::move(command), current.left_siblings_length + current.vertex.left_child().length());
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


class OrderPolicy {
  protected:
    virtual ~OrderPolicy() {}
    virtual void schedule(InspectorTask&& task) = 0;
    InspectorTask process(const InspectorTask& current_task);
    InspectorTask initial_task(const Vertex& root);
};

class Preorder : public OrderPolicy {
  protected:
    InspectorTask process(const InspectorTask& current_task) {
      schedule(InspectorTask::for_right_child(current_task, InspectorTask::Command::VISIT));
      return InspectorTask::for_left_child(current_task, InspectorTask::Command::VISIT);
    }

    InspectorTask initial_task(Vertex&& root) {
      return InspectorTask(::std::move(root), InspectorTask::Command::VISIT, LongInteger(0));
    }
};

class Inorder : public OrderPolicy {
  protected:
    InspectorTask process(const InspectorTask& current_task) {
      if (current_task.command == InspectorTask::Command::VISIT) {
        return InspectorTask::for_right_child(current_task, InspectorTask::Command::GO_LEFT);
      } else {
        schedule(InspectorTask::for_current(current_task, InspectorTask::Command::VISIT));
        return InspectorTask::for_left_child(current_task, InspectorTask::Command::GO_LEFT);
      }
    }
    InspectorTask initial_task(Vertex&& root) {
      return InspectorTask(::std::move(root), InspectorTask::Command::GO_LEFT, 0);
    }
};

class Postorder : public OrderPolicy {
  protected:
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

    InspectorTask initial_task(Vertex&& root) {
      return InspectorTask(::std::move(root), InspectorTask::Command::GO_LEFT, 0);
    }
};

template <typename Function>
struct AcceptPolicy {
  protected:

    Function functor_accepts;

    AcceptPolicy(Function&& functor_accepts) : functor_accepts(::std::move(functor_accepts)) {}

    bool is_task_accepted(const InspectorTask& task) {
      return task && task.vertex && functor_accepts(task);
    }
};

template <>
struct AcceptPolicy<std::nullptr_t> {
  protected:
    AcceptPolicy(std::nullptr_t&&) {}
    bool is_task_accepted(const InspectorTask& task) {
      return task && task.vertex;
    }
};

}//namespace inspector


template <typename OrderPolicy, typename AcceptFunctor = std::nullptr_t>
class Inspector : public OrderPolicy, public inspector::AcceptPolicy<AcceptFunctor> {
  public:

    virtual ~Inspector() {}

    Inspector(Vertex root, AcceptFunctor task_acceptor)
      : inspector::AcceptPolicy<AcceptFunctor>(::std::move(task_acceptor))
      , current_task_(this->initial_task(::std::move(root)))
    {
      if (!is_task_accepted(current_task_)) {
        current_task_ = InspectorTask();
      }

      if (current_task_.command != InspectorTask::Command::VISIT) {
        next();
      }
    }

    Inspector(Vertex root)
      : Inspector(::std::move(root), nullptr)
    { }

    Inspector()
      : inspector::AcceptPolicy<AcceptFunctor>(nullptr)
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
    virtual Inspector& next() {
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
  protected:
    typedef inspector::InspectorTask InspectorTask;
    using OrderPolicy::initial_task;
    using OrderPolicy::process;
    using inspector::AcceptPolicy<AcceptFunctor>::is_task_accepted;


    InspectorTask current_task_;

    ::std::vector<InspectorTask> scheduled_tasks_;


    void fallback_to_scheduled() {
      do {
        if (scheduled_tasks_.empty()) {
          current_task_ = InspectorTask();
          return;
        }
        current_task_ = std::move(scheduled_tasks_.back());
        scheduled_tasks_.pop_back();
      } while(!is_task_accepted(current_task_));
    }

    void process_current_task() {
      current_task_ = std::move(process(current_task_));

      if (!is_task_accepted(current_task_)) {
        fallback_to_scheduled();
      }
    }

    void schedule(InspectorTask&& task) final {//marking this as final to safely use in constructor
      scheduled_tasks_.push_back(std::move(task));
    }
};

typedef Inspector<inspector::Postorder> PostorderInspector;
typedef Inspector<inspector::Preorder> PreorderInspector;
typedef Inspector<inspector::Inorder> InorderInspector;

} //namespace slp
} //namespace crag

#endif /* CRAG_FREEGROUP_SLP_INSPECTOR_H_ */
