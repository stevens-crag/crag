/*
 * slp_inspector.h
 *
 *  Created on: Feb 11, 2013
 *      Author: dpantele
 */

#pragma once
#ifndef CRAG_FREEGROUP_SLP_INSPECTOR_H_
#define CRAG_FREEGROUP_SLP_INSPECTOR_H_

#include "slp_vertex.h"
#include <vector>

namespace crag {
namespace slp {

namespace inspector{

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
      : vertex(Vertex::Null)
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
  return (first.vertex == Vertex::Null && second.vertex == Vertex::Null) ||
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

    InspectorTask initial_task(const Vertex& root) {
      return InspectorTask(root, InspectorTask::Command::VISIT, LongInteger(0));
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
    InspectorTask initial_task(const Vertex& root) {
      return InspectorTask(root, InspectorTask::Command::GO_LEFT, 0);
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
        return InspectorTask::DO_NOTHING;
      }
    }

    InspectorTask initial_task(const Vertex& root) {
      return InspectorTask(root, InspectorTask::Command::GO_LEFT, 0);
    }
};

//We can't use virtual functions from constructors, and we have to call next() in constructor.
//So we require to pass some functor which accepts tasks while we construct inspector
class TaskAcceptor {
  public:
    virtual ~TaskAcceptor() {}
    virtual bool accept(const InspectorTask& task) {
      return task != InspectorTask::DO_NOTHING && task.vertex != Vertex::Null;
    }

    TaskAcceptor* clone() const {
      return new TaskAcceptor(*this);
    }
};


}//namespace inspector

template <typename OrderPolicy>
class Inspector : public OrderPolicy {
  public:
    Inspector()
      : current_task_(InspectorTask::DO_NOTHING)
      , task_checker_(nullptr)
    { }

    virtual ~Inspector() {
      delete task_checker_;
    }

    Inspector(const Inspector& other)
      : current_task_(other.current_task_)
      , scheduled_tasks_(other.scheduled_tasks_)
      , task_checker_(other.task_checker_ ? other.task_checker_->clone() : nullptr)
    { }

    Inspector(const Vertex& root, inspector::TaskAcceptor* task_acceptor)
      : current_task_(this->initial_task(root))
      , task_checker_(task_acceptor)
    {
      if (!task_checker_->accept(current_task_)) {
        current_task_ = InspectorTask::DO_NOTHING;
      }

      if (current_task_.command != InspectorTask::Command::VISIT) {
        next();
      }
    }

    Inspector(const Vertex& root)
      : Inspector(root, new inspector::TaskAcceptor())
    { }

    Inspector& operator++() {
      return next();
    }

    Inspector operator++(int) {
      Inspector copy(*this);
      ++*this;
      return copy;
    }

    bool operator==(const Inspector& other) const {
      return stopped() ? other.stopped() : !other.stopped() &&
              current_task_.left_siblings_length == other.current_task_.left_siblings_length &&
              vertex() == other.vertex();
    }

    bool operator!=(const Inspector& other) const {
      return !(*this==other);
    }

    virtual Inspector& next() {
      if (stopped()) {
        return *this;
      }

      do {
        process_current_task();
      } while(!stopped() && current_task_.command != InspectorTask::Command::VISIT);

      return *this;
    }

    const Vertex& operator*() const {
      return vertex();
    }

    const Vertex& vertex() const {
      return current_task_.vertex;
    }

    const LongInteger& vertex_left_siblings_length() const {
      return current_task_.left_siblings_length;
    }

    bool stopped() const {
      return !current_task_;
    }

  protected:
    typedef inspector::InspectorTask InspectorTask;
    using OrderPolicy::initial_task;
    using OrderPolicy::process;

    InspectorTask current_task_;
    inspector::TaskAcceptor* task_checker_;

    ::std::vector<InspectorTask> scheduled_tasks_;

    void fallback_to_scheduled() {
      if (scheduled_tasks_.empty()) {
        current_task_ = InspectorTask::DO_NOTHING;
        return;
      }

      current_task_ = std::move(scheduled_tasks_.back());
      scheduled_tasks_.pop_back();
    }

    void process_current_task() {
      current_task_ = std::move(process(current_task_));

      if (!task_checker_->accept(current_task_)) {
        fallback_to_scheduled();
      }
    }



    void schedule(InspectorTask&& task) final {//marking this as final to safely use in constructor
      if (task_checker_->accept(task)) {
        scheduled_tasks_.push_back(std::move(task));
      }
    }

    void initialize_if_required() {
      if (scheduled_tasks_.size() == 1 && !current_task_) { //need initialization
        init();
      }
    }

    //We would do it in constructor, but is_task_valid must be virtual, and next() depends on it also
    virtual void init() {
      if (!task_checker_->accept(scheduled_tasks_.back())) { //fist task is invalid, so we do stop
        scheduled_tasks_.pop_back();
      } else {
        current_task_ = ::std::move(scheduled_tasks_.back());
        scheduled_tasks_.pop_back(); //at this point initialize_if_required() don't call init() anymore
        if (current_task_.command != InspectorTask::Command::VISIT) {
          next();
        }
      }
    }

};

typedef Inspector<inspector::Postorder> PostorderInspector;
typedef Inspector<inspector::Preorder> PreorderInspector;
typedef Inspector<inspector::Inorder> InorderInspector;

} //namespace slp
} //namespace crag

#endif /* CRAG_FREEGROUP_SLP_INSPECTOR_H_ */
