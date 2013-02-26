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

    InspectorTask(const Vertex& vertex, Command&& command, const LongInteger& left_siblings_length)
      : vertex(vertex)
      , command(command)
      , left_siblings_length(left_siblings_length)
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


}//namespace inspector

template <typename OrderPolicy>
class Inspector : public OrderPolicy {
  public:
    Inspector()
      : current_task_(InspectorTask::DO_NOTHING)
    { }

    virtual ~Inspector() {}

    Inspector(const Vertex& root)
      : current_task_(this->initial_task(root)) {
      if (!is_task_valid(current_task_)) {
        current_task_ = InspectorTask::DO_NOTHING;
      }
      if (current_task_.command != InspectorTask::Command::VISIT) {
        next();
      }
    }

    Inspector& operator++() {
      return next();
    }

    Inspector operator++(int) const {
      return Inspector(*this).next();
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
      return current_task_ == InspectorTask::DO_NOTHING;
    }

  protected:
    typedef inspector::InspectorTask InspectorTask;
    InspectorTask current_task_;

    ::std::vector<InspectorTask> task_stack_;

    using OrderPolicy::initial_task;
    using OrderPolicy::process;

    void fallback_to_scheduled() {
      if (task_stack_.empty()) {
        current_task_ = InspectorTask::DO_NOTHING;
        return;
      }

      current_task_ = std::move(task_stack_.back());
      task_stack_.pop_back();
    }

    void process_current_task() {
      current_task_ = std::move(process(current_task_));

      if (!is_task_valid(current_task_)) {
        fallback_to_scheduled();
      }
    }

    void schedule(InspectorTask&& task) {
      if (is_task_valid(task)) {
        task_stack_.push_back(std::move(task));
      }
    }

    virtual bool is_task_valid(const InspectorTask& task) {
      return task != InspectorTask::DO_NOTHING && task.vertex != Vertex::Null;
    }
};

typedef Inspector<inspector::Postorder> PostorderInspector;
typedef Inspector<inspector::Preorder> PreorderInspector;
typedef Inspector<inspector::Inorder> InorderInspector;

} //namespace slp
} //namespace crag

#endif /* CRAG_FREEGROUP_SLP_INSPECTOR_H_ */
