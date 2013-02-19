/*
 * slp_inspector.h
 *
 *  Created on: Feb 11, 2013
 *      Author: dpantele
 */

#ifndef CRAG_FREEGROUP_SLP_INSPECTOR_H_
#define CRAG_FREEGROUP_SLP_INSPECTOR_H_

#include "slp_vertex.h"

namespace crag {
namespace slp {
namespace inspector{

struct InspectorTask {
    Vertex::Ptr vertex;
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
      , command(Command::GO_LEFT)
      , left_siblings_length()
    { }

    InspectorTask(const Vertex::Ptr& vertex, Command&& command, const LongInteger& left_siblings_length)
      : vertex(vertex)
      , command(command)
      , left_siblings_length(left_siblings_length)
    { }

};

bool operator==(const InspectorTask& first, const InspectorTask& second) {
  return (!first.vertex && !second.vertex) ||
      (*first.vertex == *second.vertex && first.command == second.command && first.left_siblings_length == second.left_siblings_length);
}

bool operator!=(const InspectorTask& first, const InspectorTask& second) {
  return !(first == second);
}


const InspectorTask InspectorTask::DO_NOTHING = InspectorTask();

class OrderPolicy {
  protected:
    void schedule(const InspectorTask& task);
    InspectorTask process(const InspectorTask& current_task);
    InspectorTask initial_task(const Vertex::Ptr& root);
};

class Preorder : public OrderPolicy {
  protected:
    InspectorTask process(const InspectorTask& current_task) {
      schedule(InspectorTask(
        current_task.vertex->right_child(),
        InspectorTask::Command::VISIT,
        current_task.left_siblings_length + current_task.vertex->left_child()->length()
      ));

      return InspectorTask(current_task.vertex->left_child(), InspectorTask::Command::VISIT, current_task.left_siblings_length);
    }

    InspectorTask initial_task(const Vertex::Ptr& root) {
      return InspectorTask(root, InspectorTask::Command::VISIT, LongInteger(0));
    }
};

class Inorder : public OrderPolicy {
  protected:
    InspectorTask process(const InspectorTask& current_task) {
      if (current_task.command == InspectorTask::Command::VISIT) {
        return InspectorTask(
          current_task.vertex->right_child(),
          InspectorTask::Command::GO_LEFT,
          current_task.left_siblings_length + current_task.vertex->left_child()->length()
        );
      } else {
        schedule(InspectorTask(current_task.vertex, InspectorTask::Command::VISIT, current_task.left_siblings_length));
        return InspectorTask(current_task.vertex->left_child(), InspectorTask::Command::GO_LEFT, current_task.left_siblings_length);
      }
    }
    InspectorTask initial_task(const Vertex::Ptr& root) {
      return InspectorTask(root, InspectorTask::Command::GO_LEFT, 0);
    }
};

class Postorder : public OrderPolicy {
  protected:
    InspectorTask process(const InspectorTask& current_task) {
      if (current_task.command == InspectorTask::Command::GO_LEFT) {
        schedule(InspectorTask(current_task.vertex, InspectorTask::Command::GO_RIGHT, current_task.left_siblings_length));

        return InspectorTask(current_task.vertex->left_child(), InspectorTask::Command::GO_LEFT, current_task.left_siblings_length);
      } else if (current_task.command == InspectorTask::Command::GO_RIGHT) {
        schedule(InspectorTask(current_task.vertex, InspectorTask::Command::VISIT, current_task.left_siblings_length));

        return InspectorTask(
          current_task.vertex->right_child(),
          InspectorTask::Command::GO_LEFT,
          current_task.left_siblings_length
        );
      } else {
        return InspectorTask::DO_NOTHING;
      }
    }

    InspectorTask initial_task(const Vertex::Ptr& root) {
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

    Inspector(const Vertex::Ptr& root)
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

    Inspector operator++() const {
      return Inspector(*this).next();
    }

    bool operator==(const Inspector& other) {
      return stopped() ? other.stopped() :
              current_task_.left_siblings_length == other.current_task_.left_siblings_length &&
              (*vertex()) == (*other.vertex());
    }

    Inspector& next() {
      do {
        process_current_task();
      } while(current_task_.command != InspectorTask::Command::VISIT);
      return *this;
    }

    const Vertex::Ptr& operator*() const {
      return vertex();
    }

    const Vertex::Ptr& vertex() const {
      return current_task_.vertex;
    }

    bool stopped() const {
      return current_task_ == InspectorTask::DO_NOTHING;
    }

  protected:
    typedef inspector::InspectorTask InspectorTask;
    InspectorTask current_task_;

    std::vector<InspectorTask> task_stack_;

    using OrderPolicy::initial_task;
    using OrderPolicy::process;

    void process_current_task() {
      current_task_ = std::move(process(current_task_));

      if (!is_task_valid(current_task_)) {
        if (task_stack_.empty()) {
          current_task_ = InspectorTask::DO_NOTHING;
          return;
        }

        current_task_ = std::move(task_stack_.back());
        task_stack_.pop_back();
      }
    }

    void schedule(InspectorTask&& task) {
      if (is_task_valid(task)) {
        task_stack_.push_back(std::move(task));
      }
    }

    bool is_task_valid(const InspectorTask& task) {
      return task != InspectorTask::DO_NOTHING && task.vertex != Vertex::Null;
    }
};

typedef Inspector<inspector::Postorder> PostorderInspector;
typedef Inspector<inspector::Preorder> PreorderInspector;
typedef Inspector<inspector::Inorder> InorderInspector;

} //namespace slp
} //namespace crag

#endif /* CRAG_FREEGROUP_SLP_INSPECTOR_H_ */
