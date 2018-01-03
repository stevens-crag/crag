/**
 * \file slp_common_prefix.h
 * \brief Algorithm to find the longest common prefix of two slps
 */

#pragma once
#ifndef CRAG_FREEGROUP_SLP_COMMON_PREFIX_H_
#define CRAG_FREEGROUP_SLP_COMMON_PREFIX_H_

#include <cassert>
#include <vector>
#include <gmpxx.h>

#include "slp_vertex.h"
#include "slp_inspector.h"
#include "slp_pattern_matching.h"

typedef mpz_class LongInteger;

namespace crag {
namespace slp {

namespace inspector {

class SegmentTracker {
  private:
    struct FollowState {
      Vertex vertex;
      LongInteger left_siblings_length;
    };

    const Vertex root_;

    std::vector<FollowState> state_history_;
    FollowState current_state_;

  public:
    SegmentTracker(const Vertex& root)
      : root_(root)
      , current_state_({root, 0})
    {}

    const Vertex& vertex() const {
      return current_state_.vertex;
    }

    const LongInteger& vertex_left_siblings() const {
      return current_state_.left_siblings_length;
    }

    void synchronize(const LongInteger& segment_begin, const LongInteger& segment_length) {
      assert(segment_begin + segment_length <= root_.length());

      LongInteger segment_end = segment_begin + segment_length;

      while (current_state_.vertex && (current_state_.left_siblings_length > segment_begin || current_state_.left_siblings_length + current_state_.vertex.length() < segment_end)) {
        assert(!state_history_.empty());

        current_state_ = std::move(state_history_.back());
        state_history_.pop_back();
      }

      while (current_state_.vertex) {//this condition should never be false
        LongInteger split_position = current_state_.vertex.split_point() + current_state_.left_siblings_length;
        if (split_position < segment_begin) {
          state_history_.push_back(std::move(current_state_));
          current_state_ = {state_history_.back().vertex.right_child(), state_history_.back().left_siblings_length + state_history_.back().vertex.split_point()};
        } else if (split_position > segment_end) {
          state_history_.push_back(std::move(current_state_));
          current_state_ = {state_history_.back().vertex.left_child(), state_history_.back().left_siblings_length};
        } else {
          return;
        }
      }
    }
};

template <typename AcceptFunctor = std::nullptr_t>
class LongestPrefixInspectorPath : public InspectorPath<AcceptFunctor> {
  public:
    LongestPrefixInspectorPath() {}
    LongestPrefixInspectorPath(Vertex other, const MatchingTable& matching_table)
      : InspectorPath<AcceptFunctor>(nullptr)
      , matching_table_(matching_table)
      , other_root_(std::move(other))
      , other_slp_inspector_(other_root_)
    { }

    using InspectorPath<AcceptFunctor>::schedule;
    using InspectorPath<AcceptFunctor>::is_task_accepted;
    using InspectorPath<AcceptFunctor>::pop_scheduled;
    using InspectorPath<AcceptFunctor>::clear_scheduled;

    InspectorTask process(const InspectorTask& current_task) {
      if (current_task.command == InspectorTask::Command::VISIT) {
        InspectorTask next_task = pop_scheduled();
        if (next_task && next_task.command == InspectorTask::Command::VISIT) {
          syncronize(next_task);
          if (!is_prefix(next_task)) {
            clear_scheduled();
            return InspectorTask();
          }
        }
        return next_task;
      }

      syncronize(current_task);
      if (matching_table_.is_calculated(current_task.vertex, other_slp_inspector_.vertex())) {
        if (is_prefix(current_task)) {
          return InspectorTask::for_current(current_task, InspectorTask::Command::VISIT);
        }
      }

      if (current_task.command == InspectorTask::Command::GO_LEFT) {
        schedule(InspectorTask::for_current(current_task, InspectorTask::Command::GO_RIGHT));
        return InspectorTask::for_left_child(current_task, InspectorTask::Command::GO_LEFT);
      } else {
        InspectorTask next_task = InspectorTask::for_right_child(current_task, InspectorTask::Command::GO_LEFT);

        if (!is_task_accepted(next_task)) {
          next_task = InspectorTask::for_current(current_task, InspectorTask::Command::VISIT);
          if (next_task && (next_task.command == InspectorTask::Command::VISIT)) {
            syncronize(next_task);
            if (!is_prefix(next_task)) {
              clear_scheduled();
              return InspectorTask();
            }
          }
          return next_task;
        }

        schedule(InspectorTask::for_current(current_task, InspectorTask::Command::VISIT));
        return next_task;
      }
    }

    InspectorTask initial_task(const Vertex& root) {
      return InspectorTask(Vertex(root), InspectorTask::Command::GO_LEFT, LongInteger());

//      while (current_task) {
//        syncronize(current_task);
//        if (matching_table_.is_calculated(current_task.vertex, other_slp_inspector_.vertex())) {
//          if (is_prefix(current_task)) {
//            return current_task;
//          }
//          current_task = InspectorTask::for_left_child(current_task, InspectorTask::Command::VISIT);
//        }
//      }
//      other_slp_inspector_.synchronize(current_task.left_siblings_length, current_task.vertex.length());
//
//      while (!matching_table_.is_calculated(current_task.vertex, other_slp_inspector_.vertex())) {
//        schedule(InspectorTask::for_right_child(current_task, InspectorTask::Command::GO_LEFT));
//        current_task = InspectorTask::for_left_child(current_task, InspectorTask::Command::VISIT);
//        syncronize(current_task);
//      }
//
//      if (!is_prefix(current_task)) {
//        current_task = InspectorTask();
//        clear_scheduled();
//      }
//
//      return current_task;
    }

    bool is_prefix(const InspectorTask& current_task) {
      auto matches = matching_table_.matches(current_task.vertex, other_vertex());
      matches.shift_right(other_vertex_left_siblings());
      return matches.contains(current_task.left_siblings_length);
    }

    const Vertex& other_vertex() const {
      return other_slp_inspector_.vertex();
    }

    const LongInteger& other_vertex_left_siblings() const {
      return other_slp_inspector_.vertex_left_siblings();
    }

  private:
    void syncronize(const InspectorTask& current_task) {
      other_slp_inspector_.synchronize(current_task.left_siblings_length, current_task.vertex.length());
    }

    MatchingTable matching_table_;
    Vertex other_root_;
    SegmentTracker other_slp_inspector_;
};
} //namespace inspector

LongInteger longest_common_prefix(const Vertex& first, const Vertex& second, MatchingTable* matching_table);

inline LongInteger longest_common_prefix(const Vertex& first, const Vertex& second) {
  MatchingTable temp;
  return longest_common_prefix(first, second, &temp);
}

} //namespace slp
} //namespace crag

#endif /* SLP_COMMON_PREFIX_H_ */
