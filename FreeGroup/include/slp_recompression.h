/*
 * slp_recompression.h
 *
 *  Created on: May 16, 2013
 *      Author: dpantele
 */

#pragma once
#ifndef CRAG_FREEGROUP_SLP_RECOMPRESSION_H_
#define CRAG_FREEGROUP_SLP_RECOMPRESSION_H_

#include <memory>
#include <cassert>
#include <forward_list>
#include <map>
#include <type_traits>
#include <list>
#include <unordered_map>
#include <unordered_set>

#include <gmpxx.h>

#include "slp_vertex.h"
#include "slp_inspector.h"

namespace crag {
namespace slp {
namespace recompression {

typedef mpz_class LetterPower;
typedef int64_t TerminalId;

class Rule;

class RuleLetter {
  public:
    typedef std::list<RuleLetter>::iterator IterRef;

    RuleLetter(const RuleLetter& other)
      : terminal_id_(other.terminal_id_)
      , terminal_power_(other.terminal_power_)
      , nonterminal_rule_(other.nonterminal_rule_)
    { }

    bool is_valid() const {
      return terminal_id_ || nonterminal_rule_;
    }

    bool is_power_of(TerminalId terminal) const {
      return terminal_id_ == terminal;
    }

    bool is_power() const {
      return terminal_power_ > 1;
    }

    bool is_nonterminal() const {
      return nonterminal_rule_;
    }

    inline bool is_empty_nonterminal() const;

    TerminalId terminal_id() const {
      return terminal_id_;
    }

    const LetterPower& terminal_power() const {
      return terminal_power_;
    }

    Rule* nonterminal_rule() {
      return nonterminal_rule_;
    }

    const Rule* nonterminal_rule() const {
      return nonterminal_rule_;
    }

    inline TerminalId first_terminal_letter_id() const;

    inline TerminalId last_terminal_letter_id() const;

    static TerminalId fresh_terminal_id_;

    static TerminalId next_fresh_terminal() {
      ++fresh_terminal_id_;
      assert(fresh_terminal_id_ > 0);
      return fresh_terminal_id_;
    }

    static TerminalId last_terminal() {
      return fresh_terminal_id_;
    }

    explicit RuleLetter(Rule* vertex_rule)
      : terminal_id_(0)
      , terminal_power_()
      , nonterminal_rule_(vertex_rule)
    {
      assert(vertex_rule);
    }

    RuleLetter(TerminalId terminal, LetterPower power)
      : terminal_id_(terminal)
      , terminal_power_(std::move(power))
      , nonterminal_rule_(nullptr)
    {
      assert(terminal != 0);
    }

    void debug_print(std::ostream* os) const;

  private:
    friend class Rule;

    inline TerminalId* first_terminal_shared();
    inline TerminalId* last_terminal_shared();

  private:
    TerminalId terminal_id_;
    LetterPower terminal_power_;

    Rule* nonterminal_rule_;
};

struct LetterPosition;

//! Definition of rule.
/**
 * Rule is stored as a double-linked list with explicit garbage removal
 */
class Rule {
  public:
    typedef std::list<RuleLetter>::iterator iterator;
    typedef std::list<RuleLetter>::const_iterator const_iterator;

    Rule(std::initializer_list<RuleLetter> letters);

    Rule(const Rule& other) = delete;
    Rule(Rule&& other) = delete;
    //it is hard to move, since
    //we should change pointers in nonterminal_index_
//      : letters_(std::move(other.letters_))
//      , first_terminal_letter_(std::move(other.first_terminal_letter_))
//      , last_terminal_letter_(std::move(other.last_terminal_letter_))
//      , nonterminal_index_(std::move(other.nonterminal_index_))
//      , debug_id(other.debug_id)
//    { }

    bool empty() const {
      return letters_.empty();
    }

//    void clear();

    size_t size() const {
      return letters_.size();
    }

    iterator begin() {
      return letters_.begin();
    }

    const_iterator begin() const {
      return letters_.begin();
    }


    iterator end() {
      return letters_.end();
    }

    const_iterator end() const {
      return letters_.end();
    }

    iterator first_nonempty() {
      auto begin = letters_.begin();
      while (begin != letters_.end() && begin->is_empty_nonterminal()) {
        ++begin;
      }

      return begin;
    }

    iterator last_nonempty() {
      if (empty()) {
        return letters_.end();
      }
      auto last = --letters_.end();
      while (last != letters_.begin() && last->is_empty_nonterminal()) {
        --last;
      }

      if (last->is_empty_nonterminal()) {
        return letters_.end();
      }

      return last;
    }

    inline TerminalId first_terminal_id() const {
      return first_terminal_letter_ ? *first_terminal_letter_ : 0;
    }

    inline TerminalId last_terminal_id() const {
      return last_terminal_letter_ ? *last_terminal_letter_ : 0;
    }

    static void collect_garbage() {
      letters_garbage().clear();
    }

    RuleLetter::IterRef delete_letter(RuleLetter::IterRef letter) {
      letter->terminal_id_ = 0;
      letter->nonterminal_rule_ = nullptr;
      letters_garbage().splice(letters_garbage().end(), letters_, letter);
      return letter;
    }

    void pop_first();
    void pop_last();

    RuleLetter::IterRef remove_empty_letter(RuleLetter::IterRef position);
    void remove_if_empty();

    void compress_power(RuleLetter::IterRef position, TerminalId new_terminal);

    void compress_pair(
        RuleLetter::IterRef first,
        RuleLetter::IterRef second,
        TerminalId new_terminal
    );

    size_t debug_id;

    void debug_print(::std::ostream* os) const;

    void debug_print_exposed(::std::ostream* os) const;

  private:
    inline void register_inclusion(Rule* rule, RuleLetter::IterRef letter);

    void copy_first_terminal();
    void copy_last_terminal();

    void set_first_terminal(TerminalId*);
    void set_last_terminal(TerminalId*);

    void pop_right_from_letter(
        RuleLetter::IterRef letter_position,
        const RuleLetter& popped_letter);

    void pop_left_from_letter(
        RuleLetter::IterRef letter_position,
        const RuleLetter& popped_letter);

    static std::list<RuleLetter>& letters_garbage() {
      static std::list<RuleLetter> letters_garbage;
      return letters_garbage;
    }

  private:
    friend class RuleLetter;

    std::list<RuleLetter> letters_;
    TerminalId* first_terminal_letter_;
    TerminalId* last_terminal_letter_;
    std::vector<LetterPosition> nonterminal_index_;
};

//void PrintTo(const Rule& rule, ::std::ostream* os) {
//  rule.debug_print(os);
//}

struct LetterPosition {
    LetterPosition(Rule* rule, RuleLetter::IterRef letter)
      : rule_(rule)
      , letter_(letter)
    { }

    Rule* rule_;
    RuleLetter::IterRef letter_;
};

TerminalId RuleLetter::first_terminal_letter_id() const {
  assert(nonterminal_rule_ == nullptr || !nonterminal_rule_->empty());
  return nonterminal_rule_ ?
      (nonterminal_rule_->first_terminal_letter_ ?
          *(nonterminal_rule_->first_terminal_letter_) :
          0):
      terminal_id_;
}

TerminalId RuleLetter::last_terminal_letter_id() const {
  return nonterminal_rule_ ?
      (nonterminal_rule_->last_terminal_letter_ ?
          *(nonterminal_rule_->last_terminal_letter_) :
          0):
      terminal_id_;
}

TerminalId* RuleLetter::first_terminal_shared() {
  return nonterminal_rule_ ?
      nonterminal_rule_->first_terminal_letter_ :
      new TerminalId(terminal_id_);
}

TerminalId* RuleLetter::last_terminal_shared() {
  return nonterminal_rule_ ?
      nonterminal_rule_->last_terminal_letter_ :
      new TerminalId(terminal_id_);
}

inline bool RuleLetter::is_empty_nonterminal() const {
  return nonterminal_rule_ && nonterminal_rule_->empty();
}

inline void Rule::register_inclusion(Rule* rule, RuleLetter::IterRef letter) {
  nonterminal_index_.push_back(LetterPosition(rule, letter));
}

struct JezRules {
    JezRules(const Vertex& slp);

    void remove_crossing_blocks();
    std::vector<LetterPosition> list_blocks();

    void compress_blocks(const std::vector<LetterPosition>& blocks);

    std::list<Rule> rules_;
    std::unordered_map<Vertex, Rule*> vertex_rules_;
    std::map<TerminalId, Vertex> terminal_vertices_;
    std::unordered_map<Vertex, TerminalId> vertex_terminals_;

    void empty_cleanup();

    void debug_print(std::ostream* os) const;

  private:
    RuleLetter get_letter(const Vertex& vertex) {
      if (vertex.height() > 1) {
        assert(vertex_rules_.count(vertex));
        return RuleLetter(vertex_rules_[vertex]);
      }

      auto rules_terminal = vertex_terminals_.insert(std::make_pair(vertex, 0));

      if (rules_terminal.second) {
        (rules_terminal.first->second) = RuleLetter::next_fresh_terminal();
        terminal_vertices_[RuleLetter::last_terminal()] = vertex;
      }

      return RuleLetter(rules_terminal.first->second, 1);
    }

    Vertex get_exponent(Vertex vertex, LetterPower power);
};

class OneStepPairs {
  public:
    OneStepPairs(JezRules* rules);

    std::tuple<std::unordered_set<TerminalId>, std::unordered_set<TerminalId>>
    greedy_pairs() const;

    TerminalId compress_pair(
        TerminalId first,
        TerminalId second,
        const std::vector<LetterPosition>& occurencies
    );

    void remove_crossing(
        const std::unordered_set<TerminalId> & left_letters,
        const std::unordered_set<TerminalId> & right_letters
    );

    void compress_pairs_from_letter_lists(
        const std::unordered_set<TerminalId> & left_letters,
        const std::unordered_set<TerminalId> & right_letters
    );


  private:
    struct LetterLeftRight {
      TerminalId id_;

      LetterLeftRight(TerminalId id)
        : id_(id)
      { }

      struct RightListEntry {
          RightListEntry(TerminalId right_id)
            : id_(right_id)
          { }

          TerminalId id_;
          std::vector<LetterPosition> occurencies;
      };

      std::list<RightListEntry> right_letters_;
      std::list<TerminalId> left_letters_;
    };

    std::list<LetterLeftRight> pairs_;
    JezRules* rules_;
};

Vertex normal_form(Vertex root);

namespace mad_sorts {
inline unsigned char reverse_bits_in_byte(unsigned char byte) {
  return ((byte * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
}

template <size_t N>
struct bits_reverse {
    template <typename T>
    inline static T reverse_N_bits(T value) {
      auto res = ((static_cast<T>(reverse_bits_in_byte(value))) << (8 * (N-1))) | (bits_reverse<N-1>::reverse_N_bits(value >> 8));
//      std::cout << static_cast<uint64_t>(value) << std::endl;
//      std::cout << static_cast<uint64_t>(reverse_bits_in_byte(value)) << std::endl;
//      std::cout << res << std::endl << std::endl;
      return res;
    }
};

template <>
struct bits_reverse<1> {
    template <typename T>
    inline static T reverse_N_bits(T value) {
//      std::cout << static_cast<uint64_t>(value) << std::endl;
//      std::cout << static_cast<uint64_t>(reverse_bits_in_byte(value)) << std::endl << std::endl;
      return reverse_bits_in_byte(value);
    }
};

template <typename T>
inline T reverse_bits(T value) {
  return bits_reverse<sizeof(T)>::reverse_N_bits(value);
}

#define SIZ(x) ((x)->_mp_size)
#define PTR(x) ((x)->_mp_d)
#define ABS(x) ((x) >= 0 ? (x) : -(x))
inline bool reverse_bit_mpz_less(const LetterPower& first, const LetterPower& second) {
  mp_size_t  usize, vsize, msize, dsize, asize;
  mpz_srcptr u = first.get_mpz_t();
  mpz_srcptr v = second.get_mpz_t();
  mp_srcptr  up, vp;
  int        cmp;

  usize = SIZ(u);
  vsize = SIZ(v);
  dsize = usize - vsize;
  msize = dsize < 0 ? usize : vsize;

  assert(msize >= 0);

  up = PTR(u);
  vp = PTR(v);

  mp_size_t  gmp_i;
  mp_limb_t  gmp_x, gmp_y;

  gmp_i = 0;
  while (gmp_i < msize) {
    gmp_x = reverse_bits((up)[gmp_i]);
    gmp_y = reverse_bits((vp)[gmp_i]);
    if (gmp_x != gmp_y)
    {
      return gmp_x < gmp_y;
    }
    ++gmp_i;
  }

  return dsize < 0;
}

#undef SIZ
#undef PTR
#undef ABS
} //namespace mad_sorts

}
}
}


#endif /* CRAG_FREEGROUP_SLP_RECOMPRESSION_H_ */
