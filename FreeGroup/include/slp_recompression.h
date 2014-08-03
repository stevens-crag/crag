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
    bool is_nonterminal() const {
      return status_ & TERMINAL_NONTERMINAL_MASK;
    }

    //RuleLetter(const RuleLetter& other)
    //  : status_(
    //{
    //  if (other.is_nonterminal()) {
    //    nonterminal_rule_ = other.nonterminal_rule_;
    //    status_ = NON_TERMINAL;
    //  } else {
    //    new (&terminal_) Terminal(other.terminal_.id, LetterPower(other.terminal_.power));
    //    status_ = TERMINAL;
    //  }
    //  assert(is_valid());
    //}

    bool is_valid() const {
      return !(status_ & FLAG_INVALID);
    }

    bool is_power_of(TerminalId terminal) const {
      return !is_nonterminal() && terminal_.id == terminal;
    }

    bool is_power() const {
      return !is_nonterminal() && terminal_.power > 1;
    }

    inline bool is_empty_nonterminal() const;

    TerminalId terminal_id() const {
      return is_nonterminal() ? 0 : terminal_.id;
    }

    const LetterPower& terminal_power() const {
      return is_nonterminal() ? zero_power() : terminal_.power;
    }

    Rule* nonterminal_rule() {
      return is_nonterminal() ? nonterminal_rule_ : nullptr;
    }

    const Rule* nonterminal_rule() const {
      return is_nonterminal() ? nonterminal_rule_ : nullptr;
    }

    inline TerminalId first_terminal_letter_id() const;

    inline TerminalId last_terminal_letter_id() const;

    explicit RuleLetter(Rule* vertex_rule)
      : status_(NON_TERMINAL)
      , terminal_(0, 0)
      , nonterminal_rule_(vertex_rule)
    {
      assert(vertex_rule);
      assert(is_nonterminal());
      assert(is_valid());
    }

    RuleLetter(TerminalId terminal, LetterPower power)
      : status_(TERMINAL)
      , terminal_(terminal, std::move(power))
    {
      assert(!is_nonterminal());
      assert(status_ == 0);
      assert(terminal != 0);
      assert(is_valid());

      //new (&terminal_) Terminal();
    }

    ~RuleLetter() {
      //if (!is_nonterminal()) {
      //  terminal_.~Terminal();
      //}
    }

    void debug_print(std::ostream* os) const;

  private:
    friend class Rule;

    inline const RuleLetter* first_terminal_ptr() const;

    inline RuleLetter* first_terminal_ptr() {
      return const_cast<RuleLetter*>(static_cast<const RuleLetter*>(this)->first_terminal_ptr());
    }

    inline const RuleLetter* last_terminal_ptr() const;

    inline RuleLetter* last_terminal_ptr() {
      return const_cast<RuleLetter*>(static_cast<const RuleLetter*>(this)->last_terminal_ptr());
    }
  private:
    static const LetterPower& zero_power() {
      static LetterPower zero = 0;
      return zero;
    }

    void make_invalid() {
      status_ |= FLAG_INVALID;
      assert(!is_valid());
    }

    struct Terminal {
        TerminalId id;
        LetterPower power;

        Terminal(TerminalId id, LetterPower&& power)
          : id(id)
          , power(std::move(power))
        { }
    };

    unsigned char status_; //bitset, #0 - is_nonterminal, #1 - is_invalid

    static const unsigned char TERMINAL = 0u;
    static const unsigned char NON_TERMINAL = 1u;
    static const unsigned char TERMINAL_NONTERMINAL_MASK = 1u;
    static const unsigned char FLAG_INVALID = (1u << 1);

    Terminal terminal_;
    Rule* nonterminal_rule_ = nullptr;

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
    //it is hard to move, since
    //we should change pointers in nonterminal_index_
    Rule(Rule&& other) = delete;

    bool empty() const {
      return letters_.empty();
    }

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

    inline TerminalId first_terminal_id() const {
      if (!first_terminal_letter_) {
        return 0;
      }

      if (!first_terminal_letter_->is_valid()) {
        first_terminal_letter_ = const_cast<RuleLetter*>(letters_.front().first_terminal_ptr());
      }

      assert(first_terminal_letter_->is_valid());

      return first_terminal_letter_->terminal_id();
    }

    inline TerminalId last_terminal_id() const {
      if (!last_terminal_letter_) {
        return 0;
      }

      if (!last_terminal_letter_->is_valid()) {
        last_terminal_letter_ = const_cast<RuleLetter*>(letters_.back().last_terminal_ptr());
      }

      assert(last_terminal_letter_->is_valid());

      return last_terminal_letter_->terminal_id();
    }

    static void collect_garbage() {
      letters_garbage().clear();
    }

    iterator pop_first_from_letter(iterator letter_position);
    iterator pop_last_from_letter(iterator letter_position);

    Rule::const_iterator delete_letter(iterator letter) {
      letter->make_invalid();
      letters_garbage().splice(letters_garbage().end(), letters_, letter);
      return letter;
    }

    std::pair<iterator, iterator> remove_empty_letter(iterator position);

    Rule::iterator compress_power(iterator position, TerminalId new_terminal);

    Rule::iterator compress_pair(
        iterator first,
        iterator second,
        TerminalId new_terminal
    );

    size_t debug_id;

    void debug_print(::std::ostream* os) const;

    void debug_print_exposed(::std::ostream* os) const;

  private:
    inline void register_inclusion(Rule* rule, iterator letter);

    void set_first_terminal(TerminalId*);
    void set_last_terminal(TerminalId*);

    void insert_popped_letter_right(
        iterator letter_position,
        const RuleLetter& popped_letter);

    void insert_popped_letter_left(
        iterator letter_position,
        const RuleLetter& popped_letter);

    static std::list<RuleLetter>& letters_garbage() {
      static std::list<RuleLetter> letters_garbage;
      return letters_garbage;
    }

  private:
    friend class RuleLetter;

    std::list<RuleLetter> letters_;
    mutable RuleLetter* first_terminal_letter_;
    mutable RuleLetter* last_terminal_letter_;
    std::vector<LetterPosition> nonterminal_index_;
};

//void PrintTo(const Rule& rule, ::std::ostream* os) {
//  rule.debug_print(os);
//}

struct LetterPosition {
    LetterPosition(Rule* rule, Rule::iterator letter)
      : rule_(rule)
      , letter_(letter)
    { }

    Rule* rule_;
    Rule::iterator letter_;
};

TerminalId RuleLetter::first_terminal_letter_id() const {
  assert(!is_nonterminal() || !nonterminal_rule_->empty());

  if (!is_nonterminal()) {
    return terminal_.id;
  }

  if (nonterminal_rule_->first_terminal_letter_) {
    return first_terminal_ptr()->terminal_id();
  } else {
    return 0;
  }
}

TerminalId RuleLetter::last_terminal_letter_id() const {
  assert(!is_nonterminal() || !nonterminal_rule_->empty());

  if (!is_nonterminal()) {
    return terminal_.id;
  }

  if (nonterminal_rule_->last_terminal_letter_) {
    return last_terminal_ptr()->terminal_id();
  } else {
    return 0;
  }
}

const RuleLetter* RuleLetter::first_terminal_ptr() const {
  if (!is_nonterminal()) {
    return this;
  }

  if (!nonterminal_rule_->first_terminal_letter_) {
    return nonterminal_rule_->first_terminal_letter_;
  }

  if (!nonterminal_rule_->first_terminal_letter_->is_valid()) {
    nonterminal_rule_->first_terminal_letter_ = nonterminal_rule_->letters_.front().first_terminal_ptr();
  }

  assert(nonterminal_rule_->first_terminal_letter_->is_valid());

  return nonterminal_rule_->first_terminal_letter_;
}

const RuleLetter* RuleLetter::last_terminal_ptr() const {
  if (!is_nonterminal()) {
    return this;
  }

  if (!nonterminal_rule_->last_terminal_letter_) {
    return nonterminal_rule_->last_terminal_letter_;
  }

  if (!nonterminal_rule_->last_terminal_letter_->is_valid()) {
    nonterminal_rule_->last_terminal_letter_ = nonterminal_rule_->letters_.back().last_terminal_ptr();
  }

  assert(nonterminal_rule_->last_terminal_letter_->is_valid());

  return nonterminal_rule_->last_terminal_letter_;
}

inline bool RuleLetter::is_empty_nonterminal() const {
  return is_nonterminal() && nonterminal_rule_->empty();
}

inline void Rule::register_inclusion(Rule* rule, Rule::iterator letter) {
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

    TerminalId next_fresh_terminal() {
      ++fresh_terminal_id_;
      assert(fresh_terminal_id_ > 0);
      return fresh_terminal_id_;
    }

    TerminalId last_terminal() const {
      return fresh_terminal_id_;
    }

  private:
    TerminalId fresh_terminal_id_;

    RuleLetter get_letter(const Vertex& vertex) {
      if (vertex.height() > 1) {
        assert(vertex_rules_.count(vertex));
        return RuleLetter(vertex_rules_[vertex]);
      }

      auto rules_terminal = vertex_terminals_.insert(std::make_pair(vertex, 0));

      if (rules_terminal.second) {
        (rules_terminal.first->second) = next_fresh_terminal();
        terminal_vertices_[rules_terminal.first->second] = vertex;
      }

      return RuleLetter(rules_terminal.first->second, 1);
    }

    Vertex get_exponent(Vertex vertex, LetterPower power);
};

class OneStepPairs {
  public:
    OneStepPairs(JezRules* rules);

    std::tuple<std::vector<unsigned char>, std::vector<unsigned char>>
    greedy_pairs() const;

    TerminalId compress_pair(
        TerminalId first,
        TerminalId second,
        const std::vector<LetterPosition>& occurencies
    );

    void remove_crossing(
        const std::vector<unsigned char> & left_letters,
        const std::vector<unsigned char> & right_letters
    );

    void compress_pairs_from_letter_lists(
        const std::vector<unsigned char> & left_letters,
        const std::vector<unsigned char> & right_letters
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
      return res;
    }
};

template <>
struct bits_reverse<1> {
    template <typename T>
    inline static T reverse_N_bits(T value) {
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
  mp_size_t  usize, vsize, msize, dsize;
  mpz_srcptr u = first.get_mpz_t();
  mpz_srcptr v = second.get_mpz_t();
  mp_srcptr  up, vp;

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
