// Copyright (C) 2005 Alexander Ushakov

#ifndef CRAG_LINKED_BRAID_STRUCTURE_H
#define CRAG_LINKED_BRAID_STRUCTURE_H

#include <list>
#include <map>
#include <ostream>
#include <set>
#include <tuple>
#include <vector>

#include "Word.h"

struct BraidNode {
public:
  BraidNode(
      int64_t num = 0, bool tp = true, BraidNode* l = nullptr, BraidNode* a = nullptr, BraidNode* r = nullptr,
      BraidNode* bl = nullptr, BraidNode* b = nullptr, BraidNode* br = nullptr)
      : left(l)
      , ahead(a)
      , right(r)
      , back_left(bl)
      , back(b)
      , back_right(br)
      , type(tp)
      , link(nullptr)
      , weight(0)
      , theNumber(num) {}

public:
  BraidNode* left;
  BraidNode* ahead;
  BraidNode* right;

  BraidNode* back_left;
  BraidNode* back;
  BraidNode* back_right;

  // shows whether a crossing positive or negative
  bool type;

  // auxiliary member, used in LinkedBraidStructure copy constructor and in
  // translateIntoWord
  mutable BraidNode* link;

  // auxiliary member, used in remove Handle functions
  mutable int64_t weight;

  // a unique number associated with the node in a LBS
  int64_t theNumber;
};

std::ostream& operator<<(std::ostream& os, const BraidNode& bn);

struct LinkedBraidStructureTransform {
  enum TRANSFORM {
    ERASED,
    ADDED,
    CHANGE_TYPE,
  };

  LinkedBraidStructureTransform(
      int64_t n, int64_t p, TRANSFORM tr, bool t = false, int64_t l = 0, int64_t a = 0, int64_t r = 0, int64_t bl = 0,
      int64_t b = 0, int64_t br = 0)
      : theTransform(tr)
      , left(l)
      , ahead(a)
      , right(r)
      , back_left(bl)
      , back(b)
      , back_right(br)
      , type(t)
      , theNumber(n)
      , thePosition(p) {}

  TRANSFORM theTransform;

  int64_t left;
  int64_t ahead;
  int64_t right;

  int64_t back_left;
  int64_t back;
  int64_t back_right;

  bool type;

  int64_t theNumber;
  int64_t thePosition;
};

class LinkedBraidStructure {
public:
  LinkedBraidStructure(size_t N);
  LinkedBraidStructure(size_t N, const Word& w);

  LinkedBraidStructure(const LinkedBraidStructure& LBS);
  LinkedBraidStructure& operator=(const LinkedBraidStructure& LBS);

public:
  // function check if the braid-word represented by *this is shortlex smaller
  // than the one given by lbs, usually the computation of the words is not
  // required, so the function is often fast
  bool operator<(const LinkedBraidStructure& lbs) const;

  size_t size() const {
    return the_nodes_.size();
  }

  void clear();

  LinkedBraidStructureTransform push_back(int g);
  LinkedBraidStructureTransform push_front(int g);

  void removeLeftHandles(std::list<LinkedBraidStructureTransform>* result = nullptr);
  void removeRightHandles(std::list<LinkedBraidStructureTransform>* result = nullptr);

  Word translateIntoWord() const;

  void undo(const std::list<LinkedBraidStructureTransform>& lbst_seq);
  void undo(const LinkedBraidStructureTransform& lbst);

private:
  using NODE = std::tuple<int64_t, int64_t, BraidNode*>;

  struct NodeCompare {
    using int64_triple = std::tuple<int64_t, int64_t, int64_t>;

    bool operator()(const NODE& lhs, const NODE& rhs) const {
      return std::less<int64_triple>()(toInt64Triple_(lhs), toInt64Triple_(rhs));
    }

    int64_triple toInt64Triple_(const NODE& n) const {
      return std::make_tuple(std::get<0>(n), std::get<1>(n), std::get<2>(n)->theNumber);
    }
  };

  using node_set = std::set<NODE, NodeCompare>;

  LinkedBraidStructureTransform make_EraseTransform(BraidNode* bn, int64_t pos) const;
  LinkedBraidStructureTransform make_AddTransform(BraidNode* bn, int64_t pos) const;
  LinkedBraidStructureTransform make_ChangeType(BraidNode* bn, int64_t pos) const;

  int64_t checkIfStartsLeftHandle(int64_t pos, BraidNode* bn);
  int64_t checkIfStartsRightHandle(int64_t pos, BraidNode* bn);

  void removeLeftHandle(NODE node, node_set& to_check, std::list<LinkedBraidStructureTransform>* lst);
  void removeRightHandle(NODE node, node_set& to_check, std::list<LinkedBraidStructureTransform>* lst);

  LinkedBraidStructureTransform removeNode(BraidNode* bn, int64_t pos);

  BraidNode* insertBackRight(BraidNode* bn, int64_t pos, bool type);
  BraidNode* insertBackLeft(BraidNode* bn, int64_t pos, bool type);
  BraidNode* insert(const LinkedBraidStructureTransform& lbst);

  void processTree(int al, BraidNode* node, Word& result) const;

  void clearLinks() const;

private:
  //! Number of generators!!!
  size_t the_index_;

  std::vector<BraidNode*> front_nodes_;
  std::vector<BraidNode*> back_nodes_;

  std::map<int64_t, BraidNode> the_nodes_;

  int64_t max_node_number_;
};

//! Compares two braids using LinkedBraidStructure.
bool areEqualBraids(size_t n, const Word& lhs, const Word& rhs);

bool isTrivialBraid(size_t n, const Word& w);

#endif // CRAG_LINKED_BRAID_STRUCTURE_H
