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
  BraidNode(int num = 0, bool tp = true, BraidNode* l = NULL, BraidNode* a = NULL, BraidNode* r = NULL,
            BraidNode* bl = NULL, BraidNode* b = NULL, BraidNode* br = NULL)
      : left(l)
      , ahead(a)
      , right(r)
      , back_left(bl)
      , back(b)
      , back_right(br)
      , type(tp)
      , link(NULL)
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
  mutable int weight;

  // a unique number associated with the node in a LBS
  int theNumber;
};

std::ostream& operator<<(std::ostream& os, const BraidNode& bn);

struct LinkedBraidStructureTransform {
  enum TRANSFORM {
    ERASED,
    ADDED,
    CHANGE_TYPE,
  };

  LinkedBraidStructureTransform(int n, int p, TRANSFORM tr, bool t = false, int l = 0, int a = 0, int r = 0, int bl = 0,
                                int b = 0, int br = 0)
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

  int left;
  int ahead;
  int right;

  int back_left;
  int back;
  int back_right;

  bool type;
  int theNumber;
  int thePosition;
};

class LinkedBraidStructure {
public:
  LinkedBraidStructure(int N);
  LinkedBraidStructure(int N, const Word& w);
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

  void removeLeftHandles(std::list<LinkedBraidStructureTransform>* result = NULL);
  void removeRightHandles(std::list<LinkedBraidStructureTransform>* result = NULL);

  Word translateIntoWord() const;

  void undo(const std::list<LinkedBraidStructureTransform>& lbst_seq);
  void undo(const LinkedBraidStructureTransform& lbst);

private:
  using NODE = std::tuple<int, int, BraidNode*>;

  LinkedBraidStructureTransform make_EraseTransform(BraidNode* bn, int pos) const;
  LinkedBraidStructureTransform make_AddTransform(BraidNode* bn, int pos) const;
  LinkedBraidStructureTransform make_ChangeType(BraidNode* bn, int pos) const;

  int checkIfStartsLeftHandle(int pos, BraidNode* bn);
  int checkIfStartsRightHandle(int pos, BraidNode* bn);

  void removeLeftHandle(NODE node, std::set<NODE>& to_check, std::list<LinkedBraidStructureTransform>* lst);
  void removeRightHandle(NODE node, std::set<NODE>& to_check, std::list<LinkedBraidStructureTransform>* lst);

  LinkedBraidStructureTransform removeNode(BraidNode* bn, int pos);

  BraidNode* insertBackRight(BraidNode* bn, int pos, bool type);
  BraidNode* insertBackLeft(BraidNode* bn, int pos, bool type);
  BraidNode* insert(const LinkedBraidStructureTransform& lbst);

  void processTree(int al, BraidNode* node, Word& result) const;

  void clearLinks() const;

private:
  //! Number of generators!!!
  int the_index_;

  std::vector<BraidNode*> front_nodes_;
  std::vector<BraidNode*> back_nodes_;
  std::map<int, BraidNode> the_nodes_;

  uint64_t max_node_number_;
};

#endif // CRAG_LINKED_BRAID_STRUCTURE_H
