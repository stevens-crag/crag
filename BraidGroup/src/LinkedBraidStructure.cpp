// Copyright (C) 2005 Alexander Ushakov

#include "LinkedBraidStructure.h"

std::ostream& operator<<(std::ostream& os, const BraidNode& bn) {
  os << &bn << "{ " << bn.type << "| " << bn.left << ", " << bn.ahead << ", " << bn.right << ", " << bn.back_left
     << ", " << bn.back << ", " << bn.back_right << " }";

  return os;
}

LinkedBraidStructureTransform LinkedBraidStructure::make_EraseTransform(BraidNode* bn, int64_t pos) const {
  return LinkedBraidStructureTransform(
      bn->theNumber, pos, LinkedBraidStructureTransform::ERASED, bn->type, bn->left ? bn->left->theNumber : -1,
      bn->ahead ? bn->ahead->theNumber : -1, bn->right ? bn->right->theNumber : -1,
      bn->back_left ? bn->back_left->theNumber : -1, bn->back ? bn->back->theNumber : -1,
      bn->back_right ? bn->back_right->theNumber : -1);
}

LinkedBraidStructureTransform LinkedBraidStructure::make_AddTransform(BraidNode* bn, int64_t pos) const {
  return LinkedBraidStructureTransform(bn->theNumber, pos, LinkedBraidStructureTransform::ADDED, bn->type);
}

LinkedBraidStructureTransform LinkedBraidStructure::make_ChangeType(BraidNode* bn, int64_t pos) const {
  return LinkedBraidStructureTransform(bn->theNumber, pos, LinkedBraidStructureTransform::CHANGE_TYPE, bn->type);
}

bool LinkedBraidStructure::operator<(const LinkedBraidStructure& lbs) const {
  if (size() < lbs.size()) {
    return true;
  }

  if (size() > lbs.size()) {
    return false;
  }

  if (the_index_ < lbs.the_index_) {
    return true;
  }

  if (the_index_ > lbs.the_index_) {
    return false;
  }

  return translateIntoWord() < lbs.translateIntoWord();
}

LinkedBraidStructure::LinkedBraidStructure(size_t N)
    : max_node_number_(0)
    , the_index_(N)
    , back_nodes_(std::vector<BraidNode*>(N, nullptr))
    , front_nodes_(std::vector<BraidNode*>(N, nullptr)) {}

LinkedBraidStructure::LinkedBraidStructure(size_t N, const Word& w)
    : max_node_number_(0)
    , the_index_(N)
    , back_nodes_(std::vector<BraidNode*>(N, nullptr))
    , front_nodes_(std::vector<BraidNode*>(N, nullptr)) {
  for (auto w_it = w.begin(); w_it != w.end(); ++w_it) {
    push_back(*w_it);
  }
}

LinkedBraidStructure& LinkedBraidStructure::operator=(const LinkedBraidStructure& LBS) {
  // cout << "   #" << endl;

  the_index_ = LBS.the_index_;
  the_nodes_ = LBS.the_nodes_;
  max_node_number_ = LBS.max_node_number_;
  back_nodes_ = LBS.back_nodes_;
  front_nodes_ = LBS.front_nodes_;

  auto n_it1 = LBS.the_nodes_.begin();
  auto n_it2 = the_nodes_.begin();

  for (; n_it1 != LBS.the_nodes_.end(); ++n_it1, ++n_it2) {
    (*n_it1).second.link = &(*n_it2).second;
  }

  n_it2 = the_nodes_.begin();

  for (; n_it2 != the_nodes_.end(); ++n_it2) {
    BraidNode& cur = (*n_it2).second;

    if (cur.left) {
      cur.left = cur.left->link;
    }

    if (cur.ahead) {
      cur.ahead = cur.ahead->link;
    }

    if (cur.right) {
      cur.right = cur.right->link;
    }

    if (cur.back_left) {
      cur.back_left = cur.back_left->link;
    }

    if (cur.back) {
      cur.back = cur.back->link;
    }

    if (cur.back_right) {
      cur.back_right = cur.back_right->link;
    }
  }

  for (size_t i = 0; i < the_index_; ++i) {
    if (back_nodes_[i]) {
      back_nodes_[i] = back_nodes_[i]->link;
    }

    if (front_nodes_[i]) {
      front_nodes_[i] = front_nodes_[i]->link;
    }
  }

  LBS.clearLinks();

  return *this;
}

LinkedBraidStructure::LinkedBraidStructure(const LinkedBraidStructure& LBS)
    : the_index_(LBS.the_index_)
    , the_nodes_(LBS.the_nodes_)
    , max_node_number_(LBS.max_node_number_)
    , back_nodes_(LBS.back_nodes_)
    , front_nodes_(LBS.front_nodes_) {
  // cout << "   &" << endl;
  auto n_it1 = LBS.the_nodes_.begin();
  auto n_it2 = the_nodes_.begin();

  for (; n_it1 != LBS.the_nodes_.end(); ++n_it1, ++n_it2) {
    (*n_it1).second.link = &(*n_it2).second;
  }

  n_it2 = the_nodes_.begin();

  for (; n_it2 != the_nodes_.end(); ++n_it2) {
    BraidNode& cur = (*n_it2).second;

    if (cur.left) {
      cur.left = cur.left->link;
    }

    if (cur.ahead) {
      cur.ahead = cur.ahead->link;
    }

    if (cur.right) {
      cur.right = cur.right->link;
    }

    if (cur.back_left) {
      cur.back_left = cur.back_left->link;
    }

    if (cur.back) {
      cur.back = cur.back->link;
    }

    if (cur.back_right) {
      cur.back_right = cur.back_right->link;
    }
  }

  for (size_t i = 0; i < the_index_; ++i) {
    if (back_nodes_[i]) {
      back_nodes_[i] = back_nodes_[i]->link;
    }

    if (front_nodes_[i]) {
      front_nodes_[i] = front_nodes_[i]->link;
    }
  }

  LBS.clearLinks();
}

LinkedBraidStructureTransform LinkedBraidStructure::push_back(int g) {
  int ag = std::abs(g);
  BraidNode* left = (ag > 1) ? back_nodes_[ag - 2] : nullptr;
  BraidNode* ahead = back_nodes_[ag - 1];
  BraidNode* right = (ag < the_index_) ? back_nodes_[ag] : nullptr;

  BraidNode& newNode = the_nodes_[max_node_number_] = BraidNode(max_node_number_, g > 0);

  max_node_number_++;

  newNode.ahead = ahead;

  if (ahead) {
    ahead->back = &newNode;
  }

  if (!ahead || ahead && ahead->back_left) {
    newNode.left = left;

    if (left) {
      left->back_right = &newNode;
    }
  }

  if (!ahead || ahead && ahead->back_right) {
    newNode.right = right;

    if (right) {
      right->back_left = &newNode;
    }
  }

  if (!front_nodes_[ag - 1]) {
    front_nodes_[ag - 1] = &newNode;
  }

  back_nodes_[ag - 1] = &newNode;

  return make_AddTransform(&newNode, ag - 1);
}

LinkedBraidStructureTransform LinkedBraidStructure::push_front(int g) {
  int ag = std::abs(g);

  BraidNode* left = (ag > 1) ? front_nodes_[ag - 2] : 0;
  BraidNode* ahead = front_nodes_[ag - 1];
  BraidNode* right = (ag < the_index_) ? front_nodes_[ag] : 0;

  BraidNode& newNode = the_nodes_[max_node_number_] = BraidNode(max_node_number_, g > 0);

  max_node_number_++;

  newNode.back = ahead;

  if (ahead) {
    ahead->ahead = &newNode;
  }

  if (!ahead || ahead && ahead->left) {
    newNode.back_left = left;

    if (left) {
      left->right = &newNode;
    }
  }

  if (!ahead || ahead && ahead->right) {
    newNode.back_right = right;

    if (right) {
      right->left = &newNode;
    }
  }

  if (!back_nodes_[ag - 1]) {
    back_nodes_[ag - 1] = &newNode;
  }

  front_nodes_[ag - 1] = &newNode;

  return make_AddTransform(&newNode, ag - 1);
}

void LinkedBraidStructure::removeLeftHandles(list<LinkedBraidStructureTransform>* result) {
  // set of nodes to check (1st value = number of uppercrossings, 2nd value =
  // -position in the Linked Structure)
  node_set to_check;

  // form the initial set of nodes to check (all nodes from the Linked
  // Structure)
  for (size_t i = 0; i < the_index_; ++i) {
    for (BraidNode* c = front_nodes_[i]; c != back_nodes_[i]; c = c->back) {
      const auto weight = checkIfStartsLeftHandle(i, c);

      if (weight != -1) {
        c->weight = weight;
        to_check.insert(NODE(weight, i, c));
      }
    }
  }

  // check the nodes one by one starting from the rightmost
  for (; !to_check.empty();) {
    // a. choose the node and remove it from to_check
    NODE cur_node = *to_check.begin();
    to_check.erase(to_check.begin());

    // b. remove the handle
    const auto weight = checkIfStartsLeftHandle(std::get<1>(cur_node), std::get<2>(cur_node));

    if (weight != -1) {
      removeLeftHandle(cur_node, to_check, result);
    }
  }
}

int64_t LinkedBraidStructure::checkIfStartsLeftHandle(int64_t pos, BraidNode* bn) {
  BraidNode* back = bn->back;

  // if there is no more x_i
  if (!back) {
    return -1;
  }

  // if there is an obstacle for handle
  if (bn->back_right) {
    return -1;
  }

  // if crossings have the same orientation
  if (bn->type == back->type) {
    return -1;
  }

  // compute the number of upper crossings
  BraidNode* l1 = bn->back_left;
  BraidNode* l2 = back->left;

  int64_t counter = 0;
  bool crossingType;

  for (; l1; l1 = l1->back) {
    if (counter == 0) {
      crossingType = l1->type;
    } else {
      if (crossingType != l1->type) {
        return -1;
      }
    }

    ++counter;

    if (l1 == l2) {
      break;
    }
  }

  return counter;
}

namespace {
struct PairCmp {
  using pair_t = std::pair<int64_t, BraidNode*>;
  using int64_pair = std::pair<int64_t, int64_t>;

  bool operator()(const pair_t& lhs, const pair_t& rhs) const {
    return std::less<int64_pair>()(toInt64Pair(lhs), toInt64Pair(rhs));
  }

  int64_pair toInt64Pair(const pair_t& p) const {
    return std::make_pair(p.first, p.second->theNumber);
  }
};
} // namespace

void LinkedBraidStructure::removeLeftHandle(
    NODE node, node_set& to_check, std::list<LinkedBraidStructureTransform>* lst) {
  const auto pos = std::get<1>(node);

  BraidNode* n1 = std::get<2>(node);
  BraidNode* n2 = n1->back;

  bool type = n1->type;

  BraidNode* l1 = n1->back_left;
  BraidNode* l2 = n2->left;

  // Removal of a handle can introduce new handles. Here we store some nodes to
  // check. A few will be added later.
  std::set<std::pair<int64_t, BraidNode*>, PairCmp> to_check2;

  if (n1->ahead) {
    to_check2.insert(std::make_pair(pos, n1->ahead));
  }

  if (n1->left) {
    to_check2.insert(std::make_pair(pos - 1, n1->left));
  }

  for (BraidNode* cn = n1; cn; cn = cn->ahead) {
    if (!cn->right) {
      continue;
    }

    to_check2.insert(std::make_pair(pos + 1, cn->right));

    break;
  }

  // A. process left nodes

  for (; l1;) {
    BraidNode* l3 = l1->back;
    BraidNode* new_node = insertBackRight(l1, pos, l1->type);
    BraidNode* new_node2 = insertBackLeft(new_node, pos - 1, type);

    if (lst) {
      lst->push_back(make_AddTransform(new_node, pos));
      lst->push_back(make_AddTransform(new_node2, pos - 1));
      lst->push_back(make_ChangeType(l1, pos));
    }

    l1->type = !type;
    to_check2.insert(std::make_pair(pos - 1, new_node2));

    if (l1 == l2) {
      to_check2.insert(std::make_pair(pos, new_node));
      break;
    }

    l1 = l3;
  }

  // B. Remove "boundary" crossings that formed a handle.
  to_check.erase(NODE(n1->weight, pos, n1));
  to_check.erase(NODE(n2->weight, pos, n2));

  if (lst) {
    lst->push_back(removeNode(n1, pos));
    lst->push_back(removeNode(n2, pos));
  } else {
    removeNode(n1, pos);
    removeNode(n2, pos);
  }

  // Check if new handles were introduced.
  for (auto c_it = to_check2.begin(); c_it != to_check2.end(); ++c_it) {
    const auto pos = c_it->first;
    BraidNode* bn = c_it->second;

    const auto weight = checkIfStartsLeftHandle(pos, bn);

    if (weight != -1) {
      // remove old version
      to_check.erase(NODE(bn->weight, pos, bn));
      // add a new
      bn->weight = weight;
      to_check.insert(NODE(weight, pos, bn));
    }
  }
}

void LinkedBraidStructure::removeRightHandles(list<LinkedBraidStructureTransform>* result) {
  // set of nodes to check (1st value = number of uppercrossings, 2nd value =
  // position in the Linked Structure)
  node_set to_check;

  // form the initial set of nodes to check (all nodes from the Linked
  // Structure)
  for (size_t i = 0; i < the_index_; ++i) {
    for (BraidNode* c = front_nodes_[i]; c != back_nodes_[i]; c = c->back) {
      const auto weight = checkIfStartsRightHandle(i, c);

      if (weight != -1) {
        c->weight = weight;
        to_check.insert(NODE(weight, i, c));
      }
    }
  }

  // check the nodes one by one starting from the rightmost
  for (; !to_check.empty();) {
    // choose the node and remove it from to_check
    NODE cur_node = *to_check.begin();
    to_check.erase(to_check.begin());

    // b. remove the handle
    const auto weight = checkIfStartsRightHandle(std::get<1>(cur_node), std::get<2>(cur_node));

    if (weight != -1) {
      removeRightHandle(cur_node, to_check, result);
    }
  }
}

int64_t LinkedBraidStructure::checkIfStartsRightHandle(int64_t pos, BraidNode* bn) {
  BraidNode* back = bn->back;

  // if there is no more x_i
  if (!back) {
    return -1;
  }

  // if there is an obstacle for handle
  if (bn->back_left) {
    return -1;
  }

  // if crossings have the same orientation
  if (bn->type == back->type) {
    return -1;
  }

  // compute the number of upper crossings
  BraidNode* l1 = bn->back_right;
  BraidNode* l2 = back->right;

  int64_t counter = 0;
  bool crossingType;

  for (; l1; l1 = l1->back) {
    if (counter == 0) {
      crossingType = l1->type;
    } else {
      if (crossingType != l1->type) {
        return -1;
      }
    }

    ++counter;

    if (l1 == l2) {
      break;
    }
  }

  return counter;
}

void LinkedBraidStructure::removeRightHandle(NODE node, node_set& to_check, list<LinkedBraidStructureTransform>* lst) {
  const auto pos = std::get<1>(node);

  BraidNode* n1 = std::get<2>(node);
  BraidNode* n2 = n1->back;

  bool type = n1->type;

  BraidNode* r1 = n1->back_right;
  BraidNode* r2 = n2->right;

  // Removal of a handle can introduce new handles. Here we store some nodes to
  // check. A few will be added later.
  std::set<std::pair<int64_t, BraidNode*>, PairCmp> to_check2;

  if (n1->ahead) {
    to_check2.insert(std::make_pair(pos, n1->ahead));
  }

  if (n1->right) {
    to_check2.insert(std::make_pair(pos + 1, n1->right));
  }

  for (BraidNode* cn = n1; cn; cn = cn->ahead) {
    if (!cn->left) {
      continue;
    }

    to_check2.insert(std::make_pair(pos - 1, cn->left));

    break;
  }

  // B. process right nodes
  for (; r1;) {
    // cout << "Process right node" << endl;
    BraidNode* r3 = r1->back;
    BraidNode* new_node = insertBackLeft(r1, pos, r1->type);
    BraidNode* new_node2 = insertBackRight(new_node, pos + 1, type);

    if (lst) {
      lst->push_back(make_AddTransform(new_node, pos));
      lst->push_back(make_AddTransform(new_node2, pos + 1));
      lst->push_back(make_ChangeType(r1, pos));
    }

    r1->type = !type;

    to_check2.insert(std::make_pair(pos + 1, new_node2));

    if (r1 == r2) {
      to_check2.insert(std::make_pair(pos, new_node));
      break;
    }

    r1 = r3;
  }

  // B. Remove "boundary" crossings that formed a handle.
  to_check.erase(NODE(n1->weight, pos, n1));
  to_check.erase(NODE(n2->weight, pos, n2));

  if (lst) {
    lst->push_back(removeNode(n1, pos));
    lst->push_back(removeNode(n2, pos));
  } else {
    removeNode(n1, pos);
    removeNode(n2, pos);
  }

  // Check if new handles were introduced.
  for (auto c_it = to_check2.begin(); c_it != to_check2.end(); ++c_it) {
    const auto pos = c_it->first;
    BraidNode* bn = c_it->second;

    const auto weight = checkIfStartsRightHandle(pos, bn);

    if (weight != -1) {
      // remove old version
      to_check.erase(NODE(bn->weight, pos, bn));
      // add a new
      bn->weight = weight;
      to_check.insert(NODE(weight, pos, bn));
    }
  }
}

LinkedBraidStructureTransform LinkedBraidStructure::removeNode(BraidNode* bn, int64_t pos) {
  LinkedBraidStructureTransform result = make_EraseTransform(bn, pos);

  // A. process left
  if (bn->left) {
    bn->left->back_right = bn->back_left ? nullptr : bn->back;
  }

  // B. process ahead
  if (bn->ahead) {
    bn->ahead->back = bn->back;

    if (!bn->ahead->back_left) {
      bn->ahead->back_left = bn->back_left;
    }

    if (!bn->ahead->back_right) {
      bn->ahead->back_right = bn->back_right;
    }
  }

  // C. process right
  if (bn->right) {
    bn->right->back_left = bn->back_right ? nullptr : bn->back;
  }

  // D. process back_left
  if (bn->back_left) {
    bn->back_left->right = bn->left ? nullptr : bn->ahead;
  }

  // E. process back
  if (bn->back) {
    bn->back->ahead = bn->ahead;

    if (!bn->back->left) {
      bn->back->left = bn->left;
    }

    if (!bn->back->right) {
      bn->back->right = bn->right;
    }
  }

  // F. process back_right
  if (bn->back_right) {
    bn->back_right->left = bn->right ? nullptr : bn->ahead;
  }

  // G. process front_nodes_ and back_nodes_
  if (front_nodes_[pos] == bn) {
    front_nodes_[pos] = bn->back;
  }

  if (back_nodes_[pos] == bn) {
    back_nodes_[pos] = bn->ahead;
  }

  // H. finally delet the memory allocated for the node
  the_nodes_.erase(bn->theNumber);

  return result;
}

BraidNode* LinkedBraidStructure::insertBackLeft(BraidNode* bn, int64_t pos, bool type) {
  // A. determine the "main" nodes around the new node
  BraidNode* back = bn->back;

  BraidNode* left = nullptr;

  for (BraidNode* c = bn; c != nullptr && !left; c = c->ahead) {
    if (c->left) {
      left = c->left;
    }
  }

  BraidNode* back_left = nullptr;

  for (BraidNode* c = bn; c != nullptr && !back_left; c = c->back) {
    if (c->back_left) {
      back_left = c->back_left;
    }
  }

  BraidNode* left_back_left;

  if (left) {
    left_back_left = left->back_left;
  } else if (back_left) {
    if (!back_left->left) {
      left_back_left = nullptr;
    } else {
      left_back_left = pos > 0 ? front_nodes_[pos - 1] : nullptr;
    }
  } else {
    left_back_left = pos > 0 ? front_nodes_[pos - 1] : nullptr;
  }

  BraidNode* back_left_right;

  if (back_left) {
    back_left_right = back_left->right != bn ? bn->back : nullptr;
  } else {
    back_left_right = back;
  }

  // B. create a new node and link it to other nodes
  BraidNode& newNode = the_nodes_[max_node_number_] =
      BraidNode(max_node_number_, type, nullptr, left, bn, left_back_left, back_left, back_left_right);

  max_node_number_++;

  bn->back_left = &newNode;

  if (left) {
    left->back = &newNode;
    left->back_left = nullptr;
  }

  if (back_left) {
    back_left->ahead = &newNode;
    back_left->right = back_left->right == bn ? nullptr : back_left->right;
  }

  if (left_back_left) {
    left_back_left->right = &newNode;
  }

  if (back_left_right) {
    back_left_right->left = &newNode;
  }

  if (!newNode.back) {
    back_nodes_[pos] = &newNode;
  }

  if (!newNode.ahead) {
    front_nodes_[pos] = &newNode;
  }

  return &newNode;
}

BraidNode* LinkedBraidStructure::insertBackRight(BraidNode* bn, int64_t pos, bool type) {
  // A. determine the "main" nodes around the new node
  BraidNode* back = bn->back;

  BraidNode* right = nullptr;

  for (BraidNode* c = bn; c != nullptr && !right; c = c->ahead) {
    if (c->right) {
      right = c->right;
    }
  }

  BraidNode* back_right = nullptr;

  for (BraidNode* c = bn; c != nullptr && !back_right; c = c->back) {
    if (c->back_right) {
      back_right = c->back_right;
    }
  }

  BraidNode* right_back_right;

  if (right) {
    right_back_right = right->back_right;
  } else if (back_right) {
    if (!back_right->right) {
      right_back_right = nullptr;
    } else {
      right_back_right = (pos + 1 < the_index_) ? front_nodes_[pos + 1] : nullptr;
    }
  } else {
    right_back_right = (pos + 1 < the_index_) ? front_nodes_[pos + 1] : nullptr;
  }

  BraidNode* back_right_left;
  if (back_right) {
    back_right_left = back_right->left != bn ? bn->back : nullptr;
  } else {
    back_right_left = back;
  }

  // B. create a new node and link it to other nodes
  BraidNode& newNode = the_nodes_[max_node_number_] =
      BraidNode(max_node_number_, type, bn, right, nullptr, back_right_left, back_right, right_back_right);

  max_node_number_++;

  bn->back_right = &newNode;

  if (right) {
    right->back = &newNode;
    right->back_right = nullptr;
  }

  if (back_right) {
    back_right->ahead = &newNode;
    back_right->left = back_right->left == bn ? nullptr : back_right->left;
  }

  if (back_right_left) {
    back_right_left->right = &newNode;
  }

  if (right_back_right) {
    right_back_right->left = &newNode;
  }

  // update back_nodes_ and front_nodes_
  if (!newNode.back) {
    back_nodes_[pos] = &newNode;
  }

  if (!newNode.ahead) {
    front_nodes_[pos] = &newNode;
  }

  return &newNode;
}

BraidNode* LinkedBraidStructure::insert(const LinkedBraidStructureTransform& lbst) {
  BraidNode& newNode = the_nodes_[lbst.theNumber] = BraidNode(
      lbst.theNumber, lbst.type, lbst.left == -1 ? nullptr : &the_nodes_[lbst.left],
      lbst.ahead == -1 ? nullptr : &the_nodes_[lbst.ahead], lbst.right == -1 ? nullptr : &the_nodes_[lbst.right],
      lbst.back_left == -1 ? nullptr : &the_nodes_[lbst.back_left], lbst.back == -1 ? nullptr : &the_nodes_[lbst.back],
      lbst.back_right == -1 ? nullptr : &the_nodes_[lbst.back_right]);

  // A. determine the "main" nodes around the new node
  if (newNode.left) {
    newNode.left->back_right = &newNode;
  }

  if (newNode.ahead) {
    newNode.ahead->back = &newNode;

    if (!newNode.left) {
      newNode.ahead->back_left = nullptr;
    }

    if (!newNode.right) {
      newNode.ahead->back_right = nullptr;
    }
  }

  if (newNode.right) {
    newNode.right->back_left = &newNode;
  }

  if (newNode.back_left) {
    newNode.back_left->right = &newNode;
  }

  if (newNode.back) {
    newNode.back->ahead = &newNode;

    if (!newNode.back_left) {
      newNode.back->left = nullptr;
    }

    if (!newNode.back_right) {
      newNode.back->right = nullptr;
    }
  }

  if (newNode.back_right) {
    newNode.back_right->left = &newNode;
  }

  // B. update back_nodes_ and front_nodes_
  if (!newNode.back) {
    back_nodes_[lbst.thePosition] = &newNode;
  }

  if (!newNode.ahead) {
    front_nodes_[lbst.thePosition] = &newNode;
  }

  return &newNode;
}

void LinkedBraidStructure::processTree(int al, BraidNode* node, Word& result) const {
  if (node->back_left && !node->back_left->link) {
    processTree(al - 1, node->back_left, result);
  }

  if (node->back_right && !node->back_right->link) {
    processTree(al + 1, node->back_right, result);
  }

  if (node->back && !node->back->link) {
    processTree(al, node->back, result);
  }

  result.push_front(node->type ? al : -al);

  node->link = node;
}

Word LinkedBraidStructure::translateIntoWord() const {
  Word result;

  for (size_t i = 0; i < the_index_; ++i) {
    if (front_nodes_[i] && !front_nodes_[i]->link) {
      processTree(i + 1, front_nodes_[i], result);
    }
  }

  clearLinks();

  return result;
}

void LinkedBraidStructure::clearLinks() const {
  auto n_it = the_nodes_.begin();

  for (; n_it != the_nodes_.end(); ++n_it) {
    (*n_it).second.link = nullptr;
  }
}

void LinkedBraidStructure::clear() {
  the_nodes_.clear();

  front_nodes_ = std::vector<BraidNode*>(the_index_, nullptr);
  back_nodes_ = std::vector<BraidNode*>(the_index_, nullptr);
}

void LinkedBraidStructure::undo(const std::list<LinkedBraidStructureTransform>& lbst_seq) {
  std::list<LinkedBraidStructureTransform>::const_iterator s_it = lbst_seq.end();

  for (; s_it != lbst_seq.begin();) {
    undo(*--s_it);
  }
}

void LinkedBraidStructure::undo(const LinkedBraidStructureTransform& lbst) {
  if (lbst.theTransform == LinkedBraidStructureTransform::ADDED) {
    removeNode(&the_nodes_[lbst.theNumber], lbst.thePosition);
  } else if (lbst.theTransform == LinkedBraidStructureTransform::ERASED) {
    insert(lbst);
  } else if (lbst.theTransform == LinkedBraidStructureTransform::CHANGE_TYPE) {
    the_nodes_[lbst.theNumber].type = lbst.type;
  }
}

bool areEqualBraids(size_t n, const Word& lhs, const Word& rhs) {
  return isTrivialBraid(n, -lhs * rhs);
}

bool isTrivialBraid(size_t n, const Word& w) {
  LinkedBraidStructure lbs(n - 1, w);
  lbs.removeLeftHandles();

  return lbs.size() == 0;
}
