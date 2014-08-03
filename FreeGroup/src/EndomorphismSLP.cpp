/*
 * \file EndomorphismSLP.cpp
 *
 * 
 */

#include "EndomorphismSLP.h"

namespace crag {

EndomorphismSLP EndomorphismSLP::inverse() const {
  if (images_.size() == 0) //identity
    return *this;
  if (images_.size() > 1)
    throw std::invalid_argument("Unsupported endomorphism with more than one non-trivial terminal image!");

  auto img = images_.begin();
  TerminalSymbol symbol = img->first;
  slp::Vertex image = img->second;
  if (image.height() > 2)
    throw std::invalid_argument("Unsupported endomorphism with unsupported slp height > 2!");

  if (image.height() == 1) {//inverter
    return *this;
  }

  TerminalVertex left(image.left_child());
  TerminalVertex right(image.right_child());

  TerminalSymbol left_symbol = left.terminal_symbol();
  TerminalSymbol right_symbol = right.terminal_symbol();

  if (! (left_symbol == symbol || right_symbol == symbol))
    throw std::invalid_argument("Unsupported endomorphism not mapping the symbol to the product of another one and itself!");

  if (left_symbol == symbol) {
    return right_multiplier(symbol, -right_symbol);
  } else {
    return left_multiplier(-left_symbol, symbol);
  }
}

bool EndomorphismSLP::operator==(const EndomorphismSLP& a) const {
  if (this == &a)
    return true;
  if (non_trivial_images_num() != a.non_trivial_images_num())
    return false;
  auto images = non_trivial_images_range();
  auto a_images = a.non_trivial_images_range();
  (void) images; (void) a_images; //TODO:check this place

  //checking that the same letters have non-trivial images
  auto img_iterator = images_.begin();
  auto a_img_iterator = a.images_.begin();
  while (img_iterator != images_.end()) {
    if (img_iterator->first != a_img_iterator->first)
      return false;
    ++img_iterator;
    ++a_img_iterator;
  }

  //checking that the images are the same
  //todo: rewrite this method if it is used
  slp::MatchingTable mt;
  img_iterator = images_.begin();
  a_img_iterator = a.images_.begin();
  while (img_iterator != images_.end()) {
    slp::VertexWord word(img_iterator->second);
    slp::VertexWord a_word(a_img_iterator->second);
    if (!word.is_equal_to(a_word, &mt))
      return false;
    ++img_iterator;
    ++a_img_iterator;
  }

  return true;
}

EndomorphismSLP& EndomorphismSLP::operator*=(const EndomorphismSLP& a) {
  std::unordered_map<slp::Vertex, slp::Vertex> new_vertices;//a's vertices to new vertices correspondence

  for (const auto& root_entry: a.images_) {//mapping vertices of #a to new ones
    slp::map_vertices(root_entry.second, &new_vertices,
                      std::bind(&EndomorphismSLP::map_vertex, *this, std::placeholders::_1, std::placeholders::_2));
  }

  //replacing roots
  std::map<TerminalSymbol, slp::Vertex> new_images;
  for (const auto& root_entry: a.images_) {
    auto new_root = new_vertices.find(root_entry.second)->second;
    new_images.insert(std::make_pair(root_entry.first, new_root));
  }
  //adding images that were not inserted
  for (const auto& root_entry: images_) {
    if (new_images.find(root_entry.first) == new_images.end())//it was not mapped by a
      new_images.insert(root_entry);
  }

  swap(images_, new_images);
  return *this;
}

EndomorphismSLP EndomorphismSLP::conjugate_with(const AutomorphismDescription<EndomorphismSLP>& conjugator) const {
  return conjugator() * (*this) * conjugator.inverse();
}


slp::Vertex EndomorphismSLP::map_vertex(const slp::Vertex& vertex, const std::unordered_map<slp::Vertex, slp::Vertex>& images) const {
  auto item = images.find(vertex.negate());
  if (item != images.end()) {//already mapped inverse
    return item->second.negate();
  }

  if (!vertex)
    return vertex;//Mapping null vertex

  if (vertex.height() == 1) {//the vertex is terminal so map it to our corresponding root
    TerminalSymbol symbol = TerminalVertex(vertex).terminal_symbol();
    bool is_positive = is_positive_terminal_symbol(symbol);
    TerminalSymbol positive_symbol = is_positive ? symbol : -symbol;
    slp::Vertex v = image(positive_symbol);
    if (TerminalVertex(positive_symbol) == v)//if id map
      return vertex;
    return is_positive ? v : v.negate();
  } else {//for a nonterminal we already processed its children because postorder inspector
    auto left_val = images.find(vertex.left_child());
    auto right_val = images.find(vertex.right_child());
    assert(left_val != images.end() && right_val != images.end());
    const slp::Vertex& left = left_val->second;
    const slp::Vertex& right = right_val->second;
    if (left == vertex.left_child() && right == vertex.right_child()) //if children were not copied, then we should not copy vertex
      return vertex;
    return slp::NonterminalVertex(left, right);
  }
}

void EndomorphismSLP::save_to(std::ostream* out) const {
  long vertex_num = 0;

  std::unordered_map<size_t, std::pair<long, long>> non_terminals;
  std::unordered_map<size_t, TerminalSymbol> terminals;

  //we save the order in wich vertices occur during postorder inspection because we want to save them in this
  //order so it would be easier to reconstruct the automorphism later
  std::vector<long> non_terminals_order;
  std::vector<long> terminals_order;

  auto processor = [&] (const slp::Vertex& vertex, std::unordered_map<slp::Vertex, long>& mapped_images) {
    auto item = mapped_images.find(vertex.negate());
    if (item != mapped_images.end()) {//inverse was already visited
      return -item->second; //negate_index
    }

    ++vertex_num;
    if (vertex.height() == 1) {//the vertex is terminal
      TerminalSymbol symbol = TerminalVertex(vertex).terminal_symbol();
      bool is_positive = symbol > 0;
      const TerminalSymbol positive_symbol = is_positive ? symbol : - symbol;
      terminals.insert(std::make_pair(vertex_num, positive_symbol));
      terminals_order.push_back(vertex_num);
      if (is_positive) {
        mapped_images.insert(std::make_pair(vertex.negate(), -vertex_num));
        return vertex_num;
      } else {
        mapped_images.insert(std::make_pair(vertex.negate(), vertex_num));
        return -vertex_num;
      }
    } else {//nonterminal
      size_t left_val = mapped_images.find(vertex.left_child())->second;
      size_t right_val = mapped_images.find(vertex.right_child())->second;
      non_terminals.insert(std::make_pair(vertex_num, std::make_pair(left_val, right_val)));
      non_terminals_order.push_back(vertex_num);
    }
    return vertex_num;
  };

  std::unordered_map<slp::Vertex, long> vertex_numbers;
  for (const auto& root_entry: images_) {
    slp::map_vertices(root_entry.second, &vertex_numbers,
                      processor);
  }

  //writing
  //"number of nontrivial images" "number of terminals" "number of non-terminals"
  *out << images_.size() << " " << terminals.size() << " " << non_terminals.size() << std::endl;

  //writing terminal symbols
  //"terminal vertex index" "terminal symbol"
  for (long terminal_vertex_index: terminals_order) {
    const auto& terminal = *(terminals.find(terminal_vertex_index));
    *out << terminal.first << " " << terminal.second << std::endl;
  }

  //writing non-terminals
  //"vertex index" "left child index" "right child index"
  for (size_t non_terminal_vertex_index: non_terminals_order) {
    const auto& non_terminal = *(non_terminals.find(non_terminal_vertex_index));
    *out << non_terminal.first << " " << non_terminal.second.first << " " << non_terminal.second.second << std::endl;
  }

  //writing roots
  //"terminal symbol" "vertex index"
  for (const auto& root_entry: images_)
    *out << root_entry.first << " " << vertex_numbers.find(root_entry.second)->second << std::endl;
}

EndomorphismSLP EndomorphismSLP::load_from(std::istream* in) {
  size_t roots_num;
  size_t terminals_num;
  size_t non_terminals_num;
  *in >> roots_num >> terminals_num >> non_terminals_num;

  std::unordered_map<long, slp::Vertex> vertices;
  for (size_t i = 0; i < terminals_num; ++i) {
    long index;
    TerminalSymbol image;
    *in >> index >> image;
    in->ignore();
    vertices.insert(std::make_pair(index, TerminalVertex(image)));
  }

  auto get_vertex = [&vertices] (long index) {
    bool is_positive = index > 0;
    long positive_index = is_positive ? index : -index;
    const slp::Vertex& v = vertices.find(positive_index)->second;
    return is_positive ? v : v.negate();
  };

  for (size_t i = 0; i < non_terminals_num; ++i) {
    long index;
    long l_index;
    long r_index;
    *in >> index >> l_index >> r_index;
    in->ignore();
    const slp::Vertex& left = get_vertex(l_index);
    const slp::Vertex& right = get_vertex(r_index);
    vertices.insert(std::make_pair(index, slp::NonterminalVertex(left, right)));
  }

  EndomorphismSLP e;
  for (size_t i = 0; i < roots_num; ++i) {
    TerminalSymbol key;
    size_t index;
    *in >> key >> index;
    in->ignore();
    const slp::Vertex& root = get_vertex(index);
    e.images_.insert(std::make_pair(key, root));
  }
  return e;
}


void EndomorphismSLP::save_graphviz(std::ostream *p_out, const std::string& name) const {
  static const char* INDENT = "\t";

  std::ostream& out = *p_out;
  out << "digraph " << name << " {" << std::endl;
  out << "node [shape=point]" << std::endl;

  long vertex_num = 0;

  std::unordered_map<size_t, std::pair<long, long>> non_terminals;
  std::unordered_map<size_t, TerminalSymbol> terminals;
  std::unordered_map<TerminalSymbol, size_t> sym_to_vertex_num;

  auto processor = [&] (const slp::Vertex& vertex, std::unordered_map<slp::Vertex, long>& mapped_images) {
    auto item = mapped_images.find(vertex.negate());
    if (item != mapped_images.end()) {//inverse was already visited
      auto negate_index = item->second;
      return -negate_index;
    }

    ++vertex_num;
    if (vertex.height() == 1) {//the vertex is terminal
      TerminalSymbol symbol = TerminalVertex(vertex).terminal_symbol();
      bool is_positive = symbol > 0;
      const TerminalSymbol positive_symbol = is_positive ? symbol : - symbol;
      terminals.insert(std::make_pair(vertex_num, positive_symbol));
      sym_to_vertex_num.insert(std::make_pair(positive_symbol, vertex_num));
      if (is_positive) {
        mapped_images.insert(std::make_pair(vertex.negate(), -vertex_num));
        return vertex_num;
      } else {
        mapped_images.insert(std::make_pair(vertex.negate(), vertex_num));
        return -vertex_num;
      }
    } else {//nonterminal
      size_t left_val = mapped_images.find(vertex.left_child())->second;
      size_t right_val = mapped_images.find(vertex.right_child())->second;
      non_terminals.insert(std::make_pair(vertex_num, std::make_pair(left_val, right_val)));
    }
    return vertex_num;
  };

  std::unordered_map<slp::Vertex, long> vertex_numbers;
  for (const auto& root_entry: images_) {
    slp::map_vertices(root_entry.second, &vertex_numbers,
                      processor);
  }

  //writing root styles
  for (const auto& root_entry: images_) {
    out << INDENT << "\"i" << root_entry.first << "\" [shape=plaintext, label=\"" << root_entry.first << "\"];" << std::endl;

    out << INDENT << "\"i" << root_entry.first << "\" -> ";
    const slp::Vertex& img = root_entry.second;
    if (img.height() <= 1) {
      TerminalSymbol symbol = TerminalVertex(img).terminal_symbol();
      bool is_positive = symbol > 0;
      TerminalSymbol positive_symbol = is_positive ? symbol : -symbol;
      const auto& terminal = sym_to_vertex_num.find(positive_symbol);
      if (terminal != sym_to_vertex_num.end()) {
        if (is_positive) {
          out << terminal->second << ";" << std::endl;
        } else {
          out << terminal->second << "[label=\"-\"];" << std::endl;
        }
      } else {
        out << "\"" << symbol << "\"[shape=plaintext];" << std::endl;
      }
    } else {
      auto non_terminal_index = vertex_numbers.find(img)->second;
      if (non_terminal_index > 0) {
        out << non_terminal_index << " [style=dotted];" << std::endl;
      } else {
        out << (-non_terminal_index) << " [style=dotted, label=\"-\"];" << std::endl;
      }
    }
  }

  //writing terminals style
  for (const auto& pair: terminals) {
    out << INDENT << pair.first << " [shape=plaintext, label=\"" << pair.second << "\"];" << std::endl;
  }

  //writing non-terminals
  for (const auto& non_terminal: non_terminals) {
    size_t non_terminal_index = non_terminal.first;

    auto left_index = non_terminal.second.first;
    auto right_index = non_terminal.second.second;
    if (left_index > 0) {
      out << INDENT << non_terminal_index << " -> " << left_index << ";" << std::endl;
    } else {
      out << INDENT << non_terminal_index << " -> " << -left_index << "[label=\"-\"];" << std::endl;
    }
    if (right_index > 0){
      out << INDENT << non_terminal_index << " -> " << right_index << "[color=red,style=dashed];" << std::endl;
    } else {
      out << INDENT << non_terminal_index << " -> " << -right_index << "[color=red,style=dashed,label=\"-\"];" << std::endl;
    }
  }


  out << "}" << std::endl;
}

void EndomorphismSLP::print(std::ostream *p_out) const {
  std::ostream& out = *p_out;

  size_t vertex_num = 0;

  std::unordered_map<size_t, std::pair<size_t, size_t>> non_terminals;
  std::unordered_map<size_t, TerminalSymbol> terminals;

  auto processor = [&] (const slp::Vertex& vertex, const std::unordered_map<slp::Vertex, size_t>& mapped_images) {
    ++vertex_num;
    if (vertex.height() == 1) {//the vertex is terminal
      const TerminalSymbol& symbol = TerminalVertex(vertex).terminal_symbol();
      terminals.insert(std::make_pair(vertex_num, symbol));
    } else {//nonterminal
      size_t left_val = mapped_images.find(vertex.left_child())->second;
      size_t right_val = mapped_images.find(vertex.right_child())->second;
      non_terminals.insert(std::make_pair(vertex_num, std::make_pair(left_val, right_val)));
    }
    return vertex_num;
  };

  std::unordered_map<slp::Vertex, size_t> vertex_numbers;
  for (auto root_entry: images_) {
    slp::map_vertices(root_entry.second, &vertex_numbers,
                      processor);
  }

  std::unordered_map<size_t, TerminalSymbol> non_terminals_to_roots;

  //writing root styles
  for (auto root_entry: images_) {
    auto img = root_entry.second;
    if (img.height() <= 1) {
      auto symbol = TerminalVertex(img).terminal_symbol();
      out << "\"" << root_entry.first << "\" -> \"" << symbol << "\"" << std::endl;
    } else {
      non_terminals_to_roots.insert(std::make_pair(vertex_numbers.find(img)->second, root_entry.first));
    }
  }

  //writing non-terminals
  for (auto non_terminal: non_terminals) {
    size_t non_terminal_index = non_terminal.first;

    size_t left_index = non_terminal.second.first;
    size_t right_index = non_terminal.second.second;
    auto left_terminal = terminals.find(left_index);
    auto right_terminal = terminals.find(right_index);

    auto root_symbol = non_terminals_to_roots.find(non_terminal_index);
    if (root_symbol == non_terminals_to_roots.end()) {
      out << non_terminal_index;
    } else {
      out << "\"" << root_symbol->second << "\"";
    }

    out << " -> ";
    if (left_terminal == terminals.end()) {
      out << left_index;
    } else {
      out << "\"" << left_terminal->second << "\"";
    }
    out << " ";
    if (right_terminal == terminals.end()) {
      out << right_index;
    } else {
      out << "\"" << right_terminal->second << "\"";
    }
    out << std::endl;
  }
}


namespace internal {
  class Packer {//todo remove class
    public:
      static slp::Vertex pack_slps_into_one(const std::vector<slp::Vertex>& v) {
        std::size_t num = v.size();
        assert (num > 0);
        if (num == 1) {
          return v[0];
        }
        std::vector<slp::Vertex> packed_v;
        unsigned int packed_num = (num / 2) + (num % 2);
        packed_v.reserve(packed_num);
        for (std::size_t i = 0; i < num - 1; i += 2) {
          packed_v.push_back(slp::NonterminalVertex(v[i], v[i + 1]));
        }
        if (num % 2 == 1) {
          packed_v.push_back(v.back());
        }
        return pack_slps_into_one(packed_v);
      }
  };
}

EndomorphismSLP EndomorphismSLP::normal_form() const {
  if (non_trivial_images_num() == 0)
    return EndomorphismSLP();
  //we rewrite all vertices into a single SLP then find normal form
  //and then split into pieces according to original vertices lengths

  //end positions in the unifying SLP
  std::vector<std::pair<slp::TerminalSymbol, LongInteger> > end_positions;

  //vertices to combine into SLP
  std::vector<slp::Vertex> vertices;
  vertices.reserve(non_trivial_images_num());

  LongInteger length_accumulator(0);
  for_each_non_trivial_image([&] (const symbol_image_pair_type& pair) {
    const TerminalSymbol& k = pair.first;
    const slp::Vertex& v = pair.second;

    //fill positions
    length_accumulator += v.length();
    end_positions.push_back(std::make_pair(k, length_accumulator));

    //adding vertices to process
    vertices.push_back(v);
  });

  slp::Vertex slp = internal::Packer::pack_slps_into_one(vertices);
  auto nf = slp::recompression::normal_form(slp);

  //extracting slps
  EndomorphismSLP result;
  LongInteger begin(0);
  for (auto key_pos_pair: end_positions) {
    auto end = key_pos_pair.second;
    auto sub_slp = slp::get_sub_slp(nf, begin, end);
    result.images_.insert(std::make_pair(key_pos_pair.first, sub_slp));
    begin = end;
  }

  return result;
}

//-------------------------------------------------------------------------
// Implementation of helper functions.


unsigned int height(const EndomorphismSLP& e) {
  unsigned int h = 0;
  auto pick_max_height = [&h] (const EndomorphismSLP::symbol_image_pair_type& v) {
    const unsigned int v_h = v.second.height();
    if (v_h > h)
      h = v_h;
  };
  e.for_each_non_trivial_image(pick_max_height);
  return h;
}


unsigned int slp_vertices_num(const EndomorphismSLP& e) {
  std::unordered_set<slp::Vertex> visited_vertices;

  auto acceptor = [&visited_vertices] (const slp::inspector::InspectorTask& task) {
    return visited_vertices.count(task.vertex) == 0
        && visited_vertices.count(task.vertex.negate()) == 0;
  };

  auto inspect_root =[&acceptor,&visited_vertices] (const EndomorphismSLP::symbol_image_pair_type& v) {
    slp::Inspector<slp::inspector::Postorder, decltype(acceptor)> inspector(v.second, acceptor);
    while (!inspector.stopped()) {
      visited_vertices.insert(inspector.vertex());
      inspector.next();
    }
  };

  e.for_each_non_trivial_image(inspect_root);
  return visited_vertices.size();
}

unsigned int slp_unique_images_length_num(const EndomorphismSLP& e) {
  std::set<LongInteger> visited_vertices;//map sum images length -> vertex

  auto acceptor = [&visited_vertices] (const slp::inspector::InspectorTask& task) {
    return visited_vertices.find(task.vertex.length()) == visited_vertices.end();
  };

  auto inspect_root =[&acceptor,&visited_vertices] (const EndomorphismSLP::symbol_image_pair_type& v) {
    slp::Inspector<slp::inspector::Postorder, decltype(acceptor)> inspector(v.second, acceptor);
    while (!inspector.stopped()) {
      visited_vertices.insert(inspector.vertex().length());
      inspector.next();
    }
  };

  e.for_each_non_trivial_image(inspect_root);
  return visited_vertices.size();
}


std::map<slp::TerminalSymbol, LongInteger> images_length(const EndomorphismSLP& e) {
  std::map<slp::TerminalSymbol, LongInteger> key_to_lengths;
  auto add_length = [&key_to_lengths] (const EndomorphismSLP::symbol_image_pair_type& pair) {
   auto key = pair.first;
   slp::Vertex v = pair.second;
   key_to_lengths.insert(std::make_pair(key, v.length()));
  };
  e.for_each_non_trivial_image(add_length);
  return key_to_lengths;
}


}  // namespace crag


