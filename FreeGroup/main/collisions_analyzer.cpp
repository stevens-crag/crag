
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "EndomorphismSLP.h"

const std::vector<crag::EndomorphismSLP> load(std::istream* p_in) {
  std::size_t num;
  *p_in >> num;
  std::vector<crag::EndomorphismSLP> morphisms;
  morphisms.reserve(num);
  for (std::size_t i = 0; i < num; ++i) {
    morphisms.push_back(crag::EndomorphismSLP::load_from(p_in));
  }
  return morphisms;
}

void to_json(std::ostream* p_out, const crag::EndomorphismSLP& e) {
  std::ostream& out = *p_out;
  out << "{" << std::endl;

  out << "repr : " << "\"";
  e.save_to(p_out);
  out << "\"," << std::endl;

  out << "dot : " << "\"";
  e.save_graphviz(p_out);
  out << "\"," << std::endl;

  out << "height : " << crag::height(e) << "," << std::endl;
  out << "vertices num : " << crag::slp_vertices_num(e) << std::endl;

  out << "}";
}

template<typename Iterator>
void iter_to_json(std::ostream *p_out, Iterator begin, Iterator end) {
  std::ostream& out = *p_out;
  out << "[" << std::endl;
  for(;begin != end; ++begin) {
    to_json(p_out, *begin);
    out << "," << std::endl;
  }
  out << "]" << std::endl;
}


int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Not enough parameters." << std::endl;
//    print_usage();
    return 1;
  }
  std::string filename(argv[1]);
  std::cout << filename << std::endl;
  std::ifstream file_in(filename);
  auto result = load(&file_in);
  if (argc == 2) {
    iter_to_json(&std::cout, result.cbegin(), result.cend());
  } else {
    std::ofstream out(std::string(argv[2]));
    iter_to_json(&out, result.cbegin(), result.cend());
  }
}
