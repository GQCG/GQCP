#include "elements.hpp"

#include <map>


namespace GQCG {
namespace elements {


/**
 *  Given a @param symbol for the name of an element, @return its atomic number
 */
size_t element_to_atomic_number(const std::string& symbol) {

    std::map<std::string, size_t> element_to_Z_map = {  // Z is the atomic number
        {"H",  1},
        {"He", 2},
        {"Li", 3},
        {"Be", 4},
        {"B",  5},
        {"C",  6},
        {"N",  7},
        {"O",  8},
        {"F",  9},
        {"Ne", 10},
        {"Na", 11},
        {"Mg", 12},
        {"Al", 13},
        {"Si", 14},
        {"P",  15},
        {"S",  16},
        {"Cl", 17},
        {"Ar", 18}
    };

    return element_to_Z_map[symbol];
};



}  // namespace elements
}  // namespace GQCG
