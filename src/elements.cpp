#include "elements.hpp"

#include <boost/bimap.hpp>


namespace GQCG {
namespace elements {


/*
 *  STATIC VARIABLES
 */

/**
 *  Create a vector of bimap values, which provides .begin() and .end() iterators
 *
 *  Adapted from https://stackoverflow.com/a/20290421/7930415
 */
std::vector<boost::bimap<std::string, size_t>::value_type> elements_list {
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


/**
 *  Use the .begin() and .end() iterators of a std::vector to construct a boost::bimap
 *
 *  Adapted from https://stackoverflow.com/a/20290421/7930415
 */
static const boost::bimap<std::string, size_t> periodic_table (elements_list.begin(), elements_list.end());



/*
 * FUNCTION IMPLEMENTATIONS
 */

/**
 *  Given a @param symbol for the name of an element, @return its atomic number
 */
size_t elementToAtomicNumber(const std::string& symbol) {

    return periodic_table.left.at(symbol);
}


/**
 *  Given an @param atomic_number, @return the symbol of the corresponding element
 */
const std::string& atomicNumberToElement(size_t atomic_number) {

    return periodic_table.right.at(atomic_number);
}


}  // namespace elements
}  // namespace GQCG
