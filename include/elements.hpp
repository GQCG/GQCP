#ifndef GQCG_ELEMENTS_HPP
#define GQCG_ELEMENTS_HPP


#include <string>


namespace GQCG {
namespace elements {


/**
 *  Given a @param symbol for the name of an element, @return its atomic number
 */
size_t element_to_atomic_number(const std::string& symbol);



}  // namespace elements
}  // namespace GQCG


#endif  // GQCG_ELEMENTS_HPP
