#ifndef GQCG_ELEMENTS_HPP
#define GQCG_ELEMENTS_HPP


#include <string>


namespace GQCG {
namespace elements {


/**
 *  Given a @param symbol for the name of an element, @return its atomic number
 */
size_t elementToAtomicNumber(const std::string& symbol);


/**
 *  Given an @param atomic_number, @return the symbol of the corresponding element
 */
const std::string& atomicNumberToElement(size_t atomic_number);



}  // namespace elements
}  // namespace GQCG


#endif  // GQCG_ELEMENTS_HPP
