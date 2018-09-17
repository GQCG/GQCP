#ifndef GQCG_ADDRESSINGSCHEME_HPP
#define GQCG_ADDRESSINGSCHEME_HPP


#include "common.hpp"



namespace GQCG {


/**
 *  An implementation of the addressing scheme from
 *      Molecular Electronic-Structure Theory (August 2000) by Trygve Helgaker, Poul Jorgensen, and Jeppe Olsen
 */
class AddressingScheme {
private:
    const size_t K;  // number of spatial orbitals
    const size_t N;  // number of electrons

    Matrixu vertex_weights;


public:
    // CONSTRUCTORS
    /**
     *  Constructor based on a number of spatial orbitals @param K and number electrons @param N
     */
    AddressingScheme(size_t K, size_t N);


    // GETTERS
    Matrixu get_vertex_weights() const { return this->vertex_weights; }
    size_t get_vertex_weights(size_t p, size_t m) const { return this->vertex_weights[p][m]; }
    size_t get_K() const { return this->K; }
    size_t get_N() const { return this->N; }
};


}  // namespace GQCG



#endif  // GQCG_ADDRESSINGSCHEME_HPP
