// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#include "Basis/Shell.hpp"

#include "utilities/miscellaneous.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param l                the angular momentum of the shell (x + y + z)
 *  @param atom             the atom on which the shell is centered
 *  @param exponents        the exponents, which are shared for every contraction
 *  @param coefficients     the contraction coefficients
 */
Shell::Shell(size_t l, const Atom& atom, const std::vector<double>& exponents, const std::vector<double>& coefficients) :
    l (l),
    atom (atom),
    exponents (exponents),
    coefficients (coefficients)
{
    if (exponents.size() != coefficients.size()) {
        throw std::invalid_argument("Shell(size_t, Atom, std::vector<double>, std::vector<double>): the exponents and contraction coefficients must match in size.");
    }
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @return the number of basis functions that are in this shell
 */
size_t Shell::numberOfBasisFunctions() const {
    return (this->l + 1) * (this->l + 2) / 2;  // Cartesian shell
}


/**
 *  @return the basis functions that are represented by this shell
 */
std::vector<BasisFunction> Shell::basisFunctions() const {

    std::vector<BasisFunction> bfs;  // basis functions
    bfs.reserve(this->numberOfBasisFunctions());


    // Generate all possible Cartesian exponents corresponding to this shell, according to its angular momentum
    std::vector<CartesianExponents> all_cartesian_exponents;
    all_cartesian_exponents.reserve(this->numberOfBasisFunctions());

    // Partition l into maximally 3 integers and then make all permutations of these partitions
    auto unique_partitions = uniquePartitions<3>(this->l);

    for (const auto& partition : unique_partitions) {
        CartesianExponents exponents (partition);
        for (const auto& permutation : exponents.allPermutations()) {
            all_cartesian_exponents.push_back(permutation);
        }
    }

    std::sort(all_cartesian_exponents.begin(), all_cartesian_exponents.end());


    // Create the explicit basis functions corresponding to the previously constructed exponents
    // The basis functions are linear combinations of CartesianGTOs
    for (const auto& cartesian_exponents : all_cartesian_exponents) {

        // Construct the 'functions' of the linear combination: CartesianGTO
        std::vector<CartesianGTO> gtos;
        gtos.reserve(this->contractionLength());

        for (size_t i = 0; i < this->contractionLength(); i++) {
            double alpha = this->exponents[i];
            gtos.emplace_back(alpha, cartesian_exponents, this->atom.position);
        }

        bfs.emplace_back(LinearCombination<double, CartesianGTO>(this->coefficients, gtos));
    }

    return bfs;
}


/**
 *  @return the length of the contraction in the shell, i.e. the number of primitives contracted in this shell
 */
size_t Shell::contractionLength() const {
    return this->coefficients.size();
}


}  // namespace GQCP
