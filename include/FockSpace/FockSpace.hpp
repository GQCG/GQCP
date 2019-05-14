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
#ifndef GQCP_FOCKSPACE_HPP
#define GQCP_FOCKSPACE_HPP


#include "FockSpace/BaseFockSpace.hpp"
#include "FockSpace/FockPermutator.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"

#include <Eigen/Sparse>


namespace GQCP {


/**
 *  The full Fock space for a number of orbitals and number of electrons
 *
 *  The ONVs and addresses are linked with a hashing function calculated with an addressing scheme. The implementation of the addressing scheme is from Molecular Electronic-Structure Theory (August 2000) by Trygve Helgaker, Poul Jorgensen, and Jeppe Olsen
 *
 */
class FockSpace: public BaseFockSpace, public FockPermutator<FockSpace> {
private:
    Matrixu vertex_weights;  // vertex_weights of the addressing scheme


    template<class T>
    class conta {
        T cont;

        conta(size_t dim, size_t memory) {
            cont = T::Zero(dim, dim);
        }

        void add(size_t i, size_t j, double value) {
            cont(i, j) = value;
        }

        T get_matrix() {
            return cont;
        }

        conta() {

        }

    };

    template<>
    class conta<Eigen::SparseMatrix<double>> {
        Eigen::SparseMatrix<double> cont;
        std::vector<Eigen::Triplet<double>> triplet_vector;
        size_t dim;
        size_t mem;
        conta(size_t dim, size_t memory) {
            cont =  Eigen::SparseMatrix<double>(dim, dim);
            mem = memory;
            triplet_vector.reserve(mem);
        }

        Eigen::SparseMatrix<double> get_matrix() {
            this->cont = Eigen::SparseMatrix<double>(dim, dim);
            this->cont.setFromTriplets(this->triplet_vector.begin(),this->triplet_vector.end());
            return cont;
        }

        void add(size_t i, size_t j, double value) {
            this->triplet_vector.emplace_back(i,j, value);
        }
    };

public:
    // CONSTRUCTORS
    /**
     *  @param K        the number of orbitals
     *  @param N        the number of electrons
     */
    FockSpace(size_t K, size_t N);


    // DESTRUCTORS
    ~FockSpace() override = default;


    // GETTERS
    size_t get_vertex_weights(size_t p, size_t m) const { return this->vertex_weights[p][m]; }
    const Matrixu& get_vertex_weights() const { return this->vertex_weights; }
    FockSpaceType get_type() const override { return FockSpaceType::FockSpace; }


    // STATIC PUBLIC METHODS
    /**
     *  @param K        the number of orbitals
     *  @param N        the number of electrons
     *
     *  @return the dimension of the Fock space
     */
    static size_t calculateDimension(size_t K, size_t N);


    // PUBLIC OVERRIDEN METHODS
    /**
     *  @param representation       a representation of an ONV
     *
     *  @return the next bitstring permutation in the Fock space
     *
     *      Examples:
     *          011 -> 101
     *          101 -> 110
     */
    size_t ulongNextPermutation(size_t representation) const override;

    /**
     *  @param representation      a representation of an ONV
     *
     *  @return the address (i.e. the ordering number) of the given ONV
     */
    size_t getAddress(size_t representation) const override;

    /**
      *  Calculate unsigned representation for a given address
      *
      *  @param address                 the address of the representation is calculated
      *
      *  @return unsigned representation of the address
      */
    size_t calculateRepresentation(size_t address) const override;

    /**
     *  @param onv       the ONV
     *
     *  @return the amount of ONVs (with a larger address) this ONV would couple with given a one electron operator
     */
    size_t countOneElectronCouplings(const ONV& onv) const override;

    /**
     *  @param onv       the ONV
     *
     *  @return the amount of ONVs (with a larger address) this ONV would couple with given a two electron operator
     */
    size_t countTwoElectronCouplings(const ONV& onv) const override;

    /**
     *  @return the amount non-zero (non-diagonal) couplings of a one electron coupling scheme in the Fock space
     */
    size_t countTotalOneElectronCouplings() const override;

    /**
     *  @return the amount non-zero (non-diagonal) couplings of a two electron coupling scheme in the Fock space
     */
    size_t countTotalTwoElectronCouplings() const override;


    // PUBLIC METHODS
    /**
     *  If we have
     *      FockSpace fock_space;
     *
     *  This makes sure that we can call
     *      fock_space.getAddress(onv);
     *  instead of the syntax
     *      fock_space.FockPermutator<FockSpace>::getAddress(onv);
     */
    using FockPermutator<FockSpace>::getAddress;

    /**
     *  Find the next unoccupied orbital in a given ONV,
     *  update the electron count, orbital index,
     *  and update the address by calculating a shift
     *  resulting from a difference between the initial vertex weights for the encountered occupied orbitals
     *  and the corrected vertex weights accounting for previously annihilated electrons
     *
     *  @tparam T        the amount of previously annihilated electrons
     *
     *  @param address   the address which is updated
     *  @param onv       the ONV for which we search the next unnocupied orbital
     *  @param q         the orbital index
     *  @param e         the electron count
     */
    template<int T>
    void shiftUntilNextUnoccupiedOrbital(const ONV& onv, size_t& address, size_t& q, size_t& e) const {

        // Test whether the current orbital index is occupied
        while (e < this->N && q == onv.get_occupation_index(e)) {

            // Take the difference of vertex weights for the encountered electron weights to that of a vertex weight path with "a" fewer electrons
            // +1 is added to the electron index, because of how the addressing scheme is arrayed.
            address += this->get_vertex_weights(q, e + 1 - T) - this->get_vertex_weights(q, e + 1);

            // move to the next electron and orbital
            e++;
            q++;
        }
    }

    /**
     *  Find the next unoccupied orbital in a given ONV,
     *  update the electron count, orbital index, sign,
     *  and update the address by calculating a shift
     *  resulting from a difference between the initial vertex weights for the encountered occupied orbitals
     *  and the corrected vertex weights accounting for previously annihilated electrons
     *
     *  @tparam T        the amount of previously annihilated electrons
     *
     *  @param address   the address which is updated
     *  @param onv       the ONV for which we search the next unnocupied orbital
     *  @param q         the orbital index
     *  @param e         the electron count
     *  @param sign      the sign which is flipped for each iteration
     */
    template<int T>
    void shiftUntilNextUnoccupiedOrbital(const ONV& onv, size_t& address, size_t& q, size_t& e, int& sign) const {

        // Test whether the current orbital index is occupied
        while (e < this->N && q == onv.get_occupation_index(e)) {

            // Take the difference of vertex weights for the encountered electron weights to that of a vertex weight path with "a" fewer electrons
            // +1 is added to the electron index, because of how the addressing scheme is arrayed.
            address += this->get_vertex_weights(q, e + 1 - T) - this->get_vertex_weights(q, e + 1);

            // move to the next electron and orbital
            e++;
            q++;
            sign *= -1;
        }
    }

    /**
     *  Find the previous unoccupied orbital in a given ONV,
     *  update the electron count, orbital index, sign,
     *  and update the address by calculating a shift
     *  resulting from a difference between the initial vertex weights for the encountered occupied orbitals
     *  and the corrected vertex weights accounting for newly created electrons
     *
     *  @tparam T        the amount of newly created electrons
     *
     *  @param address   the address which is updated
     *  @param onv       the ONV for which we search the next unoccupied orbital
     *  @param q         the orbital index
     *  @param e         the electron count
     *  @param sign      the sign which is flipped for each iteration
     */
    template<int T>
    void shiftUntilPreviousUnoccupiedOrbital(const ONV &onv, size_t &address, size_t &q, size_t &e, int &sign) const {

        // Test whether the current orbital index is occupied
        while (e != -1 && q == onv.get_occupation_index(e)) {

            int shift = static_cast<int>(this->get_vertex_weights(q, e + 1 + T)) - static_cast<int>(this->get_vertex_weights(q, e + 1));
            address += shift;

            e--;
            q--;
            sign *= -1;
        }
    }

   /**
    *  eval
    */
   template<class T>
   T evaluateOperator(const OneElectronOperator<double>& oe, bool diagonal_values = false) {

       size_t K = this->get_K();
       size_t N = this->get_N();
       size_t dim = this->get_dimension();
       size_t memory =  this->countTotalTwoElectronCouplings();

       if (diagonal_values) {
           memory += dim;
       }

       conta<T> evaluation (dim, memory);

       ONV onv = this->makeONV(0);  // onv with address 0
       for (size_t I = 0; I < dim; I++) {  // I loops over all the addresses of the onv
           for (size_t e1 = 0; e1 < N; e1++) {  // e1 (electron 1) loops over the (number of) electrons
               size_t p = onv.get_occupation_index(e1);  // retrieve the index of a given electron
               // remove the weight from the initial address I, because we annihilate
               size_t address = I - this->get_vertex_weights(p, e1 + 1);

               if(diagonal_values){
                   evaluation.add(I, I, oe(p, p));
               }

               // The e2 iteration counts the amount of encountered electrons for the creation operator
               // We only consider greater addresses than the initial one (because of symmetry)
               // Hence we only count electron after the annihilated electron (e1)
               size_t e2 = e1 + 1;
               size_t q = p + 1;

               int sign_e2 = 1;
               // perform a shift
               this->shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign_e2);

               while (q < K) {
                   size_t J = address + this->get_vertex_weights(q, e2);
                   double value = sign_e2*oe(p, q);
                   evaluation.add(I, J, value);
                   evaluation.add(J, I, value);

                   q++; // go to the next orbital

                   // perform a shift
                   this->shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign_e2);
               }  //  (creation)
           } // e1 loop (annihilation)

           // Prevent last permutation
           if (I < dim - 1) {
               this->setNextONV(onv);
           }
       }

       return evaluation.get_matrix();

   }

    /**
     *  eval
     */
    template<class T>
    T evaluateOperator(const TwoElectronOperator<double>& te, bool diagonal_values = false) {

        size_t K = this->get_K();
        size_t N = this->get_N();
        size_t dim = this->get_dimension();
        size_t memory =  this->countTotalTwoElectronCouplings();

        if (diagonal_values) {
            memory += dim;
        }

        conta<T> evaluation (dim, memory);

        // The algorithm is written for the notation in which the two-electron operators are decomposed in a two- and one-electron partition.
        auto k = te.effectiveOneElectronPartition();

        ONV onv = this->makeONV(0);  // onv with address 0
        for (size_t I = 0; I < dim; I++) {  // I loops over all addresses in the Fock space
            if (I > 0) {
                this->setNextONV(onv);
            }
            int sign1 = -1;  // start with -1 because we flip at the start of the annihilation (so we start at 1, followed by:  -1, 1, ...)
            for (size_t e1 = 0; e1 < N; e1++) {  // A1 (annihilation 1)

                sign1 *= -1;
                size_t p = onv.get_occupation_index(e1);  // retrieve the index of a given electron
                size_t address = I - this->get_vertex_weights(p, e1 + 1);

                // Strictly diagonal values
                if (diagonal_values) {
                    evaluation.add(I, I, k(p,p));
                    for (size_t q = 0; q < K; q++) {  // q loops over SOs
                        if (onv.isOccupied(q)) {  // q is in Ia
                            evaluation(I, I) += 0.5 * te(p, p, q, q);
                        } else {  // q is not in I_alpha
                            evaluation(I, I) += 0.5 * te(p, q, q, p);
                        }
                    }
                }

                size_t address1 = address;
                size_t e2 = e1;
                size_t q = p;

                /**
                 *  A1 > C1 (annihlation 1 > creation 1)
                 */
                int sign2 = sign1;
                q--;
                e2--;
                this->shiftUntilPreviousUnoccupiedOrbital<1>(onv, address1, q, e2, sign2);
                while (q != -1) {

                    size_t address2 = address1 + this->get_vertex_weights(q, e2 + 2);

                    /**
                     *  C2 > A2
                     */
                    int sign3 = sign1;
                    for (size_t e3 = e1 + 1; e3 < N; e3++) {
                        sign3 *= -1;  // initial sign3 = sign of the annhilation, with one extra electron(from crea) = *-1
                        size_t r = onv.get_occupation_index(e3);
                        size_t address3 = address2 - this->get_vertex_weights(r, e3 + 1);

                        size_t e4 = e3 + 1;
                        size_t s = r + 1;

                        int sign4 = sign3;
                        this->shiftUntilNextUnoccupiedOrbital<1>(onv, address3, s, e4, sign4);

                        while (s < K) {
                            size_t J = address3 + this->get_vertex_weights(s, e4);
                            int signev = sign1 * sign2 * sign3 * sign4;
                            double value = signev * 0.5 * (te(p, q, r, s) +
                                                           te(r, s, p, q) -
                                                           te(p, s, r, q) -
                                                           te(r, q, p, s));


                            evaluation.add(I,J, value);
                            evaluation.add(J,I, value);

                            s++;
                            this->shiftUntilNextUnoccupiedOrbital<1>(onv, address3, s, e4, sign4);
                        }
                    }
                    q--;
                    this->shiftUntilPreviousUnoccupiedOrbital<1>(onv, address1, q, e2, sign2);
                }

                e2 = e1 + 1;
                q = p + 1;
                sign2 = sign1;
                this->shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign2);

                /**
                 *  A1 < C1
                 */
                while (q < K) {

                    address1 = address + this->get_vertex_weights(q, e2);

                    /**
                     *  A2 > C1
                     */
                    int sign3 = sign2;
                    for (size_t e3 = e2; e3 < N; e3++) {
                        sign3 *= -1; // -1 cause we created electron (creation) sign of A is now the that of C *-1
                        size_t r = onv.get_occupation_index(e3);
                        size_t address3 = address1 - this->get_vertex_weights(r, e3 + 1);

                        size_t e4 = e3 + 1;
                        size_t s = r + 1;
                        int sign4 = sign3;
                        this->shiftUntilNextUnoccupiedOrbital<1>(onv, address3, s, e4, sign4);

                        while (s < K) {
                            size_t J = address3 + this->get_vertex_weights(s, e4);
                            int signev = sign1 * sign2 * sign3 * sign4;

                            double value = signev * 0.5 * (te(p, q, r, s) +
                                                           te(r, s, p, q) -
                                                           te(r, q, p, s) -
                                                           te(p, s, r, q));

                            evaluation.add(I,J, value);
                            evaluation.add(J,I, value);

                            s++;  // go to the next orbital
                            this->shiftUntilNextUnoccupiedOrbital<1>(onv, address3, s, e4, sign4);

                        }  // (creation)

                    }

                    size_t r = q;
                    sign3 = sign2;
                    size_t address1c = address1;

                    /**
                     *  A2 < C1, (A2 > A1)
                     */
                    for (size_t e3 = e2 - 1; e3 > e1; e3--) {
                        sign3 *= -1;
                        size_t e4 = e2;
                        address1c += this->get_vertex_weights(r, e3) -
                                     this->get_vertex_weights(r, e3 + 1);
                        r = onv.get_occupation_index(e3);
                        size_t address2 = address1c - this->get_vertex_weights(r, e3);
                        int sign4 = sign2;
                        size_t s = q + 1;
                        this->shiftUntilNextUnoccupiedOrbital<1>(onv, address2, s, e4, sign4);
                        while (s < K) {

                            size_t J = address2 + this->get_vertex_weights(s, e4);

                            int signev = sign1 * sign2 * sign3 * sign4;
                            double value = signev * 0.5 * (te(p, q, r, s) +
                                                           te(r, s, p, q) -
                                                           te(r, q, p, s) -
                                                           te(p, s, r, q));

                            evaluation.add(I,J, value);
                            evaluation.add(J,I, value);

                            s++;
                            this->shiftUntilNextUnoccupiedOrbital<1>(onv, address2, s, e4, sign4);

                        }
                    }

                    /**
                     *  A2 = C1
                     */
                    int signev = sign2 * sign1;

                    double value_I =  k(p, q);  // cover the one electron calculations

                    for (size_t s = 0; s < K; s++) {
                        if(!onv.isOccupied(s)){
                            value_I += 0.5 * (te(p, s, s, q));
                        } else {

                            value_I  += 0.5 *  (te(s, s, p, q) - te(s, q, p, s) + te(p, q, s, s));
                        }
                    }

                    value_I *= signev;

                    q++;

                    evaluation.add(I,address1, value_I);
                    evaluation.add(address1,I, value_I);

                    this->shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign2);
                }
            }
        }

        return evaluation.get_matrix();

    }

    /**
     *  eval
     */
    template<class T>
    T evaluateOperator(const HamiltonianParameters<double>& ham_par) {

        return evaluateOperator<T>(ham_par.get_h()) + evaluateOperator<T>(ham_par.get_g());

    }
};


}  // namespace GQCP


#endif  // GQCP_FOCKSPACE_HPP
