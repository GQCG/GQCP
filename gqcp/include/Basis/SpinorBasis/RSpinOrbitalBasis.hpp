// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#pragma once


#include "Basis/Integrals/IntegralCalculator.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/SpinorBasis/CurrentDensityMatrixElement.hpp"
#include "Basis/SpinorBasis/SimpleSpinOrbitalBasis.hpp"
#include "Basis/SpinorBasis/Spinor.hpp"
#include "Basis/Transformations/JacobiRotation.hpp"
#include "Basis/Transformations/RTransformation.hpp"
#include "Domain/MullikenDomain/RMullikenDomain.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"
#include "Operator/FirstQuantized/AngularMomentumOperator.hpp"
#include "Operator/FirstQuantized/CoulombRepulsionOperator.hpp"
#include "Operator/FirstQuantized/CurrentDensityOperator.hpp"
#include "Operator/FirstQuantized/DiamagneticOperator.hpp"
#include "Operator/FirstQuantized/ElectronicDensityOperator.hpp"
#include "Operator/FirstQuantized/ElectronicDipoleOperator.hpp"
#include "Operator/FirstQuantized/ElectronicQuadrupoleOperator.hpp"
#include "Operator/FirstQuantized/FQMolecularHamiltonian.hpp"
#include "Operator/FirstQuantized/FQMolecularMagneticHamiltonian.hpp"
#include "Operator/FirstQuantized/KineticOperator.hpp"
#include "Operator/FirstQuantized/LinearMomentumOperator.hpp"
#include "Operator/FirstQuantized/NuclearAttractionOperator.hpp"
#include "Operator/FirstQuantized/OrbitalZeemanOperator.hpp"
#include "Operator/FirstQuantized/OverlapOperator.hpp"
#include "Operator/SecondQuantized/EvaluableRSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/RSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/RSQTwoElectronOperator.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Utilities/complex.hpp"
#include "Utilities/type_traits.hpp"


namespace GQCP {


/**
 *  A restricted spin-orbital basis, i.e. a spin-orbital basis where the alpha- and beta-spinors are equal.
 *
 *  @tparam _ExpansionScalar        The scalar type used to represent an expansion coefficient of the spin-orbitals in the underlying scalar orbitals: real or complex.
 *  @tparam _Shell                  The type of shell the underlying scalar basis contains.
 */
template <typename _ExpansionScalar, typename _Shell>
class RSpinOrbitalBasis:
    public SimpleSpinOrbitalBasis<_ExpansionScalar, _Shell, RSpinOrbitalBasis<_ExpansionScalar, _Shell>> {
public:
    // The scalar type used to represent an expansion coefficient of the spinors in the underlying scalar orbitals: real or complex.
    using ExpansionScalar = _ExpansionScalar;

    // The type of shell the underlying scalar bases contain.
    using Shell = _Shell;

    // The type of the base spinor basis.
    using BaseSpinorBasis = SimpleSpinorBasis<_ExpansionScalar, RSpinOrbitalBasis<_ExpansionScalar, _Shell>>;

    // The type of transformation that is naturally related to an `RSpinOrbitalBasis`.
    using Transformation = RTransformation<ExpansionScalar>;

    // The type that is used for representing the primitive for a basis function of this spin-orbital basis' underlying AO basis.
    using Primitive = typename Shell::Primitive;

    // The type that is used for representing the underlying basis functions of this spin-orbital basis.
    using BasisFunction = typename Shell::BasisFunction;

    // The type that is used to represent a spatial orbital for this spin-orbital basis.
    using SpatialOrbital = EvaluableLinearCombination<ExpansionScalar, BasisFunction>;

    // The type that represents a density distribution for this spin-orbital basis.
    using DensityDistribution = FunctionProduct<SpatialOrbital>;

    // The type of the derivative of a primitive. The derivative of a Cartesian GTO is a linear combination of Cartesian GTOs.
    using PrimitiveDerivative = EvaluableLinearCombination<double, Primitive>;

    // The type of the derivative of a basis function.
    using BasisFunctionDerivative = EvaluableLinearCombination<double, PrimitiveDerivative>;

    // The type of the derivative of a spatial orbital.
    using SpatialOrbitalDerivative = EvaluableLinearCombination<ExpansionScalar, BasisFunctionDerivative>;

    // The type that represents a current density distribution for this spin-orbital basis.
    using CurrentDensityDistribution = EvaluableLinearCombination<complex, FunctionProduct<SpatialOrbital, SpatialOrbitalDerivative>>;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SimpleSpinOrbitalBasis`'s constructors.
    using SimpleSpinOrbitalBasis<_ExpansionScalar, _Shell, RSpinOrbitalBasis<_ExpansionScalar, _Shell>>::SimpleSpinOrbitalBasis;


    /*
     *  MARK: General information
     */

    /**
     *  @return The number of different spatial orbitals that are used in this restricted spin-orbital basis.
     */
    size_t numberOfSpatialOrbitals() const { return this->scalar_basis.numberOfBasisFunctions(); }

    /**
     *  @return The number of spinors that are described by this restricted spin-orbital basis.
     */
    size_t numberOfSpinors() const { return 2 * this->numberOfSpatialOrbitals(); /* alpha and beta spinors are equal*/ }

    /**
     *  @return The number of spin-orbitals that are described by this restricted spin-orbital basis.
     */
    size_t numberOfSpinOrbitals() const { return this->numberOfSpinors(); }


    /*
     *  MARK: Quantization of first-quantized operators (GTOShell)
     */

    /**
     *  Quantize a first-quantized one-electron operator in this restricted spin-orbital basis.
     *
     *  @param fq_one_op                            The first-quantized one-electron operator.
     *
     *  @tparam FQOneElectronOperator               The type of the first-quantized one-electron operator.
     *
     *  @return The second-quantized operator corresponding to the given first-quantized operator, i.e. expressed in/projected onto this spin-orbital basis.
     */
    template <typename FQOneElectronOperator, typename Z = Shell>
    auto quantize(const FQOneElectronOperator& fq_one_op) const -> enable_if_t<std::is_same<Z, GTOShell>::value, RSQOneElectronOperator<product_t<typename FQOneElectronOperator::Scalar, ExpansionScalar>, typename FQOneElectronOperator::Vectorizer>> {

        using ResultScalar = product_t<typename FQOneElectronOperator::Scalar, ExpansionScalar>;
        using ResultOperator = RSQOneElectronOperator<ResultScalar, typename FQOneElectronOperator::Vectorizer>;

        const auto one_op_par = IntegralCalculator::calculateLibintIntegrals(fq_one_op, this->scalarBasis());  // in AO/scalar basis

        ResultOperator op {one_op_par};   // op for 'operator'
        op.transform(this->expansion());  // now in the spatial/spin-orbital basis
        return op;
    }


    /**
     *  Quantize the Coulomb operator in this restricted spin-orbital basis.
     *
     *  @param fq_op                The first-quantized Coulomb operator.
     *
     *  @return The second-quantized operator corresponding to the Coulomb operator.
     */
    template <typename Z = Shell>
    auto quantize(const CoulombRepulsionOperator& fq_op) const -> enable_if_t<std::is_same<Z, GTOShell>::value, RSQTwoElectronOperator<product_t<CoulombRepulsionOperator::Scalar, ExpansionScalar>, CoulombRepulsionOperator::Vectorizer>> {

        using ResultScalar = product_t<CoulombRepulsionOperator::Scalar, ExpansionScalar>;
        using ResultOperator = RSQTwoElectronOperator<ResultScalar, CoulombRepulsionOperator::Vectorizer>;

        const auto g_par = IntegralCalculator::calculateLibintIntegrals(fq_op, this->scalarBasis());  // In AO/scalar basis.

        auto g = SquareRankFourTensor<ResultScalar>::Zero(g_par.dimension(0));
        for (size_t i = 0; i < g_par.dimension(0); i++) {
            for (size_t j = 0; j < g_par.dimension(1); j++) {
                for (size_t k = 0; k < g_par.dimension(2); k++) {
                    for (size_t l = 0; l < g_par.dimension(3); l++) {
                        g(i, j, k, l) = g_par(i, j, k, l);
                    }
                }
            }
        }

        ResultOperator op {g};            // 'op' for 'operator'.
        op.transform(this->expansion());  // Now in spatial/spin-orbital basis.
        return op;
    }


    /*
     *  MARK: Quantization of first-quantized operators (GTOShell)
     */

    /**
     *  Quantize the angular momentum operator in this restricted spin-orbital basis.
     *
     *  @param fq_one_op                            The first-quantized angular momentum operator.
     *
     *  @return The second-quantized angular momentum operator, i.e. expressed in/projected onto this spin-orbital basis.
     */
    template <typename Z = Shell>
    auto quantize(const AngularMomentumOperator& fq_one_op) const -> enable_if_t<std::is_same<Z, GTOShell>::value, RSQOneElectronOperator<product_t<AngularMomentumOperator::Scalar, ExpansionScalar>, typename AngularMomentumOperator::Vectorizer>> {

        using ResultScalar = product_t<AngularMomentumOperator::Scalar, ExpansionScalar>;
        using ResultOperator = RSQOneElectronOperator<ResultScalar, AngularMomentumOperator::Vectorizer>;

        auto primitive_engine = GQCP::IntegralEngine::InHouse<GTOShell>(fq_one_op);
        const auto one_op_par = IntegralCalculator::calculate(primitive_engine, this->scalarBasis().shellSet(), this->scalarBasis().shellSet());  // In AO/scalar basis.

        std::array<SquareMatrix<ResultScalar>, 3> one_op_par_square;
        for (size_t i = 0; i < 3; i++) {
            one_op_par_square[i] = SquareMatrix<ResultScalar>(one_op_par[i]);
        }

        ResultOperator op {one_op_par_square};  // 'op' for 'operator'.
        op.transform(this->expansion());        // Now in the spin-orbital basis.
        return op;
    }


    /**
     *  Quantize the linear momentum operator in this restricted spin-orbital basis.
     *
     *  @param fq_one_op                            The first-quantized linear momentum operator.
     *
     *  @return The second-quantized linear momentum operator, i.e. expressed in/projected onto this spin-orbital basis.
     */
    template <typename Z = Shell>
    auto quantize(const LinearMomentumOperator& fq_one_op) const -> enable_if_t<std::is_same<Z, GTOShell>::value, RSQOneElectronOperator<product_t<LinearMomentumOperator::Scalar, ExpansionScalar>, typename LinearMomentumOperator::Vectorizer>> {

        using ResultScalar = product_t<LinearMomentumOperator::Scalar, ExpansionScalar>;
        using ResultOperator = RSQOneElectronOperator<ResultScalar, LinearMomentumOperator::Vectorizer>;

        auto primitive_engine = GQCP::IntegralEngine::InHouse<GTOShell>(fq_one_op);
        const auto one_op_par = IntegralCalculator::calculate(primitive_engine, this->scalarBasis().shellSet(), this->scalarBasis().shellSet());  // In AO/scalar basis.

        std::array<SquareMatrix<ResultScalar>, 3> one_op_par_square;
        for (size_t i = 0; i < 3; i++) {
            one_op_par_square[i] = SquareMatrix<ResultScalar>(one_op_par[i]);
        }

        ResultOperator op {one_op_par_square};  // 'op' for 'operator'.
        op.transform(this->expansion());        // Now in the spatial/spin-orbital basis.
        return op;
    }


    /**
     *  Quantize the (one-electron) electronic density operator.
     *
     *  @param fq_density_op                    The first-quantized density operator.
     *
     *  @return The second-quantized density operator.
     */
    ScalarEvaluableRSQOneElectronOperator<DensityDistribution> quantize(const ElectronicDensityOperator& fq_density_op) const {

        // There aren't any 'integrals' to be calculated for the density operator: we can just multiply every pair of spatial orbitals.
        const auto phi = this->spatialOrbitals();
        const auto K = this->numberOfSpatialOrbitals();

        SquareMatrix<DensityDistribution> rho_par {K};
        for (size_t p = 0; p < K; p++) {
            for (size_t q = 0; q < K; q++) {
                rho_par(p, q) = phi[p] * phi[q];
            }
        }

        return ScalarEvaluableRSQOneElectronOperator<DensityDistribution> {rho_par};
    }


    /**
     *  Quantize the (one-electron) current density operator.
     *
     *  @param fq_current_density_op            The first-quantized current density operator.
     *
     *  @return The second-quantized current density operator.
     */
    VectorEvaluableRSQOneElectronOperator<CurrentDensityMatrixElement<ExpansionScalar, CartesianGTO>> quantize(const CurrentDensityOperator& fq_current_density_op) const {

        using namespace GQCP::literals;

        // There aren't any 'integrals' to be calculated for the current density operator: we can just multiply every pair of (spatial orbital, spatial orbital derivative).
        const auto K = this->numberOfSpatialOrbitals();
        const auto phi = this->spatialOrbitals();
        const auto dphi = this->spatialOrbitalGradients();

        std::vector<SquareMatrix<CurrentDensityMatrixElement<ExpansionScalar, CartesianGTO>>> j_par_vector {};
        for (size_t i = 0; i < 3; i++) {

            SquareMatrix<CurrentDensityMatrixElement<ExpansionScalar, CartesianGTO>> j_i {K};
            for (size_t p = 0; p < K; p++) {
                for (size_t q = 0; q < K; q++) {
                    j_i(p, q) = CurrentDensityMatrixElement<ExpansionScalar, CartesianGTO> {phi[p], phi[q], dphi[p](i), dphi[q](i)};
                }
            }

            j_par_vector.push_back(j_i);
        }

        return VectorEvaluableRSQOneElectronOperator<CurrentDensityMatrixElement<ExpansionScalar, CartesianGTO>> {j_par_vector};
    }


    /*
     *  MARK: Quantization of first-quantized operators (LondonGTOShell)
     */

    /**
     *  Quantize a spin-independent one-electron operator in this general spinor basis. Spin-independent one-electron operators are those whose two-component matrix operator form contains the same scalar operator in the top-left and bottom-right corner.
     *
     *  @param fq_one_op            A spin-independent first-quantized operator.
     *
     *  @return The second-quantized representation of the given operator.
     */
    template <typename FQOneElectronOperator, typename Z = Shell>
    auto quantize(const FQOneElectronOperator& fq_one_op) const -> enable_if_t<std::is_same<Z, LondonGTOShell>::value, RSQOneElectronOperator<product_t<typename FQOneElectronOperator::Scalar, ExpansionScalar>, typename FQOneElectronOperator::Vectorizer>> {

        using ResultScalar = product_t<typename FQOneElectronOperator::Scalar, ExpansionScalar>;
        using ResultOperator = RSQOneElectronOperator<ResultScalar, typename FQOneElectronOperator::Vectorizer>;
        using Vectorizer = typename FQOneElectronOperator::Vectorizer;

        const auto N = FQOneElectronOperator::NumberOfComponents;
        const auto& vectorizer = FQOneElectronOperator::vectorizer;


        // The strategy for calculating the matrix representation of the one-electron operator in this spin-orbital basis is to:
        //  1. Express the operator in the underlying scalar bases and;
        //  2. Afterwards transform them using the current coefficient matrix.
        const auto K = this->numberOfSpatialOrbitals();


        // 1. Express the operator in the underlying scalar basis.
        auto engine = GQCP::IntegralEngine::InHouse<GQCP::LondonGTOShell>(fq_one_op);
        const auto F = GQCP::IntegralCalculator::calculate(engine, this->scalarBasis().shellSet(), this->scalarBasis().shellSet());


        // For each of the components of the operator, place the scalar basis representations into the spinor basis representation.
        std::array<SquareMatrix<ResultScalar>, N> fs;
        for (size_t i = 0; i < N; i++) {
            fs[i] = SquareMatrix<ResultScalar>(F[i]);
        }


        // 2. Transform using the current coefficient matrix.
        StorageArray<SquareMatrix<ResultScalar>, Vectorizer> array {fs, vectorizer};
        ResultOperator op {array};  // 'op' for 'operator'.
        op.transform(this->expansion());
        return op;
    }


    /**
     *  Quantize the orbital Zeeman operator in this restricted spin-orbital basis.
     *
     *  @param op               The (first-quantized) orbital Zeeman operator.
     *
     *  @return The orbital Zeeman operator expressed in this restricted spin-orbital basis.
     */
    template <typename Z = Shell>
    auto quantize(const OrbitalZeemanOperator& op) const -> enable_if_t<std::is_same<Z, LondonGTOShell>::value, RSQOneElectronOperator<product_t<OrbitalZeemanOperator::Scalar, ExpansionScalar>, OrbitalZeemanOperator::Vectorizer>> {

        // Return the orbital Zeeman operator as a contraction beween the magnetic field and the angular momentum operator.
        const auto L = this->quantize(op.angularMomentum());
        const auto& B = op.magneticField().strength();
        return 0.5 * L.dot(B);
    }


    /**
     *  Quantize the quadrupole operator in this restricted spin-orbital basis.
     *
     *  @param op               The (first-quantized) electronic quadrupole operator.
     *
     *  @return The quadrupole operator expressed in this restricted spin-orbital basis.
     */
    template <typename Z = Shell>
    auto quantize(const ElectronicQuadrupoleOperator& op) const -> enable_if_t<std::is_same<Z, LondonGTOShell>::value, RSQOneElectronOperator<product_t<ElectronicQuadrupoleOperator::Scalar, ExpansionScalar>, ElectronicQuadrupoleOperator::Vectorizer>> {
        return this->quantize(op);
    }


    /**
     *  Quantize the diamagnetic operator in this restricted spin-orbital basis.
     *
     *  @param op               The (first-quantized) diamagnetic operator.
     *
     *  @return The diamagnetic operator expressed in this restricted spin-orbital basis.
     */
    template <typename Z = Shell>
    auto quantize(const DiamagneticOperator& op) const -> enable_if_t<std::is_same<Z, LondonGTOShell>::value, RSQOneElectronOperator<product_t<DiamagneticOperator::Scalar, ExpansionScalar>, DiamagneticOperator::Vectorizer>> {

        using ResultScalar = product_t<DiamagneticOperator::Scalar, ExpansionScalar>;
        using ResultOperator = RSQOneElectronOperator<ResultScalar, DiamagneticOperator::Vectorizer>;


        // Return the diamagnetic operator as a contraction beween the magnetic field and the electronic quadrupole operator.
        const auto Q = this->quantize(op.electronicQuadrupole()).allParameters();
        const auto& B = op.magneticField().strength();

        // Prepare some variables.
        const auto& B_x = B(CartesianDirection::x);
        const auto& B_y = B(CartesianDirection::y);
        const auto& B_z = B(CartesianDirection::z);

        const auto& Q_xx = Q[DyadicCartesianDirection::xx];
        const auto& Q_xy = Q[DyadicCartesianDirection::xy];
        const auto& Q_xz = Q[DyadicCartesianDirection::xz];
        const auto& Q_yy = Q[DyadicCartesianDirection::yy];
        const auto& Q_yz = Q[DyadicCartesianDirection::yz];
        const auto& Q_zz = Q[DyadicCartesianDirection::zz];


        SquareMatrix<ResultScalar> D_par = 0.125 * ((std::pow(B_y, 2) + std::pow(B_z, 2)) * Q_xx +
                                                    (std::pow(B_x, 2) + std::pow(B_z, 2)) * Q_yy +
                                                    (std::pow(B_x, 2) + std::pow(B_y, 2)) * Q_zz -
                                                    2 * B_x * B_y * Q_xy -
                                                    2 * B_x * B_z * Q_xz -
                                                    2 * B_y * B_z * Q_yz);

        return ResultOperator {D_par};
    }


    /**
     *  Quantize the Coulomb operator in this restricted spin-orbital basis.
     *
     *  @param fq_op                The first-quantized Coulomb operator.
     *
     *  @return The second-quantized operator corresponding to the Coulomb operator.
     */
    template <typename Z = Shell>
    auto quantize(const CoulombRepulsionOperator& fq_op) const -> enable_if_t<std::is_same<Z, LondonGTOShell>::value, RSQTwoElectronOperator<product_t<CoulombRepulsionOperator::Scalar, ExpansionScalar>, CoulombRepulsionOperator::Vectorizer>> {

        using ResultScalar = product_t<CoulombRepulsionOperator::Scalar, ExpansionScalar>;
        using ResultOperator = RSQTwoElectronOperator<ResultScalar, CoulombRepulsionOperator::Vectorizer>;


        auto coulomb_engine = GQCP::IntegralEngine::InHouse<GQCP::LondonGTOShell>(CoulombRepulsionOperator());
        const auto g_par = GQCP::IntegralCalculator::calculate(coulomb_engine, this->scalarBasis().shellSet(), this->scalarBasis().shellSet())[0];  // In AO basis.

        auto g = SquareRankFourTensor<ResultScalar>::Zero(g_par.dimension(0));

        for (size_t i = 0; i < g_par.dimension(0); i++) {
            for (size_t j = 0; j < g_par.dimension(1); j++) {
                for (size_t k = 0; k < g_par.dimension(2); k++) {
                    for (size_t l = 0; l < g_par.dimension(3); l++) {
                        g(i, j, k, l) = g_par(i, j, k, l);
                    }
                }
            }
        }

        ResultOperator op {g};            // 'op' for 'operator'.
        op.transform(this->expansion());  // Now in spatial/spin-orbital basis.
        return op;
    }


    /**
     *  Quantize the molecular magnetic Hamiltonian.
     *
     *  @param fq_hamiltonian           The molecular magnetic Hamiltonian.
     *
     *  @return The second-quantized molecular magnetic Hamiltonian.
     */
    template <typename Z = Shell>
    enable_if_t<std::is_same<Z, LondonGTOShell>::value, RSQHamiltonian<ExpansionScalar>> quantize(const FQMolecularMagneticHamiltonian& fq_hamiltonian) const {

        const auto T = this->quantize(fq_hamiltonian.kinetic());
        const auto OZ = this->quantize(fq_hamiltonian.orbitalZeeman());
        const auto D = this->quantize(fq_hamiltonian.diamagnetic());

        const auto V = this->quantize(fq_hamiltonian.nuclearAttraction());

        const auto g = this->quantize(fq_hamiltonian.coulombRepulsion());

        return RSQHamiltonian<ExpansionScalar> {{T, OZ, D, V}, {g}};
    }


    /*
     *  MARK: Quantization of first-quantized operators
     */

    /**
     *  Quantize the molecular Hamiltonian.
     *
     *  @param fq_hamiltonian           The molecular Hamiltonian.
     *
     *  @return The second-quantized molecular Hamiltonian.
     */
    RSQHamiltonian<ExpansionScalar> quantize(const FQMolecularHamiltonian& fq_hamiltonian) const {

        const auto T = this->quantize(fq_hamiltonian.kinetic());
        const auto V = this->quantize(fq_hamiltonian.nuclearAttraction());

        const auto g = this->quantize(fq_hamiltonian.coulombRepulsion());

        return RSQHamiltonian<ExpansionScalar> {T + V, g};
    }


    /*
     *  MARK: Orbitals
     */

    /**
     *  @return The set of spatial orbitals that is associated to this spin-orbital basis.
     */
    std::vector<SpatialOrbital> spatialOrbitals() const {

        // The spatial orbitals are a linear combination of the basis functions, where every column of the coefficient matrix describes one expansion of a spatial orbital in terms of the basis functions.
        const auto basis_functions = this->scalar_basis.basisFunctions();
        const auto& C = this->C;


        // For all spatial orbitals, proceed to calculate the contraction between the associated coefficient matrix column and the basis functions.
        std::vector<SpatialOrbital> spatial_orbitals;
        spatial_orbitals.reserve(this->numberOfSpatialOrbitals());
        for (size_t p = 0; p < this->numberOfSpatialOrbitals(); p++) {

            // Calculate the spatial orbitals as a contraction between a column of the coefficient matrix and the basis functions.
            SpatialOrbital spatial_orbital {};
            for (size_t mu = 0; mu < basis_functions.size(); mu++) {
                const auto coefficient = this->expansion().matrix().col(p)(mu);
                const auto& function = basis_functions[mu];
                spatial_orbital.append(coefficient, function);
            }

            spatial_orbitals.push_back(spatial_orbital);
        }

        return spatial_orbitals;
    }


    /**
     *  @return The gradients of each of the spatial orbitals that is associated to this spin-orbital basis.
     */
    std::vector<Vector<SpatialOrbitalDerivative, 3>> spatialOrbitalGradients() const {

        const auto K = this->numberOfSpatialOrbitals();
        const auto spatial_orbitals = this->spatialOrbitals();

        std::vector<GQCP::Vector<SpatialOrbitalDerivative, 3>> spatial_orbital_gradients {K};
        for (size_t m = 0; m < 3; m++) {
            for (size_t p = 0; p < K; p++) {
                const auto& spatial_orbital = spatial_orbitals[p];

                // A spatial orbital is a linear combination of basis functions (which are contracted GTOs).
                const auto& expansion_coefficients = spatial_orbital.coefficients();
                const auto& basis_functions = spatial_orbital.functions();

                std::vector<Vector<BasisFunctionDerivative, 3>> basis_function_gradients {K};
                for (size_t mu = 0; mu < K; mu++) {
                    const auto& expansion_coefficient = expansion_coefficients[mu];
                    const auto& basis_function = basis_functions[mu];
                    const auto contraction_length = basis_function.length();

                    // A basis function (a.k.a. contracted GTO) is a contraction (i.e. a linear combination) of Cartesian GTOs.
                    const auto& contraction_coefficients = basis_function.coefficients();
                    const auto& primitives = basis_function.functions();

                    for (size_t d = 0; d < contraction_length; d++) {
                        const auto& contraction_coefficient = contraction_coefficients[d];
                        const auto primitive_gradient = primitives[d].calculatePositionGradient();
                        basis_function_gradients[mu](m).append(contraction_coefficient, primitive_gradient(m));
                    }

                    spatial_orbital_gradients[p](m).append(expansion_coefficient, basis_function_gradients[mu](m));
                }
            }
        }

        return spatial_orbital_gradients;
    }


    /**
     *  @return The set of spin-orbitals that is associated to this spin-orbital basis.
     */
    std::vector<Spinor<ExpansionScalar, BasisFunction>> spinOrbitals() const {

        // The spin-orbitals for a restricted spin-orbital basis can be easily constructed from the spatial orbitals, by assigning a zero component once for the beta component of the spin-orbital and once for the alpha component of the spin-orbital.
        const auto spatial_orbitals = this->spatialOrbitals();

        std::vector<Spinor<ExpansionScalar, BasisFunction>> spin_orbitals;
        spin_orbitals.reserve(this->numberOfSpinors());
        for (const auto& spatial_orbital : spatial_orbitals) {

            // Add the alpha- and beta-spin-orbitals accordingly.
            const Spinor<ExpansionScalar, BasisFunction> alpha_spin_orbital {spatial_orbital, 0};  // the '0' int literal can be converted to a zero EvaluableLinearCombination
            spin_orbitals.push_back(alpha_spin_orbital);

            const Spinor<ExpansionScalar, BasisFunction> beta_spin_orbital {0, spatial_orbital};
            spin_orbitals.push_back(beta_spin_orbital);
        }

        return spin_orbitals;
    }


    /*
     *  MARK: Mulliken domain
     */

    /**
     *  Partition this set of restricted spin-orbitals according to the Mulliken partitioning scheme.
     *
     *  @param selector             A function that returns true for basis functions that should be included the Mulliken domain.
     *
     *  @return A `RMullikenDomain` for the AOs selected by the supplied selector function.
     */
    RMullikenDomain<ExpansionScalar> mullikenDomain(const std::function<bool(const BasisFunction&)>& selector) const {

        const auto ao_indices = this->scalarBasis().basisFunctionIndices(selector);
        return RMullikenDomain<ExpansionScalar> {ao_indices, ao_indices.size()};
    }


    /**
     *  Partition this set of restricted spin-orbitals according to the Mulliken partitioning scheme.
     *
     *  @param selector             A function that returns true for shells that should be included the Mulliken domain.
     *
     *  @return A `RMullikenDomain` for the AOs selected by the supplied selector function.
     */
    RMullikenDomain<ExpansionScalar> mullikenDomain(const std::function<bool(const Shell&)>& selector) const {

        const auto ao_indices = this->scalarBasis().basisFunctionIndices(selector);
        return RMullikenDomain<ExpansionScalar> {ao_indices, ao_indices.size()};
    }
};


/*
 *  MARK: Convenience aliases
 */
template <typename ExpansionScalar, typename Shell>
using RSpinorBasis = RSpinOrbitalBasis<ExpansionScalar, Shell>;


/*
 *  MARK: SpinorBasisTraits
 */

/**
 *  A type that provides compile-time information on spinor bases that is otherwise not accessible through a public class alias.
 */
template <typename _ExpansionScalar, typename _Shell>
struct SpinorBasisTraits<RSpinOrbitalBasis<_ExpansionScalar, _Shell>> {
    // The scalar type used to represent an expansion coefficient of the spinors in the underlying scalar orbitals: real or complex.
    using ExpansionScalar = _ExpansionScalar;

    // The type of transformation that is naturally related to an `RSpinOrbitalBasis`.
    using Transformation = RTransformation<ExpansionScalar>;

    // The second-quantized representation of the overlap operator related to an RSpinOrbitalBasis.
    using SQOverlapOperator = ScalarRSQOneElectronOperator<ExpansionScalar>;
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename _ExpansionScalar, typename _Shell>
struct BasisTransformableTraits<RSpinOrbitalBasis<_ExpansionScalar, _Shell>> {

    // The type of transformation that is naturally related to an `RSpinOrbitalBasis`.
    using Transformation = RTransformation<_ExpansionScalar>;
};


/*
 *  MARK: JacobiRotatableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 */
template <typename _ExpansionScalar, typename _Shell>
struct JacobiRotatableTraits<RSpinOrbitalBasis<_ExpansionScalar, _Shell>> {

    // The type of Jacobi rotation that is naturally related to an `RSpinOrbitalBasis`.
    using JacobiRotationType = JacobiRotation;
};


}  // namespace GQCP
