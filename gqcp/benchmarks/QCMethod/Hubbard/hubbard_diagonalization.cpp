/**
 *  A benchmark executable that checks the performance of finding Hubbard eigenvectors using a dense solver. The number of electron pairs is kept at 2, while the number of sites varies from 4 to 9.
 */

#include "Mathematical/Optimization/Eigenproblem/EigenproblemSolver.hpp"
#include "Operator/SecondQuantized/ModelHamiltonian/HubbardHamiltonian.hpp"
#include "QCMethod/CI/CI.hpp"
#include "QCMethod/CI/CIEnvironment.hpp"

#include <benchmark/benchmark.h>


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 4; i < 10; ++i) {  // need int instead of size_t
        b->Args({i, 2});            // sites, electron pairs
    }
}


static void diagonalizeHubbardMatrix(benchmark::State& state) {

    // Prepare a random Hubbard model Hamiltonian and an appropriate ONV basis.
    const size_t K = state.range(0);    // number of sites
    const size_t N_P = state.range(1);  // number of electron pairs
    const auto H = GQCP::HoppingMatrix<double>::Random(K);
    const GQCP::HubbardHamiltonian<double> hubbard_hamiltonian {H};

    const GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};


    // Code inside this loop is measured repeatedly.
    // Solve the Hubbard method, which is a CI method of the Hubbard model Hamiltonian in a the full spin-resolved ONV basis.
    for (auto _ : state) {

        // In the QCMethod/QCModel framework, we'll have to specify a dense solver and associated environment.
        auto environment = GQCP::CIEnvironment::Dense(hubbard_hamiltonian, onv_basis);
        auto solver = GQCP::EigenproblemSolver::Dense();

        const auto electronic_energy = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, environment).groundStateEnergy();

        benchmark::DoNotOptimize(electronic_energy);  // make sure that the variable is not optimized away by compiler
    }

    state.counters["Sites"] = K;
    state.counters["Electron pairs"] = N_P;
    state.counters["Dimension"] = onv_basis.dimension();
}


BENCHMARK(diagonalizeHubbardMatrix)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK_MAIN();
