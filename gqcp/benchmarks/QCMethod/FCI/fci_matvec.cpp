/**
 *  A benchmark executable that times the performance of one FCI matrix-vector product in a full spin-resolved ONV basis with 10 orbitals and 2 to 5 electron pairs.
 */

#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/HamiltonianBuilder/FCI.hpp"

#include <benchmark/benchmark.h>


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 2; i < 6; ++i) {  // need int instead of size_t
        b->Args({10, i});          // spatial orbitals, electron pairs
    }
}


static void matvec(benchmark::State& state) {

    const size_t K = state.range(0);    // number of spatial orbitals
    const size_t N_P = state.range(1);  // number of electron pairs


    // Set up a second-quantized Hamiltonian and a full spin-resolved ONV basis
    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Random(K);
    GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};


    GQCP::FCI fci {onv_basis};
    const auto diagonal = fci.calculateDiagonal(sq_hamiltonian);
    const auto x = onv_basis.randomExpansion();

    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        const auto matvec = fci.matrixVectorProduct(sq_hamiltonian, x, diagonal);

        benchmark::DoNotOptimize(matvec);  // make sure that the variable is not optimized away by compiler
    }

    state.counters["Spatial orbitals"] = K;
    state.counters["Electron pairs"] = N_P;
    state.counters["Dimension"] = onv_basis.dimension();
}


BENCHMARK(matvec)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK_MAIN();
