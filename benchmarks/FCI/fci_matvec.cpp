/**
 *  A benchmark executable for the FCI matvec
 */

#include <benchmark/benchmark.h>

#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "HamiltonianBuilder/FCI.hpp"


static void matvec(benchmark::State& state) {

    // Prepare the Hamiltonian
    size_t K = state.range(0);
    size_t N = state.range(1);
    GQCP::ProductFockSpace fock_space (K, N, N);
    GQCP::FCI fci (fock_space);

    GQCP::SQHamiltonian<double> sq_hamiltonian = GQCP::SQHamiltonian<double>::Random(K);
    GQCP::VectorX<double> diagonal = fci.calculateDiagonal(sq_hamiltonian);
    GQCP::VectorX<double> x = fock_space.randomExpansion();

    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        GQCP::VectorX<double> matvec = fci.matrixVectorProduct(sq_hamiltonian, x, diagonal);

        benchmark::DoNotOptimize(matvec);  // make sure the variable is not optimized away by compiler
    }

    state.counters["Orbitals"] = K;
    state.counters["Electron pairs"] = N;
    state.counters["Dimension"] = fock_space.get_dimension();
}


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 2; i < 6; ++i) {  // need int instead of size_t
        b->Args({10, i});  // orbitals, electron pairs
    }
}


// Perform the benchmarks
BENCHMARK(matvec)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK_MAIN();
