/**
 *  A benchmark executable for the construction of the FCI Hamiltonian, which is constructed for full spin-resolved SpinUnresolvedONV bases for 8 orbitals and 2 to 6 electron pairs.
 */

#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/HamiltonianBuilder/FCI.hpp"

#include <benchmark/benchmark.h>


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 2; i < 6; ++i) {  // need int instead of size_t
        b->Args({8, i});           // number of spatial orbitals, number of electron pairs
    }
}


static void constructHamiltonian(benchmark::State& state) {

    const size_t K = state.range(0);    // number of spatial orbitals
    const size_t N_P = state.range(1);  // number of electron pairs


    // Prepare the second-quantized Hamiltonian and set up a full spin-resolved ONV basis
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Random(K);
    const GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};
    const GQCP::FCI fci {onv_basis};


    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        const auto H = fci.constructHamiltonian(sq_hamiltonian);

        benchmark::DoNotOptimize(H);  // make sure that the variable is not optimized away by compiler
    }

    state.counters["Spatial orbitals"] = K;
    state.counters["Electron pairs"] = N_P;
    state.counters["Dimension"] = onv_basis.dimension();
}


BENCHMARK(constructHamiltonian)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK_MAIN();
