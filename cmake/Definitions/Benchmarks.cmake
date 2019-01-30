# List all benchmarks


# Find the benchmarks folder
set(PROJECT_BENCHMARKS_FOLDER ${CMAKE_SOURCE_DIR}/benchmarks)


# Find the source files for the benchmarks
set(PROJECT_BENCHMARK_SOURCE_FILES
        ${PROJECT_BENCHMARKS_FOLDER}/DOCI/doci_case.cpp
        ${PROJECT_BENCHMARKS_FOLDER}/DOCI/doci_matrix.cpp
        ${PROJECT_BENCHMARKS_FOLDER}/DOCI/doci_matvec.cpp

        ${PROJECT_BENCHMARKS_FOLDER}/FCI/fci_hchain.cpp
        ${PROJECT_BENCHMARKS_FOLDER}/FCI/fci_matrix.cpp
        ${PROJECT_BENCHMARKS_FOLDER}/FCI/fci_matvec.cpp

        ${PROJECT_BENCHMARKS_FOLDER}/Hubbard/hubbard_diagonalization.cpp
        ${PROJECT_BENCHMARKS_FOLDER}/Hubbard/hubbard_matrix.cpp
        ${PROJECT_BENCHMARKS_FOLDER}/Hubbard/hubbard_matvec.cpp
    )
