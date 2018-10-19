#ifndef OrbitalOptimizationOptions_hpp
#define OrbitalOptimizationOptions_hpp


namespace GQCP {


    struct OrbitalOptimizationOptions {
        double convergence_threshold = 1.0e-08;
        size_t maximum_number_of_iterations = 128;
    };


}  // namespace GQCP


#endif /* OrbitalOptimizationOptions_hpp */
