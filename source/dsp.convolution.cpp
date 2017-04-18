/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    ConvolutionOptimizer * ConvolutionOptimizer::instance = nullptr;

    int get_lg2_optimal_partition_size(GradientDescent & gd,
                                       int n_iterations,
                                       int n_channels,
                                       int n_frames,
                                       int length_impulse,
                                       bool constraint,
                                       float & min_val,
                                       int n_tests) {
        gd.setFunction( [n_frames, length_impulse, constraint, n_tests, n_channels](int lg2_partition_size, float & val){
            if(lg2_partition_size < 0) {
                return ParamState::OutOfRange;
            }
            if(lg2_partition_size > 20) {
                throw std::logic_error("Gradient descent is not working?");
            }
            auto const partition_size = pow2(lg2_partition_size);
            
            if(constraint) {
                if(n_channels * n_frames >= partition_size) {
                    return ParamState::OutOfRange;
                }
            }
            
            struct Test {
                
                Test(size_t partition_size, int length_impulse) {
                    pfftcv.set_partition_size(partition_size);
                    pfftcv.setCoefficients(cacheline_aligned_allocated::vector<float>(length_impulse));
                    for(int i=0; i<partition_size-1; ++i) {
                        pfftcv.step(0.f); // these should do next to nothing...
                        pfftcv.get();
                    }
                }
                
                void run() {
                    pfftcv.step(0.f); // ... this one should do the ffts
                    pfftcv.get();
                }
            private:
                PartitionnedFFTConvolution<float> pfftcv;
            };
            
            // prepare tests
            
            std::vector<Test> tests;
            tests.reserve(n_tests);
            for(int i=0; i<n_tests;++i) {
                tests.emplace_back(partition_size, length_impulse);
            }
            
            // right now everything might be in cache, so before running the test we need to pollute the cache.
            
            {
                std::vector<int> pollution(100000);
                auto i = 0;
                for(auto & p : pollution) {
                    p = i++;
                }
                i = 0;
                for(auto & p : pollution) {
                    i += p;
                }
                std::cout << ((i%2)?"." : ":");
            }
            
            {
                auto begin = std::chrono::high_resolution_clock::now();
                for(auto & t : tests) {
                    t.run();
                }
                auto end = std::chrono::high_resolution_clock::now();
                
                val = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
            }
            
            val /= n_tests; // val represents 'one computation'
            
            auto n_max_computes_per_callback = n_frames / partition_size;
            if(n_frames != n_max_computes_per_callback * partition_size) {
                // in the worst case, we have one more
                ++ n_max_computes_per_callback;
            }
            if(constraint) {
                if(n_max_computes_per_callback != 1) {
                    throw std::logic_error("the constraint ensures that the number of"
                                           " computes per callback is 1/n_channels on average");
                }
                // n_frames is small enough and partition_size is big enough so that
                // there is enough "room" to spread the computes of different channels over different callback calls,
                // provided we "phase" the different partitionned convolutions correctly.
                // Hence we take this advantage into account here:
                val /= n_channels;
            }
            
            val *= n_max_computes_per_callback; // val represents 'worst computation time over one callback'
            
            val /= n_frames; // val represents 'worst computation time over one callback, averaged per frame'
            
            return ParamState::Ok;
        });
        
        auto start_lg2_partition = 5;
        if(constraint) {
            // to ensure that the constraint is met in first try
            start_lg2_partition = 1 + power_of_two_exponent(n_channels * n_frames);
        }
        
        return gd.findMinimum(n_iterations,
                              start_lg2_partition,
                              min_val);
    }
}
