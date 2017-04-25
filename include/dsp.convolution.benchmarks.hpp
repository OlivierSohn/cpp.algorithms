/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    // obsolete, takes time to maintain...
    /*
    struct ConvolutionBenchmark {
        static constexpr auto count_frames_exp = 8;
        static constexpr auto count_impulse_exp = 10;
        
        static constexpr auto n_frames_base_exp = 5; // 32
        static constexpr auto l_impulse_base_exp = 10; // 1024
        
        using Minimum = std::pair<int, float>;

        ConvolutionBenchmark(int n_channels, int n_iterations) : n_channels(n_channels) {
            initialize_vectors(false, n_iterations);
            if(n_channels > 1) {
                initialize_vectors(true, n_iterations);
            }
        }
        
        Minimum get(bool with_constraint, int n_frames_index, int l_impulse_index) const {
            if(with_constraint && n_channels == 1) {
                throw std::logic_error("mono channel cannot benefit from multichannel optimization");
            }
            auto const & v = getVectors(with_constraint);
            if(n_frames_index >= v.size()) {
                std::cout << n_frames_index << " " << v.size() << std::endl;
                throw std::logic_error("out of range frame index");
            }
            auto const & w = v[n_frames_index];
            if(l_impulse_index >= w.size()) {
                std::cout << l_impulse_index << " " << w.size() << std::endl;
                throw std::logic_error("out of range impulse index");
            }
            return w[l_impulse_index];
        }

        bool hasConstraint() const { return n_channels > 1; }
    private:
        int n_channels;
        
        using Vectors = std::vector<std::vector<Minimum>> ;
        Vectors best_lg2_partition_size, best_lg2_partition_size_with_channel_spread;
        
        Vectors const & getVectors(bool constraint) const {
            return constraint? best_lg2_partition_size_with_channel_spread : best_lg2_partition_size;
        }
        
        Vectors & editVectors(bool constraint) {
            return constraint? best_lg2_partition_size_with_channel_spread : best_lg2_partition_size;
        }

        void initialize_vectors(bool constraint, int n_iterations) {
            Vectors & v = editVectors(constraint);
            v.resize(count_frames_exp);
            auto n_frames_exp = 0;
            for(auto & w : v) {
                w.resize(count_impulse_exp);
                
                auto const n_frames = pow2(n_frames_exp + n_frames_base_exp);
                
                auto l_impulse_exp = 0;
                for(auto & v : w) {
                    auto const l_impulse = pow2(l_impulse_exp + l_impulse_base_exp);
                    constexpr auto n_tests = 1;
                    GradientDescent gd;
                    v.first = get_lg2_optimal_partition_size(gd,
                                                             n_iterations,
                                                             n_channels,
                                                             n_frames,
                                                             l_impulse,
                                                             constraint,
                                                             v.second,
                                                             n_tests);
                    gd.debug(true);
                    if(v.first > 20) {
                        std::cout << v.first;
                        throw std::logic_error("out of range");
                    }
                    std::cout << "min is " << v.first << std::endl;
                    
                    ++l_impulse_exp;
                }
                ++n_frames_exp;
            }
        }
    };
    
    struct ConvolutionBenchmarks {
        ConvolutionBenchmarks(int n_iterations) :
        per_n_channel_stats
        {{
            {1, n_iterations},
            {2, n_iterations}
        }}
        {}
        
        auto get(int n_channels, bool with_constraint, int n_frames_index, int l_impulse_index) {
            return per_n_channel_stats[n_channels - min_n_channels].get(with_constraint, n_frames_index, l_impulse_index);
        }
        
        static constexpr auto min_n_channels = 1;
        static constexpr auto max_n_channels = 2;

        ConvolutionBenchmark const & getBenchmark(int n_channels) {
            return per_n_channel_stats[n_channels - min_n_channels];
        }
    private:
        
        std::array<ConvolutionBenchmark, 2> per_n_channel_stats;
    };
    
    
    struct ConvolutionOptimizer {
        ConvolutionOptimizer(int n_iterations) : benchmarks(n_iterations) {}
        
        int getLg2PartSize(int n_channels,
                           int n_audio_frames_per_callback,
                           int impulse_response_length,
                           bool can_use_multichannel_spread,
                           bool & use_multichannel_spread) {

            use_multichannel_spread = false;
            
            int l_impulse_index = static_cast<int>(power_of_two_exponent(ceil_power_of_two(impulse_response_length))) - ConvolutionBenchmark::l_impulse_base_exp;
            if(l_impulse_index < 0) {
                assert(0);
                l_impulse_index = 0;
            }
            
            std::pair<int, int> n_frames_index;
            n_frames_index.first = static_cast<int>(power_of_two_exponent(floor_power_of_two(n_audio_frames_per_callback))) - ConvolutionBenchmark::n_frames_base_exp;
            if(n_frames_index.first < 0) {
                assert(0);
                n_frames_index.first = 0;
            }
            
            n_frames_index.second = static_cast<int>(power_of_two_exponent(ceil_power_of_two(n_audio_frames_per_callback))) - ConvolutionBenchmark::n_frames_base_exp;
            if(n_frames_index.second < 0) {
                assert(0);
                n_frames_index.second = 0;
            }
            
            assert(n_frames_index.first <= n_frames_index.second);
            
            if(n_channels == 1 || !can_use_multichannel_spread) {
                auto res = benchmarks.get(n_channels, false, n_frames_index.first, l_impulse_index);
                return res.first;
            }
            
            auto normal =
            benchmarks.get(n_channels, false, n_frames_index.first, l_impulse_index);
            auto constrained =
            benchmarks.get(n_channels, true, n_frames_index.first, l_impulse_index);
            
            if(normal.second < constrained.second) {
                return normal.first;
            }
            else {
                use_multichannel_spread = true;
                return constrained.first;
            }
        }
        
    private:
        static ConvolutionOptimizer * instance;
        ConvolutionBenchmarks benchmarks;
    };*/
}
