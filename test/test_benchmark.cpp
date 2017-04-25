// notion of benchmark is obsolete and costly to maintain...
/*
namespace imajuscule {
    namespace benchmark_partconv {
        void dump(std::ofstream & f, ConvolutionBenchmark const & b, bool constraint) {
            using namespace std;
            
            if(constraint) {
                f << "with";
            }
            else {
                f << "without";
            }
            f << " multichannel spread" << endl;
            
            for(int i = 0; i < ConvolutionBenchmark::count_frames_exp; ++i) {
                for(int j = 0; j < ConvolutionBenchmark::count_impulse_exp; ++j) {
                    auto min_ = b.get(constraint, i, j);
                    f << min_.first;
                    //f << " (" << min_.second << ")";
                    f << ", ";
                }
                f << endl;
            }
            f << endl;
        }
        
        void dump(std::ofstream & f, ConvolutionBenchmark const & b) {
            dump(f, b, false);
            if(b.hasConstraint()) {
                dump(f, b, true);
            }
        }
    }
}

// reason for DISABLED:
// doesn't test anything + slow
TEST(BenchmarkPartitionnedConvolution, DISABLED_ToFile) {
    using namespace imajuscule;
    using namespace benchmark_partconv;

    // with 30 iterations we get consistent results
    constexpr auto n_iterations = 30;
    ConvolutionBenchmarks b(n_iterations);
    {
        std::ofstream f;
        f.open ("/Users/Olivier/bench_partconv_mono.csv");
        dump(f, b.getBenchmark(1));
        f.close();
    }
    {
        std::ofstream f;
        f.open ("/Users/Olivier/bench_partconv_stereo.csv");
        dump(f, b.getBenchmark(2));
        f.close();
    }
}

// reason for DISABLED:
// doesn't test anything + slow
TEST(BenchmarkPartitionnedConvolution, DISABLED_withOptimizer) {
    using namespace imajuscule;
    constexpr auto n_iterations = 1;
    ConvolutionOptimizer o(n_iterations);
    
    int n_channels = 2;
    int n_audio_frames_per_cb = 192;
    bool can_use_spread = true;
    bool use_spread;
    
    int ir_length_fact = 1000;
    for(int i=1; i<300; ++i) {
        int ir_length = ir_length_fact * i;
        auto lg2sz = o.getLg2PartSize(n_channels,
                                      n_audio_frames_per_cb,
                                      ir_length,
                                      can_use_spread,
                                      use_spread);

        std::cout << ir_length << " : " << pow2(lg2sz) << (use_spread ? " with" : " without") << " spread" << std::endl;
        
    }
}
*/
