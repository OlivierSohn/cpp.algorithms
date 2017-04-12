
namespace imajuscule {
    namespace measurereverb {
        std::ofstream myfile;
        
        template<typename ALGO>
        void test(int exponent) {
            
            using FPT = typename ALGO::FPT;
            
            auto N = pow2(exponent);
            constexpr auto n_samples = 10000000;
            auto in_seconds_mono = n_samples / static_cast<float>(44100);
            
            std::cout << "Measuring " << N << " "; COUT_TYPE(ALGO); std::cout << std::endl;
            
            ALGO algo;
            
            algo.setCoefficients(std::vector<FPT>(N));
            {
                clock_t t = clock();
                for(int i=0; i<n_samples; ++i) {
                    algo.step(0);
                    algo.get();
                }
                float time = (float)(clock() - t)/CLOCKS_PER_SEC;
                std::cout << "for " << in_seconds_mono << "s : " << time << std::endl;
                auto ratio = time / in_seconds_mono;
                std::cout << "ratio : " << ratio << std::endl;
                myfile << ratio << ", ";
            }
        }
        
        template<typename ALGO>
        void test_nocache(int exponent) {
            
            using FPT = typename ALGO::FPT;
            
            auto N = pow2(exponent);
            auto n_samples = N; // so that there is one computation only, to avoid cache effects
            auto in_seconds_mono = n_samples / static_cast<float>(44100);
            
            std::cout << "Measuring " << N << " "; COUT_TYPE(ALGO); std::cout << std::endl;
            
            ALGO algo;
            
            algo.setCoefficients(std::vector<FPT>(N));
            {
                for(int i=0; i<n_samples-1; ++i) { // measure only last iteration
                    algo.step(0);
                    algo.get();
                }
                clock_t t = clock();
                algo.step(0);
                algo.get();
                float time = (float)(clock() - t)/CLOCKS_PER_SEC;
                std::cout << "for " << in_seconds_mono << "s : " << time << std::endl;
                auto ratio = time / in_seconds_mono;
                std::cout << "ratio : " << ratio << std::endl;
                myfile << ratio << ", ";
                
                // compute the min chunk size that the os asks at each computation so that we have no problem:
                //auto n_iterations = N / chunk_sz;
                // such that 1/n_iterations == ratio;
                // because in that case, one step doesn't exceed the corresponding real time of audio.
                auto min_chunck_sz = N * ratio;
                std::cout << "min chunck size for " << N << " : " << min_chunck_sz << std::endl;
            }
        }
        
    }
}


TEST(Benchmark, reverb_algo) {
    using namespace imajuscule::measurereverb;
    myfile.open ("/Users/Olivier/reverb_msr.csv");
    
    using namespace imajuscule::measurereverb;
    
    std::array<int, 9> exponents { {10, 11, 12, 13, 14, 15, 16, 17, 18} };
    
    for(auto e : exponents) {
        
        //test<FIRFilter<float>>(e);
        //test<FIRFilter<double>>(e);
        test<FFTConvolution<float>>(e);
        test<FFTConvolution<double>>(e);
        
        myfile << std::endl;
    }
    
    myfile.close();
}

/*
 * this test shows that if we make the fft computation more asynchronous,
 * it should be ok: for long ffts (260000 items ie 6s reverb in mono), only 2% of the time is used to compute it! 
 * the problem is that these 2% are not well distributed : if we assume that at each callback call, the os asks for
 * 512 frames, we will have 507 iterations where no computation occurs, and one iteration where the computation happends
 */
TEST(Benchmark, reverb_algo_nocache) {
    using namespace imajuscule::measurereverb;
    myfile.open ("/Users/Olivier/reverb_msr_nocache.csv");
    
    using namespace imajuscule::measurereverb;
    
    std::array<int, 9> exponents { {10, 11, 12, 13, 14, 15, 16, 17, 18} };
    
    for(auto e : exponents) {
        
        //test_nocache<FIRFilter<float>>(e);
        //test_nocache<FIRFilter<double>>(e);
        test_nocache<FFTConvolution<float>>(e);
        test_nocache<FFTConvolution<double>>(e);
        
        myfile << std::endl;
    }
    
    myfile.close();
}
