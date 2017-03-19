
namespace imajuscule {
    namespace fft {
        
        // will be called with those parameters:
        //
        // 0     N
        // 1     N
        // ...
        // N/2-1 N
        //
        // 0     N/2
        // 1     N/2
        // ...
        // N/4-1 N/2
        //
        // and i / M is the same as i*2 / M*2
        //
        // we need to compute N/2-1 values
        
        
        template<typename T>
        static complex<T> make_root_of_unity(unsigned int index, unsigned int size) {
            return polar(-2. * M_PI * index / size);
        }
        
        template<typename T>
        using FFTVec = typename std::vector<complex<T>>;
        
        template<typename T>
        FFTVec<T> compute_roots_of_unity(unsigned int N) {
            assert(is_power_of_two(N));
            auto n_roots = N/2;
            FFTVec<T> res;
            res.reserve(n_roots);
            for(unsigned int i=0; i<n_roots; ++i) {
                res.push_back(make_root_of_unity<T>(i,N));
            }
            return res;
        }
        
        // https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
        
        template<typename T>
        struct Algo {
            using ROOTS_OF_UNITY = FFTVec<T>;
            
            Algo(ROOTS_OF_UNITY const & roots_of_unity) :
            roots_of_unity(roots_of_unity)
            {}
            
            template<typename ITER>
            void run(ITER it,
                     ITER result,
                     unsigned int N,
                     unsigned int stride) {
                static_assert(std::is_same<T,typename ITER::value_type::FPT>::value, "");
                assert(N/2 == roots_of_unity.size());
                do_run(it, result, N, stride);
            }
            
        private:
            ROOTS_OF_UNITY const & roots_of_unity;

            template<typename ITER>
            void do_run(ITER it,
                        ITER result,
                        unsigned int N,
                        unsigned int stride) {
                if(N==1) {
                    *result = *it;
                    return;
                }
                N /= 2;
                do_run(it, result, N, 2*stride );
                do_run(it+stride, result + N, N, 2*stride );
                
                for(unsigned int i=0; i<N; ++i) {
                    auto it_result_a = result + i;
                    auto t = *it_result_a;
                    auto it_result_b = it_result_a + N;
                    auto e = roots_of_unity[i * stride];
                    *it_result_b *= e;
                    *it_result_a += *it_result_b;
                    *it_result_b = t - *it_result_b;
                }
            }
            
        };
        
    } // NS fft
} // NS imajuscule

