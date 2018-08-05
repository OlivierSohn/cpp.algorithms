#ifndef IMJ_USE_SLOW_FFT

namespace imajuscule {
    namespace profiling {
        namespace measure_madd {

            template<typename CONTAINER>
            auto divide(CONTAINER const & v1, CONTAINER const & v2) {
                CONTAINER res;
                res = v1;
                for(int i=0; i<res.size(); ++i) {
                    res[i] /= v2[i];
                }
                return res;
            }

            template<typename CONTAINER>
            void accumulate(CONTAINER & to, CONTAINER add) {
                if(to.empty()) {
                    to = std::move(add);
                    return;
                }
                auto sz = to.size();
                assert(sz == add.size());
                for(int i=0; i<sz; ++i) {
                    to[i] += add[i];
                }
            }

            template<typename CONTAINER>
            void min_(CONTAINER & to, CONTAINER o) {
                if(to.empty()) {
                    to = std::move(o);
                    return;
                }
                auto sz = to.size();
                assert(sz == o.size());
                for(int i=0; i<sz; ++i) {
                    to[i] = std::min(to[i], o[i]);
                }
            }

            // the results could be more pronounced if using vectorize techniques here.
            template<typename ITER>
            void simple_multiply_add( ITER res, ITER i1, ITER i2, int N) {
                using T = typename ITER::value_type;
                using SC = accelerate::SplitComplex<T>;

                SC M1 {
                    &*i1,
                    (&*i1) + N/2
                };
                SC M2 {
                    &*i2,
                    (&*i2) + N/2
                };
                SC Accum {
                    &*res,
                    (&*res) + N/2
                };

                accelerate::API<T>::f_zvma(&M1, 1,
                                           &M2, 1,
                                           &Accum, 1,
                                           &Accum, 1, N/2);


//                for(int i=0; i<N; ++i, ++res, ++i1, ++i2) {
  //                  *res += *i1 * *i2;
    //            }
            }

            constexpr auto n_repeats = 5;

            using CONTAINER = a64::vector<float>;

            static inline auto testNonContiguous(int lg2BlockSize) {
                using namespace std;
                using namespace std::chrono;
                auto block_sz = pow2(lg2BlockSize);

                using VCONTAINER = std::vector<CONTAINER>; // non contiguous

                //                auto lg2_n_blocks_max = 19 - lg2BlockSize;
                auto lg2_n_blocks_max = 10;

                vector<float> times;
                times.reserve(lg2_n_blocks_max);

                CONTAINER w;
                w.resize(block_sz);

                VCONTAINER v1, v2;
                v1.reserve(pow2(lg2_n_blocks_max-1));
                v2.resize(pow2(lg2_n_blocks_max-1));

                for(int lg2_n_blocks=0; lg2_n_blocks<lg2_n_blocks_max; ++lg2_n_blocks) {
                    //auto n_blocks = pow2(lg2_n_blocks);
                    auto n_blocks = lg2_n_blocks+1;

                    v1.resize(n_blocks);
                    v2.resize(n_blocks);
                    for(auto & v : v1) {
                        v.resize(block_sz);
                    }
                    for(auto & v : v2) {
                        v.resize(block_sz);
                    }

                    auto wit = w.begin();

                    auto t = avg(measure_n<high_resolution_clock>(n_repeats, []{}, [wit,&v1,&v2]() {
                        auto n_blocks = v1.size();
                        for(int i=0; i<n_blocks; ++i) {
                            simple_multiply_add(wit, v1[i].begin(), v2[i].begin(), v1[i].size());
                        }
                    }));

                    //cout << "  for " << n_blocks <<  " blocks : " << t << endl;
                    times.push_back( t );
                }
                return std::move(times);
            }

            static inline auto testContiguous(int lg2BlockSize) {
                using namespace std;
                using namespace std::chrono;
                auto block_sz = pow2(lg2BlockSize);

                //                auto lg2_n_blocks_max = 19 - lg2BlockSize;
                auto lg2_n_blocks_max = 10;

                vector<float> times;
                times.reserve(lg2_n_blocks_max);

                CONTAINER w;
                w.resize(block_sz);

                CONTAINER v1, v2;
                v1.reserve(pow2(lg2_n_blocks_max-1) * block_sz);
                v2.reserve(pow2(lg2_n_blocks_max-1) * block_sz);

                for(int lg2_n_blocks=0; lg2_n_blocks<lg2_n_blocks_max; ++lg2_n_blocks) {
                    //auto n_blocks = pow2(lg2_n_blocks);
                    auto n_blocks = lg2_n_blocks+1;

                    v1.resize(n_blocks * block_sz);
                    v2.resize(n_blocks * block_sz);
                    auto it1 = v1.begin();
                    auto it2 = v2.begin();
                    auto wit = w.begin();

                    auto t = avg(measure_n<high_resolution_clock>(n_repeats, []{}, [wit,it1_=it1,it2_=it2,n_blocks, block_sz]() {
                        auto it1 = it1_;
                        auto it2 = it2_;
                        for(int i=0; i<n_blocks; ++i,
                            it1 += block_sz,
                            it2 += block_sz)
                        {
                            simple_multiply_add(wit, it1, it2, block_sz);
                        }
                    }));

                    //cout << "  for " << n_blocks <<  " blocks : " << t << endl;
                    times.push_back( t );
                }
                return std::move(times);
            }

            static inline void run_multiplyadd_test() {
                using namespace imajuscule::profiling::measure_madd;
                using namespace imajuscule;
                using namespace std;

                for(int i=0; i<19; ++i) {
                    std::vector<float> contiguous, noncontiguous;
                    constexpr auto n_repeats = 15;
                    for(int r = 0; r < n_repeats; ++r) {
                        //cout << "nc" << endl;
                        accumulate(noncontiguous, testNonContiguous(i));
                        //cout << "c" << endl;
                        accumulate(contiguous, testContiguous(i));
                    }
                    auto ratio = divide(contiguous, noncontiguous);
                    int index = 0;
                    cout << endl << " block size " << pow2(i) << endl;
                    for(auto r : ratio) {
                        cout << index << ": " << r << " (" << sizeof(float) * (index+1) * pow2(i) << ")" << endl;
                        ++index;
                    }
                }
            }

        }
    }
}

#endif
