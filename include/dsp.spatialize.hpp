/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    namespace audio {
        
        /*
         * Naïve unoptimized version, to test functionnality.
         *
         * Sources are spatialized, i.e. each source location provides one impulse response per ear
         */
        template<int N_EARS, typename Convolution>
        struct Spatializer {

            static constexpr auto nEars = N_EARS;
            
            using T = typename Convolution::FPT;
            using FPT = T;
            
            FPT getEpsilon() const { return convs[0][0].getEpsilon(); }
            
            void set_partition_size(int sz) {
                size_partition = sz;
            }
            
            void clear() {
                forEachEar([](auto & convolutions) {
                    convolutions.clear();
                });
            }
            
            bool empty() const { return convs.empty() || convs[0].empty() || convs[0][0].empty(); }
            
            bool isValid() const {
                return convs[0][0].isValid();
            }
            
            template<typename T>
            void applySetup(T const & setup) {
                forEachEar([&setup](auto & convolutions) {
                    for(auto & c : convolutions) {
                        c.applySetup(setup);
                    }
                });
            }
            
            int countSources() const { return convs[0].size(); }
            
            int getLatency() const {
                return convs[0][0].getLatency();
            }
            
            void addSourceLocation(std::array<a64::vector<T>, nEars> vcoeffs) {
                forEachEar(vcoeffs, [this](auto & convolutions, auto & coeffs) {
                    convolutions.emplace_back();
                    
                    auto & c = convolutions.back();
                    assert(size_partition >= 0);
                    c.set_partition_size(size_partition);
                    c.setCoefficients(std::move(coeffs));
                    
                });
            }

            void setMaxMultiplicationGroupLength() {
                forEachEar([](auto & convolutions) {
                    for(auto & c : convolutions) {
                        c.setMultiplicationGroupLength(c.getHighestValidMultiplicationsGroupSize());
                    }
                });
            }

            void setMultiplicationGroupLength(int i) {
                forEachEar([i](auto & convolutions) {
                    for(auto & c : convolutions) {
                        c.setMultiplicationGroupLength(i);
                    }
                });
            }

            void step(T const * vals) {
                forEachEar([vals](auto & convolutions) {
                    int i = 0;
                    for(auto & c : convolutions) {
                        c.step(vals[i]);
                        ++i;
                    }
                });
            }
            
            void get(T * samples) {
                forEachEar(samples, [](auto & convolutions, auto & s) {
                    s = 0;
                    for(auto & c : convolutions) {
                        s += c.get();
                    }
                });
            }
            
            void dephaseComputations(int phase) {
                int n = 0;
                forEachEar([phase_=phase, &n](auto & convolutions) {
                    for(auto & c : convolutions) {
                        auto const phase = n * phase_;
                        for(int j=0; j<phase; ++j) {
                            c.step(0);
                        }
                        ++n;
                    }
                });
            }
            
        private:
            std::array<std::vector<Convolution>, nEars> convs;
            int size_partition = -1;
            
            template<typename Input, typename F>
            void forEachEar(Input & i, F f) {
                //assert(i.size() == nEars);
                int j=0;
                for(auto & c : convs) {
                    f(c, i[j]);
                    ++j;
                }
            }

            template<typename F>
            void forEachEar(F f) {
                for(auto & c : convs) {
                    f(c);
                }
            }
        };
    }

}
