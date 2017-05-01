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
            
            FPT getEpsilon() const { return Epsilon<Convolution>::get(convs[0][0]); }
            
            void set_partition_size(int sz) {
                size_partition = sz;
            }
            
            bool isValid() const {
                return convs[0][0].isValid();
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
                    
                    // may be suboptimal
                    c.setMultiplicationGroupLength(c.getHighestValidMultiplicationsGroupSize());
                });
            }
            /*
            void setMultiplicationGroupLength(int i) {
                forEachEar([i](auto & out) {
                    for(auto & c : convs) {
                        c.setMultiplicationGroupLength(i);
                    }
                });
            }
            */
            void step(T const * vals) {
                forEachEar([vals](auto & convolutions) {
                    int i = 0;
                    for(auto & c : convolutions) {
                        c.step(vals[i]);
                        ++i;
                    }
                });
            }
            
            void get(std::array<T, nEars> & samples) {
                forEachEar(samples, [](auto & convolutions, auto & s) {
                    s = 0;
                    for(auto & c : convolutions) {
                        s += c.get();
                    }
                });
            }
            
        private:
            std::array<std::vector<Convolution>, nEars> convs;
            int size_partition = -1;
            
            template<typename Input, typename F>
            void forEachEar(Input & i, F f) {
                assert(i.size() == nEars);
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

    template<int N_EARS, typename Convolution>
    struct Epsilon< audio::Spatializer<N_EARS, Convolution> > {
        using T = typename Convolution::FPT;
        static T get(audio::Spatializer<N_EARS, Convolution> const & c) {
            return c.getEpsilon();
        }
    };

}
