/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    namespace audio {
        /*
         * Sources are spatialized, i.e. each source location provides one impulse response per ear
         *
         * TODO perform parameter optimization globally: since we know we will do
         * ** several ** convolutions per sample, we can relax optimization constraints
         * somewhat.
         */
        template<int N_EARS, typename Convolution>
        struct Spatializer {

            static constexpr auto nEars = N_EARS;

            using T = typename Convolution::FPT;
            using FPT = T;

          double getEpsilon() const { return epsilonOfNaiveSummation(convs[0]); }

            void clear() {
                forEachEar([](auto & convolutions) {
                    convolutions.clear();
                });
            }

            bool empty() const { return convs.empty() || convs[0].empty() || convs[0][0].empty(); }

            bool isValid() const {
                return convs[0][0].isValid();
            }

            int countSources() const { return convs[0].size(); }

            int getLatency() const {
                return convs[0][0].getLatency();
            }

            template<typename SetupP>
            void addSourceLocation(std::array<a64::vector<T>, nEars> vcoeffs,
                                   std::pair<int, // size_partition
                                             SetupP> const & setup) {
                forEachEar(vcoeffs, [this,&setup](auto & convolutions, auto & coeffs) {
                    convolutions.emplace_back();
                    auto & c = convolutions.back();
                    auto size_partition = setup.first;
                    assert(size_partition >= 0);
                    setPartitionSize(c,size_partition);
                    ::imajuscule::applySetup(c,setup.second);
                    c.setCoefficients(std::move(coeffs));
                });
            }

          // for tests
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
