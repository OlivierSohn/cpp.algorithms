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

          double getEpsilon() const { return epsilonOfNaiveSummation(earsConvs[0]); }

            void clear() {
                forEachEar([](auto & earConvs) {
                    earConvs.clear();
                });
            }

            bool empty() const { return earsConvs.empty() || earsConvs[0].empty() || earsConvs[0][0].empty(); }

            bool isValid() const {
                return !empty() && earsConvs[0][0].isValid();
            }

            int countSources() const {
              return earsConvs.empty() ? 0 : earsConvs[0].size();
            }

            int getLatency() const {
                return (earsConvs.empty() || earsConvs[0].empty()) ? 0 : earsConvs[0][0].getLatency();
            }

            template<typename SetupP>
            void addSourceLocation(std::array<a64::vector<T>, nEars> vcoeffs,
                                   std::pair<int, // size_partition
                                             SetupP> const & setup) {
                forEachEar(vcoeffs, [this,&setup](auto & earConvs, auto & coeffs) {
                    earConvs.emplace_back();
                    auto & c = earConvs.back();
                    auto size_partition = setup.first;
                    assert(size_partition >= 0);
                    setPartitionSize(c,size_partition);
                    ::imajuscule::applySetup(c,setup.second);
                    c.setCoefficients(std::move(coeffs));
                });
            }

          // for tests
            void setMaxMultiplicationGroupLength() {
                forEachEar([](auto & earConvs) {
                    for(auto & c : earConvs) {
                        c.setMultiplicationGroupLength(c.getHighestValidMultiplicationsGroupSize());
                    }
                });
            }

            void setMultiplicationGroupLength(int i) {
                forEachEar([i](auto & earConvs) {
                    for(auto & c : earConvs) {
                        c.setMultiplicationGroupLength(i);
                    }
                });
            }

            void step(T * inout, T const dry, T const wet) {
              Assert(dry == 0 || getLatency() == 0); // else dry and wet signals are out of sync
              
              std::array<T, nEars> in;
              for(int i=0; i<nEars; ++i) {
                in[i] = inout[i];
                inout[i] *= dry;
              }

              for(auto & earConvs : earsConvs) {

                T wetSignal = 0.;
                int i = 0;
                for(auto & c : earConvs) {
                  wetSignal += c.step(in[i]);
                  ++i;
                }
                *inout += wet * wetSignal;
                ++inout;
              }
            }

            void dephaseComputations(int phase) {
                int n = 0;
                forEachEar([phase_=phase, &n](auto & earConvs) {
                    for(auto & c : earConvs) {
                        auto const phase = n * phase_;
                        for(int j=0; j<phase; ++j) {
                            c.step(0);
                        }
                        ++n;
                    }
                });
            }

        private:
            std::array<std::vector<Convolution>, nEars> earsConvs;

            template<typename Input, typename F>
            void forEachEar(Input & i, F f) {
                //assert(i.size() == nEars);
                int j=0;
                for(auto & earConvs : earsConvs) {
                    f(earConvs, i[j]);
                    ++j;
                }
            }

            template<typename F>
            void forEachEar(F f) {
                for(auto & earConvs : earsConvs) {
                    f(earConvs);
                }
            }
        };
    }

}
