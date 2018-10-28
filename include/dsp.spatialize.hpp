/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
  template<typename C>
  int countScales(C & rev) {
    auto & lateHandler = rev.getB();
    {
      auto & inner = lateHandler.getB().getInner().getInner();
      if(inner.isZero()) {
        return 1;
      }
    }
    {
      auto & inner = lateHandler.getB().getInner().getInner().getB().getInner().getInner();
      if(inner.isZero()) {
        return 2;
      }
    }
    {
      auto & inner = lateHandler.getB().getInner().getInner().getB().getInner().getInner().getB().getInner().getInner();
      if(inner.isZero()) {
        return 3;
      }
    }
    static_assert(4==nMaxScales);
    return 4;
  }
  
  template<typename C>
  bool scalesAreValid(int n_scales, C & rev) {
    auto & lateHandler = rev.getB();
    for(int i=1; i<n_scales; ++i) {
      // the top-most will be stepped 'base_phase' times,
      // then each scale after that will be stepped by a quarter grain size.
      // phases are cumulative, so stepping a scale also steps subsequent scales.
      switch(i) {
        case 1:
        {
          auto & inner = lateHandler.getB().getInner().getInner();
          if(!inner.isValid()) {
            return false;
          }
          if(inner.isZero()) {
            return false;
          }
          break;
        }
        case 2:
        {
          auto & inner = lateHandler.getB().getInner().getInner().getB().getInner().getInner();
          if(!inner.isValid()) {
            return false;
          }
          if(inner.isZero()) {
            return false;
          }
          break;
        }
        case 3:
        {
          auto & inner = lateHandler.getB().getInner().getInner().getB().getInner().getInner().getB().getInner().getInner();
          if(!inner.isValid()) {
            return false;
          }
          if(inner.isZero()) {
            return false;
          }
          break;
        }
        default:
          throw std::logic_error("out of bound");
      }
    }
    
    return true;
  }
  
  template<typename C>
  void dephase(int phase, int n_scales, C & rev) {
    auto & lateHandler = rev.getB();
    for(int i=0; i<n_scales; ++i) {
      // the top-most will be stepped 'base_phase' times,
      // then each scale after that will be stepped by a quarter grain size.
      // phases are cumulative, so stepping a scale also steps subsequent scales.
      static_assert(nMaxScales==4);
      switch(i) {
        case 0:
          for(int j=0; j<phase; ++j) {
            lateHandler.step(0);
          }
          break;
        case 1:
        {
          auto & inner = lateHandler.getB().getInner().getInner();
          assert(inner.isValid());
          assert(!inner.isZero());
          int quarter_grain_size = inner.getA().getGranularMinPeriod() / 4;
          for(int j=0; j<quarter_grain_size; ++j) {
            inner.step(0);
          }
          break;
        }
        case 2:
        {
          auto & inner = lateHandler.getB().getInner().getInner().getB().getInner().getInner();
          assert(inner.isValid());
          assert(!inner.isZero());
          int quarter_grain_size = inner.getA().getGranularMinPeriod() / 4;
          for(int j=0; j<quarter_grain_size; ++j) {
            inner.step(0);
          }
          break;
        }
        case 3:
        {
          auto & inner = lateHandler.getB().getInner().getInner().getB().getInner().getInner().getB().getInner().getInner();
          assert(inner.isValid());
          assert(!inner.isZero());
          int quarter_grain_size = inner.getGranularMinPeriod() / 4;
          for(int j=0; j<quarter_grain_size; ++j) {
            inner.step(0);
          }
          break;
        }
        default:
          throw std::logic_error("out of bound");
      }
    }
    
  }

  
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

            bool empty() const { return earsConvs.empty() || earsConvs[0].empty() || earsConvs[0][0].isZero(); }

            bool isValid() const {
                return !empty() && earsConvs[0][0].isValid();
            }

            int countSources() const {
              return earsConvs.empty() ? 0 : earsConvs[0].size();
            }
          
          int countScales() {
            return (earsConvs.empty() || earsConvs[0].empty()) ? 0 : imajuscule::countScales(earsConvs[0][0]);
          }

            int getLatency() const {
                return (earsConvs.empty() || earsConvs[0].empty()) ? 0 : earsConvs[0][0].getLatency();
            }

            template<typename SetupP>
            void addSourceLocation(std::array<a64::vector<T>, nEars> vcoeffs,
                                   SetupP const & setup,
                                   int & n_scales,
                                   int scale_sz) {
                forEachEar(vcoeffs, [this,&n_scales, scale_sz, &setup](auto & earConvs, auto & coeffs) {
                    earConvs.emplace_back();
                    auto & c = earConvs.back();
                    prepare(setup, c, n_scales, scale_sz);
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

            void dephaseComputations(int phase, int n_scales) {
                int n = 0;
                forEachEar([phase, n_scales, &n](auto & earConvs) {
                    for(auto & c : earConvs) {
                      dephase(n * phase, n_scales, c);
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
