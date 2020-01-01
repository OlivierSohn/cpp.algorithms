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
            void flushToSilence() {
                forEachEar([](auto & earConvs) {
                    for(auto & c : earConvs) {
                        c->flushToSilence();
                    }
                });
            }

            bool empty() const { return earsConvs.empty() || earsConvs[0].empty() || earsConvs[0][0]->isZero(); }

            bool isValid() const {
                return !empty() && earsConvs[0][0]->isValid();
            }

            int countSources() const {
              return earsConvs.empty() ? 0 : earsConvs[0].size();
            }
            
            int countConvolutions() const {
                return nEars * countSources();
            }
          
          int countScales() {
              if (earsConvs.empty() || earsConvs[0].empty()) {
                  return 0;
              }
              if constexpr(Convolution::has_subsampling) {
                  return imajuscule::countScales(*earsConvs[0][0]);
              }
              else {
                  return 1;
              }
          }

            int getLatency() const {
                return (earsConvs.empty() || earsConvs[0].empty()) ? 0 : earsConvs[0][0]->getLatency();
            }

            template<typename SetupP>
            void addSourceLocation(std::array<a64::vector<T>, nEars> vcoeffs,
                                   SetupP const & setup) {
                forEachEar(vcoeffs, [&setup](auto & earConvs, auto & coeffs) {
                    earConvs.push_back(std::make_unique<Convolution>());
                    auto & c = earConvs.back();
                    c->setup(setup);
                    c->setCoefficients(std::move(coeffs));
                });
            }

          // for tests
            void setMaxMultiplicationGroupLength() {
                forEachEar([](auto & earConvs) {
                    for(auto & c : earConvs) {
                        c->setMultiplicationGroupLength(c->getHighestValidMultiplicationsGroupSize());
                    }
                });
            }

            void setMultiplicationGroupLength(int i) {
                forEachEar([i](auto & earConvs) {
                    for(auto & c : earConvs) {
                        c->setMultiplicationGroupLength(i);
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
                int n = 0;
                for(auto & c : earConvs) {
                  wetSignal += c->step(in[n]);
                  ++n;
                }
                *inout += wet * wetSignal;
                ++inout;
              }
            }
            
            void addWet(T const * const in, T * out, int stride) {
                for(auto & earConvs : earsConvs) {
                  T wetSignal = 0.;
                  int n = 0;
                  for(auto & c : earConvs) {
                    wetSignal += c->step(in[n]);
                    n += stride;
                  }
                  *out += wetSignal;
                  out += stride;
                }
            }
            
            template<typename FPT2>
            void assignWet(FPT2 const * const * const input_buffers,
                           int nInputBuffers,
                           FPT2 ** output_buffers,
                           int nOutputBuffers,
                           int nFramesToCompute)
            {
                Assert(nOutputBuffers == earsConvs.size());
                for(int i_out=0; i_out < nOutputBuffers; ++i_out) {
                    auto & earConvs = earsConvs[i_out];
                    FPT2 * out = output_buffers[i_out];
                    
                    bool assign = true;
                    Assert(nInputBuffers == earConvs.size());
                    for(int i_in = 0; i_in < nInputBuffers; ++i_in) {
                        auto & c = earConvs[i_in];
                        FPT2 const * const in = input_buffers[i_in];
                        if(assign) {
                            c->stepAssignVectorized(in,
                                                   out,
                                                   nFramesToCompute);
                        }
                        else {
                            c->stepAddVectorized(in,
                                                out,
                                                nFramesToCompute);
                        }
                        assign = false;
                    }
                }
            }
            template<typename FPT2>
            bool assignWetVectorized(FPT2 const * const * const input_buffers,
                                     int nInputBuffers,
                                     FPT2 ** output_buffers,
                                     int nOutputBuffers,
                                     int nFramesToCompute,
                                     int vectorLength)
            {
                bool success = true;
                
                Assert(nOutputBuffers == earsConvs.size());
                for(int i_out=0; i_out < nOutputBuffers; ++i_out) {
                    auto & earConvs = earsConvs[i_out];
                    FPT2 * out = output_buffers[i_out];
                    
                    bool assign = true;
                    Assert(nInputBuffers == earConvs.size());
                    for(int i_in = 0; i_in < nInputBuffers; ++i_in) {
                        auto & c = earConvs[i_in];
                        FPT2 const * const in = input_buffers[i_in];
                        for(int i=0; i<nFramesToCompute; i += vectorLength) {
                            if(assign) {
                                c->stepAssignVectorized(in + i,
                                                       out + i,
                                                       std::min(vectorLength, nFramesToCompute-i));
                            }
                            else {
                                c->stepAddVectorized(in + i,
                                                    out + i,
                                                    std::min(vectorLength, nFramesToCompute-i));
                            }
                        }
                        assign = false;
                        if constexpr (Convolution::step_can_error) {
                            success = !c->hasStepErrors() && success;
                        }
                    }
                }
                return success;
            }
            template<typename FPT2>
            void addWetVectorized(FPT2 const * const * const input_buffers,
                                  int nInputBuffers,
                                  FPT2 ** output_buffers,
                                  int nOutputBuffers,
                                  int nFramesToCompute,
                                  int vectorLength)
            {
                Assert(nOutputBuffers == earsConvs.size());
                for(int i_out=0; i_out < nOutputBuffers; ++i_out) {
                    auto & earConvs = earsConvs[i_out];
                    FPT2 * out = output_buffers[i_out];

                    Assert(nInputBuffers == earConvs.size());
                    for(int i_in = 0; i_in < nInputBuffers; ++i_in) {
                        auto & c = earConvs[i_in];
                        FPT2 const * const in = input_buffers[i_in];
                        for(int i=0; i<nFramesToCompute; i += vectorLength) {
                            c->stepAddVectorized(in + i,
                                                out + i,
                                                std::min(vectorLength, nFramesToCompute-i));
                        }
                    }
                }
            }

            template<typename FPT2>
            bool addWetInputZeroVectorized(FPT2 ** output_buffers,
                                           int nOutputBuffers,
                                           int nFramesToCompute,
                                           int vectorLength) {
                bool success = true;

                Assert(nOutputBuffers == earsConvs.size());
                for(int i_out=0; i_out < nOutputBuffers; ++i_out) {
                    auto & earConvs = earsConvs[i_out];
                    FPT2 * out = output_buffers[i_out];
                    for(int i=0; i<nFramesToCompute; i += vectorLength) {
                        for(auto & c : earConvs) {
                            c->stepAddInputZeroVectorized(out+i,
                                                          std::min(vectorLength, nFramesToCompute-i));
                            if constexpr (Convolution::step_can_error) {
                                success = !c->hasStepErrors() && success;
                            }
                        }
                    }
                }
                return success;
            }

            template<typename PS>
            void dephaseComputations(PS spec) {
                int n = 0;
                int const total = countConvolutions();
                forEachEar([spec, total, &n](auto & earConvs) {
                    for(auto & c : earConvs) {
                      dephase(total, n, spec, *c);
                      ++n;
                    }
                });
            }
            
            template<typename F>
            void foreachConvReverb(F f) const {
                for(auto const & earConvs : earsConvs) {
                    for(auto const & c : earConvs) {
                        f(*c);
                    }
                }
            }
        private:
            std::array<std::vector<std::unique_ptr<Convolution>>, nEars> earsConvs;

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
