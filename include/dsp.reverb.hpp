
namespace imajuscule
{

//////////////////////////////////////////

  template<typename C>
  void applySetup(C&c, typename C::SetupParam const & p) {
    c.applySetup(p);
  }
  
  template<typename C>
  void applySetup(Delayed<C>&c, typename Delayed<C>::SetupParam const & p) {
    c.setTheDelay(p.delay);
    applySetup(c.getInner(),p.innerParams);
  }
  template<LatencySemantic L, typename C>
  void applySetup(SubSampled<L, C>&c, typename C::SetupParam const & p) {
    applySetup(c.getInner(), p);
  }
  
  template<typename A, typename B>
  void applySetup(SplitConvolution<A,B> &c, typename SplitConvolution<A,B>::SetupParam const & p) {
    applySetup(c.getA(), p.aParams);
    applySetup(c.getB(), p.bParams);
    if(!c.getB().isValid()) { // for case where lower resolution tail is not used
      c.setSplit(noSplit);
    }
    else {
      c.setSplit( B::nCoefficientsFadeIn + c.getB().getLatency() - c.getA().getLatency() );
    }
  }


////////////////////////////////////////

  template<typename SP, typename T, typename FFTTag>
  void prepare(SP const & params,
               FinegrainedPartitionnedFFTConvolution<T,FFTTag> & rev,
               int const n_scales,
               int const scale_sz ) {    
    assert(n_scales == 1);
    applySetup(rev, params);
  }
  
  template<typename SP, typename T, typename FFTTag>
  void prepare(SP const & params,
               ZeroLatencyScaledFineGrainedPartitionnedConvolution<T,FFTTag> & rev,
               int & n_scales,
               int const scale_sz ) {
    int delay = 0;
    if(n_scales > 1) {
      delay = SameSizeScales::getDelays(scale_sz, params.partition_size);
      assert(delay > 0);
    }

    // set the delays
    
    // we disable the unused scales by setting the partition size to 0.
    auto zero = SP::makeInactive();
    static_assert(4==nMaxScales);
    std::array<SP, nMaxScales> ps {
      params,
      (n_scales >= 2)?params:zero,
      (n_scales >= 3)?params:zero,
      (n_scales >= 4)?params:zero
    };
    assert(n_scales <= nMaxScales);
    typename ZeroLatencyScaledFineGrainedPartitionnedConvolution<T,FFTTag>::SetupParam p
    {
        {{},{}},
      {
        ps[0],
        {
          (n_scales >= 2)?delay:0,
          {
            ps[1],
            {
              (n_scales >= 3)?delay:0,
              {
                ps[2],
                {
                  (n_scales >= 4)?delay:0,
                  ps[3]
                }
              }
            }
          }
        }
      }
    };
    applySetup(rev, p);
  }
  
  template<typename SP, typename T, typename FFTTag>
  void prepare(SP const & params,
               ZeroLatencyScaledAsyncConvolution<T,FFTTag> & rev,
               int & n_scales,
               int const scale_sz ) {
    assert(n_scales <= 1);
    applySetup(rev, params);
  }

  template<typename SP, typename T>
  void prepare(SP const & params,
               FIRFilter<T> & rev,
               int const n_scales,
               int const scale_sz) {
    assert(n_scales <= 1);
    applySetup(rev, params);
  }

template<typename SP, typename T, typename FFTTag>
void prepare(SP const & params,
             OptimizedFIRFilter<T, FFTTag> & rev,
             int const n_scales,
             int const scale_sz) {
    assert(n_scales <= 1);
    applySetup(rev, params);
}

  /////////////////////////////////


template<typename C>
void dephase(int phase,
             int n_scales,
             C & rev) {
    if constexpr (C::has_subsampling) {
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
    else {
        assert(n_scales <= 1);
        for(int j=0; j<phase; ++j) {
            rev.step(0);
        }
    }
}


  
}
