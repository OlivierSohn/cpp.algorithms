
namespace imajuscule
{
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
  template<typename A, typename B>
  void applySetup(SplitConvolution<A,ScaleConvolution<B>> &c,
                  typename SplitConvolution<A,ScaleConvolution<B>>::SetupParam const & p) {
    applySetup(c.getA(), p.aParams);
    applySetup(c.getB(), p.bParams);
    c.setSplit(pow2(c.getB().getEarlyDroppedConvolutions()) - 1);
  }
  
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
    std::array<SP, 4> ps {
      params,
      (n_scales >= 2)?params:zero,
      (n_scales >= 3)?params:zero,
      (n_scales >= 4)?params:zero
    };
    assert(n_scales <= 4);
    typename ZeroLatencyScaledFineGrainedPartitionnedConvolution<T,FFTTag>::SetupParam p
    {
      {},
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
  
  
  
  
}
