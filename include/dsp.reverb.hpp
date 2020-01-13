
namespace imajuscule
{

template<typename C, typename PS>
void dephase(int const total_instances,
             int const index_instance,
             PS const & ps,
             C & rev)
{
    Assert(total_instances);
    {
        float ratio = index_instance / static_cast<float>(total_instances);
        auto periodicities = rev.getComputePeriodicities();
        for(auto &p :periodicities) {
            p = static_cast<int>(p * ratio);
        }
        auto const phase_increments_late_handler = ps.getPhase();
        if(phase_increments_late_handler &&
           phase_increments_late_handler.value() &&
           !periodicities.empty())
        {
            // override late handler phase
            *periodicities.rbegin() = phase_increments_late_handler.value() * index_instance;
        }
        rev.setComputeProgresses(periodicities);

        // commented out because not true when a late handler is zero
        /*
        auto progresses = rev.getComputeProgresses();
        if(progresses != periodicities) {
            throw std::logic_error("setComputeProgresses error");
        }
         */
    }
    
    if constexpr(C::has_subsampling) {
        auto & lateHandler = rev.getB();
        int n_scales = count_scales(ps);
        for(int i=1; i<n_scales; ++i) {
            // the top-most will be stepped 'base_phase' times,
            // then each scale after that will be stepped by a quarter grain size.
            // phases are cumulative, so stepping a scale also steps subsequent scales.
            static_assert(nMaxScales==4);
            switch(i) {
                case 1:
                    if constexpr (C::has_subsampling)
                    {
                        auto & inner = lateHandler.getB().getInner().getInner();
                        assert(inner.isValid());
                        assert(!inner.isZero());
                        int quarter_grain_size = inner.getA().getGranularMinPeriod() / 4;
                        for(int j=0; j<quarter_grain_size; ++j) {
                            inner.step(0);
                        }
                    }
                    else {
                        throw std::logic_error("cannot subsample 1");
                    }
                    break;
                case 2:
                    if constexpr (C::has_subsampling)
                    {
                        auto & inner = lateHandler.getB().getInner().getInner().getB().getInner().getInner();
                        assert(inner.isValid());
                        assert(!inner.isZero());
                        int quarter_grain_size = inner.getA().getGranularMinPeriod() / 4;
                        for(int j=0; j<quarter_grain_size; ++j) {
                            inner.step(0);
                        }
                    }
                    else {
                        throw std::logic_error("cannot subsample 2");
                    }
                    break;
                case 3:
                    if constexpr (C::has_subsampling)
                    {
                        auto & inner = lateHandler.getB().getInner().getInner().getB().getInner().getInner().getB().getInner().getInner();
                        assert(inner.isValid());
                        assert(!inner.isZero());
                        int quarter_grain_size = inner.getGranularMinPeriod() / 4;
                        for(int j=0; j<quarter_grain_size; ++j) {
                            inner.step(0);
                        }
                    }
                    else {
                        throw std::logic_error("cannot subsample 3");
                    }
                    break;
                default:
                    throw std::logic_error("out of bound");
            }
        }
    }
}


  
}
