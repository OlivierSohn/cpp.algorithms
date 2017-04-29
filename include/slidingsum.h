/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    template <typename CYCLIC, typename T = typename CYCLIC::value_type>
    T computeMaxSlidingSum(CYCLIC const & cycle, int slide_width) {
        assert(slide_width > 0);
        assert(!cycle.empty());
        if(cycle.empty()) {
            return {};
        }
        int n_full = slide_width / cycle.size();
        T fulls_sum = n_full ? (n_full * std::accumulate(cycle.begin(), cycle.end(), T{})) : 0;
        
        int n_partials = slide_width - n_full * cycle.size();
        assert(n_partials >= 0);
        assert(n_partials < cycle.size());
        if(0 == n_partials) {
            return fulls_sum;
        }
        
        auto cycle_end =  cycle.end();
        auto it = cycle.begin();
        auto end = it + n_partials;
        T max_partial = std::accumulate(it, end, T{});
        T cur_sum = max_partial;
        do {
            cur_sum += *end - *it;  // optimization to not redo the whole sum at each step
            
            max_partial = std::max(max_partial,
                                   cur_sum);
            ++end;
            if(end == cycle_end) {
                end = cycle.begin();
            }
            ++it;
        }
        while(it != cycle_end);
        
        return fulls_sum + max_partial;
    }
    
    template <typename CYCLIC>
    void compute_phased_sum(CYCLIC const & cycle, int const phase, int const n_iterators, CYCLIC & result) {
        assert(phase >= 0);
        assert(n_iterators > 1);
        using T = typename CYCLIC::value_type;
        result.resize(cycle.size());
        auto it = result.begin();
        auto end = result.end();
        auto begin_cycle = cycle.begin();
        int ph = 0;
        auto sz = cycle.size();
        
        for(; it!=end; ++it, ++ph) {
            *it = {};
            auto ph_it = ph;
            for(int n = 0; n<n_iterators; ++n, ph_it += phase) {
                while(ph_it >= sz) {
                    ph_it -= sz;
                }
                assert(ph_it >= 0);
                *it += *(begin_cycle + ph_it);
            }
        }
    }
}
