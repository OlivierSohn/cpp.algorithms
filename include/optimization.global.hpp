/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

#define DEBUG_GRADIENT_DESCENT 1

namespace imajuscule
{
    /*
     * Searches a min value within an integral parameter range
     */
    template<typename Value>
    struct RangeSearch : public SearchMin<Value> {
        using Base = SearchMin<Value>;
        using FUNCTION = typename Base::FUNCTION;
        
        using Result = typename Base::Result;
        
        using Base::verify_func_exists;
        
        using Base::f;
        using Base::results;
        using Base::getMinValue;
        
        RangeSearch(range<int> admissible_params, FUNCTION f) :
        Base(f),
        admissible_params(admissible_params) {
            if(admissible_params.empty()) {
                throw std::logic_error("empty range");
            }
        }
        
        template<typename F>
        int findMinimum(F stop_search, Value & min_value) {
            
            int param;
            while(next(param)) {
                assert(results.find(param) == results.end());

                auto res = f(param, min_value);
                results[param] = {param, res, min_value};
                if(res == ParamState::Ok && stop_search(min_value)) {
                    break;
                }
            }

            return getMinValue(min_value);
        }
                
    private:
        range<int> admissible_params;
        
        bool next(int & param) const {
            int candidate_fwd = std::uniform_int_distribution<int> {
                admissible_params.getMin(),
                admissible_params.getMax()
            }(mersenne<SEEDED::No>());
            
            auto fwd = results.lower_bound(candidate_fwd);
            auto bwd = fwd;
            auto candidate_bwd = candidate_fwd;
            bool do_fwd = true;
            bool do_bwd = (bwd != results.end());

            do {
                if(do_fwd) {
                    if(fwd == results.end()) {
                        if(candidate_fwd <= admissible_params.getMax()) {
                            param = candidate_fwd;
                            return true;
                        }
                        do_fwd = false;
                    }
                    else {
                        if(fwd->first != candidate_fwd) {
                            param = candidate_fwd;
                            return true;
                        }
                        ++fwd;
                    }
                    ++candidate_fwd;
                }
                
                if(do_bwd) {
                    --candidate_bwd;
                    if(bwd == results.begin()) {
                        if(candidate_bwd >= admissible_params.getMin()) {
                            param = candidate_bwd;
                            return true;
                        }
                        do_bwd = false;
                    }
                    else {
                        bwd--;
                        if(bwd->first != candidate_bwd) {
                            param = candidate_bwd;
                            return true;
                        }
                    }
                }
            } while(do_fwd || do_bwd);
            
            return false;
        };
    };
    
    /*
     * finds the global minimum or a minimum smaller than threshold, in the given range.
     */
    template<typename F>
    int findGlobalMinimumOrSmallerThan(range<int> r, float threshold, F f, float & min_) {
        RangeSearch<float> rs(r, f);
        
        return rs.findMinimum([threshold](auto cost){
            return cost < threshold;
        }, min_);
    }
    template<typename F>
    int findGlobalMinimum(range<int> r, F f, float & min_) {
        RangeSearch<float> rs(r, f);
        
        return rs.findMinimum([](auto){
            return false;
        }, min_);
    }

}
