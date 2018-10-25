/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    
    /*
     * Naive gradient descent of a discrete function
     */
    
    template<typename Value>
    struct GradientDescent : public SearchMin<Value> {
        using Base = SearchMin<Value>;
        using Result = typename Base::Result;
        using FUNCTION = typename Base::FUNCTION;
        
        using Base::verify_func_exists;
        
        using Base::f;
        using Base::make_exhaustive;
        using Base::plot;
        using Base::results;
        using Base::getValue;
        using Base::eval;
        using Base::begin;
        using Base::end;
        
        GradientDescent() { reset(); }
        
        GradientDescent(FUNCTION f) : Base(f) {
            reset();
        }
        
        Optional<int> findLocalMinimum(int const n_iterations, int start_param, Optional<Value> & min_val) {
            verify_func_exists();
            
            for(int i=0; i<n_iterations; ++i) {
                reset(); // this->min_param is reset here but the information persists in the results map
                
                first_param = start_param;
                
                int next = first_param;
                while(1) {
                    next = do_run(next);
                    if(next == first_param) {
                        break;
                    }
                }
              if(min_value) {
                min_val = *min_value;
              }
              if(min_param) {
                start_param = *min_param;
              }
            }
            
            return min_param;
        }
        
        void debug(bool logdraw = false) {
            verify_func_exists();

            plot(logdraw); // to see result before we add points

            auto const min_index = std::min(0, begin()->first);
          auto const end_index = std::max(9, (min_param?*min_param:0) + 5);
            
            auto const ret2 = make_exhaustive(range<int>{min_index, end_index-1});
            plot(logdraw);                        
        }
        
    private:
        int decreasing_direction;
        int first_param;
        Optional<int> min_param;
        Optional<Value> min_value;

        void reset() { onFunctionChanged(); }
        
        void onFunctionChanged() override {
            Base::onFunctionChanged();
            decreasing_direction = 0;
          min_param = {};
        }
        
        int do_run(int const param) {
            Value val;
            auto res = eval(param, val);
            
            if(res == ParamState::OutOfRange) {
                // we went too far
                if(!decreasing_direction) {
                    if(first_param == param) {
                      LG(INFO,"the first parameter passed is out of range");
                      return first_param;
                    }
                    // we want to go away from where we are
                    decreasing_direction = (first_param > param) ? 1 : -1;
                    return first_param + decreasing_direction;
                }
                return first_param; // we're done, we were decreasing and we found the limit
            }
            assert(res == ParamState::Ok);
            if(!decreasing_direction) {
                if(param == first_param) {
                    min_value = val;
                    min_param = param;
                    return first_param + 1; // + 1 is arbitrary it could have been - 1 too
                }
                if(val <= min_value) {
                    // we want to continue in this direction
                    min_value = val;
                    min_param = param;
                    decreasing_direction = (first_param > param) ? -1 : 1;
                    return param + decreasing_direction;
                }
                else {
                    //we need to go the other way
                    decreasing_direction = (first_param > param) ? 1 : -1;
                    return first_param + decreasing_direction;
                }
            }
            else {
                assert(param != first_param);
                if(val <= min_value) {
                    // we want to continue
                    min_value = val;
                    min_param = param;
                    assert(decreasing_direction == (first_param > param) ? -1 : 1);
                    return param + decreasing_direction;
                }
                else {
                    // we just passed the minimum
                    return first_param;
                }
            }
        }
        
        void extend_results(Result r, int const increment, range<int> indices) {
            while(1) {
                r.index += increment;
                if(!indices.contains(r.index)) {
                    break;
                }
                r.state = f(r.index, r.val);
                assert(results.count(r.index) == 0);
                results[r.index] = r;
            }
        }
        
        void extend_results(int const increment, range<int> indices) {
            assert(!results.empty());
            
            auto res = (increment > 0) ? results.rbegin()->second : results.begin()->second;
            
            extend_results(res, increment, indices);
        }
    };
    
    /*
     * Finds the local minimum of a discrete function
     *
     * @param n_iterations : the number of times the gradient descent is redone (pass > 1 only if succeessive f calls return different values)
     * @param from : the 'x index' at which gradient descent starts
     * @param f : the function
     * @return min_value : the 'y min_value'
     *
     * @return the 'x index' at which the 'y min_value' was found
     */
    template<typename F>
    Optional<int> findLocalMinimum(int n_iterations,
                         int from,
                         F f,
                         Optional<float> & min_value)
    {
        GradientDescent<float> gd(f);
        
        return gd.findLocalMinimum(n_iterations,
                                   from,
                                   min_value);
    }
}
