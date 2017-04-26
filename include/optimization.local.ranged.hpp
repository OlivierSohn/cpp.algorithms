/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    /*
     * This algorithm was introduced for functions with randomicity,
     * when computing the derivative using close neighbours leads to big errors.
     *
     * A Naive approach would have been to do a binary search abd compute the derivative using neighbours at each step.
     * Here we first use "far away" neighbours to do this computation and progressively reduce the reach.
     * Complexity-wise, there is just a small constant factor in front of the Log N of the binary search to pay.
     */
    
    template<typename Value>
    struct RangedGradientDescent : public SearchMin<Value> {
        using Base = SearchMin<Value>;
        using Result = typename Base::Result;
        using FUNCTION = typename Base::FUNCTION;
        
        using Base::verify_func_exists;
        
        using Base::f;
        using Base::results;
        using Base::eval;
        using Base::begin;
        using Base::end;
        
        RangedGradientDescent() {}
        
        RangedGradientDescent(FUNCTION f) : Base(f) {
        }
        
        int findLocalMinimum(int const n_iterations, range<int> const r, Value & min_val) {
            verify_func_exists();
            
            bool first = true;
            Result result;
            for(int i=0; i<n_iterations; ++i) {
                std::array<Result, 5> a;
                
                bool loop = true;
                
                eval(a[0], r.getMin());
                
                if(r.getCenter() != r.getMin()) {
                    eval(a[2], r.getCenter());
                }
                else {
                    loop = false;
                }
                
                if(r.getMax() != r.getCenter()) {
                    eval(a[4], r.getMax());
                }
                else{
                    loop = false;
                }
                
                if(loop) {
                    while(1) {
                        bool stop = false;
                        if(a[2].index - a[0].index >= 2) {
                            eval(a[1], (a[0].index + a[2].index) / 2);
                        }
                        else {
                            a[1].state = ParamState::NotEvaluated;
                            stop = true;
                        }
                        
                        if(a[4].index - a[2].index >= 2) {
                            eval(a[3], (a[2].index + a[4].index) / 2);
                        }
                        else {
                            a[3].state = ParamState::NotEvaluated;
                            stop = true;
                        }

                        if(stop) {
                            break;
                        }

                        for(auto const & e:a) {
                            if(e.state == ParamState::NotEvaluated) {
                                throw std::logic_error("by now all should be evaluated");
                            }
                        }
                        
                        // find index in array of min
                        int m = 0;
                        for(int i=1; i<5; ++i) {
                            if(a[i].val < a[m].val) {
                                m = i;
                            }
                        }
                        
                        // zoom on min
                        
                        if(m <= 1) {
                            a[4] = a[2];
                            a[2] = a[1];
                            // a[0] = a[0];
                        }
                        else if(m == 2) {
                            a[4] = a[3];
                            //a[2] = a[2];
                            a[0] = a[1];
                        }
                        else if(m >= 3) {
                            a[0] = a[2];
                            a[2] = a[3];
                            //a[4] = a[4];
                        }
                        else {
                            assert(0);
                        }
                    }
                }
                
                for(auto const & r : a) {
                    if(r.state!=ParamState::Ok) {
                        continue;
                    }
                    if(first) {
                        first = false;
                        result = r;
                        continue;
                    }
                    if( static_cast<float>(r.val) < static_cast<float>(result.val)) {
                        result = r;
                    }
                }
            }
            
            if(result.state != ParamState::Ok) {
                throw std::logic_error("could not range optimize");
            }
            min_val = result.val;
            return result.index;
        }
        
    private:
        
        void eval(Result & res, int param) {
            res.index = param;
            res.state = eval(res.index, res.val);
            if(res.state != ParamState::Ok) {
                throw std::logic_error("out of bounds results are not allowed in ranged optimization");
            }
        }
        
    };
    
    /*
     * Finds the local minimum of a discrete function
     *
     * @param n_iterations : the number of times the gradient descent is redone (pass > 1 only if successive f calls return different values)
     * @param from : the 'x index' at which gradient descent starts
     * @param f : the function
     * @return min_value : the 'y min_value'
     *
     * @return the 'x index' at which the 'y min_value' was found
     */
    template<typename F>
    int findRangedLocalMinimum(int n_iterations,
                               range<int> r,
                               F f,
                               float & min_value)
    {
        RangedGradientDescent<float> gd(f);
        
        return gd.findLocalMinimum(n_iterations,
                                   r,
                                   min_value);
    }
}
