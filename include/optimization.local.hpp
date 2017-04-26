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
        using Base::results;
        using Base::getValue;
        using Base::eval;
        using Base::begin;
        using Base::end;
        
        GradientDescent() { reset(); }
        
        GradientDescent(FUNCTION f) : Base(f) {
            reset();
        }
        
        int findLocalMinimum(int const n_iterations, int start_param, Value & min_val) {
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
                min_val = this->min_value;
                start_param = min_param;
            }
            
            return min_param;
        }
        
        int debug_extend(range<int> indices) {
            verify_func_exists();
            
            assert(decreasing_direction); // call run first!
            extend_results(-1, indices);
            extend_results(+1, indices);
            
            auto min_ = min_value;
            int index_min = min_param;
            for(auto v : results) {
                auto m = v.second.get_value();
                if(std::isnan(static_cast<float>(m))) {
                    continue;
                }
                if( m < min_) {
                    min_ = m;
                    index_min = v.second.index;
                    std::cout << "found better min!" << std::endl;
                }
            }
            return index_min;
        }
        
        void debug(bool logdraw = false) {
            verify_func_exists();

            auto const min_index = std::min(0, begin()->first);
            auto const end_index = std::max(9, min_param + 5);
            
            auto const ret2 = debug_extend(range<int>{min_index, end_index-1});
            
            std::vector<Value> values;
            
            int index = begin()->first;
            assert(index == min_index);
            
            for(auto const & res : results) {
                assert(index == res.second.index); // make sure we capture every index in the range
                values.push_back(res.second.get_value());
                ++index;
            }
            assert(index == end_index);
            
            constexpr auto height = 30;
            StringPlot plot(height, values.size());
            if(logdraw) {
                constexpr auto with_zero = false;
                plot.drawLog(values, default_curve_char, with_zero);
            } else {
                plot.draw(values);
            }
            plot.log();
            
            if(min_param != ret2) {
                static auto count = 0;
                ++count;
                std::cout << "mistake " << count << " : "
                << min_param << " (" << std::log(static_cast<float>(getValue(min_param))) << ") != "
                << ret2      << " (" << std::log(static_cast<float>(getValue(ret2)))      << ")" << std::endl;
            }
        }
        
    private:
        int decreasing_direction;
        int first_param;
        int min_param;
        Value min_value;

        void reset() { onFunctionChanged(); }
        
        void onFunctionChanged() override {
            Base::onFunctionChanged();
            decreasing_direction = 0;
            min_param = -1;
        }
        
        int do_run(int param) {
            Value val;
            auto res = eval(param, val);
            
            if(res == ParamState::OutOfRange) {
                // we went too far
                if(!decreasing_direction) {
                    if(first_param == param) {
                        throw std::logic_error("the first parameter passed is out of range");
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
    int findLocalMinimum(int n_iterations,
                         int from,
                         F f,
                         float & min_value)
    {
        GradientDescent<float> gd(f);
        
        return gd.findLocalMinimum(n_iterations,
                                   from,
                                   min_value);
    }
}
