/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    enum class ParamState {
        OutOfRange,
        Ok
    };
    
    namespace detail {
        
        /*
         * Finds a local minimum of a discrete function
         */
        
        struct GradientDescent {
            
            template<typename FUNCTION>
            int run(int start_param, FUNCTION f) {
                first_param = start_param;
                int next = first_param;
                while(1) {
                    next = do_run(next, f);
                    if(next == first_param) {
                        break;
                    }
                }
                return min_param;
            }
            
        private:
            int decreasing_direction = 0;
            int first_param;
            int min_param = -1;
            float min_value;
            
            
            template<typename FUNCTION>
            int do_run(int param, FUNCTION f) {
                float val;
                auto res = f(param, val);
                if(res == ParamState::OutOfRange) {
                    // we went too far
                    if(!decreasing_direction) {
                        if(first_param == param) {
                            throw std::logic_error("the first parameter passed is out of range");
                        }
                        // we want to go away from where we are
                        decreasing_direction = (first_param > param) ? 1 : -1;
                        return first_param + decreasing_direction;
                    }
                    return first_param; // we're done, we were decreasing and we found the limit
                }
                else {
                    assert(res == ParamState::Ok);
                    if(!decreasing_direction) {
                        if(param == first_param) {
                            min_value = val;
                            min_param = param;
                            return first_param + 1;
                        }
                        if(val <= min_value) {
                            // we want to continue in this direction
                            min_value = val;
                            min_param = param;
                            decreasing_direction = (first_param > param) ? -1 : 1;
                            return param + decreasing_direction;
                        }
                        else {
                            //we need to go the other way
                            decreasing_direction = (first_param > param) ? 1 : -1;
                            return first_param + decreasing_direction;
                        }
                    }
                    else {
                        assert(param != first_param);
                        if(val <= min_value) {
                            // we want to continue
                            min_value = val;
                            min_param = param;
                            assert(decreasing_direction == (first_param > param) ? -1 : 1);
                            return param + decreasing_direction;
                        }
                        else {
                            // we just passed the minimum
                            return first_param;
                        }
                    }
                }
                throw std::runtime_error("we cannot be here");
            }
        };
    }

    template<typename F>
    int findMinimun(int from, F f) {
        detail::GradientDescent gd;
        return gd.run(from, f);
    }
}
