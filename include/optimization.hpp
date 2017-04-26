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
    
    /*
     Search over integral parameters
     */
    
    template<typename Value>
    struct SearchMin {
        using FUNCTION = std::function<ParamState(int, Value&)>;
        
        virtual ~SearchMin() = default;
        
        SearchMin() = default;
        SearchMin(FUNCTION f) : f(f) {
        }
        
        void stealFrom(SearchMin && other) {
            f = std::move(other.f);
            results = std::move(results);
        }
        
        virtual void onFunctionChanged() {
        }
        
        void verify_func_exists() {
            if(f) {
                return;
            }
            throw std::logic_error("cannot run a gradient descent when there is no function");
        }
        
        template<typename F>
        void setFunction(F func) {
            if(f) {
                throw std::logic_error("cannot change the function of a gradient descent");
            }
            f = func;
            onFunctionChanged();
        }

    protected:
        FUNCTION f;
        
        struct Result {
            Value get_value() const {
                if(state == ParamState::Ok) {
                    return val;
                }
                else {
                    return {};
                }
            }
            
            Value feed(ParamState s, Value value) {
                if(state != s) {
                    throw std::runtime_error("different param states for the same index");
                }
                if(s == ParamState::OutOfRange) {
                    return {};
                }
                val = std::min(val, value);
                return val;
            }
            
            int index;
            ParamState state;
            Value val;
        };
        std::map<int, Result> results;
    };
}
