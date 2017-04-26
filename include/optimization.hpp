/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    enum class ParamState {
        Ok,
        OutOfRange,
        NotEvaluated,
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
        
        auto begin() const { return results.begin(); }
        auto end() const { return results.end(); }

        Value getValue(int key) const {
            auto it = results.find(key);
            if( it == results.end() ) {
                throw std::logic_error("key not found in gradient descent results");
            }
            return it->second.get_value();
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
            ParamState state = ParamState::NotEvaluated;
            Value val;
        };
        std::map<int, Result> results;
        
        int getMinValue(Value & min_value) const
        {
            bool first = true;
            int res{-1}; // initialization just to silence warning
            for(auto const & p : results) {
                auto const & r = p.second;
                if(r.state != ParamState::Ok) {
                    continue;
                }
                if(first) {
                    first = false;
                    min_value = r.val;
                    res = r.index;
                    continue;
                }
                if(static_cast<float>(r.val) < static_cast<float>(min_value)) {
                    min_value = r.val;
                    res = r.index;
                }
            }
            
            if(first) {
                throw std::logic_error("no value worked");
            }
            
            return res;
        }
        
        ParamState eval(int param, Value & val) {
            auto res = f(param, val);
            
            auto it = results.find(param);
            if(it == results.end()) {
                results[param] = {param, res, val};
            }
            else {
                // use information of previous runs : take the minimum
                val = it->second.feed(res, val);
            }
            return res;
        }
    };
}
