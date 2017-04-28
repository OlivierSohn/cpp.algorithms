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
        
        void plot(bool logdraw = true) const {

            if(results.empty()) {
                std::cout << "empty" << std::endl;
                return;
            }
            
            auto range_index = range<int>{
                results.begin()->first,
                results.rbegin()->first
            };
            
            constexpr auto unevaluated_val = std::numeric_limits<float>::quiet_NaN();
            std::vector<float> values(range_index.getSpan() + 1,
                                      unevaluated_val);

            for(auto const & v : results) {
                if(v.second.state == ParamState::Ok) {
                    values[v.first - range_index.getMin()] = v.second.val;
                }
            }
            
            constexpr auto height = 30;
            StringPlot p(height, values.size());
            if(logdraw) {
                p.drawLog(std::move(values));
            } else {
                p.draw(values);
            }
            p.log();
        }

        int make_exhaustive(range<int> const & r) {
            verify_func_exists();
            
            Value min_before;
            auto index_min_before = getMinValue(min_before);

            for(auto i = r.getMin(); i<=r.getMax(); ++i) {
                auto it = results.find(i);
                if(it != results.end()) {
                    continue;
                }
                Value val;
                auto res = f(i, val);
                results[i] = {i, res, val};
            }
            
            Value min_after;
            auto index_min_after = getMinValue(min_after);
            
            if(index_min_after != index_min_before) {
                using namespace std;
                static auto count = 0;
                ++count;
                cout << "mistake " << count << " : "
                << index_min_before << " (" << log(static_cast<float>(min_before)) << ") != "
                << index_min_after  << " (" << log(static_cast<float>(min_after))  << ")" << endl;
            }
            return index_min_after;
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
