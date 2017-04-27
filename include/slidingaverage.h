/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    // TODO use advanced traits for efficient parameter passing
    
    template <class T, CyclicInitialization Init = CyclicInitialization::FIRST_FEED>
    struct slidingAverage
    {
        using ParameterType = T;
        
        slidingAverage(slidingAverage&&) = default;
        slidingAverage& operator=(slidingAverage&&) = default;

        slidingAverage(size_t size) :
        values_(size, T{})
        {}
        
        void feed(ParameterType val) {
            values_.feed(val);
        }
        
        T compute() const {
            return std::accumulate(values_.begin(), values_.end(), T{}) / safe_cast<float>(values_.size());
        }
        
        void reset() {
            values_.reset();
        }
        
        auto size() const {
            return values_.size();
        }
        
    private:
        cyclic<T, Init> values_;
    };
    
    
    /*
     * the window function parameter is between 0 (out of the window) and 1 (middle of the window)
     * the 0 value is never used
     */
    template <class T, CyclicInitialization Init = CyclicInitialization::FIRST_FEED>
    struct slidingWindowedAverage
    {
        using ParameterType = T;
        
        slidingWindowedAverage(slidingWindowedAverage&&) = default;
        slidingWindowedAverage& operator=(slidingWindowedAverage&&) = default;
        
        template<typename WINDOW>
        slidingWindowedAverage(size_t size, WINDOW f) :
        values_(size, T{}),
        window(size)
        {
            assert(std::abs(f(1.f) - 1.f) < 0.0001f);

            float n_steps;
            float mid_idx;
            if(size%2) {
                // even
                n_steps = 1.f + (size-1)/2;
                mid_idx = n_steps;
            }
            else {
                // odd
                n_steps = .5f + size/2;
                mid_idx = n_steps;
            }
            
            auto it = window.begin();
            auto end = window.end();
            auto step = 1;
            for(; it!=end; ++it, ++step) {
                float real_step = step;
                if(real_step > mid_idx) {
                    real_step = 2.f * mid_idx - real_step;
                }
                
                *it = f(real_step / n_steps);
                assert(*it <= 1.f);
            }
        }
        
        void feed(ParameterType val) {
            values_.feed(val);
        }
        
        T compute() const {
            T sum{};
            auto it = window.begin();
            values_.for_each([&](auto v) {
                sum += v * *it;
                ++it;
            });
            return sum / safe_cast<float>(values_.size());
        }
        
        void reset() {
            values_.reset();
        }
        
        auto size() const {
            return values_.size();
        }
        
        auto getWindowAverage() const {
            return std::accumulate(window.begin(), window.end(), T{}) / window.size();
        }
    private:
        cyclic<T, Init> values_;
        std::vector<T> window;
    };
}
