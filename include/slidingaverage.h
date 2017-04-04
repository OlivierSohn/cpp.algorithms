/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    // TODO use advanced traits for efficient parameter passing
    
    template <class T, int KIND = OVERWRITE_INITIAL_VALUES_WITH_FIRST_FEED>
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
        cyclic<T, KIND> values_;
    };
}
