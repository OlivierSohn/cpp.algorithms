/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    enum class CyclicInitialization
    {
        FIRST_FEED,
        INITIAL_VALUES
    };
    
    template<class T, CyclicInitialization Init = CyclicInitialization::INITIAL_VALUES>
    struct cyclic
    {
        using container = typename std::vector<T>;
        using iterator = typename container::iterator;
        using const_iterator = typename container::const_iterator;
        using value_type = T;
        using ParameterType = T;
        
        operator container & () { return buf; }
        operator container const& () const { return buf; }
        
        auto & back() { return buf.back(); }
        
        auto begin() const { return buf.begin();}
        auto end() const { return buf.end();}
        auto begin() { return buf.begin();}
        auto end() { return buf.end();}
        auto rbegin() const { return buf.rbegin();}
        auto rend() const { return buf.rend();}
        auto rbegin() { return buf.rbegin();}
        auto rend() { return buf.rend();}
        const_iterator cycleEnd() const { return it;}
        auto cycleEnd() { return it;}

        auto size() const { return buf.size(); }
        auto empty() const { return buf.empty(); }
        
        cyclic() : isFirstFeed(true) {
            it = buf.begin();
        }
        
        cyclic(size_t size, ParameterType initVals = {})
        : initialValue(initVals), isFirstFeed(true) {
            buf.resize(size, initVals);
            it = buf.begin();
        }
        
        // not copyable
        cyclic(cyclic const &) = delete;
        cyclic& operator =(cyclic const &) = delete;
 
        // movable
        cyclic(cyclic&&o) :
        initialValue(std::move(o.initialValue)),
        isFirstFeed(o.isFirstFeed)
        {
            auto d =  std::distance(o.buf.begin(), o.it);
            buf = std::move(o.buf);
            it = buf.begin() + d;
        }
        
        cyclic(container i)
        : buf(std::move(i)), isFirstFeed(false) {
            it = buf.begin();
            
        }
        
        cyclic & operator =(cyclic && o) {
            if(this != &o) {
                buf = std::move(o.buf);
                initialValue = std::move(o.initialValue);
                isFirstFeed = o.isFirstFeed;
                
                it = buf.begin() + std::distance(o.it, o.buf.begin());
            }
            return *this;
        }
        
        void resize(int sz) {
            auto dist = std::distance(buf.begin(), it);
            buf.resize(sz, initialValue);
            it = buf.begin() + dist;
        }
        
        void grow(ParameterType val) {
            auto dist = std::distance(buf.begin(), it);
            buf.push_back(val);
            it = buf.begin() + dist;
        }
        
        void feed(ParameterType val) {
            if(isFirstFeed) {
                if(Init == CyclicInitialization::FIRST_FEED) {
                    std::fill(buf.begin(), buf.end(), val);
                }
                isFirstFeed = false;
            }
            
            *it = std::move(val);
            advance();
        }
        
        void setByIndex(int i) {
            it = buf.begin() + i;
            assert(it < buf.end());
        }
        
        int getIndex() const {
            return std::distance(buf.begin(), cycleEnd());
        }
        
        void advance() {
            ++it;
            if(it == buf.end()) {
                it = buf.begin();
            }
        }
        
        void go_back() {
            if(it == buf.begin()) {
                it = buf.end();
            }
            --it;
        }
        
        void reset() {
            std::fill(buf.begin(), buf.end(), initialValue);
            it = buf.begin();
            isFirstFeed = true;
        }
        
        void erase(iterator dit) {
            if(it >= dit) {
                if(it != buf.begin()) {
                    --it;
                }
            }
            buf.erase(dit);            
        }
        
        template<typename F>
        void for_each(F f) const {
            auto start = cycleEnd();
            std::for_each(start, end(), f);
            std::for_each(begin(), start, f);
        }
        
        template<typename F>
        void for_each_bkwd(F f) const {
            auto start = std::reverse_iterator<const_iterator>(cycleEnd());
            std::for_each(start, rend(), f);
            std::for_each(rbegin(), start, f);
        }
        
        auto const & get_backward(int index_backward) const {
            int end_index = getIndex();
            auto real_index = end_index - 1 - index_backward;
            auto sz = buf.size();
            while(real_index < 0) {
                real_index += sz;
            }
            while(real_index >= sz) {
                real_index -= sz;
            }
            return buf[real_index];
        }
        
        void toBack() {
            it = buf.end()-1;
        }
        
    private:
        container buf;
        iterator it;
        T initialValue;
        bool isFirstFeed:1;
    };
    
    template<typename T>
    static std::ostream& operator<<(std::ostream& os, const cyclic<T>& c)
    {
        for(auto const & v : c) {
            os << v << " ";
        }
        return os;
    }

    template<typename C, typename T = typename C::value_type>
    std::vector<T> to_vector(C const & cyclic_) {
        std::vector<T> vec;
        vec.reserve(cyclic_.size());
        auto cycleStart = cyclic_.cycleEnd();
        auto end = cyclic_.end();
        auto begin = cyclic_.begin();

        std::move(cycleStart, end, std::back_inserter(vec));
        std::move(begin, cycleStart, std::back_inserter(vec));
        return std::move(vec);
    }
}
