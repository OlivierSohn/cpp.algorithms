
namespace imajuscule
{
    enum CyclicInitializationType
    {
        OVERWRITE_INITIAL_VALUES_WITH_FIRST_FEED,
        KEEP_INITIAL_VALUES
    };
    
    template<class T, int KIND = KEEP_INITIAL_VALUES>
    struct cyclic
    {
        using Type = T;
        using container = typename std::vector<T>;
        using iterator = typename container::iterator;
        using const_iterator = typename container::const_iterator;
        using ParameterType = T;
        
        operator container & () { return buf; }
        operator container const& () const { return buf; }
        
        const_iterator begin() const { return buf.begin();}
        iterator begin() { return buf.begin();}
        const_iterator cycleEnd() const { return it;}
        iterator cycleEnd() { return it;}
        const_iterator end() const { return buf.end();}
        iterator end() { return buf.end();}

        size_t size() const { return buf.size(); }
        
        cyclic(size_t size, ParameterType initVals)
        : initialValue(initVals), isFirstFeed(true) {
            buf.resize(size, initVals);
            it = buf.begin();
        }
        
        // not copyable
        cyclic(cyclic const &) = delete;
        cyclic& operator =(cyclic const &) = delete;
 
        // movable
        cyclic(cyclic&&o) :
        buf(std::move(o.buf)),
        initialValue(std::move(o.initialValue)),
        isFirstFeed(o.isFirstFeed)
        {
            it = buf.begin() + std::distance(o.it, o.buf.begin());
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
        
        void feed(ParameterType val) {
            *it = val;
            ++it;
            if(it == buf.end())
                it = buf.begin();
            
            if(isFirstFeed)
            {
                if(KIND == OVERWRITE_INITIAL_VALUES_WITH_FIRST_FEED)
                {
                    std::fill(buf.begin(), buf.end(), val);
                }
                isFirstFeed = false;
            }
        }
        
        void reset() {
            std::fill(buf.begin(), buf.end(), initialValue);
            it = buf.begin();
            isFirstFeed = true;
        }
        
    private:
        container buf;
        iterator it;
        T initialValue;
        bool isFirstFeed:1;
    };
    
    template<typename C, typename T = typename C::Type>
    std::vector<T> to_vector(C & cyclic_) {
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
