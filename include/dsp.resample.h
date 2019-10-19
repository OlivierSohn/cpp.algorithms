
namespace imajuscule::audio {
    
    /*
     Fixed-rate resampling using linear interpolation.
     
     Using linear interpolation for resampling audio is very fast, but also leads to very low quality results.
     Hence, this will be replaced by windowed-sinc low pass filtering.
     */
    template<typename Reader, typename Buffer>
    /*
     Reader shoud have the following methods:
     
     countChannels :: int
     ReadAsOneFloat<VAL> :: VAL
     HasMore :: bool
     */
    void resampleLinear(Reader & reader, Buffer & buf, double sample_rate)
    {
        using VAL = typename Buffer::value_type;
        
        buf.clear();
        auto const n_channels = reader.countChannels();
        double frameRatio = reader.getSampleRate() / sample_rate;
        if(!frameRatio) {
            return;
        }
        int target_size = static_cast<int>(1 + (reader.countFrames()-1) / frameRatio) * reader.countChannels();
        if(target_size <= 0) {
            return;
        }
        buf.resize(target_size);
        auto it = buf.begin();
        auto end = buf.end();
        
        // represents current and next frames
        std::array<std::vector<VAL>, 2> vecs;
        for(auto & v : vecs) {
            v.resize(n_channels, {});
        }
        auto & cur = vecs[0];
        auto & next = vecs[1];

        auto readFrames = [n_channels, &reader](auto & cur, auto & next, int64_t countFrames) {
            bool the_end = false;
            for(int64_t i=0; i<countFrames; ++i) {
                for(int64_t j=0; j<n_channels; ++j) {
                    cur[j] = next[j];
                    if(reader.HasMore()) {
                        next[j] = reader.template ReadAsOneFloat<VAL>();
                    }
                    else {
                        the_end = true;
                    }
                }
            }
            return the_end;
        };
        
        // invariants :
        // cur  has index curFrameRead
        // next has index curFrameRead+1
        int64_t curFrameRead = -2;
        int64_t frameIdx = 0;
        double where = 0.;
        int64_t countReads = 1;

        constexpr double epsilon = std::numeric_limits<double>::epsilon();
        while(it != end) {
            if(countReads > 0) {
                curFrameRead += countReads;
                if(readFrames(cur, next, countReads)) {
                    return;
                }
            }
            double const next_ratio = where - curFrameRead;
            if(next_ratio >= 1.) {
                assert(next_ratio <= 1. + epsilon);
                for(int j=0; j<n_channels; ++j, ++it) {
                    *it = next[j];
                }
            }
            else {
                assert(next_ratio <= 1.);
                assert(next_ratio >= 0.);
                double const cur_ratio = 1. - next_ratio;
                
                for(int j=0; j<n_channels; ++j, ++it) {
                    *it = cur_ratio * cur[j] + (next_ratio * next[j]);
                }
            }
            ++frameIdx;
            where = frameRatio*frameIdx;
            // we subtract epsilon to NOT request one more read if we are at 1.0000000001
            countReads = static_cast<int64_t>(where*(1-epsilon)-curFrameRead);
        }
    }
    
    struct InterlacedBuffer {
        using element_type = double;
        
        template<typename Reader>
        InterlacedBuffer(Reader & reader, double sample_rate)
        {
            resampleLinear(reader, buf, sample_rate);
            
            finalize();
        }
        
        InterlacedBuffer(std::vector<element_type> const & v, int nchannels)
        : buf(v)
        , nchannels(nchannels)
        {}
        
        std::size_t countFrames() const {
            if(nchannels == 0) {
                return 0;
            }
            return buf.size() / nchannels;
        }
        
        int countChannels() const { return nchannels; }
        
        std::vector<element_type> const & getBuffer() const { return buf; }
        
        element_type getMaxAbsValue() const { return maxAbsValue; }
        
        
        // https://stackoverflow.com/a/27216842/3940228
        std::size_t hashCode() const {
            auto hashFunc = [](std::size_t &seed, element_type val){
                std::size_t &i = reinterpret_cast<std::size_t&>(val);
                static_assert(sizeof(std::size_t) == sizeof(element_type));
                seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            };
            
            // initialize with something
            std::size_t val = buf.size() / nchannels;
            
            // mix with the channel count
            hashFunc(val, nchannels);
            
            // mix with the buffer
            for(auto const & v:buf) {
                hashFunc(val, v);
            }
            
            return val;
        }
    private:
        std::vector<element_type> buf;
        int nchannels = 0;
        element_type maxAbsValue = 0;
        
        void finalize() {
            for(auto v:buf) {
                maxAbsValue = std::max(maxAbsValue, std::abs(v));
            }
        }
    };
    
} // namespace imajuscule::audio
