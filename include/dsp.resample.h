#include <unordered_set>

namespace imajuscule
{
    struct Compare7Digits {
        static constexpr double mult = 10000000;
        constexpr bool operator()( const double& lhs, const double& rhs ) const {
            return static_cast<int64_t>(mult * lhs) < static_cast<int64_t>(mult * rhs);
        }
    };
    
    // designed to be constructed from doubles in the range [0,1)
    struct AlmostDouble {
        static constexpr double mult = 10000000.;
        
        AlmostDouble(double d)
        : key_(static_cast<int64_t>(mult * d))
        {}
        
        int64_t key() const { return key_; }
        
        bool operator == (AlmostDouble const & other) const {
            return key_ == other.key_;
        }

    private:
        int64_t key_;
    };
} // namespace imajuscule::audio

namespace std {
  template <> struct hash<imajuscule::AlmostDouble>
  {
    size_t operator()(const imajuscule::AlmostDouble & x) const
    {
        return hash<int64_t>()(x.key());
    }
  };
} // namespace std

#define DEBUG_RESAMPLING 0

namespace imajuscule::audio {

    
    /*
     [Deprecated, use resampleSinc instead]
     
     Fixed-rate resampling using linear interpolation.
     
     Using linear interpolation for resampling audio is very fast, but also leads to very low quality results.

     Hence, please use resampleSinc instead.
     */
    /*
     'Reader' is expected to have the following interface:
     Properties:
     - countChannels :: int
     - countFrames :: int
     - getSampleRate :: number
     Iteration:
     - ReadAsOneFloat<VAL> :: VAL
     - HasMore :: bool
     */
    template<typename Reader, typename Buffer>
    void resampleLinear(Reader & reader, Buffer & buf, double sample_rate)
    {
        using VAL = typename Buffer::value_type;
        
        buf.clear();
        auto const n_channels = reader.countChannels();
        double const Fs_over_FsPrime = reader.getSampleRate() / sample_rate;
        if(!Fs_over_FsPrime) {
            return;
        }
        int target_size = static_cast<int>(1 + (reader.countFrames()-1) / Fs_over_FsPrime) * n_channels;
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
                assert(next_ratio <= 1. + 2*(where*epsilon));
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
            where = Fs_over_FsPrime*frameIdx;
            // we subtract epsilon to NOT request one more read if we are at 1.0000000001
            countReads = static_cast<int64_t>(where*(1-epsilon)-curFrameRead);
        }
    }
    
    struct KaiserWindow {
        // https://en.wikipedia.org/wiki/Kaiser_window
        KaiserWindow(double shape = 16.)
        : squareShape(shape * shape)
        , oneOverDenom(1.0 / zeroethOrderBessel( shape * shape ))
        {}
        
        // K in [-1,1]
        double getAt( const double K )
        {
            const double arg = 1.0 - (K * K);
            return zeroethOrderBessel( squareShape * arg ) * oneOverDenom;
        }
    private:
        static double zeroethOrderBessel( double squareX )
        {
            constexpr double eps = 0.000001;
            
            const double squareX_Over_4 = squareX * 0.25;
            
            //  initialize the series term for m=0 and the result
            double besselValue = 0;
            double term = 1;
            double m = 0;
            
            //  accumulate terms as long as they are significant
            while(term  > eps * besselValue)
            {
                besselValue += term;
                
                //  update the term
                ++m;
                term *= squareX_Over_4 / (m*m);
            }
            
            return besselValue;
        }
        
        const double squareShape;
        const double oneOverDenom;
    };

    struct ResampleSincStats {
        std::chrono::steady_clock::duration dt_resample = {};
        std::chrono::steady_clock::duration dt_read_source = {};

        friend std::ostream& operator<<(std::ostream& os, const ResampleSincStats& dt);
        
    };
    std::ostream & operator << (std::ostream & os, const ResampleSincStats& s);

    // 0 o/n 2*o/n ...
    
    // pgcd * a = o
    // pgcd * b = n
    
    // 0 a/b
    
    //     o     n
    // 96000 41000   300
    //   320   147
    //
    // n tel que  n * 320/147 = int     n*320 = m*147
    // n=147 : ca marche
    //
    // si on utilise 4 threads, on porurait diviser [0 147) en 4 intervalles
    
    template<typename OriginalBuffer, typename ResampledBuffer>
    std::chrono::steady_clock::duration resampleSincBuffer(OriginalBuffer const & original,
                                                           int const n_channels,
                                                           double const Fs_over_FsPrime,
                                                           ResampledBuffer & resampled)
    {
        using namespace std::chrono;
        using namespace profiling;
        std::chrono::steady_clock::duration dt;
        Timer<steady_clock> t(dt);

        resampled.clear();
        if(!n_channels) {
            return dt;
        }
        auto const n_orig_frames = original.size() / n_channels;

        if(!Fs_over_FsPrime) {
            return dt;
        }
        int target_size = static_cast<int>(1 + (n_orig_frames-1) / Fs_over_FsPrime) * n_channels;
        if(target_size <= 0) {
            return dt;
        }
        if(std::abs(Fs_over_FsPrime - 1.) < 1e-8) {
            resampled = original;
            return dt;
        }
        
        // using notations found in https://ccrma.stanford.edu/~jos/resample/resample.pdf
        double const norm = std::min(1., 1./Fs_over_FsPrime);
        const auto one_over_zeroCrossingDistance = std::max(1., Fs_over_FsPrime);
        // half size, excluding the central point, hence the number of non-zero slots is:
        // 1 + 2*windowHalfSize
        constexpr auto num_zerocrossings_half_window = 512;
        constexpr double one_over_num_zerocrossings_half_window = 1./num_zerocrossings_half_window;
        const auto windowHalfSize = num_zerocrossings_half_window / one_over_zeroCrossingDistance;
        int64_t const N = static_cast<int>(windowHalfSize) + 1;
        KaiserWindow window;
        
        auto hs = [norm,
                   num_zerocrossings_half_window,
                   one_over_num_zerocrossings_half_window,
                   one_over_zeroCrossingDistance,
                   &window]
        (double t) { // zero crossings every one_over_zeroCrossingDistance
            auto sinc = [num_zerocrossings_half_window,
                         &window](double x) { // zero crossings every x
                bool inWindow = (x > -num_zerocrossings_half_window) && (x < num_zerocrossings_half_window);
                if(!inWindow) {
                    return 0.;
                }
                auto winCoeff = window.getAt(std::abs(x) * one_over_num_zerocrossings_half_window);
                double sinXOverX;
                if(!x) {
                    sinXOverX = 1.;
                }
                else {
                    double pi_x = M_PI * x;
                    sinXOverX = std::sin(pi_x)/pi_x;
                }
                return winCoeff * sinXOverX;
            };
            return norm * sinc(t * one_over_zeroCrossingDistance);
        };
        
        auto foreachResampled = [target_size, Fs_over_FsPrime](auto f){
            for(int64_t frameIdx = 0; frameIdx < target_size; ++frameIdx) {
                // we are computing frame 'frameIdx' of the resampled signal.
                
                // 'where' is the corresponding location in the original signal.
                double const where = Fs_over_FsPrime*frameIdx; // >= 0
                
                // in the original signal, 'where' is between integers 'discreteSampleBefore' and 'discreteSampleBefore'+1
                int64_t const discreteSampleBefore = static_cast<int64_t>(where);
                
                double const P = where-discreteSampleBefore; // in [0,1)

                if(!f(frameIdx, P)) {
                    return;
                }
            }
        };
        
        int64_t nDifferentLocations = 1;
        foreachResampled([&nDifferentLocations](int64_t frameIdx, double P) {
            if(frameIdx && !AlmostDouble(P).key()) {
                nDifferentLocations = frameIdx;
                return false;
            };
            return true;
        });
        
        // we memoize the results
        std::unordered_map<AlmostDouble, std::unique_ptr<std::vector<double>>> results;
        results.reserve(nDifferentLocations);
        // x between 0 and 1
        auto getResults = [&results, hs, N](double x) -> std::vector<double> const & {
            assert(x >= 0.);
            assert(x <= 1.);
            AlmostDouble key(x);
            auto it = results.find(key);
            if(it != results.end()) {
                return *it->second;
            }
            auto res = std::make_unique<std::vector<double>>();
            auto & v = *res;
            v.reserve(2+2*N);
            int64_t firstIdx = -N;
            int64_t last_idx = N+1; // included
            for(int64_t i=firstIdx; i<= last_idx; ++i) {
                v.push_back(hs(x-i));
            }
            results.insert({key, std::move(res)});
            return v;
        };
        
        resampled.resize(target_size, {});
        auto it = resampled.data();
        auto end = resampled.data() + resampled.size();
        for(int64_t frameIdx = 0; it != end; ++frameIdx, it += n_channels) {
            // we are computing frame 'frameIdx' of the resampled signal.
            
            // 'where' is the corresponding location in the original signal.
            double const where = Fs_over_FsPrime*frameIdx; // >= 0
            
            // in the original signal, 'where' is between integers 'discreteSampleBefore' and 'discreteSampleBefore'+1
            int64_t const discreteSampleBefore = static_cast<int64_t>(where);
            
            double const P = where-discreteSampleBefore; // in [0,1)
            auto const & vres = getResults(P);
            // when the resampling frequency is constant, we can cache the hs results array with key 'P'
            // this array contains results for hs(P-N-1) ... hs(P-1) hs(P) hs(P+1) ... hs(P+N) where N = (int)windowHalfSize
            // the match on the key should be done with 7 digits precision, i.e (int)(P*10000000)
            
            // resampled sample value =
            // sum over every sampleLocation in the original signal of sampleValue * hs(where-sampleLocation)
            
            assert(vres.size() == 2+2*N);
            
#if DEBUG_RESAMPLING
            std::map<int64_t, double> tmp;
#endif
            for(int64_t
                i   = std::max(static_cast<int64_t>(0),
                               N-discreteSampleBefore),
                end = std::min(static_cast<int64_t>(vres.size()),
                               static_cast<int64_t>(original.size()/n_channels)+N-discreteSampleBefore);
                i < end;
                ++i) {
                // for all j in 0 .. n_channels-1,
                //
                // (discreteSampleBefore+i-N)*n_channels + j >= 0   implies:
                // discreteSampleBefore+i-N >= 0   implies:
                // i >= N-discreteSampleBefore
                //
                // (discreteSampleBefore+i-N)*n_channels + j < original.size() implies:
                // (discreteSampleBefore+i-N)*n_channels + n_channels-1 < original.size() implies:
                // (discreteSampleBefore+i-N+1)*n_channels -1 < original.size() implies:
                // (discreteSampleBefore+i-N+1)*n_channels < 1 + original.size() implies:
                // (because the left is a multiple of n_channels, and original.size() is a multiple of n_channels)
                // (discreteSampleBefore+i-N+1)*n_channels < n_channels + original.size() implies:
                // discreteSampleBefore+i-N+1 < 1 + original.size()/n_channels implies:
                // i < N + original.size()/n_channels - discreteSampleBefore
                
                auto const vresValue = vres[i];
                auto const idxBase = (discreteSampleBefore+i-N)*n_channels;
                for(int j=0; j<n_channels; ++j) {
                    auto originalIdx = idxBase + j;
#if DEBUG_RESAMPLING
                    tmp[originalIdx] = vresValue;
#endif
                    it[j] += vresValue * original[originalIdx];
                }
            }
#if DEBUG_RESAMPLING
            std::cout << "frame " << frameIdx << " where " << where << " value " << res << std::endl;
            for(auto const & [originalIdx, vresValue] : tmp) {
                std::cout << originalIdx << " " << vresValue << std::endl;
            }
#endif
        }
        std::cout << results.size() << std::endl;
        return dt;
    }
    
    template<typename Reader, typename Buffer>
    ResampleSincStats resampleSinc(Reader & reader, Buffer & resampled, double new_sample_rate)
    {
        using namespace std::chrono;
        using namespace profiling;
        ResampleSincStats stats;
        
        using VAL = typename Buffer::value_type;
        
        resampled.clear();

        std::vector<VAL> original;
        // TODO faire cela dans un thread a part et en meme temps resampler?
        {
            Timer<steady_clock> t(stats.dt_read_source);
            auto const nOrigSamples = reader.countSamples();
            original.reserve(nOrigSamples);
            while(reader.HasMore()) {
                original.push_back(reader.template ReadAsOneFloat<VAL>());
            }
            if(nOrigSamples != original.size()) {
                throw std::logic_error("invalid reader (number of samples != number os samples read : "
                                       + std::to_string(nOrigSamples) + " " + std::to_string(original.size()));
            }
        }
        double const Fs = reader.getSampleRate();
        double const FsPrime = new_sample_rate;
        if(!Fs || !FsPrime) {
            return stats;
        }
        
        double const Fs_over_FsPrime = Fs / FsPrime;

        stats.dt_resample = resampleSincBuffer(original,
                                               reader.countChannels(),
                                               Fs_over_FsPrime,
                                               resampled);
        return stats;
    }

    struct InterlacedBuffer {
        using element_type = double;
        
        template<typename Reader>
        InterlacedBuffer(Reader & reader, double sample_rate, ResampleSincStats & stats)
        : nchannels(reader.countChannels())
        {
            stats = resampleSinc(reader, buf, sample_rate);
        }
        /*
        InterlacedBuffer(std::vector<element_type> const & v, int nchannels)
        : buf(v)
        , nchannels(nchannels)
        {}*/
        
        std::size_t countFrames() const {
            if(nchannels == 0) {
                return 0;
            }
            return buf.size() / nchannels;
        }
        
        int countChannels() const { return nchannels; }
        
        std::vector<element_type> const & getBuffer() const { return buf; }
        
        element_type getMaxAbsValue() const {
            if(!maxAbsValue) {
                element_type m = 0;
                for(auto v:buf) {
                    m = std::max(m, std::abs(v));
                }
                maxAbsValue = m;
            }
            return *maxAbsValue;
        }
        
        
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
        int nchannels;
        mutable std::optional<element_type> maxAbsValue;
    };
    
} // namespace imajuscule::audio
