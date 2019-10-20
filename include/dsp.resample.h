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
        double const frameRatio = reader.getSampleRate() / sample_rate;
        if(!frameRatio) {
            return;
        }
        int target_size = static_cast<int>(1 + (reader.countFrames()-1) / frameRatio) * n_channels;
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
            where = frameRatio*frameIdx;
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

    template<typename Reader, typename Buffer>
    void resampleSinc(Reader & reader, Buffer & resampled, double new_sample_rate)
    {
        using VAL = typename Buffer::value_type;
        
        resampled.clear();
        auto const n_channels = reader.countChannels();
        auto const original_sample_rate = reader.getSampleRate();
        auto const n_orig_frames = reader.countFrames();
        double const frameRatio = original_sample_rate / new_sample_rate;
        if(!frameRatio) {
            return;
        }
        int target_size = static_cast<int>(1 + (n_orig_frames-1) / frameRatio) * n_channels;
        if(target_size <= 0) {
            return;
        }
        resampled.resize(target_size);
        auto it = resampled.begin();
        auto end = resampled.end();

        std::vector<VAL> original;
        {
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
                
        // using notations found in https://ccrma.stanford.edu/~jos/resample/resample.pdf
        double const FsPrime = new_sample_rate;
        double const Fs = original_sample_rate;
        double const norm = std::min(1., FsPrime/Fs);
        double const minFsFsPrime = std::min(Fs, FsPrime);
        const auto zeroCrossingDistance = Fs/minFsFsPrime;
        // half size, excluding the central point, hence the number of non-zero slots is:
        // 1 + 2*windowHalfSize
        constexpr auto num_zerocrossings_half_window = 512;//10.;
        const auto windowHalfSize = num_zerocrossings_half_window * zeroCrossingDistance;
        int64_t const N = static_cast<int>(windowHalfSize);
        KaiserWindow window;
        
        /*
        auto hs = [norm, minFsFsPrime, Fs, windowHalfSize, &window] (double t, bool & inWindow) { // zero crossings every Fs/minFsFsPrime
            auto sinc = [windowHalfSize, &window](double x, bool & inWindow) { // zero crossings every x
                if(!x) {
                    inWindow = true;
                    return 1.;
                }
                inWindow = (x > -windowHalfSize) && (x < windowHalfSize);
                if(!inWindow) {
                    return 0.;
                }
                auto winCoeff = window.getAt(std::abs(x)/windowHalfSize);
                double pi_x = M_PI * x;
                return winCoeff * std::sin(pi_x)/pi_x;
            };
            return norm * sinc(t * minFsFsPrime / Fs, inWindow);
        };*/
        
        
        auto hs = [norm, minFsFsPrime, Fs, windowHalfSize, &window] (double t) { // zero crossings every Fs/minFsFsPrime
            auto sinc = [windowHalfSize, &window](double x) { // zero crossings every x
                bool inWindow = (x > -windowHalfSize) && (x < windowHalfSize);
                if(!inWindow) {
                    return 0.;
                }
                auto winCoeff = window.getAt(std::abs(x)/windowHalfSize);
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
            return norm * sinc(t * minFsFsPrime / Fs);
        };
        
        // we memoize the results
        std::unordered_map<AlmostDouble, std::unique_ptr<std::vector<double>>> results;
        results.reserve(300);
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
            int64_t firstIdx = -N-1;
            int64_t last_idx = N; // included
            for(int64_t i=firstIdx; i<= last_idx; ++i) {
                v.push_back(hs(x+i));
            }
            results.insert({key, std::move(res)});
            return v;
        };
        
        for(int64_t frameIdx = 0; it != end; ++frameIdx) {
            // we are computing frame 'frameIdx' of the resampled signal.

            // 'where' is the corresponding location in the original signal.
            double const where = frameRatio*frameIdx; // >= 0
            
            // in the original signal, 'where' is between integers 'discreteSampleBefore' and 'discreteSampleBefore'+1
            int64_t const discreteSampleBefore = static_cast<int64_t>(where);

            double const P = where-discreteSampleBefore; // in [0,1)
            auto const & vres = getResults(P);
            // when the resampling frequency is constant, we can cache the hs results array with key 'P'
            // this array contains results for hs(P-N-1) ... hs(P-1) hs(P) hs(P+1) ... hs(P+N) where N = (int)windowHalfSize
            // the match on the key should be done with 7 digits precision, i.e (int)(P*10000000)

            for(int j=0; j<n_channels; ++j, ++it) {// TODO make it inner loop
                
                // resampled sample value =
                // sum over every sampleLocation in the original signal of sampleValue * hs(where-sampleLocation)
                double res = 0.;

                assert(vres.size() == 2+2*N);
                for(int64_t
                    i   = std::max(static_cast<int64_t>(0),
                                   N+1-discreteSampleBefore),
                    end = std::min(static_cast<int64_t>(vres.size()),
                                   static_cast<int64_t>(original.size()/n_channels)+N+1-discreteSampleBefore);
                    i < end;
                    ++i) {
                    // for all j in 0 .. n_channels-1,
                    //
                    // (discreteSampleBefore+i-N-1)*n_channels + j >= 0   implies:
                    // discreteSampleBefore+i-N-1 >= 0   implies:
                    // i >= N+1-discreteSampleBefore
                    //
                    // (discreteSampleBefore+i-N-1)*n_channels + j < original.size() implies:
                    // (discreteSampleBefore+i-N-1)*n_channels + n_channels-1 < original.size() implies:
                    // (discreteSampleBefore+i-N)*n_channels -1 < original.size() implies:
                    // (discreteSampleBefore+i-N)*n_channels < 1 + original.size() implies:
                    // (because the left is a multiple of n_channels, and original.size() is a multiple of n_channels)
                    // (discreteSampleBefore+i-N)*n_channels < n_channels + original.size() implies:
                    // discreteSampleBefore+i-N < 1 + original.size()/n_channels implies:
                    // i < 1 + N + original.size()/n_channels - discreteSampleBefore

                    res += vres[i] * original[(discreteSampleBefore+i-N-1)*n_channels + j];
                }

                *it = res;
            }
        }
        std::cout << "num different keys = " << results.size() << std::endl;
    }
    
    enum class ResamplingMethod {
        LinearInterpolation,
        SincInterpolation
    };

    struct InterlacedBuffer {
        using element_type = double;
        
        template<typename Reader>
        InterlacedBuffer(Reader & reader, double sample_rate, ResamplingMethod resample)
        : nchannels(reader.countChannels())
        {
            if(resample == ResamplingMethod::LinearInterpolation) {
                resampleLinear(reader, buf, sample_rate);
            }
            else if(resample == ResamplingMethod::SincInterpolation) {
                resampleSinc(reader, buf, sample_rate);
            }
            else {
                throw std::logic_error("unhandled ResamplingMethod");
            }
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
