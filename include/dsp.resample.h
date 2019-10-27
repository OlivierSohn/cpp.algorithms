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
    
    template<typename T>
    T fSinc(T x) {
        if(likely(x)) {
            return std::sin(x)/x;
        }
        return 1.;
    }

    // We will sample sinc so as to minimize the absolute error when interpolating between 2 consecutive samples.
    // To minimize
    // We want to find inflexion point of sinc, i.e where the 2nd derivative changes sign.
    // 2nd derivative of sinc is: -(2 x cos(x) + (-2 + x^2) sin(x))/x^3
    // so we need to find zeroes of the numerator : 2 x cos(x) + (-2 + x^2) sin(x)
    // so we need the derivative of the numerator (for newton root finding method) which is : x^2 cos(x)

    template<typename T>
    T der_of_der_of_NumeratorSinc(T x) {
        return (2. * x * std::cos(x)) + ((-2. + x * x) * std::sin(x));
    }

    template<typename T>
    T der_of_der_of_der_of_NumeratorSinc(T x) {
        return x * x * std::cos(x);
    }
    
    template<typename F, typename DER, typename T = decltype(std::declval<F>()({}))>
    std::vector<T> findRootsSpacedByAtleast(T minDistanceBetweenConsecutveRoots, T const from, T const to, F f, DER der) {
        std::vector<T> res;
        if(from >= to) {
            return res;
        }
        res.reserve(static_cast<int>((to-from)*minDistanceBetweenConsecutveRoots));
        
        auto constexpr epsilon = std::numeric_limits<T>::epsilon();
        
        for(T x1 = from; x1 < to; x1+=minDistanceBetweenConsecutveRoots) {
            T x2 = x1 + minDistanceBetweenConsecutveRoots;
            if(x2 > to) {
                x2 = to;
            }
            auto y1 = f(x1);
            if(std::abs(y1) < epsilon) {
                res.push_back(x1);
                continue;
            }
            auto y2 = f(x2);
            if(std::abs(y2) < epsilon) {
                res.push_back(x2);
                continue;
            }
            
            if(y1*y2 > 0.) {
                // same signe, hence there is no root in this interval
                continue;
            }
            auto mayRoot = find_root(f, der, x1, x2);
            if(mayRoot) {
                res.push_back(*mayRoot);
            }
        }
        return res;
    }
    
    template<typename T>
    auto findSincCurvatureChanges(T from, T to) {
        auto res =  findRootsSpacedByAtleast(2.,
                                        from,
                                        to,
                                        der_of_der_of_NumeratorSinc<double>,
                                        der_of_der_of_der_of_NumeratorSinc<double>);
        auto constexpr epsilon = std::numeric_limits<T>::epsilon();
        if(res.empty() || std::abs(*res.begin() - from) > epsilon) {
            res.insert(res.begin(), from);
        }
        if(res.empty() || std::abs(*res.rbegin() - to) > epsilon) {
            res.push_back(to);
        }
        return res;
    }
    
    
    template<typename T>
    struct UniformSamplerOpt {
        template<typename F>
        UniformSamplerOpt(F f, T from, T to, int numSamples)
        : from(from)
        , to(to)
        {
            numSamples = std::max(2, numSamples);
            
            increment = (to-from)/(numSamples-1);
            if(increment) {
                inv_increment = 1. / increment;
            }
            else {
                inv_increment = 0.;
            }
            
            samples.reserve(numSamples);
            for(int i=0; i<numSamples; ++i) {
                auto y = f(from + i*increment);
                samples.push_back({y,{}});
            }
            for(int i=0; i<numSamples-1; ++i) {
                samples[i].second = samples[i+1].first - samples[i].first;
            }
            samples[numSamples-1].second = 0;
        }
        
        T getAt(T x) const {
            if(x <= from) {
                return samples.begin()->first;
            }
            if(x >= to) {
                return samples.rbegin()->first;
            }
            T distFromBegin = (x-from)*inv_increment;
            int index = static_cast<int>(distFromBegin);
            T lambda = distFromBegin-index;
            //return (1.-lambda) * samples[index] + lambda * samples[index+1];
            auto const & s = samples[index];
            return s.first + lambda * s.second;
        }
    private:
        T from, to, increment, inv_increment;
        std::vector<std::pair<T,T>> samples;
    };
    
    template<typename T>
    struct UniformSampler {
        template<typename F>
        UniformSampler(F f, T from, T to, int numSamples)
        : from(from)
        , to(to)
        {
            numSamples = std::max(2, numSamples);
            
            increment = (to-from)/(numSamples-1);
            if(increment) {
                inv_increment = 1. / increment;
            }
            else {
                inv_increment = 0.;
            }
            
            samples.reserve(numSamples);
            for(int i=0; i<numSamples; ++i) {
                auto y = f(from + i*increment);
                samples.push_back(y);
            }
        }
        
        T getAt(T x) const {
            if(x <= from) {
                return *samples.begin();
            }
            if(x >= to) {
                return *samples.rbegin();
            }
            T distFromBegin = (x-from)*inv_increment;
            int index = static_cast<int>(distFromBegin);
            T lambda = distFromBegin-index;
            return (1.-lambda) * samples[index] + lambda * samples[index+1];
        }
    private:
        T from, to, increment, inv_increment;
        std::vector<T> samples;
    };

    template<typename T>
    struct NonUniformSampler {
        struct SortableXY {
            SortableXY(T x, T y)
            : x(x)
            , y(y)
            {}
            
            T x;
            T y;
            
            bool operator < (SortableXY const & other) const {
                return x < other.x;
            }
        };

        using Container = std::vector<SortableXY>;
        using ConstIterator = typename Container::const_iterator;
        
        template<typename F>
        NonUniformSampler(F functionToSample, std::vector<T> const & curvatureChanges, T maxSamplingError)
        : maxSamplingError(maxSamplingError)
        {
            sorted.reserve(100000);
            
            static_assert(std::is_same_v<T, decltype(functionToSample({}))>);
            auto it = curvatureChanges.begin();
            auto const end = curvatureChanges.end();
            if(it == end) {
                throw std::logic_error("empty curvatureChanges");
            }
            T v1 = functionToSample(*it);
            sorted.emplace_back(*it, v1);
            for(auto it2 = it+1; it2 != end; ++it2) {
                T v2 = functionToSample(*it2);
                sampleSameDerDerSign(functionToSample, *it, *it2, v1, v2);
                sorted.emplace_back(*it2, v2);
                it = it2;
                v1 = v2;
            }
        }
        T getAt(T x) const {
            ConstIterator dummy;
            return getAtWithMin(x, sorted.begin(), dummy);
        }
        T getAtWithMin(T x, ConstIterator minIt, ConstIterator & lessIt) const {
            // first, find lower+upper bounds by doubling the interval each time.
            auto const end = sorted.end();
            ConstIterator nextMinIt = end;

            int sz = 1;
            while(true)
            {
                nextMinIt = minIt + sz;
                if(nextMinIt >= end) {
                    nextMinIt = end;
                    break;
                }
                if(nextMinIt->x > x) {
                    break;
                }
                minIt = nextMinIt;
                sz *= 2;
            }
            return getAtWithMinEnd(x, minIt, nextMinIt, lessIt);
        }
        T getAtWithMinEnd(T x, ConstIterator minIt, ConstIterator endIt, ConstIterator & lessIt) const {
            lessIt = std::upper_bound(minIt,
                                      endIt,
                                      x,
                                      [](auto const & x, auto const & other) { return x < other.x; });
            if(lessIt == minIt) {
                return lessIt->y;
            }
            if(lessIt == endIt && endIt == sorted.end()) {
                return sorted.rbegin()->y;
            }
            auto l2 = lessIt;
            --lessIt;
            auto l1 = lessIt;
            auto exactDX = l2->x - l1->x;
            if(!exactDX) {
                return l2->y;
            }
            auto ratio = (x-l1->x)/exactDX;
            return (1.-ratio) * l1->y + ratio * l2->y;
        }
        
        auto const & getExactValues() const { return sorted; }
    private:
        std::vector<SortableXY> sorted;
        T maxSamplingError;

        // betweem from and to, the 2nd derivative of function to sample has a constant sign.
        template<typename F>
        void sampleSameDerDerSign(F const & functionToSample, T from, T to, T vFrom, T vTo) {
            // by design, the samples at from and to are already sampled.
            
            T mid = 0.5 * (to+from);
            if(unlikely(mid == to || mid == from)) {
                return;
            }
            T vMid = functionToSample(mid);
            T vInterp = 0.5 * (vFrom + vTo);
            T absoluteError = std::abs(vMid - vInterp);
            if(absoluteError <= maxSamplingError) {
                // we could put the sample in the map, but we don't need to since the sampling error criteria is met.
                return;
            }
            sampleSameDerDerSign(functionToSample, from, mid, vFrom, vMid);
            sorted.emplace_back(mid, vMid);
            sampleSameDerDerSign(functionToSample, mid, to, vMid, vTo);
        }
    };
    
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

    
    template<typename OriginalBuffer, typename ResampledBuffer, typename F>
    std::chrono::steady_clock::duration resampleSincBufferVariableRate(OriginalBuffer const & original,
                                                                       int const n_channels,
                                                                       F fFs_over_FsPrime,
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
        int64_t const n_orig_frames = original.size() / n_channels;
        
        auto foreachResampled = [fFs_over_FsPrime,
                                 n_orig_frames](auto f)
        {
            constexpr double epsilon = std::numeric_limits<double>::epsilon();

            double where = 0.;
            for(int64_t frameIdx = 0;
                ;
                ++frameIdx)
            {
                // we are computing frame 'frameIdx' of the resampled signal.
                
                // 'where' is the corresponding location in the original signal.
                
                // in the original signal, 'where' is between integers 'discreteSampleBefore' and 'discreteSampleBefore'+1
                int64_t const discreteSampleBefore = static_cast<int64_t>(where);
                
                if(unlikely(discreteSampleBefore >= n_orig_frames)) {
                    return;
                }
                
                double const P = where-discreteSampleBefore; // in [0,1)
                
                if(discreteSampleBefore == n_orig_frames-1 && P > frameIdx*epsilon) {
                    return;
                }
                auto Fs_over_FsPrime = fFs_over_FsPrime(where);
                f(frameIdx, Fs_over_FsPrime, discreteSampleBefore, P);
                where += Fs_over_FsPrime;
            }
        };

        int64_t target_frame_count=0;
        foreachResampled([&target_frame_count,
                          n_orig_frames
                          ](int64_t frameIdx,
                            double Fs_over_FsPrime,
                            int64_t discreteSampleBefore,
                            double P) {
            target_frame_count = frameIdx+1;
        });
        resampled.resize(target_frame_count*n_channels, {});

        KaiserWindow window;
        constexpr auto num_zerocrossings_half_window = 512;
        constexpr auto samplingPrecision = 1e-10;
        NonUniformSampler sampledSinc(fSinc<double>,
                                      findSincCurvatureChanges(0., M_PI * num_zerocrossings_half_window),
                                      samplingPrecision);

        std::vector<double> vres;

        foreachResampled([&vres,
                          &sampledSinc,
                          &resampled,
                          &original,
                          &window,
                          num_zerocrossings_half_window,
                          n_orig_frames,
                          n_channels](int64_t frameIdx,
                                      double Fs_over_FsPrime,
                                      int64_t discreteSampleBefore,
                                      double P) {
            // we are computing frame 'frameIdx' of the resampled signal.
            
            // when the resampling frequency is constant, we can cache the hs results array with key 'P'
            // this array contains results for hs(P-N-1) ... hs(P-1) hs(P) hs(P+1) ... hs(P+N) where N = (int)windowHalfSize
            // the match on the key should be done with 7 digits precision, i.e (int)(P*10000000)
            
            // resampled sample value =
            // sum over every sampleLocation in the original signal of sampleValue * hs(where-sampleLocation)
            

            
#if DEBUG_RESAMPLING
            std::map<int64_t, double> tmp;
#endif
            
            const auto one_over_zeroCrossingDistance = std::min(1., 1./Fs_over_FsPrime);
            // half size, excluding the central point, hence the number of non-zero slots is:
            // 1 + 2*windowHalfSize
            constexpr double one_over_num_zerocrossings_half_window = 1./num_zerocrossings_half_window;
            const auto windowHalfSize = num_zerocrossings_half_window / one_over_zeroCrossingDistance;
            int64_t const N = static_cast<int>(windowHalfSize) + 1;

            auto hs = [num_zerocrossings_half_window,
                       one_over_num_zerocrossings_half_window,
                       one_over_zeroCrossingDistance,
                       &window,
                       &sampledSinc]
            (double t) { // zero crossings every one_over_zeroCrossingDistance
                assert(t >= 0.);
                double x = t * one_over_zeroCrossingDistance; // zero crossings every x
                assert(x >= 0.);
                if(x >= num_zerocrossings_half_window) {
                    return 0.;
                }
                auto winCoeff = window.getAt(x * one_over_num_zerocrossings_half_window);
                return winCoeff * one_over_zeroCrossingDistance * sampledSinc.getAt(M_PI * x); // 0 <= M_PI * x < M_PI * num_zerocrossings_half_window
            };
            
            int64_t const vres_sz = 2+2*N;

            // P between 0 and 1
            {
                assert(P >= 0.);
                assert(P <= 1.);
                vres.clear();
                vres.reserve(vres_sz);
                int64_t firstIdx = -N;
                int64_t last_idx = N+1; // included
                for(int64_t i=firstIdx; i<= last_idx; ++i) {
                    vres.push_back(hs(std::abs(P-i)));
                }
                assert(vres.size() == vres_sz);
            };
            
            int64_t const discreteSampleBefore_minus_N = discreteSampleBefore-N;

            auto const resampleIdxBase = frameIdx*n_channels;

            for(int64_t
                i   = std::max(static_cast<int64_t>(0),
                               -discreteSampleBefore_minus_N),
                end = std::min(vres_sz,
                               n_orig_frames-discreteSampleBefore_minus_N);
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
                auto const idxBase = (discreteSampleBefore_minus_N + i)*n_channels;
                for(int j=0; j<n_channels; ++j) {
                    auto const originalIdx = idxBase + j;
#if DEBUG_RESAMPLING
                    tmp[originalIdx] = vresValue;
#endif
                    auto const resampleIdx = resampleIdxBase+j;
                    assert(resampleIdx < resampled.size());
                    resampled[resampleIdx] += vresValue * original[originalIdx];
                }
            }
#if DEBUG_RESAMPLING
            std::cout << "frame " << frameIdx << " where " << discreteSampleBefore + P << " value " << resampled[resampleIdxBase] << std::endl;
            for(auto const & [originalIdx, vresValue] : tmp) {
                std::cout << originalIdx << " " << vresValue << std::endl;
            }
#endif
        });
        return dt;
    }
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
        int64_t const n_orig_frames = original.size() / n_channels;

        if(!Fs_over_FsPrime) {
            return dt;
        }
        int64_t const target_frame_count = static_cast<int>(1 + (n_orig_frames-1) / Fs_over_FsPrime);
        if(target_frame_count <= 0) {
            return dt;
        }
        if(std::abs(Fs_over_FsPrime - 1.) < 1e-8) {
            resampled = original;
            return dt;
        }
        
        // using notations found in https://ccrma.stanford.edu/~jos/resample/resample.pdf
        const auto one_over_zeroCrossingDistance = std::min(1., 1./Fs_over_FsPrime);
        // half size, excluding the central point, hence the number of non-zero slots is:
        // 1 + 2*windowHalfSize
        constexpr auto num_zerocrossings_half_window = 512;
        constexpr double one_over_num_zerocrossings_half_window = 1./num_zerocrossings_half_window;
        const auto windowHalfSize = num_zerocrossings_half_window / one_over_zeroCrossingDistance;
        int64_t const N = static_cast<int>(windowHalfSize) + 1;
        KaiserWindow window;
        
        auto hs = [num_zerocrossings_half_window,
                   one_over_num_zerocrossings_half_window,
                   one_over_zeroCrossingDistance,
                   &window]
        (double t) { // zero crossings every one_over_zeroCrossingDistance
            assert(t >= 0.);
            auto sinc = [num_zerocrossings_half_window,
                         &window](double x) { // zero crossings every x
                assert(x >= 0.);
                bool inWindow = (x > -num_zerocrossings_half_window) && (x < num_zerocrossings_half_window);
                if(!inWindow) {
                    return 0.;
                }
                auto winCoeff = window.getAt(std::abs(x) * one_over_num_zerocrossings_half_window);
                return winCoeff * fSinc(M_PI * x);
            };
            return one_over_zeroCrossingDistance * sinc(t * one_over_zeroCrossingDistance);
        };
        
        auto foreachResampled = [target_frame_count, Fs_over_FsPrime](auto f){
            for(int64_t frameIdx = 0; frameIdx < target_frame_count; ++frameIdx) {
                // we are computing frame 'frameIdx' of the resampled signal.
                
                // 'where' is the corresponding location in the original signal.
                double const where = Fs_over_FsPrime*frameIdx; // >= 0
                
                // in the original signal, 'where' is between integers 'discreteSampleBefore' and 'discreteSampleBefore'+1
                int64_t const discreteSampleBefore = static_cast<int64_t>(where);
                
                double const P = where-discreteSampleBefore; // in [0,1)

                if(!f(frameIdx, discreteSampleBefore, P)) {
                    return;
                }
            }
        };
        
        // This assumes the frame rate ratio is constant!
        int64_t nDifferentLocations = 1;
        foreachResampled([&nDifferentLocations](int64_t frameIdx, int64_t, double P) {
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
                v.push_back(hs(std::abs(x-i)));
            }
            assert(v.size() == 2+2*N);
            results.insert({key, std::move(res)});
            return v;
        };

        resampled.resize(target_frame_count*n_channels, {});
        foreachResampled([&resampled,
                          &original,
                          N,
                          n_orig_frames,
                          getResults,
                          n_channels](int64_t frameIdx,
                                      int64_t discreteSampleBefore,
                                      double P) {
            // we are computing frame 'frameIdx' of the resampled signal.

            auto const & vres = getResults(P);
            // when the resampling frequency is constant, we can cache the hs results array with key 'P'
            // this array contains results for hs(P-N-1) ... hs(P-1) hs(P) hs(P+1) ... hs(P+N) where N = (int)windowHalfSize
            // the match on the key should be done with 7 digits precision, i.e (int)(P*10000000)
            
            // resampled sample value =
            // sum over every sampleLocation in the original signal of sampleValue * hs(where-sampleLocation)
            
            
#if DEBUG_RESAMPLING
            std::map<int64_t, double> tmp;
#endif
            int64_t const discreteSampleBefore_minus_N = discreteSampleBefore-N;
            int64_t const vres_sz = 2+2*N;
            assert(vres_sz == vres.size());
            
            auto const resampleIdxBase = frameIdx*n_channels;
            
            for(int64_t
                i   = std::max(static_cast<int64_t>(0),
                               -discreteSampleBefore_minus_N),
                end = std::min(vres_sz,
                               n_orig_frames-discreteSampleBefore_minus_N);
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
                auto const idxBase = (discreteSampleBefore_minus_N + i)*n_channels;
                for(int j=0; j<n_channels; ++j) {
                    auto const originalIdx = idxBase + j;
#if DEBUG_RESAMPLING
                    tmp[originalIdx] = vresValue;
#endif
                    resampled[resampleIdxBase+j] += vresValue * original[originalIdx];
                }
            }
#if DEBUG_RESAMPLING
            std::cout << "frame " << frameIdx << " where " << discreteSampleBefore + P << " value " << resampled[resampleIdxBase] << std::endl;
            for(auto const & [originalIdx, vresValue] : tmp) {
                std::cout << originalIdx << " " << vresValue << std::endl;
            }
#endif
            return true;
        });
        return dt;
    }
    
    template<typename Reader, typename Buffer, typename F>
    ResampleSincStats resampleSincVariableRate(Reader & reader, Buffer & resampled, F f_new_sample_rate)
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
        if(!Fs) {
            return stats;
        }
        
        auto fFs_over_FsPrime = [f_new_sample_rate, Fs](double where) {
            auto const FsPrime = f_new_sample_rate(where);
            if(!FsPrime) {
                assert(0);
                return 1.;
            }
            return Fs / FsPrime;
        };

        stats.dt_resample = resampleSincBufferVariableRate(original,
                                                           reader.countChannels(),
                                                           fFs_over_FsPrime,
                                                           resampled);
        return stats;
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

        template<typename Reader, typename S>
        InterlacedBuffer(Reader & reader, S sample_rate, ResampleSincStats & stats)
        : nchannels(reader.countChannels())
        {
            if constexpr (std::is_pod<S>::value) {
                stats = resampleSinc(reader, buf, sample_rate);
            }
            else {
                stats = resampleSincVariableRate(reader, buf, sample_rate);
            }
        }

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
