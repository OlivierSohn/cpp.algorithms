

namespace imajuscule {

using FftSpecs = std::map<uint32_t, int>; // sz -> historySize

template<typename T, typename Tag>
struct XAndFFTS {
    static constexpr auto zero_n_raw = fft::RealSignal_<Tag, T>::zero_n_raw;

    using Segment = std::pair<int, int>;

    std::pair<Segment, Segment> getSegments(int const h, int const pastSize) const {
        int start = progress - h - pastSize;
        while(start < 0) {
            start += x_unpadded_size;
        }
        Assert(start >= 0);
        
        int const diff = start+pastSize-x_unpadded_size;
        if(diff <= 0)
        {
            return {
                {start, pastSize},
                {0,0}
            };
        }
        else {
            return {
                {start, pastSize-diff},
                {0    , diff}
            };
        }
    }
    
    void resize(int const sz, FftSpecs const & fftSpecs) {
        reset_states();

        x_ffts.clear();
        x_ffts.reserve(fftSpecs.size());
        fftHalfSizesBits = 0;
        fftsHalfSizesBitsToCompute = 0;
        for(auto const & [fftSz, history_sz] : fftSpecs) {
            Assert(history_sz > 0);
            Assert(is_power_of_two(fftSz));
            fftHalfSizesBits |= (fftSz/2);

            x_ffts.emplace_back(fftSz,
                                history_sz);
        }
        Assert(count_set_bits(fftHalfSizesBits) == x_ffts.size());
        int const maxHalfFFTSz = floor_power_of_two(fftHalfSizesBits);

        x_unpadded_size = std::max(maxHalfFFTSz, sz);

        x.clear();
        x.resize(x_unpadded_size + maxHalfFFTSz); // add padding for biggest fft
    }
    
    void push(T v) {
        if(unlikely(progress == x_unpadded_size)) {
            progress = 0;
            padding = 0;
        }
        x[progress] = typename RealSignal::value_type(v);
        ++progress;
        padding = std::max(padding, progress);
        if(unlikely(padding == x_unpadded_size)) {
            padding = static_cast<int>(x.size());
        }
        
        uint32_t nZeroes = count_trailing_zeroes(progress);
        uint32_t mask = (1 << (nZeroes+1)) - 1;
        fftsHalfSizesBitsToCompute = fftHalfSizesBits & mask;
        for(int i=0, nComputedFFTs = count_set_bits(fftsHalfSizesBitsToCompute);
            i<nComputedFFTs;
            ++i)
        {
            auto & f = x_ffts[i];
            int halfFftSize = f.fft_length/2;
            int xStart = progress - halfFftSize;
            Assert(xStart >= 0);

            int const xEnd = progress + halfFftSize;
            int const countPadding = xEnd - padding;
            if(countPadding > 0) {
                zero_n_raw(&x[padding], countPadding);
                padding = xEnd;
            }
            auto & oldest_fft_of_delayed_x = *f.ffts.cycleEnd();
            Assert(f.fft_length == static_cast<int>(oldest_fft_of_delayed_x.size()));
            f.fft.forward(x.begin() + xStart, oldest_fft_of_delayed_x, f.fft_length);
            f.ffts.advance();
        }
    }
    
    using RealSignal = typename fft::RealSignal_<Tag, T>::type;
    
    RealSignal x;
    int progress = 0; // last + 1
    int padding = 0;
    int x_unpadded_size = 0;
    uint32_t fftHalfSizesBits = 0; // if we use the fft of sz 2^n, the nth bit is set
    uint32_t fftsHalfSizesBitsToCompute = 0;
    
    using FBins = typename fft::RealFBins_<Tag, T>::type;
    struct FFTs {
        using Contexts = fft::Contexts_<Tag, T>;
        using FFTAlgo = typename fft::Algo_<Tag, T>;

        FFTs(int fft_length, int history_sz)
        : fft(Contexts::getInstance().getBySize(fft_length))
        , fft_length(fft_length)
        , ffts(history_sz)
        {
            ffts.for_each([fft_length](auto & v){ v.resize(fft_length); });
        }
        
        FFTAlgo fft;
        
        int fft_length;
        cyclic<FBins> ffts;
    };
    
    std::vector<FFTs> x_ffts;
    
    auto const & find_ffts(int sz) const {
        for(auto const & f : x_ffts) {
            if(sz == f.fft_length) {
                return f;
            }
        }
        throw std::runtime_error("fft not found");
    }
private:
    
    void reset_states() {
        progress = 0; // last + 1
        padding = 0; // firstNonPadded
    }
};

template<typename T, typename Tag>
struct Y {
    void resize(int const blockSz,
                int const nAnticipatedWrites) {
        blockSizeMinusOne = blockSz-1;
        int const nAnticipationBlocks = nAnticipatedWrites ? (1 + (nAnticipatedWrites-1) / blockSz) : 0;
        Assert(nAnticipationBlocks >= 0);
        int const nBlocks = 1 + nAnticipationBlocks;

        y.clear();
        ySz = blockSz * nBlocks;
        y.resize(ySz); // intialized at 0

        uProgress = 0;
        zeroed_up_to = ySz;
    }
    
    void increment() {
        ++uProgress;
        Assert(uProgress <= ySz);
        // if a block boundary was crossed, zero the block starting at 'zeroed_up_to'
        Assert(is_power_of_two(blockSizeMinusOne+1));
        if(0 == (uProgress & blockSizeMinusOne)) {
            if(uProgress == ySz) {
                uProgress = 0;
            }
            Assert(zeroed_up_to <= ySz);
            unsigned int zeroStart = (zeroed_up_to == ySz) ? 0 : zeroed_up_to;
            auto blockSize = blockSizeMinusOne + 1;
            zeroed_up_to = zeroStart + blockSize;
            Assert(zeroStart < ySz);
            Assert(zeroed_up_to <= ySz);
            fft::RealSignal_<Tag, T>::zero_n_raw(&y[zeroStart], blockSize);
        }
        Assert(uProgress < ySz);
    }

    using RealSignal = typename fft::RealSignal_<Tag, T>::type;
    RealSignal y;
    uint32_t uProgress = 0;
    uint32_t ySz=0;
    uint32_t blockSizeMinusOne = 0;
    int32_t zeroed_up_to = 0;
};

struct MinSizeRequirement {
    MinSizeRequirement() = delete;
    
    /*
     The needed size of the x buffer is deduced both from the max size of the ffts, and from this parameter.
     */
    int minXSize;
        
    int minYSize;
    // "Anticipated" writes touch the block (of size minYSize) after the current block.
    int minYAnticipatedWrites;
    
    FftSpecs xFftSizes;
    
    void mergeWith(MinSizeRequirement const & o) {
        minXSize = std::max(minXSize, o.minXSize);
        
        {
            int res = static_cast<int>(ppcm(minYSize,
                                            o.minYSize));
            
            // only because in practice we have powers of 2.
            Assert(std::max(minYSize,
                            o.minYSize) == res);

            minYSize = res;
        }
        
        minYAnticipatedWrites = std::max(minYAnticipatedWrites,
                                         o.minYAnticipatedWrites);
        
        for(auto const & [sz, historySize] : o.xFftSizes) {
            auto [it, emplaced] = xFftSizes.try_emplace(sz, historySize);
            if(!emplaced) {
                it->second = std::max(it->second, historySize);
            }
        }
    }
};

template<typename T, typename Tag>
struct DspContext {
    XAndFFTS<T, Tag> x_and_ffts;
    Y<T, Tag> y;
    
    void resize(MinSizeRequirement req) {
        x_and_ffts.resize(req.minXSize, req.xFftSizes);
        y.resize(req.minYSize, req.minYAnticipatedWrites);
    }
};

/* Used to simplify tests */
template<typename Algo, typename Tag = typename Algo::Tag>
struct Convolution {
    using State = typename Algo::State;
    using FPT = typename Algo::FPT;
    using DspContext = DspContext<FPT, Tag>;
    using SetupParam = typename Algo::SetupParam;

    void setup(SetupParam const & p) {
        algo.setup(p);
    }
    void setCoefficients(a64::vector<FPT> v) {
        Assert(v.size());
        MinSizeRequirement req = state.setCoefficients(algo, std::move(v));
        
        ctxt.resize(req);
        
        int const ySz = static_cast<int>(ctxt.y.y.size());
        int const nSteps = (ySz >= 1) ? 1 : 0;

        // set y progress such that results are written in a single chunk
        for(int i=0; i<nSteps; ++i) {
            ctxt.y.increment();
        }
    }
    
    bool isValid() const {
        return algo.isValid();
    }
    double getEpsilon() const {
        return state.getEpsilon(algo);
    }
    int getLatency() const {
        return algo.getLatency();
    }
    
    FPT step(FPT val) {
        ctxt.x_and_ffts.push(val);
        
        algo.step(state,
                  ctxt.x_and_ffts,
                  ctxt.y);
        
        FPT res = fft::RealSignal_<Tag, FPT>::get_signal(ctxt.y.y[ctxt.y.uProgress]);
        ctxt.y.increment();
        return res;
    }
    
    auto & getAlgo() {
        return algo;
    }
    auto & getAlgo() const {
        return algo;
    }
private:
    DspContext ctxt;
    State state;
    Algo algo;
};

}
