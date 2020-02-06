

namespace imajuscule {

/*
 
 Contracts:
 
 ** algo **
 
 - setup()
 
 Must be called before state.setCoefficients().
 
 ** state **
 
 - setCoefficients()
 
 puts the convolution in a state such that algo.step() has an effect

 - reset()
 
 puts the convolution in a state such that algo.step() has no effect


 */


template<typename T, template<typename> typename Allocator, typename Tag>
struct FFTs {
    using Contexts = fft::Contexts_<Tag, T>;
    using FFTAlgo = typename fft::Algo_<Tag, T>;
    using FBins = typename fft::RealFBins_<Tag, T, Allocator>::type;
    static constexpr auto zero = fft::RealFBins_<Tag, T, Allocator>::zero;

    FFTs(int fft_length, int history_sz)
    : fft(Contexts::getInstance().getBySize(fft_length))
    , fft_length(fft_length)
    , ffts(history_sz)
    {
        ffts.for_each([fft_length](auto & v){ v.resize(fft_length); });
    }
    
    FFTAlgo fft;
    
    int fft_length;
    
    auto & go_back() {
        ffts.go_back();
        return *ffts.cycleEnd();
    }
    auto * get_by_age(int age) const {
        return ffts.get_forward(age).data();
    }
    
    // This method will need to be replaced by rawindex-based iteration to support the async worker case.
    template<typename F>
    void for_some_recent(int count, F f) const {
        ffts.for_some_fwd(count, f);
    }
    
    void flushToSilence() {
        for(auto & fft : ffts) {
            zero(fft);
        }
    }
    
    int size() const {
        return ffts.size();
    }
private:
    cyclic<FBins> ffts;
};


template<typename T, template<typename> typename Allocator, typename Tag, typename WorkCplxFreqs>
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

    static int getAllocationSz_Resize(FftSpecs const & fftSpecs) {
        int sz = 0;
        for(auto const & [fftSz, history_sz] : fftSpecs) {
            if(history_sz == 0 || fftSz == 0) {
                continue;
            }
            Assert(is_power_of_two(fftSz));
            sz += fftSz * history_sz;
        }
        return sz;
    }

    void resize(int const sz, FftSpecs const & fftSpecs, int const workSz) {
        reset_states();

        x_ffts.clear();
        x_ffts.reserve(fftSpecs.size());
        fftHalfSizesBits = 0;
        fftsHalfSizesBitsToCompute = 0;
        for(auto const & [fftSz, history_sz] : fftSpecs) {
            if(history_sz == 0 || fftSz == 0) {
                continue;
            }
            Assert(is_power_of_two(fftSz));
            fftHalfSizesBits |= (fftSz/2);

            x_ffts.emplace_back(fftSz,
                                history_sz);
        }
        Assert(count_set_bits(fftHalfSizesBits) == x_ffts.size());
        int const maxHalfFFTSz = floor_power_of_two(fftHalfSizesBits);

        x_unpadded_size = std::max(maxHalfFFTSz,
                                   static_cast<int>(ceil_power_of_two(std::max(1,
                                                                               sz))));

        x.clear();
        x.resize(x_unpadded_size + maxHalfFFTSz); // add padding for biggest fft
        padding = x.size();
        
        work.reserve(workSz);
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
            auto & oldest_fft_of_delayed_x = f.go_back();
            Assert(f.fft_length == static_cast<int>(oldest_fft_of_delayed_x.size()));
            // we could omit that fft computation if:
            //   - there are only 0s in the input
            //   - all the ffts in the ring buffer are only 0s.
            f.fft.forward(x.begin() + xStart, oldest_fft_of_delayed_x.data(), f.fft_length);
        }
    }
    
    int flushToSilence() {
        for(auto & fft:x_ffts) {
            fft.flushToSilence();
        }
        if(!x.empty()) {
            zero_n_raw(&x[0], x.size());
        }
        int const progressBackup = progress;
        progress = 0;
        padding = x.size();
        fftsHalfSizesBitsToCompute = 0;
        Assert(progressBackup <= x_unpadded_size);
        return progressBackup ? (x_unpadded_size-progressBackup) : 0;
    }
    
    using RealSignal = typename fft::RealSignal_<Tag, T>::type;
    RealSignal x;
    
    int32_t progress = 0;
    int32_t padding = 0;
    int32_t x_unpadded_size = 0;
    uint32_t fftHalfSizesBits = 0; // if we use the fft of sz 2^(n+1), the nth bit is set
    uint32_t fftsHalfSizesBitsToCompute = 0;
    
    std::vector<FFTs<T, Allocator, Tag>> x_ffts;

    WorkCplxFreqs & work; // temporary work buffer
    
    auto const & find_ffts(int sz) const {
        for(auto const & f : x_ffts) {
            if(sz == f.fft_length) {
                return f;
            }
        }
        throw std::runtime_error("fft not found");
    }

    void logComputeState(std::ostream & os) const {
        os << "x size:" << x.size() <<  " unpadded:" << x_unpadded_size << " progress:" << progress << std::endl;
        
        os << "ffts: ";
        for(auto const & f:x_ffts) {
            os << f.size() << "x" << f.fft_length << " ";
        }
        os << std::endl;
    }
    
    XAndFFTS(WorkCplxFreqs & work)
    : work(work)
    {}
    
private:
    
    void reset_states() {
        progress = 0;
        padding = 0;
    }
};

template<typename T, typename Tag>
struct Y {
    static constexpr auto zero_n_raw = fft::RealSignal_<Tag, T>::zero_n_raw;
    
    void resize(int const blockSz_,
                int const nAnticipatedWrites) {
        int const blockSz = std::max(1, blockSz_);
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
            zero_n_raw(&y[zeroStart], blockSize);
        }
        Assert(uProgress < ySz);
    }
    
    void flushToSilence(int const nStepsForward) {
        zero_n_raw(&y[0], ySz);
        
        Assert(nStepsForward >= 0);
        
        for(int i=0; i<nStepsForward; ++i) {
            increment();
        }
    }

    void logComputeState(std::ostream & os) const {
        os << "y size " << ySz << std::endl;
        os << "y progress " << uProgress << std::endl;
    }
    
    using RealSignal = typename fft::RealSignal_<Tag, T>::type;
    RealSignal y;
    uint32_t uProgress = 0;
    uint32_t ySz=0;
    uint32_t blockSizeMinusOne = 0;
    int32_t zeroed_up_to = 0;
};

/* Used to simplify tests */
template<typename A>
struct XYConvolution {
    using Algo = A;

    static constexpr bool step_can_error = Algo::Desc::step_can_error;
    static constexpr bool has_subsampling = Algo::Desc::has_subsampling;
    
    using State = typename Algo::State;
    using FPT = typename Algo::FPT;
    using Tag = typename Algo::Tag;
    
    template<typename TT>
    using Allocator = typename A::template Allocator<TT>;
    
    using WorkCplxFreqs = typename fft::RealFBins_<Tag, FPT, aP::Alloc>::type;
    using SetupParam = typename Algo::SetupParam;

    XYConvolution(WorkCplxFreqs & work)
    : x_and_ffts(work)
    {}

    static int getAllocationSz_Setup(SetupParam const & p) {
        MinSizeRequirement req = p.getMinSizeRequirement();
        return XAndFFTS<FPT, Allocator, Tag, WorkCplxFreqs>::getAllocationSz_Resize(req.xFftSizes);
    }
    
    void setup(SetupParam const & p) {
        MinSizeRequirement req = p.getMinSizeRequirement();
        
        x_and_ffts.resize(req.minXSize,
                          req.xFftSizes,
                          req.minWorkSize);
        y.resize(req.minYSize,
                 req.minYAnticipatedWrites);

        int const ySz = static_cast<int>(y.y.size());
        int const nSteps = (ySz >= 1) ? 1 : 0;

        // set y progress such that results are written in a single chunk
        for(int i=0; i<nSteps; ++i) {
            y.increment();
        }
        
        phase_period.reset();
        
        auto f = [this] (auto const & innerP ){
            std::optional<float> const ph = innerP.getPhase();
            if(ph) {
                // There is at most one phase specification
                // within a set of params depending on the same context.
                Assert(!phase_period);
                phase_period = *ph;
            }
        };
        
        p.forEachUsingSameContext(f);
    }

    static int getAllocationSz_SetCoefficients(SetupParam const & p) {
        return State::getAllocationSz_SetCoefficients(p);
    }
    void setCoefficients(a64::vector<FPT> v, Algo const & algo)
    {
        state.setCoefficients(algo, std::move(v));
    }
    
    std::optional<float> getPhasePeriod() const {
        if(phase_period) {
            return phase_period;
        }
        int const maxHalfFFTSz = floor_power_of_two(x_and_ffts.fftHalfSizesBits);
        if(maxHalfFFTSz) {
            return maxHalfFFTSz;
        }
        return {};
    }
    
    void dephaseByGroupRatio(float r, Algo const & algo)
    {
        this->phase_group_ratio = r;
        
        dephase(algo);
        
        state.onContextFronteer([r](auto & s, auto const & fronteerAlgo){
            s.dephaseByGroupRatio(r, fronteerAlgo);
        });
    }

    bool isZero() const {
        return state.isZero();
    }
    double getEpsilon(Algo const & algo) const {
        return state.getEpsilon(algo);
    }
    void logComputeState(std::ostream & os, Algo const & algo) const {
        x_and_ffts.logComputeState(os);
        y.logComputeState(os);
        
        os << "phase: ";
        auto per = getPhasePeriod();
        if(per) {
            os << *per;
        }
        else {
            os << "none";
        }
        os << std::endl;
        os << "phase group ratio: " << phase_group_ratio << std::endl;
        IndentingOStreambuf i(os);
        state.logComputeState(algo, os);
    }
        
    template <typename Bool = bool>
    auto hasStepErrors() const -> std::enable_if_t<step_can_error, Bool> {
        return state.hasStepErrors();
    }

    FPT step(FPT val,
             Algo const & algo) {
        x_and_ffts.push(val);
        
        algo.step(state,
                  x_and_ffts,
                  y);
        
        FPT res = fft::RealSignal_<Tag, FPT>::get_signal(y.y[y.uProgress]);
        y.increment();
        return res;
    }
    
    template<typename FPT2>
    void stepAssignVectorized(FPT2 const * const input_buffer,
                              FPT2 * output_buffer,
                              int nSamples,
                              Algo const & algo)
    {
        for(int i=0; i<nSamples; ++i) {
            output_buffer[i] = step(input_buffer[i], algo);
        }
    }
    template<typename FPT2>
    void stepAddVectorized(FPT2 const * const input_buffer,
                           FPT2 * output_buffer,
                           int nSamples,
                           Algo const & algo)
    {
        for(int i=0; i<nSamples; ++i) {
            output_buffer[i] += step(input_buffer[i], algo);
        }
    }
    template<typename FPT2>
    void stepAddInputZeroVectorized(FPT2 * output_buffer,
                                    int nSamples,
                                    Algo const & algo)
    {
        for(int i=0; i<nSamples; ++i) {
            output_buffer[i] += step({}, algo);
        }
    }
    
    void flushToSilence(Algo const & algo) {
        algo.flushToSilence(state);

        {
            int const n_steps = x_and_ffts.flushToSilence();
            // in the future it will _not_ be needed to keep the 2 in sync.
            y.flushToSilence(n_steps);
        }
        
        dephase(algo);
    }
    
private:
    // 1 for all inputs
    XAndFFTS<FPT, Allocator, Tag, WorkCplxFreqs> x_and_ffts;
    // 1 for all outputs
    Y<FPT, Tag> y;

public:
    State state;
private:
    
    // cette phase lie 'x' a 'y', il faudra ouvrir la possibilité d'écrire dans y de facon non-contigue
    float phase_group_ratio = 0.; // in [0,1[
    std::optional<float> phase_period;

    void dephase(Algo const & algo) {
        auto per = getPhasePeriod();
        if(per) {
            const int n_steps = static_cast<int>(0.5f + phase_group_ratio * *per);

            algo.dephaseSteps(state, n_steps);

            for(int i=0; i<n_steps; ++i) {
                x_and_ffts.push(0);
                // in the future it will _not_ be needed to keep them in sync
                y.increment();
            }
        }
    }
};

// avoid false sharing (for use in async convolution)
template<typename Algo>
struct alignas(64) SelfContainedXYConvolution {
    
    using State = XYConvolution<Algo>;
    
    SelfContainedXYConvolution()
    : state(async_work)
    {}

    using WorkCplxFreqs = typename XYConvolution<Algo>::WorkCplxFreqs;
    using SetupParam = typename Algo::SetupParam;
    using FPT = typename Algo::FPT;
    using FFTTag = typename Algo::Tag;
    
    template<typename TT>
    using Allocator = typename State::template Allocator<TT>;
    
    void setupAndSetCoefficients(SetupParam const & p, a64::vector<FPT> v) {
        memory.clear();
        int const requiredSetup = State::getAllocationSz_Setup(p);
        int const requiredCoeffs = State::getAllocationSz_SetCoefficients(p);
        int const required = requiredCoeffs + requiredSetup;
        memory.reserve(required);

        {
            typename MemResource::use_type use(memory);

            algo.setup(p);
            if(algo.isValid()) {
                state.setup(p);
                state.setCoefficients(std::move(v), algo);
            }
        }
        
        if constexpr (MemResource::limited) {
            //std::cout << "mem monotonic " << memory.used() << std::endl;
            //std::cout << "mem monotonic " << memory.remaining() << std::endl;
            // verify that getAllocationSz_Setup / getAllocationSz_SetCoefficients are correct
            Assert(!algo.isValid() || required == memory.used()); // or there is padding to avoid false sharing?
        }
        
    }
    
    void dephaseByGroupRatio(float r) {
        state.dephaseByGroupRatio(r, algo);
    }

    std::optional<float> getPhasePeriod() const {
        return state.getPhasePeriod();
    }
    
    template<typename F>
    void onContextFronteer(F f) {
        f(state, algo);
    }

    bool isValid() const {
        return algo.isValid();
    }
    
    bool isZero() const {
        return state.isZero();
    }
    
    FPT step(FPT val) {
        return state.step(val, algo);
    }
    void flushToSilence()
    {
        state.flushToSilence(algo);
    }

    double getEpsilon() const {
        return state.getEpsilon(algo);
    }
    
    void logComputeState(std::ostream & os) const {
        state.logComputeState(os, algo);
    }
    
    bool handlesCoefficients() const {
        return algo.handlesCoefficients();
    }
    Latency getLatency() const {
        Assert(handlesCoefficients());
        return algo.getLatency();
    }
    
    int getBiggestScale() const {
        return algo.getBiggestScale();
    }


    using FBinsAllocator = typename fft::RealFBins_<FFTTag, FPT, Allocator>::type::allocator_type;
    using MemResource = MemResource<FBinsAllocator>;

    WorkCplxFreqs async_work;
    Algo algo;
    State state;
    typename MemResource::type memory;
    
    char padding[63];
};

}
