

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
        return static_cast<int>(ffts.size());
    }
private:
    cyclic<FBins> ffts;
};


struct XSegments {
    int start;
    int size_from_start; // from start
    int size_from_zero; // from 0
    
    // The invariant is that if size_from_zero is not 0, then size_from_zero is not 0.
};

template<typename T, template<typename> typename Allocator, typename Tag>
struct XAndFFTS {
    static constexpr auto zero_n_raw = fft::RealSignal_<Tag, T>::zero_n_raw;

    XSegments getPastSegments(int const pastSize) const {
        int start = progress + 1 - pastSize;
        while(start < 0) {
            start += x_unpadded_size;
        }
        Assert(start >= 0);
        
        int const szFromZero = std::max(0,
                                        start+pastSize-x_unpadded_size);
        return {
            start,
            pastSize-szFromZero,
            szFromZero
        };
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

    void resize(int const sz, FftSpecs const & fftSpecs) {
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
        Assert(count_set_bits(fftHalfSizesBits) == static_cast<int>(x_ffts.size()));
        int const maxHalfFFTSz = floor_power_of_two(fftHalfSizesBits);

        x_unpadded_size = std::max(maxHalfFFTSz,
                                   static_cast<int>(ceil_power_of_two(std::max(1,
                                                                               sz))));

        x.clear();
        x.resize(x_unpadded_size + maxHalfFFTSz); // add padding for biggest fft
        
        padding = static_cast<int>(x.size());
        progress = x_unpadded_size-1;
    }
    
    template<typename SetupParam>
    void setPhasePeriod(SetupParam const & p) {
        Assert(countPhaseable > 0);
        
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
    
    std::optional<float> getPhasePeriod() const {
        if(phase_period) {
            return phase_period;
        }
        int const maxHalfFFTSz = floor_power_of_two(fftHalfSizesBits);
        if(maxHalfFFTSz) {
            return maxHalfFFTSz;
        }
        return {};
    }
    
    template<typename F>
    void dephase(F f) {
        auto per = getPhasePeriod();
        if(per) {
            const int n_steps = static_cast<int>(0.5f + phase_group_ratio * *per);

            for(int i=0; i<n_steps; ++i) {
                push(0);
                f(progress);
            }
        }
    }
    
    void push(T v) {
        ++progress;
        if(unlikely(progress == x_unpadded_size)) {
            progress = 0;
            padding = 0;
        }
        x[progress] = typename RealSignal::value_type(v);
        int const next_progress = progress+1;
        padding = std::max(padding, next_progress);
        if(unlikely(padding == x_unpadded_size)) {
            padding = static_cast<int>(x.size());
        }
        
        uint32_t nZeroes = count_trailing_zeroes(next_progress);
        uint32_t mask = (1 << (nZeroes+1)) - 1;
        fftsHalfSizesBitsToCompute = fftHalfSizesBits & mask;
        for(int i=0, nComputedFFTs = count_set_bits(fftsHalfSizesBitsToCompute);
            i<nComputedFFTs;
            ++i)
        {
            auto & f = x_ffts[i];
            int halfFftSize = f.fft_length/2;
            int xStart = next_progress - halfFftSize;
            Assert(xStart >= 0);

            int const xEnd = next_progress + halfFftSize;
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
        Assert(progress < x_unpadded_size);
        int const progressBackup = progress;
        progress = x_unpadded_size-1;
        padding = x.size();
        fftsHalfSizesBitsToCompute = 0;
        Assert(progressBackup <= x_unpadded_size);
        return progress - progressBackup;
    }
    
    using RealSignal = typename fft::RealSignal_<Tag, T>::type;
    RealSignal x;
    
    int32_t progress = 0; // the last position that has been written to
    int32_t padding = 0; // the position of the first non-zero (non-padded) element
    int32_t x_unpadded_size = 0;
    uint32_t fftHalfSizesBits = 0; // if we use the fft of sz 2^(n+1), the nth bit is set
    uint32_t fftsHalfSizesBitsToCompute = 0;
    
    std::vector<FFTs<T, Allocator, Tag>> x_ffts;

    float phase_group_ratio = 0.; // in [0,1[
    std::optional<float> phase_period;

    auto const & find_ffts(int sz) const {
        for(auto const & f : x_ffts) {
            if(sz == f.fft_length) {
                return f;
            }
        }
        throw std::runtime_error("fft not found");
    }

    void logComputeState(std::ostream & os) const {
        os << "Source:" << std::endl;
        IndentingOStreambuf in(os);
        
        os << "x size:" << x.size() <<  " unpadded:" << x_unpadded_size << " progress:" << progress << std::endl;
        
        os << "ffts: ";
        for(auto const & f:x_ffts) {
            os << f.size() << "x" << f.fft_length << " ";
        }
        os << std::endl;

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
    }
};

struct YSegments {
    int size_from_progress;
    int size_from_zero;

    // The invariant is that if size_from_progress is not 0, then size_from_zero is not 0.
};

template<typename T, typename Tag>
struct Y {
    static constexpr auto zero_n_raw = fft::RealSignal_<Tag, T>::zero_n_raw;
    static constexpr auto get_signal = fft::RealSignal_<Tag, T>::get_signal;
    static constexpr auto copy = fft::RealSignal_<Tag, T>::copy;
    static constexpr auto add_assign = fft::RealSignal_<Tag, T>::add_assign;

    using RealSignal = typename fft::RealSignal_<Tag, T>::type;
    using RealValue = typename RealSignal::value_type;

    T getCurrentSignal() const
    {
        return get_signal(y[progress]);
    }

    void addAssign(int unbounded_future_progress, // >= uProgress anb maybe > ySz
                   typename RealSignal::value_type * x,
                   int N)
    {
        // zero from
        //  uProgress + write_sz
        // to
        //  futureProgress

        int const unbounded_first_old = progress + write_sz;
        int write_sz_from_future_progress = unbounded_first_old - unbounded_future_progress;
        {
            if(write_sz_from_future_progress < 0) {
                int const first_old = (unbounded_first_old < ySz) ? unbounded_first_old : (unbounded_first_old - ySz);
                Assert(first_old < ySz);
                int const diff = first_old - write_sz_from_future_progress - ySz;
                if(diff > 0) {
                    zero_n_raw(&y[first_old],
                               (-write_sz_from_future_progress)-diff);
                    zero_n_raw(&y[0],
                               diff);
                }
                else {
                    zero_n_raw(&y[first_old],
                               -write_sz_from_future_progress);
                }
                write_sz_from_future_progress = 0;
            }
        }
        
        int const future_progress = (unbounded_future_progress < ySz) ? unbounded_future_progress : (unbounded_future_progress - ySz);
        Assert(future_progress < ySz);

        auto s = getFutureSegments(future_progress, N);
        if(likely(s.size_from_progress)) {
            
            if(write_sz_from_future_progress >= s.size_from_progress) {
                add_assign(&y[future_progress],
                           &x[0],
                           s.size_from_progress);
            }
            else if(write_sz_from_future_progress == 0) {
                copy(&y[future_progress],
                     &x[0],
                     s.size_from_progress);
            }
            else {
                Assert(write_sz_from_future_progress < s.size_from_progress);
                add_assign(&y[future_progress],
                           &x[0],
                           write_sz_from_future_progress);
                copy(&y[future_progress + write_sz_from_future_progress],
                     &x[write_sz_from_future_progress],
                     s.size_from_progress-write_sz_from_future_progress);
            }
            
            if(unlikely(s.size_from_zero)) {
                if(write_sz_from_future_progress >= N) {
                    add_assign(&y[0],
                               &x[s.size_from_progress],
                               s.size_from_zero);
                }
                else {
                    int const write_sz_from_zero = std::max(0,
                                                            future_progress + write_sz_from_future_progress - ySz);
                    x += s.size_from_progress;
                    if(write_sz_from_zero == 0) {
                        copy(y.data(),
                             x,
                             s.size_from_zero);
                    }
                    else {
                        Assert(write_sz_from_future_progress < N);
                        add_assign(&y[0],
                                   x,
                                   write_sz_from_zero);
                        copy(&y[write_sz_from_zero],
                             x + write_sz_from_zero,
                             s.size_from_zero - write_sz_from_zero);
                    }
                }
            }
            int const minWriteSz = unbounded_future_progress + N - progress;
            write_sz = std::max(write_sz, minWriteSz);
        }
    }
    
    void addAssignPresent(typename RealSignal::value_type * x,
                          int N)
    {
        auto s = getFutureSegments(progress, N);
        if(likely(s.size_from_progress)) {
            
            if(write_sz >= s.size_from_progress) {
                add_assign(&y[progress],
                           &x[0],
                           s.size_from_progress);
            }
            else {
                if(write_sz == 0) {
                    copy(&y[progress],
                         &x[0],
                         s.size_from_progress);
                }
                else {
                    Assert(write_sz < s.size_from_progress);
                    add_assign(&y[progress],
                               &x[0],
                               write_sz);
                    copy(&y[progress + write_sz],
                         &x[write_sz],
                         s.size_from_progress-write_sz);
                }
            }
            
            if(unlikely(s.size_from_zero)) {
                if(write_sz >= N) {
                    add_assign(&y[0],
                               &x[s.size_from_progress],
                               s.size_from_zero);
                }
                else {
                    int const write_sz_from_zero = std::max(0,
                                                            progress + write_sz - ySz);
                    x += s.size_from_progress;
                    if(write_sz_from_zero == 0) {
                        copy(y.data(),
                             x,
                             s.size_from_zero);
                    }
                    else {
                        Assert(write_sz < N);
                        add_assign(&y[0],
                                   x,
                                   write_sz_from_zero);
                        copy(&y[write_sz_from_zero],
                             x + write_sz_from_zero,
                             s.size_from_zero - write_sz_from_zero);
                    }
                }
            }

            write_sz = std::max(write_sz, N);
        }
    }
    
private:
    YSegments getFutureSegments(int const from, int const future) const {
        int const end = from + future;
        int const szFromZero = std::max(0,
                                        end - ySz);
        return {
            future - szFromZero,
            szFromZero
        };
    }
public:
    
    void writeOne(RealValue val)
    {
        if(write_sz) {
            y[progress] += val;
        }
        else {
            y[progress] = val;
            write_sz = 1;
        }
    }

    void resize(int const sz_) {
        ySz = std::max(1, sz_);
        y.clear();
        y.resize(ySz); // intialized at 0

        progress = 0;
        write_sz = 0;
    }
    
    void increment() {
        write_sz = std::max(0, write_sz-1);
        ++progress;
        Assert(progress <= ySz);
        if(progress == ySz) {
            progress = 0;
        }
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
        os << "y progress " << progress << std::endl;
    }
private:
    RealSignal y;
public:
    int32_t progress = 0;
    int32_t ySz=0;
    int32_t write_sz = 0; // the number of locations, starting from uProgress, that are _not_ old samples.
};

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
    using WorkData = typename WorkCplxFreqs::value_type;
    
    using SetupParam = typename Algo::SetupParam;

    static int getAllocationSz_Setup(SetupParam const & p) {
        MinSizeRequirement req = p.getMinSizeRequirement();
        return XAndFFTS<FPT, Allocator, Tag>::getAllocationSz_Resize(req.xFftSizes);
    }
    
    std::optional<float> getPhasePeriod() const {
        return x_and_ffts.getPhasePeriod();
    }

    void setup(SetupParam const & p,
               WorkCplxFreqs & work)
    {
        MinSizeRequirement req = p.getMinSizeRequirement();
        
        work.reserve(req.minWorkSize);
        
        x_and_ffts.resize(req.minXSize,
                          req.xFftSizes);
        x_and_ffts.setPhasePeriod(p);
        y.resize(req.minYSize);

        // set y progress such that results are written in a single chunk
        if(y.ySz >= 1) {
            y.increment();
        }
    }

    static int getAllocationSz_SetCoefficients(SetupParam const & p) {
        return State::getAllocationSz_SetCoefficients(p);
    }
    void setCoefficients(a64::vector<FPT> v, Algo const & algo)
    {
        state.setCoefficients(algo, std::move(v));
    }
    
    void dephaseByGroupRatio(float r, Algo const & algo)
    {
        x_and_ffts.phase_group_ratio = r;
        
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
        
        IndentingOStreambuf i(os);
        state.logComputeState(algo, os);
    }
        
    template <typename Bool = bool>
    auto hasStepErrors() const -> std::enable_if_t<step_can_error, Bool> {
        return state.hasStepErrors();
    }

    FPT step(FPT val,
             Algo const & algo,
             WorkData * workData) {
        x_and_ffts.push(val);
        
        algo.step(state,
                  x_and_ffts,
                  y,
                  workData);
        
        FPT res = y.getCurrentSignal();
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
            y.flushToSilence(n_steps);
        }
        
        dephase(algo);
    }
    
private:
    // 1 for all inputs
    XAndFFTS<FPT, Allocator, Tag> x_and_ffts;
    // 1 for all outputs
    Y<FPT, Tag> y;

public:
    State state;

private:
    void dephase(Algo const & algo)
    {
        x_and_ffts.dephase([&algo, this](int x_progress) {
            algo.dephaseStep(state,
                             x_progress);
            y.increment();
        });
    }
};

// avoid false sharing (for use in async convolution)
template<typename Algo>
struct alignas(64) SelfContainedXYConvolution {
    
    using State = XYConvolution<Algo>;
    
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
                state.setup(p, work);
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
        return state.step(val,
                          algo,
                          work.data());
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

    WorkCplxFreqs work;
    Algo algo;
    State state;
    typename MemResource::type memory;
    
    char padding[63];
};

}
