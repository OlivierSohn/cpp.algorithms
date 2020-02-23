
namespace imajuscule {

/*
 * In dsp.convolution.hpp we see this algorithm:
 *
 * The impulse response h is split in parts of equal length h1, h2, ... hn
 *
 *                     FFT(h1)
 *                        v
 *       +-----+    +-----------+  +-----+  +------+
 * x +-->| FFT |-|->| cplx mult |->| Add |->| IFFT |--> y
 *       +-----+ v  +-----------+  +-----+  +------+
 *          +--------+               ^
 *          |Delay(N)| FFT(h2)       |
 *          +--------+    v          |
 *               |  +-----------+    |
 *               |->| cplx mult |----|
 *               |  +-----------+    |
 *               .         .         .
 *               .         .         .
 *               .                   .
 *               v                   |
 *          +--------+               |
 *          |Delay(N)| FFT(hn)       |
 *          +----+---+    v          |
 *               |  +-----------+    |
 *               -->| cplx mult |----|
 *                  +-----------+
 *
 * 'PartitionnedFFTConvolutionCRTP' carries this monolithic computation
 * when needed, with no buffering, with a latency of the size of a partition.
 *
 * When partition sizes are an order of magnitude bigger than the audio callback buffer
 * (it is a very common case, since large partition
 * sizes is what makes this algorithm efficient on average),
 * many successive callback calls have very little work to do, and suddenly
 * a single callback call has to perform this huge monolithic computation.
 *
 * This is problematic because this audio callback can miss its deadline,
 * and then we'll hear a loud audio crack.
 *
 * To fix this, at the cost of a longer latency (twice the size of a partition size),
 * here we split the computation in several "grains" that can be computed at different times:
 *
 * - First there is the "FFT" grain which computes the fft of a chunk of the input signal
 *     (and does a little more than that, see the code)
 * - Then there are "multiplication" grains, where we multiply some delayed FFT of the input signal
 *    by some of the the FFT of the impulse response.
 *    The (max) number of vector multiplications per grain is the "multiplication group size".
 * - Finally there is the "IFFT" grain where we sum the results of the multiplications
 *    and do its inverse fft.
 *
 * An algorithm computes the optimal parameters (to minimize the worst case cost for a single callback):
 *  - The size of the partitions
 *  - The count of multiplications per multiplication grain
 *  - The "phasing" of simultaneous convolutions to best interleave
 *      high-cost grains.
 * based on :
 *  - The callback buffer size
 *  - The count of simultaneous convolutions happening in a callback
 */

enum class GrainType {
    FFT,
    IFFT,
    MultiplicationGroup
};

struct FinegrainedSetupParam : public Cost {
    static constexpr int nCoefficientsFadeIn = 0;
    static constexpr bool has_subsampling = false;

    explicit FinegrainedSetupParam() {}
    
    FinegrainedSetupParam(int partitionSz,
                          int partition_count,
                          int multiplication_group_size,
                          int phase)
    : Cost(phase)
    , multiplication_group_size(multiplication_group_size)
    , partition_size(partitionSz)
    , partition_count(partition_count)
    {}
    
    int multiplication_group_size = 0;
    int partition_size = 0;
    int partition_count = 0;
    
    void logSubReport(std::ostream & os) const override {
        os << "Finegrained, " << partition_count << " partitions of size " << partition_size << ", mult group size: " << multiplication_group_size << std::endl;
    }
    
    bool handlesCoefficients() const {
        return partition_size > 0;
    }
    
    void adjustWork(int targetNCoeffs) {
        partition_count = countPartitions(targetNCoeffs, partition_size);
        if(!handlesCoefficients()) {
            setCost(0.);
        }
        else {
            // we could reduce the cost...
        }
    }
    
    int get_fft_length() const {
        return 2 * partition_size;
    }

    int getGranularity() const {
        return partition_size / countGrains();
    }
    
    static constexpr int countNonMultiplicativeGrains() { return 2; }

    int countMultiplicativeGrains() const
    {
        return multiplication_group_size ?
        (1 + (partition_count - 1) / multiplication_group_size) :
        0;
    }
    
    int countGrains() const {
        return countMultiplicativeGrains() + countNonMultiplicativeGrains();
    }
    
    int getLowestValidMultiplicationsGroupSize() const {
        // lowest valid mult_grp_len verifies:
        
        // countGrains() == getBlockSize()
        // 2 + 1 + (ffts_of_partitionned_h.size() - 1)/mult_grp_len == partition_size
        
        auto constexpr min_number_grains = countNonMultiplicativeGrains() + 1;
        auto diff = partition_size - min_number_grains;
        if(diff < 0) {
            // invalid configuration
            return getHighestValidMultiplicationsGroupSize();
        }
        if(partition_count == 0) {
            return 0;
        }
        for(int i=1;; ++i) {
            if( (partition_count - 1)/i <= diff) {
                return i;
            }
        }
    }
    
    int getHighestValidMultiplicationsGroupSize() const {
        return partition_count;
    }
    
    template<Overlap Mode, typename FFTAlgo>
    MinSizeRequirement getMinSizeRequirement(int const maxVectorSz) const
    {
        if(!partition_size) {
            return {0,0,{},0};
        }
        Assert(maxVectorSz);
        int const n_max_ifftgrain_per_vector = countPartitions(maxVectorSz,
                                                               partition_size);
        int const n_max_partitions_touched_by_vector = 1 + countPartitions(maxVectorSz-1,
                                                                           partition_size);

        int const biggest_y_write = [this, n_max_ifftgrain_per_vector, maxVectorSz](){
            int const blockProgressForIFFTGrain = (countMultiplicativeGrains()+1) * getGranularity();
            Assert(blockProgressForIFFTGrain <= partition_size);
            int const gap = partition_size - blockProgressForIFFTGrain ;
            return gap + n_max_ifftgrain_per_vector * (get_fft_length()/2) + [this](){
                if constexpr(Mode == Overlap::Add) {
                    return get_fft_length()/2;
                }
                else {
                    return 0;
                }
            }();
        }();
        
        // during the worst iteration, we can:
        // - write at most 'biggest_y_write' to y
        // - and _then_ read 'maxVectorSz' from y
        int const y_size = biggest_y_write + maxVectorSz;
        
        return {
            0, // x block size
            y_size, // y block size
            {
                {
                    get_fft_length(),
                    std::max(0,
                             partition_count + (n_max_partitions_touched_by_vector-1) + (multiplication_group_size - 1))
                }
            },
            FFTAlgo::inplace_dft ? 0 : get_fft_length() // work size
        };
    }
    
    static constexpr Latency getLatencyForPartitionSize(int sz) {
        return Latency(2 * sz - 1);
    }

    Latency getImpliedLatency() const {
        Assert(handlesCoefficients());
        return getLatencyForPartitionSize(partition_size);
    }
    
    int const getBiggestScale() const {
        return partition_size;
    }
    
    template<typename F>
    void forEachUsingSameContext(F f) const {
        f(*this);
    }

    static FinegrainedSetupParam makeInactive() {
        FinegrainedSetupParam res{0,0,0,0};
        res.setCost(0.f);
        return res;
    }
};

template<typename T, typename Tag>
struct FinegrainedPartitionnedFFTConvolutionSimulation {
    using SetupParam = FinegrainedSetupParam;
    
    void setup(SetupParam const & p) {
        reset();
        this->p = p;

        partition_size_minus_one = p.partition_size - 1;

        count_multiplicative_grains = p.countMultiplicativeGrains();
        granularity = p.getGranularity();
        
        costXForwardFft = fft::AlgoCosts<Tag,T>::cost_fft_forward(get_fft_length());
        costXInverseFft = fft::AlgoCosts<Tag,T>::cost_fft_inverse(get_fft_length());
        cost_mult = fft::RealFBinsCosts<Tag,T>::cost_multiply(p.partition_size);
        cost_mult_add = fft::RealFBinsCosts<Tag,T>::cost_multiply_add(p.partition_size);
        cost_add_assign = fft::RealSignalCosts<Tag,T>::cost_add_assign(get_fft_length());
    }
    
    auto const & getParam() const {
        return p;
    }
    
    void setMultiplicationGroupLength(int l) {
        auto p2 = p;
        p2.multiplication_group_size = l;
        setup(p2);
    }
    
    int getBlockSize() const {
        return 1+partition_size_minus_one;
    }

    bool isValid() const {
        if(p.multiplication_group_size == 0) {
            return p.partition_count == 0;
        }
        return p.countGrains() <= p.partition_size;
    }
    
    bool isZero() const {
        return p.partition_count == 0;
    }

    void reset() {
        count_multiplicative_grains = 0;
        granularity = 0;
        
        reset_states();
    }
    
    /*
     Dual method of FinegrainedPartitionnedFFTConvolutionSimulation::step()
     */
    double simuStep(XFFTsCostsFactors const & xFftCostFactors) {
        double cost = 0.;
        if(unlikely(isZero())) {
            return cost;
        }
        ++grain_counter;
        auto g = nextGrain();
        assert(g.first >= 0);
        if(unlikely(g.first == 0)) {
            cost += simuDoGrain(g.second,
                                xFftCostFactors);
            updatePostGrain(g.second);
            grain_counter = 0;
        }
        ++x_progress; // progress may overflow, that's ok because we use only its lower bits.

        return cost;
    }
    
private:
    int32_t grain_counter = 0;
    int32_t grain_number = 0;
    uint32_t x_progress = 0;

    SetupParam p;
    
    int32_t count_multiplicative_grains = 0;
    int32_t granularity = 0;
    uint32_t partition_size_minus_one = -1;

    double costXForwardFft = 0.;
    double costXInverseFft = 0.;
    double cost_mult = 0.;
    double cost_mult_add = 0.;
    double cost_add_assign = 0;
    
    void reset_states() {
        grain_counter = 0;
        grain_number = 0;
        x_progress = 0;
    }
    
    void updatePostGrain(GrainType g) {
        if(g==GrainType::FFT) {
            this->grain_number = 0;
        }
        else {
            ++this->grain_number;
        }
    }
    
    std::pair<int, GrainType> nextGrain() const {
        
        auto const n_mult_grains_remaining = count_multiplicative_grains - grain_number;
        
        if(unlikely(n_mult_grains_remaining < 0)) {
            int distToFFTGrain = partition_size_minus_one - (x_progress & partition_size_minus_one);
            assert(distToFFTGrain >= 0);
            return {distToFFTGrain, GrainType::FFT};
        }
        
        int const dist = granularity - grain_counter;
        
        assert(dist >= 0);
        if(unlikely(n_mult_grains_remaining == 0)) {
            return {dist, GrainType::IFFT};
        }
        else {
            Assert(n_mult_grains_remaining > 0);
            return {dist, GrainType::MultiplicationGroup};
        }
    }
    
    double simuDoGrain(GrainType g,
                       XFFTsCostsFactors const & xFftCostFactors) const
    {
        double cost {};

        auto const fft_length = get_fft_length();

        if(g == GrainType::FFT) {
            Assert(0 == ((x_progress+1) & partition_size_minus_one)); // make sure 'rythm is good'

            // the fft is done by x_and_ffts
            cost += costXForwardFft * xFftCostFactors.get(fft_length);
        }
        else {
            if(g == GrainType::MultiplicationGroup) {
                auto const M = getMultiplicationsGroupMaxSize();
                int offset = M * grain_number;
                assert(offset >= 0);
                assert(offset < p.partition_count);

                int offset_end = std::min(offset + M,
                                          p.partition_count);
                if(offset == 0) {
                    cost += cost_mult;
                    offset = 1;
                }
                int const nRemaining = std::max(0,
                                                offset_end-offset);
                cost += nRemaining * cost_mult_add;
            }
            else {
                Assert(g==GrainType::IFFT);

                cost += costXInverseFft;
                cost += cost_add_assign;
            }
        }
        return cost;
    }
    
    auto get_fft_length() const { return 2 * (1+partition_size_minus_one); }
    auto getMultiplicationsGroupMaxSize() const { return p.multiplication_group_size; }
};
    
template<typename T, typename FFTTag>
struct Simulation_<FinegrainedSetupParam, T, FFTTag> {
    using type = FinegrainedPartitionnedFFTConvolutionSimulation<T, FFTTag>;
};
  
constexpr Latency minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter {
    ScaleConvolution_::latencyForDroppedConvolutions(ScaleConvolution_::nDroppedOptimalFor_Split_Bruteforce_Fft)
};
    
template<typename SetupParam>
int constexpr getMinLg2PartitionSz() {
    int partition_sz = 1;
    for(;;partition_sz *= 2) {
        if(SetupParam::getLatencyForPartitionSize(partition_sz) >= minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter) {
            break;
        }
    }
    
    return power_of_two_exponent(partition_sz);
}

template<typename C1, typename C2>
void add_cyclic(C1 & c1, C2 const & c2) {
    // verify that there is a n >= 1 such that c..size() == n*c2.size()
    Assert((c1.size() / c2.size()) * c2.size() == c1.size());

    auto it2 = c2.begin();
    auto end2 = c2.end();
    for(auto it=c1.begin(), end = c1.end(); it!=end; ++it, ++it2) {
        if(it2 == end2) {
            it2 = c2.begin();
        }
        *it += *it2;
    }
}

/*
 input parameters :
 - sz_audio_cb, n_phaseable_channels
 - impulse response length
 
 output parameters:
 - lg(partition size)
 
 1D - Gradient descent according to cost 'max grain computation time' with variable parameters 'lg(partition_size)'
 deducing 'number of multiplications per grain' by finding the parameter that leads to grain computations time just below max(fft, ifft),
 with the constraint that one computation at most occurs per 'equivalent' audio callback.
 */
template<typename Sim, typename SetupParam = typename Sim::SetupParam, typename GradientDescent = GradientDescent<SetupParam>>
auto find_optimal_partition_size(GradientDescent & gradient_descent,
                                 int const n_iterations,
                                 int const n_phaseable_channels,
                                 int const n_scales,
                                 int const sz_audio_cb,
                                 int const zero_latency_response_size,
                                 XFFTsCostsFactors const & xFFTCostFactors,
                                 std::function<int(int, std::vector<double> & )> const & get_early_monochannel_costs,
                                 std::ostream & os) -> std::optional<SetupParam>
{
    if(n_scales != 1) {
        throw std::runtime_error("optimization with scales not implemented");
    }

    std::vector<double> frames_costs, early_monochannel_costs;
    std::vector<double> phased_sums;
    int const sz_reserve = std::min(32768,
                                    static_cast<int>(ceil_power_of_two(zero_latency_response_size)));
    early_monochannel_costs.reserve(sz_reserve);
    frames_costs.reserve(sz_reserve);
    phased_sums.reserve(sz_reserve);
    
    // inclusive upper bound
    int const max_lg2_partition_size = power_of_two_exponent(ceil_power_of_two(zero_latency_response_size));

    auto bestMultiplicationGroupSizeForPartitionSize =
    [sz_audio_cb,
     n_iterations,
     n_phaseable_channels,
     &os,
     &xFFTCostFactors,
     &frames_costs,
     &phased_sums,
     &early_monochannel_costs,
     &get_early_monochannel_costs,
     zero_latency_response_size,
     max_lg2_partition_size]
    (int const lg2_partition_size, auto & val)
    {
        using namespace profiling;
        using namespace std;
        using namespace std::chrono;
        
        if(lg2_partition_size < 0) {
            return ParamState::OutOfRange;
        }
        if(lg2_partition_size > max_lg2_partition_size) {
            return ParamState::OutOfRange;
        }
        int const partition_size = pow2(lg2_partition_size);
        //            cout << "partition size : " << partition_size << endl;
        
        Assert(partition_size <= ceil_power_of_two(zero_latency_response_size));
        int const n_early_coeffs = get_early_monochannel_costs(partition_size,
                                                               early_monochannel_costs);
        int const length_impulse = zero_latency_response_size - n_early_coeffs;
        
        // the value for multiplication group size is not very important (it will be overriden later on)
        // but needs to lead to a valid FinegrainedSetupParam. We use the highest valid number:
        int const n_partitions = countPartitions(length_impulse,
                                                 partition_size);
        
        Sim sim;
        sim.setup({
            partition_size,
            n_partitions,
            n_partitions,
            0
        });
        
        Assert(sim.getBlockSize() == partition_size);
        
        if(!sim.isValid()) {
            return ParamState::OutOfRange;
        }
        
        frames_costs.resize(partition_size);
        phased_sums.resize(partition_size);

        struct CostForPhase {
            /*
             We favor bigger multiplication groups because in reality,
             they lead to better performances due to instruction cache / data cache / prefetching
             (these effects are not taken into account in virtual simulation)
             */
            double getCost() const {
                constexpr double grp_sz_max_influence = 0.01;
                
                double const normalized_grp_sz = mult_grp_sz / static_cast<double>(std::max(1, count_partitions));
                Assert(normalized_grp_sz >= 0.);
                Assert(normalized_grp_sz <= 1.);
                return virtual_simulation_cost * (1. - grp_sz_max_influence * normalized_grp_sz);
            }
            double virtual_simulation_cost = 0.;
            double phase = 0.;
            int mult_grp_sz = 0;
            int count_partitions = 0;
        };
        
        RangedGradientDescent<CostForPhase> rgd([sz_audio_cb, n_phaseable_channels, &frames_costs, &phased_sums, &sim, &xFFTCostFactors, &early_monochannel_costs](int multiplication_group_size, auto & result)
        {
            sim.setMultiplicationGroupLength(multiplication_group_size);
            if(!sim.isValid()) {
                return ParamState::OutOfRange;
            }

            result.mult_grp_sz = multiplication_group_size;
            result.count_partitions = sim.getParam().partition_count;
            
            for(int i=0, end=frames_costs.size(); i<end; ++i) {
                frames_costs[i] = sim.simuStep(xFFTCostFactors);
            }
            add_cyclic(frames_costs, early_monochannel_costs);
            
            result.virtual_simulation_cost = computeMaxSlidingSum(frames_costs,
                                                                  sz_audio_cb);
            result.phase = 0.;
            result.virtual_simulation_cost *= n_phaseable_channels;
            // now cost is the 'phase == 0' cost

            int const granularity = sim.getParam().getGranularity();
            
            if(n_phaseable_channels >= 2) {
                int const nMinFullCbInGranularity = granularity / sz_audio_cb;
                if(nMinFullCbInGranularity >= n_phaseable_channels) {
                    // naive phasing leads to an optimal solution : with a phase of 'sz_audio_cb',
                    // during every callback call there will be _at most_ one grain among all channels
                    // that will be computed.
                    result.phase = sz_audio_cb;

                    compute_phased_sum(frames_costs,
                                       result.phase,
                                       n_phaseable_channels,
                                       phased_sums);
                    result.virtual_simulation_cost = computeMaxSlidingSum(phased_sums,
                                                                          sz_audio_cb);
                }
                else {
                    // naive phasing is not enough to have an optimal solution
                    // because the grains are too close to each other

                    // the peaks of the finegrained part have a period of "granularity"
                    // the _main_ peaks of the early part have a period of "early_monochannel_costs.size()" (aka getBiggestScale)
                    // so we use these base phases:
                    
                    float const base_phase_early = static_cast<float>(early_monochannel_costs.size()) / n_phaseable_channels;
                    float const base_phase_late = static_cast<float>(granularity) / n_phaseable_channels;
                    float const base_phase_mix = 0.5f * (base_phase_late + base_phase_early);
                    float const base_phase_zero = 0.f;
                    
                    std::set<int> base_phases {
                        static_cast<int>(0.5f + base_phase_early),
                        static_cast<int>(0.5f + base_phase_late),
                        static_cast<int>(0.5f + base_phase_mix),
                        static_cast<int>(0.5f + base_phase_zero)
                    };
                    
                    for(auto base_phase : base_phases) {
                        
                        int const maxPhase = 1 + (frames_costs.size() / n_phaseable_channels);
                        int const phaseStep = std::max(1, granularity);
                        for(int phase = base_phase;
                            phase <= maxPhase;
                            phase += phaseStep)
                        {
                            if(0 == phase) {
                                // already tested
                                continue;
                            }
                            compute_phased_sum(frames_costs,
                                               phase,
                                               n_phaseable_channels,
                                               phased_sums);
                            
                            auto phased_cost = computeMaxSlidingSum(phased_sums,
                                                                    sz_audio_cb);
                            if(phased_cost < result.virtual_simulation_cost) {
                                result.virtual_simulation_cost = phased_cost;
                                result.phase = phase;
                            }
                        }
                    }
                }
            }
            
            // the phase is now the difference between two close phaseable,
            // but we need it to be multiplied by the number of phaseable:
            result.phase *= n_phaseable_channels;
            
            // cost is now the sum of costs of each phaseable_channel
            // but cost should be per sample, not per frame, so
            // we divide by the number of channels
            result.virtual_simulation_cost /= static_cast<float>(n_phaseable_channels);

            result.virtual_simulation_cost /= sz_audio_cb;
            // cost == 'worst computation time over one callback, averaged per frame of phaseable channels'
            constexpr bool display = false;
            if constexpr (display) {
                std::cout <<
                "n_phaseable_channels:      " << n_phaseable_channels << std::endl <<
                "sz_audio_cb:     " << sz_audio_cb << std::endl <<
                "partition_sz:    " << sim.getParam().partition_size << std::endl <<
                "partition_count: " << sim.getParam().partition_count << std::endl <<
                "mult_sz:         " << multiplication_group_size << std::endl <<
                "granularity:     " << granularity << std::endl <<
                "phase:           " << result.phase << std::endl <<
                "cost:            (" << result.virtual_simulation_cost << ") " << result.getCost() << std::endl <<
                std::endl;
                {
                    StringPlot plot(20, frames_costs.size());
                    plot.draw(frames_costs, '+');
                    plot.log();
                }
                {
                    StringPlot plot(20, early_monochannel_costs.size());
                    plot.draw(early_monochannel_costs, '+');
                    plot.log();
                }
            }

            return ParamState::Ok;
        });
        
        range<int> const range_multiplication_group_length {
            sim.getParam().getLowestValidMultiplicationsGroupSize(),
            sim.getParam().getHighestValidMultiplicationsGroupSize()
        };
        
        CostForPhase best;
        val.multiplication_group_size = rgd.findLocalMinimum(n_iterations, range_multiplication_group_length, best);
        val.setCost(best.getCost());
        val.setPhase(best.phase);
        val.partition_size = partition_size;
        val.partition_count = countPartitions(length_impulse, partition_size);
        
        constexpr auto debug = false;
        if(debug) {
            std::cout << "partition_size:" << partition_size << std::endl;
            rgd.plot(true, os);
            rgd.make_exhaustive(range_multiplication_group_length, os);
            rgd.plot(true, os);
        }
        return ParamState::Ok;
    };
    
    // must yield a valid result:
    int const first_lg2_partitionsz = std::min(max_lg2_partition_size,
                                               getMinLg2PartitionSz<SetupParam>());
    gradient_descent.setFunction(bestMultiplicationGroupSizeForPartitionSize);

    Optional<SetupParam> min_val;
    auto res = gradient_descent.findLocalMinimum(n_iterations,
                                                 first_lg2_partitionsz,
                                                 min_val);
    if(res) {
        assert(min_val);
        assert(min_val->partition_size == pow2(*res));
    }
    return min_val;
}
    
template<typename T, typename Tag>
struct PartitionAlgo< FinegrainedSetupParam, T, Tag > {
    using Sim = FinegrainedPartitionnedFFTConvolutionSimulation<T, Tag>;
    using SetupParam = FinegrainedSetupParam;
    
    static std::optional<SetupParam> run(int const n_phaseable_channels,
                                         int const n_scales,
                                         int const n_audio_frames_per_cb,
                                         int const zero_latency_response_size,
                                         std::ostream & os,
                                         XFFTsCostsFactors const & xFftCostFactors,
                                         std::function<int(int, std::vector<double> &)> const & get_early_monochannel_costs)
    {
        os << "Optimization of FinegrainedPartitionnedFFTConvolution for " << n_scales << " scale(s)" << std::endl;
        IndentingOStreambuf i(os);
        
        Assert(n_phaseable_channels > 0);
        GradientDescent<SetupParam> gd;
        constexpr auto n_iterations = 1;
        std::optional<SetupParam> res = find_optimal_partition_size<Sim>(gd,
                                                                         n_iterations,
                                                                         n_phaseable_channels,
                                                                         n_scales,
                                                                         n_audio_frames_per_cb,
                                                                         zero_latency_response_size,
                                                                         xFftCostFactors,
                                                                         get_early_monochannel_costs,
                                                                         os);
        constexpr auto debug_gradient_descent = false;
        if constexpr (debug_gradient_descent) {
            os << "Gradient descent report :" << std::endl;
            gd.debug(true, os); // this takes time because it makes the search exhaustive
        }
        return res;
    }
};


static inline std::ostream& operator<<(std::ostream& os, const FinegrainedSetupParam& p)
{
    using namespace std;
    os
    << "multiplication group size : " << p.multiplication_group_size << endl;
    return os;
}

}
