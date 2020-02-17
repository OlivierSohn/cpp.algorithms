
namespace imajuscule {

template <typename Parent>
struct FFTConvolutionIntermediateSimulation : public Parent {
    using FPT = typename Parent::FPT;
    using FFTTag = typename Parent::FFTTag;
    using Parent::doSetCoefficientsCount;
    static constexpr bool has_subsampling = false;
    
    using Parent::cost_compute_convolution;
    using Parent::get_fft_length;
    using Parent::getBlockSize;
    
    void setCoefficientsCount(int64_t szCoeffs) {
        doSetCoefficientsCount(szCoeffs);
    }
    
    FFTConvolutionIntermediateSimulation()
    : minorCost(costWriteNConsecutive<FPT*>(1) +
                costReadNConsecutive<FPT>(1))
    {}
    
    double simuMinorStep() {
        return minorCost;
    }
    
    double simuMajorStep(XFFTsCostsFactors const & xFftCostFactors) {
        if(unlikely(!majorCostExceptXForwardFft)) {
            majorCostExceptXForwardFft =
            
            // assuming that cost_compute_convolution is CONSTANT
            cost_compute_convolution() +
            
            fft::AlgoCosts<FFTTag, FPT>::cost_fft_inverse(get_fft_length()) +
            
            // if get_fft_length() has the same order of magnitude than the size of y:
            // - many 'add_assign' will be replaced by 'copy'
            // - many 'add_assign' will be done in 2 separate passes
            //
            // To take these effects fully into account we would need to know
            // - the size of y
            // - the phase at which we will write in y
            // - the rythm of other convolution parts (and other channels) writing in the same y
            //
            // "This is complicated to do" and we tend to think that the magnitude of the error
            // will be small compared to the overall sum, hence we simulate using a single 'add_assign':
            
            fft::RealSignalCosts<FFTTag, FPT>::cost_add_assign(get_fft_length()) +
            
            costReadNConsecutive<FPT>(1);
        }
        if(unlikely(!costXForwardFft)) {
            costXForwardFft = fft::AlgoCosts<FFTTag, FPT>::cost_fft_forward(get_fft_length());
        }
        //LG(INFO, "*:%f fft: %f", *majorCostExceptXForwardFft, *costXForwardFft);
        return
        *majorCostExceptXForwardFft +
        *costXForwardFft * xFftCostFactors.get(get_fft_length());
    }
private:
    double minorCost;
    std::optional<double> majorCostExceptXForwardFft;
    std::optional<double> costXForwardFft;
};

struct FFTConvolutionCRTPSetupParam : public Cost
{
    static constexpr bool has_subsampling = false;
    static constexpr int nCoefficientsFadeIn = 0;

    FFTConvolutionCRTPSetupParam(int blockSize)
    : blockSize(blockSize)
    {}
    
    void logSubReport(std::ostream & os) const override {
        os << "FFTConvolutionCRTPSetupParam block size " << blockSize << std::endl;
    }
    
    bool handlesCoefficients() const {
        return blockSize > 0;
    }
    
    template<Overlap Mode>
    MinSizeRequirement getMinSizeRequirement(int const maxVectorSz) const
    {
        if(!blockSize) {
            return {0,0,{},0};
        }
        int const n_partitions_in_vec = countPartitions(maxVectorSz, blockSize);
        Assert(n_partitions_in_vec);

        int const fft_length = 2 * blockSize;
        
        // when we have vectorization, we may need to write in the future of y
        int const yAdditionalWrite = (maxVectorSz > 1) ? (fft_length/2) : 0;
        
        int const y_base_write_size = (maxVectorSz-1) + n_partitions_in_vec * (fft_length/2) + [fft_length](){
            if constexpr(Mode == Overlap::Add) {
                return fft_length/2;
            }
            else {
                return 0;
            }
        }();
        
        int const biggest_y_write = y_base_write_size + yAdditionalWrite;
        
        // during the worst iteration, we can:
        // - write at most 'biggest_y_write' to y
        // - and _then_ read 'maxVectorSz' from y
        int const y_size = biggest_y_write + maxVectorSz;
        
        return {
            0, // x block size
            y_size, // y block size
            {
                {
                    fft_length,
                    n_partitions_in_vec
                }
            },
            2*fft_length // work size
        };
    }
    
    Latency getImpliedLatency() const {
        Assert(handlesCoefficients());
        return Latency(blockSize-1);
    }

    template<typename F>
    void forEachUsingSameContext(F f) const {
        f(*this);
    }

    int blockSize;
};

template <typename T, typename Tag>
struct FFTConvolutionCRTPSimulation {
    using SetupParam = FFTConvolutionCRTPSetupParam;
    using FPT = T;
    using FFTTag = Tag;
    
    void doSetCoefficientsCount(int64_t szCoeffs) {
        fft_length = get_fft_length(szCoeffs);
        N = fft_length/2;
    }

    double cost_compute_convolution() {
        return fft::RealFBinsCosts<Tag, FPT>::cost_mult_assign(N);
    }
    
    static auto get_fft_length(int n) {
        auto N_nonzero_y = 2 * n;
        return ceil_power_of_two(N_nonzero_y);
    }

    auto get_fft_length() const {
        return fft_length;
    }
    
    auto getBlockSize() const { return N; }

    bool handlesCoefficients() const {
        return N > 0;
    }
    Latency getLatency() const {
        Assert(handlesCoefficients());
        return Latency(N-1);
    }

    int N, fft_length = 0;
};

template<typename T, typename FFTTag>
struct Simulation_<FFTConvolutionCRTPSetupParam, T, FFTTag> {
    using type = FFTConvolutionIntermediateSimulation<FFTConvolutionCRTPSimulation<T, FFTTag>>;
};

struct PartitionnedFFTConvolutionCRTPSetupParam : public Cost {
    static constexpr bool has_subsampling = false;
    static constexpr int nCoefficientsFadeIn = 0;

    PartitionnedFFTConvolutionCRTPSetupParam(int partition_size,
                                             int partition_count)
    : partition_size(partition_size)
    , partition_count(partition_count)
    {}
    
    void adjustWork(int const nTargetCoeffs) {
        partition_count = countPartitions(nTargetCoeffs,
                                          partition_size);
        if(!handlesCoefficients()) {
            setCost(0.);
        }
    }
    
    template<Overlap Mode>
    MinSizeRequirement getMinSizeRequirement(int const maxVectorSz) const
    {
        if(!partition_size) {
            return {0,0,{},0};
        }
        int const n_partitions_in_vec = countPartitions(maxVectorSz, partition_size);
        Assert(n_partitions_in_vec);

        int const fft_length = get_fft_length();
        
        // when we have vectorization, we may need to write in the future of y
        int const yAdditionalWrite = (maxVectorSz > 1) ? (fft_length/2) : 0;
        
        int const y_base_write_size = (maxVectorSz-1) + n_partitions_in_vec * (fft_length/2) + [fft_length](){
            if constexpr(Mode == Overlap::Add) {
                return fft_length/2;
            }
            else {
                return 0;
            }
        }();
        
        int const biggest_y_write = y_base_write_size + yAdditionalWrite;
        
        // during the worst iteration, we can:
        // - write at most 'biggest_y_write' to y
        // - and _then_ read 'maxVectorSz' from y
        int const y_size = biggest_y_write + maxVectorSz;
        
        return {
            0, // x block size
            y_size, // y block size
            {
                {
                    fft_length,
                    partition_count + n_partitions_in_vec - 1
                }
            },
            2*fft_length // work size
        };
    }
    
    bool handlesCoefficients() const {
        return countMaxHandledCoeffs() > 0;
    }
    int countMaxHandledCoeffs() const {
        return partition_count * partition_size;
    }
    Latency getImpliedLatency() const {
        Assert(handlesCoefficients());
        return Latency(partition_size-1);
    }

    int getBiggestScale() const {
        return partition_size;
    }

    int get_fft_length() const {
        return 2 * partition_size;
    }
    int partition_size;
    int partition_count;
    
    bool isValid() const {
        return partition_size > 0;
    }
    
    template<typename F>
    void forEachUsingSameContext(F f) const {
        f(*this);
    }

    void logSubReport(std::ostream & os) const {
        os << partition_count << " partition(s) of size " << partition_size << std::endl;
    }
};

template <typename T, typename Tag>
struct PartitionnedFFTConvolutionCRTPSimulation {
    using SetupParam = PartitionnedFFTConvolutionCRTPSetupParam;
    using FPT = T;
    using FFTTag = Tag;

    void setup(SetupParam const & p) {
      partition_size = p.partition_size;
      assert(partition_size > 0);
      assert(is_power_of_two(partition_size));
    }

    void doSetCoefficientsCount(int64_t count) {
        n_partitions = countPartitions(count,
                                       partition_size);
    }
    
    double cost_compute_convolution() const
    {
        double cost = 0.;
        if(n_partitions) {
            cost += fft::RealFBinsCosts<Tag, FPT>::cost_multiply(partition_size);
            cost += (n_partitions-1) * fft::RealFBinsCosts<Tag, FPT>::cost_multiply_add(partition_size);
        }
        return cost;
    }
    
    auto get_fft_length() const { assert(partition_size > 0); return 2 * partition_size; }
    auto getBlockSize() const { assert(partition_size > 0); return partition_size; }
    
    bool handlesCoefficients() const {
        return n_partitions * partition_size > 0;
    }
    Latency getLatency() const {
        Assert(handlesCoefficients());
        return Latency(partition_size-1);
    }
private:
    int partition_size = -1;
    int n_partitions = 0;
};

template<typename T, typename FFTTag>
struct Simulation_<PartitionnedFFTConvolutionCRTPSetupParam, T, FFTTag> {
    using type = FFTConvolutionIntermediateSimulation<PartitionnedFFTConvolutionCRTPSimulation<T, FFTTag>>;
};

}
