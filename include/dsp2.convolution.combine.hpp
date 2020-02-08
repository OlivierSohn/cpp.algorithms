
namespace imajuscule {

template<typename T, template<typename> typename Allocator, typename FFTTag>
using AlgoOptimizedFIRFilter =
/**/AlgoSplitConvolution<
/**/  AlgoFIRFilter<T, Allocator, FFTTag>,
/**/  AlgoCustomScaleConvolution<
/**/    AlgoFFTConvolutionIntermediate<
/**/      AlgoPartitionnedFFTConvolutionCRTP<T, Allocator, FFTTag>>>>;

template<typename T, template<typename> typename Allocator, typename FFTTag>
using AlgoZeroLatencyScaledFineGrainedPartitionnedConvolution =
/**/AlgoSplitConvolution<
/**/  AlgoOptimizedFIRFilter<T, Allocator, FFTTag>,
/**/  AlgoFinegrainedPartitionnedFFTConvolution<T, Allocator, FFTTag>>;

template<typename T, template<typename> typename Allocator, typename FFTTag, PolicyOnWorkerTooSlow OnWorkerTooSlow>
using AlgoZeroLatencyScaledAsyncConvolution =
/**/AlgoSplitConvolution<
/**/  AlgoOptimizedFIRFilter<T, Allocator, FFTTag>,
/**/  AlgoAsyncCPUConvolution<
/**/    AlgoCustomScaleConvolution<
/**/      AlgoFFTConvolutionIntermediate<
/**/        AlgoPartitionnedFFTConvolutionCRTP<T, Allocator, FFTTag>>>,
/**/    OnWorkerTooSlow>>;

template<typename T, template<typename> typename Allocator, typename FFTTag, PolicyOnWorkerTooSlow OnWorkerTooSlow>
using AlgoZeroLatencyScaledAsyncConvolutionOptimized =
/**/AlgoSplitConvolution<
/**/  AlgoZeroLatencyScaledFineGrainedPartitionnedConvolution<T, Allocator, FFTTag>,
/**/  AlgoAsyncCPUConvolution<
/**/    AlgoCustomScaleConvolution<
/**/      AlgoFFTConvolutionIntermediate<
/**/        AlgoPartitionnedFFTConvolutionCRTP<T, Allocator, FFTTag>>>,
/**/    OnWorkerTooSlow>>;



enum class AudioProcessing {
    Callback, // here, we want 0 latency and the smallest worst case cost per audio callback.
    Offline // here, we want the smallest averaged cost per sample.
};



namespace detail {
template<typename T, AudioProcessing P>
struct OptimalFilter_;

template<typename T>
struct OptimalFilter_<T, AudioProcessing::Callback> {
    // TODO reevaluate once we know when async is better than sync
    using type = AlgoZeroLatencyScaledFineGrainedPartitionnedConvolution<T, a64::Alloc, fft::Fastest>;
};

template<typename T>
struct OptimalFilter_<T, AudioProcessing::Offline> {
    using type = AlgoOptimizedFIRFilter<T, a64::Alloc, fft::Fastest>;
};

}

template<typename T, AudioProcessing P>
using ZeroLatencyFilter = typename detail::OptimalFilter_<T,P>::type;

}
