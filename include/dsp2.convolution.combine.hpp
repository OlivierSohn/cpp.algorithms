
namespace imajuscule {

template<typename T, template<typename> typename Allocator, typename FFTTag>
using AlgoOptimizedFIRFilter =
/**/Convolution<
/**/  AlgoSplitConvolution<
/**/    AlgoFIRFilter<T, Allocator, FFTTag>,
/**/    AlgoCustomScaleConvolution<
/**/      AlgoFFTConvolutionIntermediate<
/**/        AlgoPartitionnedFFTConvolutionCRTP<T, Allocator, FFTTag>>>>>;

template<typename T, template<typename> typename Allocator, typename FFTTag>
using AlgoZeroLatencyScaledFineGrainedPartitionnedConvolution =
/**/Convolution<
/**/  AlgoSplitConvolution<
/**/    typename AlgoOptimizedFIRFilter<T, Allocator, FFTTag>::Algo,
/**/    AlgoFinegrainedPartitionnedFFTConvolution<T, Allocator, FFTTag>>>;

template<typename T, template<typename> typename Allocator, typename FFTTag, PolicyOnWorkerTooSlow OnWorkerTooSlow>
using AlgoZeroLatencyScaledAsyncConvolution =
/**/Convolution<
/**/  AlgoSplitConvolution<
/**/    typename AlgoOptimizedFIRFilter<T, Allocator, FFTTag>::Algo,
/**/    AlgoAsyncCPUConvolution<
/**/      AlgoCustomScaleConvolution<
/**/        AlgoFFTConvolutionIntermediate<
/**/          AlgoPartitionnedFFTConvolutionCRTP<T, Allocator, FFTTag>>>,
/**/      OnWorkerTooSlow>>>;

template<typename T, template<typename> typename Allocator, typename FFTTag, PolicyOnWorkerTooSlow OnWorkerTooSlow>
using AlgoZeroLatencyScaledAsyncConvolutionOptimized =
/**/Convolution<
/**/  AlgoSplitConvolution<
/**/    typename AlgoZeroLatencyScaledFineGrainedPartitionnedConvolution<T, Allocator, FFTTag>::Algo,
/**/    AlgoAsyncCPUConvolution<
/**/      AlgoCustomScaleConvolution<
/**/        AlgoFFTConvolutionIntermediate<
/**/          AlgoPartitionnedFFTConvolutionCRTP<T, Allocator, FFTTag>>>,
/**/      OnWorkerTooSlow>>>;

}
