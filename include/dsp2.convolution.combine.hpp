
namespace imajuscule {

template<typename T, typename FFTTag>
using AlgoOptimizedFIRFilter =
/**/Convolution<
/**/  AlgoSplitConvolution<
/**/    AlgoFIRFilter<T, FFTTag>,
/**/    AlgoCustomScaleConvolution<
/**/      AlgoFFTConvolutionIntermediate<
/**/        AlgoPartitionnedFFTConvolutionCRTP<T, FFTTag>>>>>;

template<typename T, typename FFTTag>
using AlgoZeroLatencyScaledFineGrainedPartitionnedConvolution =
/**/Convolution<
/**/  AlgoSplitConvolution<
/**/    typename AlgoOptimizedFIRFilter<T, FFTTag>::Algo,
/**/    AlgoFinegrainedPartitionnedFFTConvolution<T, FFTTag>>>;

template<typename T, typename FFTTag, PolicyOnWorkerTooSlow OnWorkerTooSlow>
using AlgoZeroLatencyScaledAsyncConvolution =
/**/Convolution<
/**/  AlgoSplitConvolution<
/**/    typename AlgoOptimizedFIRFilter<T, FFTTag>::Algo,
/**/    AlgoAsyncCPUConvolution<
/**/      AlgoCustomScaleConvolution<
/**/        AlgoFFTConvolutionIntermediate<
/**/          AlgoPartitionnedFFTConvolutionCRTP<T, FFTTag>>>,
/**/      OnWorkerTooSlow>>>;

template<typename T, typename FFTTag, PolicyOnWorkerTooSlow OnWorkerTooSlow>
using AlgoZeroLatencyScaledAsyncConvolutionOptimized =
/**/Convolution<
/**/  AlgoSplitConvolution<
/**/    typename AlgoZeroLatencyScaledFineGrainedPartitionnedConvolution<T, FFTTag>::Algo,
/**/    AlgoAsyncCPUConvolution<
/**/      AlgoCustomScaleConvolution<
/**/        AlgoFFTConvolutionIntermediate<
/**/          AlgoPartitionnedFFTConvolutionCRTP<T, FFTTag>>>,
/**/      OnWorkerTooSlow>>>;

}
