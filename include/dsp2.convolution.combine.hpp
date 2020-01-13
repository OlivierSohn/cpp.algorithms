
namespace imajuscule {

template<typename T, typename FFTTag>
using AlgoOptimizedFIRFilter =
/**/Convolution<
/**/  AlgoSplitConvolution<
/**/    AlgoFIRFilter<T, FFTTag>,
/**/    AlgoCustomScaleConvolution<
/**/      AlgoFFTConvolutionIntermediate<
/**/        AlgoPartitionnedFFTConvolutionCRTP<T, FFTTag>>>>>;

}
