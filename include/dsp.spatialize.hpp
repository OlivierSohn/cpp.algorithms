/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {
template<typename C>
int countScales(C & rev) {
    auto & lateHandler = rev.getB();
    {
        auto & inner = lateHandler.getB().getInner().getInner();
        if(inner.isZero()) {
            return 1;
        }
    }
    {
        auto & inner = lateHandler.getB().getInner().getInner().getB().getInner().getInner();
        if(inner.isZero()) {
            return 2;
        }
    }
    {
        auto & inner = lateHandler.getB().getInner().getInner().getB().getInner().getInner().getB().getInner().getInner();
        if(inner.isZero()) {
            return 3;
        }
    }
    static_assert(4==nMaxScales);
    return 4;
}

template<typename C>
bool scalesAreValid(int n_scales, C & rev) {
    auto & lateHandler = rev.getB();
    for(int i=1; i<n_scales; ++i) {
        // the top-most will be stepped 'base_phase' times,
        // then each scale after that will be stepped by a quarter grain size.
        // phases are cumulative, so stepping a scale also steps subsequent scales.
        switch(i) {
            case 1:
            {
                auto & inner = lateHandler.getB().getInner().getInner();
                if(!inner.isValid()) {
                    return false;
                }
                if(inner.isZero()) {
                    return false;
                }
                break;
            }
            case 2:
            {
                auto & inner = lateHandler.getB().getInner().getInner().getB().getInner().getInner();
                if(!inner.isValid()) {
                    return false;
                }
                if(inner.isZero()) {
                    return false;
                }
                break;
            }
            case 3:
            {
                auto & inner = lateHandler.getB().getInner().getInner().getB().getInner().getInner().getB().getInner().getInner();
                if(!inner.isValid()) {
                    return false;
                }
                if(inner.isZero()) {
                    return false;
                }
                break;
            }
            default:
                throw std::logic_error("out of bound");
        }
    }
    
    return true;
}


namespace audio {
/*
 * Sources are spatialized, i.e. each source location provides one impulse response per ear
 *
 * TODO perform parameter optimization globally: since we know we will do
 * ** several ** convolutions per sample, we can relax optimization constraints
 * somewhat.
 */
template<int N_EARS, typename Convolution>
struct Spatializer {
    
    static constexpr auto nEars = N_EARS;
    
    using T = typename Convolution::FPT;
    using FPT = T;
    
    void logReport(double sampleRate, std::ostream & os)
    {
        os << "States:" << std::endl;
        IndentingOStreambuf in(os);
        
        int i=0;
        foreachConvReverb([&os, &i](auto const & r){
            ++i;
            os << "Convolution " << i << ":" << std::endl;
            IndentingOStreambuf indent(os);
            r.logComputeState(os);
        });
    }
    
    double getEpsilon() const {
        return epsilonOfNaiveSummation(convs) / nEars;
    }
    
    void clear() {
        convs.clear();
        nConvolutionsPerEar = 0;
    }
    void flushToSilence() {
        for(auto & c : convs) {
            c->flushToSilence();
        }
    }
    
    bool empty() const {
        return convs.empty();
    }
    
    bool isValid() const {
        return !empty() && convs[0]->isValid();
    }
     
    int countScales() {
        if (convs.empty()) {
            return 0;
        }
        if constexpr(Convolution::has_subsampling) {
            return imajuscule::countScales(*convs[0]);
        }
        else {
            return 1;
        }
    }
    
    bool handlesCoefficients() const {
        if(convs.empty()) {
            return false;
        }
        return convs[0]->handlesCoefficients();
    }
    
    Latency getLatency() const {
        Assert(handlesCoefficients());
        return convs.empty() ? Latency(0) : convs[0]->getLatency();
    }
    
    template<typename SetupP>
    void setSources(int n_sources,
                    std::vector<a64::vector<T>> const & deinterlaced_coeffs,
                    SetupP const & setup) {
        convs.reserve(deinterlaced_coeffs.size());
        for(auto & coeffs : deinterlaced_coeffs) {
            convs.push_back(std::make_unique<Convolution>());
            auto & c = convs.back();
            c->setup(setup);
            c->setCoefficients(std::move(coeffs));
        }
        int const nConvolutionsPerSource = convs.size() / n_sources;
        nConvolutionsPerEar = convs.size() / nEars;
        Assert((n_sources * nConvolutionsPerSource) == convs.size());
        Assert((nEars * nConvolutionsPerEar) == convs.size());
    }
    
    // for tests
    void setMaxMultiplicationGroupLength() {
        for(auto & c : convs) {
            c->setMultiplicationGroupLength(c->getHighestValidMultiplicationsGroupSize());
        }
    }
    
    void setMultiplicationGroupLength(int i) {
        for(auto & c : convs) {
            c->setMultiplicationGroupLength(i);
        }
    }
    
    template<typename FPT2>
    bool assignWetVectorized(FPT2 const * const * const input_buffers,
                             int nInputBuffers,
                             FPT2 ** output_buffers,
                             int nOutputBuffers,
                             int nFramesToCompute,
                             int vectorLength)
    {
        bool success = true;

        int i_in = 0;
        auto itConv = convs.begin();
        
        Assert(nOutputBuffers == nEars);
        for(int i_out=0; i_out < nOutputBuffers; ++i_out) {
            FPT2 * out = output_buffers[i_out];
            
            bool assign = true;
            for(auto end = i_in+nConvolutionsPerEar;
                i_in < end;
                ++i_in, ++itConv) {
                Assert(i_in < nInputBuffers);
                
                auto & c = *itConv;
                FPT2 const * const in = input_buffers[i_in];
                for(int i=0; i<nFramesToCompute; i += vectorLength) {
                    if(assign) {
                        c->stepAssignVectorized(in + i,
                                                out + i,
                                                std::min(vectorLength, nFramesToCompute-i));
                    }
                    else {
                        c->stepAddVectorized(in + i,
                                             out + i,
                                             std::min(vectorLength, nFramesToCompute-i));
                    }
                }
                assign = false;
                if constexpr (Convolution::step_can_error) {
                    success = !c->hasStepErrors() && success;
                }
            }
            if(i_in == nInputBuffers) {
                i_in = 0;
            }
        }
        return success;
    }
    
    template<typename FPT2>
    bool addWetInputZeroVectorized(FPT2 ** output_buffers,
                                   int nOutputBuffers,
                                   int nFramesToCompute,
                                   int vectorLength) {
        bool success = true;
        
        Assert(nOutputBuffers == nEars);
        for(int i_out=0; i_out < nOutputBuffers; ++i_out) {
            FPT2 * out = output_buffers[i_out];
            for(int i=0; i<nFramesToCompute; i += vectorLength) {
                for(int i_in = 0; i_in < nConvolutionsPerEar; ++i_in) {
                    auto & c = convs[nConvolutionsPerEar*i_out + i_in];
                    c->stepAddInputZeroVectorized(out+i,
                                                  std::min(vectorLength, nFramesToCompute-i));
                    if constexpr (Convolution::step_can_error) {
                        success = !c->hasStepErrors() && success;
                    }
                }
            }
        }
        return success;
    }
    
    void dephaseComputations() {
        int n = 0;
        int const total = convs.size();
        for(auto & c : convs) {
            dephase(total, n, *c);
            ++n;
        }
    }
    
    template<typename F>
    void foreachConvReverb(F f) const {
        for(auto const & c : convs) {
            f(*c);
        }
    }
private:
    int nConvolutionsPerEar = 0;
    std::vector<std::unique_ptr<Convolution>> convs;
};
}

}
