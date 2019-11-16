

static inline std::vector<double> mkCoefficientsRamp(int sz) {
    std::vector<double> res;
    res.reserve(sz);
    for(int i=0; i<sz; ++i) {
        res.push_back(1.-(i/(double)sz));
    }
    return res;
}

// the coefficients should tend to 0, to mimic real responses
// (else, with scaling, the end of the dirac's response will be wrong)
static inline std::vector<double> mkCoefficientsTriangle(int sz) {
    std::vector<double> res;
    res.reserve(sz);
    constexpr int period = 100;
    for(int i=0; i<sz; ++i) {
        auto j = i%(2*period);
        if(j < period) {
            res.push_back(1.-(j/(double)period));
        }
        else {
            res.push_back((j-period)/(double)period);
        }
    }
    
    constexpr int fadeout = 40;
    if(sz > fadeout) {
        for(int i=0; i<fadeout; ++i) {
            res[res.size()-1-i] *= (i+1) / static_cast<double>(fadeout);
        }
    }
    return res;
}

static inline void scaleVec(double coeff, std::vector<double> & v) {
    std::for_each(v.begin(), v.end(), [coeff](auto & val) { val *= coeff; });
}

namespace imajuscule {
template<ReverbType reverbType>
void testReverbDirac() {
    
    Reverbs<1, reverbType> rs;
    
    constexpr int audio_cb_size = 99;
    
    // by default, reverb is inactive.
    {
        std::vector<double> const input = mkDirac<double>(4);
        std::vector<double> output;
        output.resize(input.size(), {});
        auto prevOutput = output;
        double const * const a_inputs[1] = {input.data()};
        double * a_outputs[1] = {output.data()};
        rs.assignWetVectorized(a_inputs,
                               1,
                               a_outputs,
                               1,
                               input.size(),
                               input.size());
        ASSERT_EQ(output, prevOutput);
    }
    
    // for 0-length responses, reverb is inactive.
    {
        try {
            ResponseStructure structure;
            rs.setConvolutionReverbIR({{}, 1}, audio_cb_size, 44100., ResponseTailSubsampling::HighestAffordableResolution,
                                      std::cout, structure);
            ASSERT_TRUE(false);
        }
        catch(std::exception const &) {
        }
        std::vector<double> const input = mkDirac<double>(4);
        std::vector<double> output;
        output.resize(input.size(), {});
        auto prevOutput = output;
        double const * const a_inputs[1] = {input.data()};
        double * a_outputs[1] = {output.data()};
        rs.assignWetVectorized(a_inputs,
                               1,
                               a_outputs,
                               1,
                               input.size(),
                               input.size());
        ASSERT_EQ(output, prevOutput);
    }
    
    std::vector<int> sizes = {1,2,3,121,1476,37860,385752,2957213};
    for(auto const sz : sizes) {
        std::vector<double> const input = mkDirac<double>(sz);
        auto const coeffs = mkCoefficientsTriangle(sz);
        for(int rts_i=0; rts_i<5; ++rts_i) {
            using Convolution = typename decltype(rs)::ConvolutionReverb;
            LG(INFO, "%d, %d,", sz, rts_i);
            auto const rts = static_cast<ResponseTailSubsampling>(rts_i);

            auto const scaleRange = getScaleCountRanges<Convolution>(rts);
            
            {
                bool res = false;
                ResponseStructure structure;
                try {
                    rs.setConvolutionReverbIR({coeffs, 1}, audio_cb_size, 44100., rts,
                                              std::cout, structure);
                    res = true;
                }
                catch(std::exception const &e) {
                    if constexpr (!Convolution::has_subsampling) {
                        if(rts != ResponseTailSubsampling::FullRes &&
                           rts != ResponseTailSubsampling::HighestAffordableResolution) {
                            // an exception is expected in that case
                            continue;
                        }
                    }
                    LG(INFO, "%s", e.what());
                }
                if(scaleRange.getMin() > 1 && coeffs.size() <= minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter) {
                    ASSERT_FALSE(res);
                    continue;
                }
                if(scaleRange.getMin() > 1 && !res) {
                    continue;
                }
                if(!res) {
                    ASSERT_TRUE(res);
                }
            }
            std::vector<double> output;
            output.resize(input.size());
            std::fill(output.begin(), output.end(), double{});

            std::vector<double> diff1;
            diff1.reserve(output.size());
            
            {
                double const * const a_inputs[1] = {input.data()};
                double * a_outputs[1] = {output.data()};
                rs.assignWetVectorized(a_inputs,
                                       1,
                                       a_outputs,
                                       1,
                                       1,
                                       1);
            }
            auto scale = coeffs[0]/output[0];
            
            {
                double const * const a_inputs[1] = {input.data()+1};
                double * a_outputs[1] = {output.data()+1};
                rs.assignWetVectorized(a_inputs,
                                       1,
                                       a_outputs,
                                       1,
                                       input.size()-1,
                                       input.size());
            }

            for(int i=0, sz = output.size(); i<sz; i++) {
                output[i] *= scale;
                
                diff1.push_back(std::abs(output[i] - coeffs[i]));
            }
            
            std::vector<std::pair<double, int>> diff;
            diff.reserve(diff1.size());
            for(int i=0; i<diff1.size(); ++i) {
                diff.emplace_back(diff1[i], i);
            }
            std::sort(diff.begin(), diff.end(), std::greater<>());
            int const n_scales = rs.countScales();
            
            ASSERT_TRUE(scaleRange.contains(n_scales));
            
            //int const totalFades = (pow2(n_scales-1)-1) * scaleFadeSz::inSmallerUnits;
            static_assert(4 == nMaxScales);
            switch(n_scales) {
                case 1:
                    ASSERT_GT(0.00001, diff[0].first);
                    break;
                case 2:
                    ASSERT_GT(0.05, diff[0].first);
                    break;
                case 3:
                    ASSERT_GT(0.08, diff[0].first);
                    break;
                case 4:
                    ASSERT_GT(0.2, diff[0].first);
                    break;
                default:
                    ASSERT_TRUE(false);
                    break;
            }
        }
    }
    
    // reverb is inactive when disabled
    {
        rs.disable();
        std::vector<double> const input = mkDirac<double>(4);
        std::vector<double> output;
        output.resize(input.size(), {});
        auto prevOutput = output;
        double const * const a_inputs[1] = {input.data()};
        double * a_outputs[1] = {output.data()};
        rs.assignWetVectorized(a_inputs,
                               1,
                               a_outputs,
                               1,
                               input.size(),
                               input.size());
        ASSERT_EQ(output, prevOutput);
    }
}
} // NS imajuscule

TEST(Reverbs, dirac) {
    using namespace imajuscule;
    testReverbDirac<ReverbType::Offline>();
    testReverbDirac<ReverbType::Realtime_Synchronous>();
    testReverbDirac<ReverbType::Realtime_Asynchronous>();
}
