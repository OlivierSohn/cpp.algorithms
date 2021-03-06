

static inline std::vector<double> mkCoefficientsRamp(int sz) {
    std::vector<double> res;
    res.reserve(sz);
    for(int i=0; i<sz; ++i) {
        res.push_back(1.-(i/(double)sz));
    }
    return res;
}

static inline void scaleVec(double coeff, std::vector<double> & v) {
    std::for_each(v.begin(), v.end(), [coeff](auto & val) { val *= coeff; });
}

namespace imajuscule::audio {
template<ReverbType reverbType, int nOut, int nIns, typename ...Args>
void testReverbDirac(Args ...args) {
    using Rev = Reverbs<nOut, reverbType, PolicyOnWorkerTooSlow::Wait>;
    typename Rev::MemResource::type memory;
    Rev rs;

    using WorkCplxFreqs = typename decltype(rs)::WorkCplxFreqs;
    WorkCplxFreqs work;

    using Convolution = typename decltype(rs)::State;
    using Desc = typename Convolution::Desc;

    constexpr int audio_cb_size = 99;
    constexpr int maxVectorSz = maxVectorSizeFromBlockSizeHypothesis(audio_cb_size);
    
    std::vector<double const *> a_inputs;
    a_inputs.resize(nIns);
    std::vector<double *> a_outputs;
    a_outputs.resize(nOut);

    // by default, reverb is inactive.
    {
        std::vector<double> const input = mkDirac<double>(4);
        std::vector<std::vector<double>> output;
        output.resize(nOut);
        for(auto & o : output) {
            o.resize(input.size(), {});
        }
        auto prevOutput = output;
        
        for(int i=0; i<nIns; ++i) {
            a_inputs[i] = input.data();
        }

        for(int i=0; i<nOut; ++i) {
            a_outputs[i] = output[i].data();
        }

        rs.assignWetVectorized(a_inputs.data(),
                               nIns,
                               a_outputs.data(),
                               nOut,
                               input.size(),
                               1);
        ASSERT_EQ(output, prevOutput);
    }
    
    // for 0-length responses, reverb is inactive.
    {
        try {
            if constexpr (Desc::has_subsampling) {
                applyBestParams(rs,
                                memory,
                                1,
                                {a64::vector<double>{}, 1}, work, audio_cb_size, maxVectorSz, 44100., std::cout,
                                ResponseTailSubsampling::HighestAffordableResolution,
                                args...);
            }
            else if constexpr (reverbType == ReverbType::Offline) {
                XFFTsCostsFactors unbiasedXFftCostFactors;
                applyBestParams(rs,memory,1,
                                {a64::vector<double>{}, 1}, work, audio_cb_size, maxVectorSz, 44100., std::cout,
                                unbiasedXFftCostFactors,
                                args...);
            }
            else {
                applyBestParams(rs,memory,1,
                                {a64::vector<double>{}, 1}, work, audio_cb_size, maxVectorSz, 44100., std::cout,
                                args...);
            }
            ASSERT_TRUE(false);
        }
        catch(std::exception const &) {
        }
        
        std::vector<double> const input = mkDirac<double>(4);
        std::vector<std::vector<double>> output;
        output.resize(nOut);
        for(auto & o : output) {
            o.resize(input.size(), {});
        }
        auto prevOutput = output;

        for(int i=0; i<nIns; ++i) {
            a_inputs[i] = input.data();
        }
        
        for(int i=0; i<nOut; ++i) {
            a_outputs[i] = output[i].data();
        }
        
        rs.assignWetVectorized(a_inputs.data(),
                               nIns,
                               a_outputs.data(),
                               nOut,
                               input.size(),
                               maxVectorSz);
        ASSERT_EQ(output, prevOutput);
    }
    
    std::vector<int> sizes = {
        1,
        2,
        3,
        121,
        1476,
        8639
#ifdef NDEBUG
        , 37867
#endif
    };
    for(auto const sz : sizes) {
        /*std::cout << "sz=" << sz << std::endl;
        IndentingOStreambuf ind(std::cout);*/

        auto const coeffs_base = mkCoefficientsTriangle(sz);
        
        for(int crosstalk = 0; crosstalk<2; ++crosstalk) {
            if(!crosstalk) {
                if(nIns != nOut) {
                    continue;
                }
            }
            /*std::cout << "crosstalk=" << crosstalk << std::endl;
            IndentingOStreambuf ind(std::cout);*/

            int const nResponses = crosstalk ? (nIns*nOut) : nOut;
            int const nResponsesPerOut = nResponses / nOut;

            std::vector<a64::vector<double>> all_coeffs;
            all_coeffs.reserve(nResponses);
            for(int i=0; i<nResponses; ++i) {
                all_coeffs.push_back(coeffs_base);
                double const ratio = (i+1)/static_cast<double>(nResponses);
                for(auto & d : all_coeffs.back()) {
                    d *= ratio;
                }
            }
            
            std::vector<a64::vector<double>> expectedOutputs;
            expectedOutputs.resize(nOut);
            for(auto & e : expectedOutputs) {
                e.resize(sz);
            }
            for(int i=0; i<nResponses; ++i) {
                auto & coeffs = all_coeffs[i];
                auto & d = expectedOutputs[crosstalk ? (i/nIns) : i];
                for( int j=0; j<sz; ++j) {
                    d[j] += coeffs[j];
                }
            }
                        
            std::vector<double> const input = mkDirac<double>(sz);
            for(int i=0; i<nIns; ++i) {
                a_inputs[i] = input.data();
            }
            
            std::vector<std::vector<double>> output;
            output.resize(nOut);
            for(auto & o : output) {
                o.resize(input.size(), {});
            }
            for(int i=0; i<nOut; ++i) {
                a_outputs[i] = output[i].data();
            }
            
            for(int rts_i=0;
                rts_i<(Desc::has_subsampling?5:1);
                ++rts_i)
            {
                bool retried = false;
            retry:
                auto const rts = static_cast<ResponseTailSubsampling>(rts_i);
                
                auto const scaleRange = getScaleCountRanges<Convolution>(rts);
                
                {
                    bool res = false;
                    try {
                        if constexpr (Desc::has_subsampling) {
                            applyBestParams(rs, memory, nIns,
                                            all_coeffs, work, audio_cb_size, maxVectorSz, 44100., std::cout,
                                            rts,
                                            args...);
                        }
                        else if constexpr (reverbType == ReverbType::Offline) {
                            XFFTsCostsFactors unbiasedXFftCostFactors;
                            applyBestParams(rs, memory, nIns,
                                            all_coeffs, work, audio_cb_size, maxVectorSz, 44100., std::cout,
                                            unbiasedXFftCostFactors,
                                            args...);
                        }
                        else {
                            applyBestParams(rs, memory, nIns,
                                            all_coeffs, work, audio_cb_size, maxVectorSz, 44100., std::cout,
                                            
                                            args...);
                        }
                        res = true;
                    }
                    catch(std::exception const &e) {
                        if constexpr (!Desc::has_subsampling) {
                            if(rts != ResponseTailSubsampling::FullRes &&
                               rts != ResponseTailSubsampling::HighestAffordableResolution) {
                                // an exception is expected in that case
                                continue;
                            }
                        }
                        LG(INFO, "%s", e.what());
                    }
                    if(scaleRange.getMin() > 1 &&
                       sz <= minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter.toInteger()) {
                        ASSERT_FALSE(res);
                        continue;
                    }
                    if(scaleRange.getMin() > 1 && !res) {
                        continue;
                    }
                    if(!res) { // retry once, to allow debugging
                        if(retried) {
                            ASSERT_TRUE(res);
                        }
                        else {
                            retried = true;
                            goto retry;
                        }
                    }
                }
                
                for(int i=0; i<8; ++i) // retry after a flushToSilence()
                {
                    /*std::cout << "i=" << i << std::endl;
                    IndentingOStreambuf ind(std::cout);*/

                    rs.assignWetVectorized(a_inputs.data(),
                                           nIns,
                                           a_outputs.data(),
                                           nOut,
                                           input.size(),
                                           maxVectorSz);
                    
                    for(int k=0; k<output.size(); ++k)
                    {
                        std::vector<double> diff1;
                        diff1.reserve(output.size());
                        
                        auto & o = output[k];
                        auto const & expected_output = expectedOutputs[k];
                        auto scale = expected_output[0]/o[0];
                        for(int i=0, sz = o.size(); i<sz; i++) {
                            o[i] *= scale;
                            
                            diff1.push_back(std::abs(o[i] - expected_output[i]));
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
                        
                        auto test = [&diff, n_scales, nResponsesPerOut]() {
                            switch(n_scales) {
                                case 1:
                                    if(nResponsesPerOut*0.00001 < diff[0].first) {
                                        return false;
                                    }
                                    break;
                                case 2:
                                    if(nResponsesPerOut*0.05 < diff[0].first) {
                                        return false;
                                    }
                                    break;
                                case 3:
                                    if(nResponsesPerOut*0.08 < diff[0].first) {
                                        return false;
                                    }
                                    break;
                                case 4:
                                    if(nResponsesPerOut*0.2 < diff[0].first) {
                                        return false;
                                    }
                                    break;
                                default:
                                    return false;
                            }
                            return true;
                        };
                        
                        ASSERT_TRUE(test());
                    }
                    
                    // before flushing to silence, we feed ones:
                    
                    {
                        int const nOnes = input.size() / (i+1);
                        Assert(nOnes <= input.size());
                        if(nOnes) {
                            rs.assignWetVectorized(a_inputs.data(),
                                                   nIns,
                                                   a_outputs.data(),
                                                   nOut,
                                                   nOnes,
                                                   maxVectorSz);
                        }
                    }
                    
                    rs.flushToSilence();
                }
            }
        }
    }
    
    // reverb is inactive when disabled
    {
        rs.clear();
        std::vector<double> const input = mkDirac<double>(4);
        std::vector<std::vector<double>> output;
        output.resize(nOut);
        for(auto & o : output) {
            o.resize(input.size(), {});
        }
        auto prevOutput = output;

        for(int i=0; i<nOut; ++i) {
            a_outputs[i] = output[i].data();
        }

        for(int i=0; i<nIns; ++i) {
            a_inputs[i] = input.data();
        }

        rs.assignWetVectorized(a_inputs.data(),
                               nIns,
                               a_outputs.data(),
                               nOut,
                               input.size(),
                               maxVectorSz);
        ASSERT_EQ(output, prevOutput);
    }
}
} // NS imajuscule

template<int nOut, int nIns>
void testReverbDiracNout()
{
    using namespace imajuscule;
    using namespace imajuscule::audio;
    testReverbDirac<ReverbType::Realtime_Asynchronous_Legacy, nOut, nIns>(SimulationPhasing::no_phasing());
    testReverbDirac<ReverbType::Realtime_Synchronous, nOut, nIns>();
    testReverbDirac<ReverbType::Offline, nOut, nIns>();
    testReverbDirac<ReverbType::Realtime_Asynchronous, nOut, nIns>(SimulationPhasing::no_phasing());
    //testReverbDirac<ReverbType::Realtime_Synchronous_Subsampled, nOut, nIns>();
}

TEST(Reverbs, dirac)
{
    // nIns = nOuts = 1
    testReverbDiracNout<1,1>();

    // nIns = 1
    testReverbDiracNout<2,1>();
    testReverbDiracNout<3,1>();
    testReverbDiracNout<4,1>();

    // nOuts = 1
    testReverbDiracNout<1,2>();
    testReverbDiracNout<1,3>();
    testReverbDiracNout<1,4>();

    // nIns = nOuts
    testReverbDiracNout<2,2>();
    testReverbDiracNout<3,3>();
    testReverbDiracNout<4,4>();
}

TEST(Reverbs, reproQueueSizeGarageband) {
    using namespace imajuscule;
    using namespace imajuscule::audio;

    using Rev = Reverbs<2, ReverbType::Realtime_Asynchronous, PolicyOnWorkerTooSlow::Wait>;
    typename Rev::MemResource::type memory;
    Rev rs;
    
    using WorkCplxFreqs = typename decltype(rs)::WorkCplxFreqs;
    WorkCplxFreqs work;
    
    using Convolution = typename decltype(rs)::State;
    
    constexpr int audio_cb_size = 128;
    
    int sz = 450000;
    std::vector<double> const input = mkDirac<double>(sz);
    auto const coeffs = mkCoefficientsTriangle(sz);
    std::vector<a64::vector<double>> vcoeffs;
    vcoeffs.push_back(coeffs);
    vcoeffs.push_back(coeffs);
    
    {
        bool res = false;
        try {
            applyBestParams(rs,
                            memory,
                            1,
                            {vcoeffs},
                            work,
                            audio_cb_size,
                            maxVectorSizeFromBlockSizeHypothesis(audio_cb_size),
                            44100.,
                            std::cout,
                            SimulationPhasing::phasing_with_group_size(2));
            res = true;
        }
        catch(std::exception const &e) {
            LG(ERR, "%s", e.what());
        }
        ASSERT_TRUE(res);
    }
}
