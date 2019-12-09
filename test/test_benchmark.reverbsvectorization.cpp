
namespace imajuscule::vecto {
    class NullBuffer : public std::streambuf
    {
    public:
      int overflow(int c) { return c; }
    };
    
    /*
     The optimal vector length could depend on:
     - the type of the input / output buffers (float or double)
     - the number of inputs, the number of outputs
     - the length of the reverb (number of coefficients)
     - the audio callback size

     In this test, input data is aligned on cache line boundaries.
     */
    template<typename FPT, int nResponses, int nOut, int nIn = nOut>
    void test() {
        using namespace imajuscule;
        using namespace imajuscule::profiling;
        using namespace std::chrono;
        std::cout << "bench vectorization ["; COUT_TYPE(FPT); std::cout << "] " << nResponses << " response(s) " << nIn << " ins " << nOut << " outs." << std::endl;
        
        std::array<a64::vector<FPT>, nIn> a_inputs;
        constexpr auto inputSz = 1000000;
        for(auto & inputs:a_inputs) {
            inputs.reserve(inputSz);
            for(int i=0; i<inputSz; ++i) {
                inputs.push_back(sin(i));
            }
        }
        constexpr auto outputSz = inputSz;
        std::array<a64::vector<FPT>, nOut> a_outputs;
        for(auto & outputs:a_outputs) {
            outputs.resize(outputSz);
        }
        
        std::vector<int> coeffsSz;
        for(int powCoeffs=15; powCoeffs<19; ++powCoeffs) {
            coeffsSz.push_back(pow2(powCoeffs+1)-1);
        }
        
        std::vector<int> blocksSz;
        for(int powBlockSz=3; powBlockSz<12; ++powBlockSz) {
            blocksSz.push_back(pow2(powBlockSz));
        }
        
        for(auto nCoeffs : coeffsSz) {
            std::vector<a64::vector<double>> db;
            db.resize(nResponses);
            for(auto & coeffs:db) {
                coeffs.reserve(nCoeffs);
                for(int i=0; i<nCoeffs; ++i) {
                    coeffs.push_back(cos(static_cast<double>(i)/100.));
                }
            }
            DeinterlacedBuffers<double> deinterlaced_buffers(std::move(db));
            
            for(auto audio_cb_size : blocksSz) {
                Reverbs<nOut> rs;
                ResponseStructure structure;
                try {
                    NullBuffer null_buffer;
                    std::ostream null_stream(&null_buffer);
                    rs.setConvolutionReverbIR(1,
                                              deinterlaced_buffers,
                                              audio_cb_size,
                                              44100.,
                                              ResponseTailSubsampling::HighestAffordableResolution,
                                              null_stream,
                                              structure);
                }
                catch(std::exception const &e) {
                    LG(INFO, "skip : %s", e.what());
                    continue;
                }
                
                auto inc = [](int & v) {
                    if(v==0) {
                        v=1;
                    }
                    else {
                        v*=2;
                    }
                };
                std::cout << "  szCoeffs=" << nCoeffs << " szAudioCb=" << audio_cb_size << std::endl;

                std::map<int, double> durationsByVectorLength;
                auto plot = [](auto const & m) {
                    std::vector<double> values;
                    values.reserve(m.size());
                    for(auto const & [key, value] : m) {
                        values.push_back(value);
                    }
                    {
                        auto p = StringPlot(20,values.size());
                        p.drawLog(values, '+');
                        p.log();
                    }
                    {
                        auto p = StringPlot(20,values.size());
                        p.draw(values, '*');
                        p.log();
                    }
                };
                // skip the first iteration to warmup the cache
                for(int test = 0; test <= 3; ++test) {
                    for(int vectorLength = 0; vectorLength <= audio_cb_size; inc(vectorLength)) {
                        std::array<FPT*, nIn> input_buffers;
                        for(int i=0; i<nIn; ++i) {
                            input_buffers[i] = a_inputs[i].data();
                        }
                        std::array<FPT*, nOut> output_buffers;
                        for(int i=0; i<nOut; ++i) {
                            output_buffers[i] = a_outputs[i].data();
                        }
                        
                        constexpr int nPeriods = 12;
                        auto dt = measure_one<high_resolution_clock>([&]() {
                            for(auto per=0; per<nPeriods; ++per) {
                                // one "period" of computations has a length of structure.totalSizePadded.
                                for(int i=0; i<structure.totalSizePadded; i += audio_cb_size) {
                                    // emulate what would be done by the audio callback
                                    int numFrames = std::min(structure.totalSizePadded-i, audio_cb_size);
                                    if(vectorLength==0) {
                                        rs.assignWet(input_buffers.data(),
                                                     nIn,
                                                     output_buffers.data(),
                                                     nOut,
                                                     numFrames
                                                     );
                                    }
                                    else {
                                        rs.assignWetVectorized(input_buffers.data(),
                                                               nIn,
                                                               output_buffers.data(),
                                                               nOut,
                                                               numFrames,
                                                               vectorLength
                                                               );
                                    }
                                }
                            }
                        }).count();
                        if(test) {
                            durationsByVectorLength[vectorLength] = dt;
                        }
                    }
                    if(test) {
                        plot(durationsByVectorLength);
                    }
                }
            }
            
        }
    }
}

/*
 Tuning of the 'vectorLength' argument in
 Reverbs::assignWetVectorized
 Reverbs::addWetInputZeroVectorized
 */
TEST(BenchmarkReverbVectorization, vectorization) {
    using namespace imajuscule::vecto;
    // true stereo
    test<double, 4, 2>();
    // stereo
    test<double, 2, 2>();
    // mono
    test<double, 1, 1>();
}
