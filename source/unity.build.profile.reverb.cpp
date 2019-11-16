/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

#include "unity.build.cpp"


namespace imajuscule::bench::vecto {
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
        
        
        int nCoeffs = 500 * 1000;
        std::vector<a64::vector<double>> db;
        db.resize(nResponses);
        for(auto & coeffs:db) {
            coeffs.reserve(nCoeffs);
            for(int i=0; i<nCoeffs; ++i) {
                coeffs.push_back(cos(static_cast<double>(i)/100.));
            }
        }
        DeinterlacedBuffers<double> deinterlaced_buffers(std::move(db));
        
        int audio_cb_size = 128;
        Reverbs<nOut> rs;
        ResponseStructure structure;
        {
            NullBuffer null_buffer;
            std::ostream null_stream(&null_buffer);
            rs.setConvolutionReverbIR(deinterlaced_buffers,
                                      audio_cb_size,
                                      44100.,
                                      ResponseTailSubsampling::HighestAffordableResolution,
                                      null_stream,
                                      structure);
        }
        
        std::cout << "  szCoeffs=" << nCoeffs << " szAudioCb=" << audio_cb_size << std::endl;
        
        int vectorLength = 0;
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
    }
}


void printUsage() {
    using namespace std;
    cout << "- profile reverb usage : " << endl;
}

int main(int argc, const char * argv[]) {
    using namespace std;
    using namespace imajuscule;
    using namespace imajuscule::bench::vecto;
    
    if(argc != 1) {
        cerr << "0 argument is needed, " << argc-1 << " given" << endl;
        printUsage();
        throw;
    }
    
    test<float,4,2>();
    return 0;
}
