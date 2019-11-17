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
    template<typename FPT, ReverbType reverbType, int nResponses, int nOut, int nIn = nOut>
    void profile() {
        using namespace imajuscule;
        using namespace imajuscule::profiling;
        using namespace std::chrono;
        std::cout << "bench vectorization ["; COUT_TYPE(FPT); std::cout << "] " << toString(reverbType) << " " << nResponses << " response(s) " << nIn << " ins " << nOut << " outs." << std::endl;
        
        std::vector<a64::vector<FPT>> a_inputs;
        a_inputs.resize(nIn);
        constexpr auto inputSz = 1000000;
        for(auto & inputs:a_inputs) {
            inputs.reserve(inputSz);
            for(int i=0; i<inputSz; ++i) {
                inputs.push_back(sin(i));
            }
        }
        constexpr auto outputSz = inputSz;
        std::vector<a64::vector<FPT>> a_outputs;
        a_outputs.resize(nOut);
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
        
        int n_audio_frames_per_cb = 128;
        double frame_rate = 44100.;
        Reverbs<nOut, reverbType> rs;
        ResponseStructure structure;
        {
            std::stringstream os;
            rs.setConvolutionReverbIR(deinterlaced_buffers,
                                      n_audio_frames_per_cb,
                                      frame_rate,
                                      ResponseTailSubsampling::HighestAffordableResolution,
                                      os,
                                      structure);
        }
        
        std::this_thread::sleep_for(std::chrono::seconds(5));
        
        std::cout << "  szCoeffs=" << nCoeffs << " szAudioCb=" << n_audio_frames_per_cb << std::endl;
        
        int vectorLength = 0;
        std::array<FPT*, nIn> input_buffers;
        for(int i=0; i<nIn; ++i) {
            input_buffers[i] = a_inputs[i].data();
        }
        std::array<FPT*, nOut> output_buffers;
        for(int i=0; i<nOut; ++i) {
            output_buffers[i] = a_outputs[i].data();
        }
        
        AudioHostSimulator simu{
            frame_rate,
            n_audio_frames_per_cb,
            std::move(a_inputs),
            std::move(a_outputs)
        };
        
        constexpr int nPeriods = 12;
        // one "period" of computations has a length of structure.totalSizePadded.
        int nTotalFrames = nPeriods * structure.totalSizePadded;
        auto p = simu.simulate([&rs, &nTotalFrames, input_buffers, &output_buffers, vectorLength]
                               (std::vector<a64::vector<float>> const & inputs,
                                std::vector<a64::vector<float>> & outputs) {
            assert(inputs.size() == outputs.size());
            assert(inputs.size());
            int const n_frames = inputs[0].size();
            
            if(vectorLength==0) {
                rs.assignWet(input_buffers.data(),
                             nIn,
                             output_buffers.data(),
                             nOut,
                             n_frames
                             );
            }
            else {
                rs.assignWetVectorized(input_buffers.data(),
                                       nIn,
                                       output_buffers.data(),
                                       nOut,
                                       n_frames,
                                       vectorLength
                                       );
            }
            
            
            nTotalFrames -= n_frames;
            return nTotalFrames > 0;
        });
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
    
    //profile<float, ReverbType::Realtime_Synchronous, 4,2>(); // 72 % of multiply add, 10% dotptr, 3% forward
    //profile<float, ReverbType::Realtime_Asynchronous, 4,2>();
    //profile<float, ReverbType::Offline, 1,1>();
    profile<float, ReverbType::Realtime_Asynchronous, 1,1>();
    return 0;
}
