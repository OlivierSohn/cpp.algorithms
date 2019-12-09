/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

#include "unity.build.cpp"


namespace imajuscule::bench::vecto {
    struct CsvFile {
        CsvFile(std::string path) {
            file.open (path);
        }
        ~CsvFile() {
            file.close();
        }
        
        template<typename T>
        void push(T && v) {
            using namespace imajuscule::profiling;
            if(n) {
                file << ", ";
            }
            file << v;
            ++n;
        }
        void newline() {
            file << std::endl;
            file.flush();
            n = 0;
        }
    private:
        std::ofstream file;
        int n = 0;
    };
    
    CsvFile csv_file("/Users/Olivier/profile_reverb.csv");

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
    template<typename FPT, ReverbType reverbType, int nOut, typename F, typename ...Args>
    void profile(int const n_sources,
                 int nResponses,
                 int const nCoeffs,
                 int const n_audio_frames_per_cb,
                 double frame_rate,
                 F f,
                 Args&&... args) {
        using namespace imajuscule;
        using namespace imajuscule::profiling;
        using namespace std::chrono;
        
        if(nResponses != n_sources * nOut) {
            throw std::logic_error("inconsitent n_sources / nOut");
        }

        std::cout <<
        nResponses << " response(s), " <<
        n_sources << " sources, " <<
        nOut << " outs, " <<
        frame_rate << "Hz, " <<
        toJustifiedString(reverbType) << ", " <<
        nCoeffs << " response frames, " <<
        n_audio_frames_per_cb << " frames per cb, ";
        COUT_TYPE(FPT); std::cout << " buffers." <<
        std::endl;

        csv_file.push(nResponses);
        csv_file.push(n_sources);
        csv_file.push(nOut);
        csv_file.push(frame_rate);
        csv_file.push(toJustifiedString(reverbType));
        csv_file.push(justifyRight(12,std::to_string(nCoeffs)));
        csv_file.push(n_audio_frames_per_cb);

        std::vector<a64::vector<double>> db;
        db.resize(nResponses);
        for(auto & coeffs:db) {
            coeffs.reserve(nCoeffs);
            for(int i=0; i<nCoeffs; ++i) {
                coeffs.push_back(cos(static_cast<double>(i)/100.));
            }
        }
        DeinterlacedBuffers<double> deinterlaced_buffers(std::move(db));
        
        Reverbs<nOut, reverbType> rs;
        ResponseStructure structure;
        {
            rs.setConvolutionReverbIR(n_sources,
                                      deinterlaced_buffers,
                                      n_audio_frames_per_cb,
                                      frame_rate,
                                      std::cout,
                                      structure,
                                      args...);
        }
        if(!f(rs)) {
            return;
        }
                
        std::vector<a64::vector<FPT>> a_inputs;
        a_inputs.resize(n_sources);
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
        
        //std::this_thread::sleep_for(std::chrono::seconds(5));
                
        int vectorLength = 0;
        std::vector<FPT const *> input_buffers;
        input_buffers.reserve(a_inputs.size());
        
        std::vector<FPT*> output_buffers;
        output_buffers.reserve(a_outputs.size());

        AudioHostSimulator simu{
            frame_rate,
            n_audio_frames_per_cb,
            std::move(a_inputs),
            std::move(a_outputs)
        };

        std::chrono::steady_clock::duration simu_wall_duration;
        std::vector<profiling::CpuDuration> sync_durations;
        sync_durations.reserve(100000);
        {
            WallTimer<std::chrono::steady_clock> t(simu_wall_duration);
            
            constexpr int nPeriods = 2;
            // one "period" of computations has a length of structure.totalSizePadded.
            int nTotalFrames = nPeriods * structure.totalSizePadded;
            auto p = simu.simulate([&rs, &nTotalFrames, &input_buffers, &output_buffers, vectorLength, &sync_durations]
                                   (std::vector<a64::vector<float>> const & inputs,
                                    std::vector<a64::vector<float>> & outputs,
                                    int const frameStart,
                                    int const n_frames) {
                assert(inputs.size());
                input_buffers.clear();
                for(auto & i : inputs) {
                    input_buffers.push_back(i.data() + frameStart);
                }
                output_buffers.clear();
                for(auto & o : outputs) {
                    output_buffers.push_back(o.data() + frameStart);
                }
                
                {
                    std::optional<profiling::CpuDuration> threadCPUDuration;
                    {
                        profiling::ThreadCPUTimer tt(threadCPUDuration);
                        
                        if(vectorLength==0) {
                            rs.assignWet(input_buffers.data(),
                                         input_buffers.size(),
                                         output_buffers.data(),
                                         output_buffers.size(),
                                         n_frames
                                         );
                        }
                        else {
                            rs.assignWetVectorized(input_buffers.data(),
                                                   input_buffers.size(),
                                                   output_buffers.data(),
                                                   output_buffers.size(),
                                                   n_frames,
                                                   vectorLength
                                                   );
                        }
                    }
                    if(!threadCPUDuration) {
                        throw std::runtime_error("no time");
                    }
                    sync_durations.emplace_back(*threadCPUDuration);
                }
                
                nTotalFrames -= n_frames;
                return nTotalFrames > 0;
            });
        }

        auto printTimes = [](auto const & c, std::optional<int> i = {})
        {
            {
                constexpr auto w = 8;
                if(i) {
                    std::stringstream ss;
                    ss << std::setfill(' ') << std::setw(w) << *i;
                    std::cout << ss.str();
                }
                else {
                    std::cout << std::string(w, ' ');
                }
            }
            std::cout << "\t" << c.getGlobal();
            //std::cout << "\t" << c.getUser() << "\t" << c.getKernel();
            std::cout << std::endl;
        };

        std::cout << "simu wall:" << CpuDuration::fromDuration(simu_wall_duration) << std::endl;
        csv_file.push(CpuDuration::fromDuration(simu_wall_duration));
        
        {
            CpuDuration totalCpu;
            totalCpu = std::accumulate(sync_durations.begin(),
                                       sync_durations.end(),
                                       totalCpu);
            if constexpr (isAsync(reverbType)) {
                rs.foreachConvReverb([&totalCpu](auto const & conv) {
                    auto & async_durations = conv.getB().async_durations;
                    totalCpu = std::accumulate(async_durations.begin(),
                                               async_durations.end(),
                                               totalCpu);
                });
            }
            std::cout << "total cpu:"; printTimes(totalCpu);
            csv_file.push(totalCpu.getGlobal());
        }

        {
            auto const & durations = sync_durations;
            {
                auto M = std::max_element(durations.begin(), durations.end());
                if(M == durations.end()) {
                    std::cout << "no sync max:" << std::endl;
                    csv_file.push("");
                }
                else {
                    std::cout << "sync max :"; printTimes(*M
                                                          //, std::distance(durations.begin(), M)
                                                          );
                    csv_file.push(M->getGlobal());
                }
            }
/*
            int i=0;
            for(auto c : durations)
            {
                printTimes(c, i);
                ++i;
            }
  //*/
            
        }

        if constexpr (isAsync(reverbType))
        {
            rs.foreachConvReverb([&printTimes](auto const & conv)
                                 {
                /*
                 auto & durations = conv.getB().async_durations;
                 {
                 auto M = std::max_element(durations.begin(), durations.end());
                 if(M == durations.end()) {
                 std::cout << "no async max" << std::endl;
                 }
                 else {
                 std::cout << "async max:"; printTimes(*M, std::distance(durations.begin(), M));
                 }
                 }
                 */
                
                /*
                 int i=0;
                 for(auto c : durations)
                 {
                 printTimes(c, i);
                 ++i;
                 }*/
                {
                    int nRTWait = conv.getB().countErrorsWorkerTooSlow();
                    if(nRTWait) {
                        std::cout << "count rt waits : " << nRTWait << std::endl;
                        std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;
                        throw std::runtime_error("count rt waits > 0");
                    }
                }
                
            });
        }
    }
}


void printUsage() {
    using namespace std;
    cout << "- profile reverb usage : " << endl;
}

enum class SystemLoadModelization {
    HigherFramerate,
    MoreTracks
};

std::ostream & operator << (std::ostream & o, SystemLoadModelization s)
{
    switch(s)
    {
        case SystemLoadModelization::HigherFramerate :
            o << "HigherFramerate";
            break;
        case SystemLoadModelization::MoreTracks :
            o << "MoreTracks";
            break;
        default:
            throw std::runtime_error("unhandled SystemLoadModelization");
    }
    return o;
}

template<ReverbType ReverbType, int nOuts, typename F, typename ...Args>
void getFramesInQueue(SystemLoadModelization slm,
                     int const system_load,
                     int const n_input_channels,
                      int const response_sz,
                      int const n_frames_cb,
                     F f_foreach_conv,
                     Args&&... args) {
    using namespace imajuscule::bench::vecto;

    IndentingOStreambuf indent{std::cout};
    std::cout << "'" << slm << "' " << system_load << ", ";
    csv_file.push(slm);
    csv_file.push(system_load);
            
    double const frame_rate = 44100. * ((slm == SystemLoadModelization::HigherFramerate) ? system_load : 1.);
    int const n_responses = n_input_channels * ((slm == SystemLoadModelization::HigherFramerate) ? 1 : system_load);
    profile<float, ReverbType, nOuts>(// n_sources:
                                      n_responses/nOuts,
                                      n_responses,
                                      response_sz,
                                      n_frames_cb,
                                      frame_rate,
                                      f_foreach_conv,
                                      args...);
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
    
    csv_file.push("load_model");
    csv_file.push("load");
    csv_file.push("n_responses");
    csv_file.push("n_sources");
    csv_file.push("n_outs");
    csv_file.push("sampling_rate");
    csv_file.push("reverb_type");
    csv_file.push("response_n_frames");
    csv_file.push("cb_n_frames");
    csv_file.push("simu_wall");
    csv_file.push("cpu_total");
    csv_file.push("cpu_sync_max_per_cb");
    csv_file.push("n_frames_in_queue"); // optionel

    std::vector<int> cb_sz{1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};
    std::vector<int> responses_sizes{1, 10, 100, 1000, 10000, 100000, 1000000};
    
    std::vector<std::pair<SystemLoadModelization, int>> loads;
    {
        constexpr int max_system_load = 1;//5
        for(int sl = 0; sl<max_system_load; ++sl)
        {
            int const system_load = pow2(sl);
            loads.emplace_back(SystemLoadModelization::MoreTracks, system_load);
            if(system_load > 1) {
                loads.emplace_back(SystemLoadModelization::HigherFramerate, system_load);
            }
        }
    }

    constexpr int nOuts = 2;
    int const n_channels = 2; // we _must_ use more than one channel to witness the effect of phasing.

    for(auto n_frames_cb : cb_sz)
    {
        for(auto [slm, system_load] : loads)
        {
            for(auto response_sz: responses_sizes)
            {
                //for(int i=0; i<2; ++i)
                {
                    //bool phasing = static_cast<bool>(i);
                    bool phasing = true;
                    
                    auto const sim_phasing = phasing? SimulationPhasing::phasing_with_group_size(n_channels): SimulationPhasing::no_phasing();
                    std::cout << "phasing " << sim_phasing << std::endl;

                    try {
                        csv_file.newline();
                        std::optional<int> nFramesInQueue;
                        getFramesInQueue<ReverbType::Realtime_Asynchronous, nOuts>(slm,
                                                                                    system_load,
                                                                                    n_channels,
                                                                                   response_sz,
                                                                                   n_frames_cb,
                                                                                    [&nFramesInQueue](auto const & rs){
                            rs.foreachConvReverb([&nFramesInQueue](auto const & conv){
                                auto & l = conv.getB();
                                auto n = l.getQueueSize() * l.getSubmissionPeriod();
                                if(nFramesInQueue) {
                                    if(*nFramesInQueue != n) {
                                        throw std::logic_error("*nFramesInQueue != n");
                                    }
                                }
                                else {
                                    nFramesInQueue = n;
                                }
                            });
                            return true;
                        },
                                                                                   sim_phasing);
                        csv_file.push(*nFramesInQueue);
                    }
                    catch (std::exception const & e) {
                        std::cout << "!!!!!!!!!!!!!!!!!!!!" << std::endl;
                        std::cout << e.what() << std::endl;
                        std::cout << "!!!!!!!!!!!!!!!!!!!!" << std::endl;
                    }
            
                    /*

                    try {
                        csv_file.newline();
                        std::optional<int> nFramesInQueue;
                        getFramesInQueue<ReverbType::Realtime_Asynchronous_Legacy, nOuts>(slm,
                                                                                        system_load,
                                                                                        n_channels,
                                                                                       response_sz,
                                                                                       n_frames_cb,
                                                                                        [&nFramesInQueue](auto const & rs){
                            rs.foreachConvReverb([&nFramesInQueue](auto const & conv){
                                auto & l = conv.getB();
                                auto n = l.getQueueSize() * l.getSubmissionPeriod();
                                if(nFramesInQueue) {
                                    if(*nFramesInQueue != n) {
                                        throw std::logic_error("*nFramesInQueue != n");
                                    }
                                }
                                else {
                                    nFramesInQueue = n;
                                }
                            });
                            return true;
                        },
                                                                                   sim_phasing);
                        csv_file.push(*nFramesInQueue);
                    }
                    catch (std::exception const & e) {
                        std::cout << "!!!!!!!!!!!!!!!!!!!!" << std::endl;
                        std::cout << e.what() << std::endl;
                        std::cout << "!!!!!!!!!!!!!!!!!!!!" << std::endl;
                    }
                    */
                }
                try {
                    csv_file.newline();
                    getFramesInQueue<ReverbType::Realtime_Synchronous, nOuts>(slm,
                                                                              system_load,
                                                                              n_channels,
                                                                              response_sz,
                                                                              n_frames_cb,
                                                                              [](auto const &){ return true; });
                }
                catch (std::exception const & e) {
                    std::cout << "!!!!!!!!!!!!!!!!!!!!" << std::endl;
                    std::cout << e.what() << std::endl;
                    std::cout << "!!!!!!!!!!!!!!!!!!!!" << std::endl;
                }
                /*
                try {
                    csv_file.newline();
                    getFramesInQueue<ReverbType::Offline, nOuts>(slm,
                                                                 system_load,
                                                                 n_channels,
                                                                 response_sz,
                                                                 n_frames_cb,
                                                                 [](auto const &){ return true; });
                }
                catch (std::exception const & e) {
                    std::cout << "!!!!!!!!!!!!!!!!!!!!" << std::endl;
                    std::cout << e.what() << std::endl;
                    std::cout << "!!!!!!!!!!!!!!!!!!!!" << std::endl;
                }
                 */
            }
        }
    }

    return 0;
}
