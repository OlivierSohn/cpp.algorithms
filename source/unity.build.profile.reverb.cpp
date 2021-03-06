/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

#include "unity.build.cpp"

#include "test_utils.h"
#include "../test/test_utils.cpp"

namespace imajuscule::bench::vecto {
    CsvFile csv_file("/Users/Olivier/profile_reverb.csv");
    
    class NullBuffer : public std::streambuf
    {
    public:
        int overflow(int c) { return c; }
    };
    
    struct audio_dropout : public std::exception {
    };

    /*
     The optimal vector length could depend on:
     - the type of the input / output buffers (float or double)
     - the number of inputs, the number of outputs
     - the length of the reverb (number of coefficients)
     - the audio callback size
     
     In this test, input data is aligned on cache line boundaries.
     */
    template<ReverbType reverbType, int nOut, typename Tag, typename F, typename ...Args>
    void profile(int const n_sources,
                 int nResponses,
                 int const nCoeffs,
                 int const n_audio_frames_per_cb,
                 std::optional<int> const vector_size,
                 double frame_rate,
                 F f,
                 Args&&... args) {
        using namespace imajuscule;
        using namespace imajuscule::profiling;
        using namespace std::chrono;
        
        using Rev = Reverbs<
        nOut,
        reverbType,
        PolicyOnWorkerTooSlow::PermanentlySwitchToDry, /* because AudioHostSimulator uses real time simulation*/
        Tag
        >;
        
        using FPT = typename Rev::FPT;
        
        using WorkCplxFreqs = typename Rev::WorkCplxFreqs;
        
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
        csv_file.push(vector_size?*vector_size : 0);
        
        Rev rs;
        WorkCplxFreqs work;
        
        std::vector<a64::vector<double>> db;
        db.resize(nResponses);
        for(auto & coeffs:db) {
            coeffs.reserve(nCoeffs);
            for(int i=0; i<nCoeffs; ++i) {
                coeffs.push_back(cos(static_cast<double>(i)/100.));
            }
        }
        DeinterlacedBuffers<double> deinterlaced_buffers(std::move(db));
        auto memory = std::make_unique<typename Rev::MemResource::type>();
        
        applyBestParams(rs,
                        *memory,
                        n_sources,
                        deinterlaced_buffers,
                        work,
                        n_audio_frames_per_cb,
                        vector_size ? *vector_size : 1,
                        frame_rate,
                        std::cout,
                        args...);
        
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
        
        std::vector<FPT const *> input_buffers;
        input_buffers.reserve(a_inputs.size());
        
        std::vector<FPT*> output_buffers;
        output_buffers.reserve(a_outputs.size());
        
        AudioHostSimulator<FPT> simu{
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
            int nTotalFrames = nPeriods * rs.getBiggestScale();
            auto p = simu.simulate([&rs, &nTotalFrames, &input_buffers, &output_buffers, vector_size, &sync_durations]
                                   (std::vector<a64::vector<FPT>> const & inputs,
                                    std::vector<a64::vector<FPT>> & outputs,
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
                        
                        if(!vector_size) {
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
                                                   *vector_size
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
                if(conv.b.hasStepErrors()) {
                    throw audio_dropout();
                }
            });
        }
        
        std::cout << "simu wall:" << CpuDuration::fromDuration(simu_wall_duration) << std::endl;
        csv_file.push(CpuDuration::fromDuration(simu_wall_duration));
        
        {
            CpuDuration totalCpu;
            totalCpu = std::accumulate(sync_durations.begin(),
                                       sync_durations.end(),
                                       totalCpu);
            if constexpr (isAsync(reverbType)) {
                rs.foreachConvReverb([&totalCpu](auto const & conv) {
                    auto & async_durations = conv.b.async_durations;
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

std::ostream & operator << (std::ostream & o, SystemLoadModelization s);


template<ReverbType ReverbType, int nOuts, typename F, typename ...Args>
void getFramesInQueue(SystemLoadModelization slm,
                      int const system_load,
                      int const n_input_channels,
                      int const response_sz,
                      int const n_frames_cb,
                      std::optional<int> const vector_sz,
                      F f_foreach_conv,
                      Args&&... args) {
    using namespace imajuscule::bench::vecto;
    using namespace imajuscule;

    IndentingOStreambuf indent{std::cout};
    std::cout << "'" << slm << "' " << system_load << ", ";
    
    double const frame_rate = 44100. * ((slm == SystemLoadModelization::HigherFramerate) ? system_load : 1.);
    int const n_responses = n_input_channels * ((slm == SystemLoadModelization::HigherFramerate) ? 1 : system_load);

    for_each(std::tuple<accelerate::Tag, accelerate2::Tag>{},
             [&](auto tag)
    {
        csv_file.newline();
        csv_file.push(slm);
        csv_file.push(system_load);
        {
            typedef void(*FT)(decltype(tag));
            auto tag_str = type_to_string<FT>()(FT());
            csv_file.push(justifyLeft(12,
                                       // remove leading  'imajuscule::'
                                       // remove trailing '::Tag'
                                       std::string(tag_str.begin() + 12,
                                                   tag_str.end() - 5)));
        }
        profile<ReverbType, nOuts, decltype(tag)>(n_responses/nOuts, // n_sources
                                                  n_responses,
                                                  response_sz,
                                                  n_frames_cb,
                                                  vector_sz,
                                                  frame_rate,
                                                  f_foreach_conv,
                                                  args...);
    });
}

void testCostsReadWrites() {
    using namespace imajuscule;
    
    // mesurable à partir de 10000 écritures.
    //
    // pour chaque taille, la premiere execution est plus lente
    // (donc malloc nous renvoie probablement les memes blocs memoire
    //  qui finissent par être dans les caches)
    //
    for(int i=1; i<10000000; i*=10) {
        a64::vector<double> v;
        v.reserve(i);
        for(int k=0; k<2; ++k) {
            auto duration = imajuscule::profiling::measure_thread_cpu_one([&v, i](){
                for(int j=0; j<i; ++j) {
                    v[j] = j;
                }
            });
            std::cout << "write " << i << ":" << duration.count() << std::endl;
        }
    }
}


void testCostsReadWrites2() {
    using namespace imajuscule;
    
    // mesurable à partir de 10000 écritures.
    //
    // pour chaque taille, la premiere execution est plus lente
    // (donc malloc nous renvoie probablement les memes blocs memoire
    //  qui finissent par être dans les caches)
    //
    int constexpr Max = 10000000;
    for(int i=1; i<Max; i*=10) {
        std::vector<a64::vector<double>> v;
        v.resize(Max/i);
        for(auto & vv : v) {
            vv.resize(i);
        }
        for(int k=0; k<2; ++k) {
            auto duration = imajuscule::profiling::measure_thread_cpu_one([&v, i](){
                for(auto & vv : v) {
                    for(int j=0; j<i; ++j) {
                        vv[j] = j;
                    }
                }
            });
            std::cout << "write " << i << ":" << duration.count() / static_cast<double>(Max/i) << std::endl;
            /*
             double sum = 0;
             for(auto & vv : v) {
             for(int j=0; j<i; ++j) {
             sum += vv[j];
             }
             }
             std::cout << sum << std::endl;
             */
        }
    }
}


void testCostsReadWrites3() {
    using namespace imajuscule;
    using namespace imajuscule::profiling;
    
    mersenne<SEEDED::Yes>().seed(0);
    
    int constexpr Max = 10000000;
    using T = double;
    a64::vector<T> v;
    int constexpr nTPerCacheline = 64/sizeof(T);
    v.resize(Max);
    
    double sideEffect = 0.;
    
    std::vector<int> indexes;
    indexes.reserve(Max);
    
    for(int i=1; i<Max; i>100?(i*=10):++i) {
        // different tests sit at least 4 cachelines appart
        int offsetBetweenTests = i + 4*nTPerCacheline;
        // make it a multiple of 512
        if(offsetBetweenTests < 512) {
            offsetBetweenTests = 512;
        }
        offsetBetweenTests = (offsetBetweenTests/512) * 512;
        std::cout << "offsetBetweenTests " << offsetBetweenTests << std::endl;
        indexes.clear();
        for(int index=0;index+i <= Max; index += offsetBetweenTests) {
            indexes.push_back(index);
        }
        std::shuffle(indexes.begin(), indexes.end(), mersenne<SEEDED::Yes>());
        
        for(int k=0; k<2; ++k)
        {
            polluteCache(sideEffect);
            loadInCache(indexes.begin(), indexes.end(), sideEffect);
            
            auto duration = measure_thread_cpu_one([&v, &indexes, &sideEffect, i](){
                for(auto index : indexes) {
                    for(int j=0; j<i; ++j) {
                        v[index + j] += sideEffect+j;
                    }
                }
                sideEffect += v[0];
            });
            std::cout << "write " << i << ":" << duration.count() / static_cast<double>(indexes.size()) << " over " << indexes.size() << std::endl;
            /*
             double sum = 0;
             for(auto & vv : v) {
             for(int j=0; j<i; ++j) {
             sum += vv[j];
             }
             }
             std::cout << sum << std::endl;
             */
        }
    }
    
    std::cout << "sideEffect:" << sideEffect << std::endl;
}

namespace imajuscule::profiling {
    struct StatBucket {
        static constexpr float width = 10.f;
        
        StatBucket(int64_t value):
        r(value, value)
        {}
        
        bool tryFeed(int64_t value) {
            int64_t expectedMin = std::min(r.getMin(), value);
            int64_t expectedMax = std::max(r.getMax(), value);
            int64_t expectedSpan =  expectedMax - expectedMin;
            if(expectedSpan > width*expectedMin) {
                return false;
            }
            r.extend(value);
            ++count;
            return true;
        }
        
        auto getRange() const {
            return r;
        }
        
        auto getCount() const {
            return count;
        }
    private:
        range<int64_t> r;
        int64_t count = 1;
    };
    std::ostream & operator << (std::ostream & os, const StatBucket& s);
    
    struct Stats {
        static constexpr int nBuckets = 5;
        
        Stats() {
            v.reserve(nBuckets);
        }
        
        void start() {
            wallTimes.first = std::chrono::steady_clock::now();
            rusage ru;
            int res = my_getrusage(RUSAGE_THREAD,&ru);
            if(res) {
                throw std::runtime_error("cannot measure time");
            }
            cpuTimes.first = {ru};
        }
        void stop() {
            wallTimes.second = std::chrono::steady_clock::now();
            rusage ru;
            int res = my_getrusage(RUSAGE_THREAD,&ru);
            if(res) {
                throw std::runtime_error("cannot measure time");
            }
            cpuTimes.second = {ru};
        }
        
        template<typename Duration>
        bool feed(Duration d) {
            auto i = std::chrono::duration_cast<std::chrono::nanoseconds>(d).count();
            for(auto & s:v) {
                if(s.tryFeed(i)) {
                    return true;
                }
            }
            if(v.size() < nBuckets) {
                v.emplace_back(i);
                return true;
            }
            return false;
        }
        
        auto getSortedBuckets() const {
            auto v2 = v;
            std::sort(v2.begin(), v2.end(), [](auto const & a, auto const & b){ return a.getRange().getMin() < b.getRange().getMin(); });
            return v2;
        }
        CpuDuration getCpuDuration() const {
            return cpuTimes.second - cpuTimes.first;
        }
        std::chrono::steady_clock::duration getWallDuration() const {
            return wallTimes.second - wallTimes.first;
        }
        
    private:
        std::vector<StatBucket> v;
        std::pair<std::chrono::steady_clock::time_point, std::chrono::steady_clock::time_point> wallTimes;
        std::pair<CpuDuration, CpuDuration> cpuTimes;
    };
    std::ostream & operator << (std::ostream & os, const Stats& s);
}

void testScheduling(imajuscule::profiling::Stats & s, int64_t niterations, std::atomic<bool> const & start)
{
    while(!start) {
        std::this_thread::yield();
    }
    s.start();
    {
        imajuscule::profiling::IntervalsTimer<std::chrono::steady_clock> t;
        for(int64_t i=0;i<niterations; ++i) {
            s.feed(t.elapsedSinceLast());
        }
    }
    s.stop();
}

void testScheduling()
{
    int64_t constexpr nThreads = 2;
    int64_t constexpr nIterations = 1000000000 / nThreads;
    std::vector<imajuscule::profiling::Stats> stats;
    stats.resize(nThreads);
    {
        std::vector<std::thread> threads;
        threads.reserve(nThreads);
        std::atomic<bool> start = false;
        for(auto & s : stats) {
            threads.emplace_back([&s, &start](){testScheduling(s, nIterations, start);});
        }
        start = true;
        for(auto & t : threads) {
            t.join();
        }
    }
    
    for(auto const & s : stats) {
        std::cout << s << std::endl;
    }
    
}

struct NopEval {
    void prepare(int thread_index) {
    }
    
    void eval() {
    }
};

void testAllocationFactors2() {
    using namespace imajuscule::fft;
    using namespace imajuscule::profiling;
    for(int nThreads = 1; nThreads < 20; nThreads *= 2) {
        for(int nMicroSecs = 1000; nMicroSecs < 1000000000; nMicroSecs *= 10) {
            std::cout << "nThreads " << nThreads << "\t us " << nMicroSecs << "\t";
            {
                MaxWallTimeIncrementEval::BythreadMaxIncrements maxIncrements;
                maxIncrements.resize(nThreads);
                auto allocs = computeAllocationFactors(nThreads,
                                                       std::chrono::microseconds(nMicroSecs),
                                                       MaxWallTimeIncrementEval(maxIncrements));
                
                std::vector<std::chrono::steady_clock::duration> maxIncrementsValues;
                maxIncrementsValues.resize(maxIncrements.size());
                std::transform(maxIncrements.begin(),
                               maxIncrements.end(),
                               maxIncrementsValues.begin(),
                               [](auto & o) {
                    if(!o.value) {
                        throw std::runtime_error("no max increment");
                    }
                    return *o.value;
                });
                
                std::cout << std::chrono::duration_cast<std::chrono::microseconds>(*std::max_element(maxIncrementsValues.begin(),
                                                                                                     maxIncrementsValues.end())).count() / 1000000.;
                std::cout << "\t " << *std::min_element(allocs.begin(),
                                                        allocs.end());
            }
            std::cout << std::endl;
        }
    }
}


void testAllocationFactors() {
    std::vector<int> nmicros{
        100000,
        1000000,
        10000000,
    };
    for(int nThreads = 1; nThreads < 100; ++nThreads) {
        for(auto nMicroSecs : nmicros) {
            auto allocs = computeAllocationFactors(nThreads,
                                                   std::chrono::microseconds(nMicroSecs),
                                                   NopEval());
            
            auto totalAllocation = std::accumulate(allocs.begin(),
                                                   allocs.end(),
                                                   0.);
            std::cout << nThreads;
            std::cout << " " << totalAllocation;
            std::cout << " " << *std::min_element(allocs.begin(), allocs.end());
            std::cout << " " << *std::max_element(allocs.begin(), allocs.end());
            /*
             for(auto a : allocs) {
             std::cout << " " << a;
             }*/
            std::cout << std::endl;
        }
    }
}

template<typename T, typename Tag>
void printConvolutionCosts2() {
    using namespace imajuscule::fft;
    
    auto compute = [](auto & f){
        std::vector<double> costs, costsPerSample;
        constexpr int64_t startSz =
        //2;
        // start where we have a full cache lime of data
        cache_line_n_bytes / sizeof(T);
        for(int64_t sz = startSz;
            sz<5000000;
            sz *= 2)
        {
            double const totalCost = f(sz);
            double const costPerSample = totalCost / sz;
            costs.push_back(totalCost);
            costsPerSample.push_back(costPerSample);
            std::cout << sz << "\t" << costPerSample << "\t" << totalCost << std::endl;
        }
        {
            StringPlot plot(20, costs.size());
            plot.drawLog(costs, '*');
            plot.log();
        }
        {
            StringPlot plot(20, costsPerSample.size());
            plot.draw(costsPerSample, '+');
            plot.log();
        }
    };
    
    std::cout << std::endl << "fft forward" << std::endl << std::endl;
    compute(AlgoCosts<Tag, T>::cost_fft_forward);
    
    std::cout << std::endl << "fft inverse" << std::endl << std::endl;
    compute(AlgoCosts<Tag, T>::cost_fft_inverse);
    
    std::cout << std::endl << "sig add scalar multiply" << std::endl << std::endl;
    compute(RealSignalCosts<Tag, T>::cost_add_scalar_multiply);
    
    std::cout << std::endl << "sig add assign" << std::endl << std::endl;
    compute(RealSignalCosts<Tag, T>::cost_add_assign);
    
    std::cout << std::endl << "sig copy" << std::endl << std::endl;
    compute(RealSignalCosts<Tag, T>::cost_copy);
    
    std::cout << std::endl << "sig dotpr" << std::endl << std::endl;
    compute(RealSignalCosts<Tag, T>::cost_dotpr);
    
    std::cout << std::endl << "freq mult assign" << std::endl << std::endl;
    compute(RealFBinsCosts<Tag, T>::cost_mult_assign);
}

template<typename T>
void printConvolutionCosts() {
    for_each(fft::Tags, [](auto t) {
        using Tag = decltype(t);
        COUT_TYPE(Tag);
        std::cout << std::endl;
        printConvolutionCosts2<T, Tag>();
    });
}

/*
 Clang bug? We need these to instantiate some undefined templates
 */
void explicitInstantiations() {
    using namespace imajuscule;
    if(false) {
        fft::RealSignal_<accelerate::Tag, double>::add_assign_output<double>(nullptr, nullptr, 0);
        fft::RealSignal_<accelerate2::Tag, double>::add_assign_output<double>(nullptr, nullptr, 0);
    }
}

int main(int argc, const char * argv[]) {
    using namespace std;
    using namespace imajuscule;
    using namespace imajuscule::bench::vecto;
    explicitInstantiations();
/*
    testAllocationFactors2();
    return 0;
    testCostsReadWrites3();
    return 0;
    printConvolutionCosts<double>();
    return 0;
    testScheduling();
    return 0;
*/
    if(argc != 1) {
        cerr << "0 argument is needed, " << argc-1 << " given" << endl;
        printUsage();
        throw;
    }
    
    csv_file.push("load_model");
    csv_file.push("load");
    csv_file.push("tag");
    csv_file.push("n_responses");
    csv_file.push("n_sources");
    csv_file.push("n_outs");
    csv_file.push("sampling_rate");
    csv_file.push("reverb_type");
    csv_file.push("response_n_frames");
    csv_file.push("cb_n_frames");
    csv_file.push("vector_size");
    csv_file.push("n_frames_in_queue"); // optionnel
    csv_file.push("simu_wall");
    csv_file.push("cpu_total");
    csv_file.push("cpu_sync_max_per_cb");
    
    std::vector<int> cb_sz{1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};
    std::reverse(cb_sz.begin(), cb_sz.end());
    std::vector<int> responses_sizes{1, 10, 100, 1000, 10000, 100000, 1000000};
    std::reverse(responses_sizes.begin(), responses_sizes.end());

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
                for(int vector_sz = 0; vector_sz <= 2*n_frames_cb; )
                {
                    std::optional<int> vsz;
                    if(vector_sz) {
                        vsz = vector_sz;
                    }
                    //for(int i=0; i<2; ++i)
                    {
                        //bool phasing = static_cast<bool>(i);
                        bool phasing = true;
                        
                        auto const sim_phasing = phasing? SimulationPhasing::phasing_with_group_size(n_channels): SimulationPhasing::no_phasing();
                        std::cout << "phasing " << sim_phasing << std::endl;
                        
                        try {
                            getFramesInQueue<ReverbType::Realtime_Asynchronous, nOuts>(slm,
                                                                                       system_load,
                                                                                       n_channels,
                                                                                       response_sz,
                                                                                       n_frames_cb,
                                                                                       vsz,
                                                                                       [](auto const & rs){
                                auto & l = rs.getAlgo().b;
                                auto nFramesInQueue = l.getQueueSize() * l.getSubmissionPeriod();
                                csv_file.push(justifyRight(5, std::to_string(nFramesInQueue)));
                                return true;
                            },
                                                                                       sim_phasing);
                        }
                        catch (std::exception const & e) {
                            std::cout << "!!!!!!!!!!!!!!!!!!!!" << std::endl;
                            std::cout << e.what() << std::endl;
                            std::cout << "!!!!!!!!!!!!!!!!!!!!" << std::endl;
                        }
                        
                        /*
                         
                         try {
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
                        getFramesInQueue<ReverbType::Realtime_Synchronous, nOuts>(slm,
                                                                                  system_load,
                                                                                  n_channels,
                                                                                  response_sz,
                                                                                  n_frames_cb,
                                                                                  vsz,
                                                                                  [](auto const &){
                            csv_file.push("  N/A");
                            return true; });
                    }
                    catch (std::exception const & e) {
                        std::cout << "!!!!!!!!!!!!!!!!!!!!" << std::endl;
                        std::cout << e.what() << std::endl;
                        std::cout << "!!!!!!!!!!!!!!!!!!!!" << std::endl;
                    }
                    /*
                     try {
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
                    if(vector_sz) {
                        vector_sz *= 2;
                    }
                    else {
                        vector_sz = 1;
                    }
                }
            }
        }
    }
    
    return 0;
}

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

            namespace imajuscule::profiling {
            
            std::ostream & operator << (std::ostream & os, const StatBucket& s){
            os << "[" << s.getRange().getMin() << " " << s.getRange().getMax() << "] " << s.getCount();
            return os;
            }
            std::ostream & operator << (std::ostream & os, const Stats& s) {
            for(auto & b:s.getSortedBuckets()) {
            os << b << std::endl;
            }
            os << s.getCpuDuration().count() << " us cpu" << std::endl;
            os << std::chrono::duration_cast<std::chrono::microseconds>(s.getWallDuration()).count() << " us wall" << std::endl;
            return os;
            }
            }
