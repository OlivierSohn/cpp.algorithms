

namespace imajuscule {
namespace testdsp2conv {

template<typename T>
struct ConvolutionTraits {
    static constexpr bool supportsOddCountOfCoefficients = true;
    static constexpr bool evenIndexesAreApproximated = false;
    
    template <typename F>
    static void adaptCoefficients(F & v) {
    }
};

template<LatencySemantic L, typename T>
struct ConvolutionTraits<SubSampled<L, T>> {
    static constexpr bool supportsOddCountOfCoefficients = false;
    static constexpr bool evenIndexesAreApproximated = true;
    
};

#ifdef NDEBUG
constexpr auto end_index = 16;
#else
constexpr auto end_index = 15;
#endif

template<typename T>
a64::vector<T> makeCoefficients(int coeffs_index) {
    switch(coeffs_index) {
        case -1: return {};
        case 0: return {{ +1. }};
        case 1: return {{ -1. }};
        case 2: return {{ .9, }};
        case 3: return {{ .9,.8 }};
        case 4: return {{ .9,.8,.7 }};
        case 5: return {{ .9,.8,.7,.6 }};
        case 6: return {{ .9,.8,.7,.6,.5 }};
        case 7: return {{ .9,.8,.7,.6,.5,.4 }};
        case 8: return {{ .9,.8,.7,.6,.5,.4,.3 }};
        case 9: return {{ .9,.8,.7,.6,.5,.4,.3,.2 }};
        case 10: return {{ .9,.8,.7,.6,.5,.4,.3,.2,.1 }};
        case 11: return {{ .9,.8,.7,.6,.5,.4,.3,.2,.1,.05 }};
        case 12: return {{ .9,.8,.7,.6,.5,.4,.3,.2,.1,.05,.025 }};
        case 13: return {{ .9,.8,.7,.6,.5,.4,.3,.2,.1,.05,.025,.01 }};
        case 14: return mkTestCoeffs<T>(200);
            // this should exceed the GPU capacity to do the fft in one kernel only,
            // and force partitionning:
        case 15: return mkTestCoeffs<T>(20000);
    }
    throw std::logic_error("coeff index too big");
}

template<typename Convolution, typename Coeffs, typename Input, typename Output>
void testGeneric(Convolution & conv,
                 typename Convolution::SetupParam const & p,
                 Coeffs const & coefficients,
                 Input const & input,
                 Output const & expectedOutput,
                 int const vectorLength)
{
    using T = typename Convolution::FPT;
    using namespace fft;
        
    using Tr = ConvolutionTraits<Convolution>;
    if(!Tr::supportsOddCountOfCoefficients
       && 1 == coefficients.size() % 2) {
        return;
    }
    conv.setupAndSetCoefficients(p, coefficients);

    if(!conv.isValid())
    {
        return;
    }

    std::vector<T> output;
    output.reserve(expectedOutput.size());
    std::vector<T> inputVec;
    inputVec.reserve(expectedOutput.size());
    
    auto const period = conv.getPhasePeriod();
    int const maxPhaseInSample = (period && *period) ? (8 * *period) : 1;
    std::cout << "phases 0 to " << maxPhaseInSample-1 << std::endl;
    // change phases
    for(int phase_in_sample=0; phase_in_sample<maxPhaseInSample; ++phase_in_sample) {
        // fast forward to do less tests
        if(maxPhaseInSample > 100 && phase_in_sample == 10) {
            phase_in_sample = maxPhaseInSample-10;
            continue;
        }
        if(phase_in_sample) {
            float ratio_for_one_sample = 1.f/ *period;
            conv.dephaseByGroupRatio(ratio_for_one_sample * phase_in_sample);
        }
        // repeat after flushToSilence
        for(int i=0; i<8; ++i) {
            
            output.clear();
            inputVec.clear();
            
            auto const eps = conv.getEpsilon();
            int idxStep=0;
            auto step = [&idxStep, &input]() {
                return (idxStep < input.size()) ? input[idxStep] : ((T)0.);
            };
            
            if(conv.handlesCoefficients()) {
                for(; idxStep<conv.getLatency().toInteger(); ++idxStep) {
                    auto res = conv.step(step());
                    if(std::abs(res) > 1000*eps) {
                        LG(INFO,"");
                    }
                    ASSERT_NEAR(0.f, res, 1000*eps); // assumes that no previous signal has been fed
                }
            }
            for(;inputVec.size() != expectedOutput.size();++idxStep) {
                inputVec.push_back(step());
            }
            if(vectorLength) {
                Assert(0);
                /*
                 // initialize with zeros
                 output.resize(inputVec.size(), {});
                 
                 // then add
                 int const nFrames = inputVec.size();
                 for(int i=0; i<nFrames; i += vectorLength) {
                 conv.stepAddVectorized(inputVec.data()+i,
                 output.data()+i,
                 std::min(vectorLength, nFrames-i));
                 }*/
            }
            else {
                for(T i : inputVec) {
                    output.push_back(conv.step(i));
                }
            }
            
            using Tr = ConvolutionTraits<Convolution>;
            for(auto j=0; j<output.size(); ++j) {
                if(Tr::evenIndexesAreApproximated && (0 == j%2)) {
                    continue;
                }
                if(!areNear(expectedOutput[j], output[j], 1000*eps)) { // uses relative error
                    std::cout << std::endl << "... "; COUT_TYPE(Convolution);
                    std::cout << std::endl << "coefficient size : " << coefficients.size() << std::endl;
                    ASSERT_NEAR(expectedOutput[j], output[j], 1000*eps); // doesn't use relative error, only absolute.
                }
            }
            
            for(int i=0; i<10; ++i) {
                conv.step(1.f);
            }
            conv.flushToSilence();
        }
    }
}

template<typename Convolution, typename Coeffs, typename Input, typename Output>
void testGeneric(Convolution & conv,
                 typename Convolution::SetupParam const & p,
                 Coeffs const & coefficients, Input const & input, Output const & expectedOutput) {
    int sz = expectedOutput.size();
    
    std::set<int> vectorLengths;
    
    vectorLengths.insert(0);// no vectorization
    /*
     vectorLengths.insert(1);// minimal vectorization
     vectorLengths.insert(2);
     vectorLengths.insert(3);
     vectorLengths.insert(sz/5);
     vectorLengths.insert(sz/4);
     vectorLengths.insert(sz/4-1);
     vectorLengths.insert(sz/4+1);
     vectorLengths.insert(sz/2);
     vectorLengths.insert(sz);
     vectorLengths.insert(2*sz);*/
    
    
    for(auto vectorLength: vectorLengths) {
        if(vectorLength < 0) {
            continue;
        }
        testGeneric(conv, p, coefficients, input, expectedOutput, vectorLength);
    }
}

// we take 'coefficients' by value because we modify it in the function
template<typename Convolution, typename Coeffs>
void test(Convolution & conv,
          typename Convolution::SetupParam const & p,
          Coeffs const coefficients)
{
    using T = typename Convolution::FPT;
    using Tr = ConvolutionTraits<Convolution>;
    
    // test using a dirac
    {
        auto output = coefficients;
        if(output.empty()) {
            output.push_back({});
        }
        auto diracInput = mkDirac<T>(output.size());
        Tr::adaptCoefficients(output);
        testGeneric(conv, p, coefficients, diracInput, output);
    }
    
    // test using a random (but reproducible) signal
    if(!Tr::evenIndexesAreApproximated)
    {
        // brute-force convolution is costly (N^2) for high number of coefficients, hence we use caching
        static std::map<Coeffs, std::pair<std::vector<T>, std::vector<T>>> cache;
        
        std::vector<T> randomInput, output;
        
        auto it = cache.find(coefficients);
        if(it == cache.end()) {
            const int sz = 5*coefficients.size();
            output.reserve(sz);
            randomInput.reserve(sz);
            std::srand(0); // to have reproducible random numbers
            for(int i=0; i<sz; ++i) {
                randomInput.push_back(randf(1.f, -1.f));
            }
            // produce expected output using brute force convolution
            {
                FIRFilter<T, a64::Alloc> filter;
                filter.setCoefficients(coefficients);
                if(filter.handlesCoefficients()) {
                    EXPECT_EQ(0, filter.getLatency().toInteger());
                }
                for(int i=0; i<randomInput.size(); ++i) {
                    output.push_back(filter.step(randomInput[i]));
                }
                // "flush" the reverb
                for(int i=0; i<coefficients.size(); ++i) {
                    output.push_back(filter.step(0.));
                }
            }
            
            cache.emplace(coefficients, std::make_pair(randomInput, output));
        }
        else {
            randomInput = it->second.first;
            output = it->second.second;
        }
        
        Tr::adaptCoefficients(output);
        testGeneric(conv, p, coefficients, randomInput, output);
    }
}

template<typename T, typename F>
void testPartitionned(int coeffs_index, F f) {
    const auto coefficients = makeCoefficients<T>(coeffs_index);
    
    if(coefficients.size() < 1024) {
        for(int i=0; i<5;i++)
        {
            const auto part_size = pow2(i);
            f(part_size, coefficients);
        }
    }
    else {
        constexpr auto part_size = 256;
        f(part_size, coefficients);
    }
}

enum class TestFinegrained {
    Begin,
    
    Low = Begin,
    Med,
    High,
    
    End
};

template<typename T, template<typename> typename Allocator, typename Tag>
void testDiracFinegrainedPartitionned(int coeffs_index) {
    auto f = [](int part_size, auto & coefficients)
    {
        for(auto type = TestFinegrained::Begin;
            type != TestFinegrained::End;
            increment(type))
        {
            SelfContainedXYConvolution<AlgoFinegrainedPartitionnedFFTConvolution<T, Allocator, Tag>> conv;
            auto const n_partitions = imajuscule::countPartitions(coefficients.size(), part_size);
            
            typename decltype(conv)::SetupParam p {
                part_size,
                n_partitions,
                1000,
                0
            };
            
            conv.setupAndSetCoefficients(p, coefficients);
            if(!conv.isValid()) {
                continue;
            }
            
            range<int> r {
                p.getLowestValidMultiplicationsGroupSize(),
                p.getHighestValidMultiplicationsGroupSize()
            };
            
            switch(type) {
                case TestFinegrained::Low:
                    p.multiplication_group_size = r.getMin();
                    break;
                case TestFinegrained::High:
                    p.multiplication_group_size = r.getMax();
                    break;
                case TestFinegrained::Med:
                    p.multiplication_group_size = r.getExpCenter();
                    break;
                default:
                    throw std::logic_error("not supported");
            }
            test(conv, p, coefficients);
        }
    };
    
    testPartitionned<T>(coeffs_index, f);
}

template<typename T, template<typename> typename Allocator, typename Tag>
void testDiracPartitionned(int coeffs_index) {
    
    auto f = [](int part_size, auto & coefficients){
        SelfContainedXYConvolution<AlgoFFTConvolutionIntermediate<AlgoPartitionnedFFTConvolutionCRTP<T, Allocator, Tag>>> conv;
        
        test(conv, {
            part_size,
            countPartitions(coefficients.size(),
                            part_size)
        }, coefficients);
    };
    
    testPartitionned<T>(coeffs_index, f);
}


template<typename Convolution>
void testDirac2(int coeffs_index, Convolution & conv, typename Convolution::SetupParam const & p) {
    
    if(!conv.isValid()) {
        return;
    }
    
    using T = typename Convolution::FPT;
    
    test(conv, p, makeCoefficients<T>(coeffs_index));
}

template<typename T, template<typename> typename Allocator, typename Tag>
void testDirac() {
    using namespace fft;
    int end = end_index;
    if constexpr (!std::is_same_v<Tag, fft::Fastest>) {
        --end;
    }
    for(int i=-1; i<end; ++i) {
        LG(INFO,"index %d", i);
        int const countCoeffs = makeCoefficients<T>(i).size();
        testDiracFinegrainedPartitionned<T, Allocator, Tag>(i);
        testDiracPartitionned<T, Allocator, Tag>(i);
        
        if(i<10)
        {
            auto c = SelfContainedXYConvolution<AlgoFIRFilter<T, Allocator, Tag>>{};
            testDirac2(i, c, {countCoeffs});
        }
        {
            auto c = SelfContainedXYConvolution<AlgoFFTConvolutionIntermediate<AlgoFFTConvolutionCRTP<T, Allocator, Tag>>>{};
            testDirac2(i, c, {static_cast<int>(ceil_power_of_two(countCoeffs))});
        }
        for(int partition_sz = std::max(1, static_cast<int>(floor_power_of_two(countCoeffs/50)));
            partition_sz <= ceil_power_of_two(countCoeffs);
            partition_sz *= 2)
        {
            auto c = SelfContainedXYConvolution<AlgoFFTConvolutionIntermediate<AlgoPartitionnedFFTConvolutionCRTP<T, Allocator, Tag>>>{};
            testDirac2(i, c, {
                partition_sz,
                countPartitions(countCoeffs, partition_sz)
            });
        }
        
        for(int partition_sz = std::max(1, static_cast<int>(floor_power_of_two(countCoeffs/50)));
            partition_sz <= ceil_power_of_two(countCoeffs);
            partition_sz *= 2)
        {
            auto c = SelfContainedXYConvolution<AlgoFinegrainedPartitionnedFFTConvolution<T, Allocator, Tag>>{};
            auto const n_partitions = imajuscule::countPartitions(countCoeffs, partition_sz);
            
            testDirac2(i, c, {
                partition_sz,
                n_partitions,
                std::max(1, n_partitions/4), // size of multiplication group
                0
            });
        }
        
        for(int firstSz=1; firstSz <= 32; firstSz *= 2) {
            auto c = SelfContainedXYConvolution<AlgoCustomScaleConvolution<AlgoFFTConvolutionIntermediate<AlgoFFTConvolutionCRTP<T, Allocator, Tag>>>>{};
            using ScalingParam = typename decltype(c)::SetupParam::ScalingParam;
            std::vector<ScalingParam> params;
            int remainingCoeffs = countCoeffs;
            int sz = firstSz;
            while(remainingCoeffs > 0) {
                ScalingParam s{
                    // cannot do std::min(sz, remainingCoeffs) here else the latency would be wrong
                    sz,
                    {sz}
                };
                
                params.push_back(s);
                
                remainingCoeffs -= sz;
                sz *= 2;
            }
            testDirac2(i, c, {params});
        }
        
        // same as above but with some PartitionnedFFTConvolutionCRTP
        {
            for(int firstSz=1; firstSz<32; firstSz *= 2) {
                auto c = SelfContainedXYConvolution<AlgoCustomScaleConvolution<AlgoFFTConvolutionIntermediate<AlgoPartitionnedFFTConvolutionCRTP<T, Allocator, Tag>>>>{};
                using ScalingParam = typename decltype(c)::SetupParam::ScalingParam;
                std::vector<ScalingParam> params;
                int remainingCoeffs = countCoeffs;
                int sz = firstSz;
                while(remainingCoeffs > 0) {
                    int countCoeffs = std::min(sz, remainingCoeffs); // last coefficient group may be padded with 0s
                    ScalingParam s {
                        countCoeffs,
                        {
                            sz,
                            countPartitions(countCoeffs,sz)
                        }
                    };
                    
                    params.push_back(s);
                    
                    remainingCoeffs -= sz;
                    sz *= 2;
                }
                testDirac2(i, c, {params});
            }
        }
        // same as above but skipping one scale out of 2
        {
            for(int firstSz=1; firstSz<32; firstSz *= 2) {
                auto c = SelfContainedXYConvolution<AlgoCustomScaleConvolution<AlgoFFTConvolutionIntermediate<AlgoPartitionnedFFTConvolutionCRTP<T, Allocator, Tag>>>>{};
                using ScalingParam = typename decltype(c)::SetupParam::ScalingParam;
                std::vector<ScalingParam> params;
                int remainingCoeffs = countCoeffs;
                int sz = firstSz;
                while(remainingCoeffs > 0) {
                    int const countCoeffs = std::min(sz*3, remainingCoeffs); // last coefficient group may be padded with 0s
                    ScalingParam s {
                        countCoeffs,
                        {
                            sz,
                            countPartitions(countCoeffs,sz)
                        }
                    };
                    
                    params.push_back(s);
                    
                    remainingCoeffs -= countCoeffs;
                    sz *= 4;
                }
                testDirac2(i, c, {params});
            }
        }
        // custom scale followed by finegrained with same block size as the last in the scale
        {
            for(int firstSz=1; firstSz<32; firstSz *= 2) {
                auto c = SelfContainedXYConvolution<AlgoSplitConvolution<
                AlgoCustomScaleConvolution<AlgoFFTConvolutionIntermediate<AlgoPartitionnedFFTConvolutionCRTP<T, Allocator, Tag>>>,
                AlgoFinegrainedPartitionnedFFTConvolution<T, Allocator, Tag>
                >>{};
                using ScalingParam = typename decltype(c)::SetupParam::AParam::ScalingParam;
                std::vector<ScalingParam> paramsFull;
                int remainingCoeffs = countCoeffs;
                int sz = firstSz;
                while(remainingCoeffs > 0) {
                    int countCoeffs = std::min(sz, remainingCoeffs); // last coefficient group may be padded with 0s
                    ScalingParam s {
                        countCoeffs,
                        {
                            sz,
                            countPartitions(countCoeffs,sz)
                        }
                    };
                    
                    paramsFull.push_back(s);
                    
                    remainingCoeffs -= sz;
                    sz *= 2;
                }
                int const nParams = paramsFull.empty() ?
                0:
                std::max(1,
                         static_cast<int>(paramsFull.size()/2));
                
                std::vector<ScalingParam> params{paramsFull.begin(), paramsFull.begin()+nParams};
                std::vector<ScalingParam> discardedParams{paramsFull.begin()+nParams, paramsFull.end()};
                int countDiscardedCoeffs = 0;
                for(auto const & d : discardedParams) {
                    countDiscardedCoeffs += d.countCoeffs;
                }
                int const partition_sz = params.empty() ? 0 : params.back().setupParam.partition_size;
                auto const n_partitions = imajuscule::countPartitions(countDiscardedCoeffs, partition_sz);
                
                testDirac2(i, c, {
                    {
                        params
                    },
                    {
                        partition_sz,
                        n_partitions,
                        std::max(1, n_partitions/4), // size of multiplication group
                        0
                    }
                });
            }
        }
        {
            auto c = SelfContainedXYConvolution<AlgoAsyncCPUConvolution<AlgoFIRFilter<T, Allocator, Tag>, PolicyOnWorkerTooSlow::Wait>>{};
            std::vector<int> submisionPeriods{
                -1, // invalid
                0, // invalid
                1,
                10,
                100
            };
            std::vector<int> queueSizes{
                -1, // invalid
                0, // invalid
                1,
                10,
                100
            };
            for(auto submisionPeriod:submisionPeriods) {
                for(auto queueSize:queueSizes) {
                    testDirac2(i, c, {
                        submisionPeriod,
                        queueSize,
                        {countCoeffs}
                    });
                }
            }
        }
    }
}
}
}

TEST(Convolution2, dirac) {
    using namespace imajuscule;
    using namespace imajuscule::testdsp2conv;
    
    for_each(fft::Tags, [](auto t) {
        testDirac<float, a64::Alloc, decltype(t)>();
        testDirac<double, a64::Alloc, decltype(t)>();
    });
}
