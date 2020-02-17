
namespace imajuscule::testspatialize {
    constexpr auto end_index = 4;
    
    template<typename T>
    a64::vector<T> makeCoefficients(int coeffs_index) {
        switch(coeffs_index) {
            case 0: return {{ .9,.8,.7,.6,.3,.2,.1,0. }};
            case 1: return {{ +1. }};
            case 2: return {{ -1. }};
            case 3: {
                constexpr auto sz = 2000;
                a64::vector<T> v(sz);
                auto index = 0;
                for(auto & value: v) {
                    value = (sz - index) / static_cast<T>(sz);
                    ++index;
                }
                return std::move(v);
            }
        }
        throw std::logic_error("coeff index too big");
    }
    
    template<typename Spatialize, typename T = typename Spatialize::FPT>
    void test(Spatialize & spatialize,
              std::array<a64::vector<T>, Spatialize::nEars> const & ear_signals,
              int vectorLength,
              bool test_zero) {
        using namespace fft;
        constexpr auto nEars = Spatialize::nEars;
        
        if(!spatialize.isValid()) {
            /*
             std::cout << std::endl << "Not testing invalid setup for "; COUT_TYPE(Convolution);
             std::cout << std::endl <<
             "coefficient size : " << coefficients.size() << std::endl <<
             "partition size : " << conv.getBlockSize() << std::endl <<
             "partition count : " << conv.countPartitions() << std::endl;
             */
            return;
        }
        
        int lat = spatialize.getLatency().toInteger();
        
        std::vector<std::array<T, nEars>> results;
        
        /*
         when test_zero is not set, we feed the entire dirac using assignWetVectorized.
         when test_zero is set, we feed the first frame of the dirac using assignWetVectorized
           and subsequent frames using addWetInputZeroVectorized.
         */
        // feed a dirac
        results.emplace_back();
        std::fill(results.back().begin(), results.back().end(), 1);
        
        while(results.size() < lat + ear_signals[0].size()) {
            results.emplace_back();
            std::fill(results.back().begin(), results.back().end(), 0);
        }
        
        // process the dirac
        Assert(vectorLength);
        std::vector<T> v_output, v_input;
        v_input.reserve(nEars * results.size());
        for(int i=0; i<nEars; ++i) {
            for(auto const & r : results) {
                v_input.push_back(r[i]);
            }
        }
        {
            v_output.resize(nEars * results.size(), {});
            T const * a_inputs[nEars];
            for(int i=0; i<nEars; ++i) {
                a_inputs[i] = &v_input[i*results.size()];
            }
            T * a_outputs[nEars];
            for(int i=0; i<nEars; ++i) {
                a_outputs[i] = &v_output[i*results.size()];
            }
            spatialize.assignWetVectorized(a_inputs,
                                           nEars,
                                           a_outputs,
                                           nEars,
                                           test_zero ? 1 : results.size(),
                                           vectorLength);
        }
        if(test_zero && (results.size() > 1)) {
            T * a_outputs[nEars];
            for(int i=0; i<nEars; ++i) {
                a_outputs[i] = &v_output[i*results.size() + 1];
            }
            spatialize.addWetInputZeroVectorized(a_outputs,
                                                 nEars,
                                                 results.size() - 1,
                                                 vectorLength);
        }
        for(int j=0; j<results.size(); ++j) {
            for(int i=0; i<nEars; ++i) {
                results[j][i] = v_output[i*results.size() + j];
            }
        }
        
        auto eps = spatialize.getEpsilon();
        for(int j=0; j<lat; ++j) {
            int i=0;
            for(auto r : results[j]) {
                ASSERT_NEAR(0, r, eps);
                ++i;
            }
        }
        for(int j=lat; j<results.size(); ++j) {
            int i=0;
            for(auto r : results[j]) {
                if(!areNear(ear_signals[i][j-lat], r, eps)) { // uses relative error
                    ASSERT_NEAR(ear_signals[i][j-lat], r, eps);
                }
                ++i;
            }
        }
    }
    template<typename Spatialize, typename T = typename Spatialize::FPT>
    void test(Spatialize & spatialize,
              std::array<a64::vector<T>, Spatialize::nEars> const & ear_signals,
              int const maxVectorSz) {
        int sz = std::max(1,
                          static_cast<int>(spatialize.getBiggestScale()));

        std::set<int> vectorLengths;
        vectorLengths.insert(1);
        vectorLengths.insert(2);
        vectorLengths.insert(3);
        vectorLengths.insert(sz/4-1);
        vectorLengths.insert(sz/4+1);
        vectorLengths.insert(sz/4);
        vectorLengths.insert(sz/2);
        vectorLengths.insert(sz);
        
        for(auto l : vectorLengths) {
            if(l <= 0) {
                continue;
            }
            if(l > maxVectorSz) {
                continue;
            }
            test(spatialize, ear_signals, l, false);
            test(spatialize, ear_signals, l, true);
        }
    }
    
    template<typename T, typename F>
    void testSpatialized(int coeffs_index, F f) {
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
    
    void testDiracFinegrainedPartitionned(int coeffs_index) {
        
        auto f = [](int part_size, auto const & coefficients)
        {
            int const max_cb_size = std::max(1,
                                             static_cast<int>(floor_power_of_two(coefficients.size() / 10)));
            for(int cb_size = 1; cb_size <= max_cb_size; cb_size*=2) {
                int const maxVectorSz = maxVectorSizeFromBlockSizeHypothesis(cb_size);
                typename std::remove_const<typename std::remove_reference<decltype(coefficients)>::type>::type zero_coeffs, opposite_coeffs;
                zero_coeffs.resize(coefficients.size());
                int const n_partitions = countPartitions(coefficients.size(), part_size);
                
                XFFTsCostsFactors emptyCostFactors;
                using WorkCplxFreqs = typename fft::RealFBins_<fft::Fastest, double, aP::Alloc>::type;
                WorkCplxFreqs work;
                
                {
                    constexpr auto nOutMono = 1;
                    using Rev = Reverbs<nOutMono, ReverbType::Offline, PolicyOnWorkerTooSlow::Wait>;
                    typename Rev::MemResource::type memory;
                    Rev spatialized;
                    
                    applyBestParams(spatialized,
                                    memory,
                                    1,
                                    {{coefficients}},
                                    work,
                                    cb_size,
                                    maxVectorSz,
                                    44100.,
                                    std::cout,
                                    emptyCostFactors
                                    );
                    
                    test(spatialized, {{coefficients}}, maxVectorSz);
                }
                {
                    constexpr auto nOutStereo = 2;
                    using Rev = Reverbs<nOutStereo, ReverbType::Offline, PolicyOnWorkerTooSlow::Wait>;
                    typename Rev::MemResource::type memory;
                    Rev spatialized;
                    
                    // no cross-talk :
                    // source1 has only left component
                    // source2 has only right component
                    
                    applyBestParams(spatialized,
                                    memory,
                                    2,
                                    {{coefficients, zero_coeffs, zero_coeffs, coefficients}},
                                    work,
                                    cb_size,
                                    maxVectorSz,
                                    44100.,
                                    std::cout,
                                    emptyCostFactors
                                    );
                    
                    test(spatialized, {{coefficients, coefficients}}, maxVectorSz);
                }
                
                {
                    constexpr auto nOutStereo = 2;
                    using Rev = Reverbs<nOutStereo, ReverbType::Offline, PolicyOnWorkerTooSlow::Wait>;
                    typename Rev::MemResource::type memory;
                    Rev spatialized;
                    
                    opposite_coeffs = coefficients;
                    for(auto & o : opposite_coeffs) {
                        o = -o;
                    }
                    // extreme cross-talk :
                    // source1 right component is the opposite of source 2
                    // source2 right component is the opposite of source 1
                    
                    applyBestParams(spatialized,
                                    memory,
                                    2,
                                    {{coefficients, opposite_coeffs, opposite_coeffs, coefficients}},
                                    work,
                                    cb_size,
                                    maxVectorSz,
                                    44100.,
                                    std::cout,
                                    emptyCostFactors
                                    );
                    
                    
                    test(spatialized, {{zero_coeffs,zero_coeffs}}, maxVectorSz);
                }
                {
                    constexpr auto nOutStereo = 2;
                    using Rev = Reverbs<nOutStereo, ReverbType::Offline, PolicyOnWorkerTooSlow::Wait>;
                    typename Rev::MemResource::type memory;
                    Rev spatialized;
                    
                    opposite_coeffs = coefficients;
                    for(auto & o : opposite_coeffs) {
                        o = -o;
                    }
                    
                    applyBestParams(spatialized,
                                    memory,
                                    2,
                                    {{coefficients, opposite_coeffs, zero_coeffs, coefficients}},
                                    work,
                                    cb_size,
                                    maxVectorSz,
                                    44100.,
                                    std::cout,
                                    emptyCostFactors
                                    );
                    
                    test(spatialized, {{zero_coeffs,coefficients}}, maxVectorSz);
                }
            }
        };
        
        testSpatialized<double>(coeffs_index, f);
    }
    
    template<typename Tag, template<typename> typename Alloc>
    void testDirac() {
        using namespace fft;
        int end = end_index;
        if constexpr (!std::is_same_v<Tag, fft::Fastest>) {
            --end;
        }
        for(int i=0; i<end; ++i) {
            testDiracFinegrainedPartitionned(i);
        }
    }
}

TEST(Spatialization, dirac) {
    using namespace imajuscule;
    using namespace imajuscule::testspatialize;
    
    for_each(fft::Tags, [](auto t) {
        testDirac<decltype(t), a64::Alloc>();
    });
}

