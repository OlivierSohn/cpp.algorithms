

namespace imajuscule {
    namespace testfftfbins {
    template<typename Tag, typename T, template<typename> typename Allocator>
    void testMult() {
        using FBins = fft::RealFBins_<Tag, T, Allocator>;
        using namespace fft::slow_debug;
        
        constexpr auto N = 6;
        
        auto f1 = FBins::make({{    {2,0}, {0,+1}, {0,-1}, {4,0}, {0,-1}, {0,+1} }});
        auto f2 = FBins::make({{    {3,0}, {0,-1}, {0,+1}, {5,0}, {0,+1}, {0,-1} }});
        auto accum = FBins::make({{ {1,0}, {3,4},  {5,6}, {2,0}, {3,-4}, {5,-6} }});
        
        FBins::multiply(accum.data(),
                        f1.data(),
                        f2.data(),
                        N/2);
        
        auto const res = unwrap_frequencies<Tag>(accum, N);
        
        EXPECT_FLOAT_EQ(3*2,  res[0].real());
        EXPECT_FLOAT_EQ(0,    res[0].imag());
        
        EXPECT_FLOAT_EQ(1, res[1].real());
        EXPECT_FLOAT_EQ(0, res[1].imag());
        
        EXPECT_FLOAT_EQ(1, res[2].real());
        EXPECT_FLOAT_EQ(0, res[2].imag());
        
        EXPECT_FLOAT_EQ(4*5, res[3].real());
        EXPECT_FLOAT_EQ(0,   res[3].imag());
    }
    template<typename Tag>
    void testMult() {
        testMult<Tag, float, a64::Alloc>();
        testMult<Tag, double, a64::Alloc>();
    }

        template<typename Tag, typename T, template<typename> typename Allocator>
        void testMultAdd() {
            using FBins = fft::RealFBins_<Tag, T, Allocator>;
            using namespace fft::slow_debug;
            
            constexpr auto N = 6;
            
            auto f1 = FBins::make({{    {2,0}, {0,+1}, {0,-1}, {4,0}, {0,-1}, {0,+1} }});
            auto f2 = FBins::make({{    {3,0}, {0,-1}, {0,+1}, {5,0}, {0,+1}, {0,-1} }});
            auto accum = FBins::make({{ {1,0}, {3,4},  {5,6}, {2,0}, {3,-4}, {5,-6} }});
            
            FBins::multiply_add(accum.data(),
                                f1.data(),
                                f2.data(),
                                N/2);
            
            auto const res = unwrap_frequencies<Tag>(accum, N);
            
            EXPECT_FLOAT_EQ(1 + 3*2,  res[0].real());
            EXPECT_FLOAT_EQ(0,        res[0].imag());
            
            EXPECT_FLOAT_EQ(3 + 1, res[1].real());
            EXPECT_FLOAT_EQ(4 + 0, res[1].imag());
            
            EXPECT_FLOAT_EQ(5 + 1, res[2].real());
            EXPECT_FLOAT_EQ(6 + 0, res[2].imag());
            
            EXPECT_FLOAT_EQ(2 + 4*5, res[3].real());
            EXPECT_FLOAT_EQ(0,       res[3].imag());
        }
        
        template<typename Tag>
        void testMultAdd() {
            testMultAdd<Tag, float, a64::Alloc>();
            testMultAdd<Tag, double, a64::Alloc>();
        }
        
        template<typename Tag, typename T, template<typename> typename Allocator>
        void testMultAssign() {
            using FBins = fft::RealFBins_<Tag, T, Allocator>;
            using namespace fft::slow_debug;
            
            constexpr auto N = 6;
            
            auto f1 = FBins::make({{    {2,0}, {0,+1}, {0,-1}, {4,0}, {0,-1}, {0,+1} }});
            auto f2 = FBins::make({{    {3,0}, {0,-1}, {0,+1}, {5,0}, {0,+1}, {0,-1} }});
            
            FBins::mult_assign(f1, f2);
            
            auto const res = unwrap_frequencies<Tag>(f1, N);
            
            EXPECT_FLOAT_EQ(3*2, res[0].real());
            EXPECT_FLOAT_EQ(0,   res[0].imag());
            
            EXPECT_FLOAT_EQ(1, res[1].real());
            EXPECT_FLOAT_EQ(0, res[1].imag());
            
            EXPECT_FLOAT_EQ(1, res[2].real());
            EXPECT_FLOAT_EQ(0, res[2].imag());
            
            EXPECT_FLOAT_EQ(4*5, res[3].real());
            EXPECT_FLOAT_EQ(0,   res[3].imag());
        }
        
        template<typename Tag, template<typename> typename Allocator>
        void testMultAssign() {
            testMultAssign<Tag, float, Allocator>();
            testMultAssign<Tag, double, Allocator>();
        }
        
        template<typename Tag, typename T, template<typename> typename Allocator>
        void testZero() {
            using FBins = fft::RealFBins_<Tag, T, Allocator>;
            using namespace fft::slow_debug;
            
            constexpr auto N = 6;
            
            auto f = FBins::make({{    {2,0}, {0,+1}, {0,-1}, {4,0}, {0,-1}, {0,+1} }});
            
            FBins::zero(f);
            
            auto const res = unwrap_frequencies<Tag>(f, N);
            
            EXPECT_FLOAT_EQ(0, res[0].real());
            EXPECT_FLOAT_EQ(0, res[0].imag());
            
            EXPECT_FLOAT_EQ(0, res[1].real());
            EXPECT_FLOAT_EQ(0, res[1].imag());
            
            EXPECT_FLOAT_EQ(0, res[2].real());
            EXPECT_FLOAT_EQ(0, res[2].imag());
            
            EXPECT_FLOAT_EQ(0, res[3].real());
            EXPECT_FLOAT_EQ(0,  res[3].imag());
        }
        
        template<typename Tag, template<typename> typename Allocator>
        void testZero() {
            testZero<Tag, float, Allocator>();
            testZero<Tag, double, Allocator>();
        }
        
        template<typename Tag, typename T, template<typename> typename Allocator>
        void testSqMag() {
            using FBins = fft::RealFBins_<Tag, T, Allocator>;
            using namespace fft::slow_debug;
            using namespace fft;
            
            {
                auto f = FBins::make({{    {2,0}, {0,+1}, {0,-1}, {4,0}, {0,-1}, {0,+1} }});
                auto p = FBins::getMaxSquaredAmplitude(f);
                EXPECT_EQ(3, p.first);
            }
            {
                auto f = FBins::make({{    {1,0}, {1,1}, {0,-1}, {0,0}, {0,-1}, {0,+1} }});
                auto p = FBins::getMaxSquaredAmplitude(f);
                EXPECT_EQ(1, p.first);
            }
            {
                using Sig = fft::RealSignal_<Tag, T>;
                auto highest_sine = Sig::make({{ 1, 1, 1, 1, 1, 1, 1, 1 }});
                using Algo = Algo_<Tag, T>;
                Algo a;
                auto size_fft = highest_sine.size();
                typename FBins::type f(size_fft);
                ScopedContext_<Tag, T> sc(size_fft);
                a.setContext(sc.get());
                a.forward(highest_sine.begin(),f.data(),size_fft);
                auto p = FBins::getMaxSquaredAmplitude(f);
                EXPECT_EQ(0, p.first);
                EXPECT_EQ(1, p.second);
            }
            {
                using Sig = fft::RealSignal_<Tag, T>;
                auto highest_sine = Sig::make({{ -1, +1, -1, +1, -1, +1, -1, +1 }});
                using Algo = Algo_<Tag, T>;
                Algo a;
                auto size_fft = highest_sine.size();
                auto nyquist_bin_index = size_fft / 2;
                typename FBins::type f(size_fft);
                ScopedContext_<Tag, T> sc(size_fft);
                a.setContext(sc.get());
                a.forward(highest_sine.begin(),f.data(),size_fft);
                auto p = FBins::getMaxSquaredAmplitude(f);
                EXPECT_EQ(nyquist_bin_index, p.first);
                EXPECT_EQ(1, p.second);
            }
            {
                using Sig = fft::RealSignal_<Tag, T>;
                auto highest_sine = Sig::make({{ -1, 0, +1, 0, -1, 0, +1, 0 }});
                using Algo = Algo_<Tag, T>;
                Algo a;
                auto size_fft = highest_sine.size();
                auto nyquist_bin_index = size_fft / 2;
                typename FBins::type f(size_fft);
                ScopedContext_<Tag, T> sc(size_fft);
                a.setContext(sc.get());
                a.forward(highest_sine.begin(),f.data(),size_fft);
                auto p = FBins::getMaxSquaredAmplitude(f);
                EXPECT_EQ(nyquist_bin_index/2, p.first);
                EXPECT_EQ(0.25, p.second);
            }
            {
                using Sig = fft::RealSignal_<Tag, T>;
                auto highest_sine = Sig::make({{ 1, .9, .7, .3, -.3, -.7, -.9, -1 }});
                using Algo = Algo_<Tag, T>;
                Algo a;
                auto size_fft = highest_sine.size();
                typename FBins::type f(size_fft);
                ScopedContext_<Tag, T> sc(size_fft);
                a.setContext(sc.get());
                a.forward(highest_sine.begin(),f.data(),size_fft);
                auto p = FBins::getMaxSquaredAmplitude(f);
                EXPECT_TRUE(1 == p.first || 7 == p.first);
            }
        }
        
        template<typename Tag, template<typename> typename Allocator>
        void testSqMag() {
            testSqMag<Tag, float, Allocator>();
            testSqMag<Tag, double, Allocator>();
        }
    }
}

TEST(FFTFBins, mult) {
    using namespace imajuscule;
    using namespace imajuscule::testfftfbins;
    
    for_each(fft::Tags, [](auto t) {
        testMult<decltype(t)>();
    });
}

TEST(FFTFBins, mult_add) {
    using namespace imajuscule;
    using namespace imajuscule::testfftfbins;
    
    for_each(fft::Tags, [](auto t) {
        testMultAdd<decltype(t)>();
    });
}

TEST(FFTFBins, mult_assign) {
    using namespace imajuscule;
    using namespace imajuscule::testfftfbins;
    
    for_each(fft::Tags, [](auto t) {
        testMultAssign<decltype(t), a64::Alloc>();
    });
}

TEST(FFTFBins, zero) {
    using namespace imajuscule;
    using namespace imajuscule::testfftfbins;
    
    for_each(fft::Tags, [](auto t) {
        testZero<decltype(t), a64::Alloc>();
    });
}

TEST(FFTFBins, sqMag) {
    using namespace imajuscule;
    using namespace imajuscule::testfftfbins;
    
    for_each(fft::Tags, [](auto t) {
        testSqMag<decltype(t), a64::Alloc>();
    });
}
