

namespace imajuscule {
    namespace testfftfbins {
        
        template<typename Tag, typename T>
        void testMultAdd() {
            using FBins = fft::RealFBins_<Tag, T>;
            using namespace fft::slow_debug;
            
            constexpr auto N = 6;
            
            auto f1 = FBins::make({{    {2,0}, {0,+1}, {0,-1}, {4,0}, {0,-1}, {0,+1} }});
            auto f2 = FBins::make({{    {3,0}, {0,-1}, {0,+1}, {5,0}, {0,+1}, {0,-1} }});
            auto accum = FBins::make({{ {1,0}, {3,4},  {5,6}, {2,0}, {3,-4}, {5,-6} }});
            
            FBins::multiply_add(accum, f1, f2);
            
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
            testMultAdd<Tag, float>();
            testMultAdd<Tag, double>();
        }
        
        template<typename Tag, typename T>
        void testMultAssign() {
            using FBins = fft::RealFBins_<Tag, T>;
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
        
        template<typename Tag>
        void testMultAssign() {
            testMultAssign<Tag, float>();
            testMultAssign<Tag, double>();
        }
        
        template<typename Tag, typename T>
        void testZero() {
            using FBins = fft::RealFBins_<Tag, T>;
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
        
        template<typename Tag>
        void testZero() {
            testZero<Tag, float>();
            testZero<Tag, double>();
        }
    }
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
        testMultAssign<decltype(t)>();
    });
}

TEST(FFTFBins, zero) {
    using namespace imajuscule;
    using namespace imajuscule::testfftfbins;
    
    for_each(fft::Tags, [](auto t) {
        testZero<decltype(t)>();
    });
}
