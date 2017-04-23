

namespace imajuscule {
    namespace testfftsignal {
        template<typename Tag, typename T>
        void testGetSignal() {
            using Sig = fft::RealSignal_<Tag, T>;
            
            a_64::vector<T> const signal {{ 1, 2, 3 }};
            
            auto sig = Sig::make(signal);
            
            EXPECT_EQ(sig.size(), signal.size());
            for(int i=0; i<sig.size(); ++i) {
                EXPECT_EQ(signal[i], Sig::get_signal(sig[i]));
            }
        }
        
        template<typename Tag>
        void testGetSignal() {
            testGetSignal<Tag, float>();
            testGetSignal<Tag, double>();
        }
        
        template<typename Tag, typename T>
        void testASM() {
            using Sig = fft::RealSignal_<Tag, T>;
            using namespace fft::slow_debug;
            
            auto constexpr N = 3;
            auto add1 = Sig::make({{ 1, 2, 3 }});
            auto add2 = Sig::make({{ 4, 6, 8 }});
            auto m = 5;
            
            auto result = Sig::make({{ 0, 0, 0 }});
            
            Sig::add_scalar_multiply(result.begin(), add1.begin(), add2.begin(), m, N);
            
            auto const res = unwrap_signal<Tag>(result, N);
            
            EXPECT_FLOAT_EQ(m*(1+4), res[0].real());
            EXPECT_FLOAT_EQ(m*(2+6), res[1].real());
            EXPECT_FLOAT_EQ(m*(3+8), res[2].real());
            EXPECT_FLOAT_EQ(0, res[0].imag());
            EXPECT_FLOAT_EQ(0, res[1].imag());
            EXPECT_FLOAT_EQ(0, res[2].imag());
        }
        
        template<typename Tag>
        void testASM() {
            testASM<Tag, float>();
            testASM<Tag, double>();
        }
        
        template<typename Tag, typename T>
        void testCopy() {
            using Sig = fft::RealSignal_<Tag, T>;
            using namespace fft::slow_debug;
            
            auto constexpr N = 3;
            auto source = Sig::make({{ 1, 2, 3 }});
            auto dest = Sig::make({{ 0, 0, 0 }});
            
            Sig::copy(dest.begin(), source.begin(), N);
            
            auto const res = unwrap_signal<Tag>(dest, N);
            
            EXPECT_FLOAT_EQ(1, res[0].real());
            EXPECT_FLOAT_EQ(2, res[1].real());
            EXPECT_FLOAT_EQ(3, res[2].real());
            EXPECT_FLOAT_EQ(0, res[0].imag());
            EXPECT_FLOAT_EQ(0, res[1].imag());
            EXPECT_FLOAT_EQ(0, res[2].imag());
        }
        
        template<typename Tag>
        void testCopy() {
            testCopy<Tag, float>();
            testCopy<Tag, double>();
        }
    }
}

TEST(FFTSignal, get_signal) {
    using namespace imajuscule;
    using namespace imajuscule::testfftsignal;
    
    for_each(fft::Tags, [](auto t) {
        testGetSignal<decltype(t)>();
    });
}


TEST(FFTSignal, add_scalar_multiply) {
    using namespace imajuscule;
    using namespace imajuscule::testfftsignal;
    
    for_each(fft::Tags, [](auto t) {
        testASM<decltype(t)>();
    });
}

TEST(FFTSignal, copy) {
    using namespace imajuscule;
    using namespace imajuscule::testfftsignal;
    
    for_each(fft::Tags, [](auto t) {
        testCopy<decltype(t)>();
    });
}
