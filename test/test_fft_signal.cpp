

namespace imajuscule {
    namespace testfftsignal {
        template<typename Tag, typename T>
        void testGetSignal() {
            using Sig = fft::RealSignal_<Tag, T>;
            
            a64::vector<T> const signal {{ 1, 2, 3 }};
            
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
        void testCopy() {
            using Sig = fft::RealSignal_<Tag, T>;
            using namespace fft::slow_debug;
            
            auto constexpr N = 3;
            auto source = Sig::make({{ 1, 2, 3 }});
            auto dest = Sig::make({{ 0, 0, 0 }});
            
            Sig::copy(dest.data(), source.data(), N);
            
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
        
        template<typename Tag, typename T>
        void testZero() {
            using Sig = fft::RealSignal_<Tag, T>;
            using namespace fft::slow_debug;
 
            std::vector<int> vN{{0,1,2,3,4,5,6,7,8,80,81,82,83,84,85,86,87,88,89,90,1001,1002,1003}};
            for(auto N : vN) {
                a64::vector<T> v;
                v.resize(N, 1);
                auto sig = Sig::make(std::move(v));

                {
                    auto const res = unwrap_signal<Tag>(sig, N);
                    for(int i=0; i<N; ++i){
                        ASSERT_EQ(1, res[i].real());
                    }
                }
              
              if(!sig.empty()) {
                Sig::zero(sig);
              }
                
                {
                    auto const res = unwrap_signal<Tag>(sig, N);
                    for(int i=0; i<N; ++i){
                        ASSERT_EQ(0, res[i].real());
                    }
                }
            }
        }
        
        template<typename Tag>
        void testZero() {
            testZero<Tag, float>();
            testZero<Tag, double>();
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

TEST(FFTSignal, copy) {
    using namespace imajuscule;
    using namespace imajuscule::testfftsignal;
    
    for_each(fft::Tags, [](auto t) {
        testCopy<decltype(t)>();
    });
}

TEST(FFTSignal, zero) {
    using namespace imajuscule;
    using namespace imajuscule::testfftsignal;
    
    for_each(fft::Tags, [](auto t) {
        testZero<decltype(t)>();
    });
}
