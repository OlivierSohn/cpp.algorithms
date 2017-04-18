

namespace imajuscule {
    namespace testfftsignal {
        template<typename Tag, typename T>
        void testGetSignal() {
            using Sig = fft::RealSignal_<Tag, T>;
            
            std::vector<T> const signal {{ 1, 2, 3 }};
            
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
    }
}

TEST(FFTSignal, mult) {
    using namespace imajuscule;
    using namespace imajuscule::testfftsignal;
    
    for_each(fft::Tags, [](auto t) {
        testGetSignal<decltype(t)>();
    });
}
