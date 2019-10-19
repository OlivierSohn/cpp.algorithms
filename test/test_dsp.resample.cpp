
namespace imajuscule::audio::testresample {
    template<typename FLT>
    struct VectorReader {
        VectorReader(const VectorReader &) = delete;
        VectorReader(const VectorReader &&) = delete;
        VectorReader & operator = (const VectorReader &) const = delete;
        VectorReader & operator = (const VectorReader &&) const = delete;

        VectorReader(std::vector<FLT> const & v, int nchans, float sample_rate)
        : v(v)
        , nchans(nchans)
        , sample_rate(sample_rate)
        {
            it = this->v.begin();
        }
        
        template<typename FLT2>
        FLT2 ReadAsOneFloat() {
            assert(HasMore());
            FLT2 res = static_cast<FLT2>(*it);
            ++it;
            return res;
        }
        
        bool HasMore() const { return it != v.end(); }
        int countChannels() const { return nchans; }
        float getSampleRate() const { return sample_rate; }
        
        int countFrames() const {
            return static_cast<int>(v.size()) / nchans;
        }

    private:
        std::vector<FLT> v;
        typename std::vector<FLT>::const_iterator it;
        int nchans;
        float sample_rate;
    };

    struct FloatsVecsEqual {
        FloatsVecsEqual(double epsilon)
        : epsilon(epsilon)
        {}
        
        template<typename Vec>
        bool operator()(Vec const & v1, Vec const & v2) const {
            if(v1.size() != v2.size()) {
                return false;
            }
            
            auto it1 = v1.begin();
            auto it2 = v2.begin();
            for(;it1!=v1.end(); ++it1, ++it2) {
                auto asum = std::abs(*it1+*it2);
                auto adiff = std::abs(*it1-*it2);
                if(asum < epsilon) {
                    if(adiff > epsilon) {
                        return false;
                    }
                }
                else {
                    if(adiff / asum > epsilon) {
                        return false;
                    }
                }
            }
            return true;
        }
        
        double const epsilon;
    };
    
    template<typename FLT, typename FLT2>
    void test() {
        constexpr double epsilon = 10*
            std::max(static_cast<double>(std::numeric_limits<FLT>::epsilon()),
                     static_cast<double>(std::numeric_limits<FLT2>::epsilon()));
        
        FloatsVecsEqual floatVecsEqual = {epsilon};
        
        // 10 downsampling empty
        {
            float sample_rate_from = 10.f;
            float sample_rate_to = 100.f;
            VectorReader<FLT> reader({}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleLinear(reader, resampled, sample_rate_to);
            std::vector<FLT2> resampledExpected {};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
        // 10 upsampling empty
        {
            float sample_rate_from = 100.f;
            float sample_rate_to = 10.f;
            VectorReader<FLT> reader({}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleLinear(reader, resampled, sample_rate_to);
            std::vector<FLT2> resampledExpected {};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
        // same rate empty
        {
            float sample_rate_from = 100.f;
            float sample_rate_to = 100.f;
            VectorReader<FLT> reader({}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleLinear(reader, resampled, sample_rate_to);
            std::vector<FLT2> resampledExpected {};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
        
        // 10 downsampling singleton
        {
            float sample_rate_from = 10.f;
            float sample_rate_to = 100.f;
            VectorReader<FLT> reader({0.5}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleLinear(reader, resampled, sample_rate_to);
            std::vector<FLT2> resampledExpected {0.5};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
        // 10 upsampling singleton
        {
            float sample_rate_from = 100.f;
            float sample_rate_to = 10.f;
            VectorReader<FLT> reader({0.5}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleLinear(reader, resampled, sample_rate_to);
            std::vector<FLT2> resampledExpected {0.5};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
        // same rate singleton
        {
            float sample_rate_from = 100.f;
            float sample_rate_to = 100.f;
            VectorReader<FLT> reader({0.5}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleLinear(reader, resampled, sample_rate_to);
            std::vector<FLT2> resampledExpected {0.5};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
        
        // 10 downsampling
        {
            float sample_rate_from = 10.f;
            float sample_rate_to = 100.f;
            VectorReader<FLT> reader({0., 1.}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleLinear(reader, resampled, sample_rate_to);
            std::vector<FLT2> resampledExpected {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
        // 10 upsampling
        {
            float sample_rate_from = 100.f;
            float sample_rate_to = 10.f;
            VectorReader<FLT> reader({0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleLinear(reader, resampled, sample_rate_to);
            std::vector<FLT2> resampledExpected {0., 1.};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
        // 2 upsampling
        {
            float sample_rate_from = 100.f;
            float sample_rate_to = 50.f;
            VectorReader<FLT> reader({0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleLinear(reader, resampled, sample_rate_to);
            std::vector<FLT2> resampledExpected {0., 0.2, 0.4, 0.6, 0.8, 1.};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
        // same rate
        {
            float sample_rate_from = 100.f;
            float sample_rate_to = 100.f;
            VectorReader<FLT> reader({0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleLinear(reader, resampled, sample_rate_to);
            std::vector<FLT2> resampledExpected {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
        // 3/2 downsampling
        {
            float sample_rate_from = 100.f;
            float sample_rate_to = 100.f*3/2;
            VectorReader<FLT> reader({0., 0.5, 1.0}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleLinear(reader, resampled, sample_rate_to);
            std::vector<FLT2> resampledExpected {0./3, 1./3, 2./3, 3./3};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
        // 3/2 downsampling with reverted input
        {
            float sample_rate_from = 100.f;
            float sample_rate_to = 100.f*3/2;
            VectorReader<FLT> reader({1., 0.5, 0.0}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleLinear(reader, resampled, sample_rate_to);
            std::vector<FLT2> resampledExpected {3./3, 2./3, 1./3, 0./3};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
        // 3/2 downsampling with longer input
        {
            float sample_rate_from = 100.f;
            float sample_rate_to = 100.f*3/2;
            VectorReader<FLT> reader({1., 0.0, 1., 0.0, 1., 0.0}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleLinear(reader, resampled, sample_rate_to);
            std::vector<FLT2> resampledExpected {3./3, 1./3, 1./3, 3./3, 1./3, 1./3, 3./3, 1./3};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
    }

} // namespace imajuscule::testresample

TEST(Resample, linear) {
    using namespace imajuscule::audio::testresample;

    // same types
    test<float, float>();
    test<double, double>();
    
    // different types
    test<float, double>();
    test<double, float>();
}

