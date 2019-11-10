
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
        int countSamples() const { return v.size(); }
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
    
    template<typename T, typename T2>
    void mkSameSizesWith(std::vector<std::vector<T>> & vv, T2 fill) {
        std::size_t maxSz = 0;
        for(auto const & v:vv) {
            maxSz = std::max(maxSz, v.size());
        }
        for(auto & v:vv) {
            v.resize(maxSz, fill);
        }
    }

    template<typename FLT, typename FLT2>
    void testSincVariable() {
        constexpr double epsilon = 10*
        std::max(static_cast<double>(std::numeric_limits<FLT>::epsilon()),
                 static_cast<double>(std::numeric_limits<FLT2>::epsilon()));
        
        FloatsVecsEqual floatVecsEqual = {epsilon};
        /*
        {
            float sample_rate_from = 10.f;
            VectorReader<FLT> reader({0., 1., 0., 1., 0., 1., 0.}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleSincVariableRate(reader, resampled, [](int, double x){ return 10.f + 10.f * x; });
            std::vector<FLT2> resampledExpected {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }

        // 10 downsampling empty
        {
            float sample_rate_from = 10.f;
            float sample_rate_to = 100.f;
            VectorReader<FLT> reader({}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleSincVariableRate(reader, resampled, [sample_rate_to](int, double){ return sample_rate_to; });
            std::vector<FLT2> resampledExpected {};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
         */
        // 10 upsampling empty
        /*
         {
         float sample_rate_from = 100.f;
         float sample_rate_to = 10.f;
         VectorReader<FLT> reader({}, 1, sample_rate_from);
         std::vector<FLT2> resampled;
         resampleSinc(reader, resampled, sample_rate_to);
         std::vector<FLT2> resampledExpected {};
         ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
         }
         */
        // same rate empty
        {
            float sample_rate_from = 100.f;
            float sample_rate_to = 100.f;
            VectorReader<FLT> reader({}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleSincVariableRate(reader, resampled, [sample_rate_to](int, double){ return sample_rate_to; });
            std::vector<FLT2> resampledExpected {};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
        
        // 10 downsampling singleton
        {
            float sample_rate_from = 10.f;
            float sample_rate_to = 100.f;
            VectorReader<FLT> reader({0.5}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleSincVariableRate(reader, resampled, [sample_rate_to](int, double){ return sample_rate_to; });
            std::vector<FLT2> resampledExpected {0.5};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
        // 10 upsampling singleton
        // en commentaire car il faudrait tenir compte de la normalisation
        /*
         {
         float sample_rate_from = 100.f;
         float sample_rate_to = 10.f;
         VectorReader<FLT> reader({0.5}, 1, sample_rate_from);
         std::vector<FLT2> resampled;
         resampleSinc(reader, resampled, sample_rate_to);
         std::vector<FLT2> resampledExpected {0.5};
         ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
         }
         */
        
        {
            std::vector<FLT> buf;
            int const sz = 18;
            int const sincWindowHS = 10000;
            buf.reserve(2*sincWindowHS + sz);
            for(int i=0; i<sincWindowHS; ++i) {
                buf.push_back(1.f);
            }
            for(int i=0; i<sz; ++i) {
                buf.push_back(std::abs((sz/2)-i));
            }
            for(int i=0; i<sincWindowHS; ++i) {
                buf.push_back(1.f);
            }
            
            float sample_rate_from = 100.f;
            float sample_rate_to = 100.f;
            VectorReader<FLT> reader(buf, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleSincVariableRate(reader, resampled, [sample_rate_to](int, double){ return sample_rate_to; });
            std::vector<FLT> shortBuf(buf.begin()+(sincWindowHS),
                                      buf.begin()+(sincWindowHS) + sz);
            std::vector<FLT> shortBufRS(resampled.begin()+(sincWindowHS),
                                        resampled.begin()+(sincWindowHS) + sz);
            ASSERT_TRUE(floatVecsEqual(shortBuf, shortBufRS));
        }
        /*
        {
            KaiserWindow window;
            std::vector<double> v;
            for(int i=-100;i<=100; ++i) {
                v.push_back(window.getAt(i/100.));
            }
            std::vector<std::vector<double>> vv{v};
            write_wav("/Users/Olivier/Dev/Audiofiles", "kaiser_win.wav", vv, 100);
        }*/
        /*
        {
            std::vector<FLT> buf;
            int const sz = 400;
            for(int i=0; i<sz; ++i) {
                //buf.push_back(std::abs((sz/2)-i));
                buf.push_back(i < (sz/2) ? i : (sz-i));
            }
            float sample_rate_from = 100.f;
            float sample_rate_to = 101.f;
            VectorReader<FLT> reader(buf, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleSincVariableRate(reader, resampled, [sample_rate_to](int, double){ return sample_rate_to; });
            std::vector<FLT> resampled2(resampled.begin(), resampled.end());
            std::vector<std::vector<FLT>> vv{buf, resampled2};
            mkSameSizesWith(vv, 0.f);
            write_wav("/Users/Olivier/Dev/Audiofiles", "resampled_origin.wav", vv, sample_rate_from);
        }*/
        /*
        {
            std::vector<FLT> buf;
            int const sz = 400;
            for(int i=0; i<sz; ++i) {
                //buf.push_back(std::abs((sz/2)-i));
                buf.push_back(i < (sz/2) ? i : (sz-i));
            }
            float sample_rate_from = 100.f;
            float sample_rate_to = 99.f;
            VectorReader<FLT> reader(buf, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleSincVariableRate(reader, resampled, [sample_rate_to](int, double){ return sample_rate_to; });
            std::vector<FLT> resampled2(resampled.begin(), resampled.end());
            std::vector<std::vector<FLT>> vv{buf, resampled2};
            mkSameSizesWith(vv, 0.f);
            write_wav("/Users/Olivier/Dev/Audiofiles", "resampled_origin2.wav", vv, sample_rate_from);
        }*/
        
        // same rate singleton
        {
            float sample_rate_from = 100.f;
            float sample_rate_to = 100.f;
            VectorReader<FLT> reader({0.5}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleSincVariableRate(reader, resampled, [sample_rate_to](int, double){ return sample_rate_to; });
            std::vector<FLT2> resampledExpected {0.5};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
        // same rate
        {
            float sample_rate_from = 100.f;
            float sample_rate_to = 100.f;
            VectorReader<FLT> reader({0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleSincVariableRate(reader, resampled, [sample_rate_to](int, double){ return sample_rate_to; });
            std::vector<FLT2> resampledExpected {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
        /*
        // 10 downsampling
        {
            float sample_rate_from = 10.f;
            float sample_rate_to = 100.f;
            VectorReader<FLT> reader({0., 1., 0.}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleSincVariableRate(reader, resampled, [sample_rate_to](int, double){ return sample_rate_to; });
            std::vector<FLT2> resampledExpected {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }*/
        // 10 upsampling
        // en commentaire car il faudrait tenir compte de la normalisation
        /*
         {
         float sample_rate_from = 100.f;
         float sample_rate_to = 10.f;
         VectorReader<FLT> reader({0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}, 1, sample_rate_from);
         std::vector<FLT2> resampled;
         resampleSinc(reader, resampled, sample_rate_to);
         std::vector<FLT2> resampledExpected {0., 1.};
         ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
         }
         */
        // 2 upsampling
        // en commentaire car il faudrait tenir compte de la normalisation
        /*
         {
         float sample_rate_from = 100.f;
         float sample_rate_to = 50.f;
         VectorReader<FLT> reader({0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}, 1, sample_rate_from);
         std::vector<FLT2> resampled;
         resampleSinc(reader, resampled, sample_rate_to);
         std::vector<FLT2> resampledExpected {0., 0.2, 0.4, 0.6, 0.8, 1.};
         ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
         }
         */
        /*
        // 3/2 downsampling
        {
            float sample_rate_from = 100.f;
            float sample_rate_to = 100.f*3/2;
            VectorReader<FLT> reader({0., 0.5, 1.0}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleSincVariableRate(reader, resampled, [sample_rate_to](int, double){ return sample_rate_to; });
            std::vector<FLT2> resampledExpected {0./3, 1./3, 2./3, 3./3};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
        // 3/2 downsampling with reverted input
        {
            float sample_rate_from = 100.f;
            float sample_rate_to = 100.f*3/2;
            VectorReader<FLT> reader({1., 0.5, 0.0}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleSincVariableRate(reader, resampled, [sample_rate_to](int, double){ return sample_rate_to; });
            std::vector<FLT2> resampledExpected {3./3, 2./3, 1./3, 0./3};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
        // 3/2 downsampling with longer input
        {
            float sample_rate_from = 100.f;
            float sample_rate_to = 100.f*3/2;
            VectorReader<FLT> reader({1., 0.0, 1., 0.0, 1., 0.0}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleSincVariableRate(reader, resampled, [sample_rate_to](int, double){ return sample_rate_to; });
            std::vector<FLT2> resampledExpected {3./3, 1./3, 1./3, 3./3, 1./3, 1./3, 3./3, 1./3};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }*/
    }

    template<typename FLT, typename FLT2>
    void testSinc() {
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
            resampleSinc(reader, resampled, sample_rate_to);
            std::vector<FLT2> resampledExpected {};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
        // 10 upsampling empty
        /*
        {
            float sample_rate_from = 100.f;
            float sample_rate_to = 10.f;
            VectorReader<FLT> reader({}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleSinc(reader, resampled, sample_rate_to);
            std::vector<FLT2> resampledExpected {};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
         */
        // same rate empty
        {
            float sample_rate_from = 100.f;
            float sample_rate_to = 100.f;
            VectorReader<FLT> reader({}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleSinc(reader, resampled, sample_rate_to);
            std::vector<FLT2> resampledExpected {};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
        
        // 10 downsampling singleton
        {
            float sample_rate_from = 10.f;
            float sample_rate_to = 100.f;
            VectorReader<FLT> reader({0.5}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleSinc(reader, resampled, sample_rate_to);
            std::vector<FLT2> resampledExpected {0.5};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
        // 10 upsampling singleton
        // en commentaire car il faudrait tenir compte de la normalisation
        /*
        {
            float sample_rate_from = 100.f;
            float sample_rate_to = 10.f;
            VectorReader<FLT> reader({0.5}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleSinc(reader, resampled, sample_rate_to);
            std::vector<FLT2> resampledExpected {0.5};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
         */
        
        {
            std::vector<FLT> buf;
            int const sz = 18;
            int const sincWindowHS = 10000;
            buf.reserve(2*sincWindowHS + sz);
            for(int i=0; i<sincWindowHS; ++i) {
                buf.push_back(1.f);
            }
            for(int i=0; i<sz; ++i) {
                buf.push_back(std::abs((sz/2)-i));
            }
            for(int i=0; i<sincWindowHS; ++i) {
                buf.push_back(1.f);
            }

            float sample_rate_from = 100.f;
            float sample_rate_to = 100.f;
            VectorReader<FLT> reader(buf, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleSinc(reader, resampled, sample_rate_to);
            std::vector<FLT> shortBuf(buf.begin()+(sincWindowHS),
                                      buf.begin()+(sincWindowHS) + sz);
            std::vector<FLT> shortBufRS(resampled.begin()+(sincWindowHS),
                                        resampled.begin()+(sincWindowHS) + sz);
            ASSERT_TRUE(floatVecsEqual(shortBuf, shortBufRS));
        }
        
        {
            KaiserWindow window;
            std::vector<double> v;
            for(int i=-100;i<=100; ++i) {
                v.push_back(window.getAt(i/100.));
            }
            std::vector<std::vector<double>> vv{v};
            write_wav("/Users/Olivier/Dev/Audiofiles", "kaiser_win.wav", vv, 100);
        }
        {
            std::vector<FLT> buf;
            int const sz = 400;
            for(int i=0; i<sz; ++i) {
                //buf.push_back(std::abs((sz/2)-i));
                buf.push_back(i < (sz/2) ? i : (sz-i));
            }
            float sample_rate_from = 100.f;
            float sample_rate_to = 101.f;
            VectorReader<FLT> reader(buf, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleSinc(reader, resampled, sample_rate_to);
            std::vector<FLT> resampled2(resampled.begin(), resampled.end());
            std::vector<std::vector<FLT>> vv{buf, resampled2};
            mkSameSizesWith(vv, 0.f);
            write_wav("/Users/Olivier/Dev/Audiofiles", "resampled_origin.wav", vv, sample_rate_from);
        }
        {
            std::vector<FLT> buf;
            int const sz = 400;
            for(int i=0; i<sz; ++i) {
                //buf.push_back(std::abs((sz/2)-i));
                buf.push_back(i < (sz/2) ? i : (sz-i));
            }
            float sample_rate_from = 100.f;
            float sample_rate_to = 99.f;
            VectorReader<FLT> reader(buf, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleSinc(reader, resampled, sample_rate_to);
            std::vector<FLT> resampled2(resampled.begin(), resampled.end());
            std::vector<std::vector<FLT>> vv{buf, resampled2};
            mkSameSizesWith(vv, 0.f);
            write_wav("/Users/Olivier/Dev/Audiofiles", "resampled_origin2.wav", vv, sample_rate_from);
        }
        
        // same rate singleton
        {
            float sample_rate_from = 100.f;
            float sample_rate_to = 100.f;
            VectorReader<FLT> reader({0.5}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleSinc(reader, resampled, sample_rate_to);
            std::vector<FLT2> resampledExpected {0.5};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
        // same rate
        {
            float sample_rate_from = 100.f;
            float sample_rate_to = 100.f;
            VectorReader<FLT> reader({0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleSinc(reader, resampled, sample_rate_to);
            std::vector<FLT2> resampledExpected {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
        // 10 downsampling
        /*
        {
            float sample_rate_from = 10.f;
            float sample_rate_to = 100.f;
            VectorReader<FLT> reader({0., 1., 0.}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleSinc(reader, resampled, sample_rate_to);
            std::vector<FLT2> resampledExpected {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
         */
        // 10 upsampling
        // en commentaire car il faudrait tenir compte de la normalisation
        /*
        {
            float sample_rate_from = 100.f;
            float sample_rate_to = 10.f;
            VectorReader<FLT> reader({0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleSinc(reader, resampled, sample_rate_to);
            std::vector<FLT2> resampledExpected {0., 1.};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
         */
        // 2 upsampling
        // en commentaire car il faudrait tenir compte de la normalisation
/*
        {
            float sample_rate_from = 100.f;
            float sample_rate_to = 50.f;
            VectorReader<FLT> reader({0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleSinc(reader, resampled, sample_rate_to);
            std::vector<FLT2> resampledExpected {0., 0.2, 0.4, 0.6, 0.8, 1.};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
  */
/*
        // 3/2 downsampling
        {
            float sample_rate_from = 100.f;
            float sample_rate_to = 100.f*3/2;
            VectorReader<FLT> reader({0., 0.5, 1.0}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleSinc(reader, resampled, sample_rate_to);
            std::vector<FLT2> resampledExpected {0./3, 1./3, 2./3, 3./3};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
        // 3/2 downsampling with reverted input
        {
            float sample_rate_from = 100.f;
            float sample_rate_to = 100.f*3/2;
            VectorReader<FLT> reader({1., 0.5, 0.0}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleSinc(reader, resampled, sample_rate_to);
            std::vector<FLT2> resampledExpected {3./3, 2./3, 1./3, 0./3};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
        // 3/2 downsampling with longer input
        {
            float sample_rate_from = 100.f;
            float sample_rate_to = 100.f*3/2;
            VectorReader<FLT> reader({1., 0.0, 1., 0.0, 1., 0.0}, 1, sample_rate_from);
            std::vector<FLT2> resampled;
            resampleSinc(reader, resampled, sample_rate_to);
            std::vector<FLT2> resampledExpected {3./3, 1./3, 1./3, 3./3, 1./3, 1./3, 3./3, 1./3};
            ASSERT_TRUE(floatVecsEqual(resampled, resampledExpected));
        }
 */
    }

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


TEST(Resample, sincVariable) {
    using namespace imajuscule::audio::testresample;

    // same types
    testSincVariable<float, float>();
    testSincVariable<double, double>();
    
    // different types
    testSincVariable<float, double>();
    testSincVariable<double, float>();
}

TEST(Resample, sinc) {
    using namespace imajuscule::audio::testresample;

    // same types
    testSinc<float, float>();
    testSinc<double, double>();
    
    // different types
    testSinc<float, double>();
    testSinc<double, float>();
}

