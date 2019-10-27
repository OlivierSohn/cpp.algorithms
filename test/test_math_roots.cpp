TEST(MathRoots, find_roots) {
    using namespace imajuscule;

    // constant positive function
    {
        std::optional<double> res = find_root([](double){ return 1.; },
                                              [](double){ return 0.; },
                                              -1.,
                                              1.);
        ASSERT_FALSE(res);
    }

    {
        std::optional<double> res = find_root([](double){ return 1.; },
                                              [](double){ return 0.; },
                                              1.,
                                              -1.);
        ASSERT_FALSE(res);
    }
    
    // constant negative function
    {
        std::optional<double> res = find_root([](double){ return -1.; },
                                              [](double){ return 0.; },
                                              -1.,
                                              1.);
        ASSERT_FALSE(res);
    }

    {
        std::optional<double> res = find_root([](double){ return -1.; },
                                              [](double){ return 0.; },
                                              1.,
                                              -1.);
        ASSERT_FALSE(res);
    }

    // identity function
    {
        std::optional<double> res = find_root([](double x){ return x; },
                                              [](double){ return 1.; },
                                              1.,
                                              -1.);
        ASSERT_EQ(0, *res);
    }

    {
        std::optional<double> res = find_root([](double x){ return x; },
                                              [](double){ return 1.; },
                                              -1.,
                                              1.);
        ASSERT_EQ(0, *res);
    }

    {
        std::optional<double> res = find_root([](double x){ return x; },
                                              [](double){ return 1.; },
                                              10.,
                                              -1.);
        ASSERT_EQ(0, *res);
    }

    {
        std::optional<double> res = find_root([](double x){ return x; },
                                              [](double){ return 1.; },
                                              -10.,
                                              1.);
        ASSERT_EQ(0, *res);
    }
    
    {
        std::optional<double> res = find_root([](double x){ return x; },
                                              [](double){ return 1.; },
                                              2.,
                                              1.);
        ASSERT_FALSE(res);
    }
    
    {
        std::optional<double> res = find_root([](double x){ return x; },
                                              [](double){ return 1.; },
                                              1.,
                                              2.);
        ASSERT_FALSE(res);
    }
    
    {
        std::optional<double> res = find_root([](double x){ return x; },
                                              [](double){ return 1.; },
                                              -2.,
                                              -1.);
        ASSERT_FALSE(res);
    }
    
    {
        std::optional<double> res = find_root([](double x){ return x; },
                                              [](double){ return 1.; },
                                              -1.,
                                              -2.);
        ASSERT_FALSE(res);
    }

    // opposite identity function
    {
        std::optional<double> res = find_root([](double x){ return -x; },
                                              [](double){ return -1.; },
                                              1.,
                                              -1.);
        ASSERT_EQ(0, *res);
    }

    {
        std::optional<double> res = find_root([](double x){ return -x; },
                                              [](double){ return -1.; },
                                              -1.,
                                              1.);
        ASSERT_EQ(0, *res);
    }

    {
        std::optional<double> res = find_root([](double x){ return -x; },
                                              [](double){ return -1.; },
                                              10.,
                                              -1.);
        ASSERT_EQ(0, *res);
    }

    {
        std::optional<double> res = find_root([](double x){ return -x; },
                                              [](double){ return -1.; },
                                              -10.,
                                              1.);
        ASSERT_EQ(0, *res);
    }

    {
        std::optional<double> res = find_root([](double x){ return -x; },
                                              [](double){ return -1.; },
                                              2.,
                                              1.);
        ASSERT_FALSE(res);
    }
    
    {
        std::optional<double> res = find_root([](double x){ return -x; },
                                              [](double){ return -1.; },
                                              1.,
                                              2.);
        ASSERT_FALSE(res);
    }
    
    {
        std::optional<double> res = find_root([](double x){ return -x; },
                                              [](double){ return -1.; },
                                              -2.,
                                              -1.);
        ASSERT_FALSE(res);
    }
    
    {
        std::optional<double> res = find_root([](double x){ return -x; },
                                              [](double){ return -1.; },
                                              -1.,
                                              -2.);
        ASSERT_FALSE(res);
    }
    
    auto constexpr epsilon = std::numeric_limits<double>::epsilon();

    // x^3

    {
        auto const f = [](double x){ return x * x * x; };
        auto const fder = [](double x){ return 3. * x * x; };
        auto testCube = [f, fder] (auto a, auto b) {
            return find_root(f, fder, a, b);
        };
        {
            std::optional<double> res = testCube(-1.,2.);
            ASSERT_LT(std::abs(f(*res)), epsilon);
        }
        {
            std::optional<double> res = testCube(0.,2.);
            ASSERT_LT(std::abs(f(*res)), epsilon);
        }
        {
            std::optional<double> res = testCube(0.,-2.);
            ASSERT_LT(std::abs(f(*res)), epsilon);
        }
        {
            std::optional<double> res = testCube(1.,-2.);
            ASSERT_LT(std::abs(f(*res)), epsilon);
        }
        {
            std::optional<double> res = testCube(2., 0.);
            ASSERT_LT(std::abs(f(*res)), epsilon);
        }
        {
            std::optional<double> res = testCube(-2.,0.);
            ASSERT_LT(std::abs(f(*res)), epsilon);
        }
    }
    
    // x^2 - 2

    {
        auto const f = [](double x){ return x * x - 2.; };
        auto const fder = [](double x){ return 2. * x; };
        auto testCube = [f, fder] (auto a, auto b) {
            return find_root(f, fder, a, b);
        };
        {
            std::optional<double> res = testCube(-2.,0.);
            ASSERT_LT(std::abs(*res - (-std::sqrt(2.0))), 10.*epsilon);
            ASSERT_LT(std::abs(f(*res)), 10.*epsilon);
        }
        {
            std::optional<double> res = testCube(2.,0.);
            ASSERT_LT(std::abs(*res - (std::sqrt(2.0))), 10.*epsilon);
            ASSERT_LT(std::abs(f(*res)), 10.*epsilon);
        }
    }
    using namespace imajuscule::audio;
    {
        auto res = findSincCurvatureChanges(0., 10.);
        for(auto r:res) {
            std::cout << r << std::endl;
        }
    }
    constexpr auto numIteration = 10000000;
    constexpr auto numSamples = 10000000;
    using namespace std::chrono;
    using namespace profiling;
    {
        std::chrono::steady_clock::duration dt;
        double sum=0;
        std::cout << "sinc" << std::endl;
        {
            Timer<steady_clock> t(dt);
            for(double i=1; i<numIteration; ++i) {
                sum += fSinc(i);
            }
        }
        std::cout << sum << std::endl;
        std::cout << duration_cast<milliseconds>(dt).count() << " ms" << std::endl;
    }
    {
        std::chrono::steady_clock::duration dt;
        double sum=0;
        UniformSampler s(fSinc<double>,
                        0.,
                        1.,
                        numSamples);
        std::cout << "sinc resampled uniform" << std::endl;
        {
            Timer<steady_clock> t(dt);
            for(double i=1; i<numIteration; ++i) {
                sum += s.getAt(i);
            }
        }
        std::cout << sum << std::endl;
        std::cout << duration_cast<milliseconds>(dt).count() << " ms" << std::endl;
    }
    {
        std::chrono::steady_clock::duration dt;
        double sum=0;
        UniformSampler s(fSinc<double>,
                        0.,
                        1.,
                        numSamples);
        std::cout << "sinc resampled uniform opt" << std::endl;
        {
            Timer<steady_clock> t(dt);
            for(double i=1; i<numIteration; ++i) {
                sum += s.getAt(i);
            }
        }
        std::cout << sum << std::endl;
        std::cout << duration_cast<milliseconds>(dt).count() << " ms" << std::endl;
    }
    {
        std::chrono::steady_clock::duration dt;
        double sum=0;
        KaiserWindow window;
        std::cout << "kaiser" << std::endl;
        {
            Timer<steady_clock> t(dt);
            constexpr double last = numIteration;
            for(double i=1; i<last; ++i) {
                sum += window.getAt(i/last);
            }
        }
        std::cout << sum << std::endl;
        std::cout << duration_cast<milliseconds>(dt).count() << " ms" << std::endl;
    }
    {
        std::chrono::steady_clock::duration dt;
        double sum=0;
        KaiserWindow window;
        auto kaiser = [&window](double x) { return window.getAt(x); };
        UniformSampler s(kaiser,
                        0.,
                        1.,
                        numSamples);
        std::cout << "kaiser resampled uniform" << std::endl;
        {
            Timer<steady_clock> t(dt);
            constexpr double last = numIteration;
            for(double i=1; i<last; ++i) {
                sum += s.getAt(i/last);
            }
        }
        std::cout << sum << std::endl;
        std::cout << duration_cast<milliseconds>(dt).count() << " ms" << std::endl;
    }
    {
        std::chrono::steady_clock::duration dt;
        double sum=0;
        KaiserWindow window;
        auto kaiser = [&window](double x) { return window.getAt(x); };
        UniformSamplerOpt s(kaiser,
                        0.,
                        1.,
                        numSamples);
        std::cout << "kaiser resampled uniform opt" << std::endl;
        {
            Timer<steady_clock> t(dt);
            constexpr double last = numIteration;
            for(double i=1; i<last; ++i) {
                sum += s.getAt(i/last);
            }
        }
        std::cout << sum << std::endl;
        std::cout << duration_cast<milliseconds>(dt).count() << " ms" << std::endl;
    }
    {
        auto linear = [](double x) { return 2*x+1; };
        
        UniformSampler s(linear,
                         0.,
                         10.,
                         265);
        {
            auto v = s.getAt(-1.);
            ASSERT_EQ(linear(0.), v);
        }
        {
            auto v = s.getAt(0.);
            ASSERT_EQ(linear(0.), v);
        }
        {
            auto v = s.getAt(10.);
            ASSERT_EQ(linear(10.), v);
        }
        {
            auto v = s.getAt(11.);
            ASSERT_EQ(linear(10.), v);
        }
        {
            auto v = s.getAt(9.9);
            ASSERT_FLOAT_EQ(linear(9.9), v);
        }
        {
            auto v = s.getAt(0.1);
            ASSERT_FLOAT_EQ(linear(0.1), v);
        }
        {
            auto v = s.getAt(5);
            ASSERT_FLOAT_EQ(linear(5), v);
        }
    }
}
