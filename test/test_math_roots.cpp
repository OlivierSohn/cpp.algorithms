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
}
