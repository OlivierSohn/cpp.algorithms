namespace imajuscule {

std::set<std::vector<Scaling>> mkScalings(int const firstSz,
                                          std::vector<std::vector<int>> const & repeats) {
    std::set<std::vector<Scaling>> res;
    for(auto const & rs : repeats) {
        std::vector<Scaling> v;
        int sz = firstSz;
        for(auto r : rs) {
            if(r) {
                v.emplace_back(sz, r);
            }
            sz *= 2;
        }
        res.emplace(v);
    }
    return res;
}

template<typename T>
void analyzeSetDifferences(std::set<T> const & s1,
                           std::set<T> const & s2) {
    if(s1 == s2) {
        return;
    }
    for(auto const & e1 : s1) {
        if(s2.count(e1) == 0) {
            std::cout << "diff ";
        }
    }
    for(auto const & e2 : s2) {
        if(s1.count(e2) == 0) {
            std::cout << "diff ";
        }
    }
}

}

TEST(ConvolutionScale, iterateScales_noLastSize) {
    using namespace imajuscule;
    
    // verify that 0 coeff leads to a single empty scaling.
    {
        ScalingsIterator it{1,0};
        int count = 0;
        it.forEachScaling([&count](auto const & v){
            ASSERT_TRUE(v.empty());
            ++count;
        });
        ASSERT_EQ(1, count);
    }

    // verify that 1 coeff leads to a single scaling.
    for(int firstSz = 1; firstSz < 10; ++firstSz)
    {
        ScalingsIterator it{firstSz,1};
        int count = 0;
        it.forEachScaling([&count, firstSz](auto const & v){
            std::vector<Scaling> v2 {{firstSz,1}};
            ASSERT_EQ(v2, v);
            ++count;
        });
        ASSERT_EQ(1, count);
    }

    // verify that less coeffs than first size leads to a single scaling.
    {
        int const firstSz = 10;
        for(int nCoeffs = 0; nCoeffs <= firstSz; ++nCoeffs)
        {
            ScalingsIterator it{firstSz,nCoeffs};
            int count = 0;
            it.forEachScaling([&count, firstSz, nCoeffs](auto const & v){
                if(nCoeffs) {
                    std::vector<Scaling> v2 {{firstSz,1}};
                    ASSERT_EQ(v2, v);
                }
                else {
                    ASSERT_TRUE(v.empty());
                }
                ++count;
            });
            ASSERT_EQ(1, count);
        }
    }
    
    // verify that coeffs between firstsize+1 and 2*firstsize lead to a single scaling.
    {
        int const firstSz = 10;
        for(int nCoeffs = firstSz+1; nCoeffs <= 2*firstSz; ++nCoeffs)
        {
            ScalingsIterator it{firstSz,nCoeffs};
            int count = 0;
            it.forEachScaling([&count, firstSz](auto const & v){
                std::vector<Scaling> v2 {{firstSz,2}};
                ASSERT_EQ(v2, v);
                ++count;
            });
            ASSERT_EQ(1, count);
        }
    }
    
    
    // verify that coeffs between 2*firstsize+1 and 3*firstsize lead to 2 scalings.
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
            {3},
            {1, 1}
        });
        for(int nCoeffs = 2*firstSz+1; nCoeffs <= 3*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }

    // verify that coeffs between 3*firstsize+1 and 4*firstsize lead to 2 scalings.
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
            {4},
            {1, 2}
        });
        for(int nCoeffs = 3*firstSz+1; nCoeffs <= 4*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }
    // verify that coeffs between 4*firstsize+1 and 5*firstsize lead to 2 scalings.
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
            {5},
            {1, 2}
        });
        for(int nCoeffs = 4*firstSz+1; nCoeffs <= 5*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }
    // verify that coeffs between 5*firstsize+1 and 6*firstsize lead to 4 scalings.
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
            {6},
            {1, 3},
            {1, 1, 1},
            {3, 0, 1}
        });
        for(int nCoeffs = 5*firstSz+1; nCoeffs <= 6*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
            {7},
            {1, 3},
            {1, 1, 1},
            {3, 0, 1}
        });
        for(int nCoeffs = 6*firstSz+1; nCoeffs <= 7*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
            {8},
            {1, 4},
            {1, 1, 2},
            {3, 0, 2}
        });
        for(int nCoeffs = 7*firstSz+1; nCoeffs <= 8*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
            {9},
            {1, 4},
            {1, 1, 2},
            {3, 0, 2}
        });
        for(int nCoeffs = 8*firstSz+1; nCoeffs <= 9*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
            {10},
            {1, 5},
            {1, 1, 2},
            {3, 0, 2}
        });
        for(int nCoeffs = 9*firstSz+1; nCoeffs <= 10*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
            {11},
            {1, 5},
            {1, 1, 2},
            {3, 0, 2}
        });
        for(int nCoeffs = 10*firstSz+1; nCoeffs <= 11*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
            {12},
            {7, 0, 0, 1},
            {3, 0, 3},
            {3, 0, 1, 1},
            {1, 6},
            {1, 3, 0, 1},
            {1, 1, 3},
            {1, 1, 1, 1}
        });
        for(int nCoeffs = 11*firstSz+1; nCoeffs <= 12*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            analyzeSetDifferences(valids, results);
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
            {13},
            {7, 0, 0, 1},
            {3, 0, 3},
            {3, 0, 1, 1},
            {1, 6},
            {1, 3, 0, 1},
            {1, 1, 3},
            {1, 1, 1, 1}
        });
        for(int nCoeffs = 12*firstSz+1; nCoeffs <= 13*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            analyzeSetDifferences(valids, results);
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
            {14},
            {7, 0, 0, 1},
            {3, 0, 3},
            {3, 0, 1, 1},
            {1, 7},
            {1, 3, 0, 1},
            {1, 1, 3},
            {1, 1, 1, 1}
        });
        for(int nCoeffs = 13*firstSz+1; nCoeffs <= 14*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            analyzeSetDifferences(valids, results);
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
            {15},
            {7, 0, 0, 1},
            {3, 0, 3},
            {3, 0, 1, 1},
            {1, 7},
            {1, 3, 0, 1},
            {1, 1, 3},
            {1, 1, 1, 1}
        });
        for(int nCoeffs = 14*firstSz+1; nCoeffs <= 15*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            analyzeSetDifferences(valids, results);
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
            {16},
            {7, 0, 0, 2},
            {3, 0, 4},
            {3, 0, 1, 2},
            {1, 8},
            {1, 3, 0, 2},
            {1, 1, 4},
            {1, 1, 1, 2}
        });
        for(int nCoeffs = 15*firstSz+1; nCoeffs <= 16*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            analyzeSetDifferences(valids, results);
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
            {17},
            {7, 0, 0, 2},
            {3, 0, 4},
            {3, 0, 1, 2},
            {1, 8},
            {1, 3, 0, 2},
            {1, 1, 4},
            {1, 1, 1, 2}
        });
        for(int nCoeffs = 16*firstSz+1; nCoeffs <= 17*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            analyzeSetDifferences(valids, results);
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }
}


TEST(ConvolutionScale, iterateScales_lastSize) {
    using namespace imajuscule;
    
    // verify that 0 coeff leads to a single empty scaling.
    {
        ScalingsIterator it{1,0,2};
        it.forEachScaling([](auto const & v){
            ASSERT_TRUE(false);
        });
    }
    {
        ScalingsIterator it{1,0,1};
        it.forEachScaling([](auto const & v){
            ASSERT_TRUE(false);
        });
    }

    // verify that 1 coeff leads to a single scaling.
    for(int firstSz = 1; firstSz < 10; ++firstSz)
    {
        ScalingsIterator it{firstSz,1,firstSz};
        int count = 0;
        it.forEachScaling([&count, firstSz](auto const & v){
            std::vector<Scaling> v2 {{firstSz,1}};
            ASSERT_EQ(v2, v);
            ++count;
        });
        ASSERT_EQ(1, count);
    }
    for(int firstSz = 1; firstSz < 10; ++firstSz)
    {
        ScalingsIterator it{firstSz,1,2*firstSz};
        it.forEachScaling([](auto const & v){
            ASSERT_TRUE(false);
        });
    }

    // verify that less coeffs than first size leads to a single scaling.
    {
        int const firstSz = 10;
        for(int nCoeffs = 0; nCoeffs <= firstSz; ++nCoeffs)
        {
            ScalingsIterator it{firstSz,nCoeffs,firstSz};
            int count = 0;
            it.forEachScaling([&count, firstSz, nCoeffs](auto const & v){
                if(nCoeffs) {
                    std::vector<Scaling> v2 {{firstSz,1}};
                    ASSERT_EQ(v2, v);
                }
                else {
                    ASSERT_TRUE(v.empty());
                }
                ++count;
            });
            ASSERT_EQ(nCoeffs?1:0, count);
        }
        for(int nCoeffs = 0; nCoeffs <= firstSz; ++nCoeffs)
        {
            ScalingsIterator it{firstSz,nCoeffs, 2*firstSz};
            it.forEachScaling([](auto const & v){
                ASSERT_TRUE(false);
            });
        }
    }

    // verify that coeffs between firstsize+1 and 2*firstsize lead to a single scaling.
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
            {2}
        });
        for(int nCoeffs = firstSz+1; nCoeffs <= 2*firstSz; ++nCoeffs)
        {
            ScalingsIterator it{firstSz,nCoeffs,firstSz};
            int count = 0;
            std::set<std::vector<Scaling>> results;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
            {1, 1}
        });
        for(int nCoeffs = firstSz+1; nCoeffs <= 2*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs,2*firstSz};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
        });
        for(int nCoeffs = firstSz+1; nCoeffs <= 2*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs,4*firstSz};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }
    
    // verify that coeffs between 2*firstsize+1 and 3*firstsize lead to 2 scalings.
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
            {3}
        });
        for(int nCoeffs = 2*firstSz+1; nCoeffs <= 3*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs,firstSz};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
            {1, 1}
        });
        for(int nCoeffs = 2*firstSz+1; nCoeffs <= 3*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs,2*firstSz};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
        });
        for(int nCoeffs = 2*firstSz+1; nCoeffs <= 3*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs,4*firstSz};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }

    // verify that coeffs between 3*firstsize+1 and 4*firstsize lead to 2 scalings.
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
            {4}
        });
        for(int nCoeffs = 3*firstSz+1; nCoeffs <= 4*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs,firstSz};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
            {1, 2}
        });
        for(int nCoeffs = 3*firstSz+1; nCoeffs <= 4*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs,2*firstSz};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
            {1, 1, 1},
            {3, 0, 1}
        });
        for(int nCoeffs = 3*firstSz+1; nCoeffs <= 4*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs,4*firstSz};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
        });
        for(int nCoeffs = 3*firstSz+1; nCoeffs <= 4*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs,8*firstSz};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }

    // verify that coeffs between 4*firstsize+1 and 5*firstsize lead to 2 scalings.
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
            {5}
        });
        for(int nCoeffs = 4*firstSz+1; nCoeffs <= 5*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs,firstSz};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
            {1, 2}
        });
        for(int nCoeffs = 4*firstSz+1; nCoeffs <= 5*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs,2*firstSz};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
            {1, 1, 1},
            {3, 0, 1}
        });
        for(int nCoeffs = 4*firstSz+1; nCoeffs <= 5*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs,4*firstSz};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
        });
        for(int nCoeffs = 4*firstSz+1; nCoeffs <= 5*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs,8*firstSz};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }

    // verify that coeffs between 5*firstsize+1 and 6*firstsize lead to 4 scalings.
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
            {6}
        });
        for(int nCoeffs = 5*firstSz+1; nCoeffs <= 6*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs,firstSz};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
            {1, 3}
        });
        for(int nCoeffs = 5*firstSz+1; nCoeffs <= 6*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs,2*firstSz};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
            {1, 1, 1},
            {3, 0, 1}
        });
        for(int nCoeffs = 5*firstSz+1; nCoeffs <= 6*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs,4*firstSz};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }
    {
        int const firstSz = 10;
        auto valids = mkScalings(firstSz, {
        });
        for(int nCoeffs = 5*firstSz+1; nCoeffs <= 6*firstSz; ++nCoeffs)
        {
            std::set<std::vector<Scaling>> results;
            ScalingsIterator it{firstSz,nCoeffs,8*firstSz};
            int count = 0;
            it.forEachScaling([&count, &results](auto const & v){
                results.emplace(v);
                ++count;
            });
            ASSERT_EQ(valids.size(), count);
            ASSERT_EQ(results, valids);
        }
    }
}

TEST(ConvolutionScale, simulateBatch) {

    using namespace imajuscule;

    using C = CustomScaleConvolution<FFTConvolutionIntermediate < PartitionnedFFTConvolutionCRTP<double, fft::Fastest> >>;

    int const firstSz = 4;
    int const nCoeffs = 4545;
    ScalingsIterator it{
        firstSz,
        nCoeffs
    };
    it.forEachScaling([nCoeffs](auto const & v){
        auto sim = mkSimulation<C>(v, nCoeffs);
        auto const batchSize = sim.getBiggestScale();
        double cost{};
        for(int i=0; i<batchSize; ++i) {
            cost += sim.simuStep();
        }
        double batchCost = sim.simuBatch(batchSize);
        ASSERT_NEAR(cost, batchCost, 1e-3);
    });
}
