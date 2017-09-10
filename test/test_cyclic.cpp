
using namespace imajuscule;

TEST(Cyclic, traversal) {
    constexpr auto sz = 3;
    cyclic<int> c(sz, {});
    int i = 0;
    
    // fill-in the cycle
    for(int j=0; j<sz; ++j) {
        c.feed(i++);
    }
    
    // verify elements are traversed in the right order
    for(int j=0; j<3*sz; ++j) {
        std::vector<int> fwd, bwd;
        c.for_each([&fwd](auto v) {fwd.push_back(v);});
        c.for_each_bkwd([&bwd](auto v) {bwd.push_back(v);});
        
        ASSERT_EQ(i-3, fwd[0]);
        ASSERT_EQ(i-2, fwd[1]);
        ASSERT_EQ(i-1, fwd[2]);
        
        ASSERT_EQ(i-3, bwd[2]);
        ASSERT_EQ(i-2, bwd[1]);
        ASSERT_EQ(i-1, bwd[0]);
        
        c.feed(i++);
    }
}

TEST(Cyclic, traversalLeftRight) {
    constexpr auto sz = 10;
    cyclic<int> c(sz, {});
    int i = 0;
    
    // fill-in the cycle
    for(int j=0; j<sz; ++j) {
        c.feed(i++);
    }
    
    // verify elements are traversed in the right order
    for(int j=0; j<3*sz; ++j) {
        for(int k=0; k<sz; ++k) {
            c.advance();
            
            std::vector<int> left_right;
            c.for_each_left_and_right(j, [&left_right](auto v) {
                left_right.push_back(v);
            });
            
            ASSERT_EQ(std::min(1+2*j, sz), left_right.size());
            ASSERT_EQ(*c.cycleEnd(), left_right[0]);
        }
    }
}

TEST(Cyclic, grow) {
    for(int i=0; i<10; ++i) {
        cyclic<int> c;
        ASSERT_EQ(0, c.size());
        
        for(int j=0; j<i; ++j) {
            int v=j;
            c.grow(std::move(v));
        }
        ASSERT_EQ(i, c.size());
    }
}
TEST(Cyclic, erase) {
    {
        cyclic<int> c(2);
        c.feed(1);
        c.feed(2);
        c.advance();
        c.erase(c.cycleEnd());
        ASSERT_EQ(1, *c.cycleEnd());
    }
    {
        cyclic<int> c(2);
        c.feed(1);
        c.feed(2);
        c.erase(c.cycleEnd());
        ASSERT_EQ(2, *c.cycleEnd());
    }
}
