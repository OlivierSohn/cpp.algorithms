
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
