
using namespace imajuscule;

namespace imajuscule {
    namespace freelist {
        template<typename FREELIST>
        void test() {
            auto size = FREELIST::size;
            std::vector<typename FREELIST::value_type*> values;
            
            FREELIST l;
            for(size_t i=0; i<size; ++i) {
                auto p = l.Take();
                ASSERT_NE(nullptr, p);
                values.push_back(p);
            }
            // when the free list is full we cannot take more :
            ASSERT_EQ(nullptr, l.Take());
            
            
            // when the free list is full and one element is returned,
            // the next element taken is the returned element :
            for(int i=0; i<values.size(); ++i) {
                auto ret = values[i];
                l.Return(ret);
                auto p = l.Take();
                ASSERT_EQ(ret, p);
            }
        }
    }
}

TEST(FreeList, simple) {
    using namespace imajuscule::freelist;
    constexpr auto size = 4;
    
    test<FreeList<double, size, void*>>();
    test<FreeList<uint16_t, size, uint16_t>>();
    test<FreeList<int, size, uint16_t>>();
}


