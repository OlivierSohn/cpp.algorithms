
using namespace imajuscule;

namespace imajuscule {
    namespace freelist {
        template<typename FREELIST>
        void testNoSublists() {
            auto size = FREELIST::sizeBlock;
            std::vector<typename FREELIST::value_type*> values;

            FREELIST l;
            for(size_t i=0; i<size; ++i) {
                auto p = l.Take();
                ASSERT_NE(nullptr, p);
                values.push_back(p);
            }
            // the main free list is full, the next call to Take would allocate a sub freeList.

            // when the main free list is full and one element is returned,
            // the next element taken is the returned element :
            for(int i=0; i<values.size(); ++i) {
                auto ret = values[i];
                l.Return(ret);
                auto p = l.Take();
                ASSERT_EQ(ret, p);
            }
        }

        template<typename FREELIST>
        void testSublists() {
            auto size = FREELIST::sizeBlock;
            std::vector<typename FREELIST::value_type*> values;

            FREELIST l;
            for(size_t i=0; i<4*size; ++i) {
                auto p = l.Take();
                ASSERT_NE(nullptr, p);
                values.push_back(p);
            }

            auto prevValues = values;


            // return all values, in random order:
            Shuffle(values);
            for(int i=0; i<values.size(); ++i) {
                l.Return(values[i]);
            }

            // by now the free list is empty.

            std::vector<typename FREELIST::value_type*> newValues;

            for(size_t i=0; i<4*size; ++i) {
                newValues.push_back(l.Take());
            }
            
            StdSort(values);
            StdSort(newValues);

            // verify that taking the values again will give the same results (modulo sorting):
            EXPECT_EQ(values, newValues);
        }
    }
}

TEST(FreeList, simple) {
    using namespace imajuscule::freelist;
    constexpr auto size = 400;

    testNoSublists<FreeList<double, size, void*>>();
    testNoSublists<FreeList<uint16_t, size, uint16_t>>();
    testNoSublists<FreeList<int, size, uint16_t>>();
    testSublists<FreeList<double, size, void*>>();
    testSublists<FreeList<uint16_t, size, uint16_t>>();
    testSublists<FreeList<int, size, uint16_t>>();
}
