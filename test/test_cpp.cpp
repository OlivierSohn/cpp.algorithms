
namespace imajuscule {
    namespace test_universal_ref {
        template<typename I>
        static void increment(I&& i) {
            ++i;
        }
    }
}
TEST(UniversalReference, pod) {
    int i=0;
    imajuscule::test_universal_ref::increment(i);
    ASSERT_EQ(1, i);
}

