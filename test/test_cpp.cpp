
namespace imajuscule {
    namespace test {
        template<typename I>
        static void f(I&& i) {
            ++i;
        }
    }
}
TEST(UniversalReference, pod) {
    using namespace imajuscule::test;
    int i=0;
    imajuscule::test::f(i);
    ASSERT_EQ(1, i);
}
