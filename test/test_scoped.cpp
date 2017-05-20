
namespace imajuscule {
    struct TestScoped {
    protected:
        static thread_local int n, counter;
        
        void f() {
            ++counter;
        }
    };
    
    thread_local int TestScoped::n = 0;
    thread_local int TestScoped::counter = 0;
}


TEST(Scoped, multi_scopes) {
    using namespace imajuscule;
    using TestScoped = scoped::OnLeavingLast<TestScoped>;
    
    auto expected = 0;
    
    {
        TestScoped t;

        ASSERT_EQ(expected, TestScoped::counter);
    } ++ expected;
    ASSERT_EQ(expected, TestScoped::counter);
    
    {
        TestScoped t;
        {
            TestScoped t;
        } // expected doesn't change because t is in the parent scope


        ASSERT_EQ(expected, TestScoped::counter);
    } ++ expected;
    ASSERT_EQ(expected, TestScoped::counter);
}

