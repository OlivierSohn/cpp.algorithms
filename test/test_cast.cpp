
using namespace imajuscule;

TEST(Cast, safe_cast) {
    {
        int i = 3;
        unsigned int j = safe_cast<unsigned int>(i);
        ASSERT_EQ(3, j);
    }
    {
        struct C_ {
            virtual ~C_() = default;
        };

        struct D: public C_ {
        };
        
        C_ c;
        
        // cast to references
        {
            auto c1 = dynamic_cast<C_&>(c);
        }
        {
            auto c1 = safe_cast<C_&>(c);
        }
        {
            ASSERT_THROW(auto c_as_D = dynamic_cast<D&>(c), std::bad_cast);
#ifndef NDEBUG
            ASSERT_THROW(auto c_as_D = safe_cast<D&>(c), std::bad_cast);
#else
            // in release, safe_cast is a safe_cast
            ASSERT_NO_THROW(auto c_as_D = safe_cast<D&>(c));
#endif
        }
        
        // cast to pointers, from null
        {
            C_*null_c = nullptr;
#ifndef NDEBUG
            EXPECT_DEBUG_ASSERT(safe_cast<C_*>(null_c));
#else
            // in release, safe_cast is a safe_cast
            auto ptr_c = safe_cast<C_*>(null_c);
            ASSERT_EQ(nullptr, ptr_c);
#endif
        }
        {
            void* null_void = nullptr;
#ifndef NDEBUG
            EXPECT_DEBUG_ASSERT(safe_cast<C_*>(null_void));
#else
            // in release, safe_cast is a safe_cast
            auto ptr_c = safe_cast<C_*>(null_void);
            ASSERT_EQ(nullptr, ptr_c);
#endif

        }
        
        // cast to pointers, from non-null
        {
            auto ptr_c = dynamic_cast<C_*>(&c);
            ASSERT_NE(nullptr, ptr_c);
        }
        {
            auto ptr_c = safe_cast<C_*>(&c);
            ASSERT_NE(nullptr, ptr_c);
        }
        
        {
            auto ptr_c = dynamic_cast<D*>(&c);
            ASSERT_EQ(nullptr, ptr_c);
        }
        {
#ifndef NDEBUG
            EXPECT_DEBUG_ASSERT(safe_cast<D*>(&c));
#else
            // in release, safe_cast is a safe_cast
            auto ptr_c = safe_cast<D*>(&c);
            ASSERT_NE(nullptr, ptr_c);
#endif
        }

    }
}
