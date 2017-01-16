
using namespace imajuscule;

template<int Niter, int Nalloc, typename T>
void test() {
    struct C{
        
    };
    //    ASSERT_TRUE(std::is_nothrow_move_constructible<StackAllocator<int>>::va‌​lue);
    //ASSERT_TRUE(std::is_nothrow_move_constructible<C>::va‌​lue);

    auto & pool = AdaptiveStack::getInstance();
    ASSERT_EQ(0, pool.count());
    auto ptr = pool.GetNext(alignof(T), sizeof(T), 1);
    ASSERT_EQ(1, pool.count());
    ASSERT_TRUE(!!ptr);
    pool.Free(ptr, sizeof(T), 1);
    ASSERT_EQ(0, pool.count());
    
    ptr = pool.GetNext(alignof(T), sizeof(T), 1);
    ASSERT_TRUE(!!ptr);
    ASSERT_EQ(1, pool.count());
    
    std::vector<T*> v;
    v.push_back(new (ptr) T);
    v[0]->~T();
    v.clear();
    pool.Free(ptr, sizeof(T), 1);
    ASSERT_EQ(0, pool.count());

    auto n_items = 0;
    
    auto const size = Nalloc * sizeof(T);
    for(int i=0; i<Niter; i++) {
        ptr = pool.GetNext(alignof(T), size, Nalloc);
        if(size > pool.maxElemSize()) {
            ASSERT_EQ(nullptr, ptr);
        }
        else {
            ASSERT_NE(nullptr, ptr);
            n_items += Nalloc;
        }
    }
    ASSERT_EQ(n_items, pool.count());
    while(n_items) {
        pool.Free(nullptr, size, Nalloc);
        n_items -= Nalloc;
    }
    ASSERT_EQ(0, pool.count());
}

template<int Niter, int Nalloc>
void test () {
    test<Niter, Nalloc, int>();
    test<Niter, Nalloc, long int>();
    struct C {
        int a,b,c,d;
        bool e:1;
    };
    test<Niter, Nalloc, C>();
    struct D {
        int a,b,c;
    };
    test<Niter, Nalloc, D>();
    test<Niter, Nalloc, std::pair<const float, uint16_t>>();
}

TEST(AdaptiveStack, basic) {
    test<10000, 1>();
    test<5000, 2>();
    test<2000, 4>();
    test<200, 41>();
    test<20, 401>();
    test<1, 8001>();
}

TEST(Alignment, align) {
    std::aligned_storage_t<64, 256> t;
    EXPECT_EQ(256, alignof(t));
    
    constexpr auto cache_line_n_bytes = 64;
    static constexpr auto n_frames_per_buffer = cache_line_n_bytes / 4;
    static constexpr auto buffer_alignment = cache_line_n_bytes;
    using buffer_placeholder_t = std::aligned_storage_t<n_frames_per_buffer * sizeof(float), buffer_alignment>;
    struct A {
        bool:1;

        union {
            buffer_placeholder_t placeholder; // used to constrain alignment
            float buffer[n_frames_per_buffer];
        }u;
    };
    A a;
    EXPECT_EQ(0, reinterpret_cast<unsigned long>(a.u.buffer) % buffer_alignment);
    EXPECT_EQ(buffer_alignment, alignof(a.u.placeholder));
    EXPECT_EQ(buffer_alignment, alignof(buffer_placeholder_t));

}
