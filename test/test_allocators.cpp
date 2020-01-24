
using namespace imajuscule;

template<int Niter, int Nalloc, typename T>
void testAdaptiveStack() {
    struct C{
        
    };
    //    ASSERT_TRUE(std::is_nothrow_move_constructible<StackAllocator<int>>::va‌​lue);
    //ASSERT_TRUE(std::is_nothrow_move_constructible<C>::va‌​lue);

    auto & pool = AdaptiveStack::getThreadLocalInstance();
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
void testAdaptiveStack () {
    testAdaptiveStack<Niter, Nalloc, int>();
    testAdaptiveStack<Niter, Nalloc, long int>();
    struct C {
        int a,b,c,d;
        bool e:1;
    };
    testAdaptiveStack<Niter, Nalloc, C>();
    struct D {
        int a,b,c;
    };
    testAdaptiveStack<Niter, Nalloc, D>();
    testAdaptiveStack<Niter, Nalloc, std::pair<const float, uint16_t>>();
}

void testAdaptiveStack(void) {
    testAdaptiveStack<10000, 1>();
    testAdaptiveStack<5000, 2>();
    testAdaptiveStack<2000, 4>();
    testAdaptiveStack<200, 41>();
    testAdaptiveStack<20, 401>();
    testAdaptiveStack<1, 8001>();
}

TEST(AdaptiveStack, single_thread) {
    testAdaptiveStack();
}

TEST(AdaptiveStack, multi_thread) {
    std::vector<std::thread> threads;
    threads.reserve(100);
    for(int i=0; i<100; ++i) {
        threads.emplace_back([](){
            testAdaptiveStack();
        });
    }
    for(auto & t : threads) {
        t.join();
    }
}

TEST(AdaptiveStack, thread_local_stuff) {
    std::vector<std::thread> threads;
    std::vector<void*> locations;
    threads.reserve(100);
    locations.reserve(100);
    for(int i=0; i<100; ++i) {
        threads.emplace_back([i, &locations] () {
            locations[i] = &AdaptiveStack::getThreadLocalInstance();
        });
    }
    for(auto & t : threads) {
        t.join();
    }
    std::set<void*> s(locations.begin(), locations.end());
    ASSERT_EQ(locations.size(), s.size()); // unicity of locations
}

TEST(Alignment, align) {
    std::aligned_storage_t<64, 256> t;
    EXPECT_EQ(256, alignof(t));
    
    constexpr auto cache_line_n_bytes = 64;
    constexpr auto n_frames_per_buffer = cache_line_n_bytes / 4;
    constexpr auto buffer_alignment = cache_line_n_bytes;
    using buffer_placeholder_t = std::aligned_storage_t<n_frames_per_buffer * sizeof(float), buffer_alignment>;
    struct A {
        bool:1;

        union {
            buffer_placeholder_t placeholder; // used to constrain alignment
            float buffer[n_frames_per_buffer];
        }u;
    };
    EXPECT_EQ(buffer_alignment, alignof(buffer_placeholder_t));

    A a;
    EXPECT_EQ(0, reinterpret_cast<unsigned long>(a.u.buffer) % buffer_alignment);
    EXPECT_EQ(buffer_alignment, alignof(a.u.placeholder));
    std::vector<std::unique_ptr<A>> v;
    for(int i=0; i<100; i++) {
        auto p = std::make_unique<A>();
        EXPECT_TRUE(&p->u.buffer[0] == p->u.buffer);
        EXPECT_TRUE(static_cast<void*>(&p->u.placeholder) == static_cast<void*>(p->u.buffer));
        //std::cout << reinterpret_cast<unsigned long>(&p->u.buffer[0]) % buffer_alignment << std::endl;
        v.push_back(std::move(p));
    }
    auto p = std::make_unique<A>();
    EXPECT_TRUE(&p->u.buffer[0] == p->u.buffer);
    EXPECT_TRUE(static_cast<void*>(&p->u.placeholder) == static_cast<void*>(p->u.buffer));
    EXPECT_EQ(0, reinterpret_cast<unsigned long>(&p->u.buffer[0]) % buffer_alignment);
    EXPECT_EQ(buffer_alignment, alignof(p->u.placeholder));
    EXPECT_EQ(buffer_alignment, alignof(A));
    
}

TEST(AlignedAllocator, alignment) {
    AlignedAllocator<int, Alignment::CACHE_LINE> aligned_allocator;
    auto mem = aligned_allocator.allocate(1);
    EXPECT_EQ(0, reinterpret_cast<unsigned long>(mem) % to_underlying(Alignment::CACHE_LINE));
    
}

TEST(MonotonicAllocator, adresses)
{
    using namespace imajuscule;
    
    using base_allocator = AlignedAllocator<double, Alignment::PAGE>;

    using vector =
    monotonic::vector<base_allocator>;
    
    using monotonic_buffer =
    MonotonicBuffer<base_allocator>;
    
    using use_monotonic_buffer =
    UseMonotonicBuffer<base_allocator>;
    
    {
        monotonic_buffer mb(2);

        vector v1;
        ASSERT_THROW(v1.resize(1),
                     std::runtime_error);
        
        {
            use_monotonic_buffer u(mb);
            v1.resize(1);
        }
        vector v2;
        ASSERT_THROW(v2.resize(1),
                     std::runtime_error);
        {
            use_monotonic_buffer u(mb);
            v2.resize(1);
            
            ASSERT_THROW(v2.resize(2),
                         std::runtime_error);
        }
        ASSERT_EQ(&v1[0] + 1,
                  &v2[0]);
        ASSERT_EQ(0,
                  mb.remaining());
    }
    
    
    {
        monotonic_buffer mbA(2);
        monotonic_buffer mbB(2);

        vector v1, v2, v3, v4;
        
        {
            use_monotonic_buffer u(mbA);
            v1.resize(1);
            {
                use_monotonic_buffer u(mbB);
                v3.resize(1);
                v4.resize(1);
            }
            v2.resize(1);
        }
        ASSERT_EQ(&v1[0] + 1,
                  &v2[0]);
        ASSERT_EQ(&v3[0] + 1,
                  &v4[0]);
        ASSERT_EQ(0,
                  mbA.remaining());
        ASSERT_EQ(0,
                  mbB.remaining());
    }
}
