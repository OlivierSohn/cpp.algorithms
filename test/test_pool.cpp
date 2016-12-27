
using namespace imajuscule;

template<int Niter, int Nalloc, typename T>
void test() {
    
    auto & pool = Pool::getInstance();
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

TEST(Pool, basic) {
    test<10000, 1>();
    test<5000, 2>();
    test<2000, 4>();
    test<200, 41>();
    test<20, 401>();
    test<1, 8001>();
}
