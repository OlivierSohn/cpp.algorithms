
// we measure lookup and insertion for small number of elements
//
// we compare :
// - vector of unsorted key-value pairs
// - vector of sorted key-value pairs
// - map
// - unordered map

namespace imajuscule {
    std::ofstream myfile;

    template<typename Key, typename Value>
    struct SortedVector {
        using K = Key;
        using V = Value;
        using Pair = std::pair<K, V>;
        std::vector<Pair> v;
    };
    
    namespace stl_msr {
        template<int N, typename Key>
        struct RandomKeys {
            RandomKeys() : keys(N) {
                std::iota(keys.begin(), keys.end(), 0);
                std::shuffle(keys.begin(), keys.end(), rng::mersenne());
            }
            
            RandomKeys(RandomKeys const & o) : keys(o.keys) {
                do {
                    std::shuffle(keys.begin(), keys.end(), rng::mersenne());
                } while(keys == o.keys);
            }
            
            bool operator ==(RandomKeys const & o) const {
                return keys == o.keys;
            }
            bool operator !=(RandomKeys const & o) const {
                return keys != o.keys;
            }

            // non const versions intentionally not there
            auto begin() const { return keys.begin(); }
            auto end() const { return keys.end(); }
            
        private:
            std::vector<Key> keys;
        };
        
        namespace {
            template<typename Container, typename Key, typename Value>
            struct AddKeyValue {
                void operator()(Container & c, Key && key, Value && value) {
                    c.emplace(std::move(key), std::move(value));
                }
            };
            
            template<typename Key, typename Value>
            struct AddKeyValue<std::vector<std::pair<Key, Value>>, Key, Value> {
                void operator()(std::vector<std::pair<Key, Value>> & c, Key && key, Value && value) {
                    c.emplace_back(std::move(key), std::move(value));
                }
            };
            
            template<typename Key, typename Value>
            struct AddKeyValue<SortedVector<Key, Value>, Key, Value> {
                void operator()(SortedVector<Key, Value> & c, Key && key, Value && value) {
                    auto it = std::lower_bound(c.v.begin(), c.v.end(), std::pair<int,int>(key, value),
                                               [](auto const & v1, auto const & v2) {
                                                   return v1.first < v2.first;
                                               });
                    c.v.emplace(it, std::move(key), std::move(value));
                }
            };
            
            template<typename Container, typename Key, typename Value>
            void addKeyValue(Container & c, Key && key, Value && value) {
                AddKeyValue<Container, Key, Value> a;
                a(c, std::move(key), std::move(value));
            }
        }
        
        namespace {
            template<typename Container, typename Key, typename Value>
            struct GetValueAtKey {
                void operator()(Container & c, Key && key, Value & value) {
                    auto it = c.find(key);
                    if(it != c.end()) {
                        value = it->second;
                    }
                }
            };
            
            template<typename Key, typename Value>
            struct GetValueAtKey<std::vector<std::pair<Key, Value>>, Key, Value> {
                void operator()(std::vector<std::pair<Key, Value>> & c, Key && key, Value & value) {
                    for(auto const & p : c) {
                        if(p.first != key) {
                            continue;
                        }
                        value = p.second;
                        return;
                    }
                }
            };
            
            template<typename Key, typename Value>
            struct GetValueAtKey<SortedVector<Key, Value>, Key, Value> {
                void operator()(SortedVector<Key, Value> & c, Key && key, Value & value) {
                    auto it = std::lower_bound(c.v.begin(), c.v.end(), std::pair<int,int>(key, {}),
                                               [](auto const & v1, auto const & v2) {
                                                   return v1.first < v2.first;
                                               });
                    value = it->second;
                }
            };
            
            template<typename Container, typename Key, typename Value>
            void getValueAtKey(Container & c, Key && key, Value & value) {
                GetValueAtKey<Container, Key, Value> a;
                a(c, std::move(key), value);
            }
        }
        
        inline int lg(unsigned int v) {
            unsigned int r = 0; // r will be lg(v)
            
            while (v >>= 1) // unroll for more speed...
            {
                r++;
            }
            return r;
        }
        
        template<int N, typename Container, typename Key, typename Value>
        void DoTest() {
//            cout << "Measuring " << N << " "; COUT_TYPE(Container); cout << endl;

            static RandomKeys<N, Key> keys;
            static RandomKeys<N, Key> keys_for_retrieval(keys);
            
            if(keys == keys_for_retrieval) {
                throw;
            }
            
            auto n_iterations = 1000000 / (N*lg(N));

            std::vector<Value> sums;
            sums.reserve(n_iterations);

            std::vector<Container> vector_of_containers(n_iterations);
            
            
            clock_t t = clock();
            for(auto & c : vector_of_containers)
            {
                // fill the container up
                Value v{};
                for(auto key : keys) {
                    addKeyValue<Container, Key, Value>(c, std::move(key), v++);
                }
            }
            float time_insert = 1000000.f * (float)(clock() - t)/CLOCKS_PER_SEC;
            time_insert /= n_iterations;
            time_insert /= (float)N;
            
            t = clock();
            for(auto & c : vector_of_containers)
            {
                Value sum{};
                
                // retrieve the values by key
                for(auto key : keys_for_retrieval) {
                    Value v;
                    getValueAtKey<Container, Key, Value>(c, std::move(key), v);
                    sum += v;
                }
                
                sums.push_back(sum);
            }
            float time_get = 1000000.f * (float)(clock() - t)/CLOCKS_PER_SEC;
            time_get /= n_iterations;
            time_get /= (float)N;

            // to make sure it's not optimised away
            for(auto sum : sums) {
                if(sum != (N * (N-1))/2) {
                    throw;
                }
            }
            
//            cout<< "insert: " << time_insert << " get: " << time_get << endl;
            imajuscule::myfile <<time_insert << ","<<time_get;

        }
        
        template<int N, typename Key, typename Value>
        void TestSortedVector()
        {
            DoTest<N, SortedVector<Key, Value>, Key, Value>();
        }
        
        template<int N, typename Key, typename Value>
        void TestVector()
        {
            DoTest<N, std::vector<std::pair<Key, Value>>, Key, Value>();
        }
        
        template<int N, typename Key, typename Value>
        void TestMap()
        {
            DoTest<N, std::map<Key, Value>, Key, Value>();
        }
        
        template<int N, typename Key, typename Value>
        void TestUnorderedMap()
        {
            DoTest<N, std::unordered_map<Key, Value>, Key, Value>();
        }

        template<int N, typename Key = int, typename Value = int>
        void Test() {
            imajuscule::myfile << N << ",";
            TestSortedVector<N, Key, Value>();
            imajuscule::myfile << ",";
            TestVector<N, Key, Value>();
            imajuscule::myfile << ",";
            TestMap<N, Key, Value>();
            imajuscule::myfile << "\n";
        }
    }
}

TEST(Algorithm, stl_measure) {
    imajuscule::myfile.open ("/Users/Olivier/stl_msr.csv");
    imajuscule::myfile << "N,svec_i,svec_g,vec_i,vec_g,map_i,map_g,umap_i,umap_g,mmap_i,mmap_g\n";

    imajuscule::stl_msr::Test<5>();
    imajuscule::stl_msr::Test<10>();
    imajuscule::stl_msr::Test<20>();
    imajuscule::stl_msr::Test<40>();
    imajuscule::stl_msr::Test<80>();
    imajuscule::stl_msr::Test<160>();
    imajuscule::stl_msr::Test<380>();
    imajuscule::stl_msr::Test<760>();
    imajuscule::stl_msr::Test<1320>();
    imajuscule::stl_msr::Test<2640>();
    imajuscule::stl_msr::Test<5280>();
    
    imajuscule::myfile.close();
}


template<typename Map, typename Keys>
auto testMap(Keys const & keys, Map & map, std::vector<uint16_t> & result) {
    // sort
    uint16_t i = 0;
    for(auto k : keys) {
        map.emplace(k, i);
        ++i;
    }
    
    // traverse
    for(auto const & p : map) {
        result.push_back(p.second);
    }
    map.clear(); // included in the benchmark because I need to clear the map in the real case even if the number of items doesn't change
}

template<typename Vector>
void testVector(Vector & vec, std::vector<uint16_t> & result) {
    
    // sort
    using Pair = typename Vector::Pair;
    std::sort(vec.v.begin(), vec.v.end(), [](Pair const & p1, Pair const & p2) { return p1.first < p2.first; });
    
    // traverse
    for(auto const & p : vec.v) {
        result.push_back(p.second);
    }
}

template<int N>
void testMapSortedVector() {
    using namespace imajuscule;
    using namespace imajuscule::stl_msr;

    imajuscule::myfile << N << ",";

    std::vector<uint16_t> v0, v1, v2, v3, v4;
    v0.reserve(N);
    v1.reserve(N);
    v2.reserve(N);
    v3.reserve(N);
    v4.reserve(N);
    
    RandomKeys<N, int> keys;
    
    {
        auto t = clock();

        SortedVector<int, uint16_t> vec;
        vec.v.reserve(N);
        uint16_t i=0;
        for(auto k : keys) {
            vec.v.emplace_back(k, i);
            ++i;
        }
        
        testVector(vec, v0);
        float time = 1000000.f * (float)(clock() - t)/CLOCKS_PER_SEC;
        imajuscule::myfile << time << ",";
    }
    
    {
        std::map<int, uint16_t> map1;
        
        auto t = clock();
        testMap(keys, map1, v3);
        float time = 1000000.f * (float)(clock() - t)/CLOCKS_PER_SEC;
        imajuscule::myfile << time << ",";
    }
    
    {
        std::multimap<int, uint16_t> map1;
        
        auto t = clock();
        testMap(keys, map1, v1);
        float time = 1000000.f * (float)(clock() - t)/CLOCKS_PER_SEC;
        imajuscule::myfile << time << ",";
    }
    {
        adaptive_stack_allocated::map<int, uint16_t> map2;
        
        auto t = clock();
        testMap(keys, map2, v4);
        float time = 1000000.f * (float)(clock() - t)/CLOCKS_PER_SEC;
        imajuscule::myfile << time << ",";
    }
    
    {
        adaptive_stack_allocated::map<int, uint16_t> map2;
        
        auto t = clock();
        testMap(keys, map2, v2);
        float time = 1000000.f * (float)(clock() - t)/CLOCKS_PER_SEC;
        imajuscule::myfile << time << "\n";
    }
    
    EXPECT_EQ(v0.size(), v1.size());
    EXPECT_EQ(v1.size(), v2.size());
    EXPECT_EQ(v2.size(), v3.size());
    EXPECT_EQ(v3.size(), v4.size());
    
    for(int i=0; i<v0.size(); ++i) {
        EXPECT_EQ(v0[i], v1[i]);
        EXPECT_EQ(v1[i], v2[i]);
        EXPECT_EQ(v2[i], v3[i]);
        EXPECT_EQ(v3[i], v4[i]);
    }
}

/*
 * This test shows that a sorted vector is way more efficient than multimaps or maps
 * when the number of items don't change (no need to insert in the middle)
 */
TEST(Comparison, map_sortedvector) {
    
    imajuscule::myfile.open ("/Users/Olivier/map_sortedvec.csv");
    imajuscule::myfile << "N,svec,map,mmap,mapa,mmapa\n";

    testMapSortedVector<10>();
    testMapSortedVector<100>();
    testMapSortedVector<1000>();
    testMapSortedVector<10000>();
    testMapSortedVector<100000>();
    
    imajuscule::myfile.close();

}
