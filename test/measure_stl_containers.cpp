#include "gtest/gtest.h"

#include <map>
#include <unordered_map>
#include <vector>
#include <random>
#include <ctime>

#include "print_type.hpp"
#include <iostream>
#include <fstream>

using namespace std;

// we measure lookup and insertion for small number of elements
//
// we compare :
// - vector of unsorted key-value pairs
// - vector of sorted key-value pairs
// - map
// - unordered map

namespace imj {
    random_device rnd_device;
    mt19937 mersenne_engine(rnd_device());
    ofstream myfile;

    template<typename Key, typename Value>
    struct SortedVector {
        std::vector<std::pair<Key, Value>> v;
    };
    
    namespace stl_msr {
        template<int N, typename Key>
        struct RandomKeys {
            RandomKeys() : keys(N) {
                std::iota(keys.begin(), keys.end(), 0);
                std::shuffle(keys.begin(), keys.end(), mersenne_engine);
            }
            
            RandomKeys(RandomKeys const & o) : keys(o.keys) {
                std::shuffle(keys.begin(), keys.end(), mersenne_engine);
            }
            
            bool operator ==(RandomKeys const & o) const {
                return keys == o.keys;
            }
            bool operator !=(RandomKeys const & o) const {
                return keys != o.keys;
            }
            
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
            imj::myfile <<time_insert << ","<<time_get;

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
            imj::myfile << N << ",";
            TestSortedVector<N, Key, Value>();
            imj::myfile << ",";
            TestVector<N, Key, Value>();
            imj::myfile << ",";
            TestMap<N, Key, Value>();
            imj::myfile << ",";
            TestMyHasTable<N, Key, Value>();
            imj::myfile << "\n";
        }
    }
}

TEST(Algorithm, stl_measure) {
    imj::myfile.open ("/Users/Olivier/stl_msr.csv");
    imj::myfile << "N,svec_i,svec_g,vec_i,vec_g,map_i,map_g,umap_i,umap_g,mmap_i,mmap_g\n";

    imj::stl_msr::Test<5>();
    imj::stl_msr::Test<10>();
    imj::stl_msr::Test<20>();
    imj::stl_msr::Test<40>();
    imj::stl_msr::Test<80>();
    imj::stl_msr::Test<160>();
    imj::stl_msr::Test<380>();
    imj::stl_msr::Test<760>();
    imj::stl_msr::Test<1320>();
    imj::stl_msr::Test<2640>();
    imj::stl_msr::Test<5280>();
    imj::myfile.close();
}
