#pragma once

#include <vector>
#include <list>
#include <array>
#include <cassert>
#include <iostream>
#include <algorithm>

#include "traits.hpp"
#include "math.hpp"

namespace imajuscule {
    
    
    /*
     * Defines the range of acceptable load factors.
     *
     * If the load factor is outside the range, the hash table should be reduced or grown.
     *
     */
    
    struct LoadFactor {
        float min_, max_;

        LoadFactor()
        {
            auto const load_factor = 1.f;

            auto const hysteresis = 0.2f;
            
            min_ = (load_factor / 2.f) / (1.f + hysteresis);
            max_ = load_factor * (1.f + hysteresis);
            
            assert(min_ < max_);
            
            assert(min_ < load_factor / 2.f);
            assert(max_ > load_factor);
        }
    };
    
    

    /*
     * The hash function used in the hash table, based on :
     *
     * https://en.wikipedia.org/wiki/Universal_hashing
     *
     */
    
    struct UniversalHashFunction {
        
        void init(int lgSize) {
            prime = good_primes[lgSize];
            
            assert(prime >= 2); // else infinite loop
            
            // 'b' can be 0 ...
            b = rand() % prime;
            
            // ... but 'a' cannot
            do {
                a = rand() % prime;
            } while (a == 0);
        }
        
        size_t exec(int value, size_t table_size) const {
            assert(prime); // else not initialized
            return ((a*value + b) % prime) % table_size;
        }
        
    private:
        int prime = 0;
        int a, b;

        using primes = std::array<int, 31>;
        static primes good_primes;
    };

    
    
    /*
     * A hash table using chaining
     *
     */
    
    template< typename Value >
    struct HashTable {
        
        // 2^min_lg_size is the lower bound for the table size
        static constexpr auto min_lg_size = 3;
        
        HashTable(int lgSize_ = min_lg_size)
        : lgSize(std::max(lgSize_,
                          static_cast<int>(min_lg_size)))
        {
            assert(lgSize >= 0);
            assert(lgSize >= min_lg_size);
            
            auto size = pow2(lgSize);
            table_.resize(size);
            
            hash_function.init(lgSize);
        }
        
        void insert(Value value) {
            auto & chain_ = editChain(value);
            auto it = std::find( chain_.begin(), chain_.end(), value);
            if(it != chain_.end()) {
                return;
            }
            ++count;
            chain_.insert(chain_.begin(), value);
            
            grow_if_needed();
        }
        
        void remove(Value value) {
            auto & chain_ = editChain(value);
            auto it = std::find( chain_.begin(), chain_.end(), value);
            if(it == chain_.end()) {
                return;
            }
            --count;
            chain_.erase(it);

            reduce_if_needed();
        }

        
        bool has(Value value) const {
            auto & chain_ = getChain(value);
            
            auto it = std::find( chain_.begin(), chain_.end(), value);
            
            return it != chain_.end();
        }
        
        size_t size() const { return count; }
        
        bool empty() const { return count == 0; }
        
        template<typename F>
        void forEach(F f) {
            for(auto const & chain_ : table_ ) {
                for(auto v : chain_) {
                    f(v);
                }
            }
        }
        
        size_t max_chain_length() const {
            size_t len = 0;
            for(auto const & chain_ : table_ ) {
                len = std::max(len, chain_.size());
            }
            return len;
        }
        
        bool isConsistent() {
            auto load_factor_ = load_factor();
            
            if( can_reduce() && load_factor_ < min_load_factor() ) {
                std::cerr << "load_factor " << load_factor_ << " < " << min_load_factor() << std::endl;
                return false;
            }
            
            if( load_factor_ > max_load_factor() ) {
                std::cerr << "load_factor " << load_factor_ << " > " << max_load_factor() << std::endl;
                return false;
            }
            
            return true;
        }
        
    private:
        using chain = std::list<Value>;
        using table = std::vector<chain>;
        table table_;
        size_t count = 0;
        int lgSize = 0;

        UniversalHashFunction hash_function;
        
        static const LoadFactor loadfactor;
        
        chain & editChain(Value value) {
            auto h = hash(value);
            
            assert( h >= 0 );
            assert( h < table_.size() );
            
            return table_[h];
        }

        chain const & getChain(Value value) const  {
            auto h = hash(value);
            
            assert( h >= 0 );
            assert( h < table_.size() );
            
            return table_[h];
        }
        
        float load_factor() const { return ((float)count) / (float)table_.size(); }
        
        float max_load_factor() const { return loadfactor.max_; }
        float min_load_factor() const { return loadfactor.min_; }
        
        void reduce_if_needed() {
            if (can_reduce() && load_factor() < min_load_factor()) {
                reduce();
            }
        }
        
        void grow_if_needed() {
            if(load_factor() > max_load_factor()) {
                grow();
            }
        }
        
        void grow() {
            auto newLgSize = lgSize + 1;
            changeLgSize( newLgSize );
        }
        
        bool can_reduce() const {
            return lgSize > min_lg_size;
        }
        
        void reduce() {
            assert(can_reduce());
            auto newLgSize = lgSize - 1;
            changeLgSize( newLgSize );
        }
        
        void changeLgSize(int new_lg_size) {
            assert(new_lg_size != lgSize); // else logic error
            
            HashTable<Value> new_hash_table(new_lg_size);

            forEach([&new_hash_table](Value v) {
                new_hash_table.insert(v);
            });
            
            assert(new_hash_table.size() == size());

            this->operator=(std::move(new_hash_table));

            assert(load_factor() < max_load_factor());
            assert(load_factor() > min_load_factor());
        }
        
        size_t hash(Value value) const {
            int pre_hash = Traits<Value>::prehash(value);
            return hash_function.exec(pre_hash, table_.size());
        }
    };

    template<typename T>
    const LoadFactor HashTable<T>::loadfactor;
    

} // NS imajuscule

