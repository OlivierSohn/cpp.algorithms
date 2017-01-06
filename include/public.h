#pragma once

#include <cassert>
#include <cstddef>
#include <iostream>
#include <memory>
#include <vector>
#include <map>
#include <unordered_map>

#include "safe_cast.hpp"
#include "print_type.hpp"
#include "pool.h"
#include "allocator.hpp"

namespace imajuscule {
    namespace pool {
        template<typename T>
        using vector = std::vector<T, Allocator<T>>;
        
        template<typename K, typename V, class Compare = std::less<K>>
        using map = std::map<K, V, Compare, Allocator<std::pair<const K,V>>>;
        
        template<typename K, typename V, class Hash = std::hash<K>, class Pred = std::equal_to<K>>
        using unordered_map = std::unordered_map<K, V, Hash, Pred, Allocator<std::pair<const K,V>>>;
    }
    
    template<typename T>
    struct AssertUnchangedCapacity {
        AssertUnchangedCapacity(T const & container)
#ifndef NDEBUG
        : container(container)
#endif
        {
#ifndef NDEBUG
            capacity = container.capacity();
#endif
        }
        ~AssertUnchangedCapacity() {
#ifndef NDEBUG
            assert(capacity == container.capacity());
#endif
        }
    private:
#ifndef NDEBUG
        size_t capacity;
        T const & container;
#endif
    };

    template<typename T>
    struct StaticVector {
    private:
        using assert_unchanged_capacity = AssertUnchangedCapacity<pool::vector<T>>;
        pool::vector<T> v;
    public:
        // non copyable
        StaticVector(const StaticVector &) = delete;
        StaticVector & operator=(const StaticVector&) = delete;
        
        // movable
        StaticVector(StaticVector &&) = default;
        StaticVector& operator =(StaticVector &&) = default;

        StaticVector() = default;
        
        StaticVector(size_t size)  {
            reserve(size);
        }
        
        StaticVector( std::initializer_list<T> init ) :
        v(std::move(init)) {
#ifndef NDEBUG
            cpgo.acquire();
#endif
        }
        
        template<typename Iterator>
        StaticVector( Iterator first, Iterator last ) :
        v(std::move(first), std::move(last)) {
#ifndef NDEBUG
            cpgo.acquire();
#endif
        }
        
        ~StaticVector() {
            // first deallocate children (important for StaticVector<StaticVector<T>>)
            {
                // in this operation we don't want our vector to change its allocation ...
                assert_unchanged_capacity c(v);
                // ... unfortunately the standard doesn't guarantee that clear doesnt reallocate.
                v.clear();
            }
        }
        
        void reserve(size_t size) {
            v.reserve(size);
#ifndef NDEBUG
            cpgo.acquire();
#endif
        }
        
        void resize(size_t size) {
            assert_unchanged_capacity c(v);
            v.resize(size);
        }
        
        template<typename Iterator>
        auto erase(Iterator it) -> decltype(v.erase(it)) {
            assert_unchanged_capacity c(v);
            return v.erase(std::move(it));
        }
        
        auto begin() const -> decltype(v.begin()) { return v.begin(); }
        auto end() const -> decltype(v.end()) { return v.end(); }
        auto begin() -> decltype(v.begin()) { return v.begin(); }
        auto end() -> decltype(v.end()) { return v.end(); }
        auto rbegin() const -> decltype(v.rbegin()) { return v.rbegin(); }
        auto rend() const -> decltype(v.rend()) { return v.rend(); }
        auto rbegin() -> decltype(v.rbegin()) { return v.rbegin(); }
        auto rend() -> decltype(v.rend()) { return v.rend(); }
        
        size_t size() const { return v.size(); }
        size_t capacity() const { return v.capacity(); }
        bool empty() const { return v.empty(); }

        void clear() {
            assert_unchanged_capacity c(v);
            return v.clear();
        }
        
        void pop_back() { return v.pop_back(); }
        
        auto back() -> decltype(v.back()) { return v.back(); }
        auto back() const -> decltype(v.back()) { return v.back(); }
        
        auto front() -> decltype(v.front()) { return v.front(); }
        auto front() const -> decltype(v.front()) { return v.front(); }
        
        auto operator[](size_t n) -> decltype(v[n]) { return v[n]; };
        auto operator[](size_t n) const -> decltype(v[n]) { return v[n]; };

        template <class... Args>
        void emplace_back(Args&&... args) {
            assert_unchanged_capacity c(v);
            v.emplace_back( std::forward<Args>(args)... );
        }

        template <typename Iterator, class... Args>
        void emplace(Iterator it, Args&&... args) {
            assert_unchanged_capacity c(v);
            v.emplace(std::move(it),
                      std::forward<Args>(args)... );
        }

        void push_back(T const & val) {
            assert_unchanged_capacity c(v);
            v.push_back(val);
        }
        void push_back(T && val) {
            assert_unchanged_capacity c(v);
            v.push_back(std::move(val));
        }
        
        auto data() -> decltype(v.data()) { return v.data(); }
        
    private:
#ifndef NDEBUG
        ControlledPoolGrowthObject cpgo;
#endif
    };
}
