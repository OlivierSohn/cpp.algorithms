/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {
    
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
    struct StackVector {
    private:
        using assert_unchanged_capacity = AssertUnchangedCapacity<adaptive_stack_allocated::vector<T>>;
        adaptive_stack_allocated::vector<T> v;
    public:
        // non copyable
        StackVector(const StackVector &) = delete;
        StackVector & operator=(const StackVector&) = delete;
        
        // movable
        StackVector(StackVector &&) = default;
        StackVector& operator =(StackVector &&) = default;
        
        StackVector() = default;
        
        StackVector(size_t size) noexcept {
            reserve(size);
        }
        
        StackVector(size_t size, T val) noexcept {
            reserve(size);
            resize(size);
            fill(val);
        }
        
        StackVector( std::initializer_list<T> init ) noexcept :
        v(std::move(init)) {
#ifndef NDEBUG
            cpgo.acquire();
#endif
        }
        
        template<typename Iterator>
        StackVector( Iterator first, Iterator last ) noexcept :
        v(std::move(first), std::move(last)) {
#ifndef NDEBUG
            cpgo.acquire();
#endif
        }
        
        ~StackVector() {
            // first deallocate children (important for StackVector<StackVector<T>>)
            {
                // in this operation we don't want our vector to change its allocation ...
                // ... unfortunately the standard doesn't guarantee that clear doesnt reallocate.
                // The assert will tell us if our assumption is ok.
                assert_unchanged_capacity c(v);
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
        
        T & operator[](size_t n) { return v[n]; };
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
        StackGrowthControl cpgo;
#endif
        
        template <class... Args>
        void fill(Args&&... args) {
            std::fill(v.begin(), v.end(), std::forward<Args>(args)... );
        }
    };
    
    // todo make this work:
    //static_assert(std::is_nothrow_move_constructible<StackVector<int>>::va‌​lue, "");
    //static_assert(std::is_nothrow_move_constructible<StackGrowthControl>::va‌​lue, "");
    //static_assert(std::is_nothrow_move_constructible<adaptive_stack_allocated::vector<int>>::va‌​lue, "");
//    static_assert(std::is_nothrow_move_constructible<std::vector<int>>::va‌​lue, "");
    //static_assert(std::is_nothrow_move_constructible<int>::va‌​lue, "");
    
}
