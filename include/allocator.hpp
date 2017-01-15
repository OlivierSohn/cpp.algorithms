
namespace imajuscule {
    
    template <class T> class StackAllocator;
    
    template <> class StackAllocator<void>
    {
    public:
        typedef void * pointer;
        typedef const void* const_pointer;
        typedef void value_type;
        template <class U> struct rebind { typedef StackAllocator<U> other; };
    };
    
    template <class T> class StackAllocator
    {
    public:
        typedef T value_type;
        typedef size_t size_type;
        typedef ptrdiff_t difference_type;
        typedef T* pointer;
        typedef const T* const_pointer;
        typedef T& reference;
        typedef const T& const_reference;
        
        template <class U> struct rebind { typedef StackAllocator<U> other; };
        
        StackAllocator() noexcept
        : pool(Pool::getInstance()) {};
        
        StackAllocator( const StackAllocator& other ) noexcept
        : StackAllocator() {};
        
        template <class U>
        StackAllocator(const StackAllocator<U>&) noexcept
        : StackAllocator() {};
        
        size_type max_size() const noexcept {
            return size_type(~0) / sizeof(T);
        };
        
        pointer
        address(reference r) const noexcept
        { return std::addressof(r); };
        
        const_pointer
        address(const_reference r) const noexcept
        { return std::addressof(r); };
        
        pointer
        allocate( size_type n, StackAllocator<void>::const_pointer hint = 0 )
        {
            auto const size_request = n*sizeof(T);
            pointer return_value = reinterpret_cast<pointer>( pool.GetNext(alignof(T), size_request, n) );
            if ( return_value != 0 ) {
                return return_value;
            }
            if(!pool.used()) {
                pool.resize(2*size_request);
                return_value = reinterpret_cast<pointer>( pool.GetNext(alignof(T), size_request, n) );
                if ( return_value == 0 ) {
                    throw std::bad_alloc();
                }
                return return_value;
            }
            
            // size_request was too big to be handled by an overflow pool.
            // I think this happends when we reserve a vector that is already allocated for example.
            // Or when there is a design error, in that the allocator using this pool is currently not the only user
            throw std::bad_alloc();
        };
        
        void
        deallocate(pointer p, size_type n) noexcept
        {
            pool.Free(p, n*sizeof(T), n);
        };
        
        template <class U, class... Args>
        void
        construct(U* p, Args&&... args)
        { ::new(reinterpret_cast<void*>(p)) U(std::forward<Args>(args)...); }

        void
        destroy(pointer p)
        { p->~T(); };
        
    private:
        Pool & pool;
    };
    
    // Return true if allocators b and a can be safely interchanged. "Safely interchanged" means that b could be
    // used to deallocate storage obtained through a and vice versa.
    template <class T, class T2>
    bool
    operator == ( const StackAllocator<T>& a, const StackAllocator<T2>& b)
    { return true; }; // because they use the same pool.
    
    template <class T, class T2>
    bool
    operator != ( const StackAllocator<T>& a, const StackAllocator<T2>& b)
    { return !(operator==(a,b)); };

} // NS
