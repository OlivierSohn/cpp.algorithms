#pragma once

#pragma push_macro( "new" )
#undef new
#include <new>

namespace imajuscule {
    
    template <class T1> class Allocator;
    
    // Specialize for void
    template <> class Allocator<void>
    {
    public:
        typedef void * pointer;
        typedef const void* const_pointer;
        typedef void value_type;
        template <class U1> struct rebind { typedef Allocator<U1> other; };
    };
    
    // todo policy : either it uses the global pool, typically for a container
    // on ths stack, or it uses its own pool, for long lived containers
    template <class T1> class Allocator
    {
    public:
        typedef T1 value_type;
        typedef size_t size_type;
        typedef ptrdiff_t difference_type;
        typedef T1* pointer;
        typedef const T1* const_pointer;
        typedef T1& reference;
        typedef const T1& const_reference;
        
        // The rebind member allows a container to construct an allocator for some arbitrary type out of
        // the allocator type provided as a template parameter.
        template <class U1> struct rebind { typedef Allocator<U1> other; };
        
        // Constructors
        Allocator() : pool(Pool::getInstance()) {};
        Allocator( const Allocator& other ) : Allocator() {};
        template <class U1> Allocator(const Allocator<U1>&) : Allocator() {};
        
        // Destructor
        ~Allocator( void ) {};
        
        // Returns the address of r as a pointer type. This function and the following function are used
        // to convert references to pointers.
        pointer address(reference r) const { return &r; };
        const_pointer address(const_reference r) const { return &r; };
        
        // Allocate storage for n values of T1.
        pointer allocate( size_type n, Allocator<void>::const_pointer hint = 0 )
        {
            auto const size_request = n*sizeof(T1);
            pointer return_value = reinterpret_cast<pointer>( pool.GetNext(alignof(T1), size_request, n) );
            if ( return_value != 0 ) {
                return return_value;
            }
            if(pool.empty()) {
                pool.resize(2*size_request);
                return_value = reinterpret_cast<pointer>( pool.GetNext(alignof(T1), size_request, n) );
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
        
        // Deallocate storage obtained by a call to allocate.
        void deallocate(pointer p, size_type n)
        {
            pool.Free(p, n);
        };
        
        // Return the largest possible storage available through a call to allocate.
        size_type max_size() const
        {
            size_type return_value = 0xFFFFFFFF;
            return_value /= sizeof(T1);
            return return_value;
        };
        
        // Construct an object of type T1 at the location of ptr
        void construct(pointer ptr)
        {
            ::new (reinterpret_cast<void*>(ptr)) T1;
        };
        
        // Construct an object of type T1 at the location of ptr, using the value of U1 in the call to the
        // constructor for T1.
        template <class U1> void construct(pointer ptr, const U1& val)
        {
            ::new (reinterpret_cast<void*>(ptr)) T1(val);
        };
        
        // Construct an object of type T1 at the location of ptr, using the value of T1 in the call to the
        // constructor for T1.
        void construct(pointer ptr, const T1& val)
        {
            ::new (reinterpret_cast<void*>(ptr)) T1(val);
        };
        
        // Call the destructor on the value pointed to by p
        void destroy(pointer p)
        {
            p->T1::~T1();
        };
    private:
        Pool & pool;
    };
    
    // Return true if allocators b and a can be safely interchanged. "Safely interchanged" means that b could be
    // used to deallocate storage obtained through a and vice versa.
    template <class T1, class T2> bool operator == ( const Allocator<T1>& a, const Allocator<T2>& b)
    {
        return true;
    };
    // Return false if allocators b and a can be safely interchanged. "Safely interchanged" means that b could be
    // used to deallocate storage obtained through a and vice versa.
    template <class T1, class T2> bool operator != ( const Allocator<T1>& a, const Allocator<T2>& b)
    {
        return false;
    };

} // NS

#pragma pop_macro( "new" )
